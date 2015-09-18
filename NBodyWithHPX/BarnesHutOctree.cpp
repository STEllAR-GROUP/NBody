
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <hpx/lcos/wait_each.hpp>
#include <hpx/util/unwrapped.hpp>
#include <hpx/lcos/local/spinlock.hpp>
#include <hpx/exception.hpp>

#include <boost/timer.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

//#include <mutex>
//#include <boost/thread.hpp>

using hpx::naming::id_type;
using hpx::naming::invalid_id;

using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::lcos::wait_each;
using hpx::async;

using hpx::util::high_resolution_timer;

using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

using hpx::init;
using hpx::finalize;
using hpx::find_here;

using hpx::cout;
using hpx::flush;

using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;



//Parameters for computing force and well seperated cells for each node
#define th 60

//Data for each node
struct Body {

public:
    int ID1,parent;
    double m1;
    std::vector<double> r1, v1, force; };
std::vector<Body> body;

void initialize_body(int N, Body body1){
    std::vector<double> A4(3,0);

    body1.ID1=0; body1.m1=0; body1.parent=0;
    body1.r1=A4; body1.v1=A4; body1.force=A4;

    for(int i=0; i<N; ++i)
        body.push_back(body1);}

struct Cell {

public:
    int ID2, parent2,NumNodes=0,level,Ncell;
    std::vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2;
    std::vector<double> r2, rd, boundary;

};
std::vector<Cell> cell(100000);


void initialize_cell(int N){
    std::vector<int> A2(6,-1); std::vector<double> A4(3,0),A5(6,0);
    for(int i=0; i<cell.size(); ++i){
        cell[i].parent2=0; cell[i].NumNodes=0;
        cell[i].neighbors=A2; cell[i].r2=A4; cell[i].rd=A4;
    }
    cell[0].ID2=0; cell[0].NumNodes=N; cell[0].level=0; cell[0].boundary=A5;
    for (int i = 0; i < N; ++i)
        cell[0].members.push_back(i);
}

//Read input from txt file
void read_Input(){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("/home/zahra/Desktop/Inputs/ex10000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        body[i].ID1=atoi(temp.c_str()); body[i].m1=stof(temp1,&sz);
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}
        body[i].parent=0;}}

//Detremining max and min in the boundary cells
void apply_changes(int n){
    for(int j=0; j<cell[n].NumNodes; ++j){
        int i=cell[n].members[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//Finding out if the node is in that subcubic
bool InCube(int node, int p){
    for (int i=0; i<3; ++i)
        if(body[node].r1[i]<cell[p].boundary[2*i] || body[node].r1[i]>cell[p].boundary[2*i+1])
            return false;
    return true;}

//Determining the boundry of the each subcube
void det_boundary_subcube(int n){
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    int p=n*8;

    a1=cell[n].boundary[0]; a2=(cell[n].boundary[0]+cell[n].boundary[1])/2; a3=cell[n].boundary[1];
    b1=cell[n].boundary[2]; b2=(cell[n].boundary[2]+cell[n].boundary[3])/2; b3=cell[n].boundary[3];
    c1=cell[n].boundary[4]; c2=(cell[n].boundary[4]+cell[n].boundary[5])/2; c3=cell[n].boundary[5];

    std::vector<double> A1, A2, A3, A4, A5, A6, A7, A8; std::vector<double> B1(3,0); std::vector<int> B2(6,-1);
    A1.push_back(a1); A1.push_back(a2); A1.push_back(b1); A1.push_back(b2); A1.push_back(c1); A1.push_back(c2);
    A2.push_back(a1); A2.push_back(a2); A2.push_back(b1); A2.push_back(b2); A2.push_back(c2); A2.push_back(c3);
    A3.push_back(a1); A3.push_back(a2); A3.push_back(b2); A3.push_back(b3); A3.push_back(c2); A3.push_back(c3);
    A4.push_back(a1); A4.push_back(a2); A4.push_back(b2); A4.push_back(b3); A4.push_back(c2); A4.push_back(c3);
    A5.push_back(a2); A5.push_back(a3); A5.push_back(b1); A5.push_back(b2); A5.push_back(c1); A5.push_back(c2);
    A6.push_back(a2); A6.push_back(a3); A6.push_back(b1); A6.push_back(b2); A6.push_back(c2); A6.push_back(c3);
    A7.push_back(a2); A7.push_back(a3); A7.push_back(b2); A7.push_back(b3); A7.push_back(c1); A7.push_back(c2);
    A8.push_back(a2); A8.push_back(a3); A8.push_back(b2); A8.push_back(b3); A8.push_back(c2); A8.push_back(c3);

    cell[p+1].boundary=A1; cell[p+2].boundary=A2; cell[p+3].boundary=A3; cell[p+4].boundary=A4;
    cell[p+5].boundary=A5; cell[p+6].boundary=A6; cell[p+7].boundary=A7; cell[p+8].boundary=A8;}



//Computing Center of mass and position of cell
void root(int n){
    std::vector<float> CM(3,0); cell[n].m2=0;
    for (int j=0; j<cell[n].NumNodes; ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


/////////////////////////////////////////////////////////////////// Parellel Octree

hpx::future <void> A(int n) {

    for (int i = 1; i < 9; ++i) {
        int p = n * 8; std::vector<int> mem(cell[n].NumNodes,0);
        cell[p + i].members=mem;
        cell[p + i].ID2 = p + i;
        cell[p + i].level = cell[n].level + 1;
        cell[p + i].parent2 = n;
        cell[p + i].NumNodes=0;

        for (int j = 0; j < cell[n].NumNodes; ++j) {
            int temp=cell[n].members[j];
            if (InCube(temp, p + i) && body[temp].parent==n) {
                cell[p + i].members[cell[p + i].NumNodes] = body[temp].ID1;
                body[temp].parent=p+i;
                cell[p + i].NumNodes = cell[p + i].NumNodes + 1;}}

        if (cell[p + i].NumNodes > 1) root(p + i);

        if (cell[p + i].NumNodes >= 1 && cell[p + i].NumNodes <= th) {
            for (int j = 0; j < cell[p + i].NumNodes; ++j) {
                cell[n].child.emplace_back(cell[p + i].members[j]);}}

        if (cell[p + i].NumNodes > th) {
            cell[n].scell.emplace_back(p + i);
            det_boundary_subcube(p + i);}}}


hpx::future<void> strc3(int n){

    hpx::future<void> r=async(A,n);
    
    return r.then([=](hpx::future<void> f){
        std::vector<hpx::future<void>> v;
        v.reserve(9);
        
        for(int i=1; i<9; ++i){
            int m=n*8+i;
            if(cell[m].NumNodes>th)
                v.push_back(strc3(cell[m].ID2));}

        return hpx::when_all(v);
    });
}

///////////////////////////////////////////////////////////////////


int hpx_main(int argc, char** argv){

    int N;
    string temp;
    fstream textfile;
    textfile.open("/home/zahra/Desktop/Inputs/ex10000.txt");
    textfile >> temp;
    N = (int) atoi(temp.c_str());

    Body body1;
    initialize_body(N, body1);     //Initializing body and cell:
    read_Input();
    initialize_cell(N);
    apply_changes(0);
    root(0);
    det_boundary_subcube(0);


    ///////////////////////////////////////////////////////////

    hpx::future <void> r;
    boost::timer timer;
    hpx::evaluate_active_counters(true, "rate_1");
    r=strc3(0);
    r.wait();
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;

    ////////////////////////////////////////////////////////////

    hpx::evaluate_active_counters(true, "rate_2");
    std::cout<<endl;
    for(int i=0; i<100; ++i)
        std::cout<<cell[i].NumNodes<<","<<cell[i].ID2<<"-";


    hpx::finalize();
    return 0;
}

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}