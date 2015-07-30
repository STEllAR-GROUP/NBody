
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

#include <boost/timer.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <hpx/exception.hpp>
#include <boost/thread.hpp>

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

using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;



//Parameters for computing force and well seperated cells for each node
#define th 20

//Data for each node
struct Body {

public:
    int ID1,parent;
    double m1;
    std::vector<double> r1, v1, force;

    template<typename Ar> void serialize(Ar &ar, unsigned){
        ar &ID1,&parent,& m1,& r1,& v1,&force;
    }

};

std::vector<Body> body;

//Data for each cell
struct Cell {

public:
    int ID2, parent2,NumNodes=0,level,Ncell;
    vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2=0;
    std::vector<double> r2, rd, boundary;

    template<typename Ar> void serialize(Ar &ar, unsigned){
        ar & ID2, & parent2,& NumNodes,&level,&neighbors,&Ncell,& members, &child, &scell, &list_cell1, &list_cell2,& r2,& rd,&m2,&boundary;
    }
};

//std::vector<Cell> cell;
std::vector<int> cell_arrange;

//Read input from txt file
void read_Input(){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("/home/zahra/Desktop/Inputs/ex1000000.txt");
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
void apply_changes(int n, vector<Cell> &cell){
    for(int j=0; j<cell[n].NumNodes; ++j){
        int i=cell[n].members[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//Finding out if the node is in that subcubic
bool InCube(int node, int p, vector<Cell> &cell){
    for (int i=0; i<3; ++i)
        if(body[node].r1[i]<cell[p].boundary[2*i] || body[node].r1[i]>cell[p].boundary[2*i+1])
            return false;
    return true;}

//Determining the boundry of the each subcube
void det_boundary_subcube(int n, vector<Cell> &cell){
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    int p=n*8;

    a1=cell[n].boundary[0]; a2=(cell[n].boundary[0]+cell[n].boundary[1])/2; a3=cell[n].boundary[1];
    b1=cell[n].boundary[2]; b2=(cell[n].boundary[2]+cell[n].boundary[3])/2; b3=cell[n].boundary[3];
    c1=cell[n].boundary[4]; c2=(cell[n].boundary[4]+cell[n].boundary[5])/2; c3=cell[n].boundary[5];

    std::vector<double> A1, A2, A3, A4, A5, A6, A7, A8; std::vector<int> B2(6,-1); std::vector<double> B1(3,0);
    A1.push_back(a1); A1.push_back(a2); A1.push_back(b1); A1.push_back(b2); A1.push_back(c1); A1.push_back(c2);
    A2.push_back(a1); A2.push_back(a2); A2.push_back(b1); A2.push_back(b2); A2.push_back(c2); A2.push_back(c3);
    A3.push_back(a1); A3.push_back(a2); A3.push_back(b2); A3.push_back(b3); A3.push_back(c2); A3.push_back(c3);
    A4.push_back(a1); A4.push_back(a2); A4.push_back(b2); A4.push_back(b3); A4.push_back(c2); A4.push_back(c3);
    A5.push_back(a2); A5.push_back(a3); A5.push_back(b1); A5.push_back(b2); A5.push_back(c1); A5.push_back(c2);
    A6.push_back(a2); A6.push_back(a3); A6.push_back(b1); A6.push_back(b2); A6.push_back(c2); A6.push_back(c3);
    A7.push_back(a2); A7.push_back(a3); A7.push_back(b2); A7.push_back(b3); A7.push_back(c1); A7.push_back(c2);
    A8.push_back(a2); A8.push_back(a3); A8.push_back(b2); A8.push_back(b3); A8.push_back(c2); A8.push_back(c3);

    for(int j=1; j<9; ++j){ cell[p+j].r2=B1; cell[p+j].rd=B1; cell[p+j].neighbors=B2;}

    cell[p+1].boundary=A1; cell[p+2].boundary=A2; cell[p+3].boundary=A3; cell[p+4].boundary=A4;
    cell[p+5].boundary=A5; cell[p+6].boundary=A6; cell[p+7].boundary=A7; cell[p+8].boundary=A8;}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n, vector<Cell> &cell);
    void neighbor(vector<Cell> &cell);
};

//Computing Center of mass and position of cell
template <class type> void octree<type>::root(int n, vector<Cell> &cell){
    vector<float> CM(3,0);
    for (int j=0; j<cell[n].NumNodes; ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


/////////////////////////////////////////////////////////////////// Parellel Octree


int members(int n,int p, vector<Cell> &cell){
    int k=0;
    Cell cell3;
    cell3.boundary=cell[n].boundary;
    cell3.NumNodes=cell[n].NumNodes;
    cell3.members=cell[n].members;
    for (int j = 0; j < cell3.NumNodes; ++j) {
        if (InCube(cell3.members[j], p, cell)) {
            cell[p].members.push_back(cell3.members[j]);
            k=k+1;}}
    return k;
}



hpx::future <void> A(int n,vector<Cell> &cell) {
    octree<float> tree;
    int p = n*8;

    for (int i = 1; i < 9; ++i) {
        cell[p + i].ID2 = p + i; cell[p + i].level=cell[n].level+1;
        cell[p + i].parent2 = n;

        int B=members(n,p+i,cell);
        cell[p+i].NumNodes=B;

        if (cell[p + i].NumNodes > 1) tree.root(p + i,cell);

        if (cell[p + i].NumNodes >= 1 && cell[p + i].NumNodes <= th) {
            for (int j = 0; j < cell[p + i].NumNodes; ++j) {
                body[cell[p + i].members[j]].parent = n;
                cell[n].child.push_back(cell[p + i].members[j]);
            }}

        if(cell[p + i].NumNodes>th){
            cell[n].scell.push_back(p + i);
            det_boundary_subcube(p + i, cell);}}
}


hpx::future<void> strc3(int n,vector<Cell> &cell){

    int p=n*8;
    //future <void> r=async(A,n,cell);

    //return r.then([&](future<void> f){

    ///////////////////////////////////////////////////
    octree<float> tree;

    for (int i = 1; i < 9; ++i) {
        cell[p + i].ID2 = p + i; cell[p + i].level=cell[n].level+1;
        cell[p + i].parent2 = n;

        int B=members(n,p+i,cell);
        cell[p+i].NumNodes=B;

        if (cell[p + i].NumNodes > 1) tree.root(p + i,cell);

        if (cell[p + i].NumNodes >= 1 && cell[p + i].NumNodes <= th) {
            for (int j = 0; j < cell[p + i].NumNodes; ++j) {
                body[cell[p + i].members[j]].parent = n;
                cell[n].child.push_back(cell[p + i].members[j]);
            }}

        if(cell[p + i].NumNodes>th){
            cell[n].scell.push_back(p + i);
            det_boundary_subcube(p + i, cell);}}

    ///////////////////////////////////////////////////

        std::vector<hpx::future<void>> v;
        v.reserve(8);

        for(int i=1; i<9; ++i)
            if(cell[p + i].NumNodes>th)
                v.push_back(strc3(cell[p + i].ID2, cell));

        return when_all(v);//});
}

///////////////////////////////////////////////////////////////////


//Determining neighbors of each cell (left,right,up,down,back,front)
template <class type> void octree<type>::neighbor(vector<Cell> &cell){
    using hpx::parallel::for_each;
    using hpx::parallel::par;
    typedef boost::counting_iterator<int> iterator;
    using namespace hpx::parallel;
    for_each(par, iterator(1), iterator(100), [&](int i) {
        int temp=cell[i].ID2;
        for(int j=1; j<cell.size(); ++j)
            if(cell[j].members.size()!=0){
            if(cell[cell[j].ID2].level==cell[temp].level && cell[j].ID2!=temp){
                if(cell[cell[j].ID2].boundary[1]==cell[temp].boundary[0])
                    cell[temp].neighbors[0]=cell[cell[j].ID2].ID2;
                else if(cell[cell[j].ID2].boundary[0]==cell[temp].boundary[1])
                    cell[temp].neighbors[1]=cell[cell[j].ID2].ID2;
                else if(cell[cell[j].ID2].boundary[4]==cell[temp].boundary[5])
                    cell[temp].neighbors[2]=cell[cell[j].ID2].ID2;
                else if(cell[cell[j].ID2].boundary[5]==cell[temp].boundary[4])
                    cell[temp].neighbors[3]=cell[cell[j].ID2].ID2;
                else if(cell[cell[j].ID2].boundary[3]==cell[temp].boundary[2])
                    cell[temp].neighbors[4]=cell[cell[j].ID2].ID2;
                else if(cell[cell[j].ID2].boundary[2]==cell[temp].boundary[3])
                    cell[temp].neighbors[5]=cell[cell[j].ID2].ID2;}}});
}


int hpx_main(int argc, char** argv){

    octree<float> tree;
    int N;
    string temp;
    fstream textfile;
    textfile.open("/home/zahra/Desktop/Inputs/ex1000000.txt");

    textfile >> temp;
    N = (int) atoi(temp.c_str());

    Body body1;
    body1.ID1=0;
    body1.m1=0;
    body1.parent=0;

    for(int i=0; i<3; ++i){
        body1.r1.push_back(0);
        body1.v1.push_back(0);
        body1.force.push_back(0);
    }

    for(int i=0; i<N; ++i)
        body.push_back(body1);

    std::vector<Cell> cell;
    Cell cell1,cell2;
    read_Input();

    cell_arrange.push_back(0);
    cell1.NumNodes = N;
    cell1.ID2 = 0;
    cell1.level = 0;
    cell1.Ncell = 0;
    for(int i=0; i<3; ++i){
        cell1.r2.push_back(0);
        cell1.rd.push_back(0);}

    for(int i=0; i<6; ++i){
        cell1.boundary.push_back(0);
        cell1.neighbors.push_back(-1);}

    for (int i = 0; i < N; ++i){
        cell1.members.push_back(i);
        body[i].parent=0;}

    cell1.level=0;
    cell1.NumNodes=N;

    cell.push_back(cell1);
    apply_changes(0,cell);
    tree.root(0,cell);

    std::vector<int> B(N,0);
    //cell2.members=B;
    cell2.NumNodes=0;

    for(int i=0; i<8; ++i)
        cell.push_back(cell2);

    for(int k=cell.size(); k<1000000; ++k)
        cell.push_back(cell2);

    //cell[100].members[40]=100;

    //std::cout<<cell[100].members[40];
    det_boundary_subcube(0,cell);

    ///////////////////////////////////////////////////////////

    future <void> r;
    boost::timer timer;
    hpx::evaluate_active_counters(true, "rate_1");
    r=strc3(0,cell);
    r.wait();
    //tree.neighbor(cell);
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;
    ////////////////////////////////////////////////////////////


    hpx::evaluate_active_counters(true, "rate_2");

    std::cout<<endl;

    for(int i=0; i<100; ++i)
        std::cout<<cell[i].NumNodes<<","<<cell[i].ID2<<"-";

    std::cout<<endl;
    int max=0;
    for(int i=0; i<cell.size(); ++i)
        if(cell[i].level>max)
            max=cell[i].level;
    //std::cout<<max<<endl;*/

    hpx::finalize();
    return 0;
}

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}
