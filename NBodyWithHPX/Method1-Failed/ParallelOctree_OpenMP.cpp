
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>


#include <boost/timer.hpp>
#include <omp.h>


using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;



//Parameters for computing force and well seperated cells for each node
#define th 100
#define it 25
#define time_step 0.1
#define theta 0.5
#define G1 6.673*pow(10.0,-11.0)

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

std::vector<Cell> cell;
std::vector<int> cell_arrange;

//Read input from txt file
void read_Input(){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("/home/zahra/Desktop/Inputs/ex100000.txt");
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

    std::vector<double> A1, A2, A3, A4, A5, A6, A7, A8;
    A1.push_back(a1); A1.push_back(a2); A1.push_back(b1); A1.push_back(b2); A1.push_back(c1); A1.push_back(c2);
    A2.push_back(a1); A2.push_back(a2); A2.push_back(b1); A2.push_back(b2); A2.push_back(c2); A2.push_back(c3);
    A3.push_back(a1); A3.push_back(a2); A3.push_back(b2); A3.push_back(b3); A3.push_back(c2); A3.push_back(c3);
    A4.push_back(a1); A4.push_back(a2); A4.push_back(b2); A4.push_back(b3); A4.push_back(c2); A4.push_back(c3);
    A5.push_back(a2); A5.push_back(a3); A5.push_back(b1); A5.push_back(b2); A5.push_back(c1); A5.push_back(c2);
    A6.push_back(a2); A6.push_back(a3); A6.push_back(b1); A6.push_back(b2); A6.push_back(c2); A6.push_back(c3);
    A7.push_back(a2); A7.push_back(a3); A7.push_back(b2); A7.push_back(b3); A7.push_back(c1); A7.push_back(c2);
    A8.push_back(a2); A8.push_back(a3); A8.push_back(b2); A8.push_back(b3); A8.push_back(c2); A8.push_back(c3);

    for(int j=1; j<9; ++j){
        for(int i=0; i<3; ++i){ cell[p+j].r2.push_back(0); cell[p+j].rd.push_back(0);}
        for(int i=0; i<6; ++i) cell[p+j].neighbors.push_back(-1);}

    cell[p+1].boundary=A1; cell[p+2].boundary=A2; cell[p+3].boundary=A3; cell[p+4].boundary=A4;
    cell[p+5].boundary=A5; cell[p+6].boundary=A6; cell[p+7].boundary=A7; cell[p+8].boundary=A8;


}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n);
    void strc(int n);
    void strc2(int n);
};

//Computing Center of mass and position of cell
template <class type> void octree<type>::root(int n){
    vector<float> CM(3,0);
    for (int j=0; j<cell[n].NumNodes; ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


/////////////////////////////////////////////////////////////////// Parellel Octree

template <class type> void octree<type>::strc2(int n){
    octree<float> tree;
    int p=n*8;
    Cell cell1;
    for(int i=1; i<9; ++i){

        cell[p+i].ID2=p+i;
        for(int j=0; j<cell[n].members.size(); ++j){
            if(InCube(cell[n].members[j],p+i)){
                cell[p+i].members.push_back(cell[n].members[j]);
            }}
        cell[p+i].NumNodes=cell[p+i].members.size();
        if(cell[p+i].members.size()>1) tree.root(p+i);

        if(cell[p+i].members.size()>=1 && cell[p+i].members.size()<=th){
            for(int j=0; j<cell[p+i].members.size(); ++j){
                body[cell[p+i].members[j]].parent=n;
                cell[n].child.push_back(cell[p+i].members[j]);}}

        if(cell[p+i].members.size()>th){
            cell[n].scell.push_back(p+i);
            cell[p+i].parent2=n;
	    det_boundary_subcube(p+i);
        }}


    #pragma omp parallel for schedule(dynamic,10)
    for(int i=1; i<9; ++i){
	if(cell[n*8+i].members.size()>th)
            tree.strc2(cell[n*8+i].ID2);
    }

}

//...same as before...

int main(int argc, char** argv){

    octree<float> tree;
    int N;
    string temp;
    fstream textfile;
    textfile.open("/home/zahra/Desktop/Inputs/ex100000.txt");

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

    for (int i = 0; i < N; ++i)
        cell1.members.push_back(i);

    cell.push_back(cell1);
    apply_changes(0);
    tree.root(0);
    for(int i=0; i<8; ++i)
        cell.push_back(cell2);

    for(int k=cell.size(); k<100000; ++k)
        cell.push_back(cell2);

    boost::timer timer;
    det_boundary_subcube(0);
    tree.strc2(0);
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;

    return 0;
}

