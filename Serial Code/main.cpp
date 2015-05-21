
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>

#include <boost/timer.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/timer.hpp>


using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

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

};

std::vector<Body> body;

//Data for each cell
struct Cell {

public:
    int ID2, parent2,NumNodes=0,level,Ncell;
    vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2=0;
    std::vector<double> r2, rd, boundary;

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
    double a1,a2,a3,b1,b2,b3,c1,c2,c3; Cell cell1;

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

    cell1.boundary=A1;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A2;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A3;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A4;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A5;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A6;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A7;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);
    cell1.boundary=A8;
    for(int i=0; i<3; ++i){ cell1.r2.push_back(0); cell1.rd.push_back(0);}
    for(int i=0; i<6; ++i) cell1.neighbors.push_back(-1); cell.push_back(cell1);

}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n);
    void strc(int n);
    void neighbor();
    void cell_list(int N);
    void new_tree(int parent,int node);
    void insert_node(int parent, int node);
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

//Arranging cells from root to leaf in Octree
template <class type> void octree<type>::strc(int n){
    int i,j,p; octree<float> tree; Cell cell1;
    p=(cell[n].Ncell)*8;

    for(i=1; i<9; ++i){
        cell[i+p].ID2=i+p;
        cell[i+p].NumNodes=0;

        for(j=0; j<cell[n].members.size(); ++j){
            if(InCube(cell[n].members[j],i+p)){
                cell[i+p].members.push_back(cell[n].members[j]);
                cell[i+p].NumNodes=cell[i+p].NumNodes+1;}}

        if(cell[i+p].members.size()>1) tree.root(i+p);

        if(cell[i+p].members.size()>=1 && cell[i+p].members.size()<=th){
            for(j=0; j<cell[i+p].members.size(); ++j){
                body[cell[i+p].members[j]].parent=n;
                cell[n].child.push_back(cell[i+p].members[j]);}}

        if(cell[i+p].members.size()>th){
            cell[i+p].Ncell=(int)cell_arrange.size();
            cell_arrange.push_back(i+p);
            cell[n].scell.push_back(i+p);
            cell[i+p].parent2=n;
            det_boundary_subcube(i+p);}}

    for(i=1; i<9; ++i)
        if(cell[i+p].members.size()>th)
            tree.strc(i+p);}

//Modifying octree when position of nodes are changed (1)
template <class type> void octree<type>::insert_node(int parent, int node){
    octree<float> tree;

    cell[parent].NumNodes=cell[parent].NumNodes+1;
    cell[parent].members.push_back(node);
    body[node].parent=parent;

    if(cell[parent].NumNodes<=th)
        cell[parent].child.push_back(node);

    if(cell[parent].NumNodes==th+1){
        cell_arrange.push_back(cell[parent].ID2);
        cell[parent].Ncell=(int)cell_arrange.size();
        //tree.strc(parent);
    }

    if(cell[parent].NumNodes>th){
        int p=cell[parent].Ncell;
        for(int i=1; i<9; ++i){
            if(InCube(node,i+p)==true)
                tree.insert_node(i+p, node);}}
}


//Modifying octree when position of nodes are changed (2)
template <class type> void octree<type>::new_tree(int parent,int node){
    int k=-1; octree<float> tree;

    if (parent!=0){

        for(int i=0; i<6; ++i){
            if(cell[parent].neighbors[i]!=-1 && k==-1){
                if(InCube(node,cell[parent].neighbors[i])){
                    k=1; int temp=cell[parent].neighbors[i];
                    tree.insert_node(temp,node);}}}

        if(k==-1 && cell[parent].parent2!=0)
            tree.new_tree(cell[parent].parent2,node);
    }
}

//Determining neighbors of each cell (left,right,up,down,back,front)
template <class type> void octree<type>::neighbor(){
    for (int i=1; i<cell_arrange.size(); ++i){
        int temp=cell_arrange[i];
        for(int j=1; j<cell_arrange.size(); ++j)
            if(cell[cell_arrange[j]].level==cell[temp].level && cell_arrange[j]!=temp){
                if(cell[cell_arrange[j]].boundary[1]==cell[temp].boundary[0])
                    cell[temp].neighbors[0]=cell[cell_arrange[j]].ID2;
                else if(cell[cell_arrange[j]].boundary[0]==cell[temp].boundary[1])
                    cell[temp].neighbors[1]=cell[cell_arrange[j]].ID2;
                else if(cell[cell_arrange[j]].boundary[4]==cell[temp].boundary[5])
                    cell[temp].neighbors[2]=cell[cell_arrange[j]].ID2;
                else if(cell[cell_arrange[j]].boundary[5]==cell[temp].boundary[4])
                    cell[temp].neighbors[3]=cell[cell_arrange[j]].ID2;
                else if(cell[cell_arrange[j]].boundary[3]==cell[temp].boundary[2])
                    cell[temp].neighbors[4]=cell[cell_arrange[j]].ID2;
                else if(cell[cell_arrange[j]].boundary[2]==cell[temp].boundary[3])
                    cell[temp].neighbors[5]=cell[cell_arrange[j]].ID2;}}}

//creating interaction list for each cell (2)
void traverse_tree(int root,int node, int root2){
    int i;
    float D=sqrt(pow((body[node].r1[0]-cell[root2].r2[0]),2.0)+pow((body[node].r1[1]-cell[root2].r2[1]),2.0)+pow((body[node].r1[2]-cell[root2].r2[2]),2.0));
    float r=sqrt(pow(cell[root2].rd[0],2.0)+pow(cell[root2].rd[1],2.0)+pow(cell[root2].rd[2],2.0));
    float ratio=D/r;
    if (ratio<theta){
        if(cell[root2].child.size()>=1)
            for(i=0; i<cell[root2].child.size(); ++i)
                cell[root].list_cell1.push_back(cell[root2].child[i]);
        if(cell[root2].scell.size()>=1)
            for(i=0; i<cell[root2].scell.size(); ++i)
                traverse_tree(root,node,cell[root2].scell[i]);}
    else
        cell[root].list_cell2.push_back(root2);}

//creating interaction list for each cell (1)
template <class type> void octree<type>::cell_list(int N) {

    octree<float> tree;
    for (int i = 0; i < N; ++i)
        if (cell[body[i].parent].list_cell1.size() == 0 && cell[body[i].parent].list_cell2.size() == 0)
            traverse_tree(body[i].parent, i, 0);
}

//###############################################################################
//computing force, new position and new velocity

double foo5(int node, int j, double a1, double a2, double a3){

    double a; double t=time_step;
    a=a1/a2;
    a3=a3+body[node].v1[j]*t+0.5*a*t*t;

    return a3;}

std::vector<double> compute_step2(int parent, int node, int i, double a1, double a2, double a3, double a6, double a7, double a8, double a9, double a10, double a11){

    std::vector<double> out;

    for(int j=0; j<3; ++j)
        out.push_back((body[node].r1[j],body[cell[parent].list_cell1[i]].r1[j])*(a1*a2*a3)/pow((1+pow((a6-a7),2.0)+pow((a8-a9),2.0)+pow((a10-a11),2.0)),1.5));

    return out;}

std::vector<double> compute_step5(int parent, int node, int i, double a1, double a2, double a3, double a6, double a7, double a8, double a9, double a10, double a11){

    std::vector<double> out;

    for(int j=0; j<3; ++j)
        out.push_back((body[node].r1[j]-cell[cell[parent].list_cell2[i]].r2[j])*(a1*a2*a3)/pow((1+pow((a6-a7),2.0)+pow((a8-a9),2.0)+pow((a10-a11),2.0)),1.5));

    return out;}

void compute2(int node){

    double G=G1; octree<float> tree;

    int parent = body[node].parent;
    int N1 = (int) cell[parent].list_cell1.size();
    int N2 = (int) cell[parent].list_cell2.size();
    std::vector<double> A, B;
    std::vector<double> f1(3, 0), f2(3, 0);
    std::vector<std::vector<double>> solution2, solution1;

    solution1.resize(N1);
    solution2.resize(N2);

    if (N1 >= 1) {

        for(int i=0; i<N1; ++i)
            solution1[i] = compute_step2(body[node].parent, node, i, G, body[node].m1,
                                         body[cell[parent].list_cell1[i]].m1, body[node].r1[0],
                                         body[cell[parent].list_cell1[i]].r1[0],
                                         body[node].r1[1], body[cell[parent].list_cell1[i]].r1[1],
                                         body[node].r1[2], body[cell[parent].list_cell1[i]].r1[2]);

        for (int i = 0; i < N1; ++i) {
            B = solution1[i];
            for (int j = 0; j < 3; ++j)
                f1[j] = f1[j] + B[j];}}

    if (N2 >= 1) {

        for(int i=0; i<N2; ++i)
            solution2[i] = compute_step5(body[node].parent, node, i, G, body[node].m1,
                                         cell[cell[parent].list_cell2[i]].m2, body[node].r1[0],
                                         cell[cell[parent].list_cell2[i]].r2[0],
                                         body[node].r1[1], cell[cell[parent].list_cell2[i]].r2[1],
                                         body[node].r1[2], cell[cell[parent].list_cell2[i]].r2[2]);

        for (int i = 0; i < N2; ++i) {
            A = solution2[i];
            for (int j = 0; j < 3; ++j)
                f2[j] = f2[j] + A[j];}}

    for (int j = 0; j < 3; ++j)
        body[node].force[j] = f2[j] + f1[j];

    for(int i=0; i<3; ++i)
        body[node].r1[i] = foo5(node, i, body[node].force[i], body[node].m1, body[node].r1[i]);

    if(!InCube(node,parent)){
        cell[body[node].parent].NumNodes=cell[body[node].parent].NumNodes-1;
        if(InCube(node,0))
            tree.new_tree(parent,node); }

}

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

    Cell cell1;
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
    det_boundary_subcube(0);

    tree.strc(0);
    tree.neighbor();
    tree.cell_list(N);

    /////////////////////////////////////////////////////


    boost::timer timer;

    int N1=50000;
    
    for (int t = 0; t < it; t++) {
        
        for(int i=0; i<N; ++i)
            compute2(i);
    }


    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;

    /////////////////////////////////////////////////////

    return 0;
}
