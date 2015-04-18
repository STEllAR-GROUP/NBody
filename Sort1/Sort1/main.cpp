//
//  main.cpp
//  Struct1
//
//  Created by Zahra Khatami on 4/10/15.
//  Copyright (c) 2015 Zahra Khatami. All rights reserved.


#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>
#include <boost/timer/timer.hpp>
using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;

//Parameters for computing force and well seperated cells for each node
struct Parameters{
    int th_A=30, it=4;
    double time_step=0.1;
    double theta=0.1, G1=6.673*pow(10.0,-11.0);};

//Data for each node
struct Body {
    int ID1,parent,Hilbert_dis,num_A;
    double m1,v1[3],force[3]={0,0,0};
    unsigned int r1[3];};

//Data for each cell
struct Cell {
    int ID2, parent2,NumNodes=0,level,neighbors[6]={-1,-1,-1,-1,-1,-1},Ncell;
    vector<int> child, scell, list_cell1, list_cell2, A;
    double r2[3],rd=0,m2=0,boundary[6]={4444,0,4444,0,4444,0};};

//Read input from txt file
void read_Input(Body body[]){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("ex100000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        body[i].ID1=atoi(temp.c_str());
        body[i].m1=stof(temp1,&sz); body[i].num_A=i;
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}
        body[i].parent=0;}}

//Detremining max and min in the boundary cells
void apply_changes(vector<Cell> &cell, Body body[],int n){
    //cout<<"*"<<cell[n].A.size()<<"*";
    for(int j=0; j<cell[n].A.size(); ++j){
        int i=cell[n].A[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//Determining if node moves to this cell or not
bool InCube(Parameters params,Body body[],vector<Cell> &cell,int node, int newcell){
    int a1=0,b1=0,c1=0;
    if(cell[newcell].boundary[1]<body[node].r1[0]) a1=body[node].r1[0]-cell[newcell].boundary[0];
    if(cell[newcell].boundary[0]>body[node].r1[0]) a1=cell[newcell].boundary[1]-body[node].r1[0];
    if(cell[newcell].boundary[3]<body[node].r1[1]) b1=body[node].r1[1]-cell[newcell].boundary[2];
    if(cell[newcell].boundary[2]>body[node].r1[1]) b1=cell[newcell].boundary[3]-body[node].r1[1];
    if(cell[newcell].boundary[5]<body[node].r1[2]) c1=body[node].r1[2]-cell[newcell].boundary[4];
    if(cell[newcell].boundary[4]>body[node].r1[2]) c1=cell[newcell].boundary[5]-body[node].r1[2];
    
    if(a1 <= params.th_A  && b1 <= params.th_A && c1 <= params.th_A)
        return true;
    else
        return false;}

//Hilbert distance
void Hilbert_distance(coord_t* X, int b, int n){
    coord_t M = 1 << (b - 1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1){
        P = Q - 1;
        for (i = 0; i < n; i++)
            if ((X[i] & Q) != 0)
                X[0] ^= P; // invert
            else{
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;}}
    // Gray encode
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if ((X[n - 1] & Q)!=0)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;}

//Hilbert

void Hilbert(Body body[], int N){
    
    for(int i=0; i<N; ++i){
        coord_t X[3]={body[i].r1[0],body[i].r1[1],body[i].r1[2]};
        int m=10;
        int n=m*3; int Y[n];
        Hilbert_distance(X,m, 3); int sum=0;
        int j=n-1, M1=m;
        while(j>=0){
            M1=M1-1;
            for(int k=0; k<3; ++k){
                Y[j]=(X[k]>>M1 & 1);
                j=j-1; }}
        for(int i=0; i<n; ++i)
            sum=sum+Y[i]*pow(2,i);
        body[i].Hilbert_dis=sum;}}


//Sorting
void merge(vector<int> &A,Body body[],int, int , int ,int);
void sort_merge(vector<int> &A,Body body[], int low,int high, int N){
    int mid;
    if(low<high){
        mid=(low+high)/2;
        sort_merge(A,body,low,mid,N);
        sort_merge(A,body,mid+1,high,N);
        merge(A,body,low,mid,high,N);}}

void merge(vector<int> &A,Body body[],int low, int mid, int high, int N){
    int h,i,j,k,B[N];
    i=low; h=low; j=mid+1;
    while((h<=mid) && (j<=high)){
        if((body[A[h]].Hilbert_dis<body[A[j]].Hilbert_dis)){
            B[i]=body[A[h]].ID1;
            body[A[h]].num_A=i;
            h++; }
        else{
            B[i]=body[A[j]].ID1;
            body[A[j]].num_A=i;
            j++;}
        i=i+1;}
    
    if(h>mid)
        for(k=j;k<=high; ++k){
            B[i]=body[A[k]].ID1;
            body[A[k]].num_A=i;
            i=i+1;}
    else
        for(k=h; k<=mid; ++k){
            B[i]=body[A[k]].ID1;
            body[A[k]].num_A=i;
            i=i+1;}
    for(k=low; k<=high; ++k)
        A[k]=B[k];}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(vector<Cell> &cell, Body body[], int n);
    void strc(Parameters params,vector<Cell> &cell, Body body[], vector<int> &cell_arrange,int n);
    void cell_list(Parameters params,vector<Cell> &cell, Body body[]);
    void traverse_tree(Parameters params,vector<Cell> &cell, Body body[], int root,int node, int root2);
    void compute_force(Parameters params,vector<Cell> &cell,Body body[],int parent,int node);
    void new_tree(Parameters params,vector<Cell> &cell,Body body[],int parent,int node);
    void insert_node(Parameters params,vector<Cell> &cell,Body body[],int newcell, int node);
};

//Computing Center of mass and position of cell
template <class type> void octree<type>::root(vector<Cell> &cell, Body body[], int n){
    vector<float> CM(3,0); float RD=0;
    for(int i=0; i<6; ++i)
        RD=RD+cell[n].boundary[i];
    cell[n].rd= RD;
    if(!cell[n].A.empty()){
        for (int j=0; j<cell[n].A.size(); ++j){
            cell[n].m2=cell[n].m2+body[cell[n].A[j]].m1;
            for (int i=0; i<3; ++i)
                CM[i]=CM[i]+(body[cell[n].A[j]].r1[i] * body[cell[n].A[j]].m1);}
        for (int i=0; i<3;++i)
            cell[n].r2[i]=CM[i]/cell[n].m2;}}

//Arranging cells from root to leaf in Octree
template <class type> void octree<type>::strc(Parameters params,vector<Cell> &cell, Body body[], vector<int> &cell_arrange,int n){
    int i,j,p; float a1[9],b1[9],c1[9]; octree<float> tree;
    p=(cell[n].Ncell)*8; Cell cell1;
    
    for(i=1; i<9; ++i){
        cell[i+p].level=cell[n].level+1; cell[i+p].ID2=i+p; cell[n].scell.push_back(i+p);
        cell[i+p].parent2=cell[n].ID2; cell[i+p].NumNodes=(int)ceil((float)((cell[n].NumNodes)/8));
        int k=0;
        for(j=1; j<9; ++j) //neighbors for each cells
            if(i!=j){
                cell[i+p].neighbors[k]=cell[i+p].ID2;
                k=k+1;}
        
        int start=(int)ceil((i-1)*((float)(cell[n].A.size())/8)), en=(int)ceil((i)*((float)cell[n].A.size()/8)); a1[i]=0, b1[i]=0, c1[i]=0;
        if(start<en){
            for(int k=start; k<en; ++k){
                cell[i+p].A.push_back(cell[n].A[k]);
                body[cell[n].A[k]].parent=i+p;}}
        
        apply_changes(cell,body,i+p);
        tree.root(cell,body,i+p);
        
        a1[i]=cell[i+p].boundary[1]-cell[i+p].boundary[0], b1[i]=cell[i+p].boundary[3]-cell[i+p].boundary[2], c1[i]=cell[i+p].boundary[5]-cell[i+p].boundary[4];
            
        if(cell[i+p].A.size()>=1){
            if(a1[i] <= params.th_A && b1[i] <= params.th_A && c1[i] <= params.th_A){
                if(cell[i+p].A.size()==1){
                    cell[n].child.push_back(body[cell[i+p].A[0]].ID1);
                    cell[i+p].child.push_back(body[cell[i+p].A[0]].ID1);
                    body[cell[i+p].A[0]].parent=n;}
                else{
                    for(j=0; j<cell[i+p].A.size(); ++j)
                        cell[i+p].child.push_back(body[cell[i+p].A[j]].ID1);}}
            else{
                cell[i+p].Ncell=(int)cell_arrange.size();
                cell_arrange.push_back(cell[i+p].ID2);
                cell[i+p].parent2=n;
                for(j=1; j<9; ++j)
                    cell.push_back(cell1);}}}
    
    for(i=1; i<9; ++i)
        if(cell[i+p].A.size()>=1){
            if(a1[i] > params.th_A || b1[i] > params.th_A || c1[i] > params.th_A){
                tree.strc(params,cell,body,cell_arrange,cell[i+p].ID2);}}}

//Modifying octree when position of nodes are changed (1)
template <class type> void octree<type>::insert_node(Parameters params,vector<Cell> &cell,Body body[], int newcell, int node){
    octree<float> tree;int k=-1;
    
    if(!cell[newcell].scell.empty()){
        for(int i=0; i<8; ++i)
            if(InCube(params,body,cell,node,cell[newcell].scell[i])==true && k==-1){
                k=1; int temp=cell[newcell].scell[i];
                tree.insert_node(params,cell,body,temp,node);}}

    cell[newcell].NumNodes=cell[newcell].NumNodes+1;
    cell[newcell].A.push_back(node);
    body[node].parent=newcell;}

//Modifying octree when position of nodes are changed (2)
template <class type> void octree<type>::new_tree(Parameters params,vector<Cell> &cell,Body body[],int parent,int node){
    int k=-1; octree<float> tree;
    if (parent!=0){
        for(int j=0; j<7; ++j)
            if(cell[parent].neighbors[j]!=-1 && k==-1)
                if(InCube(params,body,cell,node,cell[parent].neighbors[j])==true){
                    k=1; int temp=cell[parent].neighbors[j];
                    tree.insert_node(params,cell,body,temp,node);}
        
        if(k==-1 && cell[parent].parent2!=0)
            tree.new_tree(params,cell, body,cell[parent].parent2,node);}}

//creating interaction list for each cell (1)
template <class type> void octree<type>::cell_list(Parameters params,vector<Cell> &cell, Body body[]){
    octree<float> tree;
    for (int i=0; i<cell.size(); ++i){
        if(!cell[i].A.empty())
        if(!cell[i].child.empty())//
            tree.traverse_tree(params,cell,body,cell[i].ID2,cell[i].child[0],0);}}

//creating interaction list for each cell (2)
template <class type> void octree<type>::traverse_tree(Parameters params,vector<Cell> &cell, Body body[],int root,int node, int root2){
    int i;
    octree<float> tree;
    float D=pow(pow((body[node].r1[0]-cell[root2].r2[0]),2.0)+pow((body[node].r1[1]-cell[root2].r2[1]),2.0)+pow((body[node].r1[2]-cell[root2].r2[2]),2.0),0.5);
    
    float ratio=D/cell[root2].rd;
    if (ratio<params.theta){
        if(cell[root2].child.size()>=1)
            for(i=0; i<cell[root2].child.size(); ++i){
                cell[root].list_cell1.push_back(cell[root2].child[i]);}
        if(cell[root2].scell.size()>=1)
            for(i=0; i<cell[root2].scell.size(); ++i)
                tree.traverse_tree(params,cell,body,root,node,cell[root2].scell[i]);}
    else
        cell[root].list_cell2.push_back(cell[root2].ID2);}

//computing force, new position and new velocity
template <class type> void octree<type>::compute_force(Parameters params,vector<Cell> &cell,Body body[],int parent,int node){
    vector<float> a(3,0); double t=params.time_step;
    double G=params.G1;
    if(cell[parent].list_cell1.size()!=0)
        for(int j=0; j<3; ++j)
            for (int i=0; i<cell[parent].list_cell1.size(); ++i)
                body[node].force[j]=body[node].force[j]+ (body[cell[parent].list_cell1[i]].r1[j]-body[node].r1[j]) *G* body[cell[parent].list_cell1[i]].m1 * body[node].m1 /pow((1+(pow((body[cell[parent].list_cell1[i]].r1[0]-body[node].r1[0]),2.0)+pow((body[cell[parent].list_cell1[i]].r1[1]-body[node].r1[1]),2.0)+pow((body[cell[parent].list_cell1[i]].r1[2]-body[node].r1[2]),2.0))),1.5);
    if(cell[parent].list_cell2.size()!=0)
        for(int j=0; j<3; ++j)
            for (int i=0; i<cell[parent].list_cell2.size()!=0; ++i)
                body[node].force[j]=body[node].force[j]+ (cell[cell[parent].list_cell2[i]].r2[j]-body[node].r1[j]) *G* cell[cell[parent].list_cell2[i]].m2 * body[node].m1 / pow((1+(pow((cell[cell[parent].list_cell2[i]].r2[0]-body[node].r1[0]),2.0)+pow((cell[cell[parent].list_cell2[i]].r2[1]-body[node].r1[1]),2.0)+pow((cell[cell[parent].list_cell2[i]].r2[2]-body[node].r1[2]),2.0))),1.5);
    for(int j=0; j<3; ++j){
        a[j]=body[node].force[j]/body[node].m1;
        body[node].r1[j]=body[node].r1[j]+body[node].v1[j]*t+0.5*a[j]*t*t;
        body[node].v1[j]=body[node].v1[j]+a[j]*t;}}


int main(int argc, const char * argv[]){
    
    int i;
    vector<int> cell_arrange; octree<float> tree;
    vector<int> m, list_cell1, list_cell2;
    int N; string temp; fstream textfile;
    textfile.open("ex100000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    Parameters params;
    Body *body=new Body[N];
    vector<Cell> cell;
    Cell cell1,cell2;
    read_Input(body);
    cell_arrange.push_back(0);
    cell1.NumNodes=N; cell1.ID2=0; cell1.level=0;
    cell1.Ncell=0;
    for(i=0; i<N; ++i)
        cell1.A.push_back(i);
    cell.push_back(cell1);
    Hilbert(body,N);
    sort_merge(cell[0].A,body,0,N,N);
    apply_changes(cell,body,0);
    for(int j=1; j<9; ++j)
        cell.push_back(cell2);
    
    tree.root(cell,body,0);
    tree.strc(params,cell,body,cell_arrange,0);
    
    boost::timer::cpu_timer timer;
    tree.cell_list(params,cell,body);
    for(int it=0 ; it<params.it; it++)
        for (i=0; i<N; ++i){
            tree.compute_force(params,cell,body,body[i].parent,i);
            if((InCube(params,body,cell,i,body[i].parent)==false)){
                cell[body[i].parent].NumNodes=cell[body[i].parent].NumNodes-1;
                if(InCube(params,body,cell,i,0)==true)
                    tree.new_tree(params,cell,body,body[i].parent,i);}}
    
    boost::timer::cpu_times elapsed = timer.elapsed();
    std::cout << " CPU TIME: " << (elapsed.user + elapsed.system) / 1e9 << " seconds"<< " WALLCLOCK TIME: " << elapsed.wall / 1e9 << " seconds"<< std::endl;
    
  
    
    return 0;
    
}