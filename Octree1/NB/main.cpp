#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>

using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;

//Read input from txt file
void read_Input(vector<vector<float>> &init_value){
    int N;
    float a;
    string temp;
    fstream textfile;
    textfile.open("Input.txt");
    textfile >> temp;
    N = atoi(temp.c_str());
    init_value.resize(N);
    for (int i = 0; i < N; ++i)
        init_value[i].resize(10);
    for (int i=0; i<N; ++i){
        for (int j=0; j<10; ++j){
            textfile >> temp;
            a = atoi(temp.c_str());
            init_value[i][j]=a;
        }}
}

//Detremining max and min in the boundary root
void apply_changes(vector<vector<float>> &init_value, vector<vector<float>> &boundary_subcube,vector<vector<float>> members1, int count, int ii){
    for(int j=0; j<count; ++j){
        int i=members1[ii][j];
        if (init_value[i][2]<boundary_subcube[ii][0]) boundary_subcube[ii][0]=init_value[i][2];
        if (init_value[i][2]>boundary_subcube[ii][1]) boundary_subcube[ii][1]=init_value[i][2];
        if (init_value[i][3]<boundary_subcube[ii][2]) boundary_subcube[ii][2]=init_value[i][3];
        if (init_value[i][3]>boundary_subcube[ii][3]) boundary_subcube[ii][3]=init_value[i][3];
        if (init_value[i][4]<boundary_subcube[ii][4]) boundary_subcube[ii][4]=init_value[i][4];
        if (init_value[i][4]>boundary_subcube[ii][5]) boundary_subcube[ii][5]=init_value[i][4];
    }
}

//Determining the boundry of the each subcube
void det_boundary_subcube(vector<vector<float>> &boundary_subcube,int ii, int p){
    float a1,a2,a3,b1,b2,b3,c1,c2,c3;
    vector<float> J(6,0);
    
    a1=boundary_subcube[ii][0]; a2=(boundary_subcube[ii][0]+boundary_subcube[ii][1])/2;
    a3=boundary_subcube[ii][1]; b1=boundary_subcube[ii][2];
    b2=(boundary_subcube[ii][2]+boundary_subcube[ii][3])/2; b3=boundary_subcube[ii][3];
    c1=boundary_subcube[ii][4]; c2=(boundary_subcube[ii][4]+boundary_subcube[ii][5])/2;
    c3=boundary_subcube[ii][5];
    
    boundary_subcube.push_back(J); boundary_subcube[1+p]={a1,a2,b1,b2,c1,c2};
    boundary_subcube.push_back(J); boundary_subcube[2+p]={a1,a2,b1,b2,c2,c3};
    boundary_subcube.push_back(J); boundary_subcube[3+p]={a1,a2,b2,b3,c1,c2};
    boundary_subcube.push_back(J); boundary_subcube[4+p]={a1,a2,b2,b3,c2,c3};
    boundary_subcube.push_back(J); boundary_subcube[5+p]={a2,a3,b1,b2,c1,c2};
    boundary_subcube.push_back(J); boundary_subcube[6+p]={a2,a3,b1,b2,c2,c3};
    boundary_subcube.push_back(J); boundary_subcube[7+p]={a2,a3,b2,b3,c1,c2};
    boundary_subcube.push_back(J); boundary_subcube[8+p]={a2,a3,b2,b3,c2,c3};
}

//Finding out if the node is in that subcubic
bool InTheCubic(vector<vector<float>> &init_value,int k,vector<vector<float>> &boundary_subcube, int ii){
    for (int i=0; i<3; ++i)
        if (init_value[k][i+2]<boundary_subcube[ii][2*i] || init_value[k][i+2]>boundary_subcube[ii][2*i+1])
            return 0;
    return 1;
}

//Raduis of each subcube
void radius(vector<vector<float>> &Center_MM, vector<vector<float>> &boundary_subcube, int ii){
    for (int i=0; i<3; ++i)
        Center_MM[ii][i+4]=boundary_subcube[ii][(2*i+1)]-boundary_subcube[ii][2*i];
}

//Push fnc needed in the code
void p_push(vector<vector<float>> &mat,vector <float> J,int it,int i, int cost){
    mat.push_back(J);
    for(int j=0; j<it;++j)
        mat.at(i).push_back(cost);
}

//Computing data tree structure
template<class type> class octree {
public:
    octree() {};
    void root(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM, vector<vector<type>> &members1, type ii);
    void strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members1,vector<vector<type>> &members2, int ii,int p);
    void cell_list(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int N);
    void traverse_tree(vector<vector<type>> init_value,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int node, int root, int z);
    float compute_force(vector<vector<type>> &init_value,vector<vector<type>> &Center_MM, vector<vector<type>>&list_cell1,vector<vector<type>>&list_cell2,type parent, int node, float force);
};

//Computing center of mass for each subcube
template <class type> void octree<type>::root(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members1, type ii){
    
    int j,i;
    vector<float> CM(3,0);
    for (j=0; j<Center_MM[ii][8]; ++j){
        Center_MM[ii][1]=Center_MM[ii][1]+init_value[members1[ii][j]][1];
        for (i=0; i<3; ++i)
            CM[i]=CM[i]+init_value[members1[ii][j]][i+2] * init_value[members1[ii][j]][1];}
    for (i=0; i<3;++i)
        Center_MM[ii][i+2]=CM[i]/Center_MM[ii][1];
}

//Travesing tree and dividing it if there is more that one node in the subcube
template <class type> void octree<type>::strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members1,vector<vector<type>> &members2, int ii , int p){
    octree<float> tree;
    vector <float> J;
    int i,j,k;
    for (i = 0; i < Center_MM[ii][8]; ++i)
        init_value[members1[ii][i]][9]=0;
    for (i=1; i<9; ++i){
        members1.push_back(J); members2.push_back(J); p_push(Center_MM,J,10,p+i,0);
        Center_MM[i+p][9]=ii; Center_MM[i+p][0]=i+p;
        for (j=0; j<Center_MM[ii][8]; ++j){
            k=members1[ii][j];
            if(init_value[k][9]==0){
                if(InTheCubic(init_value,k,boundary_subcube,i+p)==1){
                    members1.at(i+p).push_back(k);// cout<<k;
                    Center_MM[i+p][8]=Center_MM[i+p][8]+1;
                    init_value[k][8]=ii; init_value[k][9]=1;}}}}
    for(i=1;i<9;++i)
        if(Center_MM[p+i][8]>=2)
            members2.at(ii).push_back(i+p);
    for(j=0; j<members2[ii].size(); ++j){
        p=p+8;
        tree.root(init_value,Center_MM,members1,members2[ii][j]);
        radius(Center_MM,boundary_subcube,members2[ii][j]);
        det_boundary_subcube(boundary_subcube,members2[ii][j],p);
        tree.strc(init_value,boundary_subcube,Center_MM, members1,members2,members2[ii][j],p);
        p=p+7;}
}

//Traversing the tree
template <class type> void octree<type>::traverse_tree(vector<vector<type>> init_value,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int root, int node, int z){
    int i,theta=2; ///////Get from input params
    float D=sqrt(pow((init_value[node][2]-Center_MM[root][2]),2.0)+pow((init_value[node][3]-Center_MM[root][3]),2.0)+pow((init_value[node][4]-Center_MM[root][4]),2.0));
    float r=sqrt(pow(Center_MM[root][5],2.0)+pow(Center_MM[root][6],2.0)+pow(Center_MM[root][7],2.0));
    if (D<(r/theta)){
        for(i=0; i<members1[root].size(); ++i)
            if(init_value[members1[root][i]][8]==root)
                list_cell1.at(z).push_back(members1[root][i]);
        if(members2[root].size()!=0)
            for(i=0; i<members2[root].size(); ++i)
                traverse_tree(init_value,members1,members2,Center_MM,list_cell1,list_cell2,members2[root][i],node,z);}
    else
        list_cell2.at(z).push_back(Center_MM[root][0]);//??????
}

//Interaction list for each cell
template <class type> void octree<type>::cell_list(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int N){
    int i,j,k,z=-1; vector <type> J;
    for(i=0; i<Center_MM.size(); ++i){
        if (!members1[i].empty() && Center_MM[i][8]!=0){
            list_cell1.push_back(J); list_cell2.push_back(J); k=-1; z=z+1;
            list_cell1.at(z).push_back(Center_MM[i][0]); list_cell2.at(z).push_back(Center_MM[i][0]);
            for (j=0; j<Center_MM[i][8];++j)
                if(init_value[members1[i][j]][8]==i)
                    k=members1[i][j];
            if (k!=-1){
                traverse_tree(init_value,members1,members2,Center_MM,list_cell1,list_cell2,i,k,z);
            }}}
}

//Computing force
template <class type> float octree<type>::compute_force(vector<vector<type>> &init_value,vector<vector<type>> &Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,type parent,int node, float force){
    int s1=sizeof(list_cell1[parent]);//////
    int s2=list_cell2[parent][0].size();
    for (int i=0; i<s1; ++i)
        force=force+ 888* init_value[list_cell1[parent][i]][1]*init_value[node][1]/((init_value[list_cell1[parent][i]][2]-init_value[node][2])^2+(init_value[list_cell1[parent][i]][3]-init_value[node][3])^2+(init_value[list_cell1[parent][i]][4]-init_value[node][4])^2);//G=888?????????????
    for (int j=0; j<s2; ++j)
        force=force+ 888*Center_MM[list_cell2[parent][j]][0]*init_value[node][0]/((Center_MM[list_cell2[parent][j]][1]-init_value[node][2])^2+(Center_MM[list_cell2[parent][j]][2]-init_value[node][3])^2+(Center_MM[list_cell2[parent][j]][3]-init_value[node][4])^2);
    return force;
}

int main(int argc, const char * argv[]){
    int N,ii=0,j,i;
    octree<float> tree;
    vector<vector<float>> init_value,Center_MM,boundary_subcube,members1,members2,list_node1,list_node2,list_cell1,list_cell2;
    vector <float> J;
    read_Input(init_value);
    N=(int)init_value.size();
    float force[N];
    p_push(boundary_subcube,J,6,0,1);
    p_push(Center_MM,J,8,0,0);
    Center_MM.at(0).push_back(N);
    Center_MM.at(0).push_back(0);
    members1.push_back(J); members2.push_back(J);
    for(j=0; j<N;++j){
        list_node1.push_back(J); list_node2.push_back(J); force[j]=0;}
    int temp=Center_MM[0][8];
    for (j=0; j<temp; ++j){
        members1.at(ii).push_back(init_value[j][0]);
        init_value[j][8]=ii;}
    
    apply_changes(init_value,boundary_subcube,members1,N,0);
    radius(Center_MM,boundary_subcube,ii);
    tree.root(init_value,Center_MM,members1,ii);
    det_boundary_subcube(boundary_subcube,ii,0);
    tree.strc(init_value,boundary_subcube, Center_MM, members1,members2,ii, 0);
    tree.cell_list(init_value,Center_MM,members1,members2,list_cell1,list_cell2,N);
    /*for(i=0; i<Center_MM.size(); ++i){
     if(!members1[i].empty() ){
     cout<<"\n --------- \n"<<Center_MM[i][0]<<" , "<<Center_MM[i][9] <<"\n";
     for(j=0; j<Center_MM[i][8]; j++)
     cout<<members1[i][j]<<",";
     cout<<"\n";
     for(j=0; j<members2[i].size(); ++j)
     cout<<members2[i][j]<<",";
     }
     }*/
    
    /*for(i=0; i<list_cell1.size(); ++i){
     cout<<"\n -------"<<list_cell1[i][0]<<"\n";
     if(list_cell1[i].size()!=0){
     for(j=1; j<list_cell1[i].size(); ++j)
     cout<<list_cell1[i][j]<<",";}
     cout<<"\n";
     if(list_cell2[i].size()!=0){
     for(j=1; j<list_cell2[i].size();++j)
     cout<<list_cell2[i][j]<<"-";}
     }*/
    
    //computing forces on each node
    //for (it=0; it<Iteration_times; ++it) =>for iterations times
    /*for (i=0; i<N; ++i){
     int p=(int)init_value[i][8];
     force[i]=tree.compute_force(init_value,Center_MM,list_cell1,list_cell2,p,i,0);
     }*/
    
    
    return 0;
}

