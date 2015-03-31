#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <boost/timer/timer.hpp>

using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;

//Read input from txt file
void read_Input(vector<vector<float>> &init_value){
    int N;
    float a;
    string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("ex10000.txt");
    textfile >> temp;
    N = atoi(temp.c_str());
    init_value.resize(N);
    for (int i = 0; i < N; ++i)
        init_value[i].resize(10);
    for (int i=0; i<N; ++i){
        textfile >> temp>>temp1;
        init_value[i][0]=atoi(temp.c_str());
        init_value[i][1]=atoi(temp1.c_str());
        for (int j=2; j<10; ++j){
            textfile >> temp;
            a=stof (temp,&sz);
            init_value[i][j]=a;}
        init_value.at(i).push_back(0);//if leaf or not!
    }
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
    
    a1=boundary_subcube[ii][0]; a2=(boundary_subcube[ii][0]+boundary_subcube[ii][1])/2;
    a3=boundary_subcube[ii][1]; b1=boundary_subcube[ii][2];
    b2=(boundary_subcube[ii][2]+boundary_subcube[ii][3])/2; b3=boundary_subcube[ii][3];
    c1=boundary_subcube[ii][4]; c2=(boundary_subcube[ii][4]+boundary_subcube[ii][5])/2;
    c3=boundary_subcube[ii][5];
    
    boundary_subcube[1+p]={a1,a2,b1,b2,c1,c2}; boundary_subcube[2+p]={a1,a2,b1,b2,c2,c3};
    boundary_subcube[3+p]={a1,a2,b2,b3,c1,c2}; boundary_subcube[4+p]={a1,a2,b2,b3,c2,c3};
    boundary_subcube[5+p]={a2,a3,b1,b2,c1,c2}; boundary_subcube[6+p]={a2,a3,b1,b2,c2,c3};
    boundary_subcube[7+p]={a2,a3,b2,b3,c1,c2}; boundary_subcube[8+p]={a2,a3,b2,b3,c2,c3};
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
void p_push(vector<vector<float>> &mat,int it,int i, float cost){
    for(int j=0; j<it;++j)
        mat.at(i).push_back(cost);
}

//Computing data tree structure
template<class type> class octree {
public:
    octree() {};
    void root(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM, vector<vector<type>> &members1, type ii);
    void strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members1,vector<vector<type>> &members2, vector<type> &cells_arange, int ii,int p);
    void cell_list(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int N);
    void traverse_tree(vector<vector<type>> init_value,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int node, int root, int root2);
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
template <class type> void octree<type>::strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members1,vector<vector<type>> &members2,vector<type> &cells_arange, int ii , int p){
    octree<float> tree;
    vector <float> J;
    int i,j,k; long temp;
    if(Center_MM[ii][10]==0){
        Center_MM[ii][10]=1;
        for (i = 0; i < Center_MM[ii][8]; ++i)
            init_value[members1[ii][i]][9]=0;
        for (i=1; i<9; ++i){ //cout<<"\n -------"<<Center_MM[ii][0]<<" \n";
            p_push(Center_MM,11,p+i,0); p_push(boundary_subcube,6,i+p,1);
            Center_MM[p+i][9]=ii; Center_MM[p+i][0]=p+i;
            Center_MM[p+i][10]=0;
            for (j=0; j<Center_MM[ii][8]; ++j){
                k=members1[ii][j];
                if(init_value[k][9]==0 && init_value[k][10]==0){
                    if(InTheCubic(init_value,k,boundary_subcube,i+p)==1){
                        members1.at(i+p).push_back(k); //cout<<" ,"<<k<<", ";
                        Center_MM[i+p][8]=Center_MM[i+p][8]+1;
                        init_value[k][8]=ii; init_value[k][9]=1;}}}}
        
        for(i=1;i<9; ++i){
            if(!members1[i+p].empty()){
                if(Center_MM[p+i][8]>=2){
                    members2.at(ii).push_back(i+p); Center_MM[i+p][9]=ii;
                    cells_arange.push_back(p+i);  //cout<<" -"<<i+p<<"- ";
                    for(j=1; j<9; ++j){
                        members1.push_back(J); members2.push_back(J);
                        Center_MM.push_back(J);boundary_subcube.push_back(J);}
                }
                else if(Center_MM[p+i][8]==1)
                    Center_MM[i+p][9]=ii; init_value[members1[i+p][0]][10]=1;}}
        if(members2[ii].size()>=1){
            for(j=0; j<members2[ii].size(); ++j){
                tree.root(init_value,Center_MM,members1,members2[ii][j]);
                radius(Center_MM,boundary_subcube,members2[ii][j]);
                temp=(long)cells_arange.size();
                for(i=0; i<temp; ++i)
                    if(cells_arange[i]==members2[ii][j])
                        p=8*i;
                det_boundary_subcube(boundary_subcube,members2[ii][j],p);
                tree.strc(init_value,boundary_subcube,Center_MM, members1,members2,cells_arange,members2[ii][j],p);}}}
}

//Traversing the tree
template <class type> void octree<type>::traverse_tree(vector<vector<type>> init_value,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int root, int node, int root2){
    int i,theta=5; ///////Get from input params
    float D=sqrt(pow((init_value[node][2]-Center_MM[root2][2]),2.0)+pow((init_value[node][3]-Center_MM[root2][3]),2.0)+pow((init_value[node][4]-Center_MM[root2][4]),2.0));
    float r=sqrt(pow(Center_MM[root2][5],2.0)+pow(Center_MM[root2][6],2.0)+pow(Center_MM[root2][7],2.0));
    if (D<(r/theta)){
        if (members1[root2].size()!=0)
            for(i=0; i<members1[root2].size(); ++i)
                if(init_value[members1[root2][i]][8]==root2)
                    list_cell1.at(root).push_back(members1[root2][i]);
        if(members2[root2].size()!=0)
            for(i=0; i<members2[root2].size(); ++i)
                traverse_tree(init_value,members1,members2,Center_MM,list_cell1,list_cell2,root,node,members2[root2][i]);}
    else
        list_cell2.at(root).push_back(Center_MM[root2][0]);
}

//Interaction list for each cell
template <class type> void octree<type>::cell_list(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members1,vector<vector<type>> &members2,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,int N){
    int i,j,k; vector <type> J;
    for(i=0; i<Center_MM.size(); ++i){
        list_cell1.push_back(J); list_cell2.push_back(J); k=-1;
        if (!members1[i].empty() && Center_MM[i][8]!=0){
            for (j=0; j<Center_MM[i][8];++j)
                if(init_value[members1[i][j]][8]==i)
                    k=members1[i][j];
            if (k!=-1)
                traverse_tree(init_value,members1,members2,Center_MM,list_cell1,list_cell2,i,k,0);}}
}

//Computing force
template <class type> float octree<type>::compute_force(vector<vector<type>> &init_value,vector<vector<type>> &Center_MM,vector<vector<type>> &list_cell1,vector<vector<type>> &list_cell2,type parent,int node, float force){
    float G=6.673*pow(10.0,-11.0); //////Get it from input
    if(list_cell1[parent].size()!=0)
        for (int i=0; i<list_cell1[parent].size(); ++i)
            force=force+ G* init_value[list_cell1[parent][i]][1]*init_value[node][1]/(1+(pow((init_value[list_cell1[parent][i]][2]-init_value[node][2]),2.0)+pow((init_value[list_cell1[parent][i]][3]-init_value[node][3]),2.0)+pow((init_value[list_cell1[parent][i]][4]-init_value[node][4]),2.0)));
    if(list_cell2[parent].size()!=0)
        for (int j=0; j<list_cell2[parent].size(); ++j)
            force=force+ G*Center_MM[list_cell2[parent][j]][0]*init_value[node][0]/(pow((Center_MM[list_cell2 [parent][j]][1]-init_value[node][2]),2.0)+pow((Center_MM[list_cell2[parent][j]][2]-init_value[node][3]),2.0)+pow((Center_MM[list_cell2[parent][j]][3]-init_value[node][4]),2.0));
    return force;
}

int main(int argc, const char * argv[]){
    boost::timer::cpu_timer timer;
    int N,ii=0,j,i;
    octree<float> tree;
    vector<vector<float>> init_value,Center_MM,boundary_subcube,members1,members2,list_node1,list_node2,list_cell1,list_cell2;
    vector <float> J,cells_arange;
    
    // boost::timer::cpu_timer timer;
    read_Input(init_value);
    //boost::timer::cpu_times elapsed = timer.elapsed();
    // std::cout << " CPU TIME: " << (elapsed.user + elapsed.system) / 1e9 << " seconds"<< " WALLCLOCK TIME: " << elapsed.wall / 1e9 << " seconds"<< std::endl;
    
    
    N=(int)init_value.size();
    float force[N];
    boundary_subcube.push_back(J);
    p_push(boundary_subcube,6,0,1);
    Center_MM.push_back(J);
    p_push(Center_MM,8,0,0);
    Center_MM.at(0).push_back(N);
    Center_MM.at(0).push_back(0);
    members1.push_back(J); members2.push_back(J);
    for(j=0; j<N;++j){
        list_node1.push_back(J); list_node2.push_back(J); force[j]=0;}
    int temp=Center_MM[0][8];
    for (j=0; j<temp; ++j){
        members1.at(0).push_back(init_value[j][0]);
        init_value[j][8]=ii;}
    for(i=1; i<9; ++i){
        members1.push_back(J); members2.push_back(J);
        Center_MM.push_back(J);boundary_subcube.push_back(J); p_push(boundary_subcube,6,i,1);}
    cells_arange.push_back(0);
    apply_changes(init_value,boundary_subcube,members1,N,0);
    radius(Center_MM,boundary_subcube,ii);
    tree.root(init_value,Center_MM,members1,ii);
    det_boundary_subcube(boundary_subcube,ii,0);
    
    
    tree.strc(init_value,boundary_subcube, Center_MM, members1,members2,cells_arange,ii, 0);
    
    
    tree.cell_list(init_value,Center_MM,members1,members2,list_cell1,list_cell2,N);
    
    
    //for (it=0; it<Iteration_times; ++it) =>for iterations times
    for (i=0; i<N; ++i){
        force[i]=tree.compute_force(init_value,Center_MM,list_cell1,list_cell2,init_value[i][8],i,0);
    }
    boost::timer::cpu_times elapsed = timer.elapsed();
    std::cout << " CPU TIME: " << (elapsed.user + elapsed.system) / 1e9 << " seconds"<< " WALLCLOCK TIME: " << elapsed.wall / 1e9 << " seconds"<< std::endl;
    
    return 0;
}

