#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>

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
        for (int j=0; j<9; ++j){
            textfile >> temp;
            a = atoi(temp.c_str());
            init_value[i][j]=a;
    }}
}

//Detremining max and min in the boundary root
void apply_changes(vector<vector<float>> &init_value, vector<vector<float>> &boundary_subcube,vector<vector<float>> members, int count, int ii){
    for(int j=0; j<count; ++j){
        int i=members[ii][j];
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
    void root(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM, vector<vector<type>> &members, type ii);
    void strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members, int ii,int p);
};

//Computing center of mass for each subcube
template <class type> void octree<type>::root(vector<vector<type>> &init_value, vector<vector<type>> &Center_MM,vector<vector<type>> &members, type ii){
    
    int j,i;
    vector<float> CM(3,0);
    
    for (j=0; j<Center_MM[ii][8]; ++j){
        Center_MM[ii][0]=Center_MM[ii][0]+init_value[members[ii][j]][1];
        for (i=0; i<3; ++i)
            CM[i]=CM[i]+init_value[members[ii][j]][i+2] * init_value[members[ii][j]][1];}
    
    for (i=0; i<3;++i)
        Center_MM[ii][i+1]=CM[i]/Center_MM[ii][0];
}

//Travesing tree and dividing it if there is more that one node in the subcube
template <class type> void octree<type>::strc(vector<vector<type>> &init_value,vector<vector<type>> &boundary_subcube,vector<vector<type>> &Center_MM, vector<vector<type>> &members, int ii , int p){
    
    octree<float> tree;
    vector <float> J;
    int i,j,k, z;
    
    for (i = 0; i < Center_MM[ii][8]; ++i)
        init_value[members[ii][i]][9]=0;
    
    for (i=1; i<9; ++i){ //8 is for each subcube
        cout<<"---------"<<"\n"; z=0;
        members.push_back(J); p_push(Center_MM,J,8,p+i,0); Center_MM.at(i+p).push_back(0);
        members.push_back(J); p_push(Center_MM,J,8,p+i,0); Center_MM.at(i+p).push_back(0);
        int temp=Center_MM[ii][8];
        for (j=0; j<temp; ++j){
            k=members[ii][j];
            if(init_value[k][9]==0){
                if(InTheCubic(init_value,k,boundary_subcube,i+p)==1){
                    members.at(i+p).push_back(k);
                    cout<<members[i+p][z]<<"\n"; z=z+1;
                    Center_MM[i+p][8]=Center_MM[i+p][8]+1;
                    init_value[k][8]=ii; init_value[k][9]=1;
    }}}}
    
    for(i=1;i<9;++i){
        if(Center_MM[p+i][8]>=2){
            tree.root(init_value,Center_MM,members,p+i);
            radius(Center_MM,boundary_subcube,p+i);
            det_boundary_subcube(boundary_subcube,p+i,p+8);
            tree.strc(init_value,boundary_subcube,Center_MM, members,i+p,p+8);
    }}
}


int main(int argc, const char * argv[]){
    
    int N,ii=0,j;
    octree<float> tree;
    vector<vector<float>> init_value,Center_MM,boundary_subcube,members;
    vector <float> J;
    
    read_Input(init_value);
    N=(int)init_value.size();
    p_push(boundary_subcube,J,6,0,1);
    p_push(Center_MM,J,8,0,0);
    Center_MM.at(0).push_back(N);
    Center_MM.at(0).push_back(ii);
    members.push_back(J);
    int temp=Center_MM[0][8];
    for (j=0; j<temp; ++j){
        members.at(ii).push_back(init_value[j][0]);
        init_value[j][8]=ii;
    }
    
    apply_changes(init_value,boundary_subcube,members,N,0);
    radius(Center_MM,boundary_subcube,ii);
    tree.root(init_value,Center_MM,members,ii);
    det_boundary_subcube(boundary_subcube,ii,0);
    tree.strc(init_value,boundary_subcube, Center_MM, members, ii, 0);

    return 0;
}
