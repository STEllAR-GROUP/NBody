
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
#define th 10
#define th1 20

//Data for each node
struct Body {

public:
    int ID1,parent,Hilbert_dis,num_A;
    double m1;
    std::vector<double> r1, v1, force;

    /* template<typename Ar> void serialize(Ar &ar, unsigned){
         ar &ID1,&parent,&Hilbert_dis,&num_A,& m1,& r1,& v1,&force;
     }*/

};

std::vector<Body> body;

//Data for each cell
struct Cell {

public:
    int ID2, parent2,NumNodes=0,level,Ncell;
    vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2=0;
    std::vector<double> r2, rd, boundary;

    /* template<typename Ar> void serialize(Ar &ar, unsigned){
         ar & ID2, & parent2,& NumNodes,&level,&neighbors,&A,&Ncell,& members, &child, &scell, &list_cell1, &list_cell2,& r2,& rd,&m2,&boundary;
     }*/
};

//std::vector<Cell> cell;
std::vector<int> cell_arrange;

//Read input from txt file
void read_Input(vector<int> &A){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("/home/zahra/Desktop/Inputs/ex1000000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        body[i].ID1=atoi(temp.c_str()); body[i].m1=stof(temp1,&sz);
        body[i].num_A=i; body[i].parent=0; A.push_back(body[i].ID1);
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}}}

/////////////////////////////////////////////////// HilbertCurve


vector<int> to_binary(unsigned int coord){
    long rem,sum=0;
    int i=1; vector<int> y;
    while(coord>0){
        rem=coord%2;
        sum=sum+(i*rem);
        coord=coord/2;
        i=i*10;}

    while(sum!=0){
        y.push_back(sum%10);
        sum/=10;}

    return y;
}

vector<int> Hilbert_distance(vector<int> &X,vector<int> &Y,vector<int> &Z, int m)
{
int j; vector<int> H;
j=m;
while(j>0){
j=j-1;
if(j<Z.size())
H.push_back(Z[j]);
else
H.push_back(0);
if(j<Y.size())
H.push_back(Y[j]);
else
H.push_back(0);
if(j<X.size())
H.push_back(X[j]);
else
H.push_back(0);}
return H;
}

void merge(vector<int> &A,int, int , int ,int);
void sort_merge(vector<int> &A,int low,int high, int N){
int mid;
if(low<high){
mid=(low+high)/2;
sort_merge(A,low,mid,N);
sort_merge(A,mid+1,high,N);
merge(A,low,mid,high,N);
}
}
void merge(vector<int> &A,int low, int mid, int high, int N){
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


//////////////////////////////////////////////////

//Detremining max and min in the boundary cells
void apply_changes(int n,vector<Cell> &cell){
    for(int j=0; j<cell[n].members.size(); ++j){
        int i=cell[n].members[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//Finding out if the node is in that subcubic
bool InCube(int node, int p,vector<Cell> &cell){
    for (int i=0; i<3; ++i)
        if(body[node].r1[i]<cell[p].boundary[2*i] || body[node].r1[i]>cell[p].boundary[2*i+1])
            return false;
    return true;}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n,vector<Cell> &cell);
    void strc(int n,vector<Cell> &cell);
};

//Computing Center of mass and position of cell
template <class type> void octree<type>::root(int n,vector<Cell> &cell){
    vector<float> CM(3,0);
    for (int j=0; j<cell[n].members.size(); ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


//////////////////////////////////////////////////////////////////

vector<double> members(int n,vector<Cell> &cell){

    int p=n*8;
    double a1[9],b1[9],c1[9]; vector<double> r;
    octree<float> tree;
    for(int i=1; i<9; ++i){
        cell[i+p].level=cell[n].level+1; cell[i+p].ID2=i+p; cell[n].scell.push_back(i+p);
        cell[i+p].parent2=cell[n].ID2; cell[i+p].NumNodes=(int)ceil((float)((cell[n].members.size())/8));

        int start=(int)ceil((i-1)*((float)(cell[n].members.size())/8)), en=(int)ceil((i)*((float)cell[n].members.size()/8)); a1[i]=0, b1[i]=0, c1[i]=0;
        if(start<en){
            for(int k=start; k<en; ++k){
                cell[i+p].members.push_back(cell[n].members[k]);
                body[cell[n].members[k]].parent=i+p;}}

        for(int j=0; j<3; ++j){
            cell[i+p].r2.push_back(0);
            cell[i+p].rd.push_back(0);}

        for(int j=0; j<6; ++j)
            cell[i+p].boundary.push_back(0);

        apply_changes(i+p,cell);
        tree.root(i+p,cell);
        a1[i]=abs(cell[i+p].boundary[1]-cell[i+p].boundary[0]); r.push_back(a1[i]);
        b1[i]=abs(cell[i+p].boundary[3]-cell[i+p].boundary[2]); r.push_back(b1[i]);
        c1[i]=abs(cell[i+p].boundary[5]-cell[i+p].boundary[4]); r.push_back(c1[i]);


        if(cell[i+p].members.size()>=1){

            if(a1[i] <= th && b1[i] <= th && c1[i] <= th){
                if(cell[i+p].members.size()==1){
                    cell[n].child.push_back(body[cell[i+p].members[0]].ID1);
                    cell[i+p].child.push_back(body[cell[i+p].members[0]].ID1);
                    body[cell[i+p].members[0]].parent=n;}
                else{
                    for(int j=0; j<cell[i+p].members.size(); ++j)
                        cell[i+p].child.push_back(body[cell[i+p].members[j]].ID1);}}

            else if(cell[i+p].members.size()>=2 && a1[i] > th || b1[i] > th || c1[i] > th)
                    cell[i+p].parent2=n;}}

    return r;
}



hpx::future<void> strc3(int n,vector<Cell> &cell){

    int i,p=n*8;
    double a1[9],b1[9],c1[9];
    vector<double> r;

    r=members(n,cell);
    std::vector<hpx::future<void>> v;
    v.reserve(8);

    if(cell[n].level==0 || cell[n].level==1 || cell[n].level==2 || cell[n].level==3 || cell[n].level==4 || cell[n].level==5){
        for(i=1; i<9; ++i){
            a1[i]=r[(i-1)*3]; b1[i]= r[(i-1)*3+1]; c1[i]=r[(i-1)*3+2];
            if(cell[p+i].members.size()>th1){
                if(a1[i] > th || b1[i] > th || c1[i] > th)
                    v.push_back(strc3(p+i,cell));}}}
    else{
        for(i=1; i<9; ++i){
            a1[i]=r[(i-1)*3]; b1[i]= r[(i-1)*3+1]; c1[i]=r[(i-1)*3+2];
            if(cell[p+i].members.size()>th1){
                if(a1[i] > th || b1[i] > th || c1[i] > th)
                    strc3(cell[p+i].ID2,cell);}}}

    return when_all(v);
}

//////////////////////////////////////////////////////////////////

int hpx_main(int argc, char** argv){

    octree<float> tree;
    int N;
    string temp;
    fstream textfile;
    vector<int> A,H; int m;
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
        body1.force.push_back(0);}

    for(int i=0; i<N; ++i)
        body.push_back(body1);

    read_Input(A);
    vector<int>  X, Y, Z;
    for(int i=0; i<N; ++i){
        int sum=0;
        X=to_binary(body[i].r1[0]);
        Y=to_binary(body[i].r1[1]);
        Z=to_binary(body[i].r1[2]);

        m=(int)Z.size();
        if(X.size()>m)
            m=(int)X.size();
        if(Y.size()>m)
            m=(int)Y.size();

        H=Hilbert_distance(X,Y,Z,m);
        for(int j=0; j<H.size(); ++j)
            sum=sum+H[j]*pow(2,j);
        body[i].Hilbert_dis=sum;}
    sort_merge(A,0,N,N);
    vector<Cell> cell;
    Cell cell1,cell2;
    cell_arrange.push_back(0);
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
        cell1.members.push_back(body[A[i]].ID1);

    cell.push_back(cell1);
    apply_changes(0,cell);
    tree.root(0,cell);
    for(int i=0; i<8; ++i)
        cell.push_back(cell2);

    for(int k=cell.size(); k<1000000; ++k)
        cell.push_back(cell2);

    hpx::evaluate_active_counters(true, "rate_1");
    ///////////////////////////////////////////////////////
    future <void> r;
    boost::timer timer;
    r=strc3(0,cell);
    r.wait();
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;
    hpx::evaluate_active_counters(true, "rate_2");
    ///////////////////////////////////////////////////////

    /*for(int i=0; i<100; ++i)
        std::cout<<cell[i].members.size()<<",";*/

    std::cout<<endl;
    int max=0;
    for(int i=0; i<cell.size(); ++i)
        if(cell[i].level>max)
            max=cell[i].level;
    std::cout<<max<<endl;

    hpx::finalize();
    return 0;
}

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}
