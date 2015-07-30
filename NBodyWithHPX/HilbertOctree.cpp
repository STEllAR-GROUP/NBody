
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

#include <mutex>
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
#define th 15
#define th1 20

//Data for each node
struct Body {

public:
    int ID1,Hilbert_dis,num_A,parent;
    double m1;
    std::vector<double> r1, v1, force; };
std::vector<Body> body;

void initialize_body(int N, Body body1){
    vector<double> A4(3,0);

    body1.ID1=0; body1.m1=0; body1.parent=0;
    body1.r1=A4; body1.v1=A4; body1.force=A4;

    for(int i=0; i<N; ++i)
        body.push_back(body1);}


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
        body[i].num_A=i; A.push_back(body[i].ID1);
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}}}



struct Cell {

public:
    mutable hpx::lcos::local::spinlock mtx;
    int ID2, parent2,NumNodes=0,level,Ncell;
    vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2;
    vector<double> r2, rd, boundary;

    void apply_changes(int);
};
std::vector<Cell> cell(1000000);


void initialize_cell(int N,vector<int> &A){
    vector<int> A2(6,-1); vector<double> A4(3,0),A5(6,0);

    cell[0].ID2=0; cell[0].parent2=0; cell[0].NumNodes=N; cell[0].level=0; cell[0].Ncell=0;
    cell[0].neighbors=A2; cell[0].r2=A4; cell[0].rd=A4; cell[0].boundary=A5;
    for (int i = 0; i < N; ++i)
        cell[0].members.push_back(body[A[i]].ID1);}

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

    return y; }

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
    return H; }

void merge(vector<int> &A,int, int , int ,int);
void sort_merge(vector<int> &A,int low,int high, int N){
    int mid;
    if(low<high){
        mid=(low+high)/2;
        sort_merge(A,low,mid,N);
        sort_merge(A,mid+1,high,N);
        merge(A,low,mid,high,N);
    }}
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

void create_Hilbert(int N, vector<int> &A){
    vector<int>  X, Y, Z,H; int m;
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
    sort_merge(A,0,N,N);}

//Detremining max and min in the boundary cells
void Cell::apply_changes(int n){
    for(int j=0; j<cell[n].members.size(); ++j){
        int i=cell[n].members[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//creating octree and determining its parameters

void root(int n){
    vector<float> CM(3,0); cell[n].m2=0;
    for (int j=0; j<cell[n].NumNodes; ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


/////////////////////////////////////////////////////////////////// Parellel Octree

vector<double> A(int n){

    int p=n*8; Cell Foo; vector<double> r;
    double a1[9],b1[9],c1[9];

    for(int i=1; i<9; ++i){
        cell[i+p].level=cell[n].level+1; cell[i+p].ID2=i+p; cell[n].scell.push_back(i+p);
        cell[i+p].parent2=cell[n].ID2; cell[i+p].NumNodes=(int)ceil((float)((cell[n].members.size())/8));

        int start=(int)ceil((i-1)*((float)(cell[n].members.size())/8)), en=(int)ceil((i)*((float)cell[n].members.size()/8)); a1[i]=0, b1[i]=0, c1[i]=0;
        if(start<en){
            for(int k=start; k<en; ++k){
                cell[i+p].members.emplace_back(cell[n].members[k]);
                body[cell[n].members[k]].parent=i+p;}}

        for(int j=0; j<3; ++j){
            cell[i+p].r2.emplace_back(0);
            cell[i+p].rd.emplace_back(0);}

        for(int j=0; j<6; ++j)
            cell[i+p].boundary.emplace_back(0);

        Foo.apply_changes(i+p);
        root(i+p);
        a1[i]=abs(cell[i+p].boundary[1]-cell[i+p].boundary[0]); r.push_back(a1[i]);
        b1[i]=abs(cell[i+p].boundary[3]-cell[i+p].boundary[2]); r.push_back(b1[i]);
        c1[i]=abs(cell[i+p].boundary[5]-cell[i+p].boundary[4]); r.push_back(c1[i]);


        if(cell[i+p].members.size()>=1){

            if(a1[i] <= th && b1[i] <= th && c1[i] <= th){
                if(cell[i+p].members.size()==1){
                    cell[n].child.emplace_back(cell[i+p].members[0]);
                    cell[i+p].child.emplace_back(body[cell[i+p].members[0]].ID1);
                    body[cell[i+p].members[0]].parent=n;}
                else{
                    for(int j=0; j<cell[i+p].members.size(); ++j)
                        cell[i+p].child.emplace_back(cell[i+p].members[j]);}}

            else if(cell[i+p].members.size()>=2 && a1[i] > th || b1[i] > th || c1[i] > th)
                    cell[i+p].parent2=n;}}

    return r;
}


hpx::future<void> strc3(int n){

    double a1[9],b1[9],c1[9];
    vector<double> r; int i;


    hpx::future<double> r=async(A,n);

    return r.then([=](hpx::future<void> f){

        std::vector<hpx::future<void>> v;
        v.reserve(9);

        for(int i=1; i<9; ++i){
            int m=n*8+i;
            a1[i] = r[(i - 1) * 3]; b1[i] = r[(i - 1) * 3 + 1]; c1[i] = r[(i - 1) * 3 + 2]; int p = n * 8;
            if(cell[p + i].members.size()>th)
                if (a1[i] > th || b1[i] > th || c1[i] > th)
                    v.push_back(strc3(cell[m].ID2));}

        return hpx::when_all(v);
    });

}

///////////////////////////////////////////////////////////////////


int hpx_main(int argc, char** argv){

    int N; vector<int> A;
    string temp;
    fstream textfile;
    textfile.open("/home/zahra/Desktop/Inputs/ex1000000.txt");
    textfile >> temp;
    N = (int) atoi(temp.c_str());

    Body body1;
    initialize_body(N, body1);     //Initializing body and cell:
    read_Input(A);
    create_Hilbert(N,A);
    Cell Foo;
    initialize_cell(N,A);
    Foo.apply_changes(0);
    root(0);


    ///////////////////////////////////////////////////////////

    boost::timer timer;
    hpx::future<void> r;
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
