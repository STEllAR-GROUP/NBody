
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

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n);
    void strc(int n);
    void cell_list(int N);
    void new_tree(int parent,int node);
    void insert_node(int parent, int node);
};

//Computing Center of mass and position of cell
template <class type> void octree<type>::root(int n){
    vector<float> CM(3,0);
    for (int j=0; j<cell[n].members.size(); ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}

////////////////////////////////////////////////////////////////////////////

template <class type> void octree<type>::strc(int n){
    int i,j,p; float a1[9],b1[9],c1[9]; octree<float> tree;
    using hpx::parallel::for_each;
    using hpx::parallel::par;
    typedef boost::counting_iterator<int> iterator;
    p=n*8; Cell cell1;

    for(i=1; i<9; ++i){
        cell[i+p].level=cell[n].level+1; cell[i+p].ID2=i+p; cell[n].scell.push_back(i+p);
        cell[i+p].parent2=cell[n].ID2; cell[i+p].NumNodes=(int)ceil((float)((cell[n].members.size())/8));

        for(j=1; j<9; ++j) //neighbors for each cells
            if(i!=j){
                cell[i+p].neighbors.push_back(cell[i+p].ID2);}

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

        apply_changes(i+p);
        tree.root(i+p);
        a1[i]=abs(cell[i+p].boundary[1]-cell[i+p].boundary[0]), b1[i]=abs(cell[i+p].boundary[3]-cell[i+p].boundary[2]), c1[i]=abs(cell[i+p].boundary[5]-cell[i+p].boundary[4]);

        if(cell[i+p].members.size()>=1){
            if(a1[i] <= th && b1[i] <= th && c1[i] <= th){
                if(cell[i+p].members.size()==1){
                    cell[n].child.push_back(body[cell[i+p].members[0]].ID1);
                    cell[i+p].child.push_back(body[cell[i+p].members[0]].ID1);
                    body[cell[i+p].members[0]].parent=n;}
                else{
                    for(j=0; j<cell[i+p].members.size(); ++j)
                        cell[i+p].child.push_back(body[cell[i+p].members[j]].ID1);}}
            else if(cell[i+p].members.size()>=2 && a1[i] > th || b1[i] > th || c1[i] > th){
                cell[i+p].parent2=n;}}}

    
    using namespace hpx::parallel;
    for_each(par, iterator(1), iterator(9), [&](int i) {
        if(cell[p+i].members.size()>th1)
            if(a1[i] > th || b1[i] > th || c1[i] > th)
                tree.strc(i+p);
    });
}
int hpx_main(int argc, char** argv){

    //...same as before...
    
    vector<int> A,H; int m;
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
    
    hpx::evaluate_active_counters(true, "rate_1");
    for(int k=cell.size(); k<1000; ++k)
        cell.push_back(cell2);

    hpx::evaluate_active_counters(true, "rate_2");
    boost::timer timer;
    tree.strc(0);
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;
    hpx::evaluate_active_counters(true, "rate_3");
    std::cout<<cell.size()<<endl;

    //...same as before...

    hpx::finalize();
    return 0;
}

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}
