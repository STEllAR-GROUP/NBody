#include <sstream>
#include <vector>  

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/include/async.hpp>

#include <boost/shared_array.hpp>
#include <stack>


typedef unsigned int coord_t;
///////////////////////////////////////////////////////////////////////////////
int n=10000;
int th=10;
int nt=100;
int t=1;     
int th1=100;
double G=6.673*pow(10.0,-11.0);
float theta=0.7;
int k_th=200;

char const* stepper_basename = "/Dist6b/stepper/";
char const* gather_basename = "/Dist6b/gather/";

//--------------------------------------------------
struct Body
{
public:
    int ID1,parent,num_A,Hilbert_dis;
    double r1[3], v1[3], force[3], m1;
    
    friend class hpx::serialization::access;
    template<typename Ar> void serialize(Ar &ar, unsigned){
     ar &ID1 &parent &m1 &r1 &v1 &force;}
}body1;

struct test{
    public:
    int ID;
    double f[3];

    template<typename Ar> void serialize(Ar &ar, unsigned){ 
        ar &ID &f;}
};

std::vector<test> T;

struct Cell
{
public:
    int ID2, parent2,NumNodes=0,level;
    std::vector<int> members, child,scell,neighbors,list_cell1, list_cell2;
    double m2;
    std::vector<double> r2, rd, boundary;

    friend class hpx::serialization::access;  
    template<typename Ar> void serialize(Ar &ar, unsigned){
       ar &ID2 &parent2 &NumNodes &level &members &child &scell &neighbors &list_cell1 &list_cell2 &m2 &r2 &rd &boundary;}
};

std::vector<int> A;
std::vector<Body> b;
std::vector<Cell> cell(10000);

//--------------------------------------------------

template <typename T>
struct partition_allocator
{
private:
    typedef hpx::lcos::local::spinlock mutex_type;

public:
    partition_allocator(std::size_t max_size = std::size_t(-1))
            : max_size_(max_size)
    {
    }

    ~partition_allocator()
    {
        mutex_type::scoped_lock l(mtx_); 
        while (!heap_.empty())
        {
            T* p = heap_.top();
            heap_.pop();
            delete [] p;
        }
    }

    T* allocate(std::size_t n)
    {
        mutex_type::scoped_lock l(mtx_); 
        if (heap_.empty())
            return new T[n];

        T* next = heap_.top();
        heap_.pop();
        return next;
    }

    void deallocate(T* p)
    {
        mutex_type::scoped_lock l(mtx_); 
        if (max_size_ == static_cast<std::size_t>(-1) || heap_.size() < max_size_)
            heap_.push(p);
        else
            delete [] p;
    }

private:
    mutex_type mtx_;
    std::size_t max_size_;
    std::stack<T*> heap_;
};

////////////////////////////////////////////////////////////////////////////////
struct partition_data
{
private:
    typedef hpx::serialization::serialize_buffer<Body> buffer_type;

    struct hold_reference
    {
        hold_reference(buffer_type const& data)
                : data_(data)
        {}

        void operator()(Body*) {}     // no deletion necessary

        buffer_type data_;
    };

    static void deallocate(Body* p)
    {
        alloc_.deallocate(p);
    }

    static partition_allocator<Body> alloc_;

public:
    partition_data()
            : size_(0)
    {}

    partition_data(std::size_t size)
            : data_(alloc_.allocate(size), size, buffer_type::take,
                    &partition_data::deallocate),
              size_(size),
              min_index_(0)
    {}

    partition_data(std::size_t size, int cell_IDs, int flag)
            : data_(alloc_.allocate(size), size, buffer_type::take,
                    &partition_data::deallocate),
              size_(size),
              min_index_(0)
    {
        if(flag==1){
            for(std::size_t j=0; j<size; ++j) {
                data_[j]=b[cell[cell_IDs].members[j]];
                data_[j].ID1=b[cell[cell_IDs].members[j]].ID1;                
                data_[j].force[0]=1; data_[j].force[1]=1;
                data_[j].force[2]=1;
                }
            }
        
        else
          for(std::size_t j=0; j<size; ++j){
                data_[j].ID1=-1;
        }  
  }

    partition_data(partition_data const& base, std::size_t min_index)
            : data_(base.data_.data()+min_index, 1, buffer_type::reference,
                    hold_reference(base.data_)),     
                      size_(base.size()),
              min_index_(min_index)
    {
        HPX_ASSERT(min_index < base.size());
    }

    Body& operator[](std::size_t idx) { return data_[index(idx)]; }
    Body operator[](std::size_t idx) const { return data_[index(idx)]; }

    std::size_t size() const { return size_; }

private:
    std::size_t index(std::size_t idx) const
    {
        HPX_ASSERT(idx >= min_index_ && idx < size_);
        return idx - min_index_;
    }

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & data_ & size_ & min_index_ ;
    }

private:
    buffer_type data_;
    std::size_t size_;
    std::size_t min_index_;
};

partition_allocator<Body> partition_data::alloc_;

///////////////////////////////////////////////////////////////////////////////
inline std::size_t idx(std::size_t i, int dir, std::size_t size)
{

    if(i == 0 && dir == -1) return size-1;
    if(i == size-1 && dir == +1) return 0;

    HPX_ASSERT((i + dir) < size);
        return i + dir;      
}

///////////////////////////////////////////////////////////////////////////////
struct partition_server
        : hpx::components::simple_component_base<partition_server>
{
    enum partition_type
    {
        cell0
    };

    partition_server() {}

    partition_server(partition_data const& data)
            : data_(data)
    {}

    partition_server(std::size_t size, int cell_IDs, int flag)
            : data_(size, cell_IDs,flag)
    {}

    partition_data get_data(partition_type t) const
    {
        switch (t)
        {
            case cell0:
                break;

            default:
                HPX_ASSERT(false);
                break;
        }
        return data_;
    }

    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, get_data, get_data_action);

private:
    partition_data data_;
};

typedef hpx::components::simple_component<partition_server> partition_server_type;
HPX_REGISTER_COMPONENT(partition_server_type, partition_server);

typedef partition_server::get_data_action get_data_action;
HPX_REGISTER_ACTION(get_data_action);

///////////////////////////////////////////////////////////////////////////////
struct partition : hpx::components::client_base<partition, partition_server>
{
    typedef hpx::components::client_base<partition, partition_server> base_type;

    partition() {}

    partition(hpx::id_type where, std::size_t size, int cell_IDs, int flag)
            : base_type(hpx::new_<partition_server>(where,size,cell_IDs,flag))
    {}

    partition(hpx::id_type where, partition_data const& data)
            : base_type(hpx::new_<partition_server>(hpx::colocated(where), data))
    {}

    partition(hpx::future<hpx::id_type> && id)
            : base_type(std::move(id))
    {}

     partition(hpx::future<partition> && c)
            : base_type(std::move(c))
    {}

     hpx::future<partition_data> get_data(partition_server::partition_type t) const
    {
        partition_server::get_data_action act;
        return hpx::async(act, get_gid(), t);
    }
};

///////////////////////////////////////////////////////////////////////////////

struct stepper_server : hpx::components::simple_component_base<stepper_server>
{
    typedef std::vector<partition> space;

    stepper_server() {}

    stepper_server(std::size_t nl)
            : Left_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(), -1, nl))),
              Right_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(), +1, nl))),
              Parent_(2),
              U_(2),
              To_L_(1),
              To_R_(1)
              
    {
    }

    space do_work(std::size_t np, std::size_t nl);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

    void from_R(std::size_t time, partition p) {
        Right_receive_buffer_.store_received(time, std::move(p));}

    void from_L(std::size_t time, partition p){
        Left_receive_buffer_.store_received(time, std::move(p));}

    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_R, from_R_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_L, from_L_action);
   
protected:
    static partition compute_position(partition const &ce0,
                                      partition const &Left, partition const &Right);

    static partition changed_members(partition const &current_members, partition const & From_N,
                                     int cell_ID, int dim, int dir);

    static partition New_Members(partition const& ce0, partition const & From_N, int ID);
    static partition Non_Members(partition const& To_N_total, partition const & To_N);
    static partition copy_mem(partition const& ce0);
           
    
    partition receive_L(std::size_t time) { return Left_receive_buffer_.receive(time);}
    partition receive_R(std::size_t time) { return Right_receive_buffer_.receive(time);}

    void send_L(std::size_t time, partition p){ hpx::apply(from_R_action(), Left_.get(), time, p);}
    void send_R(std::size_t time, partition p){ hpx::apply(from_L_action(), Right_.get(), time, p);}

private:
    hpx::shared_future<hpx::id_type> Left_, Right_; 
    std::vector<space> U_, Parent_,To_L_, To_R_;
    hpx::lcos::local::receive_buffer<partition> Left_receive_buffer_;
    hpx::lcos::local::receive_buffer<partition> Right_receive_buffer_;      
};

typedef hpx::components::simple_component<stepper_server> stepper_server_type;
HPX_REGISTER_COMPONENT(stepper_server_type, stepper_server);

typedef stepper_server::do_work_action do_work_action;
HPX_REGISTER_ACTION(do_work_action);

///////////////////////////////////////////////////////////////////////////////

struct stepper : hpx::components::client_base<stepper, stepper_server>
{
    typedef hpx::components::client_base<stepper, stepper_server> base_type;

    stepper()
            : base_type(hpx::new_<stepper_server>(hpx::find_here(), hpx::get_num_localities_sync()))
    {
        hpx::register_id_with_basename(stepper_basename, get_gid(), hpx::get_locality_id());
    }

    stepper(hpx::future<hpx::id_type> && id)
            : base_type(std::move(id))
    {}

    hpx::future<stepper_server::space> do_work(std::size_t np, std::size_t nl) {
        return hpx::async(do_work_action(), get_gid(), np, nl);
    }
};
//////////////////////////////////////////////////////////////////////////////

void initialize_body(int N)
{    
    body1.ID1=0; body1.m1=0; body1.parent=0;
    for(int i=0; i<3; ++i){
        body1.r1[i]=1.0; body1.v1[i]=1.0; body1.force[i]=0.0;
    }
    for(int i=0; i<N; ++i)
        b.push_back(body1);
}

void read_Input()
{
    int N=n;
    for (int i=0; i<N; ++i){
        b[i].ID1=i; b[i].m1=1; b[i].parent=0;
        b[i].num_A=i; A.push_back(b[i].ID1);
        for(int j=0; j<3; ++j){
            b[i].r1[j]=float(std::rand() %100);
            b[i].v1[j]=0;}}
}

void initialize_cell(int n)
{
    cell[0].ID2=0; cell[0].parent2=0; cell[0].NumNodes=n; cell[0].level=0; 
    cell[0].neighbors.resize(2); cell[0].r2.resize(3); cell[0].rd.resize(3); cell[0].boundary.resize(6);

        for(int i=0; i<3; ++i){ cell[0].r2[i]=0; cell[0].rd[i]=0; cell[0].boundary[2*i]=0; cell[0].boundary[2*i+1]=100;}
    	for(int i=0; i<2; ++i){cell[0].neighbors[i]=-1;}      

    for (int i = 0; i < n; ++i){
        cell[0].members.push_back(b[A[i]].ID1);}
}
//////////////////////////////////////////////////////////////

std::vector<int> to_binary(unsigned int coord){
    long rem,sum=0;
    int i=1; std::vector<int> y;
    while(coord>0){
        rem=coord%2;
        sum=sum+(i*rem);
        coord=coord/2;
        i=i*10;}

    while(sum!=0){
        y.push_back(sum%10);
        sum/=10;}

    return y; }

std::vector<int> Hilbert_distance(std::vector<int> &X,std::vector<int> &Y,std::vector<int> &Z, int m)
{
    int j; std::vector<int> H;
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


void merge(int low, int mid, int high, int N){
    int h,i,j,k;
    std::vector<int> B(N,0);
    i=low; h=low; j=mid+1;

    while((h<=mid) && (j<=high)){
        if((b[A[h]].Hilbert_dis<b[A[j]].Hilbert_dis)){
            B[i]=b[A[h]].ID1;
            b[A[h]].num_A=i;
            h++; }
        else{
            B[i]=b[A[j]].ID1;
            b[A[j]].num_A=i;
            j++;}
        i=i+1;}

    if(h>mid)
        for(k=j;k<=high; ++k){
            B[i]=b[A[k]].ID1;
            b[A[k]].num_A=i;
            i=i+1;}
    else
        for(k=h; k<=mid; ++k){
            B[i]=b[A[k]].ID1;
            b[A[k]].num_A=i;
            i=i+1;}
    for(k=low; k<=high; ++k)
        A[k]=B[k];
}

hpx::future<void> sort_merge(int low,int high, int N){ 
   std::vector<hpx::future<void>> futures;
    int mid;
    if(low<high){  
        mid=(low+high)/2;  
        futures.push_back(hpx::async(sort_merge,low,mid,N));
        futures.push_back(hpx::async(sort_merge,mid+1,high,N)); 
        merge(low,mid,high,N);}

return hpx::when_all(futures);
}

void create_Hilbert(int N){
    std::vector<int>  X, Y, Z,H; int m;
    for(int i=0; i<N; ++i){
        int sum=0;
        X=to_binary(b[i].r1[0]);
        Y=to_binary(b[i].r1[1]);
        Z=to_binary(b[i].r1[2]);

        m=(int)Z.size();
        if(X.size()>m)
            m=(int)X.size();
        if(Y.size()>m)
            m=(int)Y.size();

        H=Hilbert_distance(X,Y,Z,m);
        for(int j=0; j<H.size(); ++j)
            sum=sum+H[j]*pow(2,j);
        b[i].Hilbert_dis=sum;}
   
    hpx::future<void>result;
    result=sort_merge(0,N,N);
    result.wait();
    
}
//////////////////////////////////////////////////////////////

void apply_changes(int n)
{
    for(int j=0; j<cell[n].NumNodes; ++j){
        int i=cell[n].members[j];
        if (b[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=b[i].r1[0];
        if (b[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=b[i].r1[0];
        if (b[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=b[i].r1[1];
        if (b[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=b[i].r1[1];
        if (b[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=b[i].r1[2];
        if (b[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=b[i].r1[2];}
}

bool InCube(int node, int p)
{
    for (int i=0; i<3; ++i)
       if(b[node].r1[i]<cell[p].boundary[2*i] || b[node].r1[i]>cell[p].boundary[2*i+1])
    
       return false;
    return true;
}

void root(int n)
{
    std::vector<float> CM; 
    CM.resize(3); CM[0]=0; CM[1]=0; CM[2]=0;
    cell[n].m2=0;
    for (int j=0; j<cell[n].members.size(); ++j){
        cell[n].m2=cell[n].m2+b[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+b[cell[n].members[j]].r1[i]*b[cell[n].members[j]].m1;
  }
   
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}
}

/////////////////////////////////////////////////////////////////// Parellel Octree

std::vector<double> AA(int n) {

    int p=n*8; std::vector<double> r;
    double a1[9],b1[9],c1[9];

    for(int i=1; i<9; ++i){
        cell[i+p].level=cell[n].level+1; cell[i+p].ID2=i+p; cell[n].scell.push_back(i+p);
        cell[i+p].parent2=cell[n].ID2; cell[i+p].NumNodes=(int)ceil((float)((cell[n].members.size())/8));

        int start=(int)floor((i-1)*((float)(cell[n].members.size())/8)), en=(int)floor((i)*((float)cell[n].members.size()/8)); 
        a1[i]=0, b1[i]=0, c1[i]=0;

        if(start<en){
            for(int k=start; k<en; ++k){
                cell[i+p].members.push_back(cell[n].members[k]);
                b[cell[n].members[k]].parent=i+p;
            }}

        for(int j=0; j<3; ++j){
            cell[i+p].r2.emplace_back(0);
            cell[i+p].rd.emplace_back(0);}

        for(int j=0; j<6; ++j)
            cell[i+p].boundary.emplace_back(0); 
        
        if(i==1){
            cell[i+p].neighbors.push_back(9+p); cell[i+p].neighbors.push_back(2+p);}
        else if(i==9){
            cell[i+p].neighbors.push_back(8+p); cell[i+p].neighbors.push_back(1+p);}
        else {
            cell[i+p].neighbors.push_back(i-1+p); cell[i+p].neighbors.push_back(i+1+p);}
        
        root(i+p);
        apply_changes(i+p);
        a1[i]=abs(cell[i+p].boundary[1]-cell[i+p].boundary[0]); r.push_back(a1[i]);
        b1[i]=abs(cell[i+p].boundary[3]-cell[i+p].boundary[2]); r.push_back(b1[i]);
        c1[i]=abs(cell[i+p].boundary[5]-cell[i+p].boundary[4]); r.push_back(c1[i]);

        if(a1[i] <= th && b1[i] <= th && c1[i] <= th){
              for(int j=0; j<cell[i+p].members.size(); ++j)
                 cell[n].child.push_back(cell[i+p].members[j]);}

        else if(a1[i] > th || b1[i] > th || c1[i] > th)
                cell[n].scell.push_back(i+p);
    }

    return r;
}

hpx::future<void> strc3(int n){

    double a1[9],b1[9],c1[9];
    std::vector<double> r; int i;
    
    std::vector<hpx::future<void>> futures;
    futures.reserve(9);

    r=AA(n);

        for (i = 1; i < 9; ++i) {
            a1[i] = r[(i - 1) * 3]; b1[i] = r[(i - 1) * 3 + 1]; c1[i] = r[(i - 1) * 3 + 2]; int p = n * 8;
            if (cell[p + i].members.size() > th1)
                if (a1[i] > th || b1[i] > th || c1[i] > th)
                    futures.push_back(hpx::async(strc3, p + i));
                }

    return hpx::when_all(futures);
}

///////////////////////////////////////////////////////////////////

void neighbor()
{
    for (int i=1; i<cell.size(); ++i){
        int temp=cell[i].ID2;
        for(int j=1; j<cell.size(); ++j)
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
                    cell[temp].neighbors[5]=cell[cell[j].ID2].ID2;}
    }
}

void traverse_tree(int root,int node, int root2)
{
    int i;
    float D=sqrt(pow((b[node].r1[0]-cell[root2].r2[0]),2.0)+pow((b[node].r1[1]-cell[root2].r2[1]),2.0)+pow((b[node].r1[2]-cell[root2].r2[2]),2.0));
    float r=sqrt(pow(cell[root2].rd[0],2.0)+pow(cell[root2].rd[1],2.0)+pow(cell[root2].rd[2],2.0));
    float ratio=float(D/r);
    if (ratio<theta){
        if(cell[root2].child.size()>=1)
            for(i=0; i<cell[root2].child.size(); ++i)
                cell[root].list_cell1.push_back(cell[root2].child[i]);
        if(cell[root2].scell.size()>=1)
            for(i=0; i<cell[root2].scell.size(); ++i)
                traverse_tree(root,node,cell[root2].scell[i]);}
    else
        cell[root].list_cell2.push_back(root2);
}

void cell_list(int N)
{
    for (int i = 0; i < N; ++i)
        if (cell[b[i].parent].list_cell1.size() == 0 && cell[b[i].parent].list_cell2.size() == 0)
            traverse_tree(b[i].parent, i, 0);
}

void create_Octree(){
    int N=n;

    initialize_body(N);
    read_Input();
    create_Hilbert(N); 
    initialize_cell(N);
    apply_changes(0);
    root(0);
    

    hpx::future<void> r;
    r=strc3(0);
    r.wait();

    neighbor();
    cell_list(N);
}

double* compute_r(int node){

    std::vector<float> a(3,0.0);
    int i,j,rs=-1;
    int parent=b[node].parent;
       
    if(cell[parent].list_cell1.size()>=1)
        for(j=0; j<3; ++j)
            for (i=0; i<cell[parent].list_cell1.size(); ++i)
                b[node].force[j]=b[node].force[j]+(G * b[node].m1 * b[cell[parent].list_cell1[i]].m1) * (b[node].r1[j]-b[cell[parent].list_cell1[i]].r1[j])/pow((1+pow((b[node].r1[0]-b[cell[parent].list_cell1[i]].r1[0]),2.0)+pow((b[node].r1[1]-b[cell[parent].list_cell1[i]].r1[1]),2) +pow((b[node].r1[2]-b[cell[parent].list_cell1[i]].r1[2]),2.0)),1.5);
   
    if(cell[parent].list_cell2.size()>=1)
        for(j=0; j<3; ++j)
            for (i=0; i<cell[parent].list_cell2.size(); ++i)
                b[node].force[j]=b[node].force[j]+(G * b[node].m1 * cell[cell[parent].list_cell2[i]].m2) * (b[node].r1[j]-cell[cell[parent].list_cell2[i]].r2[j])/pow((1+pow((b[node].r1[0]-cell[cell[parent].list_cell2[i]].r2[0]),2.0)+pow((b[node].r1[1]-cell[cell[parent].list_cell2[i]].r2[1]),2.0)+pow((b[node].r1[2]-cell[cell[parent].list_cell2[i]].r2[2]),2.0)),1.5);
  
    for(j=0; j<3; ++j){
        a[j]=b[node].force[j]/b[node].m1;
        b[node].r1[j]=b[node].r1[j]+b[node].v1[j] * t + 0.5*a[j]* t * t;}

  return b[node].force;

}

void insert_node(int parent, int node)
{
    cell[parent].NumNodes=cell[parent].NumNodes+1;
    cell[parent].members.push_back(node);
    b[node].parent=parent;

    if(cell[parent].members.size()<=th)
        cell[parent].child.push_back(node);

    if(cell[parent].members.size()>th)
    {
        int p=8*parent;
        for(int i=1; i<9; ++i)
            if(cell[i+p].members.size()>1)
                if(InCube(node,i+p))
                    insert_node(i + p, node);
    }
}

void new_tree(int node, int parent)
{
    int k=-1;
    for(int j=0; j<2; ++j)
        if(cell[parent].neighbors[j]!=-1 && k==-1)
        {
            if (InCube(node, cell[parent].neighbors[j]))
            {
                cell[parent].NumNodes = cell[parent].NumNodes + 1;
                k = 1;
                int temp = cell[parent].neighbors[j];
                insert_node(temp, node);
            }
        }

    if(k==-1 && cell[parent].ID2!=0)
        new_tree(node,cell[parent].parent2);
}

std::vector<int> get_cell_foreach_lc(std::size_t ID,std::size_t np,std::size_t nl){

    std::vector<int> cell_IDs;
    int local_np=np/nl;

    if(np==8){
        if(nl!=1)
            for(int i=0; i<local_np; ++i)
                cell_IDs.push_back(1+i+int((8/nl)*ID));
               
        else
            for(int i=0; i<local_np; ++i) 
                cell_IDs.push_back(ID+1+i);
                
    }

    if(np==64)
        for(int i=0; i<local_np; ++i)
            cell_IDs.push_back(9+i+int((64/nl)*ID));

    return cell_IDs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

partition stepper_server::compute_position(partition const &ce0,
                                           partition const &Left, partition const &Right)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;

    hpx::shared_future<partition_data> current_bodies = ce0.get_data(partition_server::cell0);
    hpx::future<partition_data> next_bodies = current_bodies.then(
            unwrapped(
                    [ce0](partition_data const&m) -> partition_data
                    {
                        std::size_t size = m.size();
                        partition_data next(size);
                      for (std::size_t i = 0; i < size; ++i) {
                            next[i].ID1=m[i].ID1;
                            if(next[i].ID1>=0){
                                double* p= compute_r(next[i].ID1);
                                next[i].force[0]=p[0]; next[i].force[1]=p[1]; next[i].force[2]=p[2];
                                b[next[i].ID1].force[0]=next[i].force[0];
                                b[next[i].ID1].force[1]=next[i].force[1];    
                                b[next[i].ID1].force[2]=next[i].force[2];
                            }
                       }
                        return next;
                    }
            )
    );

    return dataflow(
            hpx::launch::async,
            unwrapped(
                    [ce0,Left,Right](partition_data next, partition_data const& m,
                                                        partition_data const& L, partition_data const& R) -> partition
                    {
                        std::size_t size0=m.size();
                        int k=0; std::vector<std::size_t> size(2,0);
                        size[0] = L.size();
                        size[1] = R.size(); 

                        if(size[0]>0)
                            for(std::size_t i=0; i<size[0]; ++i)
                                if(k<k_th) {
                                    next[size0-k_th+k]=L[i]; k=k+1; }

                        if(size[1]>0)
                            for(std::size_t i=0; i<size[1]; ++i)
                                if(k<k_th) {
                                    next[size0-k_th+k]=R[i]; k=k+1; }

                       for(std::size_t i=size0-k_th; i<size0; ++i)
                            if(next[i].ID1>0)
                                new_tree(next[i].ID1,next[i].parent);

                        return partition(ce0.get_gid(), next);
                    }
            ),
            std::move(next_bodies),
            current_bodies,
            Left.get_data(partition_server::cell0),
            Right.get_data(partition_server::cell0)
    );
}
/////////////////////////////////////////////////////////////////////////////////
partition stepper_server::changed_members( partition const& ce0, partition const& From_N, int cell_ID, int dim, int dir)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;

    hpx::shared_future<partition_data> current_bodies = ce0.get_data(partition_server::cell0);

    hpx::future<partition_data> next_changes= current_bodies.then(
            unwrapped(
                    [ce0, cell_ID, dim, dir](partition_data const& m) -> partition_data
                    {
                        std::size_t size = m.size();
                        partition_data next(size);
                        for (std::size_t i = 0; i < size; ++i){
                            next[i]=m[i];
                            next[i].ID1=-1;
                            if(m[i].ID1>=0){
                                if (dir == -1) {
                                    if (m[i].r1[dim] < cell[cell_ID].boundary[2 * dim])
                                        next[i].ID1 = m[i].ID1;
                                 }

                                if(dir==1){
                                    if (m[i].r1[dim] > cell[cell_ID].boundary[2 * dim + 1])
                                        next[i].ID1 = m[i].ID1;
                                 }
                            }
                       }
                       return next;
                    }
            )
    );

    return dataflow(
            hpx::launch::async,
            unwrapped(
                    [ce0, cell_ID, dim, dir](partition_data next, partition_data const& m, partition_data const& l) -> partition
                    {
                        std::size_t size = m.size();
                        std::size_t size2 = l.size();
                        int k=0;
                        for(std::size_t i=0; i<size2; ++i)
                            if(l[i].ID1>0 && k<k_th){
                                if (dir == -1) {
                                    if (m[i].r1[dim] < cell[cell_ID].boundary[2 * dim]){
                                        next[size-k_th+k]=m[i]; k=k+1;}}

                                if(dir==1){
                                    if (m[i].r1[dim] > cell[cell_ID].boundary[2 * dim + 1]){
                                        next[size-k_th+k]=m[i]; k=k+1;}}
                            }
                        return partition(ce0.get_gid(), next);
                    }
            ),
            std::move(next_changes),
            current_bodies,
            From_N.get_data(partition_server::cell0)
    );
}

/////////////////////////////////////////////////////////////////////////////////


partition stepper_server::Non_Members(partition const& To_N_total, partition const & To_N)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;

    hpx::shared_future<partition_data> current_bodies = To_N_total.get_data(partition_server::cell0);

    hpx::future<partition_data> next_changes= current_bodies.then(
            unwrapped(
                    [To_N_total](partition_data const& m) -> partition_data
                    {
                        std::size_t size = m.size();
                        partition_data next(size);
                        for (std::size_t i = 0; i < size; ++i){
                            next[i]=m[i];
                        }
                        return next;
                    }
            )
    );

    return dataflow(
            hpx::launch::async,
            unwrapped(
                    [To_N_total, To_N](partition_data next,partition_data const& m,partition_data const& l) -> partition
                    {
                        std::size_t size=m.size();
                        std::size_t size2 = l.size();
                        int k=0;
                        for (std::size_t i = 0; i < size; ++i) {
                            if (next[i].ID1 < 0 && k==0) {
                                k=1;
                                if (i + size2 < size)
                                    for (std::size_t j = 0; j < size2; ++j)
                                        next[i + j] = l[j];
                                else if (i + size2 > size)
                                    for (std::size_t j = i; j < size; ++j)
                                        next[j] = l[j-i];
                            }
                        }
                      return partition(To_N_total.get_gid(), next);
                    }
            ),
            std::move(next_changes),
            current_bodies,
            To_N.get_data(partition_server::cell0)
    );
}

/////////////////////////////////////////////////////////////////////////////////

partition stepper_server::New_Members(partition const& ce0, partition const& From_N, int ID)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;

    hpx::shared_future<partition_data> current_bodies = From_N.get_data(partition_server::cell0);

    hpx::future<partition_data> next_changes= current_bodies.then(
            unwrapped(
                    [From_N, ID](partition_data const& m) -> partition_data
                    {
                        std::size_t size = m.size();
                        partition_data next(size);
                        for (std::size_t i = 0; i < size; ++i){
                            next[i]=m[i];
                        }
                        return next;
                    }
            )
    );

    return dataflow(
            hpx::launch::async,
            unwrapped(
                    [ce0, From_N, ID](partition_data next,partition_data const& l) -> partition
                    {
                        std::size_t size2 = l.size();
                        int k=0;
                        for (std::size_t i = 0; i < size2; ++i) {
                            if (l[i].ID1 >= 0 && k<k_th) {
                                if (InCube(l[i].ID1, ID)) {
                                    next[i] = l[i]; k=k+1; }
                            }
                        }
                        return partition(ce0.get_gid(), next);
                    }
            ),
            std::move(next_changes),
            From_N.get_data(partition_server::cell0)
    );
}

partition stepper_server::copy_mem(partition const & ce0){

        using hpx::lcos::local::dataflow; 
        using hpx::util::unwrapped; 

        hpx::shared_future<partition_data> current_bodies = ce0.get_data(partition_server::cell0);
        hpx::future<partition_data> next_changes= current_bodies.then( 
                unwrapped( 
                    [ce0](partition_data const& m) -> partition_data 
                    {
                        std::size_t size = m.size(); 
                        partition_data next(size);
                        for (std::size_t i = 0; i < size; ++i)
                            next[i]=m[i]; 
                        return next; 
                    }
            )
       );

      return dataflow(hpx::launch::async,unwrapped([ce0](partition_data next)->partition
                {
                    return partition(ce0.get_gid(), next); 
                } 
            ),std::move(next_changes));
}
///////////////////////////////////////////////////////////////////////////////

stepper_server::space stepper_server::do_work(std::size_t np, std::size_t nl)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;
   
    std::size_t local_np=np/nl;
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
   
    std::cout<<std::endl;    

    for (space& s: U_)
        s.resize(local_np);
    
    for (space& I: Parent_)
        I.resize(1);

    std::size_t ID_lc=hpx::get_locality_id();
    hpx::id_type here = hpx::find_here();
    std::vector<int> cell_IDs=get_cell_foreach_lc(ID_lc,np,nl);
    int s1=int(n/nl);

    for (std::size_t j = 0; j != local_np; ++j)
        U_[0][j]=partition(here,cell[cell_IDs[j]].members.size(),cell_IDs[j],1);    
 
    Parent_[0][0]=partition(here,s1,1,0);

    send_L(0, U_[0][0]);
    send_R(0, U_[0][0]);

  
    for (std::size_t time = 0; time !=nt ; ++time)
    {
        space const& current_members = U_[time % 2];
        space & next_members = U_[(time + 1) % 2]; 
        space & n_ = U_[time % 2];

        space & To_L_ = U_[(time + 1) % 2]; space & To_R_ = U_[(time + 1) % 2];
        
        space & L= Parent_[time % 2]; space & R= Parent_[time % 2];

        space & From_L= Parent_[(time + 1) % 2]; space & From_R= Parent_[(time + 1) % 2];

        From_L[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, L[0], receive_L(time), ID_lc);
        From_R[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, R[0], receive_R(time), ID_lc);
    
      for(std::size_t j=0; j<local_np; ++j) {
          int i = cell_IDs[j];

          next_members[j] = dataflow(hpx::launch::async, &stepper_server::compute_position, current_members[j],From_L[0], From_R[0]);

          n_[j]=dataflow(hpx::launch::async, &stepper_server::copy_mem, next_members[j]);
           

          To_R_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, n_[j], From_L[0], i, 1, 1);
          To_L_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, n_[j], From_R[0], i, 1, -1);

          From_L[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, L[0], To_L_[j]);
          From_R[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, R[0], To_R_[j]);

     }

     if (time != nt-1) {
        send_L(time + 1, From_L[0]);
        send_R(time + 1, From_R[0]);}
  
     }

     return U_[nt % 2];
}

HPX_REGISTER_GATHER(stepper_server::space, stepper_server_space_gatherer);

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm)
{
    using hpx::lcos::local::dataflow;
        
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    std::size_t nl = localities.size();                    // Number of localities

    int np; 

    if(nl<9) np=8; //in level 1

    else if(nl>15 && nl<=64) np=64; //in level 2

    create_Octree();
    
    stepper step;
    
    boost::uint64_t t1=hpx::util::high_resolution_clock::now();

    hpx::future<stepper_server::space> result = step.do_work(np, nl);
 
    if (0 == hpx::get_locality_id())
    {
    boost::uint64_t const num_worker_threads = hpx::get_num_worker_threads();
        
        hpx::future<std::vector<stepper_server::space> > overall_result =
                hpx::lcos::gather_here(gather_basename, std::move(result), nl);

        std::vector<stepper_server::space> solution = overall_result.get();
        for (std::size_t j = 0; j != nl; ++j)
        {
            stepper_server::space const& s = solution[j];    
            for (std::size_t i = 0; i != s.size(); ++i)             
                s[i].get_data(partition_server::cell0).wait(); 

           // boost::uint64_t e1=hpx::util::high_resolution_clock::now()-t1;
            //std::cout<<(boost::format("%.14g")%(e1/1e9))<<std::flush;
            //std::cout<<std::endl;
    }

    boost::uint64_t e=hpx::util::high_resolution_clock::now()-t1;
    std::cout<<(boost::format("%.14g")%(e/1e9))<<std::flush;
    }
    else
    {
        hpx::lcos::gather_there(gather_basename, std::move(result)).wait();
    }



 
  if(0==hpx::get_locality_id())
        for(int i=0; i<200; ++i)
            std::cout<<b[i].force[0]<<"," ; 


    std::cout<<"END"<<std::endl;
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    std::vector<std::string> cfg;
    cfg.push_back("hpx.run_hpx_main!=1");

    return hpx::init(argc, argv, cfg);
}
