#include <sstream>
#include <vector>  

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/runtime/serialization/serialize.hpp>

#include <boost/timer.hpp>
#include <boost/shared_array.hpp>
#include <stack>

///////////////////////////////////////////////////////////////////////////////
int th=100;
int nt=10;
int t=1;     
double G=6.673*pow(10.0,-11.0);
double theta=0.5;
std::size_t SIZE=1400;

char const* stepper_basename = "/Dist3/stepper/";
char const* gather_basename = "/Dist3/gather/";

//--------------------------------------------------
struct Body
{
public:
    int ID1,parent,m1;
    double r1[3], v1[3], force[3];
    
    friend class hpx::serialization::access;
    template<typename Ar> void serialize(Ar &ar, unsigned){
     ar &ID1 &parent &m1 &r1 &v1 &force;}
}body1;

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

std::vector<Body> b;
std::vector<Cell> cell(10000);
 
template<typename Ar> void serialize(Ar &ar, unsigned){ 
        ar & b & cell;}
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

    partition_data(std::size_t size,std::size_t size_mem, int cell_IDs, int flag)
            : data_(alloc_.allocate(size), size, buffer_type::take,
                    &partition_data::deallocate),
              size_(size),
              min_index_(0)
    {
        if(flag==1){
            for(std::size_t j=0; j<size_mem; ++j)
                data_[j]=b[cell[cell_IDs].members[j]];
                 
            for(std::size_t j=size_mem; j<size; ++j){
                data_[j]=b[0];
                data_[j].ID1=-1;  
            }
        }
        else
          for(std::size_t j=0; j<size; ++j){
                data_[j]=b[0];
                data_[j].ID1=-1;
        }  
  }

    partition_data(partition_data const& base, std::size_t min_index)
            : data_(base.data_.data()+min_index, 1, buffer_type::reference,
                    hold_reference(base.data_)),      // keep referenced partition alive
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

    if(dir==+1 || dir==-1)
    {
        if(i == 0 && dir == -1) return size-1;
        if(i == size-1 && dir == +1) return 0;

        HPX_ASSERT((i + dir) < size);
        return i + dir;
    }

    if(dir==+2 || dir==-2)
    {
        if(size==8 || size==16 || size==32){
            if(i==0 && dir==-2) return size-1; if(i==1 && dir==-2) return size-2;
            if(i==size-1 && dir==+2) return 0; if(i==size-2 && dir==+2) return 1;
        }

        if(size==64){
            if(i==0 && dir==-2) return size-1; if(i==1 && dir==-2) return size-2;
            if(i==2 && dir==-2) return size-3; if(i==3 && dir==-2) return size-4;
            if(i==size-1 && dir==+2) return 0; if(i==size-2 && dir==+2) return 1;
            if(i==size-3 && dir==+2) return 2; if(i==size-4 && dir==+2) return 3;            
        }
        int temp;
        if(size==8 || size==16 || size==32) temp=i+dir;
        if(size==64) temp=i+(2*dir);
        if(size<8) temp=i;
        HPX_ASSERT(temp < size);
        return temp;    
    }   
    
    if(dir==+4 || dir==-4){
       
        if(size==4){
            if(i==0 && dir==-4) return size-1; if(i==1 && dir==-4) return size-2;
            if(i==size-1 && dir==+4) return 0; if(i==size-2 && dir==+4) return 1;
        } 

        if(size==8){
            if(i==0 && dir==-4) return size-1; if(i==1 && dir==-4) return size-2;
            if(i==2 && dir==-4) return size-3; if(i==3 && dir==-4) return size-4;
            if(i==size-1 && dir==+4) return 0; if(i==size-2 && dir==+4) return 1;
            if(i==size-3 && dir==+4) return 2; if(i==size-4 && dir==+4) return 3;
        }

        if(size==16 || size==32){
            if(i==0 && dir==-4) return size-1; if(i==1 && dir==-4) return size-2;
            if(i==2 && dir==-4) return size-3; if(i==3 && dir==-4) return size-4; 
            if(i==4 && dir==-4) return size-5; if(i==5 && dir==-4) return size-6; 
            if(i==6 && dir==-4) return size-7; if(i==7 && dir==-4) return size-8; 
            if(i==size-1 && dir==+4) return 0; if(i==size-2 && dir==+4) return 1;
            if(i==size-3 && dir==+4) return 2; if(i==size-4 && dir==+4) return 3; 
            if(i==size-5 && dir==+4) return 4; if(i==size-6 && dir==+4) return 5; 
            if(i==size-7 && dir==+4) return 6;
        }

        int temp;
        if(size==4) temp=i+int(dir*0.5);
        if(size==8) temp=i+dir;
        if(size==16 || size==32) temp=i+(dir*2);
        if(size==64) temp=i+(dir*4);
        if(size<4) temp=i;
           
        HPX_ASSERT(temp < size);  
        return temp;
     }
        
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

    partition_server(std::size_t size,std::size_t size_mem, int cell_IDs, int flag)
            : data_(size,size_mem, cell_IDs,flag)
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

    partition(hpx::id_type where, std::size_t size,std::size_t size_mem, int cell_IDs, int flag)
            : base_type(hpx::new_<partition_server>(where,size,size_mem,cell_IDs,flag))
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
              Back_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(),-2, nl))),
              Front_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(),+2, nl))),
              Up_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(),-4, nl))),
              Down_(hpx::find_id_from_basename(stepper_basename, idx(hpx::get_locality_id(),+4, nl))),
              Parent_(6),
              U_(6),
              To_B_(1),
              To_F_(1),
              To_L_(1),
              To_R_(1),
              To_U_(1),
              To_D_(1)
    {
    }

    space do_work(std::size_t np, std::size_t nl);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

    void from_R(std::size_t time, partition p) {
        Right_receive_buffer_.store_received(time, std::move(p));}

    void from_L(std::size_t time, partition p){
        Left_receive_buffer_.store_received(time, std::move(p));}

    void from_B(std::size_t time, partition p) {
        send_B_buffer_.store_received(time, std::move(p)); }

    void from_F(std::size_t time, partition p) {
        send_F_buffer_.store_received(time, std::move(p)); }

    void from_U(std::size_t time, partition p) {
        send_U_buffer_.store_received(time, std::move(p)); }

    void from_D(std::size_t time, partition p) {
        send_D_buffer_.store_received(time, std::move(p)); }

    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_R, from_R_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_L, from_L_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_B, from_B_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_F, from_F_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_U, from_U_action);
    HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_D, from_D_action);

protected:
    static void create_Octree();
    
    static partition compute_position(partition const &ce0, partition const &Back, partition const &Front,
                                      partition const &Left, partition const &Right,
                                      partition const &Up, partition const &Down);

    static partition changed_members(partition const &current_members, partition const & From_N,
                                     int cell_ID, int dim, int dir);

    static partition New_Members(partition const& ce0, partition const & From_N, int ID);
    static partition Non_Members(partition const& To_N_total, partition const & To_N);
    
    void initialize_body(int N);
    void read_Input();
    void initialize_cell(int N);
    void apply_changes(int n);
    static bool InCube(int node, int p);
    static void det_boundary_subcube(int n);
    void root(int n);
    void A(int n);
    void strc3(int n);
    void neighbor();
    void traverse_tree(int root, int node, int root2);
    void cell_list(int N);
    static void compute_r(int node);
    static void insert_node(int parent, int node);
    static void new_tree(int node, int parent);
    
    partition receive_L(std::size_t time) { return Left_receive_buffer_.receive(time);}
    partition receive_R(std::size_t time) { return Right_receive_buffer_.receive(time);}

    partition receive_U(std::size_t time) { return send_U_buffer_.receive(time);}
    partition receive_D(std::size_t time) { return send_D_buffer_.receive(time);}
    partition receive_B(std::size_t time) { return send_B_buffer_.receive(time);}
    partition receive_F(std::size_t time) { return send_F_buffer_.receive(time);} 

    void send_L(std::size_t time, partition p){ hpx::apply(from_R_action(), Left_.get(), time, p);}
    void send_R(std::size_t time, partition p){ hpx::apply(from_L_action(), Right_.get(), time, p);}
    void send_U(std::size_t time, partition p) { hpx::apply(from_U_action(), Up_.get(), time, p);}
    void send_D(std::size_t time, partition p) { hpx::apply(from_D_action(), Down_.get(), time, p);}
    void send_B(std::size_t time, partition p) { hpx::apply(from_B_action(), Back_.get(), time, p);}
    void send_F(std::size_t time, partition p) { hpx::apply(from_F_action(), Front_.get(), time, p);}

private:
    hpx::shared_future<hpx::id_type> Left_, Right_, Back_, Front_, Up_, Down_; 
    std::vector<space> U_, Parent_, To_B_, To_F_, To_L_, To_R_, To_U_, To_D_;;
    hpx::lcos::local::receive_buffer<partition> Left_receive_buffer_;
    hpx::lcos::local::receive_buffer<partition> Right_receive_buffer_;      
    hpx::lcos::local::receive_buffer<partition> send_B_buffer_;
    hpx::lcos::local::receive_buffer<partition> send_F_buffer_;
    hpx::lcos::local::receive_buffer<partition> send_U_buffer_; 
    hpx::lcos::local::receive_buffer<partition> send_D_buffer_;
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

void stepper_server::initialize_body(int N)
{    
    body1.ID1=0; body1.m1=0; body1.parent=0;
    for(int i=0; i<3; ++i){
        body1.r1[i]=1.0; body1.v1[i]=1.0; body1.force[i]=0.0;
    }
    for(int i=0; i<N; ++i)
        b.push_back(body1);
}

void stepper_server::read_Input()
{
    int N; std::string temp,temp1;
    std::fstream textfile;
    std::string::size_type sz;
    textfile.open("/home/zkhatami/Dist2/ex10000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        b[i].ID1=atoi(temp.c_str()); b[i].m1=stof(temp1,&sz);
         for (int j=2; j<5; ++j){
            textfile>>temp; b[i].r1[j-2]=stof(temp,&sz);}
         for (int j=5; j<8; ++j){
            textfile>>temp; b[i].v1[j-5]=stof(temp,&sz);}  
    b[i].parent=0;
    }
}

void stepper_server::initialize_cell(int N)
{
    cell[0].ID2=0; cell[0].parent2=0; cell[0].NumNodes=N; cell[0].level=0;
    cell[0].neighbors.resize(6); cell[0].r2.resize(3); cell[0].rd.resize(3); cell[0].boundary.resize(6);
    
    cell[0].members.resize(N);
    for(int i=0; i<3; ++i){ cell[0].r2[i]=0; cell[0].rd[i]=0; cell[0].boundary[2*i]=0; cell[0].boundary[2*i+1]=100;}
    for(int i=0; i<6; ++i){cell[0].neighbors[i]=-1;}    
    for (int i = 0; i < N; ++i)
            cell[0].members[i]=i;  
}

void stepper_server::apply_changes(int n)
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

bool stepper_server::InCube(int node, int p)
{
    for (int i=0; i<3; ++i)
       if(b[node].r1[i]<cell[p].boundary[2*i] || b[node].r1[i]>cell[p].boundary[2*i+1])
    
       return false;
    return true;
}

void stepper_server::det_boundary_subcube(int n)
{
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    int p=n*8;

    a1=cell[n].boundary[0]; a2=(cell[n].boundary[0]+cell[n].boundary[1])/2; a3=cell[n].boundary[1];
    b1=cell[n].boundary[2]; b2=(cell[n].boundary[2]+cell[n].boundary[3])/2; b3=cell[n].boundary[3];
    c1=cell[n].boundary[4]; c2=(cell[n].boundary[4]+cell[n].boundary[5])/2; c3=cell[n].boundary[5];

    std::vector<double> A1, A2, A3, A4, A5, A6, A7, A8; std::vector<int> B2(6,-1); std::vector<double> B1(3,0);
    A1.push_back(a1); A1.push_back(a2); A1.push_back(b1); A1.push_back(b2); A1.push_back(c1); A1.push_back(c2);
    A2.push_back(a1); A2.push_back(a2); A2.push_back(b1); A2.push_back(b2); A2.push_back(c2); A2.push_back(c3);
    A3.push_back(a1); A3.push_back(a2); A3.push_back(b2); A3.push_back(b3); A3.push_back(c1); A3.push_back(c2);
    A4.push_back(a1); A4.push_back(a2); A4.push_back(b2); A4.push_back(b3); A4.push_back(c2); A4.push_back(c3);
    A5.push_back(a2); A5.push_back(a3); A5.push_back(b1); A5.push_back(b2); A5.push_back(c1); A5.push_back(c2);
    A6.push_back(a2); A6.push_back(a3); A6.push_back(b1); A6.push_back(b2); A6.push_back(c2); A6.push_back(c3);
    A7.push_back(a2); A7.push_back(a3); A7.push_back(b2); A7.push_back(b3); A7.push_back(c1); A7.push_back(c2);
    A8.push_back(a2); A8.push_back(a3); A8.push_back(b2); A8.push_back(b3); A8.push_back(c2); A8.push_back(c3);

    for(int j=1; j<9; ++j){ cell[p+j].r2=B1; cell[p+j].rd=B1; cell[p+j].neighbors=B2;}

    cell[p+1].boundary=A2; cell[p+2].boundary=A4; cell[p+3].boundary=A6; cell[p+4].boundary=A8;
    cell[p+5].boundary=A1; cell[p+6].boundary=A3; cell[p+7].boundary=A5; cell[p+8].boundary=A7;
}

void stepper_server::root(int n)
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

void stepper_server::A(int n)
{
    int p = n * 8;
    for (int i = 1; i < 9; ++i)
    {
        cell[p + i].ID2 = p + i;
        cell[p + i].NumNodes=0;
        cell[p + i].level = cell[n].level + 1;
        cell[p + i].parent2 = n;

        for (int j = 0; j < cell[n].NumNodes; ++j) 
            if (InCube(cell[n].members[j], p + i)) {
                cell[p + i].members.push_back(cell[n].members[j]);
                cell[p + i].NumNodes = cell[p + i].NumNodes + 1;
            }

       if (cell[p + i].members.size() > 1) root(p + i);
            if (cell[p + i].members.size() >= 1 && cell[p + i].members.size() <= th) {
                for (int j = 0; j < cell[p + i].members.size(); ++j) {
                    b[cell[p + i].members[j]].parent = n;
                    cell[n].child.push_back(cell[p + i].members[j]);
            }}

        if (cell[p + i].members.size() > th) {
            cell[n].scell.push_back(p + i);
            det_boundary_subcube(p + i);}
    }
}

void stepper_server::strc3(int n)
{
    A(n);
    int p=n*8;
    for(int i=1; i<9; ++i){
        if(cell[p + i].members.size()>th)
            strc3(p+i);}
}

void stepper_server::neighbor()
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

void stepper_server::traverse_tree(int root,int node, int root2)
{
    int i;
    float D=sqrt(pow((b[node].r1[0]-cell[root2].r2[0]),2.0)+pow((b[node].r1[1]-cell[root2].r2[1]),2.0)+pow((b[node].r1[2]-cell[root2].r2[2]),2.0));
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
        cell[root].list_cell2.push_back(root2);
}

void stepper_server::cell_list(int N)
{
    for (int i = 0; i < N; ++i)
        if (cell[b[i].parent].list_cell1.size() == 0 && cell[b[i].parent].list_cell2.size() == 0)
            traverse_tree(b[i].parent, i, 0);
}

void stepper_server::compute_r(int node){

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
}

void stepper_server::insert_node(int parent, int node)
{
    cell[parent].NumNodes=cell[parent].NumNodes+1;
    cell[parent].members.push_back(node);
    b[node].parent=parent;

    if(cell[parent].members.size()<=th)
        cell[parent].child.push_back(node);

    if(cell[parent].members.size()>th)
    {
        int p=8*parent;
        det_boundary_subcube(parent);
        for(int i=1; i<9; ++i)
            if(InCube(node,i+p))
                insert_node(i + p, node);
    }
}

void stepper_server::new_tree(int node, int parent)
{
    int k=-1;
    for(int j=0; j<6; ++j)
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
                cell_IDs.push_back(ID+1+i);}

    if(np==64)
        for(int i=0; i<local_np; ++i)
            cell_IDs.push_back(9+i+int((64/nl)*ID));

    return cell_IDs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stepper_server::create_Octree()
{
    stepper_server foo;
    std::string temp,temp1;
    std::fstream textfile;
    textfile.open("/home/zkhatami/Dist2/ex10000.txt");
    textfile >> temp;
    int N = (int) atoi(temp.c_str());

    foo.initialize_body(N);
    foo.read_Input();
    foo.initialize_cell(N);
    foo.apply_changes(0);
    foo.root(0);
    foo.det_boundary_subcube(0);
    
    foo.strc3(0);
    foo.neighbor();
    foo.cell_list(cell.size());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

partition stepper_server::compute_position(partition const &ce0,
                                           partition const &Back, partition const &Front,
                                           partition const &Left, partition const &Right,
                                           partition const &Up, partition const &Down)
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
                            next[i]=m[i];
                            if(m[i].ID1>0)
                                compute_r(m[i].ID1);
                        }
                        return next;
                    }
            )
    );

    return dataflow(
            hpx::launch::async,
            unwrapped(
                    [ce0,Back,Front,Left,Right,Up,Down](partition_data next, partition_data const& m,
                                                        partition_data const& B, partition_data const& F,
                                                        partition_data const& L, partition_data const& R,
                                                        partition_data const& U, partition_data const& D) -> partition
                    {
                        std::size_t size0=m.size();
                        int k=0; std::vector<std::size_t> size(6,0);
                        size[0] = B.size(); size[1] = F.size(); size[2] = L.size();
                        size[3] = R.size(); size[4] = U.size(); size[5] = D.size();

                        if(size[0]>0)
                            for(std::size_t i=0; i<size[0]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=B[i]; k=k+1;
                                    if(B[i].ID1>0) compute_r(B[i].ID1); }

                        if(size[1]>0)
                            for(std::size_t i=0; i<size[1]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=F[i]; k=k+1;
                                    if(F[i].ID1>0) compute_r(F[i].ID1); }

                        if(size[2]>0)
                            for(std::size_t i=0; i<size[2]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=L[i]; k=k+1;
                                    if(L[i].ID1>0) compute_r(L[i].ID1); }

                        if(size[3]>0)
                            for(std::size_t i=0; i<size[3]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=R[i]; k=k+1;
                                    if(R[i].ID1>0) compute_r(R[i].ID1); }

                        if(size[4]>0)
                            for(std::size_t i=0; i<size[4]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=U[i]; k=k+1;
                                    if(U[i].ID1>0) compute_r(U[i].ID1); }

                        if(size[5]>0)
                            for(std::size_t i=0; i<size[5]; ++i)
                                if(k<200) {
                                    next[size0-200+k]=D[i]; k=k+1;
                                    if(D[i].ID1>0) compute_r(D[i].ID1); }

                       for(std::size_t i=size0-200; i<size0; ++i)
                            if(next[i].ID1>0)
                                new_tree(next[i].ID1,next[i].parent);

                        return partition(ce0.get_gid(), next);
                    }
            ),
            std::move(next_bodies),
            current_bodies,
            Back.get_data(partition_server::cell0),
            Front.get_data(partition_server::cell0),
            Left.get_data(partition_server::cell0),
            Right.get_data(partition_server::cell0),
            Up.get_data(partition_server::cell0),
            Down.get_data(partition_server::cell0)
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
                                    if (b[m[i].ID1].r1[dim] < cell[cell_ID].boundary[2 * dim])
                                        next[i].ID1 = m[i].ID1; }

                                if(dir==1){
                                    if (b[m[i].ID1].r1[dim] > cell[cell_ID].boundary[2 * dim + 1])
                                        next[i].ID1 = m[i].ID1; }
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
                            if(l[i].ID1>0 && k<200){
                                if (dir == -1) {
                                    if (b[m[i].ID1].r1[dim] < cell[cell_ID].boundary[2 * dim]){
                                        next[size-200+k]=m[i]; k=k+1;}}

                                if(dir==1){
                                    if (b[m[i].ID1].r1[dim] > cell[cell_ID].boundary[2 * dim + 1]){
                                        next[size-200+k]=m[i]; k=k+1;}}
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
                            next[i].ID1=-1;
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
                            if (l[i].ID1 >= 0 && k<200) {
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

///////////////////////////////////////////////////////////////////////////////

stepper_server::space stepper_server::do_work(std::size_t np, std::size_t nl)
{
    using hpx::lcos::local::dataflow;
    using hpx::util::unwrapped;

    std::size_t local_np=np/nl;
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    create_Octree();  
    for (space& s: U_)
        s.resize(local_np);
    
    for (space& I: Parent_)
        I.resize(1);
    
    std::size_t ID_lc=hpx::get_locality_id();
    hpx::id_type here = hpx::find_here();
    std::vector<int> cell_IDs=get_cell_foreach_lc(ID_lc,np,nl);
    
    for(int i=0; i<cell_IDs.size(); ++i)
        std::cout<<cell_IDs[i]<<",";
    std::cout<<std::endl;
 

    for (std::size_t j = 0; j != local_np; ++j)
        U_[0][j]=partition(here,SIZE,cell[cell_IDs[j]].members.size(),cell_IDs[j],1);    


    Parent_[0][0]=partition(hpx::find_here(),SIZE,0,0,0);

    send_L(0, U_[0][0]);
    send_R(0, U_[0][0]);
    send_B(0, U_[0][0]); 
    send_F(0, U_[0][0]);
    send_U(0, U_[0][0]); 
    send_D(0, U_[0][0]);

    for (std::size_t time = 0; time !=nt ; ++time)
    {
        space const& current_members = U_[time % 2];
        space & next_members = U_[(time + 1) % 2]; //This is for next bodies in each locality

        space & To_B_ = U_[(time + 1) % 2]; space & To_F_ = U_[(time + 1) % 2];
        space & To_L_ = U_[(time + 1) % 2]; space & To_R_ = U_[(time + 1) % 2];
        space & To_U_ = U_[(time + 1) % 2]; space & To_D_ = U_[(time + 1) % 2];

        space & B= Parent_[time % 2]; space & F= Parent_[time % 2];
        space & L= Parent_[time % 2]; space & R= Parent_[time % 2];
        space & U= Parent_[time % 2]; space & D= Parent_[time % 2];

        space & From_B= Parent_[(time + 1) % 2]; space & From_F= Parent_[(time + 1) % 2];
        space & From_L= Parent_[(time + 1) % 2]; space & From_R= Parent_[(time + 1) % 2];
        space & From_U= Parent_[(time + 1) % 2]; space & From_D= Parent_[(time + 1) % 2];

        From_B[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, B[0], receive_B(time), ID_lc);
        From_F[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, F[0], receive_F(time), ID_lc);
        From_L[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, L[0], receive_L(time), ID_lc);
        From_R[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, R[0], receive_R(time), ID_lc);
        From_U[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, U[0], receive_U(time), ID_lc);
        From_D[0] = dataflow(hpx::launch::async, &stepper_server::New_Members, D[0], receive_D(time), ID_lc);

        for(std::size_t j=0; j<local_np; ++j) {
            int i = cell_IDs[j];

            next_members[j] = dataflow(hpx::launch::async, &stepper_server::compute_position, current_members[j],
                                       From_B[0], From_F[0], From_L[0], From_R[0], From_U[0], From_D[0]);

            To_B_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_B[0], i, 0, 1);
            To_F_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_F[0], i, 0, -1);
            To_R_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_R[0], i, 1, 1);
            To_L_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_L[0], i, 1, -1);
            To_U_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_U[0], i, 2, 1);
            To_D_[j] = dataflow(hpx::launch::async, &stepper_server::changed_members, current_members[j], From_D[0], i, 2, -1);

            From_B[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, B[0], To_B_[j]);
            From_F[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, F[0], To_F_[j]);
            From_L[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, L[0], To_L_[j]);
            From_R[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, R[0], To_R_[j]);
            From_U[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, U[0], To_U_[j]);
            From_D[0] = dataflow(hpx::launch::async, &stepper_server::Non_Members, D[0], To_D_[j]);
  
     }

     if (time != nt-1) {
        send_B(time + 1, From_B[0]);
        send_F(time + 1, From_F[0]);
        send_L(time + 1, From_L[0]);
        send_R(time + 1, From_R[0]);
        send_U(time + 1, From_U[0]);
        send_D(time + 1, From_D[0]);}
  
     }

     return U_[0 % 2];
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
    boost::timer timer;
    stepper step;
    std::cout<<nl<<"-"<<np<<std::endl;

    hpx::future<stepper_server::space> result = step.do_work(np, nl);
  
    if (0 == hpx::get_locality_id())
    {
        boost::uint64_t const num_worker_threads = hpx::get_num_worker_threads();

        hpx::future<std::vector<stepper_server::space> > overall_result =
                hpx::lcos::gather_here(gather_basename, std::move(result), nl);

        std::vector<stepper_server::space> solution = overall_result.get();
        for (std::size_t j = 0; j != nl; ++j)
        {
            stepper_server::space s = solution[j];    
            for (std::size_t i = 0; i != s.size(); ++i)             
                s[i].get_data(partition_server::cell0).wait(); 
        }

    }
    else
    {
        hpx::lcos::gather_there(gather_basename, std::move(result)).wait();
    }

    double elapsed_time = timer.elapsed();
    std::cout << std::endl << "FIRST CPU TIME: " << elapsed_time << std::endl;

    for(int i=0; i<200; ++i)
        std::cout<<b[i].r1[1]<<" , ";

    std::cout<<"END";
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    std::vector<std::string> cfg;
    cfg.push_back("hpx.run_hpx_main!=1");

    return hpx::init(argc, argv, cfg);
}
