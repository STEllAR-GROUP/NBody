//Determining the boundry of the each subcube
void det_boundary_subcube(int n){
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    int p=n*8;

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

    for(int j=1; j<9; ++j){
        for(int i=0; i<3; ++i){ cell[p+j].r2.push_back(0); cell[p+j].rd.push_back(0);}
        for(int i=0; i<6; ++i) cell[p+j].neighbors.push_back(-1);}

    cell[p+1].boundary=A1; cell[p+2].boundary=A2; cell[p+3].boundary=A3; cell[p+4].boundary=A4;
    cell[p+5].boundary=A5; cell[p+6].boundary=A6; cell[p+7].boundary=A7; cell[p+8].boundary=A8;


}

//creating octree and determining its parameters
template<class type> class octree {
public:
    octree() {};
    void root(int n);
    void strc(int n);
    void strc2(int n);
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

/////////////////////////////////////////////////////////////////// Parellel Octree

template <class type> void octree<type>::strc2(int n){
    octree<float> tree;
    using hpx::parallel::for_each;
    using hpx::parallel::par;
    typedef boost::counting_iterator<int> iterator;
    int p=n*8;
    Cell cell1;
    for(int i=1; i<9; ++i){

        cell[p+i].ID2=p+i;
        for(int j=0; j<cell[n].members.size(); ++j){
            if(InCube(cell[n].members[j],p+i)){
                cell[p+i].members.push_back(cell[n].members[j]);
            }}
        cell[p+i].NumNodes=cell[p+i].members.size();
        if(cell[p+i].members.size()>1) tree.root(p+i);

        if(cell[p+i].members.size()>=1 && cell[p+i].members.size()<=th){
            for(int j=0; j<cell[p+i].members.size(); ++j){
                body[cell[p+i].members[j]].parent=n;
                cell[n].child.push_back(cell[p+i].members[j]);}}

        if(cell[p+i].members.size()>th){
            cell[n].scell.push_back(p+i);
            cell[p+i].parent2=n;
            det_boundary_subcube(p+i);
         }}

    using namespace hpx::parallel;
    for_each(par, iterator(1), iterator(9), [&](int i) {
        if(cell[n*8+i].members.size()>th)
            tree.strc2(cell[n*8+i].ID2);
    });
}