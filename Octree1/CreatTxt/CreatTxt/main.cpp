// basic file operations
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main () {
    srand (static_cast <unsigned> (time(0)));
    ofstream myfile; string c;
    myfile.open ("example.txt");
    cout<<"\n Please insert the number of nodes in Gallaxy! \n";
    cin>>c;
    myfile << c<<"\n";
    int N=atoi(c.c_str());
    for(int i=0; i<N; ++i){
        myfile <<i<< " "  <<1;
        for(int j=0; j<6; ++j){
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            myfile<<" " ;
            myfile<< r;
            }
         myfile<< " "<<0<<" " <<0;
        myfile <<"\n";
    }
    
    
    myfile.close();
    return 0;
}