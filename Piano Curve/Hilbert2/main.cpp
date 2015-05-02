//
//  main.cpp
//  Hilbert2
//
//  Created by Zahra Khatami on 4/24/15.
//  Copyright (c) 2015 Zahra Khatami. All rights reserved.
//


#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>
using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;

//Data for each node
struct Body {
    int ID1,parent,Hilbert_dis,num_A;
    double m1,v1[3],force[3]={};
    unsigned int r1[3];
};

//Read input from txt file
void read_Input(Body body[],  vector<int> &A){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("Input.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        body[i].ID1=atoi(temp.c_str()); A.push_back(body[i].ID1);
        body[i].m1=stof(temp1,&sz); body[i].num_A=i; body[i].parent=0;
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}}}


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

int main(){
    int N; string temp; fstream textfile; vector<int> A,H; int m;
    textfile.open("Input.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    Body *body=new Body[N];
    read_Input(body,A);
    cout<<N<<endl;
    vector<int>  X, Y, Z;
    for(int i=0; i<N; ++i){
        cout<<"\n ----"<<i<<endl;
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
            cout<<H[j];
        
        int sum=0;
        for(int j=0; j<H.size(); ++j)
            sum=sum+H[j]*pow(2,j);
        body[i].Hilbert_dis=sum;}
    
    
    cout<<"\n ==== \n";
    for(int i=0; i<N; ++i)
        cout<<body[i].Hilbert_dis<<",";
    
    return 0;}

