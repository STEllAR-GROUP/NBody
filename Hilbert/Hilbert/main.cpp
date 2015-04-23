//
//  main.cpp
//  Hilbert
//
//  Created by Zahra Khatami on 4/11/15.
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
    textfile.open("Input2.txt");
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


//Hilbert distance
void Hilbert_distance(coord_t* X, int b, int n) //(this uint[] hilbertAxes, int bits)
{
    coord_t M = 1 << (b - 1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
            if ((X[i] & Q) != 0)
                X[0] ^= P; // invert
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;}
    } // exchange
    // Gray encode
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if ((X[n - 1] & Q)!=0)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;}


void merge(vector<int> &A,Body body[],int, int , int ,int);
void sort_merge(vector<int> &A,Body body[], int low,int high, int N){
    int mid;
    if(low<high){
        mid=(low+high)/2;
        sort_merge(A,body,low,mid,N);
        sort_merge(A,body,mid+1,high,N);
        merge(A,body,low,mid,high,N);
    }
}
void merge(vector<int> &A,Body body[],int low, int mid, int high, int N){
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



int main(int argc, const char * argv[]) {
    
    int N; string temp; fstream textfile; vector<int> A;
    textfile.open("Input2.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    Body *body=new Body[N];
    read_Input(body,A);
    
    
    for(int i=0; i<N; ++i){
        //cout<<"\n ------\n";
        coord_t X[3]={body[i].r1[0],body[i].r1[1],body[i].r1[2]};
        int m=10;
        int n=m*3; int Y[n];
        Hilbert_distance(X,m, 3); int sum=0;
        int j=n-1, M1=m;
        while(j>=0){
            M1=M1-1;
            for(int k=0; k<3; ++k){
                Y[j]=(X[k]>>M1 & 1);
                j=j-1; }}
        for(int i=0; i<n; ++i)
            sum=sum+Y[i]*pow(2,i);
        //cout<<sum<<",";
        body[i].Hilbert_dis=sum;}
    cout<<"\n";
    sort_merge(A,body,0,N,N);
    
    for(int i=0; i<N; ++i)
        cout<<body[A[i]].ID1<<",";

    return 0;}
