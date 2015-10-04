#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <random>

float th=0.5;
float th_cntr=0.5;

std::vector<std::vector<int>> cluster_groups(int, std::vector<std::vector<float>>&,
                                             std::vector<std::vector<float>>&, std::vector<int>&);
std::vector<std::vector<int>> cluster(int,std::vector<std::vector<float>>&,
                                      std::vector<std::vector<float>>&,std::vector<int>&,
                                      std::vector<std::vector<int>> &);
void First_Cluster(std::vector<std::vector<float>>&);

void ReadInput(std::vector<std::vector<float>>& data){
    std::ifstream myfile ("clusterdata2.csv");
    std::string line;
    std::vector<float> row(2,0.0);
    std::vector<std::string> values;
    std::string s1("0.0"), s2("0.0");
    
    while(std::getline(myfile,line)){
        std::istringstream s(line);
        std::string value;
        
        while(std::getline(s,value,',')){
            int flag=0; int k=0;
            for(int i=0; i<value.size(); ++i){
                if(value[i]!=' ' && value[i]!='\r' && value[i]!='\n' ){
                    s1[k]=value[i]; k=k+1;}
            }
            
            for(int i=0; i<value.size(); ++i){
                if(value[i]=='\r' || value[i]=='\n'){
                    k=i+1; flag=1;}
            }
            
            if(flag==1){
                for(int i=k; i<value.size(); ++i)
                    s2[i-k]=value[i];}
            
            values.push_back(s1);
            if(flag==1) values.push_back(s2);
        }
    }
    
    int i=0;
    while(i<values.size()){
        for(int j=0; j<2; ++j)
            row[j]=::atof(values[i+j].c_str());
        data.push_back(row);
        i=i+2;}
}

std::vector<float> New_Centroid(std::vector<int> & groups,
                                std::vector<std::vector<float>> & data){
    std::vector<float> ps(2,0.0);
    for(int j=0; j<groups.size(); ++j){
        ps[0]=ps[0]+data[groups[j]][0];
        ps[1]=ps[1]+data[groups[j]][1];
    }
    if(groups.size()!=0){
        ps[0]=ps[0]/groups.size();
        ps[1]=ps[1]/groups.size();}
    return ps;
}

std::vector<float> dist(int k,
                        std::vector<std::vector<float>> &centroids,
                        std::vector<std::vector<int>> &groups,
                        std::vector<std::vector<float>> &data){
    
    std::vector<float> d(k,0.0);
    
    for(int i=0; i<k; ++i)
        for(int j=0; j<k; ++j){
            float w2=float(groups[i].size()+groups[j].size())/float(data.size());
            d[i]=d[i]+(pow((centroids[j][0]-centroids[i][0]),2.0)+pow((centroids[j][1]-centroids[i][1]),2.0))*w2;
    }
    return d;
}

std::vector<float> SSE(int k,
                       std::vector<std::vector<float>> &centroids,
                       std::vector<std::vector<int>> &groups,
                       std::vector<std::vector<float>> &data){
    
    std::vector<float> error(k,0.0);
    for(int i=0; i<k; ++i)
        for(int j=0; j<groups[i].size(); ++j)
            error[i]=error[i]+pow((centroids[i][0]-data[groups[i][j]][0]),2.0)+pow((centroids[i][1]-data[groups[i][j]][1]),2.0);
    
    return error;
}

std::vector<std::vector<int>> New_Cluster(int k,
                                      std::vector<std::vector<float>>& centroids,
                                      std::vector<std::vector<float>> &data,
                                      std::vector<int> &ID_groups,
                                      std::vector<std::vector<int>> &groups){
    
    std::vector<std::vector<float>> new_centroids;
    std::vector<std::vector<int>> new_groups,new_groups1;
    std::vector<float> ps(2,0.0);
    std::vector<float> er_cntr;
    
    for(int i=0; i<k; ++i){
        er_cntr.push_back(0);
        new_centroids.push_back(ps);}
    
    int flag2=0;
    for(int i=0; i<k; ++i){
        ps=New_Centroid(groups[i],data);
        new_centroids[i]=ps;}
    
    new_groups=cluster_groups(k,new_centroids,data,ID_groups);
    
    for(int i=0; i<k; ++i){ //re-run based on error in centroids` positions
        er_cntr[i]=sqrt(pow((centroids[i][0]-new_centroids[i][0]),2.0)+
                        pow((centroids[i][1]-new_centroids[i][1]),2.0));
        if(er_cntr[i]>th_cntr)
            flag2=1;
     }
    
    if(flag2==1){ //should run again clustering with new centroids
        new_groups1=New_Cluster(k,new_centroids,data,ID_groups,new_groups);
        return new_groups1;
    }

    return new_groups;
    
}

std::vector<std::vector<int>> Split_Cluster(int k,
                                            std::vector<std::vector<float>>& centroids,
                                            std::vector<std::vector<int>>& groups,
                                            std::vector<std::vector<float>> &data){
    
    std::vector<float> sse,d;
    std::vector<int> grp, ID_groups;
    std::vector<std::vector<int>> new_groups,new_groups1,splits,new_splits;
    std::vector<std::vector<float>> new_cntr;
    std::vector<float> ps(2,0.0);
    float obj_fnc1=0, obj_fnc2=0;
 
    for(int i=0; i<k; ++i){
        d.push_back(0);
        sse.push_back(0);}
    
    for(int i=0; i<2; ++i) //two more new centroids
        new_cntr.push_back(ps);
    
    d=dist(k,centroids,groups,data); //dist[i]
    sse=SSE(k,centroids,groups,data); //sse[i]

    for(int i=0; i<k; ++i){
        float w=float(groups[i].size())/float(data.size());
        sse[i]=sse[i]*w;
        obj_fnc1=obj_fnc1+sse[i];
        obj_fnc2=obj_fnc2+0.5*d[i];
    }
    std::cout<<"\r new obj_fnc values after spliting: \r{"<<obj_fnc1<<","<<obj_fnc2<<"}"<<std::endl;
    
    float max=0; int cluster_for_split=0;
    for(int i=0; i<k; ++i)
        if(sse[i]>max && groups[i].size()>=3){
            cluster_for_split=i;
            max=sse[i];
        }

    for(int i=0; i<k; ++i){
        if(i != cluster_for_split){
            new_groups.push_back(groups[i]);
        }
    }

    for(int i=0; i<2; ++i){ //rand within the group!
        int rand_num=rand() % groups[cluster_for_split].size();
        new_cntr[i][0]=data[groups[cluster_for_split][rand_num]][0];
        new_cntr[i][1]=data[groups[cluster_for_split][rand_num]][1];
    }
    
    new_splits=cluster_groups(2,new_cntr,data,groups[cluster_for_split]);
    splits=New_Cluster(2,new_cntr,data,groups[cluster_for_split],new_splits);
    
    new_groups.push_back(splits[0]);
    new_groups.push_back(splits[1]);

    std::cout<<"\r size of new clusters is: "<<new_groups.size()<<std::endl;
    
    int code;
    std::cout<<"\r Do you want to split one more time for cluster number " << cluster_for_split<<" [1/0]?"  <<std::endl;
    std::cin>>code;
    if(code==1){
        
        for(int i=0; i<k-1; ++i) //two more new centroids
            new_cntr.push_back(ps);
        
        for(int i=0; i<k+1; ++i){
            ps=New_Centroid(new_groups[i],data);
            new_cntr[i]=ps;
        }
        new_groups1=Split_Cluster(k+1,new_cntr,new_groups,data);
        
        return new_groups1;
    }
    
    return groups;
}

std::vector<int> Merge_Study(int k,
                             std::vector<std::vector<int>>& groups,
                             std::vector<std::vector<float>>& centroids){

    float d2; std::vector<int> cell_merged(2,0);
    float min=1000;
    
    for(int i=0; i<k; ++i)
        for(int j=0; j<k; ++j){
            float w2=float(groups[i].size()+groups[j].size());
            d2=(pow((centroids[j][0]-centroids[i][0]),2.0)+pow((centroids[j][1]-centroids[i][1]),2.0))*w2;
            if (d2<min && d2!=0){
                min=d2; cell_merged[0]=i; cell_merged[1]=j;
            }
        }
    
    return cell_merged;
}

std::vector<std::vector<int>> Merge_Cluster(int k,
                                               std::vector<std::vector<float>>& centroids,
                                               std::vector<std::vector<int>>& groups,
                                               std::vector<std::vector<float>> &data){
    
    std::vector<float> sse,d;
    std::vector<int> grp, ID_groups;
    std::vector<std::vector<int>> new_groups,new_groups1,splits,new_splits;
    std::vector<std::vector<float>> new_cntr;
    std::vector<float> ps(2,0.0);
    float obj_fnc1=0, obj_fnc2=0;
    
    for(int i=0; i<k; ++i){
        d.push_back(0);
        sse.push_back(0);}
   
    d=dist(k,centroids,groups,data); //dist[i]
    sse=SSE(k,centroids,groups,data); //sse[i]
    
    for(int i=0; i<k; ++i){
        float w=float(groups[i].size())/float(data.size());
        sse[i]=sse[i]*w;
        obj_fnc1=obj_fnc1+sse[i];
        obj_fnc2=obj_fnc2+0.5*d[i];
    }
    std::cout<<"\r new obj_fnc values after merging:{"<<obj_fnc1<<","<<obj_fnc2<<"}"<<std::endl;
    
    std::vector<int> cell_merged=Merge_Study(k,groups,centroids);
    
    for(int i=0; i<k-1; ++i)
        new_cntr.push_back(ps);
    int t=0;
    for(int i=0; i<k; ++i){
        if(i != cell_merged[0] && i != cell_merged[1]){
            new_groups.push_back(groups[i]);
            new_cntr[t]=centroids[i]; t=t+1;
        }
    }
  
    for(int i=0; i<groups[cell_merged[0]].size(); ++i)
        groups[cell_merged[1]].push_back(groups[cell_merged[0]][i]);
   
    new_groups.push_back(groups[cell_merged[1]]);
    new_cntr[k-2]=New_Centroid(new_groups[k-2],data);
    
    for(int i=0; i<k-1; ++i)
        std::cout<<new_groups[i].size()<<",";
    
    std::cout<<"\r size of new clusters is: "<<new_groups.size()<<std::endl;
    
    int code;
    std::cout<<"\r Do you want to merge one more time? [1/0]"  <<std::endl;
    std::cin>>code;
    if(code==1){
        new_groups1=Merge_Cluster(k-1,new_cntr,new_groups,data);
        return new_groups1;
    }
    
    return groups;
}


std::vector<std::vector<int>> cluster_groups(int k,
                                      std::vector<std::vector<float>>& centroids,
                                      std::vector<std::vector<float>> &data,
                                      std::vector<int> &ID_groups){
    
    int size=int(ID_groups.size());
    std::vector<float> ps(2,0.0);
    std::vector<int> grp;
    std::vector<float> er,er_cntr;
    std::vector<std::vector<int>> groups;
    
    for(int i=0; i<k; ++i){
        groups.push_back(grp);
        er.push_back(0);
    }
    
    for(int i=0; i<size; ++i){
        float min=1000;
        int flag=0;
        
        for(int j=0; j<k; ++j){
            er[j]=sqrt(pow((data[ID_groups[i]][0]-centroids[j][0]),2.0)+pow((data[ID_groups[i]][1]-centroids[j][1]),2.0));
            
            if(er[j] < min)
                min=er[j];}
        
        for(int j=0; j<k; ++j)
            if(min==er[j] && flag==0){
                groups[j].push_back(ID_groups[i]);
                flag=1;}}

    return groups;
}

std::vector<std::vector<int>> cluster(int k,
                                      std::vector<std::vector<float>>& centroids,
                                      std::vector<std::vector<float>> &data,
                                      std::vector<int> &ID_groups,
                                      std::vector<std::vector<int>> &groups){

    std::vector<std::vector<float>> new_centroids;
    std::vector<std::vector<int>> new_groups, new_groups1, new_groups2;
    std::vector<float> ps(2,0.0);
    std::vector<float> er_cntr;
    
    for(int i=0; i<k; ++i){
        er_cntr.push_back(0);
        new_centroids.push_back(ps);
    }
    
    int flag2=0;
    for(int i=0; i<k; ++i){ //computing centroid of each cluster:
        ps=New_Centroid(groups[i],data);
        new_centroids[i]=ps;
        //std::cout<<"\r new centroid is: ["<<new_centroids[i][0]<<","<<new_centroids[i][1]<<"]";
    }
 
    std::cout<<std::endl;
    
    float obj_fnc1=0, obj_fnc2=0;
    std::vector<float> sse,d;
    for(int i=0; i<k; ++i)
        sse.push_back(0);
    
    d=dist(k,centroids,groups,data);
    sse=SSE(k,centroids,groups,data);
    
    for(int i=0; i<k; ++i){
        float w=float(groups[i].size())/float(data.size());
        obj_fnc1=obj_fnc1+sse[i]*w;
        obj_fnc2=obj_fnc2+0.5*d[i];
    }
    std::cout<<"{"<<obj_fnc1<<","<<obj_fnc2<<"}"<<std::endl;
    float ratio=obj_fnc2/obj_fnc1;
    
    if(ratio<0.9)//re-run based on error in sse/dist
        flag2=1;
  
    if(flag2==1){ //should run again clustering with new centroids
        new_groups=cluster_groups(k,new_centroids,data,ID_groups);
        new_groups1=cluster(k,new_centroids,data,ID_groups,new_groups);}
    
    else{
        int code;
        std::cout<<" \r For spliting more Enter 0, For merging more Enter 1";
        std::cin>>code;
        if(code==0){
            if(flag2==1)
                new_groups1=Split_Cluster(k,centroids,new_groups1,data);
            else
                new_groups1=Split_Cluster(k,centroids,groups,data);}
        else{
            if(flag2==1)
                new_groups1=Merge_Cluster(k,centroids,new_groups1,data);
            else
                new_groups1=Merge_Cluster(k,centroids,groups,data);
            }
     return new_groups1;
    
    }
    
    return groups; //will return clusters with best k-centroids

}


void First_Cluster(std::vector<std::vector<float>> &data, std::vector<int> &ID_groups){
    
    int size=int(data.size());
    std::vector<int> rand_num;
    std::vector<std::vector<int>> groups, new_groups;
    
    //Splitting method:
    int k;
    std::cout<<"Please enter the number for clustering: "<<std::endl;
    std::cin >> k;
    
    std::vector<float> ps(2,0.0);
    std::vector<std::vector<float>> centroids;
    for(int i=0; i<k; ++i)
        centroids.push_back(ps);
    
    for(int i=0; i<k; ++i){
        rand_num.push_back(rand() % size);
        centroids[i][0]=data[rand_num[i]][0];
        centroids[i][1]=data[rand_num[i]][1];
        std::cout<<"["<<data[rand_num[i]][0]<<","<<data[rand_num[i]][1]<<"]";
    }
    groups=cluster_groups(k,centroids,data,ID_groups);
    
    for(int i=0; i<k; ++i)
        std::cout<<groups[i].size()<<",";
    
    new_groups=cluster(k,centroids,data,ID_groups,groups);
    std::cout<<std::endl;

}

int main(int argc, const char * argv[]) {
    std::vector<std::vector<float>> data;
    std::vector<int> ID_groups;

    ReadInput(data);
    
    for(int i=0; i<data.size(); ++i)
        ID_groups.push_back(i);
    
    First_Cluster(data, ID_groups);

    return 0;
}





