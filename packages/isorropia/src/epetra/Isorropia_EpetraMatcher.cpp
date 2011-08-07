//============================================================================
// Name        : Isorropia_EpetraMatcher.cpp
// Author      : Arif Khan
// Email       : khan58@purdue.edu
// Copyright   : Sandia National Labs
// Description : Parallel Matching in C++/openMP
//============================================================================

#include"Isorropia_EpetraMatcher.hpp"
#include<algorithm>
using namespace std;

//#define ISORROPIA_MATCHING_STATS

#ifdef ISORROPIA_MATCHING_STATS
int minL,maxL;
vector<int> med;
#endif

Isorropia_EpetraMatcher::Isorropia_EpetraMatcher(const Epetra_CrsMatrix * matrixPtr,const Teuchos::ParameterList& paramlist)
{
    int rc=0,i;
    
    rc=matrixPtr->ExtractCrsDataPointers(CRS_pointers_,CRS_indices_,CRS_vals_);
    
    if(rc!=0)
        cout<<"Input Processing Failed"<<endl;
        
    U_=matrixPtr->NumGlobalRows();
    V_=matrixPtr->NumGlobalCols();
    E_=matrixPtr->NumGlobalNonzeros();
    
    choice_=1;
    if(paramlist.isParameter("PHKDW"))
        choice_=1;
    else 
        if(paramlist.isParameter("PHK"))
            choice_=2;
        else
            if(paramlist.isParameter("PDFS"))
                choice_=3;
            else
                if(paramlist.isParameter("PPF"))
                    choice_=4;
    
#ifdef ISORROPIA_MATCHING_STATS
    cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<endl;
    cout<<"choice: "<<choice_<<endl;
#endif
    
    finish_=false;
    avgDegU_=E_/U_+1;
    mateU_=new int[U_];
    mateV_=new int[V_];
        
    if(choice_==1 || choice_==2)
    {
        LU_=new int[U_];
        LV_=new int[V_];
        Queue_=new int[U_];
    }
    if(choice_==1 || choice_==4)
        lookahead_=new int[U_];
    
    if(choice_==3 || choice_==4)
        unmatchedU_=new int[U_];
    
    for(i=0;i<U_;i++)
    {   
        mateU_[i]=-1;
        if(choice_==1 || choice_==4)
            lookahead_[i]=CRS_pointers_[i];
        if(choice_==1 || choice_==2)
        {   
            LU_[i]=-1;
            Queue_[i]=-1;
        }
        if(choice_==3||choice_==4)
            unmatchedU_[i]=i;
    }
    
#ifdef ISORROPIA_HAVE_OMP
    scannedV_=new omp_lock_t[V_];
#endif
    parent_=new int[V_];
    
    for(i=0;i<V_;i++)
    {   
#ifdef ISORROPIA_HAVE_OMP
        omp_init_lock(&scannedV_[i]);
#endif
        mateV_[i]=-1;
        parent_[i]=-1;
        if(choice_==1 || choice_==2)
            LV_[i]=-1;
    }
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel
    numThread_=omp_get_num_threads();
#endif
}

Isorropia_EpetraMatcher::Isorropia_EpetraMatcher(Teuchos::RCP<const Epetra_CrsMatrix> matrixPtr,const Teuchos::ParameterList& paramlist)
{
    int rc=0,i;
    
    rc=matrixPtr->ExtractCrsDataPointers(CRS_pointers_,CRS_indices_,CRS_vals_);
    
    if(rc==0)
        cout<<"Input Processing Done"<<endl;
    else
        cout<<"Input Processing Failed"<<endl;
        
    U_=matrixPtr->NumGlobalRows();
    V_=matrixPtr->NumGlobalCols();
    E_=matrixPtr->NumGlobalNonzeros();
    cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<endl;
    
    choice_=1;
    if(paramlist.isParameter("PHKDW"))
        choice_=1;
    else 
        if(paramlist.isParameter("PHK"))
            choice_=2;
        else
            if(paramlist.isParameter("PDFS"))
                choice_=3;
            else
                if(paramlist.isParameter("PPF"))
                    choice_=4;
    
    cout<<"choice: "<<choice_<<endl;
    finish_=false;
    avgDegU_=E_/U_+1;
    mateU_=new int[U_];
    mateV_=new int[V_];
        
    if(choice_==1 || choice_==2)
    {
        LU_=new int[U_];
        LV_=new int[V_];
        Queue_=new int[U_];
    }
    if(choice_==1 || choice_==4)
        lookahead_=new int[U_];
    
    if(choice_==3 || choice_==4)
        unmatchedU_=new int[U_];
    
    for(i=0;i<U_;i++)
    {   
        mateU_[i]=-1;
        if(choice_==1 || choice_==4)
            lookahead_[i]=CRS_pointers_[i];
        if(choice_==1 || choice_==2)
        {   
            LU_[i]=-1;
            Queue_[i]=-1;
        }
        if(choice_==3||choice_==4)
            unmatchedU_[i]=i;
    }
    
#ifdef ISORROPIA_HAVE_OMP
    scannedV_=new omp_lock_t[V_];
#endif
    parent_=new int[V_];
    
    for(i=0;i<V_;i++)
    {   
#ifdef ISORROPIA_HAVE_OMP
        omp_init_lock(&scannedV_[i]);
#endif
        mateV_[i]=-1;
        parent_[i]=-1;
        if(choice_==1 || choice_==2)
            LV_[i]=-1;
    }
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel
    numThread_=omp_get_num_threads();
#endif
}

Isorropia_EpetraMatcher::Isorropia_EpetraMatcher(const Epetra_CrsGraph * graphPtr,const Teuchos::ParameterList& paramlist)
{
    
}

Isorropia_EpetraMatcher::Isorropia_EpetraMatcher(Teuchos::RCP<const Epetra_CrsGraph> graphPtr,const Teuchos::ParameterList& paramlist)
{
    
}


Isorropia_EpetraMatcher::~Isorropia_EpetraMatcher() 
{
    delete [] mateU_;
    delete [] mateV_;
    
    if(choice_==1 || choice_==2)
    {   
        delete [] LU_;
        delete [] LV_;
        delete [] Queue_;
    }
    if(choice_==1 ||choice_==4)
        delete [] lookahead_;
    if(choice_==3 ||choice_==4)
        delete [] unmatchedU_;
    
#ifdef ISORROPIA_HAVE_OMP
    for(int i=0;i<V_;i++)
        omp_destroy_lock(&scannedV_[i]);
#endif
    
}

void Isorropia_EpetraMatcher::extractRowPermutationCopy(int len, int& size, int* array) const
{
    const int *ptr=&mateU_[0];
    size=MIN(size,len);
    memcpy (array, ptr, size * sizeof(int));
}
void Isorropia_EpetraMatcher::extractColumnPermutationCopy(int len, int& size, int* array) const
{
    const int *ptr=&mateV_[0];
    size=MIN(size,len);
    memcpy (array, ptr, size * sizeof(int));
}

void Isorropia_EpetraMatcher::getMatchedEdges(int len,int& size,int* array) const
{
    int i,j;
    j=0;
    for(i=0;i<U_;i++)
    {   
        if(mateU_[i]!=-1)
        {   
            array[j]=i;
            array[j+1]=mateU_[i];
            j++;
        }
    }
}
int Isorropia_EpetraMatcher::getNumberOfMatchedVertices()
{
    return 2*matched_;
}
Epetra_Map* Isorropia_EpetraMatcher::getPermutedRowMap()
{
    return NULL;
}
Epetra_Map* Isorropia_EpetraMatcher::getPermutedColumnMap()
{
    return NULL;
}

void Isorropia_EpetraMatcher::filler()
{
    int i,j,rowfill,colfill,flag;
    
    vector<int> temp(U_,-1);
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<U_;i++)
        temp[i]=mateU_[i];
    
    rowfill=U_-matched_;
    
    j=0;
    while(rowfill>0)
    {
        for(i=0;i<U_;i++)
        {   
            if(mateU_[i]==-1)
            {
                flag=0;
                for(;j<V_;j++)
                {   
                    if(mateV_[j]==-1)
                    {
                        mateU_[i]=j;
                        j++;
                        rowfill--;
                        flag=1;
                        break;
                    }
                }
                if(flag==0)         
                {
                    mateU_[i]=j;
                    j++;
                    rowfill--;
                }
            }
        }
    }
    
    colfill=V_-matched_;
    
    j=0;
    while(colfill>0)
    {
        for(i=0;i<V_;i++)
        {   
            if(mateV_[i]==-1)
            {
                flag=0;
                for(;j<U_;j++)
                {   
                    if(temp[j]==-1)
                    {
                        mateV_[i]=j;
                        j++;
                        colfill--;
                        flag=1;
                        break;
                    }
                }
                if(flag==0)         
                {
                    mateV_[i]=j;
                    j++;
                    colfill--;
                }
            }
        }
    }
}

void Isorropia_EpetraMatcher::delete_matched_v()
{
#ifdef ISORROPIA_HAVE_OMP
    int i,j;
    #pragma omp parallel for private(j)
    for(i=Qst_;i<Qend_;i++)
    {   
        j=Queue_[i];
        if(LV_[j]==k_star_ && mateV_[j]!=-1)
            omp_test_lock(&scannedV_[j]);
    }
#endif
}

int Isorropia_EpetraMatcher::augment_matching(int tv)
{
    int u,v,t,lnt=1;
    v=tv;
    while(true)
    {
        u=parent_[v];
        t=mateU_[u];
        mateV_[v]=u;
        mateU_[u]=v;
        if(t==-1)
            break;
        lnt+=2;
        v=t;
    }
    
    return lnt;
}

int Isorropia_EpetraMatcher::construct_layered_graph()
{
    int k,i,j,t,tst,tend,fflag,s,tid,mem,pqind;
    Qst_=Qend_=tst=tend=k=fflag=0;
#ifdef ISORROPIA_MATCHING_STATS
    maxL=0;
    minL=U_+V_+1;
#endif
    int* Qsize=new int[numThread_];
    int** localQ=new int*[numThread_];
    int* startInd=new int[numThread_];
    
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<V_;i++)
    {   
#ifdef ISORROPIA_HAVE_OMP
        omp_unset_lock(&scannedV_[i]);
#endif
        LV_[i]=-1;
    }
    
    for(i=0;i<U_;i++)    
    {
        if(mateU_[i]==-1)
        {
            LU_[i]=0;
            Queue_[Qend_++]=i;
        }
        else
            LU_[i]=-1;
    }
    
    if((Qend_-Qst_)==0)
    {   
        finish_=true;
        return 0;
    }
    else
        BFSInd_=Qend_;
    
    while(true)
    {
        //mem=MIN(((avgDegU_*(Qend_-Qst_))/numThread_),V_);
        mem=V_;
        for(i=0;i<numThread_;i++)
        {
            startInd[i]=0;
            Qsize[i]=mem;
            localQ[i]=new int[mem];
        }
        
#ifdef ISORROPIA_HAVE_OMP
        #pragma omp parallel for private(i,t,j,pqind,tid)
#endif
        for(s=Qst_;s<Qend_;s++)
        {
#ifdef ISORROPIA_HAVE_OMP
            tid=omp_get_thread_num();
#endif
            
            pqind=startInd[tid];
                    
            i=Queue_[s];
            for(t=CRS_pointers_[i];t<CRS_pointers_[i+1];t++) 
            {
                j=CRS_indices_[t];                      
                if(mateU_[i]!=j)
                {   
                    if(mateV_[j]==-1)
                    {   
                        fflag=1;
                        LV_[j]=k+1;
                    }
                    else
                    {
                        if(fflag==0 && (LV_[j]==-1 || LV_[j]==k+1))
                        {
                            if(LV_[j]==-1)
                            {   
                                if(Qsize[tid]==pqind)
                                {
                                    int newsize=pqind+avgDegU_;
                                    Qsize[tid]=newsize;
                                    int * temp=new int[newsize];
                                    for(i=0;i<pqind;i++)
                                        temp[i]=localQ[tid][i];
                                    delete [] localQ[tid];
                                    localQ[tid]=temp;
                                        
                                }
                                    
                                localQ[tid][pqind++]=mateV_[j];
                            }
                            LU_[mateV_[j]]=k+2;
                            LV_[j]=k+1;
                        }
                    }
                }   
            }
#ifdef ISORROPIA_HAVE_OMP
            #pragma omp flush
#endif
            startInd[tid]=pqind;
        }   
        
        Qst_=Qend_;
        tst=Qst_;
        for(int ii=0;ii<numThread_;ii++)
            for(int jj=0;jj<startInd[ii];jj++)
                Queue_[tst++]=localQ[ii][jj];
        Qend_=tst;  
        
        for(i=0;i<numThread_;i++)
            delete [] localQ[i];
        
        if(fflag>0)
        {   
            k_star_=k+1;
            break;
        }
        else
            k=k+2;
        
        if((Qend_-Qst_)==0)
        {
            finish_=true;
            break;
        }
    }
    
    delete [] startInd;
    delete [] Qsize;
    return k_star_;
}


int Isorropia_EpetraMatcher::recursive_path_finder(int k, int p)
{
    int i,ind,res=0,lock=0;
    
    if(k>k_star_)
        return -1;
            
    for(i=CRS_pointers_[p];i<CRS_pointers_[p+1];i++)
    {
        ind=CRS_indices_[i];
        if(LV_[ind]==k+1)
        {
#ifdef ISORROPIA_HAVE_OMP
            lock=omp_test_lock(&scannedV_[ind]);
#endif
            if(lock>0)
            {
                parent_[ind]=p;
                if(mateV_[ind]==-1)
                    return ind;
                else
                {
                    res=recursive_path_finder(k+2,mateV_[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }
    }
    
    return -1;
}

int Isorropia_EpetraMatcher::dfs_path_finder(int u)
{
    int i,ind=-1,res=0,lock=0;
    
    if(choice_==1 || choice_==4)
    {
        for(i=lookahead_[u];i<CRS_pointers_[u+1];i++)
        {
            ind=CRS_indices_[i];
            assert(ind>=0 && ind<V_);
            if(mateV_[ind]==-1)
            {
#ifdef ISORROPIA_HAVE_OMP
                lock=omp_test_lock(&scannedV_[ind]);
#endif
                if(lock>0)
                {
                    parent_[ind]=u;
                    lookahead_[u]=i+1;
                    return ind;
                }
            }
        }
    
    
        if(icm_%2==1)   
        {   
            for(i=CRS_pointers_[u];i<CRS_pointers_[u+1];i++)
            {
                ind=CRS_indices_[i];
                assert(ind>=0 && ind<V_);
#ifdef ISORROPIA_HAVE_OMP
                lock=omp_test_lock(&scannedV_[ind]);
#endif
                if(lock>0)
                {
                    parent_[ind]=u;
                    res=dfs_path_finder(mateV_[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }
        else
        {
            for(i=CRS_pointers_[u+1]-1;i>=CRS_pointers_[u];i--)
            {
                ind=CRS_indices_[i];
                assert(ind>=0 && ind<V_);
#ifdef ISORROPIA_HAVE_OMP
                lock=omp_test_lock(&scannedV_[ind]);
#endif
                if(lock>0)
                {
                    parent_[ind]=u;
                    res=dfs_path_finder(mateV_[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }
    }
    else
    {   
        for(i=CRS_pointers_[u];i<CRS_pointers_[u+1];i++)
        {
            ind=CRS_indices_[i];
            assert (ind>=0 && ind<V_);
#ifdef ISORROPIA_HAVE_OMP
            lock=omp_test_lock(&scannedV_[ind]);
#endif
            if(lock>0)
            {
                parent_[ind]=u;
                if(mateV_[ind]==-1)
                    return ind;
                else
                {
                    res=dfs_path_finder(mateV_[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }
    }
    return -1;
}

int Isorropia_EpetraMatcher::find_set_del_M()
{
    int i,j,count=0;
    delete_matched_v();
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for private(j)
#endif
    for(i=0;i<BFSInd_;i++)
    {
        
        j=recursive_path_finder(0,Queue_[i]);
        if(j!=-1)
        {   
            int lnt=augment_matching(j);
#ifdef ISORROPIA_MATCHING_STATS
#ifdef ISORROPIA_HAVE_OMP
            #pragma omp atomic
#endif
            count++;
#ifdef ISORROPIA_HAVE_OMP
            #pragma omp critical
#endif
            {
                maxL=maxL>lnt?maxL:lnt;
                minL=minL<lnt?minL:lnt;
                med.push_back(lnt);
            }
            #endif
        }
                
    }
    return count;
}

int Isorropia_EpetraMatcher::DW_phase()
{
    int i,count=0;
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<V_;i++)
    {   
#ifdef ISORROPIA_HAVE_OMP
        omp_unset_lock(&scannedV_[i]);
#endif
    }
    
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<BFSInd_;i++)
    {
        int u=Queue_[i];
        if(mateU_[u]==-1)
        {   
            int ind=dfs_path_finder(u);
            if(ind!=-1)
            {   
                int lnt=augment_matching(ind);
#ifdef ISORROPIA_MATCHING_STATS
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp atomic
#endif
                count++;
                
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp critical
#endif
                {
                    maxL=maxL>lnt?maxL:lnt;
                    minL=minL<lnt?minL:lnt;
                    med.push_back(lnt);
                }
#endif
            }
        }
    }
    return count;
}

int Isorropia_EpetraMatcher::dfs_augment()
{

    int i,flag=0,flag1=0,count=0,totc=0,index=U_;
    icm_=0;
#ifdef ISORROPIA_MATCHING_STATS
    vector<int> med;
#endif
    while(true)
    {
        icm_++;
        
#ifdef ISORROPIA_HAVE_OMP
        #pragma omp parallel for
#endif
        for(i=0;i<V_;i++)
        {   
#ifdef ISORROPIA_HAVE_OMP
            omp_unset_lock(&scannedV_[i]);
#endif
        }
            
        
        flag=flag1=count=0;
#ifdef ISORROPIA_MATCHING_STATS
        maxL=0;
        minL=U_+V_+1;
#endif
        
#ifdef ISORROPIA_HAVE_OMP
        #pragma omp parallel for
#endif
        for(i=0;i<index;i++)
        {
            
#ifdef ISORROPIA_HAVE_OMP
            #pragma omp flush
#endif
            flag=1;
            int u=unmatchedU_[i];
            int ind=dfs_path_finder(u);
            if(ind!=-1)
            {   
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp flush
#endif
                flag1=1;
                int lnt=augment_matching(ind);
#ifdef ISORROPIA_MATCHING_STATS
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp atomic
#endif
                count++;
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp critical
#endif
                {
                    maxL=maxL>lnt?maxL:lnt;
                    minL=minL<lnt?minL:lnt;
                    med.push_back(lnt);
                }
#endif
            }
        }
            
        if(flag==0 || flag1==0)
            break;
        else
        {   
#ifdef ISORROPIA_MATCHING_STATS
            sort(med.begin(),med.end());
            totc+=count;
            //cout<<"["<<icm_<<"] unmatched="<<index<<" matched="<<count<<" size="<<totc<<" minL= "<<minL<<" maxL="<<maxL<<" medL="<<med[count/2]<<endl;
            cout<<icm_<<","<<index<<","<<count<<","<<(index*1.0)/count<<","<<minL<<","<<med[count/2]<<","<<maxL<<endl;
            med.clear();
#endif
            
            int cur=0;
            
            for(i=0;i<index;i++)
            {
                if(mateU_[unmatchedU_[i]]==-1)
                    unmatchedU_[cur++]=unmatchedU_[i];
            }
            index=cur;
            
        }
    }
    
    return totc;
}

int Isorropia_EpetraMatcher::match_dfs()
{
    int totc=0;
    icm_=0;
    double start,end;
#ifdef ISORROPIA_HAVE_OMP
    start=omp_get_wtime();
#endif
    totc=dfs_augment();
#ifdef ISORROPIA_HAVE_OMP
    end=omp_get_wtime();
#endif
    
        
#ifdef ISORROPIA_MATCHING_STATS
    cout<<"Total time: "<<(end-start)<<" seconds"<<" matching="<<totc<<endl;
#endif
    return 0;
}

int Isorropia_EpetraMatcher::match_hk()
{
    int totc=0,count=0;
    icm_=0;
    double start,end;
#ifdef ISORROPIA_HAVE_OMP
    start=omp_get_wtime();
#endif
    
    
    while(true)
    {
        icm_++;
        //cout<<"bfs"<<endl;
        construct_layered_graph();
        if(finish_)
            break;
        //cout<<"dfs"<<endl;
        count=find_set_del_M();
        if(choice_==1)
        {   
            //cout<<"dw"<<endl;
            count+=DW_phase();
        }
#ifdef ISORROPIA_MATCHING_STATS
        totc+=count;
        sort(med.begin(),med.end());
        //cout<<"["<<icm_<<"] unmatched="<<BFSInd_<<" matched="<<count<<" size="<<totc<<" minL= "<<minL<<" maxL="<<maxL<<" medL="<<med[count/2]<<endl;
        cout<<icm_<<","<<BFSInd_<<","<<count<<","<<(BFSInd_*1.0)/count<<","<<minL<<","<<med[count/2]<<","<<maxL<<endl;
        med.clear();
#endif
    }
#ifdef ISORROPIA_HAVE_OMP
    end=omp_get_wtime();
#endif
    
    
#ifdef ISORROPIA_MATCHING_STATS
    cout<<"Total time: "<<(end-start)<<" seconds"<<" matching="<<totc<<endl;
#endif
    return 0;
}

int Isorropia_EpetraMatcher::match()
{
    switch(choice_)
    {
        case 1:match_hk(); 
                 break;
        case 2:match_hk(); 
                 break;
        case 3:match_dfs(); 
                 break;
        case 4:match_dfs(); 
                 break;
        default:match_hk();
    }
    return 0;
}


