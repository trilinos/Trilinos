//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//************************************************************************
//@HEADER


#include"Isorropia_EpetraMatcher.hpp"

namespace Isorropia {
#ifdef HAVE_EPETRA
namespace Epetra {

//#define ISORROPIA_MATCHING_STATS

#ifdef ISORROPIA_MATCHING_STATS
int minL,maxL;
std::vector<int> med;
#endif

Matcher::Matcher(const Epetra_CrsMatrix * matrixPtr,const Teuchos::ParameterList& paramlist)
{
        //Matcher(Teuchos::RCP<const Epetra_CrsMatrix>(matrixPtr,false), paramlist);

        int rc=0,i;
    A_=matrixPtr;
    rc=matrixPtr->ExtractCrsDataPointers(CRS_pointers_,CRS_indices_,CRS_vals_);

    if(rc!=0)
        std::cout<<"Input Processing Failed"<<std::endl;

    U_=matrixPtr->NumGlobalRows();
    V_=matrixPtr->NumGlobalCols();
    E_=matrixPtr->NumGlobalNonzeros();

    choice_=1;
    std::string str(paramlist.get<std::string>("Matching Algorithm"));
    if(str.compare("PHKDW")==0)
        choice_=1;
    else
        if(str.compare("PHK")==0)
            choice_=2;
        else
            if(str.compare("PDFS")==0)
                choice_=3;
            else
                if(str.compare("PPF")==0)
                    choice_=4;

#ifdef ISORROPIA_MATCHING_STATS
    std::cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<std::endl;
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
    //std::cout<<"finished:"<<std::endl;
}

Matcher::Matcher(Teuchos::RCP<const Epetra_CrsMatrix> matrixPtr,const Teuchos::ParameterList& paramlist)
{

    //A_=matrixPtr.getRawPtr();
    //Matcher(A_,paramlist);

    int rc=0,i;
    rc=matrixPtr->ExtractCrsDataPointers(CRS_pointers_,CRS_indices_,CRS_vals_);
    if(rc!=0)
        std::cout<<"Input Processing Failed"<<std::endl;

    U_=matrixPtr->NumGlobalRows();
    V_=matrixPtr->NumGlobalCols();
    E_=matrixPtr->NumGlobalNonzeros();

    choice_=1;
    std::string str(paramlist.get<std::string>("Matching Algorithm"));
    if(str.compare("PHKDW")==0)
        choice_=1;
    else
        if(str.compare("PHK")==0)
            choice_=2;
        else
            if(str.compare("PDFS")==0)
                choice_=3;
            else
                if(str.compare("PPF")==0)
                    choice_=4;

#ifdef ISORROPIA_MATCHING_STATS
    std::cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<std::endl;
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

    std::cout<<"finished:"<<std::endl;
}

Matcher::Matcher(const Epetra_CrsGraph * graphPtr,const Teuchos::ParameterList& paramlist)
{
    U_=graphPtr->NumMyRows();
    V_=graphPtr->NumMyCols();
    E_=graphPtr->NumMyNonzeros();
    int i,j,num,sum, *ind;
    CRS_pointers_=new int[U_+1];
    CRS_indices_=new int[E_];
    sum=0;
    for(i=0;i<U_;i++)
    {
        graphPtr->ExtractMyRowView(i,num,ind);
        CRS_pointers_[i]=sum;
        for(j=0;j<num;j++)
            CRS_indices_[sum+j]=ind[j];
        sum=sum+num;
    }
    CRS_pointers_[U_]=sum;

    choice_=1;
    std::string str(paramlist.get<std::string>("Matching Algorithm"));
    if(str.compare("PHKDW")==0)
        choice_=1;
    else
        if(str.compare("PHK")==0)
            choice_=2;
        else
            if(str.compare("PDFS")==0)
                choice_=3;
            else
                if(str.compare("PPF")==0)
                    choice_=4;

#ifdef ISORROPIA_MATCHING_STATS
    std::cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<std::endl;
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

Matcher::Matcher(Teuchos::RCP<const Epetra_CrsGraph> graphPtr,const Teuchos::ParameterList& paramlist)
{
    U_=graphPtr->NumMyRows();
    V_=graphPtr->NumMyCols();
    E_=graphPtr->NumMyNonzeros();
    int i,j,num,sum, *ind;
    CRS_pointers_=new int[U_+1];
    CRS_indices_=new int[E_];
    sum=0;
    for(i=0;i<U_;i++)
    {
        graphPtr->ExtractMyRowView(i,num,ind);
        CRS_pointers_[i]=sum;
        for(j=0;j<num;j++)
            CRS_indices_[sum+j]=ind[j];
        sum=sum+num;
    }
    CRS_pointers_[U_]=sum;

    choice_=1;
    std::string str(paramlist.get<std::string>("Matching Algorithm"));
    if(str.compare("PHKDW")==0)
        choice_=1;
    else
        if(str.compare("PHK")==0)
            choice_=2;
        else
            if(str.compare("PDFS")==0)
                choice_=3;
            else
                if(str.compare("PPF")==0)
                    choice_=4;

#ifdef ISORROPIA_MATCHING_STATS
    std::cout<<"(U,V,E):"<<U_<<","<<V_<<","<<E_<<std::endl;
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


Matcher::~Matcher()
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

    if (CRS_indices_) { delete [] CRS_indices_; CRS_indices_ = NULL; }
    if (CRS_pointers_) { delete [] CRS_pointers_; CRS_pointers_ = NULL; }
    if (parent_) { delete [] parent_; parent_ = NULL; }

#ifdef ISORROPIA_HAVE_OMP
    for(int i=0;i<V_;i++)
        omp_destroy_lock(&scannedV_[i]);
#endif

}

void Matcher::getMatchedColumnsForRowsCopy(int len, int& size, int* array) const
{
    const int *ptr=&mateU_[0];
    size=MIN(matched_,len);
    memcpy (array, ptr, size * sizeof(int));
}

void Matcher::getMatchedRowsForColumnsCopy(int len, int& size, int* array) const
{
    const int *ptr=&mateV_[0];
    size=MIN(matched_,len);
    memcpy (array, ptr, size * sizeof(int));
}

int Matcher::getNumberOfMatchedVertices()
{
    return matched_;
}

Teuchos::RCP<Epetra_CrsMatrix> Matcher::applyRowPermutation()
{
    //int nmatch = matched_; // suppress compiler warning
    const int *mrows = &mateU_[0];
    //const int *mcols = &mateV_[0]; // suppress compiler warning

    // Create a new matrix with row permutation
    int max_entries = A_->MaxNumEntries();
    // TODO: This row map is ok, but not optimal. Should we create a new row
    // map with mrows[i] ??
    Teuchos::RCP<Epetra_CrsMatrix> perm_matrix = Teuchos::RCP<Epetra_CrsMatrix>
                    (new Epetra_CrsMatrix(Copy, A_->RowMap(), max_entries));
    int n =  A_->NumGlobalRows();

    double *values = new double[max_entries];
    int *indices = new int[max_entries];
    int num_entries;
    for (int i = 0; i < n ; i++)
    {
        // All in serial Comm so 0..n is fine
        A_->ExtractGlobalRowCopy(i, max_entries, num_entries, values,
                                    indices);
        perm_matrix->InsertGlobalValues(mrows[i], num_entries,
                            values, indices);
    }
    perm_matrix->FillComplete();

    delete[] values;
    delete[] indices;
    return (perm_matrix);
}

Teuchos::RCP<Epetra_CrsMatrix> Matcher::applyColumnPermutation()
{
    //int nmatch = matched_; // suppress compiler warning
    //const int *mrows = &mateU_[0]; // suppress compiler warning
    const int *mcols = &mateV_[0];

    // Create a new matrix with column permutation
    int max_entries = A_->MaxNumEntries();
    Teuchos::RCP<Epetra_CrsMatrix> perm_matrix = Teuchos::RCP<Epetra_CrsMatrix>
                    (new Epetra_CrsMatrix(Copy, A_->RowMap(), max_entries));
    int n =  A_->NumGlobalRows();

    double *values = new double[max_entries];
    int *indices = new int[max_entries];
    int num_entries;
    for (int i = 0; i < n ; i++)
    {
        // All in serial Comm so 0..n is fine
        A_->ExtractGlobalRowCopy(i, max_entries, num_entries, values,
                                    indices);
        for (int j = 0; j < num_entries; j++) indices[j]=mcols[indices[j]];
        perm_matrix->InsertGlobalValues(i, num_entries, values, indices);
    }
    perm_matrix->FillComplete();
    delete[] values;
    delete[] indices;
    return (perm_matrix);
}

/* // We don't have any use cases for these two functions now. We might use the
   // getPermutedRowMap()
   // when creating the matrix for applyRowPermutation(). Not exposing them for
   // now.
Epetra_Map* Matcher::getPermutedRowMap()
{
    int *ptr=new int[U_];
    complete_nonperfect_permutation();
    for(int i=0;i<U_;i++)
        ptr[i]=(A_->RowMap()).GID(mateU_[i]);

   if(A_->Comm().NumProc()==1)
   {
        Epetra_Map* map=new Epetra_Map(-1,U_,ptr,0,A_->Comm());
        //std::cout<<"Map created..!!"<<std::endl;
        return map;
   }
   else
   {
        std::cout<<"Original matrix is in distributed memory"<<std::endl;
        return NULL;
   }
}
Epetra_Map* Matcher::getPermutedColumnMap()
{
    int *ptr=new int[V_];
    complete_nonperfect_permutation();
    for(int i=0;i<V_;i++)
        ptr[i]=(A_->ColMap()).GID(mateV_[i]);

    if(A_->Comm().NumProc()==1)
    {
        Epetra_Map* map=new Epetra_Map(-1,V_,ptr,0,A_->Comm());
        return map;
    }
    else
    {
        std::cout<<"Original matrix is in distributed memory"<<std::endl;
        return NULL;
    }
}*/

void Matcher::complete_nonperfect_permutation()
{
    int i,j,rowfill,colfill,flag;

    std::vector<int> temp(U_,-1);

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

void Matcher::delete_matched_v()
{
    // This function only applicable to HK and HKDW.
    // It removes the matched vertices in the last layer of the layered graph.

#ifdef ISORROPIA_HAVE_OMP
    int i,j;
    #pragma omp parallel for private(j)
    for(i=Qst_;i<Qend_;i++)
    {
        j=Queue_[i];
        if(LV_[j]==k_star_ && mateV_[j]!=-1) // k_star is the last layer and
            omp_test_lock(&scannedV_[j]);    // mateV[j]!=-1 means matched..!
    }
#endif
}

int Matcher::augment_matching(int tv)
{
    //This function increases the matching size by one with the help of a
    //augmenting path. The input is an integer which is the id of the last
    //columns vertex of the augmenting path. We trace the whole path by
    //backtraking using parent array. while tracing back we flip the unmatched
    //edges to matched edges.
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

int Matcher::construct_layered_graph()
{
    //This function is used by the HK and HKDW. This is the BFS phase of HK/HKDW
    //where we implicitely build the layered sub graph out of the original
    //graph.

    int k,i,j,t,tst,tend,fflag,s,tid,mem,pqind;
    Qst_=Qend_=tst=tend=k=fflag=0;
#ifdef ISORROPIA_MATCHING_STATS
    maxL=0;
    minL=U_+V_+1;
#endif
    int* Qsize=new int[numThread_]; // holds current localQ size for each thread
    int** localQ=new int*[numThread_]; // localQ for each thread
    int* startInd=new int[numThread_]; // the tail of the localQ for each thread


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

    for(i=0;i<U_;i++)    // Collecting all remaining unmatched row vertices
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
        BFSInd_=Qend_; // BFSInd keeps record of the vertices in Layer 0.

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
        for(s=Qst_;s<Qend_;s++) // start of layer construction
        {
#ifdef ISORROPIA_HAVE_OMP
            tid=omp_get_thread_num();
#else
            // FIXME (mfh 07 Feb 2013) When I found this code, tid was
            // not initialized if ISORROPIA_HAVE_OMP was not defined.
            // I'm not sure what its value should be, but it's
            // reasonable to set it to zero (one thread, whose ID (tid
            // == "thread ID") is zero).
            tid = 0;
#endif
            pqind=startInd[tid];

            i=Queue_[s]; // starting with a unmatched row vertex
            for(t=CRS_pointers_[i];t<CRS_pointers_[i+1];t++) //scanning adj list
            {
                j=CRS_indices_[t];
                if(mateU_[i]!=j) // making sure that it is not already paired
                {
                    if(mateV_[j]==-1) // if unmatched then aug path found
                    {
                        fflag=1;
                        LV_[j]=k+1;
                    }
                    else // else adding it to the layer graph
                    {
                        //making sure that it has not already in the prev layer
                        if(fflag==0 && (LV_[j]==-1 || LV_[j]==k+1))
                        {
                            if(LV_[j]==-1)
                            {
                                if(Qsize[tid]==pqind)
                                {
                                    //localQ for this thread is overloaded so
                                    //add some new space.
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

        // A layer has been constructed now collect all the vertices from localQ
        // which is local to each thread to the global Queue.
        Qst_=Qend_;
        tst=Qst_;
        for(int ii=0;ii<numThread_;ii++)
            for(int jj=0;jj<startInd[ii];jj++)
                Queue_[tst++]=localQ[ii][jj];
        Qend_=tst;

        for(i=0;i<numThread_;i++)
            delete [] localQ[i];

        if(fflag>0) // means we found at least one augmenting path
        {
            k_star_=k+1;
            break;
        }
        else  // otherwise construct the next layer
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


int Matcher::recursive_path_finder(int k, int p)
{
    //This function is used by HK and HKDW. This the DFS phase of the HK/HKDW
    //where we try to find the vertex disjoint augmenting path from the just
    //created layered graph. The input is the layer number, k and the vertex id,
    //p. This function goes down layer by layer (by using k) recursively to find
    //a vertex disjoint path starting from layer 0.
    int i,ind,res=0,lock=0;

    if(k>k_star_)
        return -1;

    for(i=CRS_pointers_[p];i<CRS_pointers_[p+1];i++)
    {
        ind=CRS_indices_[i];
        if(LV_[ind]==k+1) //making sure that the vertex is in the next layer
        {
#ifdef ISORROPIA_HAVE_OMP
            lock=omp_test_lock(&scannedV_[ind]);
#endif
            if(lock>0) // unlocked means this vertex is not part of other paths
            {
                parent_[ind]=p;
                if(mateV_[ind]==-1) // we found a vertex disjoint path
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

int Matcher::dfs_path_finder(int u)
{
    //This function is almost similar to the previous function which is to find
    //a vertex disjoint path. It is used for the algorithm DFS, PPF and in HKDW. The
    //difference is that this function operates on the original graph not on the
    //layerd subgraph. It also does the incorporates the lookahead mechanism and
    //scanning adjacency list alternately from backward and from forward in
    //alternate iterations for PPF and HKDW.

    int i,ind=-1,res=0,lock=0;

    if(choice_==1 || choice_==4) // for HKDW and PPF
    {
        for(i=lookahead_[u];i<CRS_pointers_[u+1];i++) // the lookahead scheme
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


        if(icm_%2==1) // odd number iteration so scan the adj list forward dir
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
        else // even number iteration so scan from backward
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
    else // for DFS.. no lookahead and diferent directional scanning
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

int Matcher::find_set_del_M()
{
    //This function starts the BFS phase for the HK and HKDW. It starts from the
    //layer 0 vertices and tries to find a vertex disjoint path from each of
    //these vertices.

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
#else
            (void) lnt; // suppress compiler warning
#endif // ISORROPIA_MATCHING_STATS
        }

    }
    return count;
}

int Matcher::DW_phase()
{
    //This is the additional Duff and Wiberg phase for the HKDW. This function
    //does nothing but first unset the locks and then runs PPF from the
    //remaining unmatched row vertices after the BFS phase of HK.

    int i,count=0;

#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<V_;i++)
    {
#ifdef ISORROPIA_HAVE_OMP
        omp_unset_lock(&scannedV_[i]); // unsetting the locks
#endif
    }

#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for
#endif
    for(i=0;i<BFSInd_;i++) // calling the PPF
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
#else
                (void) lnt; // suppress compiler warning
#endif // ISORROPIA_MATCHING_STATS
            }
        }
    }
    return count;
}

int Matcher::dfs_augment()
{

    //This function is the starter function for PPF and DFS. It unsets the locks
    //the call dfs_path_finder and then call the augment_matching.

    int i,flag=0,flag1=0,count=0,totc=0,index=U_,cur=0;
    icm_=0;
#ifdef ISORROPIA_MATCHING_STATS
    std::vector<int> med;
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

        cur=0;
        for(i=0;i<U_;i++)
            if(mateU_[i]==-1)
                unmatchedU_[cur++]=i;
        index=cur;
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
#else
                (void) lnt; // suppress compiler warning
#endif // ISORROPIA_MATCHING_STATS
            }
        }

        if(flag==0 || flag1==0)
            break;
        else
        {
#ifdef ISORROPIA_MATCHING_STATS
            sort(med.begin(),med.end());
            totc+=count;
            //std::cout<<"["<<icm_<<"] unmatched="<<index<<" matched="<<count<<" size="<<totc<<" minL= "<<minL<<" maxL="<<maxL<<" medL="<<med[count/2]<<std::endl;
            std::cout<<icm_<<","<<index<<","<<count<<","<<(index*1.0)/count<<","<<minL<<","<<med[count/2]<<","<<maxL<<std::endl;
            med.clear();
#endif

            cur=0;
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

int Matcher::SGM()
{
    int i,j,lock,ind,count=0;
#ifdef ISORROPIA_HAVE_OMP
    #pragma omp parallel for private(j,ind,lock)
#endif
    for(i=0;i<U_;i++)
    {
        for(j=CRS_pointers_[i];j<CRS_pointers_[i+1];j++)
        {
            ind=CRS_indices_[j];
#ifdef ISORROPIA_HAVE_OMP
            lock=omp_test_lock(&scannedV_[ind]);
#else
            // mfh 07 Feb 2013: lock wasn't getting initialized if
            // ISORROPIA_HAVE_OMP was not defined.  omp_test_lock()
            // returns nonzero if the thread successfully acquired the
            // lock.  If there's only one thread, that thread will
            // always be successful, so the default value of lock
            // should be nonzero.
            lock = 1;
#endif
            if(lock>0)
            {

                mateU_[i]=ind;
                mateV_[ind]=i;
#ifdef ISORROPIA_MATCHING_STATS
#ifdef ISORROPIA_HAVE_OMP
                #pragma omp atomic
                count++;
#endif
#endif
                break;
            }
        }
    }
    return count;
}
int Matcher::match_dfs()
{
    // Forking function for DFS based algorithm
    int totc=0;
    icm_=0;

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    double start,end;
    start=omp_get_wtime();
#endif

    totc=dfs_augment();

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    end=omp_get_wtime();
#endif

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    std::cout<<"Total time: "<<(end-start)<<" seconds"<<" matching=";
#endif
    return totc;
}

int Matcher::match_hk()
{
    // Forking function for HK based algorithm
    int totc=0,count=0;

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    double start,end;
#endif

    icm_=0;

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    start=omp_get_wtime();
#endif

    while(true)
    {
        icm_++;
        //std::cout<<"bfs"<<std::endl;
        construct_layered_graph();
        if(finish_)
            break;
        //std::cout<<"dfs"<<std::endl;
        count=find_set_del_M();
        if(choice_==1)
        {
            //std::cout<<"dw"<<std::endl;
            count+=DW_phase();
        }
#ifdef ISORROPIA_MATCHING_STATS
        totc+=count;
        sort(med.begin(),med.end());
        //std::cout<<"["<<icm_<<"] unmatched="<<BFSInd_<<" matched="<<count<<" size="<<totc<<" minL= "<<minL<<" maxL="<<maxL<<" medL="<<med[count/2]<<std::endl;
        std::cout<<icm_<<","<<BFSInd_<<","<<count<<","<<(BFSInd_*1.0)/count<<","<<minL<<","<<med[count/2]<<","<<maxL<<std::endl;
        med.clear();
#endif
    }

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    end=omp_get_wtime();
#endif

#if defined(ISORROPIA_HAVE_OMP) && defined(ISORROPIA_MATCHING_STATS)
    std::cout<<"Total time: "<<(end-start)<<" seconds"<<" matching=";
#endif
    return totc;
}

int Matcher::match()
{
    // User interface function for the matching..

    //std::cout<<"mathc"<<std::endl;
    int totc=0;
    totc=SGM();
    switch(choice_)
    {
        case 1:totc+=match_hk();
                 break;
        case 2:totc+=match_hk();
                 break;
        case 3:totc+=match_dfs();
                 break;
        case 4:totc+=match_dfs();
                 break;
        default:totc+=match_hk();
    }
#ifdef ISORROPIA_MATCHING_STATS
    std::cout<<totc<<std::endl;
#endif
    matched_ = totc;
    /*std::cout<<"MateU: ";
    for(int i=0; i<5;i++)
        std::cout<<mateU_[i]<<",";
     std::cout<<std::endl;

    for(int i=0;i<5;i++)
        mateV_[mateU_[i]]=i;*/
    //std::cout<<"MateV: ";

    //for(int i=0; i<5;i++)
        //std::cout<<mateV_[i]<<",";
     //std::cout<<std::endl;
    return 0;
}
}// Epetra namespace
#endif
}// Isorropia namespace

