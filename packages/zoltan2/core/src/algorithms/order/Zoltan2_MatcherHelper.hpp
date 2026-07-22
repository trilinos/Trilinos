// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_MatcherHelper_hpp_
#define _ZOLTAN2_MatcherHelper_hpp_

//#define _ZOLTAN2_MatcherHelper_hpp_


#ifdef ZOLTAN2ND_HAVE_OMP
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <queue>
#include <cassert>

// PPF Matching code copied from Isorropia.  Siva has reservations about algorithm
// robustness.  It uses recursion, which can cause problems on the stack.
// Probably should rewrite at least some of these algorithms, eventually.


namespace Zoltan2 {

/** @ingroup matching_grp
    An implementation of the Matcher interface that operates on Epetra matrices
    and Graphs.
    
    matching algorithms provides an interface to solve the Bipartite Matching
    problem.

*/
////////////////////////////////////////////////////////////////////////////////


template <typename LO>
class Matcher 
{
private:

    LO *mCRS_rowPtrs;
    LO *mCRS_cols;
  
    // Number of vertices in u set (num row vertices)
    LO numU;

    // Number of vertices in v set (num col vertices)
    LO numV;

    // For each u vertex, the matching v vertex
    std::vector<LO> matchedVertexForU;

    // For each v vertex, the matching u vertex
    std::vector<LO> matchedVertexForV;

    std::vector<LO> queue;
    LO* LV_;
    LO* LU_;
    LO* unmatchedU_;
    LO *parent_;
    LO *lookahead_;
    bool finish_;
    LO numE,avgDegU_,k_star_,icm_,BFSInd_,numThread_,Qst_,Qend_,numMatched;

    LO *visitedV_;

    void delete_matched_v();
    LO SGM();  
    LO match_dfs();
    LO match_hk();
    LO construct_layered_graph();
    LO find_set_del_M();
    LO recursive_path_finder(LO, LO);
    LO dfs_path_finder(LO);
    LO dfs_augment();
    LO augment_matching(LO);
    LO DW_phase();

public:
    /** @ingroup matching_grp
    Constructor
    \param[in] row pointer for CRS matrix for of bipartite graph
    \param[in] cols for CRS matrix for of bipartite graph
    \param[in] Number of vertices in u set (num row vertices)
    \param[in] Number of vertices in v set (num col vertices)
    */
    Matcher(LO *_rowPtr, LO *_cols, LO _numU, LO _numV, LO _numE);



      
    /** @ingroup matching_grp
    Destructor
    */
     virtual ~Matcher();

    /* @ingroup matching_grp
    Returns the number of matched vertices from Maximum Cardinality
     Matching set
    */
    LO getNumberOfMatchedVertices()
    {
      LO count=0;
      for(unsigned int i=0;i<matchedVertexForU.size(); i++)
      {
        if(matchedVertexForU[i]!=-1)
	{
          count++;
	}
      }
      return count;
    }



    const std::vector<LO> &getVertexUMatches() {return matchedVertexForU;};
    const std::vector<LO> &getVertexVMatches() {return matchedVertexForV;};


    // Perhaps should be moved outside of class eventually or reworked to use internal variables
    void getVCfromMatching(const std::vector<LO> &bigraphCRSRowPtr,
		       std::vector<LO> &bigraphCRSCols,
                       const std::vector<LO> &vertUMatches,
		       const std::vector<LO> &vertVMatches,
		       const std::vector<LO> &bigraphVMapU,
		       const std::vector<LO> &bigraphVMapV,
		       std::vector<LO> &VC);


    /** @ingroup matching_grp
    Computes the maximum cardinality matching.
    */
    LO match();
};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
Matcher<LO>::Matcher(LO *_rowPtr, LO *_cols, LO _numU, LO _numV, LO _numE)
  :mCRS_rowPtrs(_rowPtr),mCRS_cols(_cols),numU(_numU),numV(_numV),numE(_numE)
{
    finish_=false;
    avgDegU_=numE/numU+1;
    matchedVertexForU.resize(numU);
    matchedVertexForV.resize(numV);

    lookahead_=new LO[numU];
    unmatchedU_=new LO[numU];

    for(LO i=0;i<numU;i++)
    {
        matchedVertexForU[i]=-1;

        lookahead_[i]=mCRS_rowPtrs[i];
        unmatchedU_[i]=i;
    }

    visitedV_=new LO[numV];

    parent_=new LO[numV];

    for(LO i=0;i<numV;i++)
    {
        visitedV_[i]=0;

        matchedVertexForV[i]=-1;
        parent_[i]=-1;
    }

    numThread_=1;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
Matcher<LO>::~Matcher()
{
    delete [] lookahead_;
    delete [] unmatchedU_;

    if (parent_)
    { 
      delete [] parent_; parent_ = NULL;
    }

    if(visitedV_)
      {
	delete [] visitedV_;
      }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This function increases the matching size by one with the help of a
// augmenting path. The input is an integer which is the id of the last
// columns vertex of the augmenting path. We trace the whole path by
// backtraking using parent array. while tracing back we flip the unmatched
// edges to matched edges.
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::augment_matching(LO tv)
{
    LO u,v,t,lnt=1;
    v=tv;
    while(true)
    {
        u=parent_[v];
        t=matchedVertexForU[u];
        matchedVertexForV[v]=u;
        matchedVertexForU[u]=v;
        if(t==-1)
            break;
        lnt+=2;
        v=t;
    }

    return lnt;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This function is almost similar to the previous function which is to find
// a vertex disjoint path. It is used for the algorithm DFS, PPF and in HKDW. The
// difference is that this function operates on the original graph not on the
// layerd subgraph. It also does the incorporates the lookahead mechanism and
// scanning adjacency list alternately from backward and from forward in
// alternate iterations for PPF and HKDW.
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::dfs_path_finder(LO u)
{
    LO i,ind=-1,res=0;

        for(i=lookahead_[u];i<mCRS_rowPtrs[u+1];i++) // the lookahead scheme
        {
            ind=mCRS_cols[i];
            assert(ind>=0 && ind<numV);
            if(matchedVertexForV[ind]==-1)
            {
	      LO lock2=0;
	      if (visitedV_[ind]==0)
                {
                  visitedV_[ind]=1;
                  lock2=1;
                }

                if(lock2>0)
                {
                    parent_[ind]=u;
                    lookahead_[u]=i+1;
                    return ind;
                }
            }
        }


        if(icm_%2==1) // odd number iteration so scan the adj list forward dir
        {
            for(i=mCRS_rowPtrs[u];i<mCRS_rowPtrs[u+1];i++)
            {
                ind=mCRS_cols[i];
                assert(ind>=0 && ind<numV);

	      LO lock2=0;
	      if (visitedV_[ind]==0)
                {
                  visitedV_[ind]=1;
                  lock2=1;
                }

                if(lock2>0)
                {
                    parent_[ind]=u;
                    res=dfs_path_finder(matchedVertexForV[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }
        else // even number iteration so scan from backward
        {
            for(i=mCRS_rowPtrs[u+1]-1;i>=mCRS_rowPtrs[u];i--)
            {
                ind=mCRS_cols[i];
                assert(ind>=0 && ind<numV);


	      LO lock2=0;
	      if (visitedV_[ind]==0)
                {
                  visitedV_[ind]=1;
                  lock2=1;
                }

                if(lock2>0)
                {
                    parent_[ind]=u;
                    res=dfs_path_finder(matchedVertexForV[ind]);
                    if(res!=-1)
                        return res;
                }
            }
        }


    return -1;
}
////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// // This function starts the BFS phase for the HK and HKDW. It starts from the
// // layer 0 vertices and tries to find a vertex disjoint path from each of
// // these vertices.
// ////////////////////////////////////////////////////////////////////////////////
// int Matcher::find_set_del_M()
// {

//     int i,j,count=0;
//     delete_matched_v();

// #ifdef ZOLTAN2ND_HAVE_OMP
//     #pragma omp parallel for private(j)
// #endif
//     for(i=0;i<BFSInd_;i++)
//     {

//         j=recursive_path_finder(0,queue[i]);
//         if(j!=-1)
//         {
//             augment_matching(j);

// 	    //MMW: Not 100% this is necessary, let this in just in case bug in original code
// #ifdef ZOLTAN2ND_HAVE_OMP
//             #pragma omp atomic
// #endif
//             count++;



//         }

//     }
//     return count;
// }
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// // This is the additional Duff and Wiberg phase for the HKDW. This function
// // does nothing but first unset the locks and then runs PPF from the
// // remaining unmatched row vertices after the BFS phase of HK.
// ////////////////////////////////////////////////////////////////////////////////
// int Matcher::DW_phase()
// {
//     int i,count=0;

// #ifdef ZOLTAN2ND_HAVE_OMP
//     #pragma omp parallel for
// #endif
//     for(i=0;i<numV;i++)
//     {
// #ifdef ZOLTAN2ND_HAVE_OMP
//         omp_unset_lock(&scannedV_[i]); // unsetting the locks
// #endif
//     }

// #ifdef ZOLTAN2ND_HAVE_OMP
//     #pragma omp parallel for
// #endif
//     for(i=0;i<BFSInd_;i++) // calling the PPF
//     {
//         int u=queue[i];
//         if(matchedVertexForU[u]==-1)
//         {
//             int ind=dfs_path_finder(u);
//             if(ind!=-1)
//             {
//                 augment_matching(ind);

// 	        //MMW: Not 100% this is necessary, let this in just in case bug in original code
// #ifdef ISORROPIA_HAVE_OMP
//                 #pragma omp atomic
// #endif
//                 count++;

//             }
//         }
//     }
//     return count;
// }
// ////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
    //This function is the starter function for PPF and DFS. It unsets the locks
    //the call dfs_path_finder and then call the augment_matching.
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::dfs_augment()
{
    LO i,flag=0,flag1=0,count=0,totc=0,index=numU,cur=0;
    icm_=0;

    while(true)
    {
        icm_++;

        for(i=0;i<numV;i++)
	{
          visitedV_[i]=0;
	}

        cur=0;
        for(i=0;i<numU;i++)
	{
	  if(matchedVertexForU[i]==-1)
	    unmatchedU_[cur++]=i;
	}
        index=cur;
        flag=flag1=count=0;

        for(i=0;i<index;i++)
        {
            flag=1;
            LO u=unmatchedU_[i];
            LO ind=dfs_path_finder(u);
            if(ind!=-1)
            {
                flag1=1;
                augment_matching(ind);
            }
        }

        if(flag==0 || flag1==0)
            break;
        else
        {
            cur=0;
            for(i=0;i<index;i++)
            {
                if(matchedVertexForU[unmatchedU_[i]]==-1)
                    unmatchedU_[cur++]=unmatchedU_[i];
            }
            index=cur;

        }
    }

    return totc;
}
////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
// int Matcher::SGM()
// {
//    int i,j,lock,ind,count=0;
// #ifdef ZOLTAN2ND_HAVE_OMP
//    #pragma omp parallel for private(j,ind,lock)
// #endif
//    for(i=0;i<numU;i++)
//    {
//       for(j=mCRS_rowPtrs[i];j<mCRS_rowPtrs[i+1];j++)
//       {
//          ind=mCRS_cols[j];
// #ifdef ZOLTAN2ND_HAVE_OMP
//          lock=omp_test_lock(&scannedV_[ind]);
// #else
//             // mfh 07 Feb 2013: lock wasn't getting initialized if
//             // ZOLTAN2ND_HAVE_OMP was not defined.  omp_test_lock()
//             // returns nonzero if the thread successfully acquired the
//             // lock.  If there's only one thread, that thread will
//             // always be successful, so the default value of lock
//             // should be nonzero.
//          lock = 1;
// #endif
//          if(lock>0)
//          {
//             matchedVertexForU[i]=ind;
//             matchedVertexForV[ind]=i;
//             break;
// 	 }
//       }
//     }
//     return count;
// }
// ////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::SGM()
{
  LO i,j,ind,count=0;

  for(i=0;i<numU;i++)
    {
      for(j=mCRS_rowPtrs[i];j<mCRS_rowPtrs[i+1];j++)
        {
	  ind=mCRS_cols[j];

	  LO lock2=0;
	  if (visitedV_[ind]==0)
          {
            visitedV_[ind]=1;
            lock2=1;
          }

	  if(lock2>0)
            {
	      matchedVertexForU[i]=ind;
	      matchedVertexForV[ind]=i;
	      break;
            }
        }
    }
  return count;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Forking function for DFS based algorithm
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::match_dfs()
{
    LO totc=0;
    icm_=0;

    totc=dfs_augment();

    return totc;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
LO Matcher<LO>::match()
{
  // User interface function for the matching..                                                                                                             

  LO totc=0;
  totc=SGM();
  totc+=match_dfs();

  numMatched = totc;

  return 0;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Calculate vertex cover (which is vertex separator) from matching
// VC = (U-L) union (B intersection L)
// 
// Let X be the exposed nodes in U (those without a match)
//
// E':
// Matched edges should be directed VtoU -- just use vertVMatches array
// All other edges should be directed UtoV -- use modified biGraphCRScols
//
// L is the set of vertices in (U union V, E') that can be reached from X
//
//
// Perhaps I should use internal class members in this function
// However, I decided not to do this for now since I'm changing some of these
// member variables in this function
////////////////////////////////////////////////////////////////////////////////
template <typename LO>
void Matcher<LO>::getVCfromMatching(const std::vector<LO> &bigraphCRSRowPtr,
		       std::vector<LO> &bigraphCRSCols,
                       const std::vector<LO> &vertSMatches,
		       const std::vector<LO> &vertTMatches,
		       const std::vector<LO> &bigraphVMapU,
		       const std::vector<LO> &bigraphVMapV,
		       std::vector<LO> &VC)
{
   LO numS = vertSMatches.size();
   LO numT = vertTMatches.size();

   //////////////////////////////////////////////////////////////////////
   // Build X
   //////////////////////////////////////////////////////////////////////
   std::set<LO> X;
   for(LO i=0;i<numS; i++)
   {
     if(vertSMatches[i]==-1)
     {
       X.insert(i);
     }
   }
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // Non matched edges should be directed UtoV
   // Removed matches edges from bipartite graph (set to -1)
   //   -- perhaps replace this with something more efficient/compact
   //////////////////////////////////////////////////////////////////////
   for (LO uID=0;uID<numS;uID++)
   {
     for (LO eIdx=bigraphCRSRowPtr[uID];eIdx<bigraphCRSRowPtr[uID+1];eIdx++)
     {
       LO vID = bigraphCRSCols[eIdx];

       if (vertSMatches[uID]==vID)
       {
         bigraphCRSCols[eIdx]=-1;
       }
     }

   }
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // Calculate L - set of vertices in (U union V, E') that can be reached from X
   //////////////////////////////////////////////////////////////////////
   std::set<LO> L;

   std::queue<LO> vqueue;

   std::copy(X.begin(), X.end(), std::inserter( L, L.begin() ) );



   typename std::set<LO>::const_iterator iter;
   for(iter = X.begin(); iter != X.end(); ++iter)
   {
     L.insert(bigraphVMapU[*iter]);   // L contains all vertices in X

     vqueue.push(*iter); // copy on to vqueue
   }

   LO iteration=0;
   while(vqueue.size()!=0)
   {
     //Number of active vertices on this side of bipartite graph
     LO nstarts=vqueue.size();

     for (LO i=0; i<nstarts; i++)
     {
       LO start = vqueue.front();
       vqueue.pop();

       //Traverse from U to V
       if(iteration%2==0)
       {
         //Traverse edges starting from vertex "start"
         for (LO eIdx=bigraphCRSRowPtr[start];eIdx<bigraphCRSRowPtr[start+1];eIdx++)
	 {
	   LO newV = bigraphCRSCols[eIdx];

           //Edge is in correct direction (U to V) => newV is valid
           if (newV!=-1)
	   {
             // If this vertex has not been reached
             if(L.find(bigraphVMapV[newV])==L.end())
	     {
               L.insert(bigraphVMapV[newV]);
               vqueue.push(newV);
	     }
	   }
	 }

       }

       //Traverse from V to U
       else
       {
         LO newU = vertTMatches[start];

         // If this vertex has not been reached
         if(L.find(bigraphVMapU[newU])==L.end())
	 {
           L.insert(bigraphVMapU[newU]);
           vqueue.push(newU);
	 }
       }
     } // for

     iteration++;
   } //while
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // Calc VC = (U-L) union (B intersection L)
   //////////////////////////////////////////////////////////////////////
   for(LO uID=0;uID<numS;uID++)
   {
     // if vertex not in L, it is in VC
     if(L.find(bigraphVMapU[uID])==L.end())
     {
       VC.push_back(bigraphVMapU[uID]);
     }
   }

   for(LO vID=0;vID<numT;vID++)
   {
     // if vertex is in L, it is in VC
     if(L.find(bigraphVMapV[vID])!=L.end())
     {
       VC.push_back(bigraphVMapV[vID]);
     }
   }
   //////////////////////////////////////////////////////////////////////



}
////////////////////////////////////////////////////////////////////////////////













} //Zoltan2 namespace
#endif //ifdef
