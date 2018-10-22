#ifndef SHYLUBASKER_SOLVE_RHS_KOKKOS_HPP
#define SHYLUBASKER_SOLVE_RHS_KOKKOS_HPP
//This file contains the kokkos functors needed for shylubasker_solve_rhs.hpp

//#define BASKER_DEBUG_SOLVE_RHS_INIT
//#define BASKER_DEBUG_SOLVE_RHS_LOWER
//#define BASKER_DEBUG_SOLVE_RHS_UPPER

#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_types.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

#include <assert.h>

namespace Basker
{
  //Local Lower Sparse Tri-Solve 
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int local_lower_tri_solve
  (
   Int     kid,
   BASKER_MATRIX   &L,
   INT_1DARRAY     &pinv,
   const ENTRY_2DARRAY &x,
   const ENTRY_2DARRAY &b
  )
  {
    //over all rhs
    //for(Int i=0; i < b.size(); i++)
    for(Int i = 0; i < 1; i++)
    {
      //over all columns
      for(Int k = 0; k < L.ncol; k++)
      {

        if(L.val[L.col_ptr[k]] == (Entry) 0)
        {
          printf("kid: %d L found zero pivot in column: %d \n",
              kid, k+L.scol);
          return -1;
        }

        #ifdef BASKER_DEBUG_SOLVE
        printf("kid b(1) = %f \n", b[i][1]); 
        printf("kid: %d %d x(%d) = %f b(%d) = %f\n",
            kid, i , k+L.scol,  x[i][k+L.scol], k+L.scol, b[i][k+L.scol]);
        #endif

        x[i][k+L.scol] = b[i][k+L.scol] / L.val[L.col_ptr[k]];

        //Over each nnz in column
        for(Int j = L.col_ptr[k]+1; j < L.col_ptr[k+1]; j++)
        {
          Int loc = pinv[L.row_idx[j]];

          #ifdef BASKER_DEBUG_SOLVE_RHS_LOWER
          printf("kid %d L update b: %d %d %f  value: %f %f \n",
              kid, i,loc, b[i][loc], L.val[j], x[i][k+L.scol]);
          #endif

          //Add Atomic: come back... only do atomics for shared vals
          Kokkos::atomic_fetch_sub(&(b[i][loc]),
              L.val[j]*x[i][k+L.scol]);

        }//end overall nnz in column
      }//end over all columns
    }//end over all rhs

    return 0;
  }//end local_lower_tri_solve
                              
   //Local Upper Sparse Tri-Solve
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int local_upper_tri_solve
  (
   Int kid,
   BASKER_MATRIX U,
   INT_1DARRAY pinv,
   const ENTRY_2DARRAY x,
   const ENTRY_2DARRAY b
  )
  {
    //over all rhs
    //for(Int i=0; i < b.size(); i++)
    for(Int i = 0; i < 1; i++)
    {
      //over all columns
      for(Int k = U.ncol-1; k > 0; k--)
      {

        if(U.val[U.col_ptr[k+1]-1] == (Entry) 0)
        {
          printf("kid: %d U found zero pivot in column: %d \n",
              kid, k+U.scol);
        }

        #ifdef BASKER_DEBUG_SOLVE_RHS_UPPER
        printf("kid: %d x(%d) = %f b(%d) = %f U.val: %f \n",
            kid, k+U.scol, x[i][k+U.scol], 
            k+U.scol, b[i][k+U.scol], 
            U.val[U.col_ptr[k+1]-1]);
        #endif

        x[i][k+U.scol] = b[i][k+U.scol] / U.val[U.col_ptr[k+1]-1];

        //Over each nnz in column
        for(Int j = U.col_ptr[k+1]-2; j >= U.col_ptr[k]; j--)
        {
          Int loc = U.row_idx[j];
          //No need for atomic, only over own blk
          b[i][loc] = b[i][loc] - U.val[j]*x[i][k+U.scol];
        }//end overall nnz in column

      }//end over all columns

      x[i][U.scol] = b[i][U.scol]/U.val[0];

    }//end over all rhs

    return 0;
  }//end local_upper_tri_solve()

  //Local SpMV Ax = b (will add option to use infor later)
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int local_spmv_update
  (
   Int kid,
   BASKER_MATRIX   &A,
   INT_1DARRAY     &pinv, 
   const ENTRY_2DARRAY   &inX,
   const ENTRY_2DARRAY   &outB,
   int *info
  )
  {
    if(outB.size() != inX.size())
    {return -1;}

    //over each vector (can be moved to inside for better performance
    //for(Int i=0; i < inX.size(); i++)
    for(Int i = 0; i < 1; i++)
    {
      //over each column
      for(Int k = 0; k < A.ncol; k++)
      {
        //over each nnz in a column
        for(Int j=A.col_ptr[k]; j<A.col_ptr[k+1]; j++)
        {
          #ifdef BASKER_DEBUG_SOLVE_RHS_UPPER
          printf("kid: %d, b(%d) = %f  x(%d) = %f Aval: %f \n",
              kid, A.row_idx[j], outB[i][A.row_idx[j]], 
              A.scol+k, inX[i][A.scol+k],  A.val[j]);
          #endif

          outB[i][A.row_idx[j]] -= A.val[j]*inX[i][A.scol+k];
        }//end over all nnz in a column
      }//end over each column
    }//end over each vector

    return 0;
  }//end local_spmv

  //Functor for lower tri-solve
  template <class Int, class Entry, class Exe_Space>
  struct lower_solve_funct
  {
    //Kokkos Typedefs
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                                             execution_space ;
    typedef  Kokkos::TeamPolicy<Exe_Space>                        team_policy;
    typedef  typename team_policy::member_type                     team_member;
    #endif
    //Basker Typedefs
      
    //Local Reference Variables
    MATRIX_2DARRAY L, U;
    ENTRY_2DARRAY X; //workspace
    ENTRY_2DARRAY x; //Ax = PLUx = b
    ENTRY_2DARRAY b;
    INT_2DARRAY  S;
    INT_2DARRAY gperm;
    Int nlvls;

    lower_solve_funct() {} //empty constructor
    
    lower_solve_funct(MATRIX_2DARRAY _L,
                      INT_2DARRAY _gperm,
                      ENTRY_2DARRAY _x,
                      ENTRY_2DARRAY _b,
                      ENTRY_2DARRAY _X,
                      INT_2DARRAY _S,
                      Int _nlvls)
    {
      L = _L;
      gperm = _gperm;
      x = _x;
      b = _b;
      X = _X;
      S = _S;
      nlvls = _nlvls;
    }//end lower_solve_funct()
                      
    //KOKKOS_INLINE_FUNCTION
    BASKER_INLINE
    void operator()(const team_member &thread) const
    {

      //Kokkos stuff
      Int tid = (Int) Exe_Space::hardware_thread_id();
      Int kid = (Int) (thread.league_rank() * thread.team_size()) + thread.team_rank();
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #ifdef BASKER_DEBUG_SOLVE_RHS_LOWER
      printf("\n\n league_rank: %i team_rank:  %i  league_size: %i  team_size: %i \n OMPRank: %d  KRank: %d  Team_Leader: %d  \n\n", 
             thread.league_rank(), thread.team_rank(),
             thread.league_size(), thread.team_size(), 
             tid, kid, team_leader);
      thread.team_barrier();
      #endif

      Int L_col = 0;
      Int L_row = 0;

      for(Int l=0; l <= nlvls; l++)
      {
        if(kid%((Int)pow(2,l))==0)
        {
          L_col = S[l][kid];
          L_row = 0;

          #ifdef BASKER_DEBUG_SOLVE_RHS_LOWER
          printf("\n\n----------------kid= %d l = %d------------------\n\n",
              kid, l);
          printf("kid: %d lvl: %d mod: %d L_col: %d \n", 
              kid, l, kid%(l+1),L_col);
          #endif

          BASKER_MATRIX &LL = L[L_col][L_row];
          local_lower_tri_solve<Int, Entry, Exe_Space>
            (kid, LL, gperm[0], x, b);
        }
      }//over all lvls

      thread.team_barrier();

    } //end operator()

  };//end lower_solve_funct
  

  //Functor for upper tri-solve
  template <class Int, class Entry, class Exe_Space>
  struct upper_solve_funct
  {
    //Kokkos Typedefs
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                                             execution_space ;
    typedef  Kokkos::TeamPolicy<Exe_Space>                        team_policy;
    typedef  typename team_policy::member_type                     team_member;
    #endif
    
    //Local Reference Variables
    MATRIX_2DARRAY L, U;
    ENTRY_2DARRAY X; //workspace
    ENTRY_2DARRAY x; //Ax = PLUx = b
    ENTRY_2DARRAY b;
    INT_2DARRAY    gperm;
    INT_2DARRAY   S;
    Int nlvls;

    upper_solve_funct() {} //empty constructor
    
    upper_solve_funct(MATRIX_2DARRAY _U,
                      INT_2DARRAY  _gperm,
                      ENTRY_2DARRAY _x,
                      ENTRY_2DARRAY _b,
                      ENTRY_2DARRAY _X,
                      INT_2DARRAY   _S, 
                      Int _nlvls)
    {
      U =     _U;
      gperm = _gperm;
      x =     _x;
      b =     _b;
      X =     _X;
      S =     _S;
      nlvls = _nlvls;
    }//end upper_solve_funct()
                      
    //KOKKOS_INLINE_FUNCTION
    BASKER_INLINE
    void operator()(const team_member &thread) const
    {
      //Kokkos stuff
      Int tid = (Int) Exe_Space::hardware_thread_id();
      Int kid = (Int) (thread.league_rank() * thread.team_size()) + thread.team_rank();
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
  

      #ifdef BASKER_DEBUG_SOLVER_RHS_UPPER
      printf("\n\n league_rank: %i team_rank:  %i  league_size: %i  team_size: %i \n OMPRank: %d  KRank: %d  Team_Leader: %d  \n\n ", 
             thread.league_rank(), thread.team_rank(),
             thread.league_size(), thread.team_size(), 
             tid, kid, team_leader);
      thread.team_barrier();
      #endif

      Int U_col = 0;
      Int U_row = 0;

      for(Int l=nlvls; l >= 0; l--)
      {
        U_col = S[l][kid];
        if(kid%((Int)pow(2,l)) == 0)
        {

          U_row = U[U_col].size()-1;

          #ifdef BASKER_DEBUG_SOLVER_RHS_UPPER
          printf("kid: %d lvl: %d mod: %d U_col: %d U_row: %d  \n", 
              kid, l, kid%(l+1), U_col, U_row);
          #endif

          BASKER_MATRIX &LU = U[U_col][U_row];
          local_upper_tri_solve<Int, Entry, Exe_Space>
            (kid, LU, gperm[0], x, b);
        }
        //could remove if iterate over on outside 
        //However that might cost more because of start-up
        thread.team_barrier(); 

        //Spmv update
        {
          //over each sublevel
          for(Int sl=l-1 ; sl >=0; sl--)
          { 
            if(kid%((Int)pow(2,sl)) == 0)
            {
              U_row = S[sl][kid]%(U[U_col].size());

              #ifdef BASKER_DEBUG_SOLVE_RHS_UPPER
              printf("SPMV.  kid %d sl: %d mod: %d U_col: %d U_row:  %d \n",
                  kid, sl, kid%(sl+1), U_col, U_row);
              #endif

              BASKER_MATRIX &LU = U[U_col][U_row];

              //may use spmv_info in future 
              int spmv_info = 0;
              local_spmv_update<Int,Entry,Exe_Space>
                (kid, LU, gperm[0], x, b, &spmv_info);   
            }

            thread.team_barrier();
            //Really don't need a thread barrier here
            //thread.team_barrier();
          }//end over each sublevel
        }//end spmv
      }//over all lvls

    } //end operator()

  };//end upper_solve_funct


  //Functor for init-solve
  template <class Int, class Entry, class Exe_Space>
  struct init_solve_funct
  {
    //Kokkos Typedefs
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                                             execution_space ;
    typedef  Kokkos::TeamPolicy<Exe_Space>                        team_policy;
    typedef  typename team_policy::member_type                    team_member;
    #endif
    //Basker Typedefs

    //Local Reference Variables
    MATRIX_2DARRAY U;

    Entry *x;
    Entry *b;
    ENTRY_2DARRAY sol;
    ENTRY_2DARRAY rhs;
    INT_2DARRAY S;
    Int nlvls;
    INT_2DARRAY gperm;
    
    init_solve_funct() {} //empty constructor
    
    init_solve_funct(Entry *_b,
                     Entry *_x,
                     ENTRY_2DARRAY _rhs,
                     ENTRY_2DARRAY _sol,
                     MATRIX_2DARRAY  _U,
                     INT_2DARRAY   _gperm,
                     INT_2DARRAY   _S,
                     Int           _nlvls)
    {
      x     = _x;
      b     = _b;
      rhs   = _rhs;
      sol   = _sol;
      U     = _U;
      gperm = _gperm;
      S     = _S;
      nlvls = _nlvls;
    }//end init_solve_funct
                      
    //KOKKOS_INLINE_FUNCTION
    BASKER_INLINE
    void operator()(const team_member &thread) const
    {
            //Kokkos stuff
      Int tid = (Int) Exe_Space::hardware_thread_id();
      Int kid = (Int) (thread.league_rank() * thread.team_size()) + thread.team_rank();
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #ifdef BASKER_DEBUG_SOLVE_RHS_INIT
      printf("\n\n league_rank: %i team_rank:  %i  league_size: %i  team_size: %i \n OMPRank: %d  KRank: %d  Team_Leader: %d  \n\n ", 
             thread.league_rank(), thread.team_rank(),
             thread.league_size(), thread.team_size(), 
             tid, kid, team_leader);
      thread.team_barrier();
      #endif
      
      //More column to inner loop!!!
      //Over each column of rhs/b
      //for(Int k = 0; k < rhs.size(); k++)
      for(Int k = 0; k < 1; k++)
      {
        //over every lvl
        for(Int l = 0; l <= nlvls; l++)
        {
          //if((l==0) || (kid%l == 0))
          if(kid%(l+1) == 0)
          {
            Int U_col = S[l][kid];
            Int U_row = U[U_col].size()-1;

            Int scol = U[U_col][U_row].scol;
            Int ecol = U[U_col][U_row].ncol + scol;

            #ifdef BASKER_DEBUG_SOLVE_RHS_INIT
            printf("kid: %d lvl: %d U_col: %d U_row: %d scol: %d ecol: %d \n",
                kid, l, U_col, U_row, scol, ecol);
            #endif

            for(Int j = scol; j < ecol; j++)
            {
              Int pj = gperm[0][j];
              //Int offset = k*(rhs.size(k))+j;
              Int offset = j; //Note:come back and fix

              rhs[k][pj] = x[offset]; //Ax=b flip for reuse
              sol[k][pj] = b[offset];
            } //over nnz in my section
          }//if correct KID
        }//end over lvl (l)
      }//over each column(k) 

    }//end operator();

  };//end init_solve_funct

} //end namespace Basker
#endif //end ifndef basker_solve_rhs_hpp
