#ifndef SHYLUBASKER_NFACTOR_KOKKOS_HPP
#define SHYLUBASKER_NFACTOR_KOKKOS_HPP
//This file contains the kokkos functors needed for basker_nfactor.hpp

#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

//#include <assert.h>

//#define BASKER_DEBUG_LOCAL_REACH 1
//#define BASKER_DEBUG_REC_FACTOR 1
//#define BASKER_DEBUG_BACKSOLVE 1
//#define BASKER_DEBUG_ATOMIC_BACKSOLVE 1
//#define BASKER_DEBUG_KOKKOS_REC_FACTOR 1
#define BASKER_DEBUG_MOVEA 1
#define BASKER_DEBUG_COLUMN_TRISOLVE 1
#define BASKER_DEBUG_COLUMN_FACTOR 1
#define BASKER_DEBUG_EXTEND 1

//#define BASKER_TIME_NFACTOR

//#define BASKER_TIME_DETAIL
//#define BASKER_TIME_EXTEND

/*Working Notes:
Changes:
1. fix pinv - move_A
2. Move into class structure
3. Dense reduce for each with token.
4. Sub factors for L-Sep
*/

namespace Basker
{
  //Will be moved at some point, add parameter kid
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int init_workspace(const Int b,
		     BASKER_MATRIX &L,
		     BASKER_MATRIX &U,
		     INT_1DARRAY   WS,
		     ENTRY_1DARRAY  X,
                     Int ws_size)
  {

    printf("init kid: %d ws_size: %d \n",
           b, ws_size);
    for(Int i=0; i < WS.size(); i++)
      { WS[i] = 0;}
    
    for(Int i = 0; i < X.size(); i++)
      {X[i] = 0; }

    /*Change later*/
    L.fill();
    U.fill();

    L.init_perm();
    return 0;
  }//end init_workspace

  //Will be moved at some point, add parameter kid
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int init_workspace(const Int b,
		     BASKER_MATRIX &L,
		     INT_1DARRAY &WS,
		     ENTRY_1DARRAY &X)
  {
    for(Int i=0; i < WS.size(); i++)
      { WS[i] = 0;}
    
    for(Int i = 0; i < X.size(); i++)
      {X[i] = 0; }

    /*Change later*/
    for(Int i = 0 ; i < L.ncol+1; i++)
      {
        L.col_ptr[i] = 0;
      }
    for(Int i = 0; i < L.nnz; i++)
      {
        L.row_idx[i] = 0; 
        L.val[i] = 0;
      }
    
    L.init_perm();
    return 0;
  }//end init_workspace

  //Will be moved at some point, add parameter kid
  //DEFUNCT
  
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int init_workspace(const Int b,
		     INT_1DARRAY &WS,
		     ENTRY_1DARRAY &X)
  {
    for(Int i=0; i < (WS.size()); i++)
      { WS[i] = 0;}
    
    for(Int i = 0; i < (X.size()); i++)
        {X[i] = 0; }
      return 0;
  }//end init_workspace

  
  
  //Note: Needs to be updated to MatrixView ?, add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int local_reach
  (
   Int b, //blk
   Int j, //starting point
   BASKER_MATRIX &L,
   INT_1DARRAY   ws,
   INT_1DARRAY   gperm,
   Int *top,
   Int bcol,
   Int brow, 
   Int ws_size
   )
  {
       //Int color_offset = 0;
    Int pattern_offset = ws_size;
    Int stack_offset = pattern_offset + ws_size;
    Int store_offset = stack_offset + ws_size;

    Int i, t, head, i1;
    Int start, end, done;
    
    start = -1;

    head = 0;
    ws[stack_offset + head] = j;
    
    while(head < L.max_idx)
      { 
        #ifdef BASKER_DEBUG_LOCAL_REACH
        printf("stack_offset: %d head: %d \n", stack_offset , head);
        ASSERT(head > -1);
        #endif

        j = ws[stack_offset + head];
        t = gperm[j];
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("----------DFS: %d %d -------------\n", j, t);
        #endif

        //t = L.lpinv[j];
        //t = gperm[j];
        //If not colored
        if(ws[j] == 0)
	  {
	    ws[j] = 1;  
            //if not permuted
            //if(t != L.max_idx )
            if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
	      {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("reach.... j: %d t:%d L.scol %d \n", j, t, L.scol);
                #endif
                start = L.col_ptr[t-L.scol];
	      }
            else
              {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("L.scol: %d  L.ncol: %d t: %d \n", L.scol, L.ncol,t);
                #endif
              }
	  }
	else
	  {
            start = ws[store_offset+j];
	  }
	done = 1;
	
        //if  permuted
        //if(t != L.max_idx)
        if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
	  {
	    end = L.col_ptr[t+1-L.scol];
	    for(i1 = start; i1 < end; i1++)
	      {
                i = L.row_idx[i1];
		if(ws[i] == 0)
		  {
                    head++;
                    ws[stack_offset + head] = i;
                    ws[store_offset + j] = i1+1;
		    done = 0;
		    break;
		  }
	      }
	  }
	if(done)
	  {
	    (*top)--;

            #ifdef BASKER_DEBUG_LOCAL_REACH
            printf("pattern_offset %d top %d \n", pattern_offset, *top);
            #endif

            ws[pattern_offset + *top] = j;
            ws[j] = 2;
	    if(head == 0)
	      { head = L.max_idx;}
	    else
	      {head--;}
	  }
      }//end while 
    return 0;
  }//end local_reach
  
  //Combine with atomic back solve, may need to update to MatrixView
  //try to get vectorized
  //Update to local_back_solve ?
  //Add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int back_solve(Int kid, Int k, Int top,
		 Int brow, 
		 Int xnnz, 
		 Int ws_size,
		 INT_1DARRAY ws,
		 ENTRY_1DARRAY X,
		 BASKER_MATRIX L,
		 INT_1DARRAY   gperm, 
		 const BASKER_STATS &stat
		 )		
  {
    //Int color_offset = 0;
    Int pattern_offset = ws_size;
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;
    for(pp = 0; pp < xnnz; pp++)
      {
        j = ws[pattern_offset+top1];
        top1++;
        ws[j] = 0; //zeros out color
    
        //t = L.lpinv[j];
        t = gperm[j];
        //if(t != L.max_idx)
          if(t != L.max_idx && (t >= L.scol) && (t < (L.scol+L.ncol)))
          {
           
            xj = X[j];
           
            #ifdef BASKER_DEBUG_BACKSOLVE
            printf("Updating column: %d  with %f \n", t, xj);
            #endif

            //Get rid of these temp variables
            Int local_offset = L.scol;
            p2 = L.col_ptr[t+1-local_offset];
            p = L.col_ptr[t-local_offset]+1;
            
            #ifdef BASKER_TIME_DETAIL
            //stat.nfactor_domain_flop[kid] += (p2-p);
            #endif

            for( ; p < p2; p++)
              {
                X[L.row_idx[p]] -= L.val[p] *xj;
              }//end for() over each nnz in the column
            
          }//end if() not permuted
      }//end for() over all nnz in LHS
    return 0;
  }//end back_solve

  //Move to basker_util.hpp
  //Make it for each Matrix and not LU combined
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void print_factor(BASKER_MATRIX &L,
		    BASKER_MATRIX &U)
  {
    printf("L.nnz: %d  L.ncol: %d \n", L.nnz, L.ncol);
    Int k;
    Int Lnnz = L.col_ptr[L.ncol];
    for(k = 0; k < Lnnz; k++)
      {
        //cout << "L[" << k << "] = " << L.val[k];
        printf("L[%d]=%f " ,k , L.val[k]); 
      }
    //cout << endl << endl;
    printf("\n");
    for(k = 0; k < Lnnz; k++)
      {
        //cout << "Li[" << k << "] = " << L.row_idx[k];
        printf("Li[%d]=%d ", k , L.row_idx[k]);
      }
    //cout << endl << endl;
    printf("\n");
    for(k = 0; k < L.ncol; k++)
      {
        //cout << "Lp[" << k  << "] = " << L.col_ptr[k];
        printf("Lp[%d]=%d ", k, L.col_ptr[k]);
      }
    //cout << endl << endl;
    printf("\n");


    printf("U.nnz: %d  U.ncol: %d \n", U.nnz, U.ncol);
    Int Unnz = U.col_ptr[U.ncol];
    for(k = 0; k < Unnz; k++)
      {
        //cout << "U[" << k  << "] = " << U.val[k];
        printf("U[%d]=%f ", k, U.val[k]);
      }
    printf("\n");
    //cout << endl << endl;
    for(k = 0; k < Unnz; k++)
      {
        printf("Ui[%d]=%d ", k, U.row_idx[k]);
        //cout << "Ui[" << k  << "] = " << U.row_idx[k];
      }
    //printf("\n");
    cout << endl << endl;
    for(k =0; k < U.ncol; k++)
      {
        printf("Up[%d] = %d ", k, U.col_ptr[k]);
        //cout << "Up[" << k << "] = " << U.col_ptr[k];
      }
    printf("\n\n");
    //cout << endl << endl;
  }//end print_factor

  //Update to use Matrix view of A (this will make a more uniform interface)
  //Add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int rect_factor(Int b, Int lval, Int uval,
		  Int bcol, Int brow,
		  INT_1DARRAY ws,
		  ENTRY_1DARRAY X,
		  INT_1DARRAY gperm,
		  BASKER_MATRIX A,
		  BASKER_MATRIX L,
		  BASKER_MATRIX U,
		  Int malloc_option,
		  BASKER_STATS const stats,
		  int *einfo)
  {   
    Int kid = b;

    Int i,j,k;
    INT_1DARRAY tptr, color, pattern, stack;
    Int top, top1, maxindex, t; //j1 and j2
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;
   
    Int newsize;
    Entry pivot, value;
    Entry absv, maxv;
    Int ws_size = A.nrow;
    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    Int scol = L.scol;
    Int ecol = L.ecol;
    
    Int color_offset = 0;
    Int pattern_offset = ws_size;
    Int stack_offset = pattern_offset + ws_size;

    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
    top1 = ws_size;
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_REC_FACTOR
    printf("b: %d scol: %d ecol: %d llnzz: %d uunzz: %d \n", 
           b, scol, ecol, L.nnz, U.nnz);
    #endif

    //for each column
    for(k = scol; k < ecol; k++)
        {
          #ifdef BASKER_DEBUG_REC_FACTOR
          printf("\n----------------K=%d--------------\n", k);
	  #endif
	  value = 0.0;
	  pivot = 0.0;
	  maxindex = A.ncol;  
	  lcnt = 0;
	  ucnt = 0;

          #ifdef BASKER_DEBUG_REC_FACTOR
          ASSERT(top == ws_size);
          //ASSERT entry workspace is clean
          for(i = 0 ; i < ws_size; i++){ASSERT(X[i] == 0);}
          //ASSERT int workspace is clean
	  for(i = 0; i <  ws_size; i++){ASSERT(ws[i] == 0 );}
          #endif

          //for each nnz in column
	  for(i = A.col_ptr[k]; i < A.col_ptr[k+1]; i++)
	    {
	      j = A.row_idx[i];
              X[j] = A.val[i];

              #ifdef BASKER_DEBUG_REC_FACTOR
              printf("i: %d row: %d  val: %g  top: %d \n", 
                     i, j ,A.val[i], top);
              printf("Nx in Ak %d %g %d color = %d \n",
                      j, X[j], brow,  
                     ws[color_offset + j] );
              #endif

              #ifdef BASKER_TIME_DETAIL
              Kokkos::Impl::Timer timer_reach;
              #endif

              //Search reach if not yet considered
               if(ws[color_offset + j] == 0)
               {            
                 local_reach<Int, Entry, Exe_Space>
                   (b, j, L, ws,gperm, &top,bcol, brow, 
                    ws_size);
		}

               #ifdef BASKER_TIME_DETAIL
               stats.nfactor_domain_reach_time[kid] += timer_reach.seconds();
               #endif


          }//end for() each nnz in column
	  xnnz = ws_size - top;

          #ifdef BASKER_TIME_DETAIL
          //stats.nfactor_domain_reach_flop[kid] += xnnz;
          #endif


          #ifdef BASKER_DEBUG_REC_FACTOR
          printf("xnnz: %d ws_size: %d top: %d \n", 
                 xnnz, ws_size, top);
          #endif
          

          #ifdef BASKER_TIME_DETAIL
          Kokkos::Impl::Timer timer_solve;
          #endif
          
          back_solve<Int,Entry, Exe_Space>
            ( b, k, top,
                   brow, xnnz, ws_size,
                     ws, X, L, gperm, stats); 
         

          #ifdef BASKER_TIME_DETAIL
          stats.nfactor_domain_solve_time[kid] += timer_solve.seconds();
          #endif

          //printf("b: %d lcnt: %d before \n", b, lcnt);
          //find pivot
          maxv = 0.0;
          for(i = top; i < ws_size; i++)
            {
              j = ws[pattern_offset + i];
              //t = L.lpinv[j];
              t = gperm[j];
              //value = X[j-brow];
              value = X[j];
              absv = abs(value);
              if(t == L.max_idx)
                {
                  lcnt++;
                  if(absv > maxv) 
                    {
                      maxv = absv;
                      pivot = value;
                      maxindex = j;                
                    }
                }
            }//for (i = top; i < ws_size)
          //printf("b: %d lcnt: %d after \n", b, lcnt);
       
          
          ucnt = ws_size - top - lcnt +1;
          if((maxindex == L.max_idx) || (pivot == 0))
            {
              cout << "Error: Matrix is singular" << endl;
              cout << "MaxIndex: " << maxindex << " pivot " 
                   << pivot << endl;
              *einfo = k;
              return 2;
            }          
          //L.lpinv[maxindex] = k;
          gperm[maxindex] = k;
          //printf("TAG1 r  maxindex = %d k= %d \n", maxindex, k);

          #ifdef BASKER_DEBUG_FACTOR
          if(maxindex != k)
            {
              cout << "Permuting Pivot: " << k << " as row " 
                   << maxindex << endl;
            }
          #endif
          

          //Come back to this!!!!
          if(lnnz + lcnt >= llnnz)
            {
              if(malloc_option == 1)
                {//return 1; 
                }
              newsize = lnnz * 1.1 + 2 *A.nrow + 1;
              printf("b: %d Reallocing L oldsize: %d newsize: %d \n",
                     b, llnnz, newsize);
            }
        
          if(unnz+ucnt >= uunnz)
            {
              if(malloc_option ==1)
                {//return 1;
                }
              newsize = uunnz*1.1 + 2*A.nrow+1;
              printf("b: %d Reallocing U oldsize: %d newsize: %d \n",
                     b, uunnz, newsize);
            }

          L.row_idx[lnnz] = maxindex;
          L.val[lnnz] = (Entry) 1.0;
          lnnz++;
     
          Entry lastU = (Entry) 0.0;
          for( i = top; i < ws_size; i++)
            {
              j = ws[pattern_offset + i];
              //t = L.lpinv[j];
              t = gperm[j];
            
              #ifdef BASKER_DEBUG_REC_FACTOR
              printf("j: %d t: %d \n", j, t);
              #endif            

              //if fill-in
              if(X[j] != 0)
                {
                  if(t != L.max_idx)
                    {
                      if(t < k)
                        {
                          //U.row_idx[unnz] = L.lpinv[j];
                          U.row_idx[unnz] = gperm[j];
                          U.val[unnz] = X[j];
                          unnz++;
                        }
                      else
                        {
                          lastU = X[j];
                        }
                    }
                  else if (t == L.max_idx)
                    {
                      L.row_idx[lnnz] = j;
                      L.val[lnnz] = X[j]/pivot;
                      lnnz++;
                    }
                }//end if() not 0
              
              //Note: move x[j] inside of if() not 0....extra ops this way
              #ifdef BASKER_DEBUG_REC_FACTOR
              printf("Zeroing element: %d \n", j);
              #endif
              
              X[j] = 0;
            }//end if(x[j-brow] != 0)

          //Fill in last element of U
          U.row_idx[unnz] = k;
          U.val[unnz] = lastU;
          unnz++;

          xnnz = 0;
          top = ws_size;
          
          L.col_ptr[k-bcol] = cu_ltop;
          L.col_ptr[k+1-bcol] = lnnz;
          cu_ltop = lnnz;
          
          U.col_ptr[k-bcol] = cu_utop;
          U.col_ptr[k+1-bcol] = unnz;
          cu_utop = unnz;
	}//end for() over all columns

    L.nnz = lnnz;
    U.nnz = unnz;

    #ifdef BASKER_DEBUG_REC_FACTOR
    print_factor<Int,Entry,Exe_Space>(L,U);
    #endif
   
    return 0;
  }//end rect_factor

  //Kokkos rect_factor_functor
  //May want to update to use matrix view in place of matrix of A
  template <class Int, class Entry, class Exe_Space>
  struct rect_factor_funct
  {
    //Kokkos Typedefs
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                                             execution_space ;
    typedef  Kokkos::TeamPolicy<Exe_Space>                        team_policy;
    typedef  typename team_policy::member_type                     team_member;
    #endif

    BASKER_MATRIX  A;
    MATRIX_2DARRAY L;
    MATRIX_2DARRAY U;
    INT_2DARRAY    WS;
    ENTRY_2DARRAY  X;
    INT_1DARRAY S; //schedule
    Int lvl;
    
    //Global perm
    INT_2DARRAY gperm;

    //Stats mod
    BASKER_STATS  stats;

    rect_factor_funct()
    {}

    rect_factor_funct(BASKER_MATRIX _A, 
                      MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
                      INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_1DARRAY _S, 
                      Int _lvl)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
    }
    
    rect_factor_funct(BASKER_MATRIX _A, 
                      MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
                      INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_1DARRAY _S, 
                      Int _lvl, INT_2DARRAY _gperm)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
      gperm = _gperm;
    }

    rect_factor_funct(BASKER_MATRIX _A, 
                      MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
                      INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_1DARRAY _S, 
                      Int _lvl, INT_2DARRAY _gperm, BASKER_STATS _stats)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
      gperm = _gperm;
      stats = _stats;
    }

  
    BASKER_INLINE
    void operator()(const team_member  &thread) const
    {
      #ifdef BASKER_TIME_DETAIL
      Kokkos::Impl::Timer timer;
      #endif

      //Need to wrap around Exe_Space
      int tid = Kokkos::OpenMP::hardware_thread_id();
      int i = thread.league_rank();
      Int kid = (Int) (thread.league_rank() * thread.team_size())+ thread.team_rank(); 

      //#ifdef BASKER_DEBUG_KOKKOS_REC_FACTOR
      printf("Hello World: %i %i \\ %i %i \n", 
             thread.league_rank(), thread.team_rank(),
             thread.league_size(), thread.team_size());
      printf("Kokkos: rank %d \n",  tid);
      //#endif

      Int b = S[i];
      int einfo = 0;
   
      //printf("Domain kid: %d working on blk %d  \n", kid, b);
         
      init_workspace<Int, Entry, Exe_Space>(kid, 
                         (L[b])[0], (U[b])[U[b].size()-1], 
                                            (WS[b]),  (X[b]), A.nrow);
   
      #ifdef BASKER_DEBUG_KOKKOS_REC_FACTOR
      printf("Kokkos: Factor: %d working on blk: %d \n", i, b);
      #endif

   
      rect_factor<Int, Entry,Exe_Space>
        (kid, 0, 0, 
         (L[b])[0].scol, (L[b])[0].srow, 
	 WS[b],
	 X[b],
         gperm[0], 
         A,
	 (L[b])[0], 
	 (U[b])[U[b].size()-1], 
          1, //malloc option 
         stats,
	 &einfo); 
   
           
      #ifdef BASKER_TIME_DETAIL
      stats.nfactor_domain_time[i] += timer.seconds();
      #endif


    }//end operator
  }; //end rect_factor_funct

  //helper wrapper to init a column, without permuation
  //add paramter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int populate_row_view(BASKER_MATRIX_VIEW &A,
			Int k, 
			Int prev = 0)
  {
    A.init_offset(k, prev);
    return A.offset;
  }//end populate_row_view


  //helper wrapper to init a column, with permuation
  //add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int populate_row_view(
			BASKER_MATRIX_VIEW &A,
			BASKER_MATRIX &L,
			INT_1DARRAY   &gperm,
			Int k, 
			Int prev = 0)
  {    
    A.init_perm(gperm);
    A.init_offset(k,prev);
    return A.offset;
  }//end populate_row_view
                        
  //Combine with (local)_back_solve
  //Add a selection for when we should us atomic
  //Add paramter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int back_solve_atomic(Int kid, Int k, Int top,
			Int brow,
			Int xnnz,
			Int ws_size,
			INT_1DARRAY          ws,
			ENTRY_1DARRAY        X,
			BASKER_MATRIX        L,
			INT_1DARRAY         gperm,
			const BASKER_STATS  stats
			)
		
  {
    //Int color_offset = 0;
    Int pattern_offset = ws_size;
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;

    Int Lend = L.scol + L.ncol;

    for(pp = 0; pp < xnnz; pp++)
      {
        j = ws[pattern_offset+top1];
        top1++;
        ws[j] = 0; //zeros out color
        //t = L.lpinv[j];
        t = gperm[j];
        
        //if(t != L.max_idx && (t >= L.scol) && (t < (L.scol+L.ncol)))
        if((t >= L.scol) && (t < (Lend)))
          {
            //ATOMIC : NOT NEEDED
            xj = X[j];
            Int local_offset = t-L.scol;
            //p2 = L.col_ptr[t+1-local_offset];
            //p = L.col_ptr[t-local_offset]+1;
            p2 = L.col_ptr[local_offset+1];
            p =  L.col_ptr[local_offset]+1;

            #ifdef BASKER_DEBUG_ATOMIC_BACKSOLVE
            printf("Atomic_tri_solve: Updating col: %d %d with val: %f \n",
                   j, t, xj);
            #endif
  
            for( ; p < p2; p++)
              {
                Int row_idx = L.row_idx[p];
                //if(L.row_idx[p] < (Lend))

                //printf("kid: %d row_idx: %d Lend: %d \n",
                //       kid, gperm[row_idx], Lend);
                if(gperm[row_idx] < Lend)
                  {          
                    //X[L.row_idx[p]] -= L.val[p] *xj;
                    X[row_idx] -= L.val[p]*xj;
                   
                    #ifdef BASKER_TIME_DETAIL
                    //stats.nfactor_sep_nonatomic_flop[kid] += 1;
                    #endif
                   
                  }
                else
                  {
                 
                    Kokkos::atomic_fetch_sub(&(X[L.row_idx[p]]),
                                             L.val[p]*xj);
                  
                    #ifdef BASKER_TIME_DETAIL
                    //stats.nfactor_sep_atomic_flop[kid] += 1;
                    #endif
         
                  }
              }//end for() over all nnz in a column

          }//end if() not permuted

      }//end for() over all nnz in LHS

    return 0;
  }//end back_solve_atomic

  //Combine with local_rec_factor ??
  //Add paramter for kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void column_tri_solve(Int kid, Int k, Int brow, Int ws_size,
			INT_1DARRAY ws,
			ENTRY_1DARRAY X,
			INT_1DARRAY gperm,
			BASKER_MATRIX_VIEW &A,
			BASKER_MATRIX &L,
			BASKER_MATRIX &U,
			BASKER_MATRIX &C,
			Int malloc_option,
			const BASKER_STATS  stats,
			int *enfo)
  {

    #ifdef BASKER_TIME_NFACTOR
    Kokkos::Impl::Timer timer_total;
    #endif
   
    Int bcol = U.scol;
    Int b = 0;
    Int j, xnnz, t;
    Int newsize;

    //Int scol = U.scol;
    //Int ecol = U.ecol;
       
    //Int color_offset = 0;
    Int pattern_offset = ws_size;
    //Int stack_offset = pattern_offset + ws_size;
    
    //Int cu_ltop = 0;
    //Int cu_utop = 0;
    Int top = ws_size;
    Int top1 = ws_size;
    
    Int lnnz = 0;
    Int unnz = 0;

    if(k != U.scol)
      {
        unnz = U.col_ptr[k-U.scol];
      }

    Int uunnz = U.nnz;
    Int llnnz = L.nnz;

    Int maxindex = L.max_idx;  
    Int lcnt = 0;
    Int ucnt = 0;

    #ifdef BASKER_TIME_NFACTOR
    double stat_time[6];
    for(Int i = 0; i < 6; i++)
      {
        stat_time[i] = 0;
      }
    #endif


    
    // #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
    ASSERT(top == ws_size);
    //note that X[i] may not be zero because of sharing
    //for(Int i = 0 ; i < ws_size; i++){ASSERT(X[i] == 0);}
    //for(Int i = 0; i <  ws_size; i++){ASSERT(ws[i] == 0 );}
    for(Int i = 0; i < ws_size; i++)
      {
        if(ws[i] != 0)
          {
            printf("---------------------------------------------ERROR--------------\n");
            printf("kid: %d k: %d  i: %d ws[i]=%d L.scol: %d \n", kid, k,  i, ws[i], L.scol);
            ASSERT(ws[i] == 0);
          }
      }
    
    //#endif

    //A.info();
    Int i = A.offset;    
    Int mi = A.m_offset;
    
    
 
    #ifdef BASKER_TIME_DETAIL
    Kokkos::Impl::Timer timer_reach;
    double t_i_reach = 0;
    #endif

    //Scan the whole column for entries in my 2D block
    for(i = A.offset; i < mi; i++)
      {        
        //if not in my block continue to scan
        if(A.good(i) == L.max_idx)
          {continue;}
       
        j = A.row_idx(i);
        //ATOMIC : NOT NEEDED READ
        X[j] = A.val(i);
        
        #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
        printf("kid: %d i: %d  val: %g  top: %d \n", kid,i, A.val(i), top);
        printf("kid: %d Nx in Ak %d %g %d color = %d \n",
               kid, j, X[j], brow,  ws[0 + j] );
        #endif

        
        #ifdef BASKER_TIME_NFACTOR
        Kokkos::Impl::Timer timer_l_r;
        #endif

        //if not colored find local reach
        if(ws[j] == 0)
          {      
            local_reach<Int, Entry, Exe_Space>(b, j, 
                     L, ws, gperm, &top, bcol, brow, ws_size);
          }

        #ifdef BASKER_TIME_NFACTOR
        t_i_reach += timer_l_r.seconds();
        #endif
        
      }//end for() over entries in my 2D block
    xnnz = ws_size - top;


    #ifdef BASKER_TIME_DETAIL
    stats.nfactor_sep_reach_time[kid] += timer_reach.seconds();
    //stats.nfactor_sep_reach_flop[kid] += xnnz;
    #endif


    #ifdef BASKER_TIME_NFACTOR
    printf("kid: %d k: %d reach_time: %f inner_reach_time: %f \n",
           kid, k, timer_reach.seconds(), t_i_reach);
    #endif
    


    #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
    for(i = top; i < ws_size; i++)
      {
        j = ws[pattern_offset+i];
        //t = L.lpinv[j];
        t = gperm[j];
        printf("found: %d  j: %d  t: %d \n",
               i, j , t);

      }
    #endif

    #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
    cout << "xnnz: " << xnnz << " ws_size: " << ws_size 
         << " top: " << top << endl;
    //ASSERT(xnnz <= ws_size);
    #endif



    //printf("kid: %d xnnz: %d \n", kid, xnnz);
    
    #ifdef BASKER_TIME_DETAIL
    Kokkos::Impl::Timer timer_solve;
    #endif

    //call an atomic backsolve
    back_solve_atomic<Int,Entry, Exe_Space>(kid, k, top,
                                            brow, xnnz,ws_size,ws,
                                            X,L, gperm, 
                                            stats); 
    
    #ifdef BASKER_TIME_DETAIL
    stats.nfactor_sep_solve_time[kid] += timer_solve.seconds();
    //printf("kid: %d k: %d back_solve: %f \n",
    //      kid, k, timer_solve.seconds());
    #endif

    #ifdef BASKER_TIME_NFACTOR
    Kokkos::Impl::Timer timer_book;
    #endif

    //count nnz will be need in this column
    for(i = top; i < ws_size; i++)
      {
        j = ws[pattern_offset + i];
        //t = L.lpinv[j];
        t = gperm[j];
        if(t == L.max_idx)
          {
            lcnt++;
          }
      }//for (i = top; i < ws_size)
    ucnt = ws_size - top - lcnt +1;


    //Note: go back and fix this up!!!
    //Need choices for Malloc
    if(unnz+ucnt >= uunnz)
      {
        if(malloc_option ==1)
          { //return 1;
          }
        newsize = uunnz*1.1 + 2*A.nrow+1;
        cout << "Reallocing U oldsize: " << uunnz 
             << " newsize " << newsize << endl;
      }

    
    //for nnz in LHS
    for( i = top; i < ws_size; i++)
      {
        j = ws[pattern_offset + i];
        //t = L.lpinv[j];
        t = gperm[j];

        // printf("considering loc: %d j: %d t:%d val: %e \n",
        //       pattern_offset+i, j, t, X[j]);
      
        #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
        printf("considering j: %d t:%d val: %e \n",
             j, t, X[j]);    
        #endif
        
        //if not zero
        if(X[j] != 0)
          {
            //Add to U
            //if(t != L.max_idx)
            if( (t != L.max_idx) && (t >= L.scol) && (t<(L.scol+L.ncol)))
              {
                #ifdef BASKER_DEBUG_COLUMN_TRISOLVE
                printf("kid: %d adding x[%d] to U\n", kid, j); 
                #endif

                //U.row_idx[unnz] = L.lpinv[j];
                U.row_idx[unnz] = gperm[j];
                //ATOMIC:: NOT NEEDED
                U.val[unnz] = X[j];
                unnz++;
                //ATOMIC:: NOT NEEDED
                X[j] = 0;
                C.union_bit[j] = false;
              }//if in U
            else
              {
                //New Move Type
                #ifdef MOVE_A_SMALL
                Int my_offset = 4*ws_size;
                ws[my_offset] = ws[my_offset]+1;
                printf("reduce vec, add %d at %d max: %d\n",
                       j, my_offset + ws[my_offset], ws.size()); 
                ws[my_offset + ws[my_offset]] = j;
                #endif
                C.union_bit[j] = true;
           
              }
 
          }//if nonzero
      }//over all X
    // printf("top: %d nnz: %d  \n", top, xnnz);
    ws[2*ws_size] = xnnz;
   
      U.col_ptr[k+1-bcol] = unnz;

      #ifdef BASKER_TIME_NFACTOR
      printf("kid: %d k: %d move_u: %f \n",
             kid, k, timer_book.seconds());
      #endif
      
      #ifdef BASKER_TIME_NFACTOR
      printf("kid: %d k: %d col_tri_total: %f \n",
             kid, k, timer_total.seconds());
      #endif

  }//end column_tri_solve()


  //Helper function to clear out nonzeros and reuse C
  //Note: Dimention of C needs to change
  //Note: Could possibliy be done in eithr column_trisolve or column_factor
  //Add paramter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void zero_col_matrix(BASKER_MATRIX C)
  {
    for(Int i = 0; i < C.col_ptr[1]; i++)
      {
        C.row_idx[i] = 0;
        C.val[i] = 0;
      }
  }//zero_col_matrix

  //Helper function.  
  //Does Reduce of of X and A
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void move_A(int kid, int k, int prev, int l, int lvl,
	      BASKER_MATRIX_VIEW  &A,
	      const BASKER_MATRIX &C,
	      BASKER_MATRIX &L,
	      INT_1DARRAY &gperm,
	      ENTRY_1DARRAY &X)
  {
    //populate row
    populate_row_view<Int, Entry, Exe_Space>
      (A, L, gperm, k, prev);
    
    //FIX later: need to reuse some of the temp workspace!!!!!
    Int *r_idx = new Int[X.size()]();

    //update X
    for(Int i = A.offset; i < A.m_offset; i++)
      {

        printf("good found: %d \n", A.good(i));
        if(A.good(i) != 0)
          {
            if(l == lvl)
              {
                //Move temp variable up or remove
                Int A_row = A.row_idx(i);
                
                #ifdef BASKER_DEBUG_MOVEA
                printf("Add (%d) A: %d, %f before value: %f kid: %d\n",
                       l, A_row, A.val(i), X[A_row], kid);
                #endif

                //ATOMIC:  !!!!!
                //Kokkos::atomic_fetch_add(&(X[A_row]), A.val(i));
                X[A_row] += A.val(i);
              }
          }//if not yet permute(touched)
        else
          {
            //Move temp variable up or remove
            Int A_row = A.row_idx(i);
            //Int A_row_p = L.lpinv[A_row];
            Int A_row_p = gperm[A_row];
            
            #ifdef BASKER_DEBUG_MOVEA
            printf("Considering: %d %d with start %d kid: %d\n", 
                   A_row, A_row_p, A.srow, kid);
            #endif

            if((A_row_p >= A.srow) && (A_row_p < A.srow+A.nrow))
              {
                #ifdef BASKER_DEBUG_MOVEA
                printf("add A: %d %f kid: %d \n", A_row, A.val(i), kid);
                #endif

                //Kokkos::atomic_fetch_add(&(X[A_row]), A.val(i));
                X[A_row] += A.val(i);
              }//if in my submatrix/view
          }//if permuted
      }//over all A in column

    
    //find nnz of X
    //Might want to chage this to some union op
    //Consider Kokkos bitsets
    Int nnz = 0;
    for(Int i = 0; i < X.size(); i++)
    //for(Int i = 0; i < A.nrow; i++)
      {
        printf("considering value: %f \n", X[i]);
        if(X[i] != 0)
          {
            printf("nnz found at: %d \n", i);
            //Int A_row_p = L.lpinv[i];
            Int A_row_p = gperm[i];
            //printf("nnz found at: %d %d compared to A.srow: %d \n",
            //       i, A_row_p, A.srow);
            if(A_row_p >= A.srow)
              {
                printf("taking found nnz at: %d \n", i);
                r_idx[nnz] = i;
                nnz++;
              }
          }
      }
    
    //Copy information into matrix
    C.col_ptr[0] = 0;
    C.col_ptr[1] = nnz;
    for(Int i = 0; i < nnz; i++)
      {
        C.row_idx[i] = r_idx[i];
        C.val[i] = X[r_idx[i]];
      }

    delete [] r_idx;
  }//end move_A

 
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void col_selective_sort(INT_1DARRAY ws,
                          Int ws_start,
                          Int ws_end,
                          Int e_start,
                          Int e_end,
                          INT_1DARRAY gperm
                          )
  {

    printf("select sort called: start: %d end: %d \n",
           ws_start, ws_end);

    Int new_end = ws_end;

    //compare if and change and move smallest first
    for(Int i = ws_end-1; i > ws_start; i--)
      {
        //first compare if in range
        Int elem = ws[i];
        //Note: might want to check range
        if((gperm[elem] >= e_start) && (gperm[elem]<= e_end))
          {
            new_end--;
            continue;
          }
        else
          {
            if(elem < ws[i-1])
              {
                ws[i] = ws[i-1];
                ws[i-1] = elem;
              }
          }
      }

    //insert sort rest
    for(Int i = ws_start+2; i < new_end; i++)
      {
        Int j = i;
        Int elem = ws[i];
        while(elem < ws[j-1])
          {
            ws[j] = ws[j-1];
            j--;
          }
        ws[j] = elem;
      }

  }//end col_selective_sort
    



  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void move_X_small(int kid, int k, int prev, int l, int lvl,
                    BASKER_MATRIX_VIEW &A,
                    BASKER_MATRIX &C,
                    BASKER_MATRIX &L,
                    INT_1DARRAY gperm,
                    ENTRY_1DARRAY X,
                    INT_1DARRAY ws_one,
                    Int   ws_one_offset,
                    INT_1DARRAY ws_two,
                    Int  ws_two_offset)
  {

    populate_row_view<Int, Entry, Exe_Space>
      (A, L, gperm, k, prev);

    Int A_Start = A.scol;
    Int A_End = A.scol+A.ncol;

    Int nnz = 0;
    Int my_ws_one_offset = 4*ws_one_offset;
    Int my_ws_one_nnz = ws_one[my_ws_one_offset];
    Int my_ws_two_offset = 4*ws_two_offset;
    Int my_ws_two_nnz = ws_two[my_ws_two_offset];
    printf("ws_start: %d \n", ws_one_offset+1);
    printf("ws_one_nnz: %d \n", my_ws_one_nnz);
    printf("ws_two_nnz: %d \n", my_ws_two_nnz);


    //sort lists
    col_selective_sort<Int,Entry,Exe_Space>
      (ws_one, my_ws_one_offset+1, 
                       my_ws_one_offset+my_ws_one_nnz+2,
                       A_Start, A_End, gperm);
    col_selective_sort<Int,Entry,Exe_Space>
      (ws_two, my_ws_two_offset+1, 
                       my_ws_two_offset+my_ws_two_nnz+2,
                       A_Start, A_End, gperm);
    
    
    //Debug
    printf("sorted ws_one \n");
    for(Int i = my_ws_one_offset+1; i < my_ws_one_offset+my_ws_one_nnz; i++)
      {
        printf("%d, ", ws_one[i]);
      }
    printf("\n");
    printf("sorted ws_two \n");
    for(Int i = my_ws_two_offset+1; i < my_ws_two_offset+my_ws_two_nnz; i++)
      {
        printf("%d, ", ws_two[i]);
      }
    printf("\n");


 
    //merge ws_one and ws_two
    Int ptr_m = 0;
    Int ptr_1 = my_ws_one_offset+1;
    Int ptr_2 = my_ws_two_offset+1;
    for(Int k = 0 ; k < (my_ws_one_nnz+my_ws_two_nnz); k++)
      {
        if(ptr_1 == my_ws_one_nnz)
          {
            C.row_idx[ptr_m] = ws_one[ptr_1];
            C.val[ptr_m] = X[ws_one[ptr_1]];
            X[ws_one[ptr_1]] = 0;
            ws_one[ptr_1] = 0;
            ptr_m++;
            ptr_1++;
            continue;
          }
        if(ptr_2 == my_ws_two_nnz)
          {
            C.row_idx[ptr_m] = ws_two[ptr_2];
            C.val[ptr_m] = X[ws_two[ptr_2]];
            X[ws_two[ptr_2]] = 0;
            ws_two[ptr_2] = 0;
            ptr_m++;
            ptr_2++;
            continue;
          }
        //Cases
        if(ws_one[ptr_1] == ws_two[ptr_2])
          {
            C.row_idx[ptr_m] = ws_one[ptr_1];
            C.val[ptr_m] = X[ws_one[ptr_1]];
            X[ws_one[ptr_1]] = 0;
            ws_one[ptr_1] = 0;
            ws_two[ptr_2] = 0;
            ptr_m++;
            ptr_1++;
            ptr_2++;
          }
        else if(ws_one[ptr_1] < ws_two[ptr_2])
          {
            C.row_idx[ptr_m] = ws_one[ptr_1];
            C.val[ptr_m] = X[ws_one[ptr_1]];
            X[ws_one[ptr_1]] = 0;
            ws_one[ptr_1] = 0;
            ptr_m++;
            ptr_1++;
          }
        else if(ws_one[ptr_1] > ws_two[ptr_2])
          {
            C.row_idx[ptr_m] = ws_two[ptr_2];
            C.val[ptr_m] = X[ws_two[ptr_2]];
            X[ws_two[ptr_2]] = 0;
            ws_two[ptr_2] = 0;
            ptr_m++;
            ptr_1++;
          }
      }
    Int C_nnz = ptr_m;

    
    Int Annz = 0;
    for(Int i = A.offset; i < A.m_offset; i++)
      {

        if(gperm[A.row_idx(i)] == L.max_idx )
          {
            if(l == lvl)
              {
                //Move temp variable up or remove
                Int A_row = A.row_idx(i);
                
                #ifdef BASKER_DEBUG_MOVEA
                printf("Add (%d) A: %d, %f before value: %f kid: %d\n",
                       l, A_row, A.val(i), X[A_row], kid);
                #endif

                X[A_row] += A.val(i);

                ws_one[my_ws_one_offset+Annz] = A_row;
                C.union_bit[A_row] = false;
                Annz++;


              }
          }//if not yet permute(touched)
        else
          {
            //Move temp variable up or remove
            Int A_row = A.row_idx(i);
            //Int A_row_p = L.lpinv[A_row];
            Int A_row_p = gperm[A_row];
            
            #ifdef BASKER_DEBUG_MOVEA
            printf("Considering: %d %d with start %d kid: %d\n", 
                   A_row, A_row_p, A.srow, kid);
            #endif

            if((A_row_p >= A.srow) && (A_row_p < A.srow+A.nrow))
              {
                #ifdef BASKER_DEBUG_MOVEA
                printf("add A: %d %f kid: %d \n", A_row, A.val(i), kid);
                #endif

                X[A_row] += A.val(i);
              
                ws_one[my_ws_one_offset+Annz] = A_row;
                C.union_bit[A_row] = false;
                Annz++;

                //printf("adding A 2: %d \n", A_row);
                
              }//if in my submatrix/view
          }//if permuted
      }//over all A in column
    

    printf("before adding A\n");
    C.print();



    //sort ws_one
    col_selective_sort<Int,Entry,Exe_Space>
      (ws_one, my_ws_one_offset, 
                       my_ws_one_offset+Annz,
                       A_Start, A_End, gperm);


    //Special Merge
    ptr_1 = 0; //ptr to C
    ptr_2 = my_ws_one_offset; //ptr to other
    for(Int k = 0; k < (C_nnz+Annz); k++)
      {
        if(ptr_2 == (my_ws_one_offset+Annz))
          {
            //all done
            break;
          }
        if(ptr_1 == C_nnz)
          {
            //just add the remainder
            C.row_idx[ptr_m] = ws_one[ptr_2];
            C.val[ptr_m] += X[ws_one[ptr_2]];
            X[ws_one[ptr_2]] = 0;
            ws_one[ptr_2] = 0;
            ptr_m++;
            ptr_2++;
          }
        //Cases
        if(C.row_idx[ptr_1] == ws_one[ptr_2])
          {
            //already been added just upate value
            C.val[ptr_1] += X[ws_one[ptr_2]];
            ptr_1++;
            ptr_2++;
          }
        if(C.row_idx[ptr_1] < ws_one[ptr_2])
          {
            //Advance C to correct point
            ptr_1++;
          }
        if(C.row_idx[ptr_1] > ws_one[ptr_2])
          {
            //not been added yet
            C.row_idx[ptr_m] = ws_one[ptr_2];
            C.val[ptr_m] += X[ws_one[ptr_2]];
            X[ws_one[ptr_2]] = 0;
            ws_one[ptr_2] = 0;
            ptr_m++;
            ptr_2++;
          }


      }
   
    C.col_ptr[0] = 0;
    C.col_ptr[1] = ptr_m;

    ws_one[my_ws_one_offset] = 0;
    ws_two[my_ws_two_offset] = 0;

  }//end move_A_small()


  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void move_A_new(int kid, int k, int prev, int l, int lvl,
		  BASKER_MATRIX_VIEW A,
		  const BASKER_MATRIX C,
		  BASKER_MATRIX L,
		  INT_1DARRAY gperm,
		  ENTRY_1DARRAY X,
		  INT_1DARRAY ws_one,
		  Int ws_one_offset,
		  INT_1DARRAY ws_two,
		  Int ws_two_offset)
  {
    Int nnz = 0;
    
    //populate row
    populate_row_view<Int, Entry, Exe_Space>
      (A, L, gperm, k, prev);
    
    //update X
    for(Int i = A.offset; i < A.m_offset; i++)
      {

        if(gperm[A.row_idx(i)] == L.max_idx )
          {
            if(l == lvl)
              {
                //Move temp variable up or remove
                Int A_row = A.row_idx(i);
                
                #ifdef BASKER_DEBUG_MOVEA
                printf("Add (%d) A: %d, %f before value: %f kid: %d\n",
                       l, A_row, A.val(i), X[A_row], kid);
                #endif

                //ATOMIC:  !!!!!
                //Kokkos::atomic_fetch_add(&(X[A_row]), A.val(i));
                X[A_row] += A.val(i);
                //Add right away to matrix
                C.row_idx[nnz] = A_row;
                C.val[nnz] = X[A_row];
                X[A_row] = 0;
                C.union_bit[A_row] = false;
                nnz++;

                //printf("adding A 1 : %d \n", A_row);

              }
          }//if not yet permute(touched)
        else
          {
            //Move temp variable up or remove
            Int A_row = A.row_idx(i);
            //Int A_row_p = L.lpinv[A_row];
            Int A_row_p = gperm[A_row];
            
            #ifdef BASKER_DEBUG_MOVEA
            printf("Considering: %d %d with start %d kid: %d\n", 
                   A_row, A_row_p, A.srow, kid);
            #endif

            if((A_row_p >= A.srow) && (A_row_p < A.srow+A.nrow))
              {
                #ifdef BASKER_DEBUG_MOVEA
                printf("add A: %d %f kid: %d \n", A_row, A.val(i), kid);
                #endif

                //Kokkos::atomic_fetch_add(&(X[A_row]), A.val(i));
                X[A_row] += A.val(i);
                //Add right away to matrix
                C.row_idx[nnz] = A_row;
                C.val[nnz] = X[A_row];
		X[A_row] = 0;
                C.union_bit[A_row] = false;
                nnz++;

                //printf("adding A 2: %d \n", A_row);
                

              }//if in my submatrix/view
          }//if permuted
      }//over all A in column

    
    //Add ones for ws_onw
    for(Int i = ws_one[ws_one_offset]; i > 0; i--)
      {
        Int row_idx = ws_one[ws_one_offset-i];
        //printf("ws one consider: %d %d  \n", ws_one_offset-i, row_idx);
        ws_one[ws_one_offset-i] = 0;
        if(C.union_bit[row_idx] == true)
          {
            printf("ws one adding: %d \n", row_idx);
            C.row_idx[nnz] = row_idx;
            C.val[nnz] = X[row_idx];
            
            C.union_bit[row_idx] = false;
            nnz++;
          }
        X[row_idx] = 0;
      }
    ws_one[ws_one_offset] = 0;
    
    //Add ones from ws_two
    for(Int i = ws_two[ws_two_offset] ; i > 0; i--)
      {
       
        Int row_idx = ws_two[ws_two_offset-i];
        //printf("ws two consider: %d %d \n", ws_one_offset-i, row_idx);
        ws_two[ws_two_offset-i] = 0;
        if(C.union_bit[row_idx] == true)
          {
            printf("ws two adding: %d \n", row_idx);
            C.row_idx[nnz] = row_idx;
            C.val[nnz] = X[row_idx];
            C.union_bit[row_idx] = false;
            nnz++;
          }
        X[row_idx] = 0;
      }
    ws_two[ws_two_offset] = 0;
    C.col_ptr[1] = nnz;
    
    //Clear bit vector
    for(Int i = 0; i < C.col_ptr[1]; i++)
      {
        Int r_idx = C.row_idx[i];
        C.union_bit[r_idx] = false;
      }

  }//end move_A_new()

  //Helper function, clear out entry work space
  //Need to come up with a trick since we can't clear in reduce as 
  //muliple threads will be accessing the X in the reduce
  //Add kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void clear_x(ENTRY_1DARRAY &X, Int size)
  {
    for(Int i = 0; i < size; i++)
      {
        X[i] = 0;
      }
  }//end clear_x

 


  //Factors a single column based on a matrix view
  //Maybe can be integrated into rec_factor as a special case
  //Add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int col_factor(Int b, Int k, Int lval, Int uval,
		 Int bcol, Int brow,
		 INT_1DARRAY ws,
		 ENTRY_1DARRAY X,
		 INT_1DARRAY &gperm,
		 BASKER_MATRIX_VIEW A,
		 BASKER_MATRIX &L,
		 BASKER_MATRIX &U,
		 Int malloc_option,
		 const BASKER_STATS stats,
		 int *einfo)
  {
    Int i,j;
    INT_1DARRAY tptr, color, pattern, stack;
    Int top, top1, maxindex, t; 
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;
   
    Int newsize;
    Entry pivot, value;
    Entry absv, maxv;
    Int ws_size = A.base->nrow;
    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    //Int scol = L.scol;
    //Int ecol = L.ecol;
    
    Int color_offset = 0;
    Int pattern_offset = ws_size;
    //Int stack_offset = pattern_offset + ws_size;

    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
    top1 = ws_size;
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_COLUMN_FACTOR
    cout << "col_rec_fact" << endl;
    //cout << "scol: " << scol << " ecol: " << ecol << endl;
    cout << "L.nnz " << llnnz << " U.nnz " << uunnz << endl;
    #endif
  
        {
          #ifdef BASKER_DEBUG_COLUMN_FACTOR 
	  cout << "-------------------------- k = " 
               << k << " -----------------"
	       << endl;
           #endif

	  value = 0.0;
	  pivot = 0.0;
	  maxindex = A.base->max_idx;  
	  lcnt = 0;
	  ucnt = 0;

          #ifdef BASKER_DEBUG_COLUMN_FACTOR
          ASSERT(top == ws_size);
          for(i = 0 ; i < ws_size; i++){ASSERT(X[i] == 0);}
	  for(i = 0; i <  ws_size; i++){ASSERT(ws[i] == 0 );}
          #endif


          #ifdef BASKER_DEBUG_COLUMN_FACTOR
          printf("--------------PAST ASSERT----------- \n");
          printf("A.col(k): %d  A.col(k+1): %d \n", 
                 A.col_ptr(k), A.col_ptr(k+1));
          #endif


          #ifdef BASKER_TIME_DETAIL
          Kokkos::Impl::Timer timer_reach;
          #endif

	  for(i = A.col_ptr(k); i < A.col_ptr(k+1); i++)
	    {
              #ifdef BASKER_DEBUG_COLUMN_FACTOR
              A.good(i);
              #endif
            
	      j = A.row_idx(i);
	  
              #ifdef BASKER_DEBUG_COLUMN_FACTOR
              printf("j: %d i: %d \n", j, i);
              #endif

              X[j] = A.val(i);

              #ifdef BASKER_DEBUG_COLUMN_FACTOR
              printf("i: %d j: %d  val: %g  top: %d \n", 
                     i, gperm[j],A.val(i), top);
              printf("Nx in Ak %d %g %d color = %d \n",
                      j, X[j], brow,  
                     ws[color_offset + j ] );
               #endif
             
               if(ws[color_offset + j ] == 0)
               {            
                 local_reach<Int, Entry, Exe_Space>
                   (b, j, L, ws, gperm, &top,bcol, brow, 
                    ws_size);
		}
          }//end reachable (for i = A.col)
	  xnnz = ws_size - top;

          #ifdef BASKER_TIME_DETAIL
          stats.nfactor_sep_reach_time[b] += timer_reach.seconds();
          //stats.nfactor_sep_reach_flop[b] += xnnz;
          #endif

          #ifdef BASKER_DEBUG_COLUMN_FACTOR
          cout << "xnnz: " << xnnz << " ws_size: " 
               << ws_size << " top: " << top << endl;
          //ASSERT(xnnz <= ws_size);
          #endif
          
          #ifdef BASKER_TIME_DETAIL
          Kokkos::Impl::Timer timer_solve;
          #endif

          /*
          back_solve_atomic<Int,Entry, Exe_Space>
            ( b, k, top,
                   brow, xnnz, ws_size,
              ws, X,  L, gperm, stats); 
          */
          
          back_solve<Int,Entry, Exe_Space>
            ( b, k, top,
                   brow, xnnz, ws_size,
              ws, X,  L, gperm, stats); 
          
          #ifdef BASKER_TIME_DETAIL
          stats.nfactor_sep_solve_time[b] += timer_solve.seconds();
          #endif

          #ifdef BASKER_DEBUG_COLUMN_FACTOR
          cout << " top " << top << " ws_size " << ws_size << endl;
          #endif          

          //find pivot
          maxv = 0.0;
          for(i = top; i < ws_size; i++)
            {
              j = ws[pattern_offset + i];
              //t = L.lpinv[j];
              t = gperm[j];
              value = X[j];
              absv = abs(value);
              if(t == L.max_idx)
                {
                  lcnt++;
                  if(absv > maxv) 
                    {
                      maxv = absv;
                      pivot = value;
                      maxindex = j;                
                    }
                }
            }//for (i = top; i < ws_size)          
          ucnt = ws_size - top - lcnt +1;


          printf("pivot found: %f \n", pivot);
          
          if((maxindex == L.max_idx) || (pivot == 0))
            {
              cout << "Error: Matrix is singular" << endl;
              cout << "MaxIndex: " << maxindex << " pivot " 
                   << pivot << endl;
              *einfo = k;
              return 2;
            }          
          //L.lpinv[maxindex] = k;
          gperm[maxindex] = k;
          //printf("TAG1:  maxindex = %d k = %d \n", maxindex, k);
          //printf("TAG2: \n");

          #ifdef BASKER_DEBUG_COLUMNFACTOR
          if(maxindex != k)
            {
              cout << "Permuting Pivot: " << k << " as row " 
                   << maxindex << endl;
            }
          #endif
          
          //NOTE: Fix:  need to fix how we do malloc
          //Come back to
          if(lnnz + lcnt >= llnnz)
            {
              if(malloc_option == 1)
                {//return 1; 
                }
              newsize = lnnz * 1.1 + 2 *A.nrow + 1;
              cout << "Reallocing L oldsize: " << llnnz 
                   << " newsize: " << newsize << endl;
            }
          if(unnz+ucnt >= uunnz)
            {
              if(malloc_option ==1)
                {//return 1;
                }
              newsize = uunnz*1.1 + 2*A.nrow+1;
              cout << "Reallocing U oldsize: " << uunnz 
                   << " newsize " << newsize << endl;
            }
          
          L.row_idx[lnnz] = maxindex;
          L.val[lnnz] = (Entry) 1.0;
          lnnz++;
     
          Entry lastU = (Entry) 0.0;
  
          //For every nnz in LHS
          for( i = top; i < ws_size; i++)
            {
              j = ws[pattern_offset + i];
              //t = L.lpinv[j];
              t = gperm[j];
              
              #ifdef BASKER_DEBUG_COLUMN_FACTOR
              cout << "j " << j << " t " << t << endl;
              #endif
              
              //if not zero
              if(X[j] != 0)
                {
                  #ifdef BASKER_DEBUG_COLUMN_FACTOR
                  cout << "found value " << X[j] << " at " << j << endl;
                  #endif
                  
                  if(t != L.max_idx)
                    {
                      if(t < k)
                        {
                          printf("U insert: %f at %d \n",
                                 X[j], gperm[j]);
                          //U.row_idx[unnz] = L.lpinv[j];
                          U.row_idx[unnz] = gperm[j];
                          U.val[unnz] = X[j];
                          unnz++;
                        }
                      else
                        {
                          lastU = X[j];
                        }
                    }
                  else if (t == L.max_idx)
                    {
                      #ifdef BASKER_DEBUG_COLUMN_FACTOR
                      printf("inserting %f at %d into %d \n", 
                             X[j]/pivot, j, lnnz ); 
                      #endif

                      L.row_idx[lnnz] = j;
                      L.val[lnnz] = X[j]/pivot;
                      lnnz++;
                    }
                }//end if() not zero             
              //move inside if not zero..... extra set ops not needed
              X[j] = 0;
            }//if(x[j-brow] != 0)

          //Fill in last element of U
          U.row_idx[unnz] = k;
          U.val[unnz] = lastU;
          unnz++;

          xnnz = 0;
          top = ws_size;
          

          printf("setting col: %d %d %d %d\n",
                 k-bcol, cu_ltop, k+1-bcol, lnnz);
          L.col_ptr[k-bcol] = cu_ltop;
          L.col_ptr[k+1-bcol] = lnnz;
          cu_ltop = lnnz;
          
          U.col_ptr[k-bcol] = cu_utop;
          U.col_ptr[k+1-bcol] = unnz;
          cu_utop = unnz;

          #ifdef BASKER_DEBUG_COLUMN_FACTOR
          printf("col_fact k: %d Unnz: %d   Lnnz: %d \n",
                 k, unnz, lnnz);
          #endif
	}//over all columns (only 1 column)

    #ifdef BASKER_DEBUG_COLUMN_FACTOR
	print_factor<Int,Entry,Exe_Space>(L,U);
    #endif
   
    return 0;
  }// end column_factor

  //Helper function.  Does a reduce on perm every level
  //Reduces does not need to happen every level.... need to figure this out
  //Add parameter kid
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void reduce_perm(BASKER_MATRIX P,
		   BASKER_MATRIX C)
  {

    for(Int i = 0; i < P.lpinv.size(); i++)
      {
        //May want to think if should use atomic is odd cases
        if(C.lpinv[i] < C.max_idx)
          {
            P.lpinv[i] = C.lpinv[i];
          }
      }
  }//end reduce_perm()
                   
  //Kokkos local_extend funct
  //Need to change name to something that matches
  //Need to move more operations out of operator function
  template <class Int, class Entry, class Exe_Space>
  struct local_extend_funct
  {
    //Kokkos Typedefs
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                                        execution_space ;
    typedef  Kokkos::TeamPolicy<Exe_Space>                   team_policy;
    typedef  typename team_policy::member_type               team_member;
    #endif
    
    MATRIX_VIEW_2DARRAY A;
    MATRIX_2DARRAY L;
    MATRIX_2DARRAY U;
    INT_2DARRAY WS;
    ENTRY_2DARRAY X;
    INT_2DARRAY S; //schedule
    Int lvl;
    Int slvl;
    BASKER_TREE tree;
    MATRIX_2DARRAY C;

    //New global perm
    INT_2DARRAY gperm;

    //Stats
    BASKER_STATS stats;
    
    local_extend_funct()
    {}

    local_extend_funct(MATRIX_VIEW_2DARRAY  _A, 
		       MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
                       INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_2DARRAY _S, 
                       Int _lvl, Int _slvl, BASKER_TREE _tree, MATRIX_2DARRAY _C)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
      slvl =_slvl;
      tree = _tree;
      C = _C;
    }

     local_extend_funct(MATRIX_VIEW_2DARRAY  _A, 
                        MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
                        INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_2DARRAY _S, 
                        Int _lvl, Int _slvl, BASKER_TREE _tree, MATRIX_2DARRAY _C, 
                        INT_2DARRAY _gperm)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
      slvl =_slvl;
      tree = _tree;
      C = _C;
      gperm = _gperm;
    }

    local_extend_funct(MATRIX_VIEW_2DARRAY  _A, 
		       MATRIX_2DARRAY _L, MATRIX_2DARRAY _U,  
		       INT_2DARRAY _WS, ENTRY_2DARRAY _X, INT_2DARRAY _S, 
		       Int _lvl, Int _slvl, BASKER_TREE _tree, MATRIX_2DARRAY _C, 
		       INT_2DARRAY _gperm, BASKER_STATS _stats)
    {
      A = _A;
      L = _L;
      U = _U;
      WS = _WS;
      X = _X;
      S = _S;
      lvl = _lvl;  
      slvl =_slvl;
      tree = _tree;
      C = _C;
      gperm = _gperm;
      stats = _stats;
    }

    //KOKKOS_INLINE_FUNCTION
    BASKER_INLINE
    void operator()(const team_member &thread) const
    {
      #ifdef BASKER_TIME_DETAIL
      Kokkos::Impl::Timer timer_sep;
      #endif

      //#ifdef BASKER_TIME
      double stat_time[6];
      double stat_tri_time[6];
      for(Int i = 0; i < 6; i++)
        {
          stat_time[i] = 0;
          stat_tri_time[i] = 0;
        }
      //#endif

      //need to add based on Exe_Space
      //int tid = Kokkos::OpenMP::hardware_thread_id();
      Int tid = (Int) Exe_Space::hardware_thread_id();
      Int kid = (Int) (thread.league_rank() * thread.team_size())+ thread.team_rank();
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());

      //#ifdef BASKER_DEBUG_EXTEND
      printf("\n\n league_rank: %i team_rank:  %i  league_size: %i  team_size: %i \n OMPRank: %d  KRank: %d  Team_Leader: %d  \n\n ", 
             thread.league_rank(), thread.team_rank(),
             thread.league_size(), thread.team_size(), 
             tid, kid, team_leader);
      thread.team_barrier(); //EXP
      //#endif

      //Need to try to reduce local variables
      Int L_col = 0;
      Int L_row = 0;
      
      Int X_col = 0;
      
      Int U_col = S[lvl+1][kid];
      Int U_row = 0;

      MATRIX_1DARRAY lol= U[U_col];

      Int scol = U[U_col][U_row].scol;
      Int ecol = U[U_col][U_row].ecol;

      bool roll_back = false;
      BASKER_MATRIX *base;
      BASKER_MATRIX col_matrix;


      {
        //over each column
        Int myend = ecol;
        //  if(lvl == 1)
   
        for(Int k = scol; k< myend; k++)
        {


          /*
          printf("\n\n\n");
          printf("---------------------------- start k=%d -------------\n",
                 k);
          printf("\n\n");
          */

          ENTRY_1DARRAY x = X[team_leader]; //Comeback and fix so all not using
          //for each level down do a tri-solve
          for(Int l = 0 ; l <= lvl; l++)
            {

              //I can do work because I am the correct thread
              //should I do a trisolve
              #ifdef BASKER_DEBUG_EXTEND
              //thread.team_barrier();
              #endif
              
              //--------------------TRI SOLVE-----------------------//
              //check module for all level p
               #ifdef BASKER_TIME_EXTEND
               Kokkos::Impl::Timer  timer_tri_solve;
               #endif
              
              //if((kid%(l+1)) == 0)
              if(kid%((Int)pow(2,l)) ==0)
                {
                  
                  #ifdef BASKER_DEBUG_EXTEND
                  printf("--------------STARTING TRI-SOLVE: %d %d ----------\n", 
                         kid, l);
                  #endif

                  #ifdef BASKER_TIME_DETAIL
                  Kokkos::Impl::Timer  timer_tri_solve_start;
                  #endif
             

                  //Get local offsets
                  L_col = S[l][kid];
                  L_row = 0;

                  X_col = S[l][kid];

                  U_col = U_col; //fixed
                  //fix this for lvls higher than 2
                  if(lvl == 0)
                     {
                       U_row = thread.team_rank();
                     }
                  else
                    {
                      U_row  = S[l][kid]%(U[U_col].size());
                    }
                  

                  #ifdef BASKER_DEBUG_EXTEND
                  printf("kid: %d u_col: %d u_row: %d \n",
                         kid, U_col, U_row);
                
                  #endif

                  BASKER_MATRIX_VIEW &LA = A[U_col][U_row];
                  BASKER_MATRIX      &LL = L[L_col][L_row];
                  BASKER_MATRIX      &LU = U[U_col][U_row];
                  INT_1DARRAY         ws = WS[L_col];
                 
                  
                  Int brow = LL.srow;
                  Int ws_size = L[0][0].nrow; 
                  Int malloc_option = 0;
                  int einfo = 0;


                  #ifdef BASKER_DEBUG_EXTEND
		  /*
                  printf("kid: %d  lvl: %d  l: %d L %d %d  U %d %d Lsize %d %d  Usize %d %d \n", 
                         kid, lvl, l, L_col, L_row, U_col, U_row,
                         LL.col_ptr.size(), LL.val.size(), 
                         LU.col_ptr.size(), LU.val.size());
		  */
                  #endif

                  
                  //get my column view (why did I call row???)
                  populate_row_view<Int, Entry, Exe_Space>
                    (LA, LL, gperm[0], k, 0);
                 
                 #ifdef BASKER_TIME_DETAIL
                  stat_tri_time[0] += timer_tri_solve_start.seconds();
                 #endif
             

                  #ifdef BASKER_TIME_DETAIL
                  Kokkos::Impl::Timer  timer_tri_solve_fill;
                  #endif
             
                  
                  //if(kid == 0)
                    {

                      #ifdef BASKER_DEBUG_EXTEND
                      //LL.print();
                      LA.info();
                      LA.base->info();
                      //LA.base->print();
                      #endif


                      LL.fill(); //will only fill if not yet been filled
                      LU.fill();
                      //do my local tri_solve
                  

                      BASKER_MATRIX &myC = C[kid/2][0];

                      #ifdef BASKER_TIME_DETAIL
                      stat_tri_time[1] += timer_tri_solve_fill.seconds();
                      #endif

                      #ifdef BASKER_TIME_DETAIL
                      Kokkos::Impl::Timer  timer_tri_solve_col;
                      #endif
              
                      column_tri_solve<Int, Entry, Exe_Space>
                        (kid,k,  brow, ws_size,
                         ws,
                         x,
                         gperm[0],
                         LA,
                         LL,
                         LU,
                         myC,
                         malloc_option,
                         stats,
                         &einfo);

                      #ifdef BASKER_TIME_DETAIL
                      stat_tri_time[2] += timer_tri_solve_col.seconds();
                      #endif

                    }
                    //set my view offset to zero.... we don't know the future
                    LA.offset  = 0 ;

                    //Do a reduction on perm
                    //We only need to do this once!!!!!
                    //Need to move this so only done once!!!!
                    //We should do this staticly before this functor call
                    Int p_L_col = S[l+1][kid];
                    BASKER_MATRIX &LP = L[p_L_col][L_row];

                    //Shoudl not need now have global perm
                    //reduce_perm<Int, Entry, Exe_Space>(LP, LL);
                    
                    //update NNZ (Done with all factoring
                    //No --- we would want to keep this around to know zeros
                    
               
                }//end if() kid... I should do try solve

                         
              ////////// --- NEEDED BEFORE REDUCE----/////////////////
              //////////--------------------------/////////////////
              thread.team_barrier();   //EXP
              ///////////------------------------///////////////


              #ifdef BASKER_TIME_EXTEND
              stat_time[0] += timer_tri_solve.seconds();
              //printf("kid: %d lvl: %d Tri-Solve-Time: %f \n",
              //       kid, l, timer_tri_solve.seconds());
              #endif

              #ifdef BASKER_TIME_EXTEND
              Kokkos::Impl::Timer timer_reduce;
              #endif

              //-----------------REDUCE---------------------///
              //Should I do the local reduce
              Int max_level_cut = pow(2,l+1);
              if((kid%(max_level_cut)) == 0)
                {
                  
                  //Get Parameters
                  U_col = U_col; //fixed
                  if(lvl == 0)
                    {
                      //U_row = thread.team_rank();
                      U_row = 2;
                    }
                  else
                    {
                      //U_row  = S[l][kid];
                      U_row = S[l+1][kid]%(U[U_col].size());
                    }
                  
                  #ifdef BASKER_DEBUG_EXTEND
                  printf("\n----------------------CHANGED POINTER kid: %d-----------------------\n A: %d %d \n\n", kid, U_col, U_row);
                  #endif
                  
                  BASKER_MATRIX_VIEW &UA = A[U_col][U_row];
                  Int p_L_col = S[l+1][kid];
                  BASKER_MATRIX      &LL = L[p_L_col][L_row];   
                  

		  Int temp_temp_idx = kid+pow(2,l);
                  Int ws_index = S[l][temp_temp_idx];
                  //printf("ws_two_index: %d \n" , ws_index);

                  INT_1DARRAY &ws_one = WS[L_col];
                  INT_1DARRAY &ws_two = WS[ws_index];
                  
                  Int ws_one_offset = 2*L[0][0].nrow;
                  Int ws_two_offset = 2*L[0][0].nrow;

                  
                  //NOTE: THIS INDEX NEEDS TO BE FIXED
                  //printf("myc %d \n", kid);
                  //MATRIX &myC = C[kid][0];
                 
                  BASKER_MATRIX &myC = C[kid/2][0];
                                
                  //merge X and A
                
                  /*
                  move_A<Int, Entry, Exe_Space>
                    (kid, k, 0, l, lvl, 
                     UA, myC, LL, gperm[0],
                     x);
		  
                  */

		  /*
                     //Faster but need to work out bug
                  move_A_new<Int,Entry,Exe_Space>
                    (kid, k, 0, l, lvl,
                     UA, myC, LL, gperm[0],
                     x, 
                     ws_one, ws_one_offset,
                     ws_two, ws_two_offset);
		  */
                  
                  #ifdef MOVE_A_SMALL
                  printf("move small called \n");
                  move_X_small<Int,Entry,Exe_Space>
                    (kid, k, 0, l, lvl,
                     UA, myC, LL, gperm[0],
                     x, 
                     ws_one, L[0][0].nrow,
                     ws_two, L[0][0].nrow);
                  #else
                  printf("move_A called \n");
                  move_A<Int, Entry, Exe_Space>
                    (kid, k, 0, l, lvl, 
                     UA, myC, LL, gperm[0],
                     x);
                  #endif
                  
                  
                  //set matrix to tempC matrix
                  base = UA.base;
                  //UA.base = &(C[kid][0]);
                  UA.base =  &(C[kid/2][0]);
                  //UA.k_offset = UA.scol;
                  UA.k_offset = k;
                  UA.offset = 0;
                  roll_back = true; //tell us that it we need to roll back to old A
                  
                }//end should I do reduce

             
              ////---------NEEDED BARRIER BEFORE FACTOR----------/////////////
              //thread.team_barrier();
              //////////--------------------------------------//////////////
             
              #ifdef BASKER_TIME_EXTEND
              stat_time[1] += timer_reduce.seconds();
              //printf("kid: %d lvl: %d time_reduce: %f \n",
              //       kid, l, timer_reduce.seconds());
              #endif

              #ifdef BASKER_TIME_EXTEND
              Kokkos::Impl::Timer timer_factor;
              #endif

              ///-------------Col Factor--------------------//
              //ed lif(kid == 100)
              {
 
                if((l== lvl) && (thread.team_rank() == 0))
                  {
                    //clear out any old updates that we will not need anymore
                                      
                    //Local offsets
                    L_col = S[lvl+1][kid];
                    L_row = 0;

                    #ifdef BASKER_DEBUG_EXTEND
                    printf("Factor Test.  id: %d kid: %d \n", 
                           thread.team_rank(), kid);
                    printf("\n----------------------FACTOR  kid: %d-----------------------\n A: %d %d L %d %d  \n\n", kid, U_col, U_row, L_col, L_row);
                    #endif
                 
                    BASKER_MATRIX_VIEW &UA = A[U_col][U_row];

                    BASKER_MATRIX      &LL = L[L_col][L_row];
                    BASKER_MATRIX      &LU = U[U_col][U_row];
                    INT_1DARRAY &ws = WS[L_col];
                    //ENTRY_1DARRAY &xx = X[L_col];


                      if(thread.team_rank() ==0)
                      { 
                        #ifdef BASKER_DEBUG_EXTEND
                        printf("clear x called by kid: %d on x: %d \n",
                               kid, team_leader);
                        #endif

                        #ifdef MOVE_A_SMALL

                        #else
                        clear_x<Int,Entry,Exe_Space>(x, L[0][0].nrow);
                        //clear_x<Int,Entry,Exe_Space>(x, A.nrow);
                        #endif
                      }
                  

                    #ifdef BASKER_DEBUG_EXTEND
                    UA.info();
                    UA.base->info();
                    #endif

                     //only fill if not yest called
                    LL.fill();
                    LU.fill(); 

                    populate_row_view<Int, Entry, Exe_Space>
                      (UA, k, 0);
                   
                    #ifdef BASKER_DEBUG_EXTEND
                    //LL.print();
                    //LU.print();
                    #endif

                    Int scol = LL.scol;
                    Int L_nnz_already = LL.col_ptr[k - scol];
                    Int U_nnz_already = LU.col_ptr[k - scol];
                    
                    Int malloc_option = 0;
                    int einfo = 0;

                    #ifdef BASKER_DEBUG_EXTEND
                    printf("HERE col_fact, scol: %d lnnz: %d unnz: %d\n", 
                           scol, L_nnz_already, U_nnz_already); 
                    UA.base->print();
                    #endif


                    //Init workspace for L.... needs to be move to symbolic!!!
  

                    //Note: I don't think we want to do this here!!!!
                    // Remthink
                    //init_workspace<Int, Entry, Exe_Space>(L_col,LL,ws,x);
                    init_workspace<Int,Entry,Exe_Space>(L_col,ws,x);


                    //if(kid == 0)
                    {
                      //Note:changes ((Int) 0 -> kid)
                      col_factor<Int,Entry,Exe_Space>
			( kid,  k, L_nnz_already, U_nnz_already,  
                                  scol, LL.srow,
                                  ws,
                                  x, gperm[0],
                                  UA,
                                  LL,
                                  LU,
                                  malloc_option,
                                  stats,
                                  &einfo);
                      
                      //printf("Zeroing column matrix: %d %d \n", kid, 0);
                      zero_col_matrix<Int,Entry,Exe_Space>(C[kid/2][0]);
                    }

                  }//end if() should factor??
                  
              }//test nothing
                                     
              //maybe get rid of this one later        
              thread.team_barrier(); //EXP
           
               #ifdef BASKER_TIME_EXTEND
              stat_time[2] += timer_factor.seconds();
              //printf("kid: %d lvl: %d time_factor: %f \n",
              //       kid, l, timer_factor.seconds());
              #endif
             
            }//for each sublevel

          //thread.team_barrier();
          //Now test if I am the thread to do the factor
          //thread.team_rank() == 0 is king
          
          #ifdef BASKER_DEBUG_EXTEND
          printf("\n-----Done with column: %d   KID %d -------\n",
                 k, kid);
          #endif

          
          #ifdef BASKER_TIME_EXTEND
          Kokkos::Impl::Timer timer_rollback;
          #endif

          if(roll_back)
            {
              #ifdef BASKER_DEBUG_EXTEND
              printf("Roll back called  KID: %d \n", kid);
              #endif

              //NOTE:FIX!!!!!!
              for(Int rv = 0 ; rv < A[U_col].size(); rv++)
                {
                  
                  A[U_col][rv].base = base;
                  A[U_col][rv].k_offset = 0;
                }
            }//end if() roll_back
       
          /////-------------NEEDED BARRIER----------///
          thread.team_barrier(); //EXP
        
          //-------------CHANGE OF COLUMNS----------.//
          #ifdef BASKER_TIMER_EXTEND
          stat_time[3] += timer_rollback.seconds();
          //printf("kid: %d lvl: %d time_rollback: %f \n",
          //       kid, l, timer_rollback.seconds());
          #endif

        }//end for() each column

        #ifdef BASKER_DEBUG_EXTEND
        //print_factor(U[U_col][0], U[U_col][1]);
        #endif

      }//end temp test


      #ifdef BASKER_TIME_DETAIL
      stats.nfactor_sep_time[kid] += timer_sep.seconds();
      #endif

      
      #ifdef BASKER_TIME_EXTEND
      printf("kid: %d Total-Time  Tri-Solve: %f  Reduce: %f Factor: %f Roll-Back: %f \n",
             kid, stat_time[0], stat_time[1], stat_time[2], stat_time[3]);
      //printf("kid: %d Tri-Solve %f %f %f \n",
      //     kid, stat_tri_time[0], stat_tri_time[1], stat_tri_time[2]);

       #endif
      
    }//end operator

  }; //end rect_factor_funct
} //end namespace Basker

#endif //end ifndef


