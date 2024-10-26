// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_DIAG_HPP
#define SHYLUBASKER_NFACTOR_DIAG_HPP

#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"
//#include "Teuchos_LAPACK.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#endif 

//#define BASKER_TIMER
//#define BASKER_DEBUG_NFACTOR_DIAG

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_diag
  {
    //Comeback and cleanup kokkos
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    Basker<Int,Entry,Exe_Space> *basker;

    kokkos_nfactor_diag()
    {}

    kokkos_nfactor_diag(Basker<Int,Entry,Exe_Space> *_basker)
    {
      basker = _basker;
    }

    //comeback and fix kokkos stuff

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {
      Int kid = (Int)(thread.league_rank()*
          thread.team_size()+
          thread.team_rank());
      //Int total_threads = thread.league_size()*thread.team_size(); //Not used


      //================OLD EQ BLK=================//

      //divide equally and launch
      //Int nchunks = (basker->btf_nblks - 
      //               basker->btf_tabs_offset);
      //Int chunk_size = nchunks/total_threads;
      //Int chunks_left_over = nchunks - 
      //                       (chunk_size*total_threads);
      //Int chunk_start = kid*chunk_size;

      //Divid into chunks based on work

      //printf("kid: %d c_size: %d numc: %d schunk: %d \n", kid, chunk_size, nchunks, chunk_start);

      //basker->t_nfactor_diag(kid, chunk_start, 
      //                       chunk_size);

      //extra
      //if(kid == 0)
      //{
      //  Int extra_start=chunk_size*total_threads;
      //  basker->t_nfactor_diag(kid, extra_start,
      //                         chunks_left_over);
      //}


      //============= NEW EQ  =============//

      Int chunk_start = basker->btf_schedule(kid);
      Int chunk_size  = basker->btf_schedule(kid+1) - basker->btf_schedule(kid);


      //printf(" > kokkos_nfactor_diag: Chunk start = %d, size = %d\n", chunk_start, chunk_size);
      int info = BASKER_SUCCESS;
      if(chunk_size > 0)
      {
        info = basker->t_nfactor_diag(kid, chunk_start, chunk_size);
      }

#if defined(BASKER_SPLIT_A) 
      if (info == BASKER_SUCCESS && basker->btf_top_nblks > 0) {
        Int tot_num_ranks = (Int)(thread.league_size()*thread.team_size());
        chunk_size = basker->btf_top_nblks / tot_num_ranks;
        chunk_start = kid * chunk_size;
        if (kid == tot_num_ranks-1) {
          chunk_size = basker->btf_top_nblks - chunk_start;
        }

        //printf( " Top D factor:t_nfactor_diag(kid=%d: chunk_start=%d, chunk_size=%d) with tab=%d top_tab=%d\n",kid,chunk_start,chunk_size,
        //        basker->btf_tabs_offset,basker->btf_top_tabs_offset );
        basker->t_nfactor_diag(kid, chunk_start, chunk_size);
      }
#endif

    }//end operator()
    
  };//end kokkos_nfactor_diag


  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_diag_remalloc
  {
    //Comeback and cleanup kokkos
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    Basker<Int,Entry,Exe_Space> *basker;
    
    INT_1DARRAY thread_start_top;
    INT_1DARRAY thread_start;

    kokkos_nfactor_diag_remalloc()
    {}

    kokkos_nfactor_diag_remalloc
    (
     Basker<Int,Entry,Exe_Space> *_basker,
     INT_1DARRAY                 _thread_start_top,
     INT_1DARRAY                 _thread_start
     )
    {
      basker           = _basker;
      thread_start_top = _thread_start_top;
      thread_start     = _thread_start;
    }

    //comeback and fix kokkos stuff

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {
      Int kid = (Int)(thread.league_rank()*
          thread.team_size()+
          thread.team_rank());
      //Int total_threads = thread.league_size()*
      //thread.team_size(); //Not used

      int info = BASKER_SUCCESS;
      Int chunk_start = thread_start(kid);
      //printf(" > kokkos_nfactor_diag_remalloc: %d %d: start = %d (max=%d) \n", kid, thread_start(kid), chunk_start, BASKER_MAX_IDX);
      if(chunk_start != BASKER_MAX_IDX)
      {
        Int chunk_size  = basker->btf_schedule(kid+1) - chunk_start;

        //printf("Chunk start: %d size: %d \n", chunk_start, chunk_size);
        if (chunk_size > 0) {
          info = basker->t_nfactor_diag(kid, chunk_start, chunk_size);
          if (info == BASKER_SUCCESS) {
            thread_start(kid) = basker->btf_schedule(kid+1);
          }
        }
      }//end if

#if defined(BASKER_SPLIT_A) 
      if (info == BASKER_SUCCESS && basker->btf_top_nblks > 0) {
        chunk_start = thread_start_top(kid);

        if(chunk_start != BASKER_MAX_IDX) {
          Int tot_num_ranks = (Int)(thread.league_size()*thread.team_size());
          Int chunk_size = basker->btf_top_nblks / tot_num_ranks;
          if (kid == tot_num_ranks-1) {
            chunk_size = basker->btf_top_nblks - chunk_start;
          }

          basker->t_nfactor_diag(kid, chunk_start, chunk_size);
        }
      }
#endif

    }//end operator()
    
  };//end kokkos_nfactor_diag_remalloc


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_nfactor_diag
  (
   Int kid,
   Int schunk,
   Int nchunk
  )
  {
    #ifdef BASKER_TIMER
    Kokkos::Timer timer_nfactor;
    timer_nfactor.reset();
    #endif

    //for(Int c = schunk+btf_tabs_offset;
    //    c < schunk+nchunk+btf_tabs_offset; ++c)
    //printf( "\n >> t_nfactor_diag( c = %d .. %d) <<\n",schunk,schunk+nchunk-1 );
    for(Int c = schunk; c < schunk+nchunk; ++c)
    {
      Int c_size = btf_tabs(c+1)-btf_tabs(c);

      if(Options.verbose == BASKER_TRUE)
      {
        printf(" > kid: %ld factoring_diag current_chunk: %ld size: %ld start: %ld\n",
            (long)kid, (long)c, (long)c_size, (long)btf_tabs(c));
      }

      Int err = BASKER_SUCCESS;
      if(c_size == 1)
      {
        //call single
        err = t_single_nfactor(kid, c);
      }
      else
      {
        //call GP alg.
        err = t_blk_nfactor(kid, c);
      }
      if(err != BASKER_SUCCESS)
      {
        return BASKER_ERROR;
      }
    }//over all chunks 
    #ifdef BASKER_TIMER
    double nfactor_time = timer_nfactor.seconds();
    std::cout << " > Basker nfactor_diag(" << schunk << ":" << schunk+nchunk-1 << ") : time: " << nfactor_time
              << std::endl << std::endl;
    #endif
    return BASKER_SUCCESS;
  }//end t_nfactor_diag


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_single_nfactor
  (
   Int kid,
   Int c
  )
  {
    const Entry zero (0.0);
    const Entry one (1.0);

    Int btab = btf_tabs_offset;

    BASKER_MATRIX  &M = (c >= btab ? BTF_C        : BTF_D);
    BASKER_MATRIX  &U = (c >= btab ? UBTF(c-btab) : U_D(c));
    BASKER_MATRIX  &L = (c >= btab ? LBTF(c-btab) : L_D(c));

    Int k = btf_tabs(c);
    Int bcol = M.scol;

    if (Options.verbose == BASKER_TRUE) {
      printf( " t_single_nfactor(c=%d)\n", (int)c );
      if (c >= btab) {
        printf( " U = UBTF(%d) and L = LBTF(%d)\n",int(c-btab),int(c-btab) );
      } else {
        printf( " U = U_D(%d) and L = L_D(%d)\n",int(c),int(c) );
      }
    }
    //Int j = M.col_ptr(k+1-bcol)-1; // was assuming the column is sorted in the ascending order of row indexes
    Entry pivot = zero;
    for (Int j = M.col_ptr(k-bcol); j < M.col_ptr(k-bcol+1); j++) {
      //printf( "%d: %d %d %e, %d\n",k-bcol, M.row_idx(j), k, M.val(j), j);
      if (M.row_idx(j) == k-bcol) pivot = M.val(j);
    }
    //Int j = M.row_idx[i];
    //will have to make a c and c'
    //printf("kid: %d chunk: %d col: %d bcol: %d j: %d\n",
    //       kid, c, k, bcol, j);
    //printf("Single blk slv, kid: %d val:%f idx:%d %d \n",
    //       kid, M.val[j], M.row_idx[j], M.srow);

    if(pivot == zero || pivot != pivot)
    {
      if (Options.verbose == BASKER_TRUE) {
        if(pivot == zero) {
          printf("Error: zero diag in single factor\n");
        } else {
          printf("Error: NaN diag in single factor\n");
        }
      }
      thread_array(kid).error_type = BASKER_ERROR_SINGULAR;
      thread_array(kid).error_blk  = c;
      thread_array(kid).error_info = k;
      return BASKER_ERROR;
    }

    U.val(0)     = pivot;
    //M already has local idxing
    //U.row_idx(0) = M.row_idx(j);
    U.row_idx(0) = 0;
    L.val(0)     = one;
    //L.row_idx(0) = M.row_idx(j);
    L.row_idx(0) = 0;
    U.col_ptr(0) = 0;
    U.col_ptr(1) = 1;
    L.col_ptr(0) = 0;
    L.col_ptr(1) = 1;

    gperm(k) = k;
    gpermi(k)= k;

    return BASKER_SUCCESS;
  }//end t_single_factor


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_blk_nfactor
  (
   Int kid, 
   Int c
  )
  {
    using STS = Teuchos::ScalarTraits<Entry>;
    using Mag = typename STS::magnitudeType;

    const Entry zero (0.0);
    const Entry one (1.0);

    Int btab = btf_tabs_offset;

    BASKER_MATRIX  &M = (c >= btab ? BTF_C : BTF_D);
    BASKER_MATRIX  &U = (c >= btab ? UBTF(c-btab) : U_D(c));
    BASKER_MATRIX  &L = (c >= btab ? LBTF(c-btab) : L_D(c));

    Int bcol = M.scol;
    //JDB: brow hack: fix.
    Int brow2 = L.srow;
    if (c >= btab) {
      brow2 -= btf_tabs(btab);
    }

    #ifdef BASKER_DEBUG_NFACTOR_DIAG
    if (Options.verbose == BASKER_TRUE) {
      printf("CURRENT BLK: %ld \n", (long)c);
      printf("btab: %ld %ld\n", (long)btab, (long)c-btab);
      printf("brow2: %ld \n", (long)brow2);
      printf("col_start: %ld \n", (long)btf_tabs(c));
      printf("BLK size %ld \n",
          (long)btf_tabs(c+1)-(long)btf_tabs(c));
    }
    #endif
    #ifdef BASKER_TIMER
    double initi_time = 0.0;
    double pivot_time = 0.0;
    double solve_time = 0.0;
    double scale_time = 0.0;
    double prune_time = 0.0;
    Kokkos::Timer timer_nfactor;
    Kokkos::Timer timer_nfactor_tot;
    timer_nfactor_tot.reset();
    #endif

    //Teuchos::LAPACK<int, Entry> lapack;
    //Mag rmin_ = lapack.LAMCH('E');
    //Mag rmin_ = lapack.LAMCH('U');
    Mag rmin_ (0.0);

    //workspace
    Int ws_size       = thread_array(kid).iws_size;
    INT_1DARRAY    ws = thread_array(kid).iws;
    ENTRY_1DARRAY   X = thread_array(kid).ews;

    Int *color    = &(ws(0));
    Int *pattern  = &(color[ws_size]);

    //local variable 
    Int top = ws_size;
    Int xnnz, lnnz, unnz, lcnt, ucnt;

    Entry pivot, value;
    //    Entry absv, maxc; //NDE - warning: unsed maxc
    // Max number of nnz allowed/allocated
    Int llnnz  = L.mnnz;
    Int uunnz  = U.mnnz;
    Mag absv = Mag (0.0);
    Mag maxv = Mag (0.0);
    Mag digv = Mag (0.0);

    Int maxindex = 0;
    Int digindex = 0;

    Int i,j, t;
    Int cu_ltop = 0;
    Int cu_utop = 0;

    Int newsize;

    unnz = 0;
    lnnz = 0;

    //#define BASKER_DEBUG_NFACTOR_DIAG
    #ifdef BASKER_DEBUG_NFACTOR_DIAG
    if (Options.verbose == BASKER_TRUE && c == 4) {
      printf( " kid=%d: c: %d bcol=%d, brow2=%d, wsize=%d\n", (int)kid, (int)c, (int)bcol, (int)brow2, (int)ws_size );
      //if (c >= btab) {
      //  printf( " BTF_C -> UBTF(%d) & LBTF(%d) \n",(int)(c-btab),(int)(c-btab) );
      //} else {
      //  printf( " BTF_D -> U_D(%d) & L_D(%d)\n",(int)c, (int)c );
      //}
      //std::cout << " K" << c << " = [" << std::endl;

      char filename[250];
      sprintf(filename,"A_%d.dat",(int)c);
      FILE *fp = fopen(filename,"w");
      for(Int k = btf_tabs(c); k < btf_tabs(c+1); ++k) {
        for( i = M.col_ptr(k-bcol); i < M.col_ptr(k-bcol+1); ++i) {
          if (M.row_idx(i) >= brow2) 
            fprintf(fp, " %d %d %d %e, %d %d %d\n", (int)i, (int)(M.row_idx(i)-brow2), (int)(k-btf_tabs(c)),
                    std::real(M.val(i)), (int)(M.row_idx(i)), (int)(k), (int)(i));
        }
      }
      fclose(fp);

      //std::cout << "];" << std::endl << std::endl << std::flush;
    }
    #endif

    // compute one-norm of diagonal block, if needed
    Mag normA_blk = Mag (0.0);
    if (Options.replace_zero_pivot) {
      for (j = btf_tabs(c); j < btf_tabs(c+1); ++j) {
        Mag anorm_j = Mag(0.0);
        for (i = M.col_ptr(j-bcol); i < M.col_ptr(j-bcol+1); ++i) {
          if (M.row_idx(i) >= brow2) {
            anorm_j += abs(M.val(i));
          }
        }
        if (anorm_j > normA_blk) {
          normA_blk = anorm_j;
        }
      }
    }

    // initialize perm vector
    // gperm_array(i) = k means that the current i-th row is the k-th row before pivot
    #ifdef BASKER_CHECK_WITH_DIAG_AFTER_PIVOT
    for(Int k = btf_tabs(c); k < btf_tabs(c+1); ++k)
    {
      gperm_array(k) = k;
      gpermi_array(k) = k;
    }
    #endif

    //for each column
    for(Int k = btf_tabs(c); k < btf_tabs(c+1); ++k)
    {
      #ifdef BASKER_DEBUG_NFACTOR_DIAG
      if (Options.verbose == BASKER_TRUE) {
        printf("\n------------ K=%d (%d) -------------\n", (int)(k-btf_tabs(c)), (int)k);
        BASKER_ASSERT(top == ws_size, "nfactor dig, top");
        for( i = 0; i < ws_size; ++i)
        {
          if(X(i) != zero)
          {
            printf("x(i) error: %d  %e k: %d \n", (int)i , std::real(X(i)), (int)k);
          }
          if(ws(i) != 0)
          {
            printf("ws(i) error: %d K: %d \n", (int)i, (int)k);
          }
          BASKER_ASSERT(X(i)  == zero, "X!=0");
          BASKER_ASSERT(ws(i) == 0, "ws!=0");
        }
      }
      #endif

      value = zero;
      pivot = zero;
      lcnt  = 0;
      ucnt  = 0;

      //for each nnz in column
      #ifdef BASKER_TIMER
      timer_nfactor.reset();
      #endif
      #ifdef BASKER_DEBUG_NFACTOR_DIAG
      if (Options.verbose == BASKER_TRUE) {
        printf(" ptr[%d] = %d:%d", (int)(k-bcol), (int)M.col_ptr(k-bcol), (int)M.col_ptr(k-bcol+1));
      }
      #endif
      for( i = M.col_ptr(k-bcol); i < M.col_ptr(k-bcol+1); ++i)
      {
        j = M.row_idx(i);
        j = j - brow2;
        #ifdef BASKER_DEBUG_NFACTOR_DIAG
        if (Options.verbose == BASKER_TRUE) {
          printf("\n Input(i=%d): %d(%d-%d) %d(%d-%d) %e", (int)i, (int)j, (int)M.row_idx(i), (int)brow2, (int)(k-bcol), (int)k, (int)bcol, std::real(M.val(i)));
          if (j < 0) {
            printf( "\n" );
          } else {
            printf(" (color[%d]=%d, gperm[%d+%d]=%d)\n", (int)j, (int)color[j], (int)j, (int)L.srow, (int)gperm(j+L.srow));
          }
        }
        #endif
        if(j < 0)
        {
          // j is in off-diagonal
          continue;
        }

        if (M.val(i) != zero ||
           (i+1 == M.col_ptr(k-bcol+1) && top == ws_size)) // the last element, and have not found non-zero entry, then go ahead and process this zero entry 
        {
          X(j) = M.val(i);
          if(color[j] == 0)
          {
            if(gperm(j+L.srow) != BASKER_MAX_IDX)
            { // j is in upper
              t_local_reach_btf(kid, L, j, top, c, k);
            }
            else
            { // j is in lower
              t_local_reach_short_btf(kid, j, top);
            }
          }
        }
      }//end over all nnz in column
      xnnz = ws_size - top;
      #ifdef BASKER_TIMER
      initi_time += timer_nfactor.seconds();
      #endif

      #ifdef BASKER_DEBUG_NFACTOR_DIAG
      if (Options.verbose == BASKER_TRUE) {
        printf("xnnz: %d ws_size: %d top: %d \n", 
               (int)xnnz, (int)ws_size, (int)top);
      }
      #endif

      //add custom local_solve
      #ifdef BASKER_TIMER
      timer_nfactor.reset();
      #endif
      t_back_solve(kid, c, L, 0, 0, k, top, xnnz);
      #ifdef BASKER_TIMER
      solve_time += timer_nfactor.seconds();
      #endif

      //Future add
      //t_locate_pivot(kid, top)  
      //find pivot
      maxv = abs(zero);
      digv = abs(zero);
      maxindex = BASKER_MAX_IDX;
      digindex = BASKER_MAX_IDX;
      #ifdef BASKER_TIMER
      timer_nfactor.reset();
      #endif
      for(i = top; i < ws_size; i++)
      {
        j = pattern[i];
        t = gperm(j+L.srow);

        value = X(j);
        if (value != value) {
          // NaN
          if (Options.verbose == BASKER_TRUE)
          {
            cout << endl;
            cout << "---------------------------" << endl;
            cout << "Error: NaN found blk: "
              << c 
              << " Column: "
              << k << std::endl;
          }
          thread_array(kid).error_type = BASKER_ERROR_NAN;
          thread_array(kid).error_blk = c;
          thread_array(kid).error_info = k;
          return BASKER_ERROR;
        }
        absv = abs(value);

        #ifdef BASKER_DEBUG_NFACTOR_DIAG
        if (Options.verbose == BASKER_TRUE) {
          #ifdef BASKER_CHECK_WITH_DIAG_AFTER_PIVOT
          printf("\nconsider X(%d)=%e -> %e, c=%d, i=%d, k=%d->%d: j=%d t=%d: perm=%d, permi=%d: val=%e,max=%e (off=%d, %d, %d)\n",
                 (int)j,X(j),absv, (int)c, (int)i, (int)k,(int)(k-L.scol), (int)j, (int)t,
                 (int)gperm_array(j+L.srow),(int)gpermi_array(j+L.srow), value,maxv, (int)L.srow,(int)L.scol,(int)btf_tabs(c));
          #else
          printf("\nconsider X(%d)=%e -> %e, c=%d, i=%d, k=%d->%d: j=%d t=%d: val=%e,max=%e (off=%d, %d, %d)\n",
                 (int)j,std::real(X(j)),absv, (int)c, (int)i, (int)k,(int)(k-L.scol), (int)j, (int)t,
                 std::real(value),maxv, (int)L.srow,(int)L.scol,(int)btf_tabs(c));
          #endif
        }
        #endif

        if(t == BASKER_MAX_IDX)
        {
          ++lcnt;
          // looking for the pivot candidate (largest off-diagonal)
          // will be compared with diagonal
          if(absv > maxv || maxindex == BASKER_MAX_IDX) 
          {
            maxv     = absv;
            pivot    = value;
            maxindex = j;
            #ifdef BASKER_DEBUG_NFACTOR_DIAG
            if (Options.verbose == BASKER_TRUE) {
              printf( " pivot=%e (k=%d, j=%d -> %d)\n", std::real(pivot), (int)k, (int)j, (int)(j+L.srow));
            }
            #endif
          }
          #ifdef BASKER_CHECK_WITH_DIAG_AFTER_PIVOT
          if (gpermi_array(j+L.srow) == k)
          #else
          if (j+L.srow == k)
          #endif
          {
            digv = absv;
            digindex = j;
            #ifdef BASKER_DEBUG_NFACTOR_DIAG
            if (Options.verbose == BASKER_TRUE) {
              #ifdef BASKER_CHECK_WITH_DIAG_AFTER_PIVOT
              printf( " digv=%e (k=%d, permii(%d) = %d, brow=%d,bcol=%d)\n",absv,k,j+L.srow,gpermi_array(j+L.srow),M.scol,M.srow);
              #else
              printf( " digv=%e (k=%d, brow=%d,bcol=%d)\n", absv, (int)(k), (int)(M.scol), (int)(M.srow));
              #endif
            }
            #endif
          }
        }
      } //for (i = top; i < ws_size)
      //printf("b: %d lcnt: %d after \n", b, lcnt);
      #ifdef BASKER_DEBUG_NFACTOR_DIAG
      if (Options.verbose == BASKER_TRUE) {
        printf(" >> pivot=%e, maxindex=%d, k=%d (digv = %e, pivot_tol = %g)\n",
               std::real(pivot), (int)maxindex, (int)k, digv, Options.pivot_tol);
      }
      #endif
      //printf( "c=%d k=%d (%d) maxindex=%d pivot=%e maxv=%e, diag=%e tol=%e (nopivot=%d)\n", c, k, k-btf_tabs(c), maxindex, pivot, maxv, digv, Options.pivot_tol,Options.no_pivot);

      // check if diagonal is relatively large enough
      if(digindex == BASKER_MAX_IDX) {
        #if 0 // diagonal may be zero (and not in CSC) in the original matrix
        if (Options.verbose == BASKER_TRUE)
        {
          cout << "----------------------------" <<endl;
          cout << "  Failed to find diagonal" << std::endl;
          cout << "----------------------------" <<endl;
        }
        #endif
      } else if(Options.no_pivot == BASKER_TRUE || digv > maxv * Options.pivot_tol)
      {
        maxindex = digindex;
        pivot    = X(maxindex);
        #ifdef BASKER_DEBUG_NFACTOR_DIAG
        if (Options.verbose == BASKER_TRUE) {
          printf( " pivot -> %d %e (using diagonal)\n",maxindex,pivot );
        }
        #endif
      }

      #ifdef BASKER_DEBUG_NFACTOR_DEBUG
      if (Options.verbose == BASKER_TRUE) {
        printf("pivot: %g maxindex: %d K:%d \n", pivot, maxindex, k);
      }
      #endif

      ucnt = ws_size - top - lcnt +1;

      if((maxindex == BASKER_MAX_IDX) || (pivot == zero))
      {
        const Mag eps = STS::eps ();
        const Mag normA = A.gnorm;
        if (Options.verbose == BASKER_TRUE)
        {
          cout << endl;
          cout << "---------------------------"
            << endl;
          cout << "Error: Diag Matrix is singular, blk: "
            << c 
            << " Column: "
            << k
            << " Size: "
            << btf_tabs(c+1)-btf_tabs(c)
            << endl;
          cout << "MaxIndex: " << maxindex 
            << " pivot " 
            << pivot << endl;
          if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
            cout << "  replace zero pivot with " << normA_blk * eps
                 << " (normA = " << normA << ", normA_blk = " << normA_blk << ", eps = " << eps << ")" << endl;
          } else {
            if (!Options.replace_zero_pivot) {
              cout << " replace_zero_pivot disabled, ";
            }
            if (normA_blk <= abs(zero)) {
              cout << " empty block, ";
            }
            if (maxindex == BASKER_MAX_IDX) {
              cout << " empty column, ";
            }
            cout << endl;
          }
        }
        if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          pivot = normA_blk * eps;
          X(maxindex) = pivot;
        } else {
          thread_array(kid).error_type = BASKER_ERROR_SINGULAR;
          thread_array(kid).error_blk  = c;
          thread_array(kid).error_info = k;
          return BASKER_ERROR;
        }
      }

      //printf("----------------PIVOT------------blk: %d %d \n", 
      //      c, btf_tabs(c+1)-btf_tabs(c));
      // store pivot
      gperm(maxindex+L.scol) = k;
      gpermi(k)              = maxindex + L.srow;
      //if (c == 4) std::cout << " k = " << k << " : gperm( " << maxindex << " + " << L.scol << " = " << maxindex+L.scol << " ) = " << k << "(" << k-btf_tabs(c) << ")" 
      //                      << (k == maxindex+L.scol ? "" : " pivot") << std::endl;
      #ifdef BASKER_TIMER
      pivot_time += timer_nfactor.seconds();
      #endif
      // > maxindex is in the original row (before pivot)
      #ifdef BASKER_CHECK_WITH_DIAG_AFTER_PIVOT
      int pivot_index = gpermi_array(maxindex+L.srow);
      if (k != pivot_index) {
        // update global perm vector for figuring out diagonal entry
        //
        // swap perm
        int pivot_row = gperm_array(k);
        gperm_array(k) = gperm_array(pivot_index);
        gperm_array(pivot_index) = pivot_row;

        // swap iperm
        int row1 = gperm_array(k);
        int row2 = gperm_array(pivot_index);

        pivot_row = gpermi_array(row1);
        gpermi_array(row1) = gpermi_array(row2);
        gpermi_array(row2) = pivot_row;
        //printf( " > swap(%d, %d)\n",row1,row2 );
      }
      #endif
      #ifdef BASKER_DEBUG_NFACTOR_DIAG
      if (Options.verbose == BASKER_TRUE) {
        printf( " > gperm(%d+%d = %d) = %d\n", (int)maxindex, (int)L.scol, (int)(maxindex+L.scol), (int)k );
        //for(Int ii = btf_tabs(c); ii < btf_tabs(c+1); ++ii) {
        //  printf( "gperm_array(%d) = %d, gpermi_array(%d) = %d\n",ii,gperm_array(ii), ii,gpermi_array(ii));
        //}
        if((maxindex+L.scol) != k)
        {
          cout << " >> Permuting Pivot: " << k
               << " as row " 
               << maxindex+L.scol << endl;
        }
      }
      #endif

      //Note: Come back to this!!!!
      if(lnnz + lcnt > llnnz)
      {
        newsize = lnnz * 1.1 + 2 *L.nrow + 1;

        if (Options.verbose == BASKER_TRUE)
        {
          printf("Diag blk: %ld Reallocing L oldsize: %ld current: %ld count: %ld newsize: %ld (n = %ld)\n",
              (long)c, (long)llnnz, (long)lnnz, (long)lcnt, (long)newsize, (long)btf_tabs(c+1)-btf_tabs(c));
          printf("Columns in blks: %ld %ld %ld \n",
              (long)btf_tabs(c), (long)btf_tabs(c+1), (long)(btf_tabs(c+1)-btf_tabs(c)));
        }

        thread_array(kid).error_blk = c;
        if(Options.realloc == BASKER_FALSE)
        {
          thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
          return BASKER_ERROR;
        }
        else
        {
          thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
          thread_array(kid).error_info = newsize;
          return BASKER_ERROR;
        }
      }
      if(unnz+ucnt > uunnz)
      {
        newsize = uunnz*1.1 + 2*L.nrow+1;

        if(Options.verbose == BASKER_TRUE)
        {
          printf("Diag blk: %ld Reallocing U oldsize: %ld %ld newsize: %ld \n",
              (long)c, (long)uunnz, (long)unnz+ucnt, (long)newsize);
          printf("blk: %ld column: %ld \n", (long)c, (long)k);
        }

        thread_array(kid).error_blk = c;
        if(Options.realloc == BASKER_FALSE)
        {
          thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
          return BASKER_ERROR;
        }
        else
        {
          thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
          thread_array(kid).error_info = newsize;
          return BASKER_ERROR;
        }
      }

      #ifdef BASKER_TIMER
      timer_nfactor.reset();
      #endif
      L.row_idx(lnnz) = maxindex; // ???
      L.val(lnnz)     = one;
      ++lnnz;

      //printf("Pivot: %f \n", pivot);

      Entry lastU = zero;
      for( i = top; i < ws_size; i++)
      {
        j = pattern[i];
        t = gperm(j+L.srow);

        #ifdef BASKER_DEBUG_NFACTOR_DIAG
        if (Options.verbose == BASKER_TRUE) {
          printf("j = %d, t = %d, k = %d, x = %g\n",
                  (int)j, (int)t, (int)k, std::real(X(j)));
        }
        #endif            

        if(t != BASKER_MAX_IDX)
        { // U
          if (Options.prune || abs(X(j)) > rmin_)
          {
            if(t < k)
            {
              U.row_idx(unnz) = t - L.srow;
              #ifdef BASKER_2DL
              U.val(unnz) = X(j);
              #else
              U.val[unnz] = X[j];
              #endif

              #ifdef BASKER_DEBUG_NFACTOR_DIAG
              if (Options.verbose == BASKER_TRUE) {
                printf("%d: U(%d,%d) = %f \n", (int)unnz, (int)U.row_idx(unnz), (int)(k-L.srow), std::real(U.val(unnz)));
              }
              #endif
              ++unnz;
            }
            else
            {
              #ifdef BASKER_2DL
              lastU = X(j);
              #else
              lastU = X[j];
              #endif
            }
          }
        }
        else if (t == BASKER_MAX_IDX)
        { // L
          if (Options.prune || abs(X(j)) > rmin_*abs(pivot))
          {
            L.row_idx(lnnz) = j;
            #ifdef BASKER_2DL
            L.val(lnnz) = X(j)/pivot;
            #else
            L.val[lnnz] = X[j]/pivot;
            #endif

            #ifdef BASKER_DEBUG_NFACTOR_DIAG
            if (Options.verbose == BASKER_TRUE) {
              printf("%d: L(%d,%d) = %f \n", (int)lnnz, (int)L.row_idx(lnnz), (int)(k-L.srow), std::real(L.val(lnnz)));
            }
            #endif
            ++lnnz;
          }
        }//end if() not 0
        else
        {
          BASKER_ASSERT(0==1, "T not be selected\n");
        }
        //Note: move x[j] inside of if() not 0..
        //..extra ops this way
        #ifdef BASKER_DEBUG_NFACTOR_DIAG
        if (Options.verbose == BASKER_TRUE) {
          printf("Zeroing element: %d \n", (int)j);
        }
        #endif

        #ifdef BASKER_2DL
        X(j) = zero;
        #else
        X[j] = zero;
        #endif
      }
      #ifdef BASKER_TIMER
      scale_time += timer_nfactor.seconds();
      #endif

      //Fill in last element of U
      U.row_idx(unnz) = k - L.scol;
      U.val(unnz)     = lastU; //NDE: lastU init to 0 line 563; only updated if 't < k' and t not BASKER_MAX_IDX; t defined by gperm i.e. pivoting
      if(lastU == zero)
      {
        printf("Basker t_blk_nfactor: diag btf zero, error, c: %ld k: %ld \n", (long)c, (long)k);
      }
      ++unnz;

      xnnz = 0;
      top = ws_size;

      L.col_ptr(k-L.srow) = cu_ltop;
      L.col_ptr(k+1-L.srow) = lnnz;
      cu_ltop = lnnz;

      U.col_ptr(k-L.srow) = cu_utop;
      U.col_ptr(k+1-L.srow) = unnz;
      cu_utop = unnz;

      if(Options.prune && L.ncol > Options.btf_prune_size)
      {
        #ifdef BASKER_TIMER
        timer_nfactor.reset();
        #endif
        t_prune_btf(kid, L, U, k-L.srow, maxindex);
        #ifdef BASKER_TIMER
        prune_time += timer_nfactor.seconds();
        #endif
      }

    }//over each column
    L.nnz = lnnz;
    U.nnz = unnz;

    L.col_ptr(L.ncol) = lnnz;
    U.col_ptr(U.ncol) = unnz;

    #if 0
    printf("L=[\n");
    for (int j = 0; j < L.ncol; j++) {
      for (int k = L.col_ptr[j]; k < L.col_ptr[j+1]; k++) {
        printf( "%d %d %e\n",(int)L.row_idx[k],j,L.val[k]);
      }
    }
    printf("];\n");
    printf("U=[\n");
    for (int j = 0; j < U.ncol; j++) {
      for (int k = U.col_ptr[j]; k < U.col_ptr[j+1]; k++) {
        printf( "%d %d %e\n",(int)U.row_idx[k],j,U.val[k]);
      }
    }
    printf("];\n");
    #endif

    #ifdef BASKER_TIMER
    double total_time = timer_nfactor_tot.seconds();
    std::cout << " ++ Basker nfactor_diag(" << c << ") : n=" << btf_tabs(c+1)-btf_tabs(c)
              << " nnz(L)=" << lnnz << " nnz(U)= " << unnz << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : init   time: " << initi_time << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : update time: " << solve_time << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : pivot  time: " << pivot_time << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : scale  time: " << scale_time << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : prune  time: " << prune_time << std::endl;
    std::cout << " ++ Basker nfactor_diag(" << c << ") : total  time: " << total_time << std::endl << std::endl;
    #endif
    return 0;
  }//end t_factor_diag()


  template <class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  void Basker<Int,Entry,Exe_Space>::t_local_reach_short_btf
  (
   const Int kid,
   const Int j,
   Int &top
  )
  {
    //printf("=======LOCAL REACH BTF SHORT CALLED (pattern[top=%d - 1] = %d) =====\n",(int)top, (int)j);
    INT_1DARRAY    ws  = thread_array(kid).iws;
    Int        ws_size = thread_array(kid).iws_size;

    Int *color   = &(ws(0));
    Int *pattern = &(ws(ws_size));

    color[j]       = 2;
    pattern[--top] = j;
  }//end t_locak_reach_short_btf

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach_btf
  (
   Int kid, 
   BASKER_MATRIX &L,
   Int root_j, 
   Int &top,
   Int c, Int k
  )
  {
    //printf("=======LOCAL REACH BTF CALLED =====\n");

    INT_1DARRAY    ws  = thread_array(kid).iws;
    Int        ws_size = thread_array(kid).iws_size;
 
    /*{
      printf("ws_size: %d \n", ws_size);
      printf("total size: %d \n", thread_array(kid).iws_size*thread_array(kid).iws_mult);
    }*/
 
    Int brow        = L.srow;
    Int *pattern    = &(ws(ws_size));
    Int *stack      = &(pattern[ws_size]);
    Int *store      = &(stack[ws_size]);

    Int start = -1;
    Int head  = 0;
    stack[0]  = root_j;

    //printf( " t_reach (L.srow = %d, L.scol = %d)\n",L.srow,L.scol );
    while (head != BASKER_MAX_IDX) // BASKER_MAX_IDX = -1
    { 
      Int j = stack[head];
      Int t = gperm(j+brow);

    #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("%d: stack[%d]=%d, gmer[%d]=%d, head=%d \n", kid, head,j ,j+brow,gperm(j+brow), head);
      //BASKER_ASSERT(head > -1,"local head val\n");
      //printf("----------DFS: %d %d -------------\n", j, t);
      printf( "%d: x ws(%d) = %d, gperm(%d + %d) = %d -> %d\n",kid, j,ws(j), j,brow,t,t-L.scol );
    #endif

      if(ws(j) == 0)
      {
        //Color
        ws(j) = 1;

        if (t != BASKER_MAX_IDX)
        {
          //Backward walk
          //start = L.col_ptr(t+1-L.scol);
          //Add pend here
          start = 
            (L.pend(t-L.scol)==BASKER_MAX_IDX) ? L.col_ptr(t+1-L.scol) : L.pend(t-L.scol);
        #ifdef BASKER_DEBUG_LOCAL_REACH
          printf("%d: reach.. j: %d t:%d L.scol %d\n", kid, j, t, L.scol);
          printf("%d: start = (pend(%d) = %d ? ptr(%d) : ptr(%d)) = %d\n", kid, t-L.scol,L.pend(t-L.scol), t+1-L.scol,t-L.scol, start);
          //printf("dfs. col: %d prune: %d \n", L.col_ptr(t+1-L.scol), L.pend(t-L.scol));
        #endif
        }
        else
        {
        #ifdef BASKER_DEBUG_LOCAL_REACH
          printf("%d: L.scol: %d L.ncol: %d t: %d \n", kid, L.scol, L.ncol,t);
        #endif
        }
      }
      else
      {
        #ifdef BASKER_DEBUG_LOCAL_REACH
          printf("%d: start = store[%d] = %d\n", kid, j,store[j]);
        #endif
        start = store[j];
      }

      Int i1;
      Int end = L.col_ptr(t-L.scol);
      #ifdef BASKER_DEBUG_LOCAL_REACH
        printf( "%d: > start = %d, end = %d (t = %d)\n",kid, start,end,t );
      #endif
      for(i1 = --start; i1 >= end; --i1)
      {
        Int i = L.row_idx(i1);

        #ifdef BASKER_DEBUG_LOCAL_REACH
          printf( "%d: > row_idx(%d) = %d\n",kid, i1,i );
        #endif
        if(ws(i) != 0)
        {
          continue;
        }
        else
        {
          if (gperm(i+brow) >= 0) // previous col -> belongs to upper-tri
          { // i is in upper
            #ifdef BASKER_DEBUG_LOCAL_REACH
              printf( "%d: > store[%d] = %d, stack[%d] = %d\n",kid, j,i1, head,i );
            #endif
            store[j] = i1;
            stack[++head] = i;
            break;
          }
          else
          { // i is in lower
            #ifdef BASKER_DEBUG_LOCAL_REACH
              printf( "%d: > ws(%d) = 1, pattern[%d - 1] = %d\n",kid, i, top,i );
            #endif
            ws(i) = 1;
            pattern[--top] = i;
          } //end if(gperm)
        } //end if/else ws(i)!=0
      } //end for start:-1:end

      //if we have used up the whole column
      if (i1 < end)
      {
        #ifdef BASKER_DEBUG_LOCAL_REACH
          printf( "%d: + ws(%d) = 1, pattern[%d - 1] = %d\n",kid, j, top,j );
        #endif
        head--;
        ws(j)          = 2;
        pattern[--top] = j;
      }
    } //end while 
    return 0;
  }//end t_local_reach_btf

    
  //Note: need to fix indexing for color
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach_old
  (
   Int kid, 
   BASKER_MATRIX &L,
   Int lvl, 
   Int l,
   Int j, 
   Int *top
  )
  {
    
    INT_1DARRAY    ws  = thread_array(kid).iws;
    Int        ws_size = thread_array(kid).iws_size;
 
    /*
    printf("ws_size: %d \n", ws_size);
    printf("total size: %d \n", 
           thread_array(kid).iws_size*thread_array(kid).iws_mult);
    */
 
    Int brow        = L.srow;
    
    Int *color       = &(ws(0));
    Int *pattern     = &(ws(ws_size));
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);

    Int i, t, head, i1;
    Int start, end, done;
    
    start = -1;

    head = 0;
    
    stack[head] = j;

    while(head != BASKER_MAX_IDX)
    { 
    #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("stack_offset: %d head: %d \n", stack_offset , head);
      BASKER_ASSERT(head > -1,"local head val\n");
    #endif

      j = stack[head];
      //t = gperm[j];
      t = gperm(j+brow);

    #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("----------DFS: %d %d -------------\n", j, t);
    #endif

      #ifdef BASKER_2DL
      //if(color[j-brow] == 0)
      if(color[j] == 0)
      #else
      if(color[j] == 0)
      #endif
        {    
        #ifdef BASKER_2DL
          //color[j-brow] = 1;
          color[j] = 1;
        #else
          color[j] = 1;
        #endif

          //if not permuted, might be able to make this smaller now
          //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
          //if(t!=L.max_idx)
          if(t!=BASKER_MAX_IDX)
          {
          #ifdef BASKER_DEBUG_LOCAL_REACH
            printf("reach.. j: %d t:%d L.scol %d\n",
                j, t, L.scol);
          #endif
            //start = L.col_ptr[t-L.scol];
            start = L.col_ptr(t-L.scol);
          }
          else
          {
          #ifdef BASKER_DEBUG_LOCAL_REACH
            printf("L.scol: %d L.ncol: %d t: %d \n",
                L.scol, L.ncol,t);
          #endif
          }
        }
        else
        {
          //start = store[j-brow];
          start = store[j];
        }
      done = 1;

      //if  permuted
      if(t!=BASKER_MAX_IDX)
      {

        //end = L.col_ptr[t+1-L.scol];
        end = L.col_ptr(t+1-L.scol);
        for(i1 = start; i1 < end; ++i1)
        {
          //i = L.row_idx[i1];
          i = L.row_idx(i1);
          #ifdef BASKER_2DL
          //if(color[i-brow] == 0)
          if(color[i] == 0)
          #else
          if(color[i] == 0)
          #endif
            {
              ++head;
              stack[head] = i;
              //store[j-brow] = i1+1;
              store[j] = i1+1;
              done = 0;

              break;
            }
        }
      }
      if(done)
      {
        //sprintf("%d \n", *top);
        pattern[--*top] = j;
        #ifdef BASKER_2DL
        //color[j-brow] = 2;
        color[j] = 2;
        #else
        color[j] = 2;
        #endif
        if(head == 0)
        {head = BASKER_MAX_IDX;}
        else
        {--head;}

      }
    }//end while 
    return 0;
  }//end t_local_reach()    
    
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::t_back_solve
  (
   Int kid,
   Int c,
   BASKER_MATRIX &L,
   Int lvl,
   Int l,
   Int k, 
   Int top,
   Int xnnz
  )
  {
    const Entry zero (0.0);
    
    INT_1DARRAY   ws = thread_array(kid).iws;
    ENTRY_1DARRAY  X = thread_array(kid).ews;
    Int      ws_size = thread_array(kid).iws_size;
    
    Int brow     = L.srow;
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
 
    //printf("back solve called. nnz: %d \n",
    //       xnnz);

    for(Int top1 = top; top1 < top+xnnz; ++top1)
    {
      Int j = pattern[top1];
      Int t = gperm(j+brow);

      #ifdef BASKER_2DL
      color[j] = 0;
      #else
      color[j] = 0;
      #endif

      #ifdef BASKER_2DL
        Entry xj = X(j);
      #else
        Entry xj = X[j];
      #endif
      if(t != BASKER_MAX_IDX && xj != zero)
      { // j is original nonzero in upper-triangullar part of A
        // p starts witth ptr(k)+1, skipping diagonal element
        Int k_i = t - L.scol;
        for(Int p = L.col_ptr(k_i)+1; p < L.col_ptr(k_i+1); ++p)
        {
          #ifdef BASKER_DEBUG_NFACTOR_DIAG
          if (Options.verbose == BASKER_TRUE) {
            #ifdef BASKER_2DL
            printf("Updating block %d row: %d with value: %f - %f * %f -> %f (j=%d, p=%d)\n",
                   (int)c, (int)L.row_idx[p],
                   std::real(X[L.row_idx[p]]), std::real(L.val[p]), std::real(xj),
                   std::real(X(L.row_idx(p)) - L.val(p)*xj), (int)j, (int)p);
            #else
            printf("Updating block %d row: %d with value: %f %f -> %f (j=%d)\n",
                   (int)c, (int)L.row_idx(p),
                   std::real(X[L.row_idx(p)]), std::real(L.val[p]*xj), std::real(X[L.row_idx[p]] - L.val[p]*xj), (int)j);
            #endif
          }
          #endif

          #ifdef BASKER_2DL
          X(L.row_idx(p)) -= L.val(p)*xj;
          #else
          X[L.row_idx[p]] -= L.val[p] *xj;
          #endif
        }//end for() over each nnz in the column
      }//end if() not permuted
    }//end for() over all nnz in LHS
    return 0;
  }//end t_back_solve()


  template <class Int, class Entry, class Exe_Space>
  inline
  void Basker<Int, Entry, Exe_Space>::t_prune_btf
  (
   const Int kid,
   const BASKER_MATRIX &L,
   const BASKER_MATRIX &U,
   const Int k,
   const Int pivotrow
  )
  {
    //printf("prune_btf called. kid: %d k: %d pr: %d \n",
    //       kid, k, pivotrow);

    for(Int ui = U.col_ptr(k); ui < U.col_ptr(k+1)-1; ++ui)
    {
      Int j = U.row_idx(ui);
      BASKER_ASSERT(j<k, "Pruning, j not less than k");

      if(L.pend(j) == BASKER_MAX_IDX)
      {

        for(Int li = L.col_ptr(j); li < L.col_ptr(j+1); ++li)
        {
          //if we can prune?
          if(pivotrow == L.row_idx(li))
          {
            //Bubble up to right position
            Int phead = L.col_ptr(j);
            Int ptail = L.col_ptr(j+1);

            while(phead < ptail)
            {
              Int i = L.row_idx(phead);
              if(gperm(i+L.srow) >= 0)
              {
                //Advance head
                phead++;
              }
              else
              {
                //flip head and tail
                ptail--;
                L.row_idx(phead) = L.row_idx(ptail);
                L.row_idx(ptail) = i;
                Entry x = L.val(phead);
                L.val(phead) = L.val(ptail);
                L.val(ptail) = x;
              }//end flip
            }//end while phead < ptail
            L.pend(j) = ptail;
            break;
          }//if pr == L.row
        }//end for over all nnz in L
      }//end if L.pend == MAXIDX
    }//end for over all nnz in U
  }//end t_prune_btf


}//end namespace BaskerNS

#undef BASKER_TIMER
#endif//end BASKER_NFACTOR_DIAG_HPP
