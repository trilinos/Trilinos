#ifndef SHYLUBASKER_NFACTOR_BLK_HPP
#define SHYLUBASKER_NFACTOR_BLK_HPP

#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"

#include <string>

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

//#define BASKER_DEBUG_NFACTOR_BLK

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space> *basker;

    kokkos_nfactor_domain()
    {}

    kokkos_nfactor_domain(Basker<Int,Entry,Exe_Space> *_basker)
    { basker = _basker;}

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const
    #endif
    {
      #ifdef BASKER_KOKKOS
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //	      thread.team_rank());
      Int kid = basker->t_get_kid(thread);
      #endif


      #ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("\n-----------BLK---Kid: %d -------------\n",
	     kid);
      #endif

      //No longer needed, added in sfactor
      //printf("before workspace init\n");
      //basker->t_init_workspace(kid);
      //printf("after workspace init\n");
      
      //if(kid == 1)
	{
      basker->t_nfactor_blk(kid);
	}
    }//end operator
  };//end kokkos_nfactor_domain struct


  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain_remalloc
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space> *basker;
    INT_1DARRAY                 thread_start;

    kokkos_nfactor_domain_remalloc()
    {}

    kokkos_nfactor_domain_remalloc
    (
     Basker<Int,Entry,Exe_Space> *_basker, 
     INT_1DARRAY                 _thread_start
     )
    {
      basker       = _basker;
      thread_start = _thread_start;
    }

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const
    #endif
    {
      #ifdef BASKER_KOKKOS
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //	      thread.team_rank());
      Int kid = basker->t_get_kid(thread);
      #endif


      if(thread_start(kid) != BASKER_MAX_IDX)
	{

      #ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("\n-----------BLK---Kid: %d -------------\n",
	     kid);
      #endif

      //No longer needed, added in sfactor
      //printf("before workspace init\n");
      //basker->t_init_workspace(kid);
      //printf("after workspace init\n");
      
      //if(kid == 8)
	{
      basker->t_nfactor_blk(kid);
	}


	}

    }//end operator
  };//end kokkos_nfactor_domain_remalloc struct
 

  //use local number on local blks (Crazy idea)
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_nfactor_blk(Int kid)
  {
    Int b = S[0][kid]; //Which blk from schedule
    BASKER_MATRIX &L   = LL(b)(0);
    BASKER_MATRIX &U   = LU(b)(LU_size(b)-1);
    BASKER_MATRIX &M   = ALM(b)(0); //A->blk
#ifdef BASKER_2DL
    //printf("Accessing blk: %d kid: %d  \n", b, kid);
    INT_1DARRAY   ws   = LL(b)(0).iws;
    ENTRY_1DARRAY X    = LL(b)(0).ews;
    Int        ws_size = LL(b)(0).iws_size;
#else  //else if BASKER_2DL
    INT_1DARRAY   ws   = thread_array[kid].iws;
    ENTRY_1DARRAY X    = thread_array[kid].ews;
    Int       ws_size  = thread_array[kid].iws_size;
#endif

    //Int          bcol  = L.scol;  //begining col //NOT UD
    Int          brow  = L.srow;  //begining row //Note: move out in future
    Int          lval  = 0;
    Int          uval  = 0;


    Int i,j,k;
    //Int top, top1, maxindex, t; //NDE - warning: top1 set but not used
    Int top, maxindex, t;
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;

    Int newsize;
    Entry pivot, value;
    double absv, maxv;
    //Entry absv, maxv;

    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    //Int scol  = L.scol; //Note: this seems like over kill --clean up variables
    //Int ecol  = L.ecol; //Not used

    //Why did we need this?
    Int col_idx_offset = M.nnz;

    //printf("test one ws_size: %d \n", ws_size);

    //Note:
    Int *color    = &(ws(0));

    Int *pattern  = &(color[ws_size]);


    maxindex = BASKER_MAX_IDX;
    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
    //    top1 = ws_size; //NDE - warning: top1 set but not used

    lnnz = lval;
    unnz = uval;

#ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("b: %d  ls: %d us: %d llnzz: %d uunzz: %d \n", 
        b, lnnz, unnz, L.nnz, U.nnz);
    printf("b: %d gperm: %d \n", b, gperm(L.srow));
#endif

    //return 0 ;

    //---TEMP--- DEBUG ---
    //ecol = 5;

    //for each column
    //for(k = scol; k < ecol; k++)
    //Note: might want to add const trick for vectorize,
    //though this loop should really not vectorize

#ifdef BASKER_DEBUG_NFACTOR_BLK
    for(i = 0; i < M.ncol; i++)
    {

      if(L.pend(i) != BASKER_MAX_IDX)
      {
        printf("pend error: i: %d b: %d p: %d \n",
            i, b, L.pend(i));

      }
      BASKER_ASSERT(L.pend(i) == BASKER_MAX_IDX, "pend");

    }
#endif

    for(k = 0; k < M.ncol; ++k)
    {
#ifdef BASKER_DEBUG_NFACTOR_BLK
      //if((k+M.scol) == 183)
      {
        printf("\n----------------K=%d--------------\n", 
            k+M.scol);
      }
#endif
      value = 0.0;
      pivot = 0.0;
      //Why do we need?
      //maxindex = M.ncol;  
      lcnt = 0;
      ucnt = 0;

#ifdef BASKER_DEBUG_NFACTOR_BLK
      ASSERT(top == ws_size);
      //ASSERT entry workspace is clean
      for(i = 0 ; i < ws_size; i++)
      {
        BASKER_ASSERT(X(i) == 0, "Xerror");
      }
      //ASSERT int workspace is clean
      for(i = 0; i <  ws_size; i++)
      {
        BASKER_ASSERT(ws(i) == 0, "wserror");
      }

#endif

      //for each nnz in column
      //Want to change this to local blk anyway
      for(i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {

        j = M.row_idx(i);

#ifdef BASKER_2D
        //Do we need this anymore ?? Don't think
        if(j >= ecol)
        {
#ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("col_break, kid: %d idx: %d \n",
              kid, i);
#endif
          col_idx_offset = i;
          break;
        }
#endif

#ifdef BASKER_2DL
        X(j) = M.val(i);
#else
        X[j] = M.val[i];
#endif

#ifdef BASKER_DEBUG_NFACTOR_BLK
        if((k+M.srow) == 183)
        {
          printf("i: %d row: %d %d  val: %g  top: %d \n", 
              i, j, j+M.srow ,M.val(i), top);
#ifdef BASKER_2DL
          if(k == 120)
          {
            printf("Nx in Ak %d %g %d color = %d \n",
                j, X[j], brow,  
                color[j] );
          }
#else
          printf("Nx in Ak %d %g %d color = %d \n",
              j, X[j], brow,  
              color[j] );
#endif
        }
#endif

        //NOTE:  Need a quick skip of dfs if 
        //j i not pivotal (KLU)	      
        if(color[j] == 0)
        {
          //we want to skip the call if we can
          if(gperm(j+brow) != BASKER_MAX_IDX)
          {
            //printf("local_reach\n");
            t_local_reach(kid,0,0,j,top);
          }
          else
          {
            //printf("short\n");
            t_local_reach_short(kid,0,0,j,top);	 
          }
        }


      }//end for() each nnz in column
      xnnz = ws_size - top;


      //Debug
      //printf("TEST  x(%d) = %f \n",k , X(k));


#ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("xnnz: %d ws_size: %d top: %d \n", 
          xnnz, ws_size, top);
#endif

      //t_back_solve_selective(kid, 0, 0, k, top, xnnz);
      t_back_solve(kid, 0,0,  k, top, xnnz);

      //Future add
      //t_locate_pivot(kid, top)	  
      //find pivot
      maxv = 0.0;
      //Entry maxp = 0.0;
      double maxp = 0.0;
      for(i = top; i < ws_size; i++)
      {
        //j = pattern[i];
        j = pattern[i];
        //t = gperm[j];
        t = gperm(j+brow);
#ifdef BASKER_2DL
        //value = X[j-brow];
        value = X(j);
#else
        value = X[j];
#endif

#ifdef BASKER_DEBUG_NFACTOR_BLK
        //if((k+M.srow) == 183)
        {
          printf("k: %d consider: %d %d %g maxv %g %g  \n",
              k, j+M.srow, t, value, maxv, pivot);
        }
#endif

        absv = EntryOP::approxABS(value);

        if(j == k)
        {
          maxp = absv;
        }

        if(t == BASKER_MAX_IDX)
        {
          lcnt++;

          if(EntryOP::gt(absv,maxv))
          {
            maxv     = absv;
            pivot    = value;
            maxindex = j;
          }
        }
      }//for (i = top; i < ws_size)
      //printf("b: %d lcnt: %d after \n", b, lcnt);

      //Need a BIAS towards the diagonl
      //if(maxv < (Options.pivot_bias*maxp))
      if(maxv < (1.00001*maxp))
      {
        //printf("Close: %d %f %d %f \n",
        //     k, maxp, 
        //     maxindex, maxv);
        if(gperm(k+brow) == BASKER_MAX_IDX)
        {
          //  printf("using diag\n");
          pivot     = X(k);
          maxindex = k;
        }
      }

      if(Options.no_pivot == BASKER_TRUE)
      {
        maxindex = k;
        pivot    = X(k);
      }

      ucnt = ws_size - top - lcnt +1;
      //if((maxindex == L.max_idx) || (pivot == 0))
      if((maxindex == BASKER_MAX_IDX) || (pivot == (Entry)(0)) )
      {
        if (Options.verbose == BASKER_TRUE)
        {
          cout << endl << endl;
          cout << "---------------------------"<<endl;
          cout << "Error: Matrix is singular, blk" 
            << endl;
          cout << "k: " << k << " MaxIndex: " 
            << maxindex 
            << " pivot " 
            << pivot << endl;
          cout << "lcnt: " << lcnt << endl;
        }
        thread_array(kid).error_type =
          BASKER_ERROR_SINGULAR;
        thread_array(kid).error_blk   = b;
        thread_array(kid).error_subblk = 0; 
        thread_array(kid).error_info  = k;
        return BASKER_ERROR;
      }          

      gperm(maxindex+brow) = k+brow;
      gpermi(k+brow) = maxindex + brow;

#ifdef BASKER_DEBUG_NFACTOR
      //if(maxindex != k)
      //  {
      //    cout << "Permuting Pivot: " << k << " as row " 
      //         << maxindex << endl;
      //  }
#endif

      //Note: Come back to this!!!!
      if(lnnz + lcnt > llnnz)
      {
        //printf("\n\n");
        // printf("----------------------\n");

        newsize = lnnz * 1.1 + 2 *M.nrow + 1;

        if (Options.verbose == BASKER_TRUE)
        {
          printf("b: %ld Reallocing L oldsize: %ld current: %ld count: %ld newsize: %ld \n",
              (long)b, (long)llnnz, (long)lnnz, (long)lcnt, (long)newsize);
        }

        if(Options.realloc == BASKER_FALSE)
        {
          thread_array(kid).error_type =
            BASKER_ERROR_NOMALLOC;
          return BASKER_ERROR;
        }
        else
        {
          thread_array(kid).error_type =
            BASKER_ERROR_REMALLOC;
          thread_array(kid).error_blk    = b;
          thread_array(kid).error_subblk = 0;
          thread_array(kid).error_info   = newsize;
          return BASKER_ERROR;
        }

      }
      if(unnz+ucnt > uunnz)
      {

        // printf("\n\n");
        //printf("-------------------\n");

        newsize = uunnz*1.1 + 2*M.nrow+1;

        if (Options.verbose == BASKER_TRUE)
        {
          printf("b: %ld Reallocing U oldsize: %ld newsize: %ld  k: %ld \n",
              (long)b, (long)uunnz, (long)unnz+ucnt, (long)k);
        }

        if(Options.realloc == BASKER_FALSE)
        {
          thread_array(kid).error_type =
            BASKER_ERROR_NOMALLOC;
          return BASKER_ERROR;
        }
        else
        {
          thread_array(kid).error_type =
            BASKER_ERROR_REMALLOC;
          thread_array(kid).error_blk    = b;
          thread_array(kid).error_subblk = -1;
          thread_array(kid).error_info   = newsize;
          return BASKER_ERROR;
        }

      }

      L.row_idx(lnnz) = maxindex;
      L.val(lnnz)     = (Entry) 1.0;
      lnnz++;

      Entry lastU = (Entry) 0.0;
      for( i = top; i < ws_size; i++)
      {
        //j = pattern[i];
        j = pattern[i];
        //t = gperm[j];
        t = gperm(j+brow);

#ifdef BASKER_DEBUG_NFACTOR_BLK
        printf("j: %d t: %d x: %g \n", j, t, X(j));
#endif            

        //Note can not excludude numeric cancel in prun
        //if fill-in
#ifdef BASKER_2DL
        //if(X[j-brow] != 0)
        //if(X(j) != 0)
#else
        //if(X[j] != 0)
#endif
        {
          //if(t != L.max_idx)
          if(t != BASKER_MAX_IDX)
          {
            if(t < (k+brow))
            {

              //U.row_idx[unnz] = gperm[j];
              //U.row_idx(unnz) = gperm(j+brow);
              U.row_idx(unnz) = t-brow;
#ifdef BASKER_2DL
              //U.val[unnz] = X[j-brow];
              U.val(unnz) = X(j);
#else
              U.val[unnz] = X[j];
#endif

#ifdef BASKER_DEBUG_NFACTOR_BLK
              //debug
              if((k+M.srow)==183)
              {
                printf("U(%d,%d): %f \n",
                    gperm(j+brow), k+bcol,
                    X(j));
              }
#endif

              unnz++;
            }
            else
            {
#ifdef BASKER_2DL
              //lastU = X[j-brow];
              lastU = X(j);
#else
              lastU = X[j];
#endif
            }
          }
          else if (t == BASKER_MAX_IDX)
          {
            L.row_idx(lnnz) = j;
#ifdef BASKER_2DL
            //L.val(lnnz) = X(j)/pivot;
            L.val(lnnz) = EntryOP::divide(X(j),pivot);
#else
            //L.val[lnnz] = X[j]/pivot;
            L.val(lnnz) = EntryOP::divde(X(j),pivot);
#endif

#ifdef BASKER_DEBUG_NFACTOR_BLK
            if((k+M.srow)==183)
            {
              printf("L(%d,%d): %f \n",
                  j, k+bcol,
                  L.val(lnnz));
            }
#endif
            //Need to comeback for local convert
            //#ifdef BASKER_INC_LVL
            //L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
            //#endif

            lnnz++;

          }
        }//end if() not 0

        //Note: move x[j] inside of if() not 0....extra ops this way
#ifdef BASKER_DEBUG_NFACTOR_BLK
        printf("Zeroing element: %d \n", j);
#endif

#ifdef BASKER_2DL
        //X[j-brow] = 0;
        X(j) = 0;
#else
        X[j] = 0;
#endif
      }//end if(x[i] != 0)

      //Fill in last element of U
      U.row_idx(unnz) = k;
      U.val(unnz) = lastU;
      unnz++;

      xnnz = 0;
      top = ws_size;

      L.col_ptr(k) = cu_ltop;
      L.col_ptr(k+1) = lnnz;
      cu_ltop = lnnz;

      U.col_ptr(k) = cu_utop;
      U.col_ptr(k+1) = unnz;
      cu_utop = unnz;

      //printf("U_col: %d \n",
      //	 U.col_ptr(k));

#ifdef BASKER_2DL
      //-----------------------Update offdiag-------------//
      for(Int blk_row = 1; blk_row < LL_size(b); ++blk_row)
      {
        //Do back solve of off-diag blocks
#ifdef BASKER_INC_LVL
        //t_back_solve_offdiag_selective(kid,
        //		   b, blk_row,
        //		   b, blk_row,
        //		   k, col_idx_offset,
        //		   U.val, U.row_idx,
        //       U.col_ptr(k-bcol+1)-U.col_ptr(k-bcol),
        //		  U.col_ptr(k-bcol),
        //		   BASKER_TRUE);
#else

        //t_back_solve_offdiag(kid,
        //		   b, blk_row,
        //		   b, blk_row,
        //		   k, col_idx_offset,
        //		   U.val, U.row_idx,
        //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
        //		   U.col_ptr[k-bcol],
        //		   BASKER_TRUE);


        t_back_solve_offdiag(kid,
            b, blk_row,
            b, blk_row,
            //Note, different meaning
            k, col_idx_offset,
            U.val, U.row_idx,
            //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
            U.col_ptr(k+1)-U.col_ptr(k),
            // U.col_ptr[k-bcol],
            U.col_ptr(k),
            BASKER_TRUE);

#endif
        //Move these factors into Local Ls
        Int move_error = 
          t_move_offdiag_L(kid,
              b, blk_row,
              b, blk_row,
              k, pivot);

        if(move_error == BASKER_ERROR)
        {
          return BASKER_ERROR;
        }
      }//end over all diag
#endif

      //Why?
      col_idx_offset = A.nnz;

      t_prune(kid,0,0,k,maxindex);

    }//end for() over all columns

    L.nnz = lnnz;
    U.nnz = unnz;

#ifdef BASKER_DEBUG_NFACTOR_BLK
    //print_factor(L,U);
#endif
    return 0;
  }//end t_nfactor_blk()


 
  template<class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  void Basker<Int,Entry,Exe_Space>::t_local_reach_short
  (
   const Int kid, 
   const Int lvl,
   const Int l,
   const Int j,
   Int &top
   )
  {
    //Setup variables
//    const Int b      = S(lvl)(kid); //NDE - warning: unused 
    const Int wsb    = S(0)(kid);
    //BASKER_MATRIX &L = LL(b)(0); //NDE - warning: unused L

    INT_1DARRAY  ws   = LL(wsb)(l).iws;
    const Int ws_size = LL(wsb)(l).iws_size;

    Int *color       = &(ws(0));
    Int *pattern     = &(ws(ws_size));

    color[j]       = 2;
    pattern[--top] = j;

  }//end t_local_reach short()
  

  template <class Int, class Entry, class Exe_Space>
  inline
  void Basker<Int,Entry,Exe_Space>::t_prune
  (
   const Int kid,
   const Int lvl, 
   const Int l, 
   const Int k, 
   const Int pivotrow 
   )
  {
    const Int b      = S(lvl)(kid);
    //const Int wsb    = S(0)(kid);
    BASKER_MATRIX &L = LL(b)(0);
    const Int U_col  = S(lvl)(kid);
    Int U_row        = LU_size(U_col)-1;
    if(lvl > 0)
      {
	//U_row = (lvl==1)?(kid%2):S(l)(kid)%LU_size(U_col);
      }
    BASKER_MATRIX &U = LU(U_col)(U_row);
    //const Int brow   = L.srow;


    //printf("t_prune,k: %d  L: %d %d U: %d %d pivotrow: %d\n", 
    //	   k,b,0, U_col, U_row, pivotrow);

    //Scan U and check if we any columns we can prune :)
    //don't include the last entry
    for(Int ui = U.col_ptr(k); ui < U.col_ptr(k+1)-1; ++ui)
      {
	
	Int j = U.row_idx(ui);

	//printf("prune, j: %d %d k: %d \n", 
	//     j, j+brow, k);


	//if(j>=k)
	  {
	    //printf("Error: j: %d k: %d U: %d %d ui: %d %d \n",
	    //j, k, U_col, U_row, U.col_ptr(k), U.col_ptr(k+1));
	BASKER_ASSERT(j < k, "Pruning, j not less than k") ;
	  }
	
	//Check if this column has not been pruned
	if(L.pend(j)==BASKER_MAX_IDX)
	  {
	    //printf("This column can be pruned\n");

	    for(Int li = L.col_ptr(j); li < L.col_ptr(j+1); ++li)
	      {
		//printf("considering row: %d %d \n",
		//     L.row_idx(li),
		//     L.row_idx(li)+L.srow);

		//if we can prune ?
		if(pivotrow == L.row_idx(li))
		  {
		    
		    //printf("This can be pruned\n");

		    //order L 
		    //partition those [above pivotrow, below]
		    Int phead = L.col_ptr(j);
		    Int ptail = L.col_ptr(j+1);
		    
		    //printf("phead: %d ptail: %d \n",
		    //	   phead, ptail);

		    //partion and bubble up
		    while(phead < ptail)
		      {
			Int i = L.row_idx(phead);
			if(gperm(i+L.srow) >= 0)
			  {
			    //advance top ptr
			    //printf("adv head(%d) %d %d \n",
			    //phead+1, 
			    //	   i+L.srow,
			    //	   gperm(i+L.srow));
			    phead++;
			  }
			else
			  {
			    //flip head and tail
			    //printf("flipping head:%d tail:%d \n",
			    //	   phead, ptail);
			    //printf("head: %d %d %f \n",
			    //	   phead, 
			    //	   L.row_idx(phead)+L.srow,
			    //	   L.val(phead));
			    //printf("tail: %d %d %f \n",
			    //	   ptail-1,
			    //	   L.row_idx(ptail-1)+L.srow,
			    //	   L.val(ptail-1));

			    ptail--;
			    L.row_idx(phead) = L.row_idx(ptail);
			    L.row_idx(ptail) = i;
			    Entry x = L.val(phead);
			    L.val(phead) = L.val(ptail);
			    L.val(ptail) = x;
			  }//End Flip
		      }//end -- while(for flip)
		    
		    //printf("prune:, ptail: %d \n",
		    //	   ptail);

		    //if(j ==39)
		    // {
		    //	printf("bad prune.  k: %d %d \n",
		    //	       k, k+L.scol);
		    // }
		    //printf("pruned, j:%d ptail: %d\n",
		    //	   j, ptail);
		    L.pend(j) = ptail;
		    break;
		    
		  }//if pivotrow == L.row_idx(..)
	      }//end scan over all Li
	  }//if L.pend == enmpty
      }//over all Ui
  }//end t_prune()
 

  //updated local reach based
  //Note that this function needs to be tight
  template <class Int, class Entry, class Exe_Space>
  inline
  void Basker<Int,Entry,Exe_Space>::t_local_reach
  (
   const Int kid, 
   const Int lvl,
   const Int l,
   Int j, 
   Int &top
   )
  {
    
    //Setup variables
    const Int      b   = S(lvl)(kid);
    const Int     wsb  = S(0)(kid);
    BASKER_MATRIX  &L  = LL(b)(0);
    const Int     brow = L.srow;

    INT_1DARRAY    ws  = LL(wsb)(l).iws;
    const Int  ws_size = LL(wsb)(l).iws_size;
   
    //Int *color       = &(ws[0]);
    Int *pattern     = &(ws(ws_size));
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);
    

    Int i, t, head, i1;
//    Int start, end, done; //NDE - warning: done set but not used
    Int start, end; //NDE - warning: done set but not used
    

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d j: %d, kid: %d \n",
	   b, 0, wsb, l, j, kid);
    #endif

    start    = -1;
    head     = 0;
    stack[0] = j;

    while(head != BASKER_MAX_IDX)
      { 
        #ifdef BASKER_DEBUG_LOCAL_REACH
        printf("stack_offset: %d head: %d \n", 
	       stack_offset , head);
        BASKER_ASSERT(head > -1, " ");
        #endif

	j = stack[head];
	t = gperm(j+brow);
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("--------DFS: %d %d %d -------------\n",
	       j, j+brow, t);
        #endif

	if(ws(j) == 0)
	  {	    

	    ws(j) = 1;
  
	    if(t!=BASKER_MAX_IDX)
	      {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("reach.... j: %d t:%d L.scol %d \n",
		       j, t, L.scol);
		printf("colend: %d pend(%d): %d %d \n",
		       L.col_ptr(t+1-L.scol),
		       t-L.scol,
		       L.pend(j),
		       L.pend(t-L.scol));
                #endif
		//start = L.col_ptr(t+1-L.scol);
		start = 
		(L.pend(t-L.scol)==BASKER_MAX_IDX)?L.col_ptr(t+1-L.scol):L.pend(t-L.scol); 

	      }
	  }
	else
	  {
	    start = store[j];
	  }
//	done = 1; //NDE - warning: done set but not used
	

	//We want to go backwards through this
	//This little insight can save time
	//end = L.col_ptr(t+1-L.scol);
	//printf("t-L %d \n", t-L.scol);
	end = L.col_ptr(t-L.scol);
	//printf("t: %d start: %d end: %d \n",t, start, end);
	for(i1 = --start; i1 >= end; --i1)
	  {
	    i = L.row_idx(i1);
	   
	    //printf("Search i1: %d  i: %d %d %d \n",
	    //	   i1, i, i+L.scol, gperm(i+L.scol));

	    if(ws(i) != 0)
	      {
		//printf("continue called\n");
		continue;
	      }
	    else
	    //if(ws(i) == 0)
	      {
		if(gperm(i+brow) >= 0)
		  {
		    //store[j] = i1+1;
		    store[j] = i1;
		    stack[++head] = i;
		    break;
		  }
		else
		  {
		    //color[i] = 2;
		    ws(i) = 2;
		    pattern[--top] = i;
		    //printf("Adding idx: %d %d  to pattern at location: %d \n",i, i+L.scol, top);

		  }
	      }
	   


	  }
	
	//if we have used up the whole column
	if(i1 < end)
	  {
	    head--;
	    //color[j] = 2;
	    ws(j) = 2;
	    pattern[--top] = j;
    
	  }
	//printf("head: %d \n", head);
	/*
	if(head == 0)
	  {head = BASKER_MAX_IDX;}
	else
	  {--head;}
	*/
	//printf("Adding idx: %d to pattern at location: %d \n",j, *top);
      	//printf("head2: %d \n", head);
      }//end while 

  }//end t_local_reach (new)

   //Note: need to fix indexing for color
  //Uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach_old
  (Int kid, Int lvl, Int l,
   Int j, Int *top)
  {

    //Setup variables
    const Int      b   = S(lvl)(kid);
    const Int     wsb  = S(0)(kid);
    BASKER_MATRIX  &L  = LL(b)(0);
    #ifdef BASKER_2DL
    INT_1DARRAY    ws  = LL(wsb)(l).iws;
    const Int  ws_size = LL(wsb)(l).iws_size;
    #else
    INT_1DARRAY    ws  = thread_array[kid].iws;
    Int        ws_size = thread_array[kid].iws_size;
    #endif

    const Int brow   = L.srow;

    Int *color       = &(ws[0]);
    Int *pattern     = &(ws[ws_size]);
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);
    


    /*
    auto color = Kokkos::subview(ws, 
				 std::make_pair((Int)0,
						(Int)ws_size+1));
    auto pattern = Kokkos::subview(ws,
				   std::make_pair((Int)ws_size,
					      (Int)(2*ws_size)-1));
    auto stack   = Kokkos::subview(ws,
			     std::make_pair((Int)2*ws_size, 
					    (Int)(3*ws_size)-1));
    auto store   = Kokkos::subview(ws,
				 std::make_pair((Int)3*ws_size,
					   (Int)(4*ws_size)-1));
    */

    Int i, t, head, i1;
    Int start, end, done;
    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d, kid: %d \n",
	   b, 0, wsb, l, kid);
    #endif

    start = -1;

    head = 0;
    

    //stack(head) = j;
    stack[head] = j;

    //while(head < L.max_idx)
    while(head != BASKER_MAX_IDX)
      { 
        #ifdef BASKER_DEBUG_LOCAL_REACH
        printf("stack_offset: %d head: %d \n", stack_offset , head);
        ASSERT(head > -1);
        #endif

	j = stack[head];
	//j = stack(head);
        //t = gperm[j];
	t = gperm(j+brow);
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("----------DFS: %d %d -------------\n", j, t);
        #endif

	#ifdef BASKER_2DL
	if(color[j] == 0)
	//if(color(j) == 0)
	#else
        if(color[j] == 0)
	#endif
	  {	    
	    #ifdef BASKER_2DL
	    color[j] = 1;
	    //color(j) = 1;
	    #else
	    color[j] = 1;
	    #endif
  
            //if not permuted, might be able to make this smaller now
            //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
	    //if(t!=L.max_idx)
	    if(t!=BASKER_MAX_IDX)
	      {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("reach.... j: %d t:%d L.scol %d \n", j, t, L.scol);
                #endif
                //start = L.col_ptr[t-L.scol];
		start = L.col_ptr(t-L.scol);
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
	    start = store[j];
	    //start = store(j);
	  }
	done = 1;
	
        //if  permuted
        //if(t != L.max_idx), might be able to make this smaller now
        //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
	//if(t!=L.max_idx)
	if(t!=BASKER_MAX_IDX)
	  {
	   
	    //We want to go through this backward
	    //That way our first search will end fast
	    //our later searchs will use our first searchs
	    //end = L.col_ptr[t+1-L.scol];
	    end = L.col_ptr(t+1-L.scol);
	    for(i1 = start; i1 < end; i1++)
	    //for(il = 
	      {
                i = L.row_idx(i1);
		#ifdef BASKER_2DL
		if(color[i] == 0)
		//if(color(i) == 0)
		#else
		if(color[i] == 0)
		#endif
		  {
                    head++;
		    stack[head] = i;
		    //stack(head) = i;
		    store[j] = i1+1;
		    //store(j) = i1+1;
		    done = 0;
		    
		    #ifdef BASKER_INC_LVL
		    inc_lvl++;
		    #endif

		    break;
		  }
	      }
	  }
	if(done)
	  {
            pattern[--*top] = j;
	    //pattern(--*top) = j;
	    #ifdef BASKER_2DL
	    color[j] = 2;
	    //color(j) = 2;
	    #else
	    color[j] = 2;
	    #endif
	    if(head == 0)
	      {head = BASKER_MAX_IDX;}
	      //{ head = L.max_idx;}
	    else
	      {head--;}

	    //printf("Adding idx: %d to pattern at location: %d \n",j, *top);


	  }
      }//end while 
    return 0;
  }//end t_local_reach()




  //local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_FINLINE
  int Basker<Int, Entry,Exe_Space>::t_back_solve
  (
   Int kid, 
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz)
  {

    const Int      b = S(lvl)(kid);
    const Int    wsb = S(0)(kid);
    BASKER_MATRIX &L = LL(b)(0);
    #ifdef BASKER_2DL
    INT_1DARRAY   ws = LL(wsb)(l).iws;
    ENTRY_1DARRAY X  = LL(wsb)(l).ews;
    Int      ws_size = LL(wsb)(l).iws_size;
    #else
    INT_1DARRAY   ws = thread_array[kid].iws;
    ENTRY_1DARRAY  X = thread_array[kid].ews;
    Int      ws_size = thread_array[kid].iws_size;
    #endif
    
    Int *color   = &(ws[0]);
    Int *pattern = &(color[ws_size]);
    
    Int brow = L.srow;

    Int top1 = top;
    Int j,t;
    
    for(Int pp = 0; pp < xnnz; ++pp)
      {
	j = pattern[top1];
	color[j] = 0;
	++top1;
      }
    
    top1 = top;

    for(Int pp = 0; pp < xnnz; ++pp)
      {
	j = pattern[top1];
	t = gperm(j+brow);

	//color[j] = 0;
	
	if(t!=BASKER_MAX_IDX)
          {
            const Int local_offset = L.scol; //L.srow
	    const Int pend = L.col_ptr(t+1-local_offset);
	    const Entry xj = X(j);

	    //We could use kokkos local memory,
	    //to get rid of this pragma
	    //
	    //What we want to think about is how to 
	    //do this in a dense way
	    //#pragma ivdep (may be slower)
	    for(Int p = L.col_ptr(t-local_offset)+1;
		p < pend; ++p)
              {
		const Int row_idx = L.row_idx(p);
		const Entry update_val = L.val(p)*xj;
		
		X(row_idx) -= update_val;

              }//end for() over each nnz in the column
          }//end if() not permuted
	top1++;
      }//end for() over all nnz in LHS

    return 0;
  }//end t_back_solve()


   //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_dense_move_offdiag_L
  (Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot)
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
    //const Int   ws_size = LL(X_col)(X_row).iws_size;
    //const Int   p_size  = LL(X_col)(X_row).p_size;
   

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
    #endif

   
    Int *color   = &(ws(0));
    //Int *pattern = &(color[ws_size]);

    //const Int    brow  = L.srow;
    //const Int    bcol  = L.scol;
    //const Int    llnnz = L.nnz;
    Int    lnnz  = L.col_ptr(k);
   
    /*
    if((p_size) > (llnnz-lnnz))
      {
	printf("-Warning, Need to remalloc L: %d %d kid: %d current size: %d used_size: %d  addition: %d \n",
	       blkcol, blkrow, kid, llnnz,lnnz,p_size  );
	
      }
    */

    ///for(Int i = 0; i < p_size; i++)
    for(Int j = 0; j < L.nrow; ++j)
      {
	//Int j = pattern[i];
	//Int t = gperm(j);
	if(X(j) != (Entry)(0) )
	  {

            //Int t = gperm(j+brow);
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, kid: %d j: %d val: %f lnnz: %d \n",
	     kid, j, X[j]/pivot, lnnz);
	#endif

	color[j] = 0;
	
	L.row_idx(lnnz) = j;

	//L.val(lnnz) = X(j)/pivot;
	L.val(lnnz) = EntryOP::divide(X(j),pivot);

	X(j) = 0;

	//Note come back to fix
	#ifdef BASKER_INC_LVL
	L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
	#endif

	lnnz++;
	  }
      }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //	   kid, k-bcol+1, lnnz);
    
    //L.col_ptr[k-bcol+1] = lnnz;
    L.col_ptr(k+1) = lnnz;

    //LL[X_col][X_row].p_size = 0;
    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_dense_offdiag_mov_L()


  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_move_offdiag_L
  (Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot)
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
    const Int   ws_size = LL(X_col)(X_row).iws_size;
    const Int   p_size  = LL(X_col)(X_row).p_size;
   

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 8 )
      {
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
      }
    #endif

   
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);

    //const Int    brow  = L.srow;
    //const Int    bcol  = L.scol;
    const Int    llnnz = L.nnz;
          Int    lnnz  = L.col_ptr(k);

    if((p_size) > (llnnz-lnnz))
      {

	Int newsize = llnnz*1.2 + L.ncol;

	 if(Options.realloc == BASKER_FALSE)
	   {
	     thread_array(kid).error_type =
	       BASKER_ERROR_NOMALLOC;
	     return BASKER_ERROR;
	   }
	 else
	   {
	     thread_array(kid).error_type =
	       BASKER_ERROR_REMALLOC;
	     thread_array(kid).error_blk    = blkcol;
	     thread_array(kid).error_subblk = blkrow;
	     thread_array(kid).error_info   = newsize;
	     return BASKER_ERROR;
	   }

	 if (Options.verbose == BASKER_TRUE)
	   {
	printf("-Warning, Need to remalloc L: %ld %ld kid: %ld current size: %ld used_size: %ld  addition: %ld \n",
	       (long)blkcol, (long)blkrow, (long)kid, (long)llnnz, (long)lnnz, (long)p_size  );
	   }
	//BASKER_ASSERT(0==1, "REALLOC LOWER BLOCK\n");
	
      }

    for(Int i = 0; i < p_size; i++)
      {
	Int j = pattern[i];
	//Int t = gperm(j+brow);
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, kid: %d j: %d val: %f lnnz: %d \n",
	     kid, j, X[j]/pivot, lnnz);
	#endif

	color[j] = 0;
	
	L.row_idx(lnnz) = j;

	//L.val(lnnz) = X(j)/pivot;
	L.val(lnnz) = EntryOP::divide(X(j),pivot);

	X(j) = 0;

	//Note come back to fix
	#ifdef BASKER_INC_LVL
	L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
	#endif

	lnnz++;
      }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //	   kid, k-bcol+1, lnnz);
    
    L.col_ptr(k+1) = lnnz;

    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_offdiag_mov_L()

  
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_dense_back_solve_offdiag
  (
   Int kid, 
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option)
  {
    //Note:  need to add support for offdiag permuation

    BASKER_MATRIX &L            = LL(blkcol)(blkrow);
    BASKER_MATRIX &B            = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws            = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X             = LL(X_col)(X_row).ews;
    //Int         ws_size         = LL(X_col)(X_row).iws_size;
    
    Int    nnz            = LL(X_col)(X_row).p_size;
    //const Int    brow           = L.srow;
    //const Int    bcol           = L.scol;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("\n\n");
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    printf("\n\n");
    #endif
    // B.info();
    //B.print();

    //Int *color =   &(ws(0));
    //Int *pattern = &(color[ws_size]);
    
    //Preload with A
    if(A_option == BASKER_TRUE)
      {
	#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	#endif
	//for(Int i = view_offset; i < B.m_offset; i++)
	//printf("t_b_s_off debug, kid: %d k: %d bcol: %d col_ptr: %d \n",
	//     kid, k, bcol, B.col_ptr[k-bcol]);
	for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
	  {
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //Bgood(remove)
	    //printf("t_back_solve_diag, kid: %d i: %d g: %d\n",
	    //	   kid, i, B.good(i));
	    #endif

	    const Int j = B.row_idx(i);

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	    printf("t_back_solve_diag, kid: %d x: %f %f \n",
		   kid, X(j), X(j)+B.val(i));
	    #endif
	    	   
	    X(j) = X(j)+ B.val(i);
	    
	  }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
	   kid, x_size);
    #endif

 

    //printf("x_offset: %d \n", x_offset);

    for(Int i = 0 ; i < x_size; ++i)
      {
	//const Int k    =   x_idx[i+x_offset];
	const Int k = x_idx[i+x_offset];
	//const Entry xj =   x[i+x_offset];
	const Entry xj = x(i+x_offset);
	//printf("kid: %d bcol: %d k: %d \n",
	//   kid, bcol, k);
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	//printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	//   kid, k, L.col_ptr[k-bcol], L.col_ptr[k-bcol+1]);
	printf("t_back_solve_diag, kid: %d  k: %d %g  x_size: %d [%d %d] \n",
	       kid, k, xj, x_size,  L.col_ptr[k], L.col_ptr[k+1]);
	#endif
	
	//for(Int j = L.col_ptr[k-bcol]; 
	//  j < L.col_ptr[k-bcol+1]; j++)
	//printf("ERROR, kid: %d k: %d \n", kid, k);
	//printf("ERROR, kid: %d blkcol: %d blkrow: %d \n",
	//     kid, blkcol, blkrow);
	//printf("ERROR, kid: %d k: %d max: %d \n",
	//     kid, k, L.ncol);
	for(Int j = L.col_ptr(k); 
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //printf("t_b_solve_d, kid: %d j: %d color: %d \n",
	    //	   kid, jj, color[jj-brow]);
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    #endif

	    //if(color[jj-brow] != 1)
	    //Do not need this in dense
	    /*
	    if(color[jj] != 1)
	      {
		//color[jj-brow] = 1;
		color[jj] = 1;
		#ifdef BASKER_DEBUG_NFACTOR_BLK
		printf("pattern index: %d kid: %d \n",
		       nnz, kid);
		#endif
		pattern[nnz++] = jj;
		//	printf("t_b_solve_d, id: %d nnz: %d\n",
		//     kid, nnz);
		//	printf("-----------PATTERN UPDATE kid: %d L: %d %d pattern(%d) = %d brow: %d \n", 
		//     kid, X_col, X_row, nnz-1, pattern[nnz-1], brow);

	      }
	    */

	    #ifdef BASKER_DEBUG_NFACTOR_BLK

	    // printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
	    //kid, jj,X[jj-brow], L.val[j], xj);
	     printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
		    kid, jj,X[jj], L.val[j], xj);
	     #endif 

	     X(jj) -= L.val(j)*xj;
	     //X[jj-brow] -= L.val[j]*xj;
	  }


	/*
	printf("\n");
	for(Int i = 0; i < L.nrow; ++i)
      {
	Int jj = i;
	//Int jj = pattern[i];
	//printf("X[%d] = %f , kid: %d  \n",
	//     jj, X[jj-brow], kid);
	printf("k: %d X[%d](%d) = %f , kid: %d  \n",
	       k, jj,jj+B.srow, X[jj], kid);

      }
	printf("\n");
	*/
      }//over all nonzero in left

    ///Just scan for pattern
    //We will want to remove this in the future
    /*
    for(Int i = 0; i < L.nrow; ++i)
      {
	if(X(i) != (Entry)0)
	  {
	    pattern[nnz++] = i;
	    color[i]       = 1;
	  }
      }
    */

    #ifdef BASKER_2DL
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    //LL[X_col][X_row].p_size = nnz;
    LL(X_col)(X_row).p_size = nnz;
    #endif


    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("kid: %d all values: \n", kid);
    for(Int i = 0; i < L.nrow; ++i)
      {
	Int jj = i;
	//Int jj = pattern[i];
	//printf("X[%d] = %f , kid: %d  \n",
	//     jj, X[jj-brow], kid);
	printf("X[%d](%d) = %f , kid: %d  \n",
	       jj,jj+B.srow, X[jj], kid);

      }
    printf("\n\n");
    #endif

    return 0;
  }//end t_offdiag_dense_back_solve();


  //local idx for local blks
  //Bgood done.
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag
  (
   Int kid, 
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option)
  {
    //Note:  need to add support for offdiag permuation

    BASKER_MATRIX &L            = LL(blkcol)(blkrow);
    BASKER_MATRIX &B            = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws            = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X             = LL(X_col)(X_row).ews;
    Int         ws_size         = LL(X_col)(X_row).iws_size;
    
    Int    nnz            = LL(X_col)(X_row).p_size;
    //const Int    brow           = L.srow;
    //const Int    bcol           = L.scol;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 8)
      {
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    printf("t_back_solve_diag, kid: %d row: %d %d col: %d %d \n",
	   kid, L.srow, L.srow+L.nrow, L.scol, L.scol+L.ncol);


      }
    #endif
    // B.info();
    //B.print();

    Int *color =   &(ws(0));
    Int *pattern = &(color[ws_size]);
    
    //Preload with A
    if(A_option == BASKER_TRUE)
      {
	//#ifdef BASKER_DEBUG_NFACTROR_BLK
	//printf("t_back_solve, A_OPTION TRUE \n");
	//#endif
	//for(Int i = view_offset; i < B.m_offset; i++)
	//printf("t_b_s_off debug, kid: %d k: %d bcol: %d col_ptr: %d \n",
	//     kid, k, bcol, B.col_ptr[k-bcol]);
	for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
	  {
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //Bgood(remove)
	    //printf("t_back_solve_diag, kid: %d i: %d g: %d\n",
	    //	   kid, i, B.good(i));
	    #endif


            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    if(kid == 8)
	      {
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	      }
	    #endif
	    const Int j = B.row_idx(i);
	    //color[j-brow] = 1;
	    color[j] = 1;
	    //X[j-brow] = B.val(i);
	    X(j) = B.val(i);
	    pattern[nnz++] = j;
	  }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 8)
      {
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
	   kid, x_size);
      }
    #endif
    for(Int i = 0 ; i < x_size; ++i)
      {
	//const Int k    =   x_idx[i+x_offset];
	const Int k = x_idx[i+x_offset];
	//const Entry xj =   x[i+x_offset];
	const Entry xj = x(i+x_offset);
	//printf("kid: %d bcol: %d k: %d \n",
	//   kid, bcol, k);

	#ifdef BASKER_DEBUG_NFACTOR_BLK
	if(kid == 8)
	  {
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	     kid, k, L.col_ptr[k], L.col_ptr[k+1]);
	  }
	#endif
	
	//for(Int j = L.col_ptr[k-bcol]; 
	//  j < L.col_ptr[k-bcol+1]; j++)
	for(Int j = L.col_ptr(k); 
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    if(kid == 8)
	      {

	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	      }
	    #endif

	    //if(color[jj-brow] != 1)
	    if(color[jj] != 1)
	      {
		//color[jj-brow] = 1;
		color[jj] = 1;
		#ifdef BASKER_DEBUG_NFACTOR_BLK
		if(kid == 8)
		  {
		printf("pattern index: %d kid: %d \n",
		       nnz, kid);
		  }
		#endif
		pattern[nnz++] = jj;
		//	printf("t_b_solve_d, id: %d nnz: %d\n",
		//     kid, nnz);
		//	printf("-----------PATTERN UPDATE kid: %d L: %d %d pattern(%d) = %d brow: %d \n", 
		//     kid, X_col, X_row, nnz-1, pattern[nnz-1], brow);

	      }

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    if(kid == 8)
	      {
	     printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
	    kid, jj,X[jj], L.val[j], xj);
	      }
	    #endif 

	     X(jj) -= L.val(j)*xj;
	     //X[jj-brow] -= L.val[j]*xj;
	  }
      }//over all nonzero in left


    #ifdef BASKER_2DL
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    //LL[X_col][X_row].p_size = nnz;
    LL(X_col)(X_row).p_size = nnz;
    #endif


    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 8)
      {
    printf("kid: %d all values: \n", kid);
    for(Int i = 0; i < nnz; ++i)
      {
	Int jj = pattern[i];
	//printf("X[%d] = %f , kid: %d  \n",
	//     jj, X[jj-brow], kid);
	printf("X[%d](%d) = %f , kid: %d  \n",
	       jj,jj+B.srow, X[jj], kid);

      }
    printf("\n\n");
      }
    #endif
      
    return 0;
  }//end t_offdiag_back_solve();





 

}//end namespace baskerNS -- contains functions

#endif//end ifndef basker_nfactor_blk
