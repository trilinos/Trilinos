#ifndef BASKER_NFACTOR_DIAG_HPP
#define BASKER_NFACTOR_DIAG_HPP


#include "basker_matrix_decl.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"
#include "basker_types.hpp"
#include "basker_stats.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

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
      Int total_threads = thread.league_size()*
	thread.team_size();


      //================OLD EQ BLK=================//

      //divide equally and launch
      //Int nchunks = (basker->btf_nblks - 
      //			   basker->btf_tabs_offset);
      //Int chunk_size = nchunks/total_threads;
      //	    Int chunks_left_over = nchunks - 
      //	      (chunk_size*total_threads);
      //	    Int chunk_start = kid*chunk_size;

      //Divid into chunks based on work
     
 	    
      //printf("kid: %d c_size: %d numc: %d schunk: %d \n", kid, chunk_size, nchunks, chunk_start);
      
      //basker->t_nfactor_diag(kid, chunk_start, 
      //			     chunk_size);
      
      //extra
      //if(kid == 0)
      //{
      //	  Int extra_start=chunk_size*total_threads;
      //	  basker->t_nfactor_diag(kid, extra_start,
      //				 chunks_left_over);
	  
      //	}


      //============= NEW EQ  =============//

      Int chunk_start = basker->btf_schedule(kid);
      Int chunk_size  = basker->btf_schedule(kid+1) - 
	basker->btf_schedule(kid);
      

      printf("Chunk start: %d size: %d \n", 
	     chunk_start, chunk_size);
      basker->t_nfactor_diag(kid, chunk_start,
			     chunk_size);




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
    
    INT_1DARRAY thread_start;

    kokkos_nfactor_diag_remalloc()
    {}

    kokkos_nfactor_diag_remalloc
    (
     Basker<Int,Entry,Exe_Space> *_basker,
     INT_1DARRAY                 _thread_start
     )
    {
      basker       = _basker;
      thread_start = _thread_start;
    }

    //comeback and fix kokkos stuff

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {
      Int kid = (Int)(thread.league_rank()*
		      thread.team_size()+
		      thread.team_rank());
      Int total_threads = thread.league_size()*
	thread.team_size();


      Int chunk_start = thread_start(kid);
      if(chunk_start != BASKER_MAX_IDX)
	{
	  Int chunk_size  = basker->btf_schedule(kid+1) - 
	    chunk_start;
      
	  printf("Chunk start: %d size: %d \n", 
		 chunk_start, chunk_size);
	  basker->t_nfactor_diag(kid, chunk_start,
				 chunk_size);

	}//end if

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

    for(Int c = schunk+btf_tabs_offset;
	c < schunk+nchunk+btf_tabs_offset; ++c)
      {
	Int c_size = btf_tabs(c+1)-btf_tabs(c);
	//printf("kid: %d current_chunk: %d size: %d \n",
	//     kid, c, c_size);

	if(c_size == 1)
	  {
	    //call single
	    t_single_nfactor(kid, c);
	  }
	else
	  {
	    //call GP alg.
	    t_blk_nfactor(kid, c);
	  }

      }//over all chunks 

    return 0;
  }//end t_nfactor_diag

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_single_nfactor
  (
   Int kid,
   Int c
   )
  {  
    Int bcol = BTF_C.scol;
    Int brow = BTF_C.srow;
    Int btab = btf_tabs_offset;
    
    BASKER_MATRIX  &M = BTF_C;
    BASKER_MATRIX  &U = UBTF(c-btab);
    BASKER_MATRIX  &L = LBTF(c-btab);
    
    Int k = btf_tabs(c);
    Int j = M.col_ptr(k+1-bcol)-1;
    //Int j = M.row_idx[i];
    

    //will have to make a c and c'
    


    //printf("kid: %d chunk: %d col: %d bcol: %d j: %d\n",
    //	   kid, c, k, bcol, j);
    //printf("Single blk slv, kid: %d val:%f idx:%d \n",
    //	   kid, M.val[j], M.row_idx[j]);
    
    
    U.val(0)     = M.val(j);
    //M already has local idxing
    U.row_idx(0) = M.row_idx(j);
    L.val(0)     = (Entry) 1.0;
    L.row_idx(0) = M.row_idx(j);
    U.col_ptr(0) = 0;
    U.col_ptr(1) = 1;
    L.col_ptr(0) = 0;
    L.col_ptr(1) = 1;

   
    gperm(k) = k;
    gpermi(k)= k;

    
    return 0;
  }//end t_single_factor

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_blk_nfactor
  (
   Int kid, 
   Int c
   )
  {
    Int bcol = BTF_C.scol;
    Int brow = BTF_C.srow;
    Int btab = btf_tabs_offset;
    
    BASKER_MATRIX  &M = BTF_C;
    BASKER_MATRIX  &U = UBTF(c-btab);
    BASKER_MATRIX  &L = LBTF(c-btab);
    
    //JDB: brow hack: fix.
    Int brow2 = L.srow;

    /*
    printf("CURRENT BLK: %d \n", c);
    printf("col_start: %d \n", btf_tabs(c));
    printf("BLK size %d \n",
	   btf_tabs(c+1)-btf_tabs(c));
    */

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
    Entry absv, maxc;

    Int llnnz  = L.nnz;
    Int uunnz  = U.nnz;
    Entry maxv;

    Int maxindex;

    Int i,j, t;
    Int cu_ltop = 0;
    Int cu_utop = 0;
    
    Int newsize;

    unnz = 0;
    lnnz = 0;
    
    //printf("c: %d kid: %d ws_size: %d \n", c, kid, ws_size);

    //for each column
    for(Int k = btf_tabs(c); k < btf_tabs(c+1); ++k)
      {
	#ifdef BASKER_DEBUG_NFACTOR_DIAG
	printf("\n------------K=%d-------------\n", k);
	BASKER_ASSERT(top == ws_size, "nfactor dig, top");
	for( i = 0; i < ws_size; ++i)
	  {
	    if(X[i] != 0)
	      {
		printf("x(i) error: %d \n", i);
	      }
            if(ws[i] != 0)
              {
                printf("ws(i) error: %d \n", i);
              }
	    BASKER_ASSERT(X[i] == 0, "X!=0");
	    BASKER_ASSERT(ws[i] == 0, "ws!=0");
	  }
	#endif

	value = 0.0;
	pivot = 0.0;
	lcnt  = 0;
	ucnt  = 0;

	//for each nnz in column
	for( i = M.col_ptr(k-bcol); 
	     i < M.col_ptr(k-bcol+1); ++i)
	  {

	   

	    j = M.row_idx(i);
	    j = j - brow2;
	    
	    
            //printf("(One) In put j: %d %d \n", 
            //	   M.row_idx(i), j);

	    
	     if(j < 0)
	      {
		continue;
	      }


	    //X[j-brow] = M.val[i];
	    X(j) = M.val(i);
	    
	   
	    //printf("In put j: %d %d \n", 
            //	   M.row_idx(i), j);
	    
            #ifdef BASKER_DEBUG_NFACTOR_DIAG
	    printf("i: %d row: %d  val: %g  top: %d \n", 
                     i, j ,M.val[i], top);
	    printf("Nx in Ak %d %g %d color = %d \n",
                      j, X[j], brow,  
                     color[j] );
	    #endif


	    //if(color[j-brow] == 0)
	    if(color[j] == 0)
	      {
		//add custom local_reach
		t_local_reach(kid, L, 0, 0, j, &top);
	      }
	  }//end over all nnz in column
	xnnz = ws_size - top;

	#ifdef BASKER_DEBUG_NFACTOR_DIAG
	printf("xnnz: %d ws_size: %d top: %d \n", 
	       xnnz, ws_size, top);
        #endif
             
	//add custom local_solve
	t_back_solve(kid, L, 0,0,  k, top, xnnz);
	
	//Future add
	//t_locate_pivot(kid, top)	  
	//find pivot
	maxv = 0.0;
	for(i = top; i < ws_size; i++)
	  {
	    j = pattern[i];
	    //t = gperm[j];
	    //t = gperm(j+brow);
	    t = gperm(j+brow2);

	    //value = X[j-brow];
	    value = X(j);
	    
	    absv = abs(value);
	    //if(t == L.max_idx)
	    if(t == BASKER_MAX_IDX)
	      {
		++lcnt;
		if(absv > maxv) 
		  {
		    maxv     = absv;
		    pivot    = value;
		    maxindex = j;                
		  }
	      }
	  }//for (i = top; i < ws_size)
          //printf("b: %d lcnt: %d after \n", b, lcnt);

	  if(Options.no_pivot == BASKER_TRUE)
	    {
	      maxindex = k - brow2;
	      pivot    = X(maxindex);
	    }
	  
          ucnt = ws_size - top - lcnt +1;
         
	  if((maxindex == BASKER_MAX_IDX) || (pivot == 0))
            {
	      cout << endl << endl;
	      cout << "---------------------------"
		   <<endl;
              cout << "Error: Matrix is singular, blk"
		   << endl;
              cout << "MaxIndex: " << maxindex 
		   << " pivot " 
                   << pivot << endl;
              return 2;
            }          

          //printf("----------------PIVOT------------blk: %d %d \n", 
          //      c, btf_tabs(c+1)-btf_tabs(c));
	  gperm(maxindex+brow2) = k;
	  gpermi(k)             = maxindex+brow2;


	  //printf("TAG1 r  maxindex = %d k= %d \n", maxindex, k);
          #ifdef BASKER_DEBUG_NFACTOR_DIAG
          if((maxindex+brow) != k)
            {
              cout << "Permuting Pivot: " << k
		   << " as row " 
                   << maxindex << endl;
            }
          #endif
          
          //Note: Come back to this!!!!
          if(lnnz + lcnt > llnnz)
            {
	      printf("\n\n");
	      printf("----------------------\n");

              newsize = lnnz * 1.1 + 2 *A.nrow + 1;
              printf("Diag blk: %d Reallocing L oldsize: %d current: %d count: %d newsize: %d \n",
                     c,
		     llnnz, lnnz, lcnt, newsize);
	      printf("Columns in blks: %d %d %d \n",
		     btf_tabs(c), btf_tabs(c+1), 
		     (btf_tabs(c+1)-btf_tabs(c)));
            }
          if(unnz+ucnt > uunnz)
            {

	      printf("\n\n");
	      printf("-------------------\n");

              newsize = uunnz*1.1 + 2*A.nrow+1;
              printf("Diag blk: %d Reallocing U oldsize: %d newsize: %d \n",
                     c, uunnz, newsize);
            }

          //L.row_idx[lnnz] = maxindex;
	  L.row_idx(lnnz) = maxindex;
          //L.val[lnnz] = (Entry) 1.0;
	  L.val(lnnz)     = (Entry) 1.0;
          ++lnnz;
     
	  //printf("Pivot: %f \n", pivot);

          Entry lastU = (Entry) 0.0;
          for( i = top; i < ws_size; i++)
            {
	      j = pattern[i];
              //t = gperm[j];
	      //t = gperm(j+brow);
	      t = gperm(j+brow2);
            
              #ifdef BASKER_DEBUG_NFACTOR_DIAG
              printf("j: %d t: %d \n", j, t);
              #endif            

              //if fill-in
	      #ifdef BASKER_2DL
	      //if(X[j-brow] != 0)
	      if(X(j) != 0)
	      #else
              if(X[j] != 0)
	      #endif
                {
                  //if(t != L.max_idx)
		  if(t != BASKER_MAX_IDX)
                    {
                      if(t < k)
                        {
                          //U.row_idx[unnz] = gperm[j];
			  //U.row_idx(unnz)   = t-brow;
			  U.row_idx(unnz) = t - brow2;
			  #ifdef BASKER_2DL
			  //U.val[unnz] = X[j-brow];
			  U.val(unnz)       = X(j);
			  #else
                          U.val[unnz] = X[j];
			  #endif
                          ++unnz;
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

                      //L.row_idx[lnnz] = j;
		      L.row_idx(lnnz) = j;
		      #ifdef BASKER_2DL
		      //L.val[lnnz] = X[j-brow]/pivot;
		      L.val(lnnz) = X(j)/pivot;
		      #else
                      L.val[lnnz] = X[j]/pivot;
		      #endif
		      /*
		      printf("L(%d,%d) = %f \n",
			     lnnz,L.row_idx(lnnz),
			     L.val(lnnz));
		      */
                      ++lnnz;

                    }
                }//end if() not 0
              
              //Note: move x[j] inside of if() not 0..
	      //..extra ops this way
              #ifdef BASKER_DEBUG_NFACTOR_DIAG
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
          //U.row_idx[unnz] = k;
	  //U.row_idx(unnz)   = k-brow;
	  U.row_idx(unnz) = k - brow2;
          //U.val[unnz]     = lastU;
	  U.val(unnz)       = lastU;
          ++unnz;

          xnnz = 0;
          top = ws_size;
          
	  /*
	  
	  printf("Lcol(%d) = %d  \n",
		 k-brow2, cu_ltop);
	  printf("Lcol(%d) = %d \n",
		 k+1-brow2, lnnz);
	  */
		 
          //L.col_ptr[k-bcol] = cu_ltop;
	  L.col_ptr(k-brow2) = cu_ltop;
          //L.col_ptr[k+1-bcol] = lnnz;
	  L.col_ptr(k+1-brow2) = lnnz;
          cu_ltop = lnnz;
          
          //U.col_ptr[k-bcol] = cu_utop;
	  U.col_ptr(k-brow2)   = cu_utop;
          //U.col_ptr[k+1-bcol] = unnz;
	  U.col_ptr(k+1-brow2)      = unnz;
          cu_utop = unnz;


      }//over each column
    L.nnz = lnnz;
    U.nnz = unnz;

    L.col_ptr(L.ncol) = lnnz;
    U.col_ptr(U.ncol) = unnz;

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    //print_factor(L,U);
    #endif
    return 0;

  }//end t_factor_daig()

    
  //Note: need to fix indexing for color
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach
  (
   Int kid, BASKER_MATRIX &L,
   Int lvl, Int l,
   Int j, Int *top
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
        printf("stack_offset: %d head: %d \n", 
	       stack_offset , head);
        BASKER_ASSERT(head > -1,"local head val\n");
        #endif

	j = stack[head];
        //t = gperm[j];
	t = gperm(j+brow);
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("----------DFS: %d %d -------------\n",
	       j, t);
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
   BASKER_MATRIX &L,
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz
   )
  {
    
    INT_1DARRAY   ws = thread_array(kid).iws;
    ENTRY_1DARRAY  X = thread_array(kid).ews;
    Int      ws_size = thread_array(kid).iws_size;
    
    Int brow     = L.srow;

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
 
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;
    for(pp = 0; pp < xnnz; ++pp)
      {
	j = pattern[top1];

	#ifdef BASKER_2DL
	//color[j-brow] = 0;
	color[j] = 0;
	#else
	color[j] = 0;
	#endif
        //t = gperm[j];
	t = gperm(j+brow);
	
	if(t!=BASKER_MAX_IDX)
          {

	    #ifdef BASKER_2DL
	    //xj = X[j-brow];
	    xj = X(j);
	    #else
            xj = X[j];
	    #endif
	   
            //Get rid of these temp variables
            const Int local_offset = L.scol; //L.srow
            p2 = L.col_ptr(t+1-local_offset);
            p  = L.col_ptr(t-local_offset)+1;
            
            for( ; p < p2; ++p)
              {


		#ifdef BASKER_DEBUG_NFACTOR_DIAG
		#ifdef BASKER_2DL
		printf("Updating row: %d  with value: %f %f \n",
		       L.row_idx[p], X[L.row_idx[p]], 
		       L.val[p]*xj);
		#else
		printf("Updateing with value: %f %f \n",
		       X[L.row_idx[p]], L.val[p]*xj);
		#endif
		#endif
		#ifdef BASKER_2DL
		//X[L.row_idx[p]-brow] -= L.val[p]*xj;
		X(L.row_idx(p)) -= L.val(p)*xj;
		#else
                X[L.row_idx[p]] -= L.val[p] *xj;
		#endif
              }//end for() over each nnz in the column
          }//end if() not permuted
	++top1;
      }//end for() over all nnz in LHS
    return 0;
  }//end t_back_solve()


}//end namespace BaskerNS
#endif//end BASKER_NFACTOR_DIAG_HPP
