#ifndef SHYLUBASKER_NFACTOR_BLK_INC_HPP
#define SHYLUBASKER_NFACTOR_BLK_INC_HPP

#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain_inc_lvl
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space> *basker;

    kokkos_nfactor_domain_inc_lvl()
    {}

    kokkos_nfactor_domain_inc_lvl(Basker<Int,Entry,Exe_Space> *_basker)
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
      
      //if(kid ==0)
      {
	basker->t_nfactor_blk_inc_lvl(kid);
      }
    }//end operator
  };//end kokkos_nfactor_domain struct


  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain_remalloc_inc_lvl
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space> *basker;
    INT_1DARRAY                 thread_start;

    kokkos_nfactor_domain_remalloc_inc_lvl()
    {}

    kokkos_nfactor_domain_remalloc_inc_lvl
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

            
      {
	basker->t_nfactor_blk_inc_lvl(kid);
      }

	}//if-thread !BASKER_MAX_IDX
      
    }//end operator
  };//end kokkos_nfactor_domain_remalloc_inc_lvl struct
 
  //use local number on local blks (Crazy idea)
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_nfactor_blk_inc_lvl(Int kid)
  {
    Int b              = S(0)(kid); //Which blk from schedule
    BASKER_MATRIX &L   = LL(b)(0);
    BASKER_MATRIX &U   = LU(b)(LU_size(b)-1);
    BASKER_MATRIX &M   = ALM(b)(0); //A->blk

    INT_1DARRAY   ws   = LL(b)(0).iws;
    ENTRY_1DARRAY X    = LL(b)(0).ews;
    Int        ws_size = LL(b)(0).iws_size;

    Int          brow  = L.srow;  //begining row 
    Int          lval  = 0;
    Int          uval  = 0;

    Int i,j,k;
//    Int top, top1, maxindex, t;  //NDE - warning: top1 set but not used
    Int top, maxindex, t; 
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;
   
    Int newsize;
    Entry pivot, value;
    Entry absv, maxv;

    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    
    //Why did we need this?
    Int col_idx_offset = M.nnz;

    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("=======NFACTOR BLK INC LVL %d========\n", kid);
    #endif
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
    printf("b: %d scol: %d ecol: %d llnzz: %d uunzz: %d \n", 
           b, scol, ecol, L.nnz, U.nnz);
    #endif

    //Debug time
    //Entry dfs_time = 0;
    
    //for each column    
    for(k = 0; k < M.ncol; ++k)
        {
          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("\n----------------K=%d %d--------------\n", k+M.scol, kid);
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
          for(i = 0 ; i < ws_size; i++){ASSERT(X[i] == 0);}
          //ASSERT int workspace is clean
	  for(i = 0; i <  ws_size; i++){ASSERT(ws[i] == 0 );}
          #endif

	 
          //for each nnz in column
	  //Wnat to change this to local blk anyway
	  //Backwards is a little faster because less flips
	  //for(i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)

	  //Debug of dim
	  //Kokkos::Impl::Timer timer;
	  if(Options.same_pattern == BASKER_FALSE)
	    {
	  for(i = M.col_ptr(k+1)-1; i >= M.col_ptr(k); --i)
	    {
	      j = M.row_idx(i);
	      X(j) = M.val(i);
	           
        #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	      printf("i: %d row: %d  val: %g  top: %d \n", 
		     i, j ,M.val(i), top);
	      printf("Nx in Ak %d %g %d color = %d \n",
		     j, X[j], brow,  
                     color[j] );
        #endif

	 

	      //Note we need to check all 
	      //LVL = min max{}+1
	      //Need to search all "posible" paths
	      if(Options.incomplete_type ==
		 BASKER_INCOMPLETE_LVL)
		{
		  //printf("reach one called\n");
		  t_local_reach_inc_lvl(kid,0,0,j,&top);
		}
	      else
		{
		  //printf("reach two called \n");
		  if(gperm(j+brow) != BASKER_MAX_IDX)
		    {
		      t_local_reach_inc_rlvl(kid,0,0,j,top);
		    }
		  else
		    {
		      t_local_reach_short_inc_rlvl(kid,0,
						   0,j,top);
		    }
		}
	    
          }//end for() each nnz in column
	  }//if--same_pattern
	  else
	    {

	      //Populate X(j) with values of New A
	      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); i++)
		{
		  j = M.row_idx(i);
		  X(j) = M.val(i);
		}

	      //Get L
	      for(Int i = L.col_ptr(k+1)-1; i >= L.col_ptr(k);
		  i--)
		{
		  j = L.row_idx(i);
		  color[j] = 2;
		  
		  //Lower does not matter if pivot or not
		  pattern[--top] = j;  
		}

	      //Get U pattern
	      for(Int i = U.col_ptr(k+1)-2; i >= U.col_ptr(k);
		  i--)
		{
		  
		  j = U.row_idx(i);
		  color[j] = 2;
		 
		  if(Options.no_pivot == BASKER_TRUE)
		    {
		      pattern[--top] = j;
		      //Note we do not store this for domains
		      //INC_TEMP_LVL(j) = U.inc_lvl(i);
		    }
		  else
		    {
		      BASKER_ASSERT(0==1, 
				    "Currently not supported");
		      pattern[--top] = j;
		    }
		}
	      

	    }//end --- same pattern
	  xnnz = ws_size - top;
	  //dfs_time += timer.seconds();

          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("xnnz: %d ws_size: %d top: %d \n", 
                 xnnz, ws_size, top);
          #endif
	 
          
	  Entry relax_value = 0;
	  if(Options.incomplete_type == 
	     BASKER_INCOMPLETE_LVL)
	    {
	      //printf("backsolve one called\n");
	      t_back_solve_inc_lvl(kid,0,0,k,top,xnnz);
	    }
	  else
	    {
	      //printf("backsolve two called \n");
	      t_back_solve_inc_rlvl(kid,0,0,k,top,xnnz,
				    relax_value);
	    }


	  //Select Pivot 
	  //Get upper and lower count
	  //Note: upper and lower count should be moved to
	  //back_solve to speed this up when not pivoting!!!
	  //Most, likely case first, 
	  //Which is no_pivot in incomplete factorization
          maxv = 0.0;
	  ucnt = 0;
	  for(i = top; i < ws_size; i++)
	    {
	      j = pattern[i];
	      t = gperm(j+brow);

	      //printf("Consider: %d %d inc: %d \n",
	      //     j, t, INC_LVL_TEMP(j+brow));

	      //Only consider those that can be included
	      if(INC_LVL_TEMP(j+brow) <= Options.inc_lvl)
		{
		  
		  value = X(j);
		  absv = EntryOP::approxABS(value);
		  if(t == BASKER_MAX_IDX)
		    {
		      lcnt++;
		      //We can short circuit the gt 
		      //so it is not used all the time
		      if((Options.no_pivot != BASKER_TRUE) &&
			 EntryOP::gt(absv,maxv))
			{
			  maxv     = absv;
			  pivot    = value;
			  maxindex = j;                
			}
		    }//if lower-half of column
		  else
		    {
		      ucnt++;
		    }
		}//if lvl < option
	    }//for over all elements

	  if(Options.no_pivot == BASKER_TRUE)
	    {
	      maxindex = k;
	      pivot    = X(k);
	    }
          //ucnt = ws_size - top - lcnt +1;
    
    if((maxindex == BASKER_MAX_IDX) || (pivot == (Entry)(0)) )
            {
	      if(Options.verbose == BASKER_TRUE)
		{
	      cout << endl << endl;
	      cout << "---------------------------"<<endl;
              cout << "Error: Matrix is singular, blk" << endl;
              cout << "MaxIndex: " << maxindex << " pivot " 
                   << pivot << endl;
              cout << "lcnt: " << lcnt << endl;
		}
	      thread_array(kid).error_type =
		BASKER_ERROR_SINGULAR;
	      thread_array(kid).error_blk  = b;
	      thread_array(kid).error_info = k;
	      return BASKER_ERROR;
            }          

	  gperm(maxindex+brow) = k+brow;
	  gpermi(k+brow) = maxindex + brow;
     
          //Note: Come back to this!!!!
          if(lnnz + lcnt > llnnz)
            {
       
              newsize = lnnz * 1.1 + 2 *M.nrow + 1;

	      if(Options.verbose == BASKER_TRUE)
		{
	      printf("\n\n");
	      printf("----------------------\n");
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

              newsize = uunnz*1.1 + 2*M.nrow+1;
	      
	      if(Options.verbose == BASKER_TRUE)
		{
	      printf("\n\n");
	      printf("-------------------\n");

              printf("b: %ld Reallocing U oldsize: %ld newsize: %ld \n",
                     (long)b, (long)uunnz, (long)newsize);
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
		  thread_array(kid).error_info  = newsize;
		  return BASKER_ERROR;
		}

            }
	  
	
          L.row_idx(lnnz) = maxindex;
          L.val(lnnz)     = (Entry) 1.0;
	  if(Options.same_pattern == BASKER_FALSE)
	    {
	      L.inc_lvl(lnnz) = INC_LVL_TEMP(maxindex+brow);
	      INC_LVL_TEMP(maxindex+brow) = BASKER_MAX_IDX;
	    }
	  #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	  printf("add L(%d): 1.0 lvl: %d \n", 
		 maxindex, L.inc_lvl(lnnz));
	  #endif
          lnnz++;
     
          Entry lastU = (Entry) 0.0;
          for( i = top; i < ws_size; i++)
            {

	      j = pattern[i];
	      pattern[i] = 0;
	      t = gperm(j+brow);
            
              #ifdef BASKER_DEBUG_NFACTOR_BLK
              printf("j: %d t: %d lvl: %d\n", 
		     j, t, INC_LVL_TEMP(j+brow));
              #endif         
	      
	      if((Options.same_pattern == BASKER_TRUE) ||
		 (INC_LVL_TEMP(j+brow) <= Options.inc_lvl))
	      {
		if(t != BASKER_MAX_IDX)
		  {
		    if(t < (k+brow))
		      {
			U.row_idx(unnz) = t-brow;
			U.val(unnz) = X(j);
			unnz++;
    		
			#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
			printf("add U(%d): %g lvl: %d \n",
			       U.row_idx(unnz-1),
			       U.val(unnz-1), 
			       INC_LVL_TEMP(j+brow));
			#endif
		      }
		    else
		      {
			lastU = X(j);
		      }
		  }
		  else if (t == BASKER_MAX_IDX)
                    {
		      L.row_idx(lnnz) = j;
		      L.val(lnnz) = EntryOP::divide(X(j),pivot);
		     
		      #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		      printf("add L(%d): %g lvl: %d \n",
			     j, L.val(lnnz), 
			     INC_LVL_TEMP(j+brow));
		      #endif
		      
		      if(Options.same_pattern == BASKER_FALSE)
			{
			  L.inc_lvl(lnnz) = INC_LVL_TEMP(j+brow);
			  INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
			}
		      
                      lnnz++;
                    }
                }//end if() lvl <
	      else
		{
		  INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
		  //Add relaxation here!
		}

              #ifdef BASKER_DEBUG_NFACTOR_BLK
	      printf("Zeroing element: %d \n", j);
              #endif

	      X(j) = 0;	  
            }//end if(x[i] != 0)

          //Fill in last element of U
	  U.row_idx(unnz) = k;
	  U.val(unnz)     = lastU;
          unnz++;

	  //printf("Add U(%d) %g \n",
	  //	 k, U.val(unnz-1));

          xnnz = 0;
          top = ws_size;
          
	  L.col_ptr(k)   = cu_ltop;
	  L.col_ptr(k+1) = lnnz;
          cu_ltop        = lnnz;
          
	  U.col_ptr(k)   = cu_utop;
	  U.col_ptr(k+1) = unnz;
          cu_utop        = unnz;

	  
	  #ifdef BASKER_2DL
	  //-----------------------Update offdiag-------------//
	  for(Int blk_row = 1; blk_row < LL_size(b); ++blk_row)
	    {
	      //Do back solve of off-diag blocks
	      // printf("before offdiag blk_row: %d kid: %d \n", 
	      //     blk_row, kid);
              
              
	      if(Options.same_pattern == BASKER_FALSE)
		{

		
              t_dom_lower_col_offdiag_find_fill(kid, L.srow,
                                                b, blk_row,
                                                b, blk_row,
                                                k,
                                                U.row_idx,
                                  U.col_ptr(k+1)-U.col_ptr(k),
                                                U.col_ptr(k),
                                                BASKER_TRUE);
              
	      
	     
	      
	      t_back_solve_offdiag_inc_lvl(kid, L.srow,
				   b, blk_row,
				   b, blk_row,
				   k, col_idx_offset,
				   U.val, U.row_idx,
		       U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k),
				   BASKER_TRUE);
	      
		}
	      else
		{
		  
	
		  t_back_solve_offdiag_same_pattern_inc_lvl(kid,
						  L.srow,
						  b,blk_row,
						  b,blk_row,
					      k, col_idx_offset,
 	      				    U.val, U.row_idx,
				    U.col_ptr(k+1)-U.col_ptr(k),
						  U.col_ptr(k),
						   BASKER_TRUE);


		}
	      //Move these factors into Local Ls
	      Int move_error = 
	      t_move_offdiag_L_inc_lvl(kid,
			       b, blk_row,
			       b, blk_row,
			       k, pivot);
	      
	      if(move_error == BASKER_ERROR)
		{
		  return BASKER_ERROR;
		}
		  
	    }//end over all diag
	 
	  //This can be made faster by only looping over used
	  //Note that this will need to be changed with pivoting
	 
	  if(Options.same_pattern == BASKER_FALSE)
	    {
	      for(Int i = U.col_ptr(k); i < U.col_ptr(k+1); i++)
		{
		  const Int j = U.row_idx(i);
		  INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
		}
	    }
	    
	  #endif


	}//end for() over all columns

    L.nnz = lnnz;
    U.nnz = unnz;

    //Debug//
    //printf("DFS time: %f \n", dfs_time);

    return 0;
  }//end t_nfactor_blk_inc_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void
  Basker<Int,Entry,Exe_Space>::t_local_reach_inc_rlvl
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
//    Int start, end, done; //NDE - warning: done set but unused
    Int start, end;


    Int inc_lvl = 0;
    
    //printf("reach blks: %d %d kid: %d  lvl: %d %d  \n",
    //	   b,0, kid, lvl, l);

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d j: %d %d, kid: %d \n",
	   b, 0, wsb, l, j+brow, brow, kid);
    #endif

    start    = -1;
    head     = 0;
    stack[0] = j;

    if(lvl!=0)
      {
	//printf("kid: %d j: %d  using inc: %d %d\n",
	//		   kid , j+brow, INC_LVL_TEMP(j+brow),
	//		   INC_LVL_TEMP(j+brow)+inc_lvl);
	inc_lvl =  INC_LVL_TEMP(j+brow);
      }
   


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
	printf("--------DFS: %d %d %d kid: %d -------------\n",
	       j, j+brow, t, kid);
        #endif

	if((ws(j) == 0) || (ws(j) == 2))
	  {	    

	    if(ws(j) == 0)
	      {
		ws(j) = 1;
	      }
	    else
	      {
		ws(j) = 3;
	      }
  
	    if(t!=BASKER_MAX_IDX)
	      {
		start = L.col_ptr(t+1-L.scol);
	      }//if not permuted
	    //printf("kid: %d j: %d start: %d inc_lvl: %d\n",
	    //	   kid, j+brow, start, inc_lvl);
	  }
	else
	  {
	    start = store[j];
	    inc_lvl = inc_lvl - L.inc_lvl(start) -1;

	    //printf("kid: %d j: %d start2: %d inc_lvl: %d \n",
	    //	   kid, j+brow, start, inc_lvl);

	    BASKER_ASSERT(inc_lvl >= 0, 
			  "inc_lvl too small 2");
	  }
//	done = 1; //NDE - warning: set but unused
	

	if(inc_lvl <= Options.inc_lvl)
	  {
	    end = L.col_ptr(t-L.scol);
	  }
	else
	  {
	    end = start;
	  }
	//printf("end: %d \n", end);
	for(i1 = --start; i1 >= end; --i1)
	  {
	    i = L.row_idx(i1);
	   
	    ///printf("Search i1: %d  i: %d %d %d lvl: %d color: %d\n",
	    //	   i1, i, i+L.scol, gperm(i+L.scol), 
	    //	   INC_LVL_TEMP(i+brow),ws(i));
	    
	    //if explored
	    if((ws(i) != 0))
	      {

		//if ws(i) ==1 || 3, we already seen it
		// and the path must be less
		

		//if ws(i) ==2, we might have not
		//printf("continue called\n");
		if(ws(i) == 2)
		  {
		if(L.inc_lvl(i1)+inc_lvl < Options.inc_lvl)
		  {
		    if(gperm(i+brow) != BASKER_MAX_IDX)
		      {
			store[j] = i1;
			stack[++head] = i;
			inc_lvl = inc_lvl + L.inc_lvl(i1)+1;
			//printf("kid: %d j: %d seen but exploring inc: %d\n",
			//     kid, i+brow, inc_lvl);
			break;
		      }
		    else
		      {
		
			INC_LVL_TEMP(i+brow) =
			  min(inc_lvl+L.inc_lvl(i1)+1,
			      INC_LVL_TEMP(i+brow));
		
		      


		      }
		  }
		else
		  {
		    continue;
		  }
		  }
	      }
	    else
	    //if(ws(i) == 0)
	      {
		//if(not explored)
		if(gperm(i+brow) >= 0)
		  {
		    //store[j] = i1+1;
		    store[j] = i1;
		    stack[++head] = i;
		    inc_lvl = inc_lvl + L.inc_lvl(i1) + 1;
		    //printf("kid: %d j: %d  break inc: %d \n",
		    //	   kid, j, inc_lvl);
		    break;
		  }
		else
		  {
		    //color[i] = 2;
		    ws(i) = 2;
		    pattern[--top] = i;

		    if(i == 4)
		      {
			//printf("\n\n\n Add 4 leaf \n \n");
		      }


		    if(INC_LVL_TEMP(i+brow) == BASKER_MAX_IDX)
		      {
			INC_LVL_TEMP(i+brow) = 
			  inc_lvl + L.inc_lvl(i1)+1;
			//printf("kid: %d setone: %d \n",
			//   kid, INC_LVL_TEMP(i+brow));
		      }
		    else
		      {
			INC_LVL_TEMP(i+brow) =
			  min(inc_lvl+L.inc_lvl(i1)+1,
			      INC_LVL_TEMP(i+brow));
			//printf("kid: %d settwo: %d \n",
			       //kid, INC_LVL_TEMP(i+brow));
		      }

		    //printf("Adding idx: %d %d  to pattern at location: %d inc_lvl: %d  \n",i, i+L.scol, top, INC_LVL_TEMP(i+brow));

		  }
	      }
	   


	  }
	
	//if we have used up the whole column
	if(i1 < end)
	  {
	    
	       if(INC_LVL_TEMP(j+brow) == BASKER_MAX_IDX)
	      {
		INC_LVL_TEMP(j+brow) = inc_lvl;
		//printf("kid: %d j: %d set three: %d \n",
		//   kid, j+brow, inc_lvl);
	      }
	    else
	      {
		INC_LVL_TEMP(j+brow) =
		  min(inc_lvl, INC_LVL_TEMP(j+brow));
		//printf("kid: %d j: %d set four: %d \n",
		//   kid, j+brow, INC_LVL_TEMP(j+brow));
	      }

	       //printf("consider adding: %d color %d lvl: %d \n",
	       //  j+brow, ws(j), INC_LVL_TEMP(j+brow));

	    head--;
	    if((ws(j) != 2) &&
	       (ws(j) != 3))
	      {
		//color[j] = 2;
		ws(j) = 2;
		pattern[--top] = j;
		//printf("add: %d \n",
		//j+brow);

		//if(j == 4)
		// {
		    // printf("\n\n\n ADD 4 root \n\n");
		//}


	      }
	    if(ws(j) == 3)
	      {
		ws(j) = 2;
	      }
	       
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

  }//end t_local_reach_inc_rlvl

  template <class Int, class Entry, class Exe_Space>
  inline
  void 
  Basker<Int,Entry,Exe_Space>::t_local_reach_short_inc_rlvl
  (
   const Int kid, 
   const Int lvl,
   const Int l,
   const Int j,
   Int &top
   )
  {
    //Setup variables
    const Int b      = S(lvl)(kid);
    const Int wsb    = S(0)(kid);
    BASKER_MATRIX &L = LL(b)(0);

    INT_1DARRAY  ws   = LL(wsb)(l).iws;
    const Int ws_size = LL(wsb)(l).iws_size;

    Int *color       = &(ws(0));
    Int *pattern     = &(ws(ws_size));


    if(color[j] != 2)
      {
    color[j]       = 2;
    pattern[--top] = j;
      }
    if(Options.incomplete == BASKER_TRUE)
      {
	//printf("short reach: j: %d inc: %d \n",
	//     j+L.srow, INC_LVL_TEMP(j+L.srow));
	if(INC_LVL_TEMP(j+L.srow) == BASKER_MAX_IDX)
	  {
	    INC_LVL_TEMP(j+L.srow) = 0;
	  }//if incomplete lvl not set
      }//end--if incomplete

    //printf("kid: %d short_reach: %d \n",
    //	   kid, j);
    return;
  }//end t_local_reach short_inc_rlvl()
  
  //Note: need to fix indexing for color
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int 
  Basker<Int,Entry,Exe_Space>::t_local_reach_inc_lvl
  (
   Int kid, Int lvl, Int l,
   Int j, Int *top
   )
  {

    //Notes: 
    //We cannot use symmetric pruning with our 
    //incomplete level appoarch.
    //Moreover, we use forward scan and not backward
    //Will want to make this backward in the future

    //Setup variables
    const Int      b   = S(lvl)(kid);
    const Int     wsb  = S(0)(kid);
    BASKER_MATRIX  &L  = LL(b)(0);
    const Int     brow = L.srow;

    INT_1DARRAY    ws  = LL(wsb)(l).iws;
    const Int  ws_size = LL(wsb)(l).iws_size;
 
    Int *color       = &(ws(0));
    Int *pattern     = &(ws(ws_size));
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);

    Int i, t, head, i1;
    Int start, end;

    BASKER_BOOL done;
    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach_inc_lvl, L: %d %d  X: %d %d, kid: %d \n",
	   b, 0, wsb, l, kid);
    #endif

    start    = -1;
    head     = 0;
    stack[0] = j;

    Int inc_lvl = 0;
    Int pop_top = *top-1; //Where we want to add


      if(lvl!=0)
      {
	//printf("kid: %d j: %d  using inc: %d %d\n",
	//		   kid , j+brow, INC_LVL_TEMP(j+brow),
	//		   INC_LVL_TEMP(j+brow)+inc_lvl);
	inc_lvl =  INC_LVL_TEMP(j+brow);
      }
   


    //Because we cut short the tree
    //We nee to find a merge point
    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("begin sample top: %d p:%d %d\n",
	   pop_top, pattern[pop_top], pattern[pop_top+1] );
    #endif
    while(((pop_top+1)!=ws_size) && (pattern[pop_top+1] < j))
      {
        #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	printf("sample top: %d p:%d \n",
	       pop_top, pattern[pop_top]);
	#endif
	pop_top++;
      }
    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    if(kid == 0)
      {
    printf("end sample top: %d p:%d %d \n",
	   pop_top, pattern[pop_top], pattern[pop_top+1]);
      }
    #endif

    //printf("short test: %d %d kid: %d \n",
    //	   j+brow, gperm(j+brow), kid);
    //Short circuit
    if(gperm(j+brow) == BASKER_MAX_IDX)
      {
        #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	printf("==DFS Short circuit: %d %d== \n",
	       j, gperm(j+brow));
	#endif

	//If we have not seen 
	if(color[j] == 0)
	  {
	    color[j] = 2;
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("test top: %d pop_top: %d \n",
		   *top-1, pop_top-1);
	    #endif
	    for(i = (*top-1);  i < (pop_top); i++)
	      {
		#ifdef BASKER_NFACTOR_BLK_INC
		printf("moving: %d %d to %d %d\n",
		       i+1, pattern[i+1],
		       i,  pattern[i]);
		#endif
		pattern[i] = pattern[i+1];
	      }
	    (*top)--;
	    pattern[pop_top] = j;
	    pop_top--;
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("add j: %d to pattern at %d \n",
		   j ,  pop_top);
	    printf("xsize: %d \n",
		   ws_size - *top);
	    #endif
	  }
       
	if(INC_LVL_TEMP(j+L.srow) == BASKER_MAX_IDX)
	  {
	    INC_LVL_TEMP(j+L.srow) = 0;
	  }
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	printf("===========LEAVING REACH ======\n");
	printf("Leave top: %d leave pop: %d \n",
	   *top, pop_top);
	#endif
	return 0;
      }//end if gperm short circuit

    


    while(head != BASKER_MAX_IDX)
      { 
        #ifdef BASKER_DEBUG_LOCAL_REACH
        printf("stack_offset: %d head: %d \n", 
	       stack_offset , head);
        BASKER_ASSERT(head > -1, "head !> -1 ");
        #endif

	j = stack[head];
        t = gperm(j+brow);
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("--------DFS: %d %d %d %d-----------\n", 
	       j, t, color[j], inc_lvl);
        #endif

	if((ws(j) == 0)||(ws(j)==2))
	  {

	    //Note:
	    //Use 4 colors
	    // 0 - new
	    // 1 - iter on new
	    // 2 - found
	    // 3 - found iter

	    if(ws(j) == 0)
	      {
		ws(j) = 1;
	      }
	    else
	      {
		ws(j) = 3;
	      }
  

	    if(t!=BASKER_MAX_IDX)
	      {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("reach j: %d t:%d L.scol%d\n",
		       j, t, L.scol);
                #endif
		
		start = L.col_ptr(j)+1;
	      }
            else
              {
                #ifdef BASKER_DEBUG_LOCAL_REACH
                printf("L.scol: %d  L.ncol: %d t: %d \n",
		       L.scol, L.ncol,t);
                #endif
              }
	  }
	else
	  {
	    if(t != BASKER_MAX_IDX)
	      {
		start = store[j];
                #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		printf("subtracting: %d %d %d \n", 
		       inc_lvl, L.inc_lvl(start-1),
		       inc_lvl-L.inc_lvl(start-1)-1);
		#endif
		inc_lvl = inc_lvl-L.inc_lvl(start-1)-1;
		
		BASKER_ASSERT(inc_lvl >=0, "inc_lvl too small\n");
	      }
	  }
	done = BASKER_TRUE;

	#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	printf("t: %d %d inc_lvl: %d \n",
	       t, BASKER_MAX_IDX, inc_lvl);
	#endif
	if(t!=BASKER_MAX_IDX)
	  {
	    end = L.col_ptr(t+1-L.scol);
	    store[j] = start;
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("inc_lvl = %d options=%d \n",
		   inc_lvl, Options.inc_lvl);
	    #endif
	    if((inc_lvl < Options.inc_lvl))
	      {
		for(i1 = start; i1 < end; ++i1)
		  {
		    i = L.row_idx(i1);
		    
		    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		    printf("considering i: %d color: %d touched: %d  \n",
			   i, color[i], L.inc_lvl(i1));
		    printf("considering2 i: %d color: %d inc_lvl: %d \n",
			       i, L.inc_lvl(i1), inc_lvl);
		    #endif
		   
		    if((L.inc_lvl(i1)+inc_lvl) < Options.inc_lvl)
		    {
		
			  { 
			    head++;
			    stack[head] = i;
			    store[j] = i1+1;
			    done = BASKER_FALSE;
		    
			    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
			    printf("additng: %d %d %d \n",
				   inc_lvl, L.inc_lvl(i1),
				   inc_lvl+L.inc_lvl(i1)+1);
			    #endif
			   
			    inc_lvl = inc_lvl + L.inc_lvl(i1)+1;
			   
			    break;
			    
			  }

		    }//check on fill level
		  }//over all nnz in the column
	      }//check fill level
	  }//end ifend not perm
	
	if(done == BASKER_TRUE)
	  {

	    if((color[j] == 0)||(color[j] == 1))
	      {
	    
		BASKER_ASSERT((*top-1) >= 0, "Top pass pattern\n");
		if((*top-1)<0)
		  {
		printf("blk: %ld kid: %ld \n", (long)b, (long)kid);
		  }

		color[j] = 2;
	       		 
	     //loop forward
	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("begin refine top: %d p:%d %d \n",
		    pop_top, pattern[pop_top], pattern[pop_top+1]);
	     #endif
	     while((pop_top+1!=ws_size)&&(pattern[pop_top+1] < j))
	       {
		 #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		 printf("sample refine top: %d p:%d %d \n",
			pop_top, pattern[pop_top], pattern[pop_top+1]);
		 #endif
		 pop_top++;
	       }
	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("end refine top: %d p:%d \n",
		    pop_top, pattern[pop_top]);
	     #endif

	    
	     //Note that this can be done with a pop
	     //In future rewrites
	     //loop backwards
	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("begin bb refine top: %d p:%d %d \n",
		    pop_top, pattern[pop_top], pattern[pop_top-1]);
	     #endif
	     while((pop_top-1!=ws_size)&&(pattern[pop_top] > j) )
	       {
		 #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		 printf("sample bb refine top: %d p:%d %d \n",
			pop_top, pattern[pop_top], pattern[pop_top-1]);
		 #endif
		 pop_top--;
	       }
	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("end bb refine top: %d p:%d %d \n",
		    pop_top, pattern[pop_top], pattern[pop_top-1]);
	     #endif
	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("test top: %d pop_top: %d \n",
		   *top-1, pop_top-1);
	     #endif
	    for(i = (*top-1);  i < (pop_top); i++)
	      {
		#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		printf("moving: %d %d to %d %d\n",
		       i+1, pattern[i+1],
		       i,  pattern[i]);
		#endif
		pattern[i] = pattern[i+1];
	      }
	    (*top)--;

	
	    pattern[pop_top] = j;
	    pop_top--;
	
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("Adding idx: %d to pattern at location: %d \n",j, pop_top+1);
	    printf("total size: %d \n",
		   ws_size-*top);
	    #endif
	      }//if colored
	    if(head == 0)
	      {head = BASKER_MAX_IDX;}
	    else
	      {head--;}
	  
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("compare(%d) inc_lvl: %d INC_LVL_TEMP: %d \n",
		   j+brow, inc_lvl, INC_LVL_TEMP(j+brow));
	    #endif
	    
	    if(INC_LVL_TEMP(j+brow) == BASKER_MAX_IDX)
	      {
		INC_LVL_TEMP(j+brow) = inc_lvl;
	      }
	    else
	      {
		INC_LVL_TEMP(j+brow) =
		  min(inc_lvl, INC_LVL_TEMP(j+brow));
	      }
	    
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("setting min level: %d to %d \n",
		   j+brow, INC_LVL_TEMP(j+brow));
	    #endif
	    
	    color[j] = 2;

	  }//if done
      }//end while head

    
    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("===========LEAVING REACH ======\n");
    printf("Leave top: %d leave pop: %d \n",
	   *top, pop_top);
    #endif
    return 0;
  }//end t_local_reach_inv_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int 
  Basker<Int, Entry,Exe_Space>::t_back_solve_inc_lvl
  (Int kid, 
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz
   )
  {

    //We note that this can be fixed to be faster

    const Int       b = S(lvl)(kid);
    const Int     wsb = S(0)(kid);
    BASKER_MATRIX  &L = LL(b)(0);
    INT_1DARRAY    ws = LL(wsb)(l).iws;
    ENTRY_1DARRAY   X = LL(wsb)(l).ews;
    const Int ws_size = LL(wsb)(l).iws_size;

    Int brow = L.srow;

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
 
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;
    //printf("back solve: nnz: %d \n",
    //	   xnnz);
    for(pp = 0; pp < xnnz; pp++)
      {
	
	j = pattern[top1];
	t = gperm(j+brow);

	//printf("backsolve column: %d %d\n",
	//  j,t);

	if(t!= BASKER_MAX_IDX)
          {
	    xj = X(j);
	   
            //Get rid of these temp variables
            Int local_offset = L.scol;
            p2 = L.col_ptr(t+1-local_offset);
            p = L.col_ptr(t-local_offset)+1;
            
            for( ; p < p2; p++)
              {
		if(color[L.row_idx(p)] == 0)
		  {
		    //printf("Continue color: %d \n",
		    //   L.row_idx(p));
		    continue;
		  }
	
		
		//if(L.row_idx(p) == 96)
		// {
		//printf("Updateing row: %d with value: %f %f \n",
		// L.row_idx(p),X(L.row_idx(p)), L.val(p)*xj);
		// }

                X(L.row_idx(p)) -= L.val(p) *xj;

              }//end for() over each nnz in the column
          }//end if() not permuted
	top1++;
      }//end for() over all nnz in LHS

    //Now we can clear the colors
    top1 = top;
    for(pp = 0; pp < xnnz; pp++)
      {

	j = pattern[top1];
	color[j] = 0;
	top1++;
       
      }//end for over all nnz_pattern


    return 0;
  }//end t_back_solve_inc_lvl()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int 
  Basker<Int, Entry,Exe_Space>::t_back_solve_inc_rlvl
  (Int kid, 
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz,
   Entry &relax_value
   )
  {

    //We note that this can be fixed to be faster
    const Int       b = S(lvl)(kid);
    const Int     wsb = S(0)(kid);
    BASKER_MATRIX  &L = LL(b)(0);
    INT_1DARRAY    ws = LL(wsb)(l).iws;
    ENTRY_1DARRAY   X = LL(wsb)(l).ews;
    const Int ws_size = LL(wsb)(l).iws_size;

    Int brow     = L.srow;
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
 
    Int top1 = top;
    //printf("back solve top: %d %d \n", top, top1);
    Int j,t,pp, p, p2;
    Entry xj = 0;
    for(pp = 0; pp < xnnz; pp++)
      {
	j = pattern[top1];
	t = gperm(j+brow);

	//printf("kid: %d backsolve %d %d \n",
	// kid, j, t);

	//if(Options.incomplete_type ==
	//  BASKER_INCOMPLETE_RLVL_LIMITED)
	  {
	    //Going to short circuit
	    //If not the right level
	    //printf("kid: %d limite, j: %d inc_lvl: %d \n",
	    //	   kid, j+brow, INC_LVL_TEMP(j+brow));
	    if(INC_LVL_TEMP(j+brow) > Options.inc_lvl)
	      {
		
		//printf("limited continue j: %d kid: %d \n", 
		//   j, kid);
		top1++;
		continue;
	      }
	  }//end if- rlvl limited

	if(t!= BASKER_MAX_IDX)
          {
	    xj = X(j);
	    //Get rid of these temp variables
            Int local_offset = L.scol;
            p2 = L.col_ptr(t+1-local_offset);
            p = L.col_ptr(t-local_offset)+1;
            
            for( ; p < p2; p++)
              {
		if(color[L.row_idx(p)] == 0)
		  {
		    // printf("Color continue %d \n", 
		    //	   L.row_idx(p));
		    //top1++;
		    continue;
		  }
		
		/*
		if((L.row_idx(p) + L.srow) == 93)
		  {
		    
		printf("j: %d X(%d): %f Lval: %f xj: %f end: %f \n",
		       t,
		       L.row_idx(p),
		       X(L.row_idx(p)),
		       L.val(p),
		       xj, 
		       X(L.row_idx(p))-L.val(p)*xj);
		  }
		*/
		
		

                X(L.row_idx(p)) -= L.val(p) *xj;

              }//end for() over each nnz in the column
          }//end if() not permuted
	//printf("end of not permuted top: %d\n", top1);
	top1++;
      }//end for() over all nnz in LHS

    //Now we can clear the colors
    //We find this to be slightly faster than 
    //clearing the color in the compute loop for some reason
    top1 = top;
    for(pp = 0; pp < xnnz; pp++)
      {
	j = pattern[top1];
	color[j] = 0;
	//printf("CLEARING COLOR: %d \n", 
	//     j);
	top1++;
	
      }//end for over all nnz_pattern

    return 0;
  }//end t_back_solve_inc_rlvl()

  

  template <class Int, class Entry, class Exe_Space>
  void
  Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag_same_pattern_inc_lvl
 (Int kid, Int pbrow,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    BASKER_MATRIX &B     = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int         ws_size  = LL(X_col)(X_row).iws_size;
    
    Int          nnz     = LL(X_col)(X_row).p_size;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
//    Int *stack   = &(pattern[ws_size]); //NDE - warning: unused
   
    //Fill Colog with correct index marks from already L
    for(Int i = L.col_ptr(k); i < L.col_ptr(k+1); i++)
      {
	const Int j = L.row_idx(i);
	color[j] = 1;
	//Copying the pattern is over kill, 
	//Should not hurt that much though
	//printf("kid: %d add pattern [%d %d] \n",
	//     kid, k+L.scol, j+L.srow);
	pattern[nnz++] = j; 
      }
    

    if(A_option == BASKER_TRUE)
      {
	#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	printf("k: %d size: %d \n",
	       k , B.col_ptr(k+1)-B.col_ptr(k));
	#endif

	for(Int i = B.col_ptr(k); 
	    i < B.col_ptr(k+1); i++)
	  {

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	    printf("t_back_solve_diag, kid: %d inc[%d] = 0\n",
		   kid, B.row_idx(i));
	    #endif
	    
	    //printf("Add to X [%d %d] %f kid: %d \n",
	    //	   k+L.scol, B.row_idx(i)+L.srow,
	    //	   B.val(i),
	    //	   kid);

	    const Int j = B.row_idx(i);
	    X(j) = B.val(i);
	  }//over all nnz in subview
      }//end if preload

    //printf("done with reload\n\n");
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %d \n",
	   kid, x_size);
    #endif
    for(Int i = 0 ; i < x_size-1; ++i)
      {
	
	const Int k    = x_idx[i+x_offset];
	const Entry xj = x(i+x_offset);
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	       kid, k+pbrow, L.col_ptr[k], L.col_ptr[k+1]);
	#endif

        #ifdef BASKER_DEBUG_NFACTOR_INC_LVL
	//A multiple of a value at lvl l results in l+1 fill-in
	printf("LVL_TEMP[%d] = %d, %d kid: %d continue? \n", 
	       k+pbrow, INC_LVL_TEMP(k+pbrow), Options.inc_lvl, kid); 
	#endif

	
	//printf("considering K: %d kid: %d\n",
	//     k+L.scol, kid);

	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);

            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    #endif

	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("L.inc_lvl[%d]: %d %d kid: %d \n",
		   j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow),kid);
	    #endif

	    
	    //if((stack[jj]) > Options.inc_lvl)
	    if(color[jj] != 1)
	      {
		#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		printf("continue, already use Linc(%d): %d %d\n", 
		       j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow));
		#endif
		continue;
	      }
	    
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("jj: %d j: %d \n",
		   jj, j);
	    printf("VALUE: before: %f %f %f AFTER: %f kid: %d jj: %d k: %d\n",
		   X(jj), L.val(j), xj, X(jj)-L.val(j)*xj, 
		   kid, jj+L.srow, k+L.scol);
	    #endif

	    X(jj) -= L.val(j)*xj;
	    
	  }
	//INC_LVL_TEMP(k+pbrow) = BASKER_MAX_IDX;
      }//over all nonzero in left


    #ifdef BASKER_2DL
    #ifdef BAKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    LL(X_col)(X_row).p_size = nnz;
    #endif

     return;
  }//end t_offdiag_back_solve_same_pattern_inc_lvl();


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int 
  Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag_inc_lvl
  (Int kid, Int pbrow,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    BASKER_MATRIX &B     = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int         ws_size  = LL(X_col)(X_row).iws_size;
    
    Int          nnz     = LL(X_col)(X_row).p_size;
    //Int          brow    = L.srow;
    //Int          bcol    = L.scol;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);
   
    if(A_option == BASKER_TRUE)
      {
	#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	printf("k: %d size: %d \n",
	       k , B.col_ptr(k+1)-B.col_ptr(k));
	#endif

	for(Int i = B.col_ptr(k); 
	    i < B.col_ptr(k+1); i++)
	  {

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	    printf("t_back_solve_diag, kid: %d inc[%d] = 0\n",
		   kid, B.row_idx(i));
	    #endif
	    const Int j = B.row_idx(i);
	    color[j] = 1;
	    X(j) = B.val(i);
	    pattern[nnz++] = j;
	  }//over all nnz in subview
      }//end if preload

    //printf("done with reload\n\n");
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %d \n",
	   kid, x_size);
    #endif
    for(Int i = 0 ; i < x_size; ++i)
      {
	
	const Int k    = x_idx[i+x_offset];
	const Entry xj = x(i+x_offset);
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	       kid, k+pbrow, L.col_ptr[k], L.col_ptr[k+1]);
	#endif

        #ifdef BASKER_DEBUG_NFACTOR_INC_LVL
	//A multiple of a value at lvl l results in l+1 fill-in
	printf("LVL_TEMP[%d] = %d, %d kid: %d continue? \n", 
	       k+pbrow, INC_LVL_TEMP(k+pbrow), Options.inc_lvl, kid); 
	#endif

	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);

            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    #endif

	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("L.inc_lvl[%d]: %d %d kid: %d \n",
		   j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow),kid);
	    #endif

	    if((stack[jj]) > Options.inc_lvl)
	      {
		#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		printf("continue, already use Linc(%d): %d %d\n", 
		       j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow));
		#endif
		continue;
	      }
	    
	    if(color[jj] != 1)
	      {
		color[jj] = 1;
		#ifdef BASKER_DEBUG_NFACTOR_BLK
		printf("pattern index: %d kid: %d \n",
		       nnz, kid);
		#endif
		pattern[nnz++] = jj;
	      }
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("jj: %d j: %d \n",
		   jj, j);
	    printf("VALUE: before: %f %f %f AFTER: %f kid: %d jj: %d k: %d\n",
		   X(jj), L.val(j), xj, X(jj)-L.val(j)*xj, 
		   kid, jj+L.srow, k+L.scol);
	    #endif

	    X(jj) -= L.val(j)*xj;
	    
	  }
	//INC_LVL_TEMP(k+pbrow) = BASKER_MAX_IDX;
      }//over all nonzero in left


    #ifdef BASKER_2DL
    #ifdef BAKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    LL[X_col][X_row].p_size = nnz;
    #endif

     return 0;
  }//end t_offdiag_back_solve_inc_lvl();




  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int 
  Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag_inc_lvl_old
  (Int kid, Int pbrow,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    BASKER_MATRIX &B     = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int         ws_size  = LL(X_col)(X_row).iws_size;
    
    Int          nnz     = LL(X_col)(X_row).p_size;
    Int          brow    = L.srow;
    Int          bcol    = L.scol;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);
   
    //need to make this so not every column 
    for(Int i = 0 ; i < L.ncol; i++)
      {
	stack[i] = BASKER_MAX_IDX;
      }
    //printf("done filling\n");
    
    //Preload with A
    Int preload_fill = 0;
    if(A_option == BASKER_TRUE)
      {
	#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	printf("k: %d size: %d \n",
	       k , B.col_ptr(k+1)-B.col_ptr(k));
	#endif

	for(Int i = B.col_ptr(k); 
	    i < B.col_ptr(k+1); i++)
	  {

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	    printf("t_back_solve_diag, kid: %d inc[%d] = 0\n",
		   kid, B.row_idx(i));
	    #endif
	    const Int j = B.row_idx(i);
	    color[j] = 1;
	    X(j) = B.val(i);
	    pattern[nnz++] = j;
	    preload_fill++;
	    stack[j] = 0;

	  }//over all nnz in subview
      }//end if preload

    //printf("done with reload\n\n");
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %d \n",
	   kid, x_size);
    #endif
    for(Int i = 0 ; i < x_size; ++i)
      {
	
	const Int k    = x_idx[i+x_offset];
	const Entry xj = x(i+x_offset);
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	       kid, k+pbrow, L.col_ptr[k], L.col_ptr[k+1]);
	#endif

        #ifdef BASKER_DEBUG_NFACTOR_INC_LVL
	//A multiple of a value at lvl l results in l+1 fill-in
	printf("LVL_TEMP[%d] = %d, %d kid: %d continue? \n", 
	       k+pbrow, INC_LVL_TEMP(k+pbrow), Options.inc_lvl, kid); 
	#endif

	//We can't do this because might be one added
	//We could do this if we put a counter in the reload
	//Come back for this optimization
	/*
	if(preload_fill == 0)
	  {
	if(INC_LVL_TEMP(k+pbrow)+1 > Options.inc_lvl)
	  {
	    #ifdef BASKER_DEBUG_NFACTOR_INC_LVL
	    printf("spmv continue column: %d fill-lvl: %d \n",
		   k+brow, INC_LVL_TEMP(k+pbrow));
	    #endif
	    continue;
	  }
	  }
	*/
       

	//printf("==before col== size: %d\n",
	//     L.col_ptr(k+1)-L.col_ptr(k));
	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);

            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    #endif

	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("L.inc_lvl[%d]: %d %d kid: %d \n",
		   j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow),kid);
	    #endif

	    Int temp_cal = L.inc_lvl(j)+INC_LVL_TEMP(k+pbrow)+1; 
	    if(stack[jj] == BASKER_MAX_IDX)
	      {
		stack[jj] = temp_cal;
	      }
	    else
	      {
		stack[jj] = min(temp_cal, stack[jj]);
	      }


	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("Assigned inc_lvl(%d) = %d \n",
		   jj, stack[jj]);
	    #endif
	   

	    if((stack[jj]) > Options.inc_lvl)
	      {
		#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		printf("continue, already use Linc(%d): %d %d\n", 
		       j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow));
		#endif
		continue;
	      }
	    
	    if(color[jj] != 1)
	      {
		color[jj] = 1;
		#ifdef BASKER_DEBUG_NFACTOR_BLK
		printf("pattern index: %d kid: %d \n",
		       nnz, kid);
		#endif
		pattern[nnz++] = jj;
	      }
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("jj: %d j: %d \n",
		   jj, j);
	    printf("VALUE: before: %f %f %f AFTER: %f\n",
		   X(jj), L.val(j), xj, X(jj)-L.val(j)*xj );
	    #endif

	    X(jj) -= L.val(j)*xj;
	    
	  }
	//INC_LVL_TEMP(k+pbrow) = BASKER_MAX_IDX;
      }//over all nonzero in left


    #ifdef BASKER_2DL
    #ifdef BAKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    LL[X_col][X_row].p_size = nnz;
    #endif

     return 0;
  }//end t_offdiag_back_solve_inc_lvl();



  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_move_offdiag_L_inc_lvl
  (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot
   )
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
    const Int   ws_size = LL(X_col)(X_row).iws_size;
    const Int   p_size  = LL(X_col)(X_row).p_size;
   

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    //if(kid == 8 )
      {
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
      }
     #endif

   
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);

    //const Int    brow  = L.srow;
    //const Int    bcol  = L.scol;
    const Int    llnnz = L.nnz;
          Int    lnnz  = L.col_ptr(k);

    if((p_size) > (llnnz-lnnz))
      {
	Int newsize = llnnz*1.2 + L.ncol;

	if(Options.verbose == BASKER_TRUE)
	  {
	printf("-Warning, Need to remalloc L: %ld %ld kid: %ld current size: %ld used_size: %ld  addition: %ld newsize: %ld  \n",
	       (long)blkcol, (long)blkrow, (long)kid,
	       (long)llnnz, (long)lnnz, (long)p_size, (long)newsize  );
	  }
	BASKER_ASSERT(0==1, "REALLOC LOWER BLOCK\n");

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


       
      }

    for(Int i = 0; i < p_size; i++)
      {
	Int j = pattern[i];
	//Int t = gperm(j+brow);
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, kid: %d j: %d val: %f lnnz: %d inc: %d \n",
	       kid, j+L.srow, X[j]/pivot, lnnz, stack[j]);
	#endif

	

	color[j] = 0;
	L.row_idx(lnnz) = j;
	if(Options.same_pattern == BASKER_FALSE)
	  {
	    L.inc_lvl(lnnz) = stack[j];
	    stack[j] = BASKER_MAX_IDX;
	  }

	L.val(lnnz) = EntryOP::divide(X(j),pivot);
	X(j) = 0;

	lnnz++;
      }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //	   kid, k-bcol+1, lnnz);
    
    //Fix later 
    if(Options.same_pattern == BASKER_FALSE)
      {
    for(Int i = 0; i < LL(X_col)(X_row).nrow; i++)
      {
	stack[i] = BASKER_MAX_IDX;
      }
      }

    L.col_ptr(k+1) = lnnz;
    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_offdiag_mov_L_inc_lvl()



  //--------------------TEMP LOCATION-----------------//

  //Uses global indexing for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_nfactor_blk_old(Int kid)
  {
    Int b = S[0][kid]; //Which blk from schedule
    BASKER_MATRIX &L   = LL[b][0];
    BASKER_MATRIX &U   = LU[b][LU_size[b]-1];
    #ifdef BASKER_2DL
    printf("Accessing blk: %d \n", b);
    INT_1DARRAY   ws   = LL[b][0].iws;
    ENTRY_1DARRAY X    = LL[b][0].ews;
    Int        ws_size = LL[b][0].iws_size;
    #else  //else if BASKER_2DL
    INT_1DARRAY   ws   = thread_array[kid].iws;
    ENTRY_1DARRAY X    = thread_array[kid].ews;
    Int       ws_size  = thread_array[kid].iws_size;
    #endif
   
    Int          bcol  = L.scol;  //begining col
    Int          brow  = L.srow;  //begining row //Note: move out in future
    Int          lval  = 0;
    Int          uval  = 0;

    Int i,j,k;
    Int top, top1, maxindex, t; 
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;
   
    Int newsize;
    Entry pivot, value;
    Entry absv, maxv;

    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    Int scol  = L.scol; //Note: this seems like over kill --clean up variables
    Int ecol  = L.ecol;
    
    Int col_idx_offset = A.nnz;

    Int *color    = &(ws[0]);
    Int *pattern  = &(color[ws_size]);
    
    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
    top1 = ws_size;
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("b: %d scol: %d ecol: %d llnzz: %d uunzz: %d \n", 
           b, scol, ecol, L.nnz, U.nnz);
    #endif
    
    //return 0 ;

    //---TEMP--- DEBUG ---
    //ecol = 5;

    //for each column
    for(k = scol; k < ecol; k++)
        {
          #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
          printf("\n----------------K=%d--------------\n", k);
	  #endif
	  value = 0.0;
	  pivot = 0.0;
	  maxindex = A.ncol;  
	  lcnt = 0;
	  ucnt = 0;

          #ifdef BASKER_DEBUG_NFACTOR_BLK
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
	      
	      #ifdef BASKER_2DL
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
	      X[j-brow] = A.val[i];
	      #else
              X[j] = A.val[i];
	      #endif
      
              #ifdef BASKER_DEBUG_NFACTOR_BLK
              printf("i: %d row: %d  val: %g  top: %d \n", 
                     i, j ,A.val[i], top);
	      #ifdef BASKER_2DL
	      printf("Nx in Ak %d %g %d color = %d \n",
                      j, X[j-brow], brow,  
                     color[j-brow] );
	      #else
              printf("Nx in Ak %d %g %d color = %d \n",
                      j, X[j], brow,  
                     color[j] );
	      #endif
              #endif

              //Search reach if not yet considered
	      #ifdef BASKER_2DL
	      if(color[j-brow] == 0)
	      #else
	      if(color[j] == 0)
	      #endif
		{        

		  #ifdef BASKER_INC_LVL
		  t_local_reach_selective(kid,0,0,j, &top);
                  #else
		  //t_local_reach(kid, 0,0, j, &top);
		  #endif
		}//if not colored	      
          }//end for() each nnz in column
	  xnnz = ws_size - top;

          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("xnnz: %d ws_size: %d top: %d \n", 
                 xnnz, ws_size, top);
          #endif
             

	  #ifdef BASKER_INC_LVL
	  t_back_solve_selective(kid, 0, 0, k, top, xnnz);
	  #else
          t_back_solve(kid, 0,0,  k, top, xnnz);
	  #endif

	  //Future add
	  //t_locate_pivot(kid, top)	  
          //find pivot
          maxv = 0.0;
          for(i = top; i < ws_size; i++)
            {
	      j = pattern[i];
              t = gperm[j];
	      #ifdef BASKER_2DL
	      value = X[j-brow];
	      #else
              value = X[j];
	      #endif

              absv = abs(value);
              //if(t == L.max_idx)
	      if(t == BASKER_MAX_IDX)
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

	  if(Options.no_pivot == BASKER_TRUE)
	    {
	      maxindex = k;
	      pivot = X[k-brow];
	    }
      
          ucnt = ws_size - top - lcnt +1;
          //if((maxindex == L.max_idx) || (pivot == 0))
	  if((maxindex == BASKER_MAX_IDX) || (pivot == 0))
            {

	      if(Options.verbose == BASKER_TRUE)
		{
	      cout << endl << endl;
	      cout << "---------------------------"
		   <<endl;
	     
              cout << "Error: Matrix is singular, blk" 
		   << endl;
              cout << "MaxIndex: " 
		   << maxindex << " pivot " 
                   << pivot << endl;
		}
              return 2;
            }          

          gperm[maxindex] = k;
	  gpermi[k] = maxindex;
          //printf("TAG1 r  maxindex = %d k= %d \n", maxindex, k);
          #ifdef BASKER_DEBUG_NFACTOR
          if(maxindex != k)
            {
              cout << "Permuting Pivot: " << k << " as row " 
                   << maxindex << endl;
            }
          #endif
          
          //Note: Come back to this!!!!
          if(lnnz + lcnt > llnnz)
            {
	      printf("\n\n");
	      printf("----------------------\n");

              newsize = lnnz * 1.1 + 2 *A.nrow + 1;
              printf("b: %ld Reallocing L oldsize: %ld current: %ld count: %ld newsize: %ld \n",
                     (long)b, (long)llnnz, (long)lnnz, (long)lcnt, (long)newsize);
            }
          if(unnz+ucnt > uunnz)
            {

	      printf("\n\n");
	      printf("-------------------\n");

              newsize = uunnz*1.1 + 2*A.nrow+1;
              printf("b: %ld Reallocing U oldsize: %ld newsize: %ld \n",
                     (long)b, (long)uunnz, (long)newsize);
            }

          L.row_idx[lnnz] = maxindex;
          L.val[lnnz] = (Entry) 1.0;
          lnnz++;
     
          Entry lastU = (Entry) 0.0;
          for( i = top; i < ws_size; i++)
            {
	      j = pattern[i];
              t = gperm[j];
            
              #ifdef BASKER_DEBUG_NFACTOR_BLK
              printf("j: %d t: %d \n", j, t);
              #endif            

	      //temp
	      /*
	      if(k<6)
		{
		  printf("k: %d j: %d t: %d  val: %f count: %d \n", k, j, t, X[j-brow], unnz);  
		  cout << "t: " << t << endl;
		}
	      */

              //if fill-in
	      #ifdef BASKER_2DL
	      if(X[j-brow] != 0)
	      #else
              if(X[j] != 0)
	      #endif
                {
                  //if(t != L.max_idx)
		  if(t != BASKER_MAX_IDX)
                    {
                      if(t < k)
                        {
                          U.row_idx[unnz] = gperm[j];
			  #ifdef BASKER_2DL
			  U.val[unnz] = X[j-brow];
			  #else
                          U.val[unnz] = X[j];
			  #endif
                          unnz++;
                        }
                      else
                        {
			  #ifdef BASKER_2DL
			  lastU = X[j-brow];
			  #else
                          lastU = X[j];
			  #endif
                        }
                    }
                  //else if (t == L.max_idx)
		  else if (t == BASKER_MAX_IDX)
                    {

		      /*
		      if(k < 6)
			{
			  printf("Idx: %d x: %f pivot: %f result: %f \n", j, X[j-brow], pivot, X[j-brow]/pivot);
			  
			}
		      */

                      L.row_idx[lnnz] = j;
		      #ifdef BASKER_2DL
		      L.val[lnnz] = X[j-brow]/pivot;
		      #else
                      L.val[lnnz] = X[j]/pivot;
		      #endif

		      #ifdef BASKER_INC_LVL
		      L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
		      #endif

                      lnnz++;

		      /*
		      if(k < 6)
			{
		      printf("L, placing.  Idx: %d nnz: %d val: %f \n", L.row_idx[lnnz-1], lnnz-1, L.val[lnnz-1]);
			}
		      */
                    }
                }//end if() not 0
              
              //Note: move x[j] inside of if() not 0....extra ops this way
              #ifdef BASKER_DEBUG_NFACTOR_BLK
              printf("Zeroing element: %d \n", j);
              #endif

	      #ifdef BASKER_2DL
	      X[j-brow] = 0;
	      #else
              X[j] = 0;
	      #endif
            }//end if(x[i] != 0)

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


	  #ifdef BASKER_2DL
	  //-----------------------Update offdiag-------------//
	  for(Int blk_row = 1; blk_row < LL_size[b]; blk_row++)
	    {
	      //Do back solve of off-diag blocks
	      #ifdef BASKER_INC_LVL
	      t_back_solve_offdiag_selective(kid,
				   b, blk_row,
				   b, blk_row,
				   k, col_idx_offset,
				   U.val, U.row_idx,
	       U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
				   U.col_ptr[k-bcol],
				   BASKER_TRUE);
	      #else
	      t_back_solve_offdiag(kid,
				   b, blk_row,
				   b, blk_row,
				   k, col_idx_offset,
				   U.val, U.row_idx,
	       U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
				   U.col_ptr[k-bcol],
				   BASKER_TRUE);
	      #endif
	      //Move these factors into Local Ls
	      t_move_offdiag_L(kid,
			       b, blk_row,
			       b, blk_row,
			       k, pivot);
	    }//end over all diag
	  #endif

	  col_idx_offset = A.nnz;
	}//end for() over all columns

    L.nnz = lnnz;
    U.nnz = unnz;

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    //print_factor(L,U);
    #endif
    return 0;
  }//end t_nfactor_blk()


  
  //Uses global idx for local blks
  //Note: need to fix indexing for color
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach_old_old
  (Int kid, Int lvl, Int l,
   Int j, Int *top)
  {

    //Setup variables
    const Int      b   = S[lvl][kid];
    const Int     wsb  = S[0][kid];
    BASKER_MATRIX  &L  = LL[b][0];
    #ifdef BASKER_2DL
    INT_1DARRAY    ws  = LL[wsb][l].iws;
    Int        ws_size = LL[wsb][l].iws_size;
    #else
    INT_1DARRAY    ws  = thread_array[kid].iws;
    Int        ws_size = thread_array[kid].iws_size;
    #endif

    const Int brow    = L.srow;

    Int *color       = &(ws[0]);
    Int *pattern     = &(ws[ws_size]);
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);

    Int i, t, head, i1;
    Int start, end, done;
    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d, kid: %d \n",
	   b, 0, wsb, l, kid);
    #endif

    start = -1;

    head = 0;
    //ws[stack_offset + head] = j;
    stack[head] = j;

    //while(head < L.max_idx)
    while(head != BASKER_MAX_IDX)
      { 
        #ifdef BASKER_DEBUG_LOCAL_REACH
        printf("stack_offset: %d head: %d \n", stack_offset , head);
        ASSERT(head > -1);
        #endif

	j = stack[head];
        t = gperm[j];
        
        #ifdef BASKER_DEBUG_LOCAL_REACH
	printf("----------DFS: %d %d -------------\n", j, t);
        #endif

	#ifdef BASKER_2DL
	if(color[j-brow] == 0)
	#else
        if(color[j] == 0)
	#endif
	  {	    
	    #ifdef BASKER_2DL
	    color[j-brow] = 1;
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
	    start = store[j-brow];
	  }
	done = 1;
	
        //if  permuted
        //if(t != L.max_idx), might be able to make this smaller now
        //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
	//if(t!=L.max_idx)
	if(t!=BASKER_MAX_IDX)
	  {
	   
	    end = L.col_ptr[t+1-L.scol];
	    for(i1 = start; i1 < end; i1++)
	      {
                i = L.row_idx[i1];
		#ifdef BASKER_2DL
		if(color[i-brow] == 0)
		#else
		if(color[i] == 0)
		#endif
		  {
                    head++;
		    stack[head] = i;
		    store[j-brow] = i1+1;
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
	    #ifdef BASKER_2DL
	    color[j-brow] = 2;
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


  //Note that this can be combined with t_back_solve_offdiage_same
  template <class Int, class Entry, class Exe_Space>
  void 
  Basker<Int,Entry,Exe_Space>::t_same_pattern_back_solve_offdiag_inc_lvl
    (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int UP_col, Int UP_row,
   Int LP_col, Int LP_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   INT_1DARRAY   x_fill, 
   Int x_size, Int x_offset,
   BASKER_BOOL A_option
     )
  {
    BASKER_MATRIX &L            = LL(blkcol)(blkrow);
    BASKER_MATRIX &B            = ALM(blkcol)(blkrow);


    /*
    printf("TEST, UP PATTERN %d %d kid: %d \n",
	   UP_col, UP_row, kid);
    printf("TEST, LP PATTERN %d %d kid: %d \n",
	   LP_col, LP_row, kid);
    */

    BASKER_MATRIX *UPP = &LU(UP_col)(0);
    if(UP_row != BASKER_MAX_IDX)
      {
	UPP = &(LU(UP_col)(UP_row));
      }
    BASKER_MATRIX &UP = *(UPP);

    BASKER_MATRIX *LPP = &LU(LP_col)(0);
    if(LP_row != BASKER_MAX_IDX)
      {
	LPP = &(LL(LP_col)(LP_row));
      }
    BASKER_MATRIX &LP = *(LPP);
    


    INT_1DARRAY   ws            = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X             = LL(X_col)(X_row).ews;
    Int         ws_size         = LL(X_col)(X_row).iws_size;
    
    Int    nnz            = LL(X_col)(X_row).p_size;
 

   

   
    Int *color =   &(ws(0));
    Int *pattern = &(color[ws_size]);
//    Int *stack   = &(pattern[ws_size]); //Temp store the fill-in //NDE - warning: unused
    
    //For Debug
    /*
    auto color = Kokkos::subview(ws, std::make_pair((Int)0,ws_size));
    auto pattern = Kokkos::subview(ws, std::make_pair(ws_size,2*ws_size));
    auto stack = Kokkos::subview(ws, std::make_pair(2*ws_size, 3*ws_size));
    */

    #ifdef BASKER_DEBUG_NFACTOR_BlK_INC
    printf("\n\n===============DEBUG FILL AT SPMV========\n\n");
    for(Int i = 0 ; i < L.nrow; i++)
      {
	printf("k: %d i: %d fill: %d  kid: %d\n",
	       k ,i+L.srow, stack[i], kid);
      }
    printf("\n");
    #endif
    

    //printf("-----------Before loading pattern----- kid: %d \n",
    //	   kid);
    //Preload pattern of U into color
    if(UP_row != BASKER_MAX_IDX)
      {
	for(Int i = UP.col_ptr(k); i < UP.col_ptr(k+1); i++)
	  {
	    const Int j = UP.row_idx(i);
	    if(color[j] != 2)
	      {
		/*
	    printf("Added Upattern: %d %d nnz: %d colors: %d psize: %d kid: %d \n",
		   j, j+UP.srow,nnz,
		   color.extent(0),
		   pattern.extent(0),
		   kid);
		*/
	    color[j] = 2;
	    pattern[nnz++] = j;
	      }
	  }
      }
    //Preload pattern of L
    if(LP_row != BASKER_MAX_IDX)
      {
	for(Int i = LP.col_ptr(k); i < LP.col_ptr(k+1); i++)
	  {
	    const Int j = LP.row_idx(i);
	    /*
	    printf("Added Lpattern: %d %d nnz: %d kid: %d \n",
		   j, j+LP.srow, nnz, kid);
	    */

	    if(color[j] != 2)
	      {
		color[j] = 2;
		pattern[nnz++] = j;
	      }

	  }
      }

    //   printf("===========Done loading pattern kid: %d \n",
    //	   kid);

    /*
    //DEBUG test x is filled
    for(Int i = 0; i < L.nrow; i++)
      {
	if( X(i) != 0 )
	  {
	    printf("\n\n   ERROR X(%d) %f kid: %d \n\n",
		   i+L.srow, X(i), kid);

	  }
      }
    */



    //Preload values of with A
    if(A_option == BASKER_TRUE)
      {
        #ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	#endif

	for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
	  {
	    const Int j = B.row_idx(i);

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i)+B.srow, B.val(i));
	    printf("t_back_solve_diag, kid: %d x: %f %f \n",
		   kid, X(j), X(j)+B.val(i));
	    #endif
	    	   
	    // printf("Add X[%d %d] %f %f kid: %d \n",
	    //	   k+L.scol, j, X(j), B.val(i), kid);
	    X(j) = X(j)+ B.val(i);
	  }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if((kid == 2)||(kid==3))
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
	   kid, x_size);
    #endif

    
    for(Int i = 0 ; i < x_size; ++i)
      {
	const Int kk    = x_idx(i+x_offset);
	const Entry xj = x(i+x_offset);
	if((A_option == BASKER_TRUE) &&
	   (kk == k))
	  {
	    //printf("CONT called %d kid: %d \n", 
	    //	   kk, kid);
	    continue;
	  }
	

        #ifdef BASKER_DEBUG_NFACTOR_BLK
	if((kid == 2)||(kid==3))
	printf("t_back_solve_diag, kid: %d  k: %d %g  x_size: %d [%d %d] \n",
	       kid, k, xj, x_size, L.col_ptr[kk], L.col_ptr[kk+1]);
	#endif
       
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	if((kid == 2)||(kid==3))
	printf("L_size: %d k: %d kid: %d \n",
	       L.col_ptr(kk+1)-L.col_ptr(kk), k, kid);
	#endif

	for(Int j = L.col_ptr(kk); 
	    j < L.col_ptr(kk+1); j++)
	  {
	    const Int jj = L.row_idx(j);
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //if(kid ==0)
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj+L.srow, color[jj]);
	    #endif


	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    if((kid == 2)||(kid==3))
	    printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
		    kid, jj+L.srow,X[jj], L.val[j], xj);
	    #endif 
	     
	 
	  
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	      //if((kid == 2)||(kid==3))
	    {
	      printf("VALUE: before %f %f %f AFTER: %f kid: %d jj: %d K: %d \n",
		     X(jj), L.val(j), xj, X(jj)-L.val(j)*xj, 
		     kid, jj+L.srow, kk+L.scol);
	    }
            #endif
	       
	    if(color[jj]  == 2)
	      {  


		X(jj) -= L.val(j)*xj;
	      }
	  }//over all nnz in column
	//REset moved from lower to here  ... move to caller
	//INC_LVL_TEMP(k+LL(blkcol)(0).srow) = BASKER_MAX_IDX;
      }//over all nonzero in left


    LL(X_col)(X_row).p_size = nnz;

    return;

  }//end t_same_pattern_back_solve_offdiage_inc_lvl

  
  template <class Int, class Entry, class Exe_Space>
  int 
  Basker<Int,Entry,Exe_Space>::t_dense_back_solve_offdiag_inc_lvl
  (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   INT_1DARRAY   x_fill, 
   Int x_size, Int x_offset,
   BASKER_BOOL A_option)
  {
    BASKER_MATRIX &L            = LL(blkcol)(blkrow);
    BASKER_MATRIX &B            = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws            = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X             = LL(X_col)(X_row).ews;
    Int         ws_size         = LL(X_col)(X_row).iws_size;
    
    Int    nnz            = LL(X_col)(X_row).p_size;
    //const Int    brow           = L.srow;
    //const Int    bcol           = L.scol;
  

    //printf("=============RIGHT ONE CALLED ============\n");


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

    Int *color =   &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]); //Temp store the fill-in


    #ifdef BASKER_DEBUG_NFACTOR_BlK_INC
    printf("\n\n===============DEBUG FILL AT SPMV========\n\n");
    for(Int i = 0 ; i < L.nrow; i++)
      {
	printf("k: %d i: %d fill: %d  kid: %d\n",
	       k ,i+L.srow, stack[i], kid);
      }
    printf("\n");
    #endif
    
    //Preload with A
    if(A_option == BASKER_TRUE)
      {
        #ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	#endif

	for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
	  {
	    const Int j = B.row_idx(i);

	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i)+B.srow, B.val(i));
	    printf("t_back_solve_diag, kid: %d x: %f %f \n",
		   kid, X(j), X(j)+B.val(i));
	    #endif
	    	   
	    X(j) = X(j)+ B.val(i);
	  }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if((kid == 2)||(kid==3))
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
	   kid, x_size);
    #endif

    for(Int i = 0 ; i < x_size; ++i)
      {
	const Int k    = x_idx(i+x_offset);
	const Entry xj = x(i+x_offset);

        #ifdef BASKER_DEBUG_NFACTOR_BLK
	if((kid == 2)||(kid==3))
	printf("t_back_solve_diag, kid: %d  k: %d %g  x_size: %d [%d %d] \n",
	       kid, k, xj, x_size, L.col_ptr[k], L.col_ptr[k+1]);
	#endif
       
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	if((kid == 2)||(kid==3))
	printf("L_size: %d k: %d kid: %d \n",
	       L.col_ptr(k+1)-L.col_ptr(k), k, kid);
	#endif

	for(Int j = L.col_ptr(k); 
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //if(kid ==0)
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj+L.srow, color[jj]);
	    #endif


	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    if((kid == 2)||(kid==3))
	    printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
		    kid, jj+L.srow,X[jj], L.val[j], xj);
	    #endif 
	     
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    //if((kid == 2)||(kid==3))
	     printf("L.inc_lvl[%d %d]: %d %d kid: %d \n",
		    jj+L.srow, j, L.inc_lvl(j), stack[jj], kid);
	    #endif

	     #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    if((kid == 2)||(kid==3))
	     printf("Assigned inc_lvl(%d) = %d \n",
		    jj, stack[jj]);
	     #endif





	    //TEST--- CAN BE REMOVED LATER
	    if(A_option == BASKER_TRUE)
	      {
		/* 
	    if(INC_LVL_TEMP(k+LL(blkcol)(0).srow) == BASKER_MAX_IDX)
	      {
		printf("ERROR INC_TEMP RESET LOWER BACKSOLVE\n");

	      }
		*/
		

	    Int temp = INC_LVL_TEMP(k+LL(blkcol)(0).srow)+ L.inc_lvl(j) + 1;
	   
	    /*
	    printf("lower row: %d kid: %d inc: %d %d %d j: %d \n",
		   jj+L.srow,
		   kid, stack[jj],
		   INC_LVL_TEMP(k+LL(blkcol)(0).srow),
		   temp,
		   k+LL(blkcol)(0).srow);
	    */	   
	    
				


	    //stack[jj] = min(stack[jj], 
	    //		    INC_LVL_TEMP(k+LL(blkcol)(0).srow)+ L.inc_lvl(j) + 1);
	       
	    if(stack[jj] != BASKER_MAX_IDX)
	      {
	    stack[jj] = min(stack[jj], temp);
	      }
	    else
	      {
		stack[jj] = temp;
	      }

	    /*
	     printf("lower update  row: %d kid: %d inc: %d \n", 
	    	   jj+L.srow,
	    	   kid, 
	    
		   stack[jj]);
	    */

	      }

	    /*
	    if(Options.incomplete_type == 
	       BASKER_INCOMPLETE_LVL)
	      {
	     
	     if((stack[jj] == BASKER_MAX_IDX) ||
		(stack[jj] > Options.inc_lvl))
	       {
                 //printf("continued compute not done \n");

		 continue;
	       }
	      }
	    */
	    

	  
	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	      //if((kid == 2)||(kid==3))
	       {
	     printf("VALUE: before %f %f %f AFTER: %f kid: %d jj: %d K: %d \n",
		    X(jj), L.val(j), xj, X(jj)-L.val(j)*xj, 
		    kid, jj+L.srow, k+L.scol);
	       }
	      #endif
	     
	     X(jj) -= L.val(j)*xj;

	  }//over all nnz in column
	//REset moved from lower to here  ... move to caller
	//INC_LVL_TEMP(k+LL(blkcol)(0).srow) = BASKER_MAX_IDX;
      }//over all nonzero in left

    #ifdef BASKER_2DL
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
	   kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
	   nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    LL(X_col)(X_row).p_size = nnz;
    #endif

    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if((kid==2)||(kid==3))
      {
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
      }
    #endif

    return 0;
  }//end t_dense_back_solve_offdiag_inc_lvl();


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_dense_move_offdiag_L_inc_lvl
  (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot
   )
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
    const Int   ws_size = LL(X_col)(X_row).iws_size;
    //const Int   p_size  = LL(X_col)(X_row).p_size; //NDE - warning: unused


   
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);

    //const Int    brow  = L.srow;
    //const Int    bcol  = L.scol;
    const Int    llnnz = L.nnz;
    Int    lnnz  = L.col_ptr(k);
   
    if(L.nrow > (llnnz - lnnz))
      {
	printf("no enough memory in dense \n");
//	printf("kid: %ld llnnz: %ld lnnz: %ld \n",kid, llnnz, lnnz);
  std::cout << "kid: " << kid
            << " llnz: " << llnnz
            << " lnz: " << lnnz << std::endl;

      }
    

    for(Int j = 0; j < L.nrow; ++j)
      {

	//printf("Consider L-move (%d %d) %f kid: %d inc: %d \n",
	//   k+L.scol, j+L.srow, X(j), kid, stack[j]);

	//if((X(j)!=0))
	if((X(j)!=(Entry)(0)))
	  {

	    if((stack[j] != BASKER_MAX_IDX)&&
	       stack[j] <= Options.inc_lvl)
	      {
	    
	   #ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, [%d %d] kid: %d j: %d val: %f lnnz: %d inc: %d L: %d %d ws: %d %d\n",
	       k+L.scol, j+L.srow,
	       kid, j, X[j]/pivot, lnnz, stack[j], 
	       blkcol, blkrow,
	       X_col, X_row);
	#endif	

	//printf("lnnz: %d %d kid: %d \n",
	//     lnnz, j, kid); 
	L.row_idx(lnnz) = j;
	L.val(lnnz) = EntryOP::divide(X(j),pivot);
        L.inc_lvl(lnnz) = stack[j];
	lnnz++;
	      }
	
	    //X(j) = 0;
	    //color[j] = 0;
	
	//Clear stack too ?
	  }//if not zero
	color[j] = 0;
	X(j)  = 0;
	//printf("clear stack %d kid: %d \n", 
	//     j+L.srow, kid);
	stack[j] = BASKER_MAX_IDX;
      }

    L.col_ptr(k+1) = lnnz;
    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_dense_offdiag_mov_L_inv_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::t_dom_lower_col_offdiag_find_fill
  (
   const Int kid, const Int pbrow,
   const Int blkcol, const Int blkrow,
   const Int X_col, const Int X_row,
   const Int k,
   INT_1DARRAY x_idx,
   const Int x_size, 
   const Int  x_offset,
   const BASKER_BOOL A_option
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    BASKER_MATRIX &B     = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int         ws_size  = LL(X_col)(X_row).iws_size;
    
    //Int          nnz     = LL(X_col)(X_row).p_size;
    //Int          brow    = L.srow;
    //Int          bcol    = L.scol;
  
    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);
   
    
    //need to make this so not every column 
    for(Int i = 0 ; i < ws_size; i++)
      {
	stack[i] = BASKER_MAX_IDX;
	//printf("clearing stack: %d \n",
	//     i+LL(X_col)(X_row).srow);
      }
    
  
    //Preload with A
    if(A_option == BASKER_TRUE)
      {
	#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	printf("k: %d size: %d \n",
	       k , B.col_ptr(k+1)-B.col_ptr(k));
	#endif

	for(Int i = B.col_ptr(k); 
	    i < B.col_ptr(k+1); i++)
	  {
	    const Int j = B.row_idx(i);
	    stack[j] = 0;
	  }//over all nnz in column
      }//end if preload

    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %d \n",
	   kid, x_size);
    #endif
    for(Int i = 0 ; i < x_size; ++i)
      {
	
	const Int k    = x_idx[i+x_offset];
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	       kid, k+pbrow, L.col_ptr[k], L.col_ptr[k+1]);
	#endif

        #ifdef BASKER_DEBUG_NFACTOR_INC_LVL
	printf("LVL_TEMP[%d] = %d, %d kid: %d continue? \n", 
	       k+pbrow, INC_LVL_TEMP(k+pbrow), Options.inc_lvl, kid); 
	#endif
	
	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);

            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    #endif

	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("L.inc_lvl[%d]: %d %d kid: %d \n",
		   j, L.inc_lvl(j), INC_LVL_TEMP(k+pbrow),kid);
	    #endif

	    Int temp_cal = L.inc_lvl(j)+INC_LVL_TEMP(k+pbrow)+1; 
	    if(stack[jj] == BASKER_MAX_IDX)
	      {
		stack[jj] = temp_cal;
	      }
	    else
	      {
		stack[jj] = min(temp_cal, stack[jj]);
	      }


	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("Assigned inc_lvl(%d) = %d \n",
		   jj, stack[jj]);
	    #endif
	    
          }//over the column
	//can not clear here
	//INC_LVL_TEMP(k+pbrow) = BASKER_MAX_IDX;
      }//over all nonzero in left

  }//end t_dom_lower_col_offdiag_find_fill
   

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_lower_col_offdiag_find_fill
  (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k,
   ENTRY_1DARRAY x,
   INT_1DARRAY x_idx,
   INT_1DARRAY x_fill,
   Int x_size, Int x_offset
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    
    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int        ws_size   = LL(X_col)(X_row).iws_size;

    //Int nnz              = LL(X_col)(X_row).p_size;
    //const Int brow       = L.srow; //Not used
    //const Int bcol       = L.scol; //Not used

    
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);

    // printf("find fill kid: %d L: %d %d X: %d %d xsize: %d \n",
    //	   kid, blkcol, blkrow, X_col, X_row, x_size);

    //Need to find away not to prefill everytime
    for(Int i = 0; i < L.nrow; i++)
      {
	stack[i] = BASKER_MAX_IDX;
	//printf("clearing stack2: %d \n",
	//     i+LL(X_col)(X_row).srow);
      }

    //Simulate SPMV to get fill-in pattern
    for(Int i = 0; i < x_size; i++)
      {
	const Int k     = x_idx(i+x_offset);
	const Int flvl  = x_fill(i+x_offset);

	/*
	if(kid == 1)
	  {
	    printf("consider U_fill: %d %d \n",
		   k+L.scol, flvl);
	  }
	*/

	//NDE Unnecessary...
	if(flvl+1 > Options.inc_lvl)
	  {
	    //printf("Continued skip because too large\n");
	    continue;
	  }
	
	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj    = L.row_idx(j);
	    const Int nflvl = L.inc_lvl(j) + flvl + 1;
	    

	    /*
	    if((jj+L.srow) == 40)
	      {
		
		printf("fill-32.  kid: %d col: %d xj-lvl: %d col-lvl: %d %d\n",
		       kid, k+L.scol, 
		       flvl,
		       L.inc_lvl(j),
		       nflvl);
		       
	      }
	    */
	    
	    
	    /*
	     if((jj+L.srow) == 20)
	       {
			
		printf("fill-20.  kid: %d col: %d xj-lvl: %d col-lvl: %d %d \n",
		       kid, k+L.scol, 
		       flvl,
		       L.inc_lvl(j),
		       nflvl);
		       
	      }
	    */




            // if((kid==2)||(kid==3))
	    //  {
	    //printf("%d -- flvl: %d nflvl: %d %d %d kid: %d \n",
	    //	   jj+L.srow, stack[jj], 
	    //	   nflvl, L.inc_lvl(j), flvl,
	    //	   kid);
            // }
	    if(stack[jj] == BASKER_MAX_IDX)
	      {
		stack[jj] = nflvl;
	      }
	    else
	      {
		stack[jj] = min(nflvl, stack[jj]);
	      }

            //if((kid==2)||(kid==3))
	    //  {
            //printf("Assigned inc_lvl(%d) = %d, kid: %d \n",
	    //	   jj+L.srow, stack[jj], kid);
            // }

	  }//end for-all nnz in L
      }//over all nnz

    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("Fill-in after test kid: %d blk:  %d %d \n",
           kid, blkcol, blkrow);
    for(Int i = 0; i < L.nrow; i++)
      {
        printf("kid: %d fi: %d %d \n",
               kid, i+L.srow, stack[i]);
      }
    #endif

    return 0;
  }//end t_lower_col_diag_find_fill()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_lower_col_offdiag_find_fill_rlvl
  (
   Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k,
   ENTRY_1DARRAY x,
   INT_1DARRAY x_idx,
   INT_1DARRAY x_fill,
   Int x_size, Int x_offset
   )
  {
    BASKER_MATRIX &L     = LL(blkcol)(blkrow);
    
    INT_1DARRAY   ws     = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X      = LL(X_col)(X_row).ews;
    Int        ws_size   = LL(X_col)(X_row).iws_size;

    //Int nnz              = LL(X_col)(X_row).p_size;
    //const Int brow       = L.srow; //Not used
    //const Int bcol       = L.scol; //Not used

    
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);

    /*
    printf("find fill rlvl kid: %d L: %d %d X: %d %d xsize: %d \n",
    	   kid, blkcol, blkrow, X_col, X_row, x_size);

    if(X_row == 3)
      {
	printf("TEST kid: %d stack[%d] = %d \n",
	       kid, 30, stack[0]);
      }
    */

    //Need to find away not to prefill everytime
    /*
    for(Int i = 0; i < L.nrow; i++)
      {
	//This will over write, need to fix
	//stack[i] = BASKER_MAX_IDX;
	//printf("clearing stack2: %d \n",
	//     i+LL(X_col)(X_row).srow);
	

	//Check if z
	
	if(stack[i] != BASKER_MAX_IDX)
	  {
	    printf("EERRO kid: %d stack not cleared: %d \n",
		   kid, i+L.srow);
	  }
	
      }
    */
    

    //Simulate SPMV to get fill-in pattern
    for(Int i = 0; i < x_size; i++)
      {
	const Int k     = x_idx(i+x_offset);
	const Int flvl  = x_fill(i+x_offset);

	if(kid == 1)
	  {
	    // printf("consider U_fill: %d %d \n",
	    //	   k+L.scol, flvl);
	  }

	//In RLVL we can not skip out
		
	for(Int j = L.col_ptr(k);
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj    = L.row_idx(j);
	    const Int nflvl = L.inc_lvl(j) + flvl + 1;
	    
	    /*
	    if((jj+L.srow) == 90)
	      {

		printf("fill-90.  kid: %d col: %d xj-lvl: %d col-lvl: %d %d\n",
		       kid, k+L.scol, 
		       flvl,
		       L.inc_lvl(j),
		       nflvl);
		       
	      }
	    */
	    
	    
	    /*
	     if((jj+L.srow) == 14)
	       {
			
		printf("fill-14.  kid: %d col: %d xj-lvl: %d col-lvl: %d %d \n",
		       kid, k+L.scol, 
		       flvl,
		       L.inc_lvl(j),
		       nflvl);
		       
	      }
	    */
	    




            // if((kid==2)||(kid==3))
	    //  {
	    //printf("%d -- flvl: %d nflvl: %d %d %d kid: %d \n",
	    //	   jj+L.srow, stack[jj], 
	    //	   nflvl, L.inc_lvl(j), flvl,
	    //	   kid);
            // }
	    if(stack[jj] == BASKER_MAX_IDX)
	      {
		stack[jj] = nflvl;
	      }
	    else
	      {
		stack[jj] = min(nflvl, stack[jj]);
	      }
	    
	    /*
	    if((jj+L.srow) == 14)
	      {

		printf("fill2-14.  kid: %d col: %d set: %d\n",
		       kid, k+L.scol, 
		       stack[jj]);
	   
		       
	      }
	    */
	    



            //if((kid==2)||(kid==3))
	    //  {
            //printf("Assigned inc_lvl(%d) = %d, kid: %d \n",
	    //	   jj+L.srow, stack[jj], kid);
            // }

	  }//end for-all nnz in L
      }//over all nnz

    //printf("TEST 30 kid: %d stack: %d \n",
    //	   kid, stack[0]);

    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("Fill-in after test kid: %d blk:  %d %d \n",
           kid, blkcol, blkrow);
    for(Int i = 0; i < L.nrow; i++)
      {
        printf("kid: %d fi: %d %d \n",
               kid, i+L.srow, stack[i]);
      }
    #endif

    return 0;
  }//end t_lower_col_diag_find_fill_rlvl()


 
  //Notes on t_popluate_col_fill
  //Assume pair is done is done finding local fill-in
  //This will be done in parallel down a column length
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_populate_col_fill
  (
   const Int kid,
   const Int blkcol, const Int blkrow,
   const Int X_col,  const Int X_row,
   const Int k, const BASKER_BOOL lower
   )
  {

    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    //if((kid == 2) || (kid == 3))
      {
    printf("============t_pop_col_fill called kid: %d============\n", kid);
    printf("populating k: %d B: %d %d X: %d %d kid: %d \n",
	   k, blkcol,blkrow,X_col,X_row, kid);
      }
     #endif
    //Note X_col is always the leader
    BASKER_MATRIX *B;
    if(lower == BASKER_TRUE)
      {
	B = &(ALM(blkcol)(blkrow));
      }
    else
      {
	B = &(AVM(blkcol)(blkrow));
      }
    BASKER_MATRIX &M = *B;
    //BASKER_MATRIX &M = ALM(blkcol)(blkrow);
    INT_1DARRAY  ws   = LL(X_col)(X_row).iws;
    const Int ws_size = LL(X_col)(X_row).iws_size;
    
    Int *color    = &(ws(0));
    Int *pattern  = &(color[ws_size]);
    Int *stack    = &(pattern[ws_size]);


    //printf("Crazy test: %d %d \n", 
    //	   X_col, X_row);
    // for(Int i = 0; i < M.nrow; i++)
    // {
    //	printf("scrazy: %d %d \n",
    //	       i+M.srow, stack[i]);
    // }

    //printf("\n");
    //Scan the column to see for nnz
    for(Int jj = M.col_ptr(k); jj < M.col_ptr(k+1); jj++)
      {


	const Int j = M.row_idx(jj);

        //printf("test pop: %d %d %d kid: %d \n",
        //     jj, j, M.srow, kid);
	//#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	//if((kid==2)||(kid==3))
	//  {
	//printf("Populating original(%d) = %d to %d kid: %d \n",
        //      j+M.srow, stack[j], 0, kid);
        // }
	//#endif

	stack[j] = 0;

	
	//printf("Populating original: k: %d j: %d %f \n",
	//   k+LU(blkcol)(blkrow).scol, 
	//   j+LU(blkcol)(blkrow).srow, M.val(jj));
       
    
      }//end-for over all nnz in column

    //A debug test
    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    //if((kid==2)||(kid == 3))
      {
    printf("DEBUG fill-in pattern:\n");
    for(Int i = 0; i < M.nrow; i++)
      {
	printf("pop i: %d %d kid: %d \n",
               i+M.srow, stack[i], kid);
      }
    printf("\n");
      }
      #endif



  }//end t_populate_col_fill()

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_reduce_col_fill
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k,
   const BASKER_BOOL lower
   )
  {

    const Int my_idx      = S(0)(kid);
    const Int team_leader = find_leader(kid,sl);
    const Int leader_idx  = S(0)(team_leader);
    //Int loop_col_idx      = S(l)(kid);

    //printf("Reduce col fill called, kid: %d leader: %d \n",
    //	   kid, team_leader);


    //If I am not the leader, then get sparsity pattern of B
    if(kid != team_leader)
      {
	Int endblk = (lower)?(LL_size(my_idx)):(l+2);
	
	for(Int blk = l+1; blk < endblk; ++blk)
	  {
//	    ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews; //NDE - warning: unused
	    INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws;
	    //Int      p_sizeL  = LL(leader_idx)(blk).p_size;
	    Int      ws_sizeL = LL(leader_idx)(blk).iws_size;
//	    ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews; //NDE - warning: unused
	    INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
	    const Int ws_size = LL(my_idx)(blk).iws_size;
	    //Int       p_size  = LL(my_idx)(blk).p_size;
	    Int       *color  = &(ws[0]);
	    Int     *pattern  = &(color[ws_size]); 
	    Int     *stack    = &(pattern[ws_size]); //used for fill
	    //Int      brow     = LL(my_idx)(blk).srow;
	    //Int      browL    = LL(leader_idx)(blk).srow;
	    
	    Int *colorL   = &(wsL(0));
	    Int *patternL = &(colorL[ws_sizeL]);
	    Int *stackL   = &(patternL[ws_sizeL]);
	    	    
	    	//over all nnnz found
	    for(Int jj = 0; jj < LL(my_idx)(blk).nrow; ++jj)
	      {
		//if(kid==3)
                // {
		//printf("checking: %d %d %d kid: %d lkid: %d \n",
		// jj+LL(my_idx)(blk).srow, stack[jj], 
		//	 stackL[jj], kid, team_leader);
                // }
		if((stack[jj] != BASKER_MAX_IDX)&&
		   (stackL[jj] != BASKER_MAX_IDX))
		  {
		    //if((kid==2)||(kid==3))
                    // {
		    //     printf("Copy fillto leader: %d %d kid: %d\n",
                    //	   stack[jj], stackL[jj], kid);
		    
                    // }
		    stackL[jj] = min(stackL[jj], stack[jj]);
		    stack[jj]  = stackL[jj];
		    
		  }//if stack(j) != BASKER_MAX_IDX
		else if(stack[jj] != BASKER_MAX_IDX)
		  {
		    //if((kid==2)||(kid==3))
		    //  {
		    //printf("Copy2 fillto leader: %d %d kid: %d\n",
                    //	   stack[jj], stackL[jj], kid);
		    
                    // }

		    stackL[jj] = stack[jj];
		  }
		else
		  {
		    //if((kid==2)||(kid==3))
                    // {
		    //printf("Copy3 fillto leader: %d %d kid: %d\n",
                    //	   stack[jj], stackL[jj], kid);		    
                    //}


		    stack[jj] = stackL[jj];
		  }

	      }//end-for over all nnz
	  }//end-for over all blks
      }//end-if not team_leader
  }//end t_reduce_col_fill()

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_add_orig_fill
  (
   const Int kid,
   const Int lvl,
   const Int l,
   const Int k,
   const BASKER_BOOL lower
   )
  {
  

    if(lower == BASKER_TRUE)
      {
	//printf("===========T ADD ORIG FILL CALLED\n");
	const Int leader_id  = find_leader(kid, l);
	const Int lteam_size = pow(2,l+1);
	const Int L_col      = S(lvl)(leader_id);
	Int L_row      = 0;
	//const Int U_col      = S(lvl)(leader_id); 
	//Int U_row            = LU_size(U_col)-1;
	//Int X_col            = S(0)(leader_id);
        Int X_col            = S(0)(kid);
	Int X_row            = l+1;
    

	//printf("kid: %d leader_kid: %d \n",
	//	   kid, leader_id);
	

	L_row += (kid-leader_id);
	X_row += (kid-leader_id);
	//printf("L_row: %d X_row: %d kid: %d \n", 
	//	   L_row, X_row, kid);
	for( ;
	     L_row < LL_size(L_col);
	     X_row+=(lteam_size), L_row+=(lteam_size))
	  {
	    t_populate_col_fill(kid, 
				L_col, L_row,
				X_col, X_row,
				k, lower);
	  }//for over all row in column
      }//end if(BASKER_TRUE)
    else
      {
	//printf("===========T ADD ORIG FILL CALLED\n");
	const Int leader_id  = find_leader(kid, l);
	//const Int lteam_size = pow(2,l+1);
//	const Int L_col      = S(lvl-1)(leader_id); //NDE - warning: unused
	//Int L_row             = 0;
	//const Int U_col      = S(lvl)(leader_id); 
	//Int U_row            = LU_size(U_col)-1;
	Int X_col            = S(0)(leader_id);
	Int X_row            = l+1;
    
       	//printf("=***== fill MY ID: %d LEADER ID: %d ===** \n",
	// kid, leader_id);

	if(kid == leader_id)
	  {
	    
	    Int bl = l+1;
	    Int A_col = S(lvl)(kid);

	    /*
	    printf("leader_id: %d kid: %d lvl: %d l: %d blk: %d %d \n",
		   leader_id, kid, lvl, l,
		   S(bl)(kid), find_leader(kid,lvl-1));
	    */
	    Int my_row_leader = find_leader(kid, lvl-1);
	    Int my_new_row = 
	      S(bl)(kid) - S(0)(my_row_leader);


	    Int A_row = (lvl==l)?(2):S(bl)(kid)%(LU_size(A_col));
	    if((S(bl)(kid)>14) &&
	       (S(bl)(kid)>LU_size(A_col)) &&
	       (lvl != 1))
	      {
		Int tm = (S(bl)(kid)+1)/16;
		A_row = ((S(bl)(kid)+1)-(tm*16))%LU_size(A_col);
	      }
	   
	    /*
	    printf("TEST---ADD kid: %d A: %d %d new: %d leader: %d %d lvl: %d %d\n",
	    kid, A_col, A_row, my_new_row, my_row_leader, L_col,lvl, bl);
	    */

	    A_row = my_new_row;


	    BASKER_ASSERT((A_row!=(LU_size(A_col)-1)),
			  "ERROR: Selected Lower");
	   
	     t_populate_col_fill(kid, 
				A_col, A_row,
				X_col, X_row,
				 k, lower);

	  }//if not leader
      }
  }//end t_add_orig_fill

}//end namespace BakserNS

#endif //end BASKER_NFACTOR_INC_HPP
