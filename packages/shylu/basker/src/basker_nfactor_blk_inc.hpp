#ifndef BASKER_NFACTOR_BLK_INC_HPP
#define BASKER_NFACTOR_BLK_INC_HPP


//#include "basker_decl.hpp"
#include "basker_matrix_decl.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"
#include "basker_types.hpp"
#include "basker_stats.hpp"

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

      //No longer needed, added in sfactor
      //printf("before workspace init\n");
      //basker->t_init_workspace(kid);
      //printf("after workspace init\n");
      
      //if(kid == 1)
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
 
    //printf("Accessing blk: %d kid: %d  \n", b, kid);
    INT_1DARRAY   ws   = LL(b)(0).iws;
    ENTRY_1DARRAY X    = LL(b)(0).ews;
    Int        ws_size = LL(b)(0).iws_size;
 
   
    Int          bcol  = L.scol;  //begining col
    Int          brow  = L.srow;  //begining row 
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
    
    //Why did we need this?
    Int col_idx_offset = M.nnz;

    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
    printf("=======NFACTOR BLK INC LVL========\n");
    #endif
    //printf("test one ws_size: %d \n", ws_size);
    
    //Note:
    Int *color    = &(ws(0));
    Int *pattern  = &(color[ws_size]);
        
    maxindex = BASKER_MAX_IDX;
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
    
    //for each column    
    for(k = 0; k < M.ncol; ++k)
    //for(k = 0; k < 36; ++k)
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
	      t_local_reach_inc_lvl(kid,0,0,j,&top);
	      
          }//end for() each nnz in column
	  xnnz = ws_size - top;

          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("xnnz: %d ws_size: %d top: %d \n", 
                 xnnz, ws_size, top);
          #endif
             
	  t_back_solve_inc_lvl(kid,0,0,k,top,xnnz);

          maxv = 0.0;
          for(i = top; i < ws_size; i++)
            {
	      j = pattern[i];
	      t = gperm(j+brow);
	   
	      value = X(j);
	   
	      absv = EntryOP::approxABS(value);
    
	      if(t == BASKER_MAX_IDX)
                {
                  lcnt++;
                  //if(absv > maxv)
		  if(EntryOP::gt(absv,maxv))
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
	      maxindex = k;
	      pivot    = X(k);
	      //printf("sym pivot: %f \n", pivot);
	    }
      
          ucnt = ws_size - top - lcnt +1;
    
	  if((maxindex == BASKER_MAX_IDX) || (pivot == 0))
            {
	      cout << endl << endl;
	      cout << "---------------------------"<<endl;
              cout << "Error: Matrix is singular, blk" << endl;
              cout << "MaxIndex: " << maxindex << " pivot " 
                   << pivot << endl;
              cout << "lcnt: " << lcnt << endl;
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
	      printf("\n\n");
	      printf("----------------------\n");

              newsize = lnnz * 1.1 + 2 *M.nrow + 1;
              printf("b: %d Reallocing L oldsize: %d current: %d count: %d newsize: %d \n",
                     b, llnnz, lnnz, lcnt, newsize);

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
		  thread_array(kid).error_blk = b;
		  thread_array(kid).error_info = newsize;
		}

            }
          if(unnz+ucnt > uunnz)
            {

	      printf("\n\n");
	      printf("-------------------\n");

              newsize = uunnz*1.1 + 2*M.nrow+1;
              printf("b: %d Reallocing U oldsize: %d newsize: %d \n",
                     b, uunnz, newsize);

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
		  thread_array(kid).error_blk = b;
		  thread_array(kid).error_info = newsize;
		}

            }

          L.row_idx(lnnz) = maxindex;
          L.val(lnnz)     = (Entry) 1.0;
	  L.inc_lvl(lnnz) = INC_LVL_TEMP(maxindex+brow);
	  INC_LVL_TEMP(maxindex+brow) = BASKER_MAX_IDX;
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
              printf("j: %d t: %d \n", j, t);
              #endif            

	      
	      //For very ill problems, will need to clear
	      //if(X(j) != 0)
                {
		  if(t != BASKER_MAX_IDX)
                    {
                      if(t < (k+brow))
                        {
			  
			  U.row_idx(unnz) = t-brow;
			  U.val(unnz) = X(j);
			  unnz++;
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
		      
		      L.inc_lvl(lnnz) = INC_LVL_TEMP(j+brow);
		      INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
                      lnnz++;
                    }
                }//end if() not 0

              #ifdef BASKER_DEBUG_NFACTOR_BLK
              printf("Zeroing element: %d \n", j);
              #endif

	      X(j) = 0;	  
            }//end if(x[i] != 0)

          //Fill in last element of U
	  U.row_idx(unnz) = k;
	  U.val(unnz)     = lastU;
          unnz++;

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
	      //printf("before offdiag \n");
	      t_back_solve_offdiag_inc_lvl(kid, L.srow,
				   b, blk_row,
				   b, blk_row,
				   k, col_idx_offset,
				   U.val, U.row_idx,
		       U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k),
				   BASKER_TRUE);
	    
	      //Move these factors into Local Ls
	      //printf("before t_move_L\n");
	      t_move_offdiag_L_inc_lvl(kid,
			       b, blk_row,
			       b, blk_row,
			       k, pivot);
	         
	    }//end over all diag
	  #endif

	  
	  //NOT PRUNE
	  //t_prune(kid,0,0,k,maxindex);

	}//end for() over all columns

    L.nnz = lnnz;
    U.nnz = unnz;

    return 0;
  }//end t_nfactor_blk_inc_lvl()

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
    
    //#ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach_inc_lvl, L: %d %d  X: %d %d, kid: %d \n",
	   b, 0, wsb, l, kid);
    //#endif

    start    = -1;
    head     = 0;
    stack[0] = j;

    Int inc_lvl = 0;
    Int pop_top = *top-1; //Where we want to add

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
    printf("end sample top: %d p:%d %d \n",
	   pop_top, pattern[pop_top], pattern[pop_top+1]);
    #endif

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
       
	INC_LVL_TEMP(j+L.srow) = 0;
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
	    
		BASKER_ASSERT((*top-1) > 0, "Top pass pattern\n");
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

    //printf("\n\t_back_solve_selective,kid: %d \n",
    //kid);

    Int brow = L.srow;

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
 
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;
    for(pp = 0; pp < xnnz; pp++)
      {
	
	j = pattern[top1];
	t = gperm(j+brow);

	//printf("backsolve column: %d \n",
	//     j,t);

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
		    continue;
		  }
	
		
		//printf("Updateing row: %d with value: %f %f \n",
		//     L.row_idx(p),X[L.row_idx[p]], L.val[p]*xj);


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

    //B.info();
    //B.print();

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
    if(kid == 8 )
      {
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
      }
    #endif

   
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);

    const Int    brow  = L.srow;
    const Int    bcol  = L.scol;
    const Int    llnnz = L.nnz;
          Int    lnnz  = L.col_ptr(k);

    if((p_size) > (llnnz-lnnz))
      {
	printf("-Warning, Need to remalloc L: %d %d kid: %d current size: %d used_size: %d  addition: %d \n",
	       blkcol, blkrow, kid, llnnz,lnnz,p_size  );
	BASKER_ASSERT(0==1, "REALLOC LOWER BLOCK\n");
	
      }

    for(Int i = 0; i < p_size; i++)
      {
	Int j = pattern[i];
	Int t = gperm(j+brow);
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, kid: %d j: %d val: %f lnnz: %d inc: %d \n",
	       kid, j, X[j]/pivot, lnnz, stack[j]);
	#endif

	color[j] = 0;
	
	L.row_idx(lnnz) = j;
	L.inc_lvl(lnnz) = stack[j];
	stack[j] = BASKER_MAX_IDX;
	L.val(lnnz) = EntryOP::divide(X(j),pivot);
	X(j) = 0;

	lnnz++;
      }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //	   kid, k-bcol+1, lnnz);
    
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
	      cout << endl << endl;

	      cout << "---------------------------"<<endl;
	     
              cout << "Error: Matrix is singular, blk" << endl;
              cout << "MaxIndex: " << maxindex << " pivot " 
                   << pivot << endl;
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
              printf("b: %d Reallocing L oldsize: %d current: %d count: %d newsize: %d \n",
                     b, llnnz, lnnz, lcnt, newsize);
            }
          if(unnz+ucnt > uunnz)
            {

	      printf("\n\n");
	      printf("-------------------\n");

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

  //Uses glb idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::t_back_solve_old
  (Int kid, 
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz)
  {

    Int            b = S[lvl][kid];
    Int          wsb = S[0][kid];
    BASKER_MATRIX &L = LL[b][0];
    #ifdef BASKER_2DL
    INT_1DARRAY   ws = LL[wsb][l].iws;
    ENTRY_1DARRAY X  = LL[wsb][l].ews;
    Int      ws_size = LL[wsb][l].iws_size;
    #else
    INT_1DARRAY   ws = thread_array[kid].iws;
    ENTRY_1DARRAY  X = thread_array[kid].ews;
    Int      ws_size = thread_array[kid].iws_size;
    #endif
   
    //#ifdef BASKER_DEBUG_NFACTOR_BLK
    /*
    if(k < 6)
      {
    printf("t_back_solve, kid: %d L: %d %d X: %d ws: %d \n",
	   kid, b, 0, kid, wsb);
      }
    */
    //#endif

    Int brow = L.srow;

    Int *color   = &(ws[0]);
    Int *pattern = &(color[ws_size]);
 
    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;
    for(pp = 0; pp < xnnz; pp++)
      {
	j = pattern[top1];

	#ifdef BASKER_2DL
	color[j-brow] = 0;
	#else
	color[j] = 0;
	#endif
        t = gperm[j];
	
	//#ifdef BASKER_DEBUG_NFACTOR_BLK
	/*
	if(k < 6)
	  {
	printf("t_back_solve, kid: %d j: %d t: %d \n",
	       kid, j, t);
	printf("t_back_solve, kid: %d j: %d xj: %f \n",
	       kid, j, X[j-brow]);
	  }
	*/
	//#endif
	

	//NOTE:::::: SHOULD NOT NEED --- IF WE DO NEED THAN TOO COSTLY
	//Could we make this check smaller in the 2dl case??
	/*
	if(t != L.max_idx && (t >= L.scol) && 
	   (t < (L.scol+L.ncol)))
	*/
	//if(t!=L.max_idx)
	if(t!=BASKER_MAX_IDX)
          {

	    #ifdef BASKER_2DL
	    xj = X[j-brow];
	    #else
            xj = X[j];
	    #endif
	   
	    // #ifdef BASKER_DEBUG_NFACTOR_BLK
	    /*
	    if(k < 6)
	      {
		printf("Updating column: %d  with %f \n", t, xj);
	      }
	    */
            //#endif

            //Get rid of these temp variables
            Int local_offset = L.scol; //L.srow
            p2 = L.col_ptr[t+1-local_offset];
            p = L.col_ptr[t-local_offset]+1;
            
            for( ; p < p2; p++)
              {

		/*
		if(k < 6)
		  {
		    printf("Updateing with value: %f %f ri: %d Lval: %f \n",
			   X[L.row_idx[p]-brow], L.val[p]*xj, L.row_idx[p], L.val[p] );
		  }
		*/

		#ifdef BASKER_DEBUG_NFACTOR_BLK
		#ifdef BASKER_2DL
		printf("Updating row: %d  with value: %f %f \n",
		       L.row_idx[p], X[L.row_idx[p]-brow], L.val[p]*xj);
		#else
		printf("Updateing with value: %f %f \n",
		       X[L.row_idx[p]], L.val[p]*xj);
		#endif
		#endif
		#ifdef BASKER_2DL
		X[L.row_idx[p]-brow] -= L.val[p]*xj;
		#else
                X[L.row_idx[p]] -= L.val[p] *xj;
		#endif
              }//end for() over each nnz in the column
          }//end if() not permuted
	top1++;
      }//end for() over all nnz in LHS
    return 0;
  }//end t_back_solve()


  //used global idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_move_offdiag_L_old
  (Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot)
  {
    BASKER_MATRIX &L    = LL[blkcol][blkrow];
    #ifdef BASKER_2DL
    INT_1DARRAY   ws    = LL[X_col][X_row].iws;
    ENTRY_1DARRAY X     = LL[X_col][X_row].ews;
    Int         ws_size = LL[X_col][X_row].iws_size;
    Int         p_size  = LL[X_col][X_row].p_size;
    #else
    INT_1DARRAY  ws = thread_array[kid].iws;
    ENTRY_1DARRAY X = thread_array[kid].ews;
    Int     ws_size = thread_array[kid].iws_size;
    Int     p_size  = 3;
    #endif

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
    #endif

   
    Int *color   = &(ws[0]);
    Int *pattern = &(color[ws_size]);

    Int    brow  = L.srow;
    Int    bcol  = L.scol;
    Int    llnnz = L.nnz;
    Int    lnnz  = L.col_ptr[k-bcol];

    if((p_size) > (llnnz-lnnz))
      {
	printf("-Warning, Need to remalloc L: %d %d kid: %d current size: %d used_size: %d  addition: %d \n",
	       blkcol, blkrow, kid, llnnz,lnnz,p_size  );
	
      }

    for(Int i = 0; i < p_size; i++)
      {
	Int j = pattern[i];
	Int t = gperm[j];
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("L-Moving, kid: %d j: %d val: %f lnnz: %d \n",
	       kid, j, X[j-brow]/pivot, lnnz);
	#endif

	color[j-brow] = 0;

	L.row_idx[lnnz] = j;
	L.val[lnnz] = X[j-brow]/pivot;
	X[j-brow] = 0;

	#ifdef BASKER_INC_LVL
	L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
	#endif

	lnnz++;
      }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //	   kid, k-bcol+1, lnnz);
    L.col_ptr[k-bcol+1] = lnnz;

    LL[X_col][X_row].p_size = 0;

    return 0;
  }//end t_offdiag_mov_L()

  //old global idx for local blks
  //Bgood done.
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag_old
  (Int kid, 
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k , Int &view_offset,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option)
  {
    //Note:  need to add support for offdiag permuation

    BASKER_MATRIX &L            = LL[blkcol][blkrow];
    //BASKER_MATRIX_VIEW &B       = AL[blkcol][blkrow];
    //B.init_offset(k, view_offset);
    BASKER_MATRIX &B            = ALM[blkcol][blkrow];

   


    #ifdef BASKER_2DL
    INT_1DARRAY   ws            = LL[X_col][X_row].iws;
    ENTRY_1DARRAY X             = LL[X_col][X_row].ews;
    Int         ws_size         = LL[X_col][X_row].iws_size;
    #else
    INT_1DARRAY   ws            = thread_array[kid].iws;
    ENTRY_1DARRAY X             = thread_array[kid].ews;
    Int         ws_size         = thread_array[kid].iws_size;
    #endif

    Int          nnz            = LL[X_col][X_row].p_size;
    Int          brow           = L.srow;
    Int          bcol           = L.scol;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_diag, kid: %d blkcol: %d blkrow: %d \n",
	   kid, blkcol, blkrow);
    printf("t_back_solve_diag, kid: %d Xcol: %d Xrow: %d \n",
	   kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d ws: %d starting psize: %d \n",
	   kid,ws_size, nnz);
    #endif
    // B.info();
    //B.print();

    Int *color =   &(ws[0]);
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
	for(Int i = B.col_ptr[k-bcol]; i < B.col_ptr[k+1-bcol]; i++)
	  {
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    //Bgood(remove)
	    //printf("t_back_solve_diag, kid: %d i: %d g: %d\n",
	    //	   kid, i, B.good(i));
	    #endif


	    //Bgood(remove)
	    //Note: Come back and slow down
	    //if(B.good(i)==L.max_idx)
	    //if(B.good(i) == BASKER_MAX_IDX)
	    //  {
	    //view_offset = i;
	    //	break;
	    // }
	    #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d, add A, kid: %d psize:%d \n",
		   kid, nnz);
	    printf("t_back_solve_diag, kid: %d A(%d) %f \n",
		   kid, B.row_idx(i), B.val(i));
	    #endif
	    Int j = B.row_idx(i);
	    color[j-brow] = 1;
	    X[j-brow] = B.val(i);
	    pattern[nnz++] = j;
	  }//over all nnz in subview
      }//end if preload

    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %d \n",
	   kid, x_size);
    #endif
    for(Int i = 0 ; i < x_size; i++)
      {
	Int k =   x_idx[i+x_offset];
	Entry xj = x[i+x_offset];
	//printf("kid: %d bcol: %d k: %d \n",
	//   kid, bcol, k);
	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
	     kid, k, L.col_ptr[k-bcol], L.col_ptr[k-bcol+1]);
	#endif
	for(Int j = L.col_ptr[k-bcol]; 
	    j < L.col_ptr[k-bcol+1]; j++)
	  {
	    Int jj = L.row_idx[j];
            #ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj-brow]);
	    #endif

	    if(color[jj-brow] != 1)
	      {
		color[jj-brow] = 1;
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

	    #ifdef BASKER_DEBUG_NFACTOR_BLK

	    printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
	    kid, jj,X[jj-brow], L.val[j], xj);
	    #endif 

	    X[jj-brow] -= L.val[j]*xj;
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
    LL[X_col][X_row].p_size = nnz;
    #endif


    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("kid: %d all values: \n", kid);
    for(Int i = 0; i < nnz; i++)
      {
	Int jj = pattern[i];
	printf("X[%d] = %f , kid: %d  \n",
	       jj, X[jj-brow], kid);

      }
    printf("\n\n");
    #endif



    return 0;
  }//end t_offdiag_back_solve();

  
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_dense_back_solve_offdiag_inc_lvl
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
    const Int    brow           = L.srow;
    const Int    bcol           = L.scol;
  
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

    //Need to find a way not to prefill everytime
    for(Int i = 0; i < B.nrow; i++)
      {
	stack[i] = BASKER_MAX_IDX;
      }

    Int preload_fill = 0;
    //Preload with A
    if(A_option == BASKER_TRUE)
      {
	//#ifdef BASKER_DEBUG_NFACTROR_BLK
	printf("t_back_solve, A_OPTION TRUE \n");
	//#endif

	for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
	  {
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
       
	    stack[j] = 0;
	    
	    preload_fill++;

	  }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
	   kid, x_size);
    #endif

    for(Int i = 0 ; i < x_size; ++i)
      {
	const Int k    = x_idx(i+x_offset);
	const Entry xj = x(i+x_offset);

	#ifdef BASKER_DEBUG_NFACTOR_BLK
	printf("t_back_solve_diag, kid: %d  k: %d %g  x_size: %d [%d %d] \n",
	       kid, k, xj, x_size,  L.col_ptr[k], L.col_ptr[k+1]);
	#endif
	
	
	#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	printf("LVL_TEMP[%d] = %d, %d kid: %d continue? \n",
	       k+pbcol, x_fill(i+x_offset), 
	       Options.inc_lvl, kid);
	#endif

	if(preload_fill != 0)
	  {
	    if(x_fill(i+x_offset)+1 > Options.inc_lvl)
	  {

	    #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	    printf("continued skip because too large\n");
	    #endif
	    continue;
	  }
	  }
	

	for(Int j = L.col_ptr(k); 
	    j < L.col_ptr(k+1); j++)
	  {
	    const Int jj = L.row_idx(j);
            //#ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_b_solve_d, kid: %d j: %d color: %d \n",
		   kid, jj, color[jj]);
	    //#endif


	    //#ifdef BASKER_DEBUG_NFACTOR_BLK
	    printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
		    kid, jj,X[jj], L.val[j], xj);
	    // #endif 
	     
	     //#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("L.inc_lvl[%d]: %d %d kid: %d \n",
		    j, L.inc_lvl(j), x_fill(i+x_offset), kid);
	     //#endif

	     Int temp_cal = L.inc_lvl(j) + x_fill(i+x_offset)+1;
	     if(stack[jj] = BASKER_MAX_IDX)
	       {
		 stack[jj] = temp_cal;
	       }
	     else
	       {
		 stack[jj] = min(temp_cal, stack[jj]);
	       }
	     
	     //#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("Assigned inc_lvl(%d) = %d \n",
		    jj, stack[jj]);
	     //#endif

	     if((stack[jj]) > Options.inc_lvl)
	       {
		 //#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
		 printf("continue,already used Linc(%d):%d %d\n",
			j, L.inc_lvl(j), x_fill(i+x_offset));
		 //#endif
		 continue;
	       }

	     //#ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	     printf("VALUE: before %f %f %f AFTER: %f\n",
		    X(jj), L.val(j), xj, X(jj)-L.val(j)*xj);
	     //#endif
	     
	     X(jj) -= L.val(j)*xj;

	  }//over all nnz in column
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
    const Int   p_size  = LL(X_col)(X_row).p_size;
   
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_move_offdiag_L, kid: %d L %d % X %d %d p_size: %d \n",
	   kid, blkcol,blkrow, X_col, blkrow,  p_size);
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);

    const Int    brow  = L.srow;
    const Int    bcol  = L.scol;
    const Int    llnnz = L.nnz;
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
	if(X(j)!=0)
	  {

	Int t = gperm(j+brow);
	
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
  }//end t_dense_offdiag_mov_L_inv_lvl()


}//end namespace BakserNS

#endif //end BASKER_NFACTOR_INC_HPP
