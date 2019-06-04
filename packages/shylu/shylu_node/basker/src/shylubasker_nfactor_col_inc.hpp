#ifndef SHYLUBASKER_NFACTOR_COL_INC_HPP
#define SHYLUBASKER_NFACTOR_COL_INC_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_thread.hpp"
#include "shylubasker_util.hpp"

#include "shylubasker_nfactor_blk_inc.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_sep2_inc_lvl
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif
    
    Basker<Int,Entry,Exe_Space> *basker;
    Int lvl;
    
    kokkos_nfactor_sep2_inc_lvl()
    {}
    
    kokkos_nfactor_sep2_inc_lvl(
                Basker<Int,Entry,Exe_Space>* _basker, Int _lvl)
    {
      basker = _basker;
      lvl = _lvl;
    }

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const  
    #endif
    {
      #ifdef BASKER_KOKKOS
      Int kid = basker->t_get_kid(thread);
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #else
      Int team_leader = 0; //Note: come back and fix
      #endif
      
      //if(kid < 8)
      //if(kid > 11 && kid < 16)
	{
      #ifdef BASKER_KOKKOS
       basker->t_nfactor_sep2_inc_lvl(kid, lvl, team_leader, thread);
      #else
      
      #endif
	}
    }//end operator ()
  };//end col_factor_funct

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_nfactor_sep2_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int team_leader,
   const TeamMember &thread
   )
  {

    const Int U_col =  S(lvl)(kid);
    Int U_row       =  0;
	  
    //const Int scol = LU(U_col)(U_row).scol;
    //const Int ecol = LU(U_col)(U_row).ecol;

    int upper_error = BASKER_MAX_IDX;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n  LVL=%d  ----- kid: %d --\n\n",
	   lvl, kid);
    #endif
    
    //Do all domains (old sublevel 0)
    //If this works well we can move this for into the function


    //for(Int k = 0; k < 1; ++k)
     for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("UPPER, kid: %d k: %d \n",
	       kid, k);
	#endif


	upper_error =
	t_upper_col_factor_inc_lvl(kid, team_leader, 
				   lvl, 0, 
				   k,
				   BASKER_FALSE);

	if(upper_error == BASKER_ERROR)
	  {
	    //Need a way to say all should stop!
	    //How to do that nice, maybe an atomic
	    //Cannot return because the communication will hang
	    //printf("kid: %d break error \n", kid);
	    break;
	  }
	
      }//over all columns / domains / old sublevel 0

     #ifdef BASKER_DEBUG_NFACTOR_COL2
     printf("\n\n\n done with UPPER, kid: %d \n\n\n", kid);
     #endif

     
    //------Need because extend does not 
    //-------------Barrier--Between Domains-------------
    Int error_leader = find_leader(kid, lvl-1);
    //printf("kid: %d eleader: %d test: %d \n", 
    //	   kid, error_leader, find_leader(kid,lvl-1));
    if(upper_error == BASKER_ERROR)
      {
	//printf("kid: %d set error leader: %d  \n", 
	//     kid, error_leader);
	basker_barrier.ExitSet(error_leader, BASKER_TRUE);
      }
    Int my_leader = find_leader(kid, 0);
    Int b_size    = pow(2,1);
    //barrier k = 0 usedl1
   
    t_basker_barrier_inc_lvl(thread,kid,my_leader,
		     b_size, 0, LU(U_col)(U_row).scol, 0);
    //printf("1 kid: %d  error_leader: %d lvl: %d  \n", kid, error_leader, lvl);
    BASKER_BOOL error_flag = BASKER_FALSE;
    basker_barrier.ExitGet(error_leader, error_flag);
    if(error_flag == BASKER_TRUE)
      {
	return;
      }
    //printf("2 kid: %d \n", kid);
    
    //DEBUG
    /*
    if(lvl == 2)
      {
	return;
      }
    */
    //return;

    //----------------Sep level upper tri-------------
    for(Int l = 1; l < (lvl); ++l)
      {
	
	//for(Int k = 2; k < 3; ++k)
	for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
	  {
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("\n\nSep, upper update, kid: %d k=%d %d\n\n",
		   kid, k,k+LU(U_col)(U_row).scol);
	    #endif
	    
	    t_add_extend_inc_lvl(thread, kid,lvl,l-1, k, 
			 LU(U_col)(U_row).scol, 
			 BASKER_FALSE);
	    
	    //where to start again
	    if(kid%((Int)pow(2,l)) == 0)
	      {
		my_leader = find_leader_inc_lvl(kid,l);
		b_size    = pow(2,l+1);
		
    #ifdef BASKER_DEBUG_NFACTOR_COL2
		printf("\n\n\n SEP UPPER, kid: %d \n\n",
		       kid);
		#endif
		
		//if((lvl != 2) || 
		//  ((kid == 0) && (k < 4)))
		  {
		Int upper_error =
		t_upper_col_factor_inc_lvl(kid, team_leader, 
				   lvl, l, 
				   k,
				   BASKER_FALSE);
		
		if(upper_error == BASKER_ERROR)
		  {
		    //Again need a nice way to exit here
		    printf("ShyLU Basker Error: BREAK CALLED\n");
		    break;
		  }
		  }
	      }//if correct kid to do this sublevels upperfactor
	    
	  }//over all columns
      }//for - over all sublevel 1...lvl-2


    //DEBUG
    /*
    if(lvl == 2)
      {
	return;
      }
    */

    //printf("Done with upper seps kid: %d \n", kid);
   
    //---------Lower Factor (old sublevel lvl-1)-------
    
    my_leader = find_leader_inc_lvl(kid, lvl-1);
    b_size    = pow(2,lvl);
    // printf("[3] barrier test, kid: %d leader: %d b_size: %d lvl: %d \n",
    //	   kid,  my_leader, b_size, lvl);
    t_basker_barrier_inc_lvl(thread, kid, my_leader,
		     b_size, 7, LU(U_col)(U_row).scol, 0);

    #ifdef BASKER_DEBUG_NFACTOR_COL_INC
    if(kid == 0)
      {
    printf("\n\n======= LOWER, KID: %d [%d %d] ======= \n\n", 
	   kid, LU(U_col)(U_row).scol,
	   LU(U_col)(U_row).scol + LU(U_col)(U_row).ncol);
      }
    #endif
    //printf("\n\n lvl: %d kid: %d  \n\n", lvl, kid);
    //if(lvl < 2)
      {
	//for(Int k=0; k < 1; ++k)
     for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("\n*******lower_update, kid: %d k: %d \n",
	       kid, k+LU(U_col)(U_row).scol);
	#endif

	//printf("test: %d \n", LU(U_col)(U_row).scol);
       
	t_add_extend_inc_lvl(thread, kid,lvl,lvl-1, k,
		     LU(U_col)(U_row).scol,
		     BASKER_TRUE);
	Entry pivot = 0;
	if((kid%(Int)(pow(2,lvl))) == 0)
	  {

	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("\n lower factor, kid: %d k: %d \n",
		   kid, k+LU(U_col)(U_row).scol);
	    #endif
	    
	    t_lower_col_factor_inc_lvl(kid, team_leader, 
			       lvl, lvl-1, 
			       k, pivot);
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("\n lower factor return, kid: %d k: %d \n",
		   kid, k+LU(U_col)(U_row).scol);
	    #endif

	  }
	
	
	//nee barrier if multiple thread uppdate
	//thread.team_barrier();
	my_leader = find_leader_inc_lvl(kid, lvl-1);
	b_size    = pow(2,lvl);
	//printf("barrier test, leader 4: %d b_size: %d lvl: %d \n",
	//     my_leader, b_size, lvl);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 8, k, lvl-1);
	
	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("lower diag factor, kid: %d k: %d \n",
	       kid, k);
	#endif
	

	//printf("\n\n===========AFTER LOWER COL FACTOR =======\n\n");

        //if((kid == 0)||(kid==1))
	t_lower_col_factor_offdiag2_inc_lvl(kid, lvl, lvl-1, k, pivot);

	//printf("\n\n===========AFTER LOWER UPDATE FACTOR =======\n\n");

	//thread.team_barrier();
	my_leader = find_leader_inc_lvl(kid, lvl-1);
	b_size    = pow(2,lvl);
	//printf("barrier test 5, leader: %d b_size: %d lvl: %d \n",
	//     my_leader, b_size, lvl);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
		     b_size, 9, k, lvl-1);


	//Need to clean INC_LVL_HERE
	if(Options.same_pattern == BASKER_FALSE)
	  {
	t_lower_col_factor_offdiag2_cleanup_inc_lvl(kid, lvl, lvl-1, k);
	  }

      }
      }
    
   
  }//end t_nfactor_sep2_inc_lvl


  //This is the most complex part of the the code
  //It is also the life-blood that makes it work!
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_add_extend_inc_lvl
  (
   const TeamMember &thread,
   const Int kid, 
   const Int lvl,
   const Int l,
   const Int k,
   const Int k_offset,
   const BASKER_BOOL lower
   )
  {
    Int my_leader = find_leader_inc_lvl(kid,l);
    Int b_size    = pow(2,l+1);

    //===============Symbolic SPMV=================//
    //loop over all sublevels to perform
    //off-diag smultiple, reduce, etc
    if(Options.same_pattern == BASKER_FALSE)
      {
    for(Int sl = 0; sl <=l; ++sl)
      {
		
	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_FALSE)
	  {
	printf("extend up, kid: %d sl: %d l: %d lvl: %d \n",
	       kid, sl, l, lvl);
	  }
	else
	  {
	printf("extend low, kid: %d sl: %d l: %d lvl: %d \n", 
	       kid, sl, l, lvl);
	  }
	#endif

	//printf("sl=%d kid: %d\n", sl, kid);
	
		
	//This will do the correct symb-spmv
	//printf("=======SFACTOR KID:%d ======\n", kid);
	t_upper_col_ffactor_offdiag2_inc_lvl(kid, lvl, sl,
					     l, k, lower);


	//if(sl == l)
	  {
	    //printf("========FILL IN KID: %d =======\n", kid);
	    t_add_orig_fill(kid, lvl, l, k, lower); 
	    //printf("=======AFFERT FILL_IN KID: %d ====\n", kid);
	  }
	 //printf("\n\n");
	
	//Barrier--Start
	my_leader = find_leader_inc_lvl(kid,sl);
	b_size    = pow(2,sl+1);
	//printf("Barrier1 kid: %d Leader: %d size: %d \n",
	//     kid, my_leader, b_size);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 1, k+k_offset, sl);
	//Barrier--End

        //Upward scan
	//Copy to leader
	//printf("\n\n");
	if(kid%((Int)pow(2,sl))==0)
	  {
	    //printf("=========REDUCE KID: %d \n", kid);
	    t_reduce_col_fill(kid,
			      lvl, sl, l, k, lower);
	  }
	
	//printf("\n\n\n");
	//printf("Done with sreduce. kid: %d \n", kid);

	//Barrier--Start
	//printf("Barrier2 kid: %d Leader: %d size: %d \n",
	//     kid, my_leader, b_size);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 2, k+k_offset, sl);
	
	//Barrier--End

      }//end-for sl (sublevel)
    //Backwards scan to reduce fill-in pattern
    for(Int sl = l-1; sl >=0; --sl)
      {
        if(kid%((Int)pow(2,sl))==0)
	  {
	    //printf("=========SCAN DOWN REDUCE KID: %d \n", kid);
	    t_reduce_col_fill(kid,
			      lvl, sl, l, k, lower);
	  }

        //Barrier--Start
        t_basker_barrier_inc_lvl(thread, kid, my_leader,
                                 b_size, 3, k+k_offset, sl);
	//Barrier--Ent
      }
      }//if same_pattern == false


    //================Numeric SPMV=================//
    //loop over sublevels to perform 
    //off-dig multiple, reduce, etc
    for(Int sl = 0; sl <= l; ++sl)
      {
	
	
	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_FALSE)
	  {
	printf("extend up, kid: %d sl: %d l: %d lvl: %d \n",
	       kid, sl, l, lvl);
	  }
	else
	  {
	printf("extend low, kid: %d sl: %d l: %d lvl: %d \n", 
	       kid, sl, l, lvl);
	  }
	#endif

	//Add barrier here that was not here before
	//Barrier--Start
	t_basker_barrier_inc_lvl(thread,kid, my_leader,
				 b_size, 4, k+k_offset,sl);

		
	//This will do the correct spmv
	if(Options.same_pattern == BASKER_FALSE)
	  {
	t_upper_col_factor_offdiag2_inc_lvl(kid, lvl, sl,l, k, lower);
	  }
	else
	  {
	    
	    //printf("t_upper_col_factor, kid: %d \n", kid);
	 t_upper_col_factor_offdiag2_same_pattern_inc_lvl(kid, lvl, sl,l, k, lower);
	  }

	
	//Barrier--Start
	my_leader = find_leader_inc_lvl(kid,sl);
	b_size    = pow(2,sl+1);
	//printf("[1] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 5, k+k_offset, sl);

	//printf("DONE SPMV KID: %d \n", kid);
	//Barrier--End

	if(kid%((Int)pow(2,sl))==0)
	  {
	    //printf("======COPY ATOMIC kid: %d ========\n", kid);
	    if(Options.same_pattern == BASKER_FALSE)
	      {
	    t_dense_blk_col_copy_atomic2_inc_lvl(kid, my_leader,
					 lvl, sl, l, k, lower);
	      }
	    else
	      {
		//printf("t_col_copy, kid: %d \n", kid);
		t_same_pattern_col_copy_inc_lvl(kid, lvl, 
						sl,l,k, lower);
	      }
	  }

	//Barrier--Start
	//printf("[2] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 6, k+k_offset, sl);
	
      }//over all sublevels

    //printf("TEST--LEADER me: %d leader: %d \n", 
    //	   kid, my_leader);
    //JDB TESTED
    if(kid == my_leader)
      {
	if(Options.same_pattern == BASKER_FALSE)
	  {
    t_dense_copy_update_matrix2_inc_lvl(kid, my_leader, lvl, l, k);
	  }
	else
	  {
	    //printf("t_copy matrix, kid: %d \n", kid);
	    t_same_pattern_update_matrix_inc_lvl(kid,my_leader,lvl,l,k);
	    
	  }
      }


  }//end t_add_add_inc_lvl

 


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_basker_barrier_inc_lvl
  (
   const TeamMember &thread,
   const Int my_kid,
   const Int leader_kid, 
   const Int size,
   const Int function_n,
   const Int k, 
   const Int l
   )
  {
    //printf("before call. lkid: %d kid: %d task: %d size: %d k: %d \n",
    //	       leader_kid, my_kid, function_n, size, k);
    
     if(size <= thread.team_size())
      {
	thread.team_barrier();
      }
    else
      {
		
	basker_barrier.BarrierDomain
	  (
	   leader_kid,
	   my_kid, 
	   function_n,
	   size, 
	   k, l
	   );
	
	
      }
  }//end t_basker_barrier_inc_lvl

  //We need an easier and faster way to do this.  
  //Could get very big
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int Basker<Int,Entry,Exe_Space>::find_leader_inc_lvl
  (
   Int kid, 
   Int l
   )
  {
    l = l+1;
    Int my_token = S(l)(kid);
    Int my_loc = kid;
    while((my_loc > 0))
      {
	my_loc--;
	if(S(l)(my_loc) != my_token)
	  {
	    my_loc++;
	    break;
	  }
      }

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("find_leader, kid: %d l: %d leader: %d \n",
	   kid, l, my_loc);
    #endif
    return my_loc;

  }//end find_leader_inc_lvl()

  //local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_upper_col_factor_inc_lvl
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l,
   Int k, 
   BASKER_BOOL sep_flg
   )
  {

    //printf("UPPER KID: %d K: %d \n",
    //	   kid, k);


    //Get needed variables
    const Int L_col = S(l)(kid);
//    const Int L_row = 0; //NDE - warning: unused 
    const Int U_col = S(lvl)(kid);
    
    Int my_row_leader = find_leader(kid,lvl-1);
    //Int my_new_row = 
    //  L_col - S(0)(my_row_leader);
    Int U_row = L_col - S(0)(my_row_leader);
  
    /*
    Int U_row = (lvl==1)?(kid%2):S(l)(kid)%LU_size(U_col);  
    if((L_col > 14) &&
       (L_col > LU_size(U_col)) &&
       (lvl != 1))
      {
	Int tm = (L_col+1)/16;
	U_row = ((L_col+1)-(tm*16))%LU_size(U_col);
      }
    */
    //printf("U-upper.  kid: %d U: %d %d  leader: %d %d new: %d lvl: %d l %d \n",
    //	   kid, U_col, U_row, my_row_leader, L_col,my_new_row, lvl, l);
    
    //JDB TEST
    //PASS TEST
    //U_row = my_new_row;


    const Int X_col = S(0)(kid);
    const Int X_row = l; //X_row = lower(L)
    //const Int col_idx_offset = 0; //we might be able to remove
  
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
    printf("kid %d, upper using L: %d %d  U: %d %d  X %d %d\n",
	   kid, L_col, L_row, U_col, U_row, X_col, X_row);
    #endif
    //end get needed variables//

    //BASKER_MATRIX        &L = LL(L_col)(L_row); //NDE - warning: unused L
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    
    //Ask C++ guru if this is ok
    BASKER_MATRIX        *Bp;
    if(l == 0)
      {
        Bp = &(AVM(U_col)(U_row));
      }
    else
      {
	Bp = &(thread_array[kid].C);
      }
    BASKER_MATRIX    &B = *Bp;
    //if(kid ==0)
    // {
      //	printf("11====KID ===11\n");
      //B.print();
    // }
    //B.print();

    INT_1DARRAY ws     = LL(X_col)(X_row).iws;
    const Int ws_size  = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY X    = LL(X_col)(X_row).ews;
  
    const Int brow = U.srow;
    //const Int bcol = U.scol;


    
    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);
    Int *stack     = &(pattern[ws_size]);
    

    /*
    auto color = Kokkos::subview(ws, std::make_pair((Int)0,ws_size));
    auto pattern = Kokkos::subview(ws, std::make_pair(ws_size,2*ws_size));
    auto stack = Kokkos::subview(ws, std::make_pair(2*ws_size, 3*ws_size));
    */


    Int j, t, xnnz;
    Int top = ws_size;
    Int unnz = 0 ; 

    if(k!=0)
      {unnz = U.col_ptr(k);}

    Int uunnz = U.nnz;
   
    Int lcnt = 0;
    Int ucnt = 0; 
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
    printf("kid: %d col_start Unnz: %d %d \n", kid, unnz, uunnz);
    #endif


    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("---------Checking Workspace, kid: %d-------\n", 
	   kid);
    BASKER_ASSERT(top == ws_size, "top == ws_size");
    for(Int i = 0 ; i < ws_size; i++)
      {
	if(ws[i] != 0)
	  {
	    printf("--------------ERROR---------");
	    printf("kid: %d k: %d i: %d ws[i]=%d L.scol: %d \n",
		   kid, k, i, ws[i], L.scol);
	    ASSERT(ws[i] == 0);
	  }
      }
    #endif

    Int k_offset = 0;
    if(l != 0)
      {
	//printf("sep flag true, offset: %d \n",
	//     B.scol);
	k_offset = k;
      }

   

    //Might be better to go through backwards
    if(Options.same_pattern == BASKER_FALSE)
      {
    for(Int i = B.col_ptr(k-k_offset); 
	i < B.col_ptr(k-k_offset+1); ++i)
      {
	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid>=0)
	printf("kid: %d index: %d %d\n", 
	       kid, i, B.row_idx(i));
	#endif

	j = B.row_idx(i);
	if((l != 0)&&
	   (INC_LVL_TEMP(j+brow)==BASKER_MAX_IDX))
	  {
	    continue;
	  }
	

	X(j) = B.val(i);

	//INC_LVL_TEMP(j+brow) = 0;

	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid==0)
        printf("kid: %d i: %d  val: %g  top: %d \n", 
	       kid,i, B.val(i), top);
	if(kid==0)
	  printf("kid: %d Nx in Ak %d %g color = %d \n",
	    kid, j, X[j],  ws[j] );
        #endif

	//Note we need to check all
	//LVL = min max{}+1
	//Need to search all "possible" paths
	if(Options.incomplete_type == 
	   BASKER_INCOMPLETE_LVL)
	  {
	   t_local_reach_inc_lvl(kid, l, l, j, &top);
	  }
	else
	  {
	    if(gperm(j+brow) != BASKER_MAX_IDX)
	      {
		t_local_reach_inc_rlvl(kid,l,l,j,top);
	      }
	    else
	      {
		t_local_reach_short_inc_rlvl(kid,l,l,j,top);
	      }
	  }


      }//end over each nnz in column
      }//end same_pattern?
    else
      {
	
	//Since A may be over filled, 
	//WAnt to look at pattern first
	
	/*
	//if(kid == 0)
	  {
	    printf("Kid: %d add pattern \n",
		   kid);
	    printf("Ucol : %d %d kid: %d \n",
		   U.col_ptr(k+1), 
		   U.col_ptr(k),
		   kid);
	  }
	*/

	//Get U pattern
	for(Int i = U.col_ptr(k+1)-1; i >=  U.col_ptr(k);
	    i--)
	  {
	    j = U.row_idx(i);
	    // printf("Add U pattern k: %d row: %d kid: %d \n",
	    //	   k, j, kid);
	    color[j] = 2;
	    
	    if(Options.no_pivot == BASKER_TRUE)
	      {
		pattern[--top] = j;
	      }
	    else
	      {
		BASKER_ASSERT(0==1,
			      "Currently not supported");
		pattern[--top] = j;
	      }

	  }//end -- over all elements of U


	
	//Now populate X with new values of A
	 for(Int i = B.col_ptr(k-k_offset); 
	     i < B.col_ptr(k-k_offset+1); ++i)
	   {

	     j = B.row_idx(i);
	     //printf("Add X pattern k: %d row: %d val: %f kid: %d \n",
	     //	    k, j, B.val(i), kid);
	     
	     //printf("color size: %d \n",
	     //	    color.extent(0));
     
	     if(color[j] != 2)
	       {
		 continue;
	       }
	      
	     X(j) = B.val(i);
	     
	   }//end populate from B
	
      }//end if same pattern
    xnnz = ws_size - top;
   
    //return 0;
 
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid>=0)
    printf("xnnz: %d ws_size: %d top: %d , kid: %d\n",
	   xnnz, ws_size, top, kid);
    #endif


    //WE SHOUD DO A UNNZ COUNT
    //count number of nnz
    ucnt = 0;
    for(Int i=top; i < ws_size; ++i)
      {
	 j = pattern[i];
	 t = gperm(j+brow);
	 //printf("Checking pattern(%d)=%d gperm(%d) %d kid:%d \n",
	 //	i, j, j+brow, t, kid);
	 if(INC_LVL_TEMP(j+brow) <= Options.inc_lvl)
	   {
	     if(t == BASKER_MAX_IDX)
	       {
		 lcnt++;
	       }
	     else
	       {
		 ucnt++;
	       }
	   }
       }
     //Note: This +1 causes some trouble
     //ucnt = ws_size - top - lcnt +1;
     
     #ifdef BASKER_DEBUG_NFACTOR_COL
     if(kid>=0)
       printf("lcnt: %d ucnt: %d , kid: %d \n", lcnt, ucnt, kid);
     #endif

     if(unnz+ucnt-1 > uunnz)
       {
	 
	 if(Options.verbose == BASKER_TRUE)
	   {
//	 printf("kid: %ld col: %ld need to realloc, unnz: %ld ucnt: %ld uunnz: %ld U_col: %ld U_row: %ld \n", kid, k, unnz, ucnt, uunnz, U_col, U_row);
   std::cout << "kid: "  << kid
             << " col: " << k
             << " need to realloc, unnz: " << unnz
             << " ucnt: " << ucnt
             << " uunnz: " << uunnz
             << " U_col: " << U_col
             << " U_row: " << U_row << std::endl;
	   }
	 BASKER_ASSERT(0==1, "USIZE\n");
	 

	 Int newsize = (unnz+U.nrow) * 1.2  ;
	 
	 if(Options.realloc == BASKER_FALSE)
	   {
	     thread_array(kid).error_type =
	       BASKER_ERROR_NOMALLOC;
	     return BASKER_ERROR;
	   }
	 else
	   {
	     //printf("HERE\n");
	     thread_array(kid).error_type =
	       BASKER_ERROR_REMALLOC;
	     thread_array(kid).error_blk    = U_col;
	     thread_array(kid).error_subblk = U_row;
	     thread_array(kid).error_info   = newsize;
	     return BASKER_ERROR;
	   }//if/else realloc
       }//if need to realloc


     Entry relax_value = 0;
     if(false)
     //if(Options.incomplete_type == 
     // BASKER_INCOMPLETE_LVL)
       {
	 t_back_solve_inc_lvl(kid, l,l, k, top, xnnz);
       }
     else
       {
	 t_back_solve_inc_rlvl(kid,l,l,k,top,xnnz,
			       relax_value);
       }

     //move over nnz to U 
     for(Int i = top; i < ws_size; ++i)
       {
	 j = pattern[i];
	 pattern[i] = 0;
	 t = gperm(j+brow);
	 
	 #ifdef BASKER_DEBUG_NFACTOR_COL
	 //if(kid>=0)
	   {
	     printf("HERE considering j: %d t:%d lvl %d val: %e, kid: %d \n",
		    j, t, INC_LVL_TEMP(j+brow), X[j], kid);
	   }
	  #endif

	 //old zero checek
	   if((Options.same_pattern == BASKER_TRUE) ||
	      (INC_LVL_TEMP(j+brow) <= Options.inc_lvl))
	  {	    
	    //Note, if we remove this test, 
	    //we might get this to unroll!!!
	    if(t != BASKER_MAX_IDX)
              {
                #ifdef BASKER_DEBUG_NFACTOR_COL
		if(kid==0)
		{
		printf("U add kid: %d adding x[%d %d] %g to U inc_lvl: %d\n", kid, k+U.scol, j+U.srow, X(j), INC_LVL_TEMP(j+brow));
		}
                #endif

 		U.row_idx(unnz) = t-brow;
		U.val(unnz) = X(j);
		if(Options.same_pattern == BASKER_FALSE)
		  {
		    U.inc_lvl(unnz) = INC_LVL_TEMP(j+brow);
		    INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
		  }
                unnz++;
		X(j) = 0;

              }//if in U
            else
              {
		//printf("----ERROR--- KID: %d EXTRA L[%d]=%f \n",
		//     kid, j, X[j-brow]);
		BASKER_ASSERT(0==1, " "); 
              }//LOWER
	  }//END NOT 0
	 else
	   {
	     //Relax stuff here 
	     X(j) = 0;
	     if(Options.same_pattern == BASKER_FALSE)
	       {
		 INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
	       }
	   }
       }//OVER ALL X

     U.col_ptr(k+1) = unnz;
  
     #ifndef BASKER_MULTIPLE_UPPER
     //This is if we want the thread to update offdig.
     //The default will be to us multple upper
     //----------------------UPDATE OFFDIAG-----------------//
     x_row++;
     l_row++;
     for(; x_row < LL_size(x_col); x_row++, l_row++)
       {
	 
	 /*
	 t_back_solve_offdiag_inc_lvl(kid,
			      L_col, L_row,
			      X_col, X_row,
			      k, col_idx_offset,
			      U.val,
			      U.row_idx,
		     U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			      U.col_ptr[k-bcol],
			      A_option);
	 */
       }//end for over all offdiag
     #endif
     

     //Bad fix later
     if(Options.same_pattern == BASKER_FALSE)
       {
     if(l>0)
       {
	 for(Int si = 0; si < U.nrow; si++)
	   {
	     stack[si] = BASKER_MAX_IDX;
	   }
       }
       }

     return 0;
  }//end t_upper_col_factor_inc_lvl()


  //symbolic(fill) spmv
  template <class Int, class Entry, class Exe_Space>
  void 
  Basker<Int,Entry,Exe_Space>::t_upper_col_ffactor_offdiag2_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k, 
   const BASKER_BOOL lower
   )
  {

    const Int my_leader = 
      (sl==0)?(kid):find_leader_inc_lvl(kid,sl-1);
    
    if(kid != my_leader)
      {
	/*
	if(lower == BASKER_TRUE)
	  {
	printf("offd, kid: %d my_leader: %d \n",
	       kid, my_leader);
	  }
	*/
	return;
      }
    
    //if(kid==24)
      {
	// printf("ffactor kid: %d lvl: %d sl: %d l: %d \n",
	// kid, lvl, sl, l);
      }

    const Int L_col     = S(sl)(my_leader);
    Int L_row           = l-sl+1; //Might have to think about th
    const Int U_col     = S(lvl)(kid);
    
    Int my_row_leader = find_leader(kid,lvl-1);
    Int my_new_row = 
      L_col - S(0)(my_row_leader);
    // Int U_row = my_new_row;
    
  
    Int U_row     =
      (lvl==1)?(kid%2):S(sl)(kid)%LU_size(U_col);
    if((S(sl)(kid) > 14) &&
       (S(sl)(kid) > LU_size(U_col)) &&
       (lvl != 1))
      {
	//printf("lower offdiag new num, %d %d \n",
	//     S(sl)(kid), LU_size(U_col));
	Int tm = (S(sl)(kid)+1)/16;
	U_row = ((S(sl)(kid)+1) - (tm*16))%LU_size(U_col);
      }
  
    //printf("UFF kid:%d U: %d %d new: %d leader: %d %d lvl: %d l: %d sl: %d \n",
    //	   kid, U_col, U_row, my_new_row,
    //	   my_row_leader, L_col, 
    //	   lvl, l, sl);
    
    //JDB PASS TEST
    U_row = my_new_row;

    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
  
    //Int col_idx_offset  = 0;
    
    BASKER_MATRIX     &U   = LU(U_col)(U_row);
    //const Int         bcol = U.scol;
         
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    if(L_row >= LL_size(L_col))
      {
	printf("assert error. Lrow:%d L_col: %d size: %d kid: %d\n",
	       L_row, L_col, LL_size(L_col), kid);
      }
    BASKER_ASSERT(L_row < LL_size(L_col), "upper-off, Lrow >= size");
    BASKER_ASSERT(X_row < LL_size(X_col), "upper-off, Xrow >=size"); 
    #endif
    

    #ifdef BASKER_DEBUG_NFACTOR_COL2
    //if(lower == BASKER_TRUE)
      {
    printf("Upper_ffact_offdiag, kid: %d leader: %d l: %d lvl: %d works_size: %d X: %d %d L: %d %d U: %d %d k: %d \n",
	   kid, my_leader, l, lvl, LL_size[X_col], X_col, X_row, L_col, L_row, U_col, U_row,  k+U.scol);
      }
     #endif

      #ifdef BASKER_DEBUG_NFACTOR_COL2
    //if(lower == BASKER_TRUE)
      {
    printf("OFF-DIAG, kid: %d, l: %d  X: %d %d L: %d %d \n",
	   kid, l,X_col, X_row, L_col, L_row);
      }
      #endif

      //const BASKER_BOOL A_option = BASKER_FALSE;


    //U.info();
    //JDB TEST
    //printf("kid: %d k: %d \n",
    //	   kid, k);
    //printf("kid: %d k: %d testpoint: %d %d \n",
    //	   kid, k, U.col_ptr(k+1), U.col_ptr(k));

      //     if(Options.incomplete_type == 
      // BASKER_INCOMPLETE_LVL)
      if(false)
	{
	  t_lower_col_offdiag_find_fill(kid,
				  L_col, L_row,
				  X_col, X_row,
				  k,
				  U.val, //can be removed
				  U.row_idx,
				  U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k));
	}
      else
	{
	   t_lower_col_offdiag_find_fill_rlvl(kid,
				  L_col, L_row,
				  X_col, X_row,
				  k,
				  U.val, //can be removed
				  U.row_idx,
				  U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k));
	}
    
    //if lower, finish off the updates
    if(lower == BASKER_TRUE)
      {
	X_row++;
	L_row++;
	for(; X_row < LL_size(X_col);++X_row, ++L_row)
	  {
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("LLL OFF-DIAG,kid:%d, l: %d X: %d %d L: %d %d U: %d %d \n",
		   kid, l,X_col, X_row, L_col, L_row, U_col, U_row);
	    #endif

	 
	    // if(Options.incomplete_type ==
	    // BASKER_INCOMPLETE_LVL)
	    if(false)
	      {
	    t_lower_col_offdiag_find_fill(kid,
					  L_col, L_row,
					  X_col, X_row,
					  k,
					  U.val,
					  U.row_idx,
					  U.inc_lvl,
				   U.col_ptr(k+1)-U.col_ptr(k),
					  U.col_ptr(k));
	      }
	    else
	      {

		 t_lower_col_offdiag_find_fill_rlvl(kid,
					  L_col, L_row,
					  X_col, X_row,
					  k,
					  U.val,
					  U.row_idx,
					  U.inc_lvl,
				   U.col_ptr(k+1)-U.col_ptr(k),
					  U.col_ptr(k));

	      }

	  }//for --over remainder blks
      }//if --lower  
  }//end t_upper_col_ffactor_offdiag2_inc_lvl()

  
  template <class Int, class Entry, class Exe_Space>
  void
  Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag2_same_pattern_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k, 
   const BASKER_BOOL lower
   )
  {

    const Int my_leader = 
      (sl==0)?(kid):find_leader_inc_lvl(kid,sl-1);
    
    if(kid != my_leader)
      {
	return;
      }
    
    const Int L_col     = S(sl)(my_leader);
    Int L_row           = l-sl+1; //Might have to think about th
    const Int U_col     = S(lvl)(kid);

    Int my_row_leader = find_leader(kid,lvl-1);
    Int my_new_row = 
      L_col - S(0)(my_row_leader);
    Int U_row  = 0;
    U_row = my_new_row;


    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
  
    Int col_idx_offset  = 0;
    
    BASKER_MATRIX     &U   = LU(U_col)(U_row);
  

    //Need to give them the output pattern
    Int U_pattern_col = S(lvl)(kid);
    Int my_pattern_leader = 
      find_leader_inc_lvl(kid,l);
    Int U_pattern_row = S(l+1)(my_pattern_leader) - 
      S(0)(my_row_leader);

    /*
    printf("Test mypleader: %d myrowleader: %d kid: %d\n", 
	   my_pattern_leader, my_row_leader, kid);
    printf("Test one: %d two: %d kid: %d \n",
	   S(l+1)(my_pattern_leader),
	   S(0)(my_row_leader),
	   kid);
    */
   
    
    Int L_pattern_col = S(lvl)(kid);
    Int L_pattern_row = BASKER_MAX_IDX;
    if(lower == BASKER_TRUE)
      {
	L_pattern_row = 0;
      }
         
    /*
    printf("DEBUG U Pattern: %d %d kid: %d \n", 
	   U_pattern_col, U_pattern_row, kid);
    printf("DEBUG L Pattern: %d %d kid: %d \n",
	   L_pattern_col, L_pattern_row, kid);
    printf("DEBUG L: %d %d kid: %d \n",
	   L_col, L_row, kid);
    */

    t_same_pattern_back_solve_offdiag_inc_lvl(kid,
				      L_col, L_row,
				      X_col, X_row,
				U_pattern_col, U_pattern_row,
				L_pattern_col, L_pattern_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				       U.inc_lvl,
				 U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),

				        BASKER_FALSE);
					      


      /*
    t_dense_back_solve_offdiag_inc_lvl(kid,
			       L_col, L_row,
			       X_col, X_row,
			       k, col_idx_offset,
			       U.val,
			       U.row_idx,
			       U.inc_lvl,       
			       U.col_ptr(k+1)-U.col_ptr(k),
			       U.col_ptr(k),
			       BASKER_FALSE);
      */
    
    //if lower, finish off the updates
    if(lower == BASKER_TRUE)
      {
	X_row++;
	L_row++;
	
	U_pattern_row = BASKER_MAX_IDX;
	L_pattern_row = 1;

	for(; X_row < LL_size(X_col);
	    ++X_row, ++L_row, ++L_pattern_row)
	  {
	   
	    
	    
    t_same_pattern_back_solve_offdiag_inc_lvl(kid,
				      L_col, L_row,
				      X_col, X_row,
				U_pattern_col, U_pattern_row,
				L_pattern_col, L_pattern_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				       U.inc_lvl,
				 U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),
				        BASKER_FALSE);
 



	    
	    /*
	    t_dense_back_solve_offdiag_inc_lvl(kid,
				       L_col, L_row,
				       X_col, X_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				       U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),
				       BASKER_FALSE
				       );
	    */

	  }//for --over remainder blks
      }//if --lower
  
  }//end t_upper_col_factor_offdiag2_same_pattern_inc_lvl()


  //uses local idx for local blks X
  template <class Int, class Entry, class Exe_Space>
  void 
  Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag2_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k, 
   const BASKER_BOOL lower
   )
  {

    const Int my_leader = 
      (sl==0)?(kid):find_leader_inc_lvl(kid,sl-1);
    
    if(kid != my_leader)
      {
	/*
	if(lower == BASKER_TRUE)
	  {
	printf("offd, kid: %d my_leader: %d \n",
	       kid, my_leader);
	  }
	*/
	return;
      }
    
    const Int L_col     = S(sl)(my_leader);
    Int L_row           = l-sl+1; //Might have to think about th
    const Int U_col     = S(lvl)(kid);



    Int my_row_leader = find_leader(kid,lvl-1);
    Int my_new_row = 
      L_col - S(0)(my_row_leader);
    // Int U_row = my_new_row;
   

    Int U_row     =
      (lvl==1)?(kid%2):S(sl)(kid)%LU_size(U_col);
    if((S(sl)(kid) > 14) &&
       (S(sl)(kid) > LU_size(U_col)) &&
       (lvl != 1))
      {
	Int tm = (S(sl)(kid)+1)/16;
	U_row = ((S(sl)(kid)+1) - (tm*16))%LU_size(U_col);
      }

    // printf("lowerspmv kid: %d U: %d %d new %d leader: %d %d lvl: %d %d %d \n",
    //	   kid, U_col, U_row, my_new_row, my_row_leader, 
    //	   lvl, l, sl);

    //JDB PASSED TEST
    U_row = my_new_row;


    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
  
    Int col_idx_offset  = 0;
    
    BASKER_MATRIX     &U   = LU(U_col)(U_row);
    //const Int         bcol = U.scol;
         
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    if(L_row >= LL_size(L_col))
      {
	printf("assert error. Lrow:%d L_col: %d size: %d kid: %d\n",
	       L_row, L_col, LL_size(L_col), kid);
      }
    BASKER_ASSERT(L_row < LL_size(L_col), "upper-off, Lrow >= size");
    BASKER_ASSERT(X_row < LL_size(X_col), "upper-off, Xrow >=size"); 
    #endif
    

    #ifdef BASKER_DEBUG_NFACTOR_COL2
    //if(lower == BASKER_TRUE)
      {
    printf("Upper_fact_offdiag, kid: %d leader: %d l: %d lvl: %d works_size: %d X: %d %d L: %d %d U: %d %d k: %d \n",
	   kid, my_leader, l, lvl, LL_size[X_col], X_col, X_row, L_col, L_row, U_col, U_row,  k+U.scol);
      }
     #endif

     #ifdef BASKER_DEBUG_NFACTOR_COL2
    //if(lower == BASKER_TRUE)
      {
    printf("OFF-DIAG, kid: %d, l: %d  X: %d %d L: %d %d \n",
	   kid, l,X_col, X_row, L_col, L_row);
      }
      #endif

      //const BASKER_BOOL A_option = BASKER_FALSE;
	    
    t_dense_back_solve_offdiag_inc_lvl(kid,
			       L_col, L_row,
			       X_col, X_row,
			       k, col_idx_offset,
			       U.val,
			       U.row_idx,
			       U.inc_lvl,       
			       U.col_ptr(k+1)-U.col_ptr(k),
			       U.col_ptr(k),
			       BASKER_FALSE);
    
    //if lower, finish off the updates
    if(lower == BASKER_TRUE)
      {
	X_row++;
	L_row++;
	for(; X_row < LL_size(X_col);++X_row, ++L_row)
	  {
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("LLL OFF-DIAG,kid:%d, l: %d X: %d %d L: %d %d U: %d %d \n",
		   kid, l,X_col, X_row, L_col, L_row, U_col, U_row);
	    #endif

	    t_dense_back_solve_offdiag_inc_lvl(kid,
				       L_col, L_row,
				       X_col, X_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				       U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),
				       BASKER_FALSE
				       );

	  }//for --over remainder blks
      }//if --lower
  
  }//end t_upper_col_factor_offdiag2_inc_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void
  Basker<Int,Entry,Exe_Space>::t_same_pattern_update_matrix_inc_lvl
  (
   const Int kid,
   const Int team_leader,
   const Int lvl,
   const Int l, 
   const Int k
   )
  {
    
    const Int leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;
//    Int gbrow = 0; //NDE - warning: unused

    //Over each blk    
//    Int last_blk = l+2; //NDE - warning: unused
    
    {
      //Copy B -> C
      Int bl = l+1;
      Int A_col = S(lvl)(kid);
      
      Int my_row_leader = find_leader(kid,lvl-1);
      Int my_new_row = 
	S(bl)(kid) - S(0)(my_row_leader);
      Int A_row = 0;
      A_row = my_new_row;

      
      BASKER_MATRIX  *Bp;
      if(A_row != (LU_size(A_col)-1))
	{
	  //printf("upper picked, kid: %d \n", kid);
	  //printf("up: %d %d kid: %d \n",
	  //	   A_col, A_row, kid);
	  Bp = &(AVM(A_col)(A_row));
	}
      else
	{
	  //printf("lower picked, kid: %d\n", kid);
	  Bp = &(ALM(A_col)(0));
	}  
      BASKER_MATRIX   &B  = *Bp;
      //printf("ADDING UPDATES TO B\n");
      //B.info();
      //B.print();
//      gbrow = B.srow; //NDE - warning: unused
      
      //return;
      
      //Int team_leader   = find_leader(kid, l);  //Not used
      ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
      INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
      Int *color = &(ws(0));
      LL(leader_idx)(bl).p_size = 0;

      //Get the columns pattern
      Int U_pattern_col = A_col;
      Int U_pattern_row = A_row;
      Int L_pattern_col = A_col;
      Int L_pattern_row = BASKER_MAX_IDX;
      if(A_row == (LU_size(A_col)-1))
	{
	  L_pattern_row = 0;
	}

        
      //Copy in X
      for(Int i = B.col_ptr(k);
	  i < B.col_ptr(k+1); ++i)
	{
	  Int B_row = B.row_idx(i);	
	  X(B_row) += B.val(i);
	}

     
      //Copy into C
      BASKER_MATRIX &Up = LU(U_pattern_col)(U_pattern_row);
      for(Int i = Up.col_ptr(k); i < Up.col_ptr(k+1); i++)
	{
	  const Int j = Up.row_idx(i);
	  C.row_idx(nnz) = j;
	  C.val(nnz)     = X(j);
	  nnz++;
	  X(j)            = 0;
	  color[j] = 0;
	}

      //if there is a L
      if(L_pattern_row != BASKER_MAX_IDX)
	{
	  BASKER_MATRIX &Lp = LL(L_pattern_col)(L_pattern_row);
	  for(Int i = Lp.col_ptr(k)+1; i < Lp.col_ptr(k+1);i++)
	    {
	      const Int j = Lp.row_idx(i);
	      C.row_idx(nnz) = j;
	      C.val(nnz)     = X(j);
	      nnz++;
	      X(j)            = 0;
	      color[j] = 0;
	    }//loop over L pattern
	}
	    
    }//DUMMY (Do we need this?)

    C.col_ptr(0) = 0;
    C.col_ptr(1) = nnz; 

  }//end t_update_matrix_inc_lvl

  
  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  void 
  Basker<Int,Entry,Exe_Space>::t_dense_copy_update_matrix2_inc_lvl
  (
   const Int kid,
   const Int team_leader,
   const Int lvl,
   const Int l, 
   const Int k
   )
  {
  
    const Int leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;
    Int gbrow = 0;

    //Over each blk    
//    Int last_blk = l+2; //NDE - warning: unused
    /*
    if(lvl ==(l+1))
      {
	last_blk = LL_size(leader_idx);
      }
    */
    
    //do we need this temp loop anymore ??
    for(Int TEMP=0; TEMP < 1; TEMP++)
      {
    //Copy B -> C
    Int bl = l+1;
    Int A_col = S(lvl)(kid);

    Int my_row_leader = find_leader(kid,lvl-1);
    Int my_new_row = 
      S(bl)(kid) - S(0)(my_row_leader);
    //Int A_row = my_new_row;
    

    Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
    if((S(bl)(kid) > 14) &&
       (S(bl)(kid) > LU_size(A_col)) &&
       (lvl != 1))
      {
	//printf("test cm %d %d %d \n",
	//     kid, S(bl)(kid), LU_size(A_col));

	Int tm = (S(bl)(kid)+1)/16;
	A_row  = ((S(bl)(kid)+1) - (tm*16))%LU_size(A_col);
      } 
     

  
    //printf("Dense TEST kid: %d A: %d %d new: %d leader: %d %d lvl: %d %d \n", 
    //	   kid, A_col, A_row, my_new_row, my_row_leader, S(bl)(kid), lvl, bl);

    //JDB PASS
    A_row = my_new_row;

    //Int CM_idx = kid; //Not used
    
    BASKER_MATRIX  *Bp;
    if(A_row != (LU_size(A_col)-1))
      {
	//printf("upper picked, kid: %d \n", kid);
	//printf("up: %d %d kid: %d \n",
	//	   A_col, A_row, kid);
	Bp = &(AVM(A_col)(A_row));
      }
    else
      {
	//printf("lower picked, kid: %d\n", kid);
	Bp = &(ALM(A_col)(0));
      }  
    BASKER_MATRIX   &B  = *Bp;
    //printf("ADDING UPDATES TO B\n");
    //B.info();
    //B.print();
    gbrow = B.srow;

    //return;

    //Int team_leader   = find_leader(kid, l);  //Not used
    ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
    //const Int brow    = LL(leader_idx)(bl).srow;
    //const Int nrow    = LL(leader_idx)(bl).nrow;
    //Int p_size        = LL(leader_idx)(bl).p_size;
    //const Int ws_size = LL(leader_idx)(bl).iws_size;
    //Int *color        = &(ws(0));
    //Int *pattern      = &(color[ws_size]);
    //Int *stack        = &(pattern[ws_size]); //used for fill

    #ifdef BASKER_DEBUG_NFACTOR_COL2    
    printf("copy, kid: %d bl: %d  A: %d %d \n", 
	   kid, bl, A_col, A_row);
    #endif

    //const Int bbcol = B.scol;
   
    for(Int i = B.col_ptr(k);
	i < B.col_ptr(k+1); ++i)
      {
	Int B_row = B.row_idx(i);	
	//Int j = gperm(B_row+B.srow); //Not used

     	X(B_row) += B.val(i);
	
	//printf("Added B(%d) %g %g \n",
	//     B_row, B.val(i), X(B_row));
	

      }//end for over all nnz
      }//DUMMY (Do we need this?)
        
    //-------------move into C------------------- 
    //(Right now have to do dense but like to do sparse)

//    last_blk = LL_size(leader_idx); //NDE - warning: unused
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    
    //Do I need this ?
    for(Int TEMP = 0; TEMP < 1; ++TEMP)
      {
	Int bl = l+1;
	//Int A_col = S(lvl)(kid);
	
	/* CAN REMOVE
	Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
	//maybe no???
	if((S(bl)(kid) > 14) &&
	   (S(bl)(kid) > LU_size(A_col)) &&
	   (lvl != 1))
      {
	//printf("test cm %d %d %d \n",
	//     kid, S(bl)(kid), LU_size(A_col));

	Int tm = (S(bl)(kid)+1)/16;
	A_row  = ((S(bl)(kid)+1) - (tm*16))%LU_size(A_col);
      } 
	*/
	
     // printf("kid: %d leader_idx: %d bl: %d \n",
     //	    kid, leader_idx, bl);


        //Int CM_idx = kid;
    ENTRY_1DARRAY   X   = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws   = LL(leader_idx)(bl).iws;
    const Int   ws_size = LL(leader_idx)(bl).ews_size;
//    const Int      brow = LL(leader_idx)(bl).srow; //NU //NDE - warning: unused
    const Int      nrow = LL(leader_idx)(bl).nrow;
    //Int p_size          = LL(leader_idx)(bl).p_size;

    //For recounting patterns in dense blk
    //Need better sparse update
    //Int p_count  =0 ; 
    
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);
    
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    //printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
    //	   kid, A_col, A_row, team_leader, bl,p_size);
    #endif
      
    //over all dim(S)
    //for(Int jj=brow; jj < (brow+nrow); jj++)
    for(Int jj=0; jj < nrow; ++jj)
      {
	//Int j = pattern[jj];
	Int j = jj;
        #ifdef BASKER_DEBUG_NFACTOR_COL22
	printf("considering: %d %d %f, kid: %d\n",
	       j,brow,X[j], kid);
        #endif
	
	
	//if(X(j) != 0)
	if(X(j) != (Entry)(0))
	  {
	    if(bl == l+1)
	      {
                #ifdef BASKER_DEBUG_NFACTOR_COL22
		printf("moving X[%d] %f kid: %d nnz: %d Csize: %d \n",
		       j+brow, X(j), kid, nnz, C.nrow);
                #endif
		C.row_idx(nnz) = j;
		C.val(nnz)     = X(j);
		nnz++;
		X(j)            = 0;
		color[j] = 0;
		
		INC_LVL_TEMP(j+gbrow) = stack[j];
		//printf("To C: index: %d lvl: %d kid: %d \n",
		//j+gbrow, stack[j], kid);

		//Clear stack
		//stack[j] = BASKER_MAX_IDX;

	      }
	    else
	      {
		#ifdef BASKER_DEBUG_NFACTOR_COL22
		printf("counting [%d] %f kid: %d \n",
		       j+brow, X[j], kid);
                #endif
		BASKER_ASSERT("1==0", "copy, should never call");
	      }
	  }//if X(j) != 0
	//printf("kid: %d clearing stack[%d] \n",
	//     kid, j+brow);
	stack[j] = BASKER_MAX_IDX;
      }//for -- length of x

      }//DUMMY	
    C.col_ptr(0) = 0;
    C.col_ptr(1) = nnz;

    //printf("Done with move, kid: %d found nnz: %d \n",
    //   kid, nnz);
    
    //C.info();
    //C.print();

  }//end t_dense_copy_update_matrix2_inc_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::t_lower_col_factor_inc_lvl
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k, 
   Entry &opivot
   )
  {
    //Get needed variables
    const Int L_col = S(lvl)(kid);
    const Int L_row = 0;
    const Int U_col = S(lvl)(kid);
    const Int U_row = LU_size(U_col)-1;
    
    const Int X_col = S(0)(kid);
    //Int col_idx_offset = 0; //can we get rid of now?
    

    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid == 0)
      {
    printf("LOWER_COL_FACTOR kid: %d \n", kid);
    printf("kid %d using L: %d %d  U: %d %d  X %d \n",
	   kid, L_col, L_row, U_col, U_row, X_col);
      }
    #endif
    //end get needed variables

    BASKER_MATRIX        &L = LL(L_col)(L_row);
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
   
    BASKER_MATRIX        &B = thread_array(kid).C;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
      {
	printf("---------------LOWER FACTOR KID: %d ------\n",
	       kid);
	//B.base->info();
	//B.base->print();
	B.info();
	B.print();
	printf("After matrix print \n");
      }
    #endif
    /*
    if(kid == 0)
      {
	B.print();
      }
    */

    INT_1DARRAY  ws       = LL(X_col)(l+1).iws;
    const Int     ws_size = LL(X_col)(l+1).iws_size;
    ENTRY_1DARRAY X       = LL(X_col)(l+1).ews;

    const Int brow     = U.srow;
    //const Int bcol     = U.scol;

    const Int lval  = L.col_ptr(k);
    const Int uval  = U.col_ptr(k);
    
    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);
    Int *stack     = &(pattern[ws_size]);

    Int i,j;
//    Int top, top1, maxindex, t; //NDE - warning: top1 set but unused
    Int top, maxindex, t;
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;

    Int newsize;
    Entry pivot, value;
    Entry absv, maxv;
    Int llnnz = L.nnz;
    Int uunnz = U.nnz;
    
    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
//    top1 = ws_size; //NDE - warning: top1 set but unused
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("-------------lower_col_factor-------\n");
    printf("JB ecol: %d L.nnz: %d U.nnz: %d brow: %d ws_size: %d  kid: %d\n",
	   U.scol+U.ncol, llnnz, uunnz, brow, ws_size, kid);
    #endif

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("\n----------------- K = %d --------------\n", 
	   k+U.scol);
    #endif
    	  
    value = 0.0;
    pivot = 0.0;

    maxindex = BASKER_MAX_IDX;
    lcnt = 0;
    ucnt = 0;
    
   #ifdef BASKER_DEBUG_NFACTOR_COL
   ASSERT(top == ws_size);
   for(i = 0 ; i < ws_size; i++){
     if(X[i] !=0)
       {
	 printf("--Error, kid: %d X[%d] = %f \n", kid, i,X[i]); 
       }
     ASSERT(X[i] == 0);}
   for(i = 0; i <  ws_size; i++){ASSERT(ws[i] == 0 );}
   #endif


   #ifdef BASKER_DEBUG_NFACTOR_COL
   printf("--------------PAST ASSERT----------- \n");
   printf("B.col(k): %d  B.col(k+1): %d \n", 
	  B.col_ptr[0], B.col_ptr[1]);
   #endif

   if(Options.same_pattern == BASKER_FALSE)
     {
   for(i = B.col_ptr(0); i < B.col_ptr(1); ++i)
     {
             
       j = B.row_idx(i);
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid==0)
       printf("j: %d i: %d \n", j, i);
       #endif
              
       if(INC_LVL_TEMP(j+brow) == BASKER_MAX_IDX)
	 {
	   continue;
	 }

       X(j) = B.val(i);
        
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid==0)
	 {
	   printf("i: %ld  j: %ld %ld  val: %g  top: %d \n", 
		  i, j, gperm(j+brow), B.val(i), top);
	 }
       if(kid==0)
       printf("Nxk in Ak %d %g color = %d  given inc: %d\n",
	      j+brow, X[j],  
	      ws[j ],
	      INC_LVL_TEMP(j+brow));
       #endif
       
       if(Options.incomplete_type == 
	  BASKER_INCOMPLETE_LVL)
	 {
	   t_local_reach_inc_lvl(kid,lvl,l+1,j,&top);
	 }
       else
	 {
	   //printf("lower reach 2 called \n");
	   if(gperm(j+brow) != BASKER_MAX_IDX)
	     {
	       t_local_reach_inc_rlvl(kid,lvl,l+1,j,top);
	     }
	   else
	     {
	       t_local_reach_short_inc_rlvl(kid,lvl,l+1,j,top);
	     }
	 }
       
     }//over each nnz in the column
     }//if--same pattern
   else
     {
       //Get L
       for(i = L.col_ptr(k+1)-1; i >= L.col_ptr(k);
	   i--)
	 {
	   j = L.row_idx(i);
	   color[j] = 2;
	   pattern[--top] = j;
	 }//end get L

       //Get U pattern
       for(i = U.col_ptr(k+1)-2; i >= U.col_ptr(k);
	   i--)
	 {
	   j = U.row_idx(i);
	   color[j] = 2;	   
	   //Pivot options here
	   pattern[--top] = j;
	 }// end get U
     
       //Fill in the needed X
        for(i = B.col_ptr(0); i < B.col_ptr(1); ++i)
	  {
	    j = B.row_idx(i);
	    if(color[j] == 2)
	      {
		X(j) = B.val(i);
	      }
	  }//for-B
     }//if same pattern
   xnnz = ws_size - top;
   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid==0)
     {
       printf("==DEBUG===\n");
       printf("xnnz: %d ws_size: %d top: %d \n",
	      xnnz, ws_size, top);
       for(Int z = top; z < ws_size; z++)
	 {
	   printf("pattern[%d]: %d %f\n",
		  z, pattern[z], X(pattern[z]));
	 }
	   
     }
   #endif
  

   Entry relax_value = 0;
   // if(Options.incomplete_type ==
   // BASKER_INCOMPLETE_LVL)
   if(false)
   {
     // printf("lower back solve 1 called \n");
       t_back_solve_inc_lvl(kid, lvl,l+1, k, top, xnnz); 
     }
   else
     {
       // printf("lower back solve 2 called \n");
       t_back_solve_inc_rlvl(kid,lvl,l+1,k,top,xnnz,
			     relax_value);
     }


   //find pivot
   maxv = 0.0;
   for(i = top; i < ws_size; i++)
     { 
       j = pattern[i];
       t = gperm(j+brow);
       
       //printf("considering: %d %f lvl: %d \n", 
       //     j, X(j), INC_LVL_TEMP(j+brow));

     
       if(INC_LVL_TEMP(j+brow) <= Options.inc_lvl)
	 {

	   value = X(j);
	   absv = EntryOP::approxABS(value);
	   //if(t == L.max_idx)
	   if(t == BASKER_MAX_IDX)
	     {
	       lcnt++;
	       if((Options.no_pivot != BASKER_TRUE) &&
		  EntryOP::gt(absv,maxv))
		 {
		   maxv     = absv;
		   pivot    = value;
		   maxindex = j;                
		 }
	     }//if lower-half
	   else
	     {
	       ucnt++;
	     }//if upper-half
	 }//if lvl <options  
   }//for (i = top; i < ws_size)
            
   //ucnt = ws_size - top - lcnt +1;
  
   //----------------------------Sym-----
   //SYM
   if(Options.no_pivot == BASKER_TRUE)
     {
       maxindex = k;
       pivot = X(k);
     }

   U.tpivot = pivot;
  
   //printf("lower pivot: %f k: %d kid: %d \n",
   //	  U.tpivot, k, kid);

   //local only
   //opivot = pivot;

   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid>=0)
   printf("pivot found: %f , kid: %d \n", pivot, kid);
   #endif

   //if((maxindex == L.max_idx) || (pivot == 0)
   if((maxindex == BASKER_MAX_IDX) || (pivot == (Entry)(0)) )
     {
       cout << "Error: Matrix is singular, col, lvl: " << l <<endl;
       cout << "MaxIndex: " << maxindex << " pivot " 
	    << pivot << endl;
       return 2;
     }          
  
   gperm(maxindex+brow) = k+brow; 
   gpermi(k+brow) = maxindex+brow;
   
   //printf("set pivot: %d %d \n",
   //	  maxindex+brow, 
   //	  k+brow);

   if(lnnz + lcnt > llnnz)
     {
       //Note: comeback
       newsize = lnnz * 1.1 + 2 *L.nrow + 1;
     
       if(Options.verbose == BASKER_TRUE)
	 {
       cout << "Lower Col Reallocing L oldsize: " 
	    << llnnz 
	    << " newsize: " << newsize << endl;
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
	   thread_array(kid).error_blk    = L_col;
	   thread_array(kid).error_subblk = -1;
	   thread_array(kid).error_info   = newsize;
	   return BASKER_ERROR;
	 }
     }
   if(unnz+ucnt > uunnz)
     {
       //Note: comeback
       newsize = uunnz*1.1 + 2*U.nrow+1;

       if(Options.verbose == BASKER_TRUE)
	 {
       cout << "Lower Col Reallocing U oldsize: " 
	    << uunnz 
	    << " newsize " << newsize << endl;
	 }

       if(Options.realloc == BASKER_FALSE)
	 {
	   thread_array(kid).error_type = 
	     BASKER_ERROR_NOMALLOC;
	 }
       else
	 {
	   thread_array(kid).error_type = 
	     BASKER_ERROR_REMALLOC;
	   thread_array(kid).error_blk    = U_col;
	   thread_array(kid).error_subblk = U_row;
	   thread_array(kid).error_info   = newsize;
	   return BASKER_ERROR;
	 }
     }
  
   L.row_idx(lnnz) = maxindex;
   L.val(lnnz)     = (Entry) 1.0;
   L.inc_lvl(lnnz) = INC_LVL_TEMP(maxindex+brow);
   INC_LVL_TEMP(maxindex+brow) = BASKER_MAX_IDX;
   lnnz++;

   #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
   printf("add L(%d) 1.0 lvl: %d \n",
	  maxindex+brow, L.inc_lvl(lnnz));
   #endif

   Entry lastU = (Entry) 0.0;
  
   //For every nnz in LHS
   for( i = top; i < ws_size; i++)
     {
       j = pattern[i];
       pattern[i] = 0;
       t = gperm(j+brow);
       
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(k>=0)
	 printf("j %d  t %d val: %g , kid: %d \n", j, t, kid, 
		X(j));
       #endif
       
       
       if((Options.same_pattern == BASKER_TRUE) ||
	  (INC_LVL_TEMP(j+brow) <= Options.inc_lvl))
	 {
           #ifdef BASKER_DEBUG_NFACTOR_COL
	   #ifdef BASKER_2DL
	   if(kid>=0)
	     {
	       //printf("found value: %f at %d, kid: %d \n",
	       //     X[j-brow], j, kid);
	       printf("found value: %f at %d, kid: %d \n",
	            X[j], j, kid);
	     }
	   #else
	   if(kid>=0)
	     printf("found value: %f at %d, kid: %d \n",
		    X[j], j, kid);
	   #endif
           #endif
	   
	   if(t != BASKER_MAX_IDX)
	     {
	       if(t < (k+brow))
		 {
		 
		   U.row_idx(unnz) = t-brow;
		   U.val(unnz) = X(j);
		   unnz++;
		   
		   #ifdef BASKER_DEBUG_NFACTOR_COL_INC
		   printf("add U(%d) [%d %d]: %g lvl: %d \n",
			  U.row_idx(unnz-1),
			  k+brow,
			  U.row_idx(unnz-1)+brow,
			  U.val(unnz-1),
			  INC_LVL_TEMP(j+brow));
		   #endif

		 }
	       else
		 {
		   
		   //printf("setting last: %g \n", X(j));
		   lastU = X(j);
		 }
	     }
	   else if (t == BASKER_MAX_IDX)
	     {
               #ifdef BASKER_DEBUG_NFACTOR_COL
	       if(kid>=0)
		 {
		   printf("inserting %f at %d into %d \n", 
			  X[j]/pivot, j, lnnz );
		 }
               #endif

	       L.row_idx(lnnz) = j;
	       L.val(lnnz) = EntryOP::divide(X(j), pivot);
	       
               #ifdef BASKER_DEBUG_NFACTOR_BLK_INC
	       printf("add L(%d) [%d %d]: %g lvl %d \n",
		      j+brow,
		      k+L.scol, 
		      j+brow,
		      L.val(lnnz), INC_LVL_TEMP(j+brow));
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
	     //Add relaxation here !
	      INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
	   }
	 //move inside if not zero..... extra set ops not needed
	 //INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
	 X(j) = 0;
       
     }//if(x[j-brow] != 0)
   
   //Fill in last element of U
   U.row_idx(unnz) = k;
   U.val(unnz) = lastU;
   unnz++;
   
   //printf("Add U(%d) %g \n",
   //	  k, U.val(unnz-1));

   xnnz = 0;
   top = ws_size;
   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   printf("setting col: %d %d %d %d\n",
	  k, cu_ltop, k+1-bcol, lnnz);
   #endif
   L.col_ptr(k) = cu_ltop;
   L.col_ptr(k+1) = lnnz;
   cu_ltop = lnnz;
   

   U.col_ptr(k) = cu_utop;
   U.col_ptr(k+1) = unnz;
   cu_utop = unnz;

   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid>=0)
   printf("col_fact k: %d Unnz: %d   Lnnz: %d \n",
	  k, unnz, lnnz);
   #endif
  
   
   #ifndef BASKER_MULTIPLE_LOWER
   //DEFAULT IS MULTIPLE LOWER
   //COME BACK AND DELETE
   BASKER_ASSERT(0==1, "MULTIPLE LOWER\N");
   #ifdef BASKER_2DL
   //----Update offdiag (this can be done in parallel l>0)---//
   //-Need to fix
   #ifdef BASKER_DEBUG_NFACTOR_COL
   printf("k: %d -----TTTEST(blk_row-- %d %d \n",
     kid, lvl+1, LL_size[L_col]);
   printf("k: %d -----TTTEST(x_row)--- %d %d \n",
       kid, l+2, LL_size[X_col]);
   #endif
   for(Int blk_row = L_row+1, x_row = l+2; 
       blk_row < LL_size(L_col); blk_row++, x_row++)

     { 
       t_dense_back_solve_offdiag(kid,
			    L_col, blk_row,
			    X_col, x_row,
			    k, col_idx_offset,
			    U.val, U.row_idx,
			    U.col_ptr(k+1)-U.col_ptr(k),
			    U.col_ptr(k),
			    BASKER_TRUE);

       t_dense_move_offdiag_L(kid, 
			L_col, blk_row,
			X_col, x_row,
			k, pivot);
     }//end for over all offdiag blks
   #endif
   #endif

 

   //FIX LATER (these will now be reset by lower)
   
   if(Options.same_pattern == BASKER_FALSE)
     {
   for(Int si = 0; si < L.nrow; si++)
     {
       stack[si] = BASKER_MAX_IDX;
     }
     }
   

   #ifdef BASKER_DEBUG_NFACTOR_DEBUG
   print_factor(L,U);
   #endif
   
   return 0;
  }//end t_lower_col_fact_inc_lvl()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void 
  Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag2_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int l,
   const Int k,
   Entry pivot
   )
  {
  
    #ifndef BASKER_MULTIPLE_LOWER
    BASKER_ASSERT(0==1,"BASKER_MULTIPLE_LOWER ERROR");
    return 0;
    #endif


    const Int leader_id   = find_leader(kid, l);
    const Int lteam_size  = pow(2,l+1);
    const Int L_col       = S(lvl)(leader_id);
    Int L_row             = 0;
    const Int U_col       = S(lvl)(leader_id);
    Int U_row             = LU_size(U_col)-1;
    Int X_col             = S(0)(leader_id);
    Int X_row             = l+1;
    Int col_idx_offset    = 0;  //can get rid of?
   
    //BASKER_MATRIX        &L = LL(L_col)(L_row); //NDE - warning: unused L
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    
    INT_1DARRAY     ws = LL(X_col)(X_row).iws;
    //const Int  ws_size = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY    X = LL(X_col)(X_row).ews;

    //const Int brow     = U.srow;
    //const Int bcol     = U.scol;

    pivot = U.tpivot;
        
    Int U_pattern_row = BASKER_MAX_IDX;
    Int U_pattern_col = L_col;
    //Int L_pattern_row = L_row; //NDE - warning: unused 
    //Int L_pattern_col = L_col; //NDE - warning: unused 


    //printf("OFF_DIAG_LOWER, kid: %d leaderid: %d t_size: %d \n"
    //   kid, leader_id, lteam_size);
    
    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; 
	 L_row < LL_size(L_col);
	 X_row+=(lteam_size), L_row+=(lteam_size))
    
     { 
       
       //Using leader_id is a bad hack!!!!

       if(Options.same_pattern == BASKER_TRUE)
	 {

	       t_same_pattern_back_solve_offdiag_inc_lvl(leader_id,
				      L_col, L_row,
				      X_col, X_row,
				U_pattern_col, U_pattern_row,
				L_col, L_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				       U.inc_lvl,
				 U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),
				        BASKER_TRUE);
 
	       //printf("Calling t_move_offdiag_L");
	       t_move_offdiag_L_inc_lvl(leader_id, 
			      L_col, L_row,
			      X_col, X_row,
			      k, pivot);


	  

	 }
       else
	 {

       t_dense_back_solve_offdiag_inc_lvl(leader_id,
				  L_col, L_row,
				  X_col, X_row,
				  k, col_idx_offset,
				  U.val, U.row_idx,
				  U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k),
				  BASKER_TRUE);


       
       t_dense_move_offdiag_L_inc_lvl(leader_id, 
			      L_col, L_row,
			      X_col, X_row,
			      k, pivot);

	 }
       
     }//end for over all offdiag blks
    
    //Now we can clear the inc-array
    //Race here
    /*
    for(Int i = 0 ; i < U.nrow; i++)
      {
	INC_LVL_TEMP(i+U.srow) = BASKER_MAX_IDX;
      }
    */

  }//end t_lower_col_factor_offdiag2_inc_lvl()


   template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void 
  Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag2_cleanup_inc_lvl
  (
   const Int kid,
   const Int lvl,
   const Int l,
   const Int k
   )
  {
  
    #ifndef BASKER_MULTIPLE_LOWER
    BASKER_ASSERT(0==1,"BASKER_MULTIPLE_LOWER ERROR");
    return 0;
    #endif


    const Int leader_id   = find_leader(kid, l);
    //const Int lteam_size  = pow(2,l+1); //NDE - warning: unused
//    const Int L_col       = S(lvl)(leader_id); //NDE - warning: unused 
//    Int L_row             = 0; //NDE - warning: unused 
    const Int U_col       = S(lvl)(leader_id);
    Int U_row             = LU_size(U_col)-1;
    Int X_col             = S(0)(leader_id);
    Int X_row             = l+1;
    //Int col_idx_offset    = 0;  //can get rid of?//NDE - warning: unused 
   
    //BASKER_MATRIX        &L = LL(L_col)(L_row); //NDE - warning: unused
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    
    INT_1DARRAY     ws = LL(X_col)(X_row).iws;
    //const Int  ws_size = LL(X_col)(X_row).iws_size; 
    ENTRY_1DARRAY    X = LL(X_col)(X_row).ews;

    if(kid == leader_id)
      {
    for(Int i = 0 ; i < U.nrow; i++)
      {
	INC_LVL_TEMP(i+U.srow) = BASKER_MAX_IDX;
      }
      }
  

  }//end t_lower_col_factor_offdiag2_cleanup_inc_lvl()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_same_pattern_col_copy_inc_lvl
   (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k,
   const BASKER_BOOL lower
   )
  {
    const Int my_idx     = S(0)(kid);

    //should remove either as a paramter or here
    Int team_leader      = find_leader(kid, sl);
    const Int leader_idx = S(0)(team_leader);
      
    //If I an not a leader, then need to copy over
    if(kid != team_leader)
      {

	Int endblk = (lower)?(LL_size(my_idx)):(l+2);
	for(Int blk = l+1; blk < endblk; ++blk)
	  {

	    //const Int blk = l+1;  
	    ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews;
//	    INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws; //NDE - warning: unused
//	    Int      p_sizeL  = LL(leader_idx)(blk).p_size; //NDE - warning: unused
//	    Int      ws_sizeL = LL(leader_idx)(blk).iws_size; //NDE - warning: unused
	    ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
	    INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
//	    const Int ws_size = LL(my_idx)(blk).iws_size; //NDE - warning: unused
	    //Int       p_size  = LL(my_idx)(blk).p_size;
	    LL(my_idx)(blk).p_size = 0;
	    Int       *color  = &(ws[0]);
//	    Int     *pattern  = &(color[ws_size]);  //NDE - warning: unused 
//	    Int     *stack    = &(pattern[ws_size]); //NDE - warning: unused
	    //used for fill
//	    Int      brow     = LL(my_idx)(blk).srow; //NU //NDE - warning: unused
	    //Int      browL    = LL(leader_idx)(blk).srow; //NU
	    
//	    Int *colorL   = &(wsL(0)); //NDE - warning: unused
	    //This may be in error in the bigger code
//	    Int *patternL = &(colorL[ws_sizeL]); //NDE - warning: unused 
//	    Int *stackL   = &(patternL[ws_sizeL]); //NDE - warning: unused
	    


	    //What to get pattern from U/L
//	    Int my_pattern_leader = find_leader(kid,l); //NDE - warning: unused
	    /*
	    printf("Copy, my pattern leader: %d kid: %d %d \n",
		   my_pattern_leader, find_leader(kid,l), 
		   kid);

	    printf("Options: %d %d kid: %d \n",
		   S(0)(0), S(l+1)(my_pattern_leader),kid);

	    */

	    /*
	    printf("NUMBER: kid: %d lvl: %d sl: %d l: %d col_leader: %d %d my_blk: %d %d U_col: %d \n",
		   kid, lvl, sl, l, 
		   find_leader(kid,lvl-1),S(0)(find_leader(kid,lvl-1)),
		   kid, S(l+1)(kid), 
		   S(l+1)(kid) - S(0)(find_leader(kid,lvl-1)));
	    */



	    
	    Int U_pattern_col = S(lvl)(kid);
	    Int U_pattern_row = BASKER_MAX_IDX;
	    
	    if(blk == l+1)
	      {
		//U_pattern_row = S(0)(my_pattern_leader) -
		//S(0)(find_leader(kid,lvl));
		//U_pattern_row = S(l+1)(kid) - 
		//S(0)(my_pattern_leader);
		U_pattern_row = S(l+1)(kid) - 
		  S(0)(find_leader(kid,lvl-1));
	      }

	    Int L_pattern_col = S(lvl)(kid);
	    Int L_pattern_row = BASKER_MAX_IDX;
	    if(lower == BASKER_TRUE)
	      {
		L_pattern_row = blk-(l+1);
	      }

	    /*
	    printf("Debug copy, Up: %d %d kid: %d\n",
		   U_pattern_col, U_pattern_row, 
		   kid);

	    printf("Debug copy, Lp: %d %d kid: %d \n",
		   L_pattern_col, L_pattern_row, 
		   kid);
	    */
	   
	    //copy over
	    if(U_pattern_row != BASKER_MAX_IDX)
	      {

		BASKER_MATRIX &UP = LU(U_pattern_col)(U_pattern_row);
	    
		for(Int jj = UP.col_ptr(k);
		    jj < UP.col_ptr(k+1);
		    jj++)
		  {
		    const Int jjj = UP.row_idx(jj);
		    color[jjj] = 0; 
		    XL(jjj) += X(jjj);
		    X(jjj)   = 0;
		  }
	      }//if UPattern
	    if(L_pattern_row != BASKER_MAX_IDX)
	      {
		BASKER_MATRIX &LP = LL(L_pattern_col)(L_pattern_row);
		for(Int jj = LP.col_ptr(k);
		    jj < LP.col_ptr(k+1);
		    jj++)
		  {
		    const Int jjj = LP.row_idx(jj);
		    color[jjj] = 0;
		    XL(jjj) += X(jjj);
		    X(jjj) = 0;
		  }
	      }//if LPattern
	  }//over blks
      }//not leader
  }//end t_same_pattern_col_copy_inc_lvl()

 
 template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry, Exe_Space>::t_dense_blk_col_copy_atomic2_inc_lvl
  (
   const Int kid,
   const Int NOT_USED,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k,
   const BASKER_BOOL lower
   )
  {
    //Setup
    //Int A_col = S(lvl)(kid);
    //Int A_row = (lvl==1)?(2):S(l+1)(kid)%(LU_size(A_col));
   
    //printf("DEBUG, kid: %d k: %d A_col: %d A_row: %d \n", 
    //	   kid, k, A_col, A_row);


    //BASKER_MATRIX &B    = AVM(A_col)(A_col);

    const Int my_idx     = S(0)(kid);

    //should remove either as a paramter or here
    Int team_leader      = find_leader(kid, sl);
    const Int leader_idx = S(0)(team_leader);
    //Int loop_col_idx     = S(l)(kid); NU

    //#ifdef BASKER_DEBUG_NFACTOR_COL2
    if(lower == BASKER_TRUE)
      {
	
	#ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("Called t_blk_col_copy_atomic kid: %d \n " , kid);
    printf("Copying col, kid: %d  k: %d lvl: %d l: %d \n", 
	   kid,k,lvl, l);
    #endif
    //printf("Copying Col, kid: %d k:%d  A: %d %d to tl: %d li: %d\n",
    //	   kid, k, A_col, A_row, team_leader, leader_idx);
    
      }
    //#endif
   
    //If I an not a leader, then need to copy over
    if(kid != team_leader)
      {

	Int endblk = (lower)?(LL_size(my_idx)):(l+2);
	//Int endblk = l+2;

	/*
	printf("l+2: %d endblk: %d \n",
	       l+2, endblk);
	*/

	for(Int blk = l+1; blk < endblk; ++blk)
	  {

	    //const Int blk = l+1;  
	ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews;
//	INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws; //NDE - warning: unused
	Int      p_sizeL  = LL(leader_idx)(blk).p_size;
//	Int      ws_sizeL = LL(leader_idx)(blk).iws_size; //NDE - warning: unused
	ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
	INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
	const Int ws_size = LL(my_idx)(blk).iws_size;
	//Int       p_size  = LL(my_idx)(blk).p_size; 
	Int       *color  = &(ws[0]);
	Int     *pattern  = &(color[ws_size]); 
	Int     *stack    = &(pattern[ws_size]); //used for fill
//	Int      brow     = LL(my_idx)(blk).srow; //NU //NDE - warning: unused
	//Int      browL    = LL(leader_idx)(blk).srow; //NU

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_TRUE)
	  {
	printf("kid: %d  COPY INDEX %d %d to %d %d \n",
	       kid, my_idx, blk, leader_idx, blk);
	  }
	#endif

//	Int *colorL   = &(wsL(0)); //NDE - warning: unused
	//This may be in error in the bigger code
//	Int *patternL = &(colorL[ws_sizeL]); //NDE - warning: unused 
//	Int *stackL   = &(patternL[ws_sizeL]); //NDE - warning: unused

        #ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_TRUE)
	  {
	printf("t_b_col_copy, kid: %d wsize: %d \n", 
	       kid, ws_size);
	printf("t_b_col_copy,kid:%d ps: %d XL: %d %d X: %d %d\n",
	       kid, p_size, leader_idx, blk, my_idx, blk);
	  }
	#endif

	//over all nnnz found
	for(Int jj = 0; jj < LL(my_idx)(blk).nrow; ++jj)
	  {
	    color[jj] = 0;
		
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    if(lower == BASKER_TRUE)
	      {
		printf("kid: %d jj: %d nrows: %d xsize: %d %d \n",
		       kid, jj, LL(my_idx)(blk).nrow, 
		       LL(my_idx)(blk).iws_size,
		       LL(leader_idx)(blk).iws_size);

	    printf("Adding X(%d)%f to XL(%d) %f, leader: %d kid: %d\n",
		  jj+brow, X[jj], jj, XL[jj], team_leader, kid);

	      }
	    #endif
	
//	    if(X(jj) !=0)
	    if(X(jj) != (Entry)(0))
	      {
		#ifdef BASKER_DEBUG_NFACTOR_COL2
		//if(lower == BASKER_TRUE)
		  {
		printf("Atomic Adding X(%d) %f to XL(%d) %f, kid: %d\n",
			       jj+brow, X[jj], jj, XL[jj], kid);
		  }
		#endif

		XL(jj) += X(jj);
		X(jj)   = 0;
		
                #ifdef BASKER_DEBUG_NFACTOR_COL_INC
		printf("Copy over mfill: %d lfill:%d \n",
		       stack[jj], stackL[jj]);
		#endif
	
		//if((stack[j] != BASKER_MAX_IDX)
	   
		//stackL[jj] = min(stackL[jj], stack[jj]);
		//Clear my stack
		
		    		    			
		
	      }//if X(j) != 0
	    //clear my stack
	    stack[jj] = BASKER_MAX_IDX;
	  }//for - over all nnz
	   
	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_TRUE)
	  {
	printf("----DEBUG--- X %d %d brow: %d kid: %d ws_size: %d\n",
	       my_idx, blk, brow, kid, ws_size);
	for(Int j = 0; j< ws_size; ++j)
	  {
	    printf("X[%d,%d] = %f , ptern: %d color: %d kid: %d \n",
		   j, j+brow, X[j], pattern[j], color[j],kid);
	  }
	  }
	 #endif

	//This can be removed in the future
	if(kid != team_leader)
	  {
	    LL(my_idx)(blk).p_size = 0;
	  }
	else
	  {
           #ifdef BASKER_DEBUG_NFACTOR_COL22
	    printf("SETTING PS: %d L:%d %d kid: %d\n",
		   p_sizeL, leader_idx, blk, kid);
	    #endif
	    LL(leader_idx)(blk).p_size = p_sizeL;
	    //p_size = 0; NOT USED
	  }//over all blks
	  }
      }//if not team_leader
   
  }//end t_dense_blk_col_copy_atomic2_inc_lvl()

}//end namespace BaskerNS

#endif //end def BASKER_NFACTOR_COL_INC_HPP
