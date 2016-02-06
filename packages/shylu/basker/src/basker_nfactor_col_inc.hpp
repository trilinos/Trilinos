#ifndef BASKER_NFACTOR_COL_INC_HPP
#define BASKER_NFACTOR_COL_INC_HPP

#include "basker_types.hpp"
#include "basker_stats.hpp"
#include "basker_thread.hpp"

#include "basker_nfactor_blk_inc.hpp"

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

      //if(kid==12 || kid==13 || kid==14 || kid==15)
      //if(kid == 2 || kid ==3 )
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
	  
    const Int scol = LU(U_col)(U_row).scol;
    const Int ecol = LU(U_col)(U_row).ecol;

    
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n  LVL=%d  ----- kid: %d --\n\n",
	   lvl, kid);
    #endif
    
    //Do all domains (old sublevel 0)
    //If this works well we can move this for into the function


    //for(Int k = 0; k < 2; ++k)
     for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("UPPER, kid: %d k: %d \n",
	       kid, k);
	#endif
	
	t_upper_col_factor_inc_lvl(kid, team_leader, 
				   lvl, 0, 
				   k,
				   BASKER_FALSE);


      }//over all columns / domains / old sublevel 0

    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n\n done with UPPER, kid: %d \n\n\n", kid);
    #endif

   
  
    //------Need because extend does not 
    //-------------Barrier--Between Domains-------------
    Int my_leader = find_leader(kid, 0);
    Int b_size    = pow(2,1);
    //barrier k = 0 usedl1
    t_basker_barrier_inc_lvl(thread,kid,my_leader,
		     b_size, 0, LU(U_col)(U_row).scol, 0);
   
    
    
    //----------------Sep level upper tri-------------
    for(Int l = 1; l < (lvl); ++l)
      {
	
	//for(Int k = 2; k < 3; ++k)
	for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
	  {
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("\n\nSep, upper update, kid: %d k=%d \n\n",
		   kid, k+LU(U_col)(U_row).scol);
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
		
		t_upper_col_factor_inc_lvl(kid, team_leader, 
				   lvl, l, 
				   k,
				   BASKER_FALSE);
		

	      }//if correct kid to do this sublevels upperfactor
	    
	  }//over all columns
      }//for - over all sublevel 1...lvl-2


    
   
    //---------Lower Factor (old sublevel lvl-1)-------
    
    my_leader = find_leader_inc_lvl(kid, lvl-1);
    b_size    = pow(2,lvl);
    // printf("[3] barrier test, kid: %d leader: %d b_size: %d lvl: %d \n",
    //	   kid,  my_leader, b_size, lvl);
    t_basker_barrier_inc_lvl(thread, kid, my_leader,
		     b_size, 7, LU(U_col)(U_row).scol, 0);

    //printf("\n\n======= LOWER, KID: %d ======= \n\n", kid);
    
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
        
        //if((kid == 0)||(kid==1))
	t_lower_col_factor_offdiag2_inc_lvl(kid, lvl, lvl-1, k, pivot);
	//thread.team_barrier();
	my_leader = find_leader_inc_lvl(kid, lvl-1);
	b_size    = pow(2,lvl);
	//printf("barrier test 5, leader: %d b_size: %d lvl: %d \n",
	//     my_leader, b_size, lvl);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
		     b_size, 9, k, lvl-1);
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


        //printf("========FILL IN KID: %d =======\n", kid);
	 t_add_orig_fill(kid, lvl, l, k, lower); 
	 //printf("=======AFFERT FILL_IN KID: %d ====\n", kid);
	 //printf("\n\n");
	
	//Barrier--Start
	my_leader = find_leader_inc_lvl(kid,sl);
	b_size    = pow(2,sl+1);
	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 1, k+k_offset, sl);
	//Barrier--End

        //Upward scan
	//Copy to leader
	if(kid%((Int)pow(2,sl))==0)
	  {
	    //printf("=========REDUCE KID: %d \n", kid);
	    t_reduce_col_fill(kid,
			      lvl, sl, l, k, lower);
	  }
	
	//printf("Done with sreduce. kid: %d \n", kid);

	//Barrier--Start
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
	t_upper_col_factor_offdiag2_inc_lvl(kid, lvl, sl,l, k, lower);
	
	//Barrier--Start
	my_leader = find_leader_inc_lvl(kid,sl);
	b_size    = pow(2,sl+1);
	//printf("[1] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 5, k+k_offset, sl);
	//Barrier--End

	if(kid%((Int)pow(2,sl))==0)
	  {
	    t_dense_blk_col_copy_atomic2_inc_lvl(kid, my_leader,
					 lvl, sl, l, k, lower);
	  }

	//Barrier--Start
	//printf("[2] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier_inc_lvl(thread, kid, my_leader,
			 b_size, 6, k+k_offset, sl);
	
      }//over all sublevels

    t_dense_copy_update_matrix2_inc_lvl(kid, my_leader, lvl, l, k);


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
    Int my_token = S[l][kid];
    Int my_loc = kid;
    while((my_loc > 0))
      {
	my_loc--;
	if(S[l][my_loc] != my_token)
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
  BASKER_INLINE
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
    //Get needed variables
    const Int L_col = S(l)(kid);
    const Int L_row = 0;
    const Int U_col = S(lvl)(kid);
    Int U_row = (lvl==1)?(kid%2):S(l)(kid)%LU_size(U_col);
  
    if((L_col > 14) &&
       (L_col > LU_size(U_col)) &&
       (lvl != 1))
      {
	//printf("modify urow, %d %d \n",
	//     L_col, LU_size(U_col));
	
	Int tm = (L_col+1)/16;
	U_row = ((L_col+1)-(tm*16))%LU_size(U_col);

      }


    const Int X_col = S(0)(kid);
    const Int X_row = l; //X_row = lower(L)
    const Int col_idx_offset = 0; //we might be able to remove
  
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
    printf("kid %d, upper using L: %d %d  U: %d %d  X %d %d\n",
	   kid, L_col, L_row, U_col, U_row, X_col, X_row);
    #endif
    //end get needed variables//

    BASKER_MATRIX        &L = LL(L_col)(L_row);
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    
    //Ask C++ guru if this is ok
    BASKER_MATRIX        *Bp;
    if(l == 0 )
      {
        Bp = &(AVM(U_col)(U_row));
      }
    else
      {
	Bp = &(thread_array[kid].C);
      }
    BASKER_MATRIX    &B = *Bp;
    // if(kid ==2)
    // {
    //	printf("2====KID === 2\n");
    // B.print();
    //  }

    INT_1DARRAY ws     = LL(X_col)(X_row).iws;
    const Int ws_size  = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY X    = LL(X_col)(X_row).ews;
  
    const Int brow = U.srow;
    const Int bcol = U.scol;

    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);
    
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
    for(Int i = B.col_ptr(k-k_offset); 
	i < B.col_ptr(k-k_offset+1); ++i)
      {
	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid>=0)
	printf("kid: %d index: %d %d\n", 
	       kid, i, B.row_idx(i));
	#endif

	j = B.row_idx(i);
	X(j) = B.val(i);

	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid>=0)
        printf("kid: %d i: %d  val: %g  top: %d \n", 
	       kid,i, B.val(i), top);
	if(kid>=0)
	  printf("kid: %d Nx in Ak %d %g color = %d \n",
	    kid, j, X[j],  ws[j] );
        #endif

	//Note we need to check all
	//LVL = min max{}+1
	//Need to search all "possible" paths
	t_local_reach_inc_lvl(kid, l, l, j, &top);

      }//end over each nnz in column
    xnnz = ws_size - top;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid>=0)
    printf("xnnz: %d ws_size: %d top: %d , kid: %d\n",
	   xnnz, ws_size, top, kid);
    #endif


    //WE SHOUD DO A UNNZ COUNT
    //count number of nnz
     for(Int i=top; i < ws_size; ++i)
       {
	 j = pattern[i];
	 t = gperm(j+brow);
	 if(t == BASKER_MAX_IDX)
	   {
	     lcnt++;
	   }
       }
     //Note: This +1 causes some trouble
     ucnt = ws_size - top - lcnt +1;
     
     #ifdef BASKER_DEBUG_NFACTOR_COL
     if(kid>=0)
       printf("lcnt: %d ucnt: %d , kid: %d \n", lcnt, ucnt, kid);
     #endif

     if(unnz+ucnt-1 > uunnz)
       {

	 printf("kid: %d col: %d need to realloc, unnz: %d ucnt: %d uunnz: %d U_col: %d U_row: %d \n", kid, k, unnz, ucnt, uunnz, U_col, U_row);
	 BASKER_ASSERT(0==1, "USIZE\n");
       }

     t_back_solve_inc_lvl(kid, l,l, k, top, xnnz);
     
     //move over nnz to U 
     for(Int i = top; i < ws_size; ++i)
       {
	 j = pattern[i];
	 pattern[i] = 0;
	 t = gperm(j+brow);
	 
	 #ifdef BASKER_DEBUG_NFACTOR_COL
	 if(kid>=0)
	   {
	     printf("considering j: %d t:%d val: %e, kid: %d \n",
		  j, t, X[j], kid);
	   }
         #endif

	 //old zero checek
	  {	    
	    //Note, if we remove this test, 
	    //we might get this to unroll!!!
	    if(t != BASKER_MAX_IDX)
              {
                #ifdef BASKER_DEBUG_NFACTOR_COL
		if(kid==2)
		{
		printf("U add kid: %d adding x[%d %d] %g to U inc_lvl: %d\n", kid, k+U.scol, j+U.srow, X(j), INC_LVL_TEMP(j+brow));
		}
                #endif

 		U.row_idx(unnz) = t-brow;
		U.val(unnz) = X(j);
		U.inc_lvl(unnz) = INC_LVL_TEMP(j+brow);
		INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
                unnz++;
		X(j) = 0;

              }//if in U
            else
              {
		printf("----ERROR--- KID: %d EXTRA L[%d]=%f \n",
		       kid, j, X[j-brow]);
		BASKER_ASSERT(0==1, " "); 
              }//LOWER
	  }//END NOT 0
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
    
    const Int L_col     = S(sl)(my_leader);
    Int L_row           = l-sl+1; //Might have to think about th
    const Int U_col     = S(lvl)(kid);
    Int U_row     =
      (lvl==1)?(kid%2):S(sl)(kid)%LU_size(U_col);

    //printf("test \n");
    if((S(sl)(kid) > 14) &&
       (S(sl)(kid) > LU_size(U_col)) &&
       (lvl != 1))
      {
	//printf("lower offdiag new num, %d %d \n",
	//     S(sl)(kid), LU_size(U_col));
	Int tm = (S(sl)(kid)+1)/16;
	U_row = ((S(sl)(kid)+1) - (tm*16))%LU_size(U_col);
      }

    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
  
    Int col_idx_offset  = 0;
    
    BASKER_MATRIX     &U   = LU(U_col)(U_row);
    const Int         bcol = U.scol;
         
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

    const BASKER_BOOL A_option = BASKER_FALSE;
	    

    t_lower_col_offdiag_find_fill(kid,
				  L_col, L_row,
				  X_col, X_row,
				  k,
				  U.val, //can be removed
				  U.row_idx,
				  U.inc_lvl,
				  U.col_ptr(k+1)-U.col_ptr(k),
				  U.col_ptr(k));
    
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

	    
	    t_lower_col_offdiag_find_fill(kid,
					  L_col, L_row,
					  X_col, X_row,
					  k,
					  U.val,
					  U.row_idx,
					  U.inc_lvl,
				   U.col_ptr(k+1)-U.col_ptr(k),
					  U.col_ptr(k));

	  }//for --over remainder blks
      }//if --lower  
  }//end t_upper_col_ffactor_offdiag2_inc_lvl()


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
    Int U_row     =
      (lvl==1)?(kid%2):S(sl)(kid)%LU_size(U_col);


    //printf("test \n");
    if((S(sl)(kid) > 14) &&
       (S(sl)(kid) > LU_size(U_col)) &&
       (lvl != 1))
      {
	//printf("lower offdiag new num, %d %d \n",
	//     S(sl)(kid), LU_size(U_col));
	Int tm = (S(sl)(kid)+1)/16;
	U_row = ((S(sl)(kid)+1) - (tm*16))%LU_size(U_col);
      }

    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
  
    Int col_idx_offset  = 0;
    
    BASKER_MATRIX     &U   = LU(U_col)(U_row);
    const Int         bcol = U.scol;
         
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

    const BASKER_BOOL A_option = BASKER_FALSE;
	    
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
    Int last_blk = l+2;   
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

    Int CM_idx = kid;
    
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


    Int team_leader   = find_leader(kid, l);
    ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
    const Int brow    = LL(leader_idx)(bl).srow;
    const Int nrow    = LL(leader_idx)(bl).nrow;
    Int p_size        = LL(leader_idx)(bl).p_size;
    const Int ws_size = LL(leader_idx)(bl).iws_size;
    Int *color        = &(ws(0));
    Int *pattern      = &(color[ws_size]);
    Int *stack        = &(pattern[ws_size]); //used for fill

    #ifdef BASKER_DEBUG_NFACTOR_COL2    
    printf("copy, kid: %d bl: %d  A: %d %d \n", 
	   kid, bl, A_col, A_row);
    #endif

    const Int bbcol = B.scol;
   
    for(Int i = B.col_ptr(k);
	i < B.col_ptr(k+1); ++i)
      {
	Int B_row = B.row_idx(i);	
	Int j = gperm(B_row+B.srow);

     	X(B_row) += B.val(i);

	//printf("Added B(%d) %g fill-in: %d \n",
	//     B_row, X(B_row), stack[B_row]);
	

      }//end for over all nnz
      }//DUMMY (Do we need this?)
        
    //-------------move into C------------------- 
    //(Right now have to do dense but like to do sparse)

    last_blk = LL_size(leader_idx);
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    
    //Do I need this ?
    for(Int TEMP = 0; TEMP < 1; ++TEMP)
      {
	Int bl = l+1;
	Int A_col = S(lvl)(kid);
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

     // printf("kid: %d leader_idx: %d bl: %d \n",
     //	    kid, leader_idx, bl);


    Int CM_idx = kid;
    ENTRY_1DARRAY   X   = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws   = LL(leader_idx)(bl).iws;
    const Int   ws_size = LL(leader_idx)(bl).ews_size;
    const Int      brow = LL(leader_idx)(bl).srow;
    const Int      nrow = LL(leader_idx)(bl).nrow;
    Int p_size          = LL[leader_idx][bl].p_size;

    //For recounting patterns in dense blk
    //Need better sparse update
    Int p_count  =0 ; 
    
    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);
    Int *stack   = &(pattern[ws_size]);
    
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
	   kid, A_col, A_row, team_leader, bl,p_size);
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
	
	
	if(X(j) != 0)
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
		//     j+gbrow, stack[j], kid);

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
    Int col_idx_offset = 0; //can we get rid of now?
    

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("LOWER_COL_FACTOR kid: %d \n", kid);
    printf("kid %d using L: %d %d  U: %d %d  X %d \n",
	   kid, L_col, L_row, U_col, U_row, X_col);
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
    //B.print();


    INT_1DARRAY  ws       = LL(X_col)(l+1).iws;
    const Int     ws_size = LL(X_col)(l+1).iws_size;
    ENTRY_1DARRAY X       = LL(X_col)(l+1).ews;

    const Int brow     = U.srow;
    const Int bcol     = U.scol;

    const Int lval  = L.col_ptr(k);
    const Int uval  = U.col_ptr(k);
    
    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);

    Int i,j;
    Int top, top1, maxindex, t;
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
    top1 = ws_size;
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("-------------lower_col_factor-------\n");
    printf("JB scol: %d ecol: %d L.nnz: %d U.nnz: %d brow: %d ws_size: %d  kid: %d\n",
	   bcol, U.scol+U.ncol, llnnz, uunnz, brow, ws_size, kid);
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

   for(i = B.col_ptr(0); i < B.col_ptr(1); ++i)
     {
             
       j = B.row_idx(i);
     
       
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid>=0)
       printf("j: %d i: %d \n", j, i);
       #endif

       X(j) = B.val(i);
        
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid>=0)
	 {
	   printf("i: %ld  j: %ld %ld  val: %g  top: %d \n", 
		  i, j, gperm(j+brow), B.val(i), top);
	 }
       if(kid>=0)
       printf("Nxk in Ak %d %g color = %d \n",
	      j, X[j],  
	      ws[j ] );
       #endif

       t_local_reach_inc_lvl(kid,lvl,l+1,j,&top);
       
     }//over each nnz in the column
   xnnz = ws_size - top;
   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid>=0)
     {
       printf("xnnz: %d ws_size: %d top: %d \n",
	      xnnz, ws_size, top);
     }
   #endif
  
   t_back_solve_inc_lvl(kid, lvl,l+1, k, top, xnnz); // note: l not lvl given

   //find pivot
   maxv = 0.0;
   for(i = top; i < ws_size; i++)
     { 
       j = pattern[i];
       t = gperm(j+brow);

       value = X(j);
       
       //printf("considering: %d %f \n", 
       //      j, value);

       //absv = abs(value);
       absv = EntryOP::approxABS(value);
       //if(t == L.max_idx)
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
   ucnt = ws_size - top - lcnt +1;
   

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
   if((maxindex == BASKER_MAX_IDX) || (pivot == 0))
     {
       cout << "Error: Matrix is singular, col, lvl: " << l <<endl;
       cout << "MaxIndex: " << maxindex << " pivot " 
	    << pivot << endl;
       return 2;
     }          
  
   //gperm[maxindex] = k;
   gperm(maxindex+brow) = k+brow; 
   //gpermi[k] = maxindex;
   gpermi(k+brow) = maxindex+brow;
   
   if(lnnz + lcnt > llnnz)
     {
       //Note: comeback
       newsize = lnnz * 1.1 + 2 *A.nrow + 1;
       cout << "Lower Col Reallocing L oldsize: " << llnnz 
	    << " newsize: " << newsize << endl;
     }
   if(unnz+ucnt > uunnz)
     {
       //Note: comeback
       newsize = uunnz*1.1 + 2*A.nrow+1;
       cout << "Lower Col Reallocing U oldsize: " << uunnz 
	    << " newsize " << newsize << endl;
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
	 printf("j %d  t %d , kid: %d \n", j, t, kid);
       #endif
       
       
       //if(X(j) != 0)
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
		   #ifdef BASKER_DEBUG_NFACTOR_COL
		   #ifdef BASKER_2DL
		   if(kid>=0)
		   printf("U insert: %f at %d \n",
			  X[j], t);
		   #else
		   if(kid>=0)
		   printf("U insert: %f at %d \n",
			  X[j], gperm[j]);
		   #endif
		   #endif

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
	       printf("add L(%d): %g lvl %d \n",
		      j, L.val(lnnz), INC_LVL_TEMP(j+brow));
	       #endif

	       L.inc_lvl(lnnz) = INC_LVL_TEMP(j+brow);
	       INC_LVL_TEMP(j+brow) = BASKER_MAX_IDX;
	       lnnz++;
	     }
	 }//end if() not zero             
       //move inside if not zero..... extra set ops not needed
  
       X(j) = 0;
       
     }//if(x[j-brow] != 0)
   
   //Fill in last element of U
   U.row_idx(unnz) = k;
   U.val(unnz) = lastU;
   unnz++;
   
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
   
    BASKER_MATRIX        &L = LL(L_col)(L_row);
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    
    INT_1DARRAY     ws = LL(X_col)(X_row).iws;
    const Int  ws_size = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY    X = LL(X_col)(X_row).ews;

    const Int brow     = U.srow;
    const Int bcol     = U.scol;

    pivot = U.tpivot;
        
    //printf("OFF_DIAG_LOWER, kid: %d leaderid: %d t_size: %d \n",
    //	   kid, leader_id, lteam_size);
    
    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; 
	 L_row < LL_size(L_col);
	 X_row+=(lteam_size), L_row+=(lteam_size))
    
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
       
     }//end for over all offdiag blks

  }//end t_lower_col_factor_offdiag2_inc_lvl()



 template <class Int, class Entry, class Exe_Space>
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
    Int loop_col_idx     = S(l)(kid);

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
	INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws;
	Int      p_sizeL  = LL(leader_idx)(blk).p_size;
	Int      ws_sizeL = LL(leader_idx)(blk).iws_size;
	ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
	INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
	const Int ws_size = LL(my_idx)(blk).iws_size;
	Int       p_size  = LL(my_idx)(blk).p_size;
	Int       *color  = &(ws[0]);
	Int     *pattern  = &(color[ws_size]); 
	Int     *stack    = &(pattern[ws_size]); //used for fill
	Int      brow     = LL(my_idx)(blk).srow;
	Int      browL    = LL(leader_idx)(blk).srow;

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_TRUE)
	  {
	printf("kid: %d  COPY INDEX %d %d to %d %d \n",
	       kid, my_idx, blk, leader_idx, blk);
	  }
	#endif

	Int *colorL   = &(wsL(0));
	//This may be in error in the bigger code
	Int *patternL = &(colorL[ws_sizeL]);
	Int *stackL   = &(patternL[ws_sizeL]);

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
	
	    if(X(jj) !=0)
	      {
		#ifdef BASKER_DEBUG_NFACTOR_COL2
		if(lower == BASKER_TRUE)
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
	
		stackL[jj] = min(stackL[jj], stack[jj]);
		    		    			
		
	      }//if X(j) != 0
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
	    p_size = 0;
	  }//over all blks
	  }
      }//if not team_leader
   
  }//end t_dense_blk_col_copy_atomic2_inc_lvl()

}//end namespace BaskerNS

#endif //end def BASKER_NFACTOR_COL_INC_HPP
