#ifndef SHYLUBASKER_NFACTOR_COL2_HPP
#define SHYLUBASKER_NFACTOR_COL2_HPP

//#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_thread.hpp"

#include "shylubasker_nfactor_blk.hpp"
#include "shylubasker_nfactor_col.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

#ifdef BASKER_DEBUG
#include <assert.h>
#endif

#ifdef HAVE_VTUNE
#include "ittnotify.h"
#endif

//#define BASKER_DEBUG_NFACTOR_COL2
//#define BASKER_DEBUG_TIME
//#define BASKER_COUNT_OPS

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_sep2
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif
    
    Basker<Int,Entry,Exe_Space> *basker;
    Int lvl;
    
    kokkos_nfactor_sep2()
    {}
    
    kokkos_nfactor_sep2(
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
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //	      thread.team_rank());
      Int kid = basker->t_get_kid(thread);

      //team leader is not using
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #else
      Int team_leader = 0; //Note: come back and fix
      #endif


      #ifdef HAVE_VTUNE
      __itt_pause();
      #endif

      //if(kid==12 || kid==13 || kid==14 || kid==15)
      //if(kid ==8 || kid ==9 || kid ==10 || kid ==11)
      //if(kid==4 || kid ==5 || kid==6 || kid==7)
	{
      #ifdef BASKER_KOKKOS
       basker->t_nfactor_sep2(kid, lvl, team_leader, thread);
      #else
      
      #endif
	}

	#ifdef HAVE_VTUNE
	__itt_resume();
	#endif
            
    }//end operator ()
  };//end col_factor_funct


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_nfactor_sep2
  (
   const Int kid,
   const Int lvl,
   const Int team_leader,
   const TeamMember &thread
   )
  {

    const Int U_col =  S(lvl)(kid);
    Int U_row       =  0;
	  
    //const Int scol = LU(U_col)(U_row).scol; //Not used
    //const Int ecol = LU(U_col)(U_row).ecol; //Not used


    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n  LVL=%d  ----- kid: %d --\n\n",
	   lvl, kid);
    #endif
    


    //Do all domains (old sublevel 0)
    //If this works well we can move this for into the function
    for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("UPPER, kid: %d k: %d \n",
	       kid, k);
	#endif
	
	t_upper_col_factor(kid, team_leader, 
			   lvl, 0, 
			   k,
			   BASKER_FALSE);


      }//over all columns / domains / old sublevel 0

    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n\n done with UPPER, kid: %d \n\n\n", kid);
    #endif

    
    // return;


    //------Need because extend does not 
    //-------------Barrier--Between Domains-------------
    Int my_leader = find_leader(kid, 0);
    Int b_size    = pow(2,1);
    //barrier k = 0 usedl1
    t_basker_barrier(thread,kid,my_leader,
		     b_size, 0, LU(U_col)(U_row).scol, 0);
    
   

    
    

    //----------------Sep level upper tri-------------
    for(Int l = 1; l < (lvl); ++l)
      {
	
	//for(Int k = 0; k < 1; ++k)
	for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
	  {
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("\n\nSep, upper update, kid: %d k=%d \n\n",
		   kid, k+LU(U_col)(U_row).scol);
	    #endif
	    
	    t_add_extend(thread, kid,lvl,l-1, k, 
			 LU(U_col)(U_row).scol, 
			 BASKER_FALSE);

	    if(kid%((Int)pow(2,l)) == 0)
	      {
		my_leader = find_leader(kid,l);
		b_size    = pow(2,l+1);

		#ifdef BASKER_DEBUG_NFACTOR_COL2
		printf("\n\n\n SEP UPPER, kid: %d \n\n",
		       kid);
		#endif

		t_upper_col_factor(kid, team_leader, 
				   lvl, l, 
				   k,
				   BASKER_FALSE);
		

	      }//if correct kid to do this sublevels upperfactor
	  }//over all columns
	//Trick to just add 1 in case of sep size == 1
	//t_basker_barrier(thread, kid, kid, 
	//		 1, 1, LU(U_col)(U_row).ncol+1, l-1);
	//t_basker_barrier(thread, kid, kid,
	//		 1, 2, LU(U_col)(U_row).ncol+1, l-1);

      }//for - over all sublevel 1...lvl-2
    
    

    //---------Lower Factor (old sublevel lvl-1)-------
    
    //printf("\n\n");
    //printf("lower team size: %d \n", thread.team_size());
    //thread.team_barrier();
    my_leader = find_leader(kid, lvl-1);
    b_size    = pow(2,lvl);
    // printf("[3] barrier test, kid: %d leader: %d b_size: %d lvl: %d \n",
    //	   kid,  my_leader, b_size, lvl);
    t_basker_barrier(thread, kid, my_leader,
		     b_size, 3, LU(U_col)(U_row).scol, 0);


    


    //printf("\n\n======= LOWER, KID: %d ======= \n\n", kid);


   
    //return;
    
    //if(lvl < 2)
      {
	//for(Int k=0; k < 1; ++k)
	for(Int k = 0; k < LU(U_col)(U_row).ncol; ++k)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("lower_update, kid: %d k: %d \n",
	       kid, k);
	#endif

	//printf("test: %d \n", LU(U_col)(U_row).scol);
       
	t_add_extend(thread, kid,lvl,lvl-1, k,
		     LU(U_col)(U_row).scol,
		     BASKER_TRUE);
	Entry pivot = 0;
	if((kid%(Int)(pow(2,lvl))) == 0)
	  {

	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL2
	    printf("lower factor, kid: %d k: %d \n",
		   kid, k);
	    #endif
	    
	    t_lower_col_factor(kid, team_leader, 
			       lvl, lvl-1, 
			       k, pivot);
	  }

	
	//nee barrier if multiple thread uppdate
	//thread.team_barrier();
	my_leader = find_leader(kid, lvl-1);
	b_size    = pow(2,lvl);
	//printf("barrier test, leader 4: %d b_size: %d lvl: %d \n",
	//     my_leader, b_size, lvl);
	t_basker_barrier(thread, kid, my_leader,
			 b_size, 4, k, lvl-1);
	
	#ifdef BASKER_DEBUG_NFACTOR_COL2
	printf("lower diag factor, kid: %d k: %d \n",
	       kid, k);
	#endif
	
	t_lower_col_factor_offdiag2(kid, lvl, lvl-1, k, pivot);
	//thread.team_barrier();
	my_leader = find_leader(kid, lvl-1);
	b_size    = pow(2,lvl);
	//printf("barrier test 5, leader: %d b_size: %d lvl: %d \n",
	//     my_leader, b_size, lvl);
	t_basker_barrier(thread, kid, my_leader,
		     b_size, 5, k, lvl-1);
      }
      //Trick for sep = 1
      //t_basker_barrier(thread, kid, kid,
      //	       1, 1, LU(U_col)(U_row).ncol+1, lvl-1);
      //t_basker_barrier(thread, kid, kid,
      //	       1, 1, LU(U_col)(U_row).ncol+1, lvl-1);

      }
    

    
  }//end t_nfactor_sep2

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_add_extend
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
    Int my_leader = find_leader(kid,l);
    Int b_size    = pow(2,l+1);

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

	//This will do the correct spmv
	t_upper_col_factor_offdiag2(kid, lvl, sl,l, k, lower);
	
	//Barrier--Start
	my_leader = find_leader(kid,sl);
	b_size    = pow(2,sl+1);
	//printf("[1] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier(thread, kid, my_leader,
			 b_size, 1, k+k_offset, sl);
	//Barrier--End

	if(kid%((Int)pow(2,sl))==0)
	  {
	    t_dense_blk_col_copy_atomic2(kid, my_leader,
					 lvl, sl, l, k, lower);
	  }

	//Barrier--Start
	//printf("[2] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
	//     kid, my_leader, k, sl);

	t_basker_barrier(thread, kid, my_leader,
			 b_size, 2, k+k_offset, sl);
	
      }//over all sublevels

    t_dense_copy_update_matrix2(kid, my_leader, lvl, l, k);


  }//end t_add_add

  //uses local idx for local blks X
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag2
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k, 
   const BASKER_BOOL lower
   )
  {

    const Int my_leader = (sl==0)?(kid):find_leader(kid,sl-1);
    
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

    Int my_row_leader  = S(0)(find_leader(kid,lvl-1));
    //Int my_new_row = 
    // L_col - my_row_leader;
    Int U_row     = L_col-my_row_leader;
   

    const Int X_col     = S(0)(my_leader);
    Int X_row     = l+1; //this will change for us 
    //const Int X_col = S(l)(my_leader);
    //Int x_row = 

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
	    
    t_dense_back_solve_offdiag(kid,
			       L_col, L_row,
			       X_col, X_row,
			       k, col_idx_offset,
			       U.val,
			       U.row_idx,
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

	    t_dense_back_solve_offdiag(kid,
				       L_col, L_row,
				       X_col, X_row,
				       k, col_idx_offset,
				       U.val,
				       U.row_idx,
				  U.col_ptr(k+1)-U.col_ptr(k),
				       U.col_ptr(k),
				       BASKER_FALSE
				       );

	  }//for --over remainder blks
      }//if --lower
  
  }//end t_upper_col_factor_offdiag2()

 
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry, Exe_Space>::t_dense_blk_col_copy_atomic2
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
    //Int loop_col_idx     = S(l)(kid);

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
	//INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws; //NDE - warning: unused 
	Int      p_sizeL  = LL(leader_idx)(blk).p_size;
	ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
	INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
	//const Int ws_size = LL(my_idx)(blk).iws_size;
	//Int       p_size  = LL(my_idx)(blk).p_size;
	Int       *color  = &(ws[0]);
	//Int     *pattern  = &(color[ws_size]); 
	//Int      brow     = LL(my_idx)(blk).srow;
	//Int      browL    = LL(leader_idx)(blk).srow;

	#ifdef BASKER_DEBUG_NFACTOR_COL2
	if(lower == BASKER_TRUE)
	  {
	printf("kid: %d  COPY INDEX %d %d to %d %d \n",
	       kid, my_idx, blk, leader_idx, blk);
	  }
	#endif

	//Int *colorL   = &(wsL(0));
	//Int *patternL = &(colorL[ws_size]);

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
	    //Int jj = pattern[j];
	    //color[jj-brow] = 0;
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
	
	    if(X(jj) != (Entry)(0) )
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
	    //p_size = 0; //not needed
	  }//over all blks
	  }
      }//if not team_leader
   
  }//end t_dense_blk_col_copy_atomic2()

  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_dense_copy_update_matrix2
  (
   const Int kid,
   const Int team_leader,
   const Int lvl,
   const Int l, 
   const Int k
   )
  {
    //printf("\n\n\n\n");
    //printf("-----------------copy_update_matrx----------");
    //printf("\n\n\n\n");

    const Int leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;

    //Over each blk    
//    Int last_blk = l+2;//NDE - warning: unused 
    /*
    if(lvl ==(l+1))
      {
	last_blk = LL_size(leader_idx);
      }
    */
    


    for(Int TEMP=0; TEMP < 1; TEMP++)
      {
    //Copy B -> C
    Int bl = l+1;
    Int A_col = S(lvl)(kid);

    Int my_row_leader = S(0)(find_leader(kid,lvl-1));
    //Int my_new_row = 
    // S(bl)(kid) - my_row_leader;
    Int A_row = S(bl)(kid) - my_row_leader;
    


    //Int CM_idx = kid;
    
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
    
    //Int team_leader   = find_leader(kid, l);
    ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
    //const Int brow    = LL(leader_idx)(bl).srow;
    //const Int nrow    = LL(leader_idx)(bl).nrow;
    //Int p_size        = LL(leader_idx)(bl).p_size;
    //const Int ws_size = LL(leader_idx)(bl).iws_size;
    //Int *color        = &(ws(0));
    //Int *pattern      = &(color[ws_size]);


    #ifdef BASKER_DEBUG_NFACTOR_COL2
    
    printf("copy, kid: %d bl: %d  A: %d %d \n", 
	   kid, bl, A_col, A_row);
    #endif

    //const Int bbcol = B.scol;
   
    for(Int i = B.col_ptr(k);
	i < B.col_ptr(k+1); ++i)
      {
	Int B_row = B.row_idx(i);	
	//Int j = gperm(B_row+B.srow);

        #ifdef BASKER_DEBUG_NFACTOR_COL22
	printf("Scanning_2 A: %d %d lvl %d l: %d bl:%d brow: % d %d K: %d \n",
	       B_row, j, lvl, l, bl, brow, B.srow, kid);
	   
	printf("Adding Aval: %f to xval: %f \n", 
	       X[B_row], B.val(i));
        #endif
	 
	X(B_row) += B.val(i);

      }//end for over all nnz
      }//DUMMY
        
    //-------------move into C------------------- 
    //(Right now have to do dense but like to do sparse)

//    last_blk = LL_size(leader_idx); //NDE - warning: unused 
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    
    for(Int TEMP = 0; TEMP < 1; ++TEMP)
      {

	Int bl = l+1;
        //Int A_col = S(lvl)(kid);
  
    /*
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
     */

    //Int CM_idx = kid;
    ENTRY_1DARRAY   X   = LL(leader_idx)(bl).ews;
    INT_1DARRAY    ws   = LL(leader_idx)(bl).iws;
    //const Int   ws_size = LL(leader_idx)(bl).ews_size;
    //const Int      brow = LL(leader_idx)(bl).srow;
    const Int      nrow = LL(leader_idx)(bl).nrow;
    //Int p_size          = LL[leader_idx][bl].p_size;

    //For recounting patterns in dense blk
    //Need better sparse update
    //Int p_count  =0 ; 
    
    Int *color   = &(ws(0));
    //Int *pattern = &(color[ws_size]);
    
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
	
	
	if(X(j) != (Entry)(0) )
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

  }//end t_dense_copy_update_matrix2()


  //local idx and lock blk
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag2
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
    
    
    
    //printf("OFF_DIAG_LOWER, kid: %d leaderid: %d t_size: %d \n",
    //	   kid, leader_id, lteam_size);
    
    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; 
	 L_row < LL_size(L_col);
	 X_row+=(lteam_size), L_row+=(lteam_size))
    
     { 


       //printf("OFF_DIAG_LOWER. kid: %d k: %d  U: %d %d L: %d %d X: %d %d pivot: %f \n", kid, k, U_col, U_row, L_col, L_row, X_col, X_row, pivot);

       /*old
       t_back_solve_offdiag(leader_id,
			    L_col, L_row,
			    X_col, X_row,
			    k, col_idx_offset,
			    U.val, U.row_idx,
			    //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			    U.col_ptr(k+1)-U.col_ptr(k),
			    //U.col_ptr[k-bcol],
			    U.col_ptr(k),
			    BASKER_TRUE);
       */

       //We might still have to do sparse here
              t_dense_back_solve_offdiag(leader_id,
			    L_col, L_row,
			    X_col, X_row,
			    k, col_idx_offset,
			    U.val, U.row_idx,
			    U.col_ptr(k+1)-U.col_ptr(k),
			    U.col_ptr(k),
			    BASKER_TRUE);


	      //printf("MOVING OFF, kid: %d k: %d L: %d %d X %d %d \n",
	      //     kid, k, L_col, L_row, X_col, X_row);

	      t_dense_move_offdiag_L(leader_id, 
			L_col, L_row,
			X_col, X_row,
			k, pivot);
       
     }//end for over all offdiag blks

  }//end t_lower_col_factor_offdiag2()


}//end namespace BaskerNS

#endif //end ifndef
