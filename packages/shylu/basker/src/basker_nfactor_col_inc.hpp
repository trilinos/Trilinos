#ifndef BASKER_NFACTOR_COL_INC_HPP
#define BASKER_NFACTOR_COL_INC_HPP

#include "basker_decl.hpp"
#include "basker_matrix_decl.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"
#include "basker_types.hpp"
#include "basker_stats.hpp"
#include "basker_thread.hpp"

#include "basker_nfactor_blk.hpp"
#include "basker_nfactor_blk_inc.hpp"


#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif 

#ifdef BASKER_DEBUG
#include <assert.h>
#endif

namespace BaskerNS
{

  //move laster
    //old using global index for local blk
  #ifdef BASKER_KOKKOS
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int, Entry,Exe_Space>::t_nfactor_sep_old
  (Int kid,
   Int lvl, 
   Int team_leader,
   const TeamMember &thread)
  #else
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_nfactor_sep_old
    (Int kid,
     Int lvl,
     Int team_leader)
  #endif
  {
  
    #ifdef BASKER_DEBUG_TIME
    double upper_time = 0; 
    double upper_diag_time = 0;
    double lower_time = 0;
    double lower_diag_time = 0;
    double reduce_time = 0;
    double copy_time = 0;
    double barrier_time = 0;
    #endif
    
    Int U_col = S[lvl][kid];
    Int U_row = 0;
    
    Int scol = LU[U_col][U_row].scol;
    Int ecol = LU[U_col][U_row].ecol;

    /*
    if(lvl == 2)
      {
	 ecol = scol+1;
      }
    */
    for(Int k = scol; k < ecol; k++)
      {

	#ifdef BASKER_DEBUG_NFACTOR_COL
	printf("-------------------k=%d--------------------\n",
	       k);
	#endif
	
       
	for(Int l = 0; l < lvl; l++) //sublevel
	  // for(Int l = 0; l < 1; l++)
	  {	    
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL
	    printf("\n\n--------sl=%d-------------\n\n", l);
	    #endif

	    //Used for barrier
	    Int my_leader = find_leader(kid, l);
	    Int b_size = pow(2,l+1);


	    //Can remove in future, don't need anymore
	    BASKER_BOOL sep_lvl_flg = BASKER_FALSE;
	  

	    #ifdef BASKER_DEBUG_TIME
	    Kokkos::Impl::Timer timer;
	    #endif
	    
	    //-----------Upper Col Factor------------//
	    if(kid%((Int)pow(2,l)) == 0)
	      {
                #ifdef BASKER_DEBUG_NFACTOR_COL
		printf("------STARTING TRI-SOLVE: %d %d -----\n", 
		   kid, l);
                #endif
		
		//if((kid==0)||(kid==1)||(kid==2)||(kid==3))
		t_upper_col_factor(kid, team_leader, lvl, l, k, sep_lvl_flg);
		
	      }//end if kid ... I should do an upper_col_factor
	    

	    #ifdef BASKER_DEBUG_TIME
	    upper_time += timer.seconds();
	    timer.reset();
	    #endif

	    
	    #ifdef BASKER_MULTIPLE_UPPER
	    //Upper_col_factor_updates
	    my_leader = (l==0)?kid:find_leader(kid,l-1);
	    b_size = pow(2,l);
	    
	    //printf("t_upper barrier, kid: %d b_size: %d \n",
	    //	   kid, b_size);
	    
	    //printf("Barrier0. kid: %d lkid: %d bsize: %d l: %d\n",
	    //	   kid, my_leader, b_size, l);

	    //t_basker_barrier(thread,my_leader,l,0,b_size);
	    t_basker_barrier(thread,kid,my_leader,b_size,0,k,l);

	    //printf("Leave Barrier0 kid: %d \n", kid);

	    #ifdef BASKER_DEBUG_TIME
	    barrier_time += timer.seconds();
	    timer.reset();
	    #endif


	    t_upper_col_factor_offdiag(kid,lvl,l,k);
	    my_leader = find_leader(kid,l);
	    b_size = pow(2,l+1);


	    #ifdef BASKER_DEBUG_TIME
	    upper_diag_time += timer.seconds();
	    timer.reset();
	    #endif

	    #else

	    my_leader = find_leader(kid,l);
	    b_size = pow(2,l+1);



	    #endif
	    

	    

	    #ifdef BASKER_OLD_BARRIER
	    thread.team_barrier();
	    #else
	    //Int my_leader = find_leader(kid, l);
	    //Int b_size = pow(2,l+1);
	    //Debug


	    
	    //printf("Barrier1, kid: %d lkid: %d bsize: %d l: %d \n",
	    //	   kid, my_leader, b_size, l);
	    //t_basker_barrier(thread, my_leader, l, 1, b_size);
	    t_basker_barrier(thread, kid, my_leader, b_size, 1,k, l);

	    //printf("Leave Barrier1: %d \n", kid);
	    #endif


	    
	    #ifdef BASKER_DEBUG_TIME
	    barrier_time += timer.seconds();
	    timer.reset();
	    #endif


	    #ifdef BASKER_2DL
	    if(kid%((Int)pow(2,l))==0)
	      {

		//Rename---not really atomic anymore
		//if((kid==0)||(kid==1)||(kid==2)||(kid==3))
		t_blk_col_copy_atomic(kid, team_leader, lvl, l, k);
	      }
	    #endif


	    #ifdef BASKER_DEBUG_TIME
	    reduce_time += timer.seconds();
	    timer.reset();
	    #endif


	    #ifdef BASKER_OLD_BARRIER
	    thread.team_barrier();
	    #else
	    //printf("Barrier2, kid: %d lkid: %d bsize: %d l: %d\n",
	    //	   kid, my_leader, b_size, l);
	    //t_basker_barrier(thread, my_leader, l, 2,b_size);
	    t_basker_barrier(thread, kid, my_leader, b_size, 2, k,l);

	    //printf("Leave Barrier2 kid: %d \n", kid);
	    #endif


	    #ifdef BASKER_DEBUG_TIME
	    barrier_time += timer.seconds();
	    timer.reset();
	    #endif
	    

            #ifdef BASKER_ATOMIC_2 //O(dim(S)*p)
	    //t_n_col_copy_atomic(kid, team_leader, lvl, l, k);
	    #endif

	    //#ifdef BASKER_KOKKOS
	    //thread.team_barrier();//might have to be added back
	    //#else
	    //t_col_barrier(kid,lvl,l);
	    //#endif

	    if(kid%((Int)pow((double)2,(double)l+1)) == 0)
	      {  
		#ifdef BASKER_ATOMIC
		//t_col_copy_atomic(kid, team_leader, lvl, l, k);
		//if((kid==0)||(kid==1)||(kid==2)||(kid==3))
		t_copy_update_matrix(kid, team_leader,lvl, l, k);
		#else

		#endif
	      }//end if(kid) should I do reduce
	    
	    
	    #ifdef BASKER_DEBUG_TIME
	    copy_time += timer.seconds();
	    timer.reset();
	    #endif

	    //Used if O(dim(S)*num_threads) Atomic
	    #ifdef BASKER_ATOMIC_2
	    //thread.team_barrier();
	    #endif

	    //-----------------lower factor-----------------//  
	    Entry pivot = 0;
	    if(((l+1)==lvl)&&
	       (kid%((Int)pow((double)2,(double)l+1)) == 0))
	      {
		//if((kid==0)||(kid==1)||(kid==2)||(kid==3))
		t_lower_col_factor(kid, team_leader,
				   lvl, l, k, pivot);
	      }//end if I should perfom factor
	   

	    #ifdef BASKER_DEBUG_TIME
	    lower_time += timer.seconds();
	    timer.reset();
	    #endif


	    //printf("Barrier3, kid: %d lkid: %d bsize: %d l: %d\n", kid, my_leader, b_size, l);

		t_basker_barrier(thread,kid, my_leader, b_size,
				 3, k,l);

		//printf("Leave Barrier3, kid: %d \n", kid);

	    #ifdef BASKER_MULTIPLE_LOWER
	    if((l+1)==lvl)
	      {
	    //need a barrier to make sure upper_lower done
	    //printf("Barrier: lower_diag off-diag \n");
	    //t_basker_barrier(thread, my_leader, l, 3,
		//		     b_size);

	

		

	    #ifdef BASKER_DEBUG_TIME
	    barrier_time += timer.seconds();
	    timer.reset();
	    #endif


	    t_lower_col_factor_offdiag(kid,lvl,l,k, pivot);


	    #ifdef BASKER_DEBUG_TIME
	    lower_diag_time += timer.seconds();
	    timer.reset();
	    #endif

	      }
	    #endif
 
	    

	    //We might be able to remove this in the future
	    #ifdef BASKER_OLD_BARRIER
	    thread.team_barrier();
	    #else
	    //my_leader = 0;
	    //b_size = num_threads;
	    //printf("Barrier4, kid: %d lkid: %d bsize: %d l: %d\n",
	    //	   kid, my_leader, b_size, l);
	    //t_basker_barrier(thread, my_leader, l,4, b_size);
	    t_basker_barrier(thread,kid, my_leader,b_size,4,k,l);

	    //printf("Leave Barrier4 kid: %d \n", kid);
	    #endif
	   


	    #ifdef BASKER_DEBUG_TIME
	    barrier_time += timer.seconds();
	    timer.reset();
	    #endif
	   
	  }// end over all sublvl
      }//over each column


    #ifdef BASKER_DEBUG_TIME
    printf("COL-TIME KID: %d  UPPER: %f UPPER-DIAG: %f LOWER: %f LOWER_DIAG: %f REDUCE: %f COPY: %f BARRIER: %f \n", 
	   kid, upper_time, upper_diag_time, lower_time, lower_diag_time, reduce_time, copy_time, barrier_time);
    #endif



    #ifdef BASKER_OPS_COUNT


    for(Int l = 0; l < lvl; l++)
      {

	printf("OPS.   KID : %d  LVL: %d OPS : %d \n",
	       kid, l, thread_array[kid].ops_counts[l][0]);
	thread_array[kid].ops_count[1][0] = 0;

      }

    #endif

    return 0;
  }//end t_nfactor_sep()


    //global idx for local number
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_upper_col_factor_old
  (Int kid,
   Int team_leader,
   Int lvl,
   Int l,
   Int k, 
   BASKER_BOOL sep_flg)
  {

    
    //Get needed variables
    Int L_col = S[l][kid];
    Int L_row = 0;
    Int U_col = S[lvl][kid];
    Int U_row = (lvl==1)?(kid%2):S[l][kid]%LU_size[U_col];
    #ifdef BASKER_2DL
    //Int X_col = L_col;
    Int X_col = S[0][kid];
    Int X_row = l; //X_row = lower(L)
    Int col_idx_offset = 0;
    #else
    #ifdef BASKER_ATOMIC_2
    Int X_col = kid;
    #else
    Int X_col = team_leader;
    #endif
    #endif

    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
    printf("kid %d, upper using L: %d %d  U: %d %d  X %d %d\n",
	   kid, L_col, L_row, U_col, U_row, X_col, X_row);
    #endif
    //end get needed variables//

    BASKER_MATRIX        &L = LL[L_col][L_row];
    BASKER_MATRIX        &U = LU[U_col][U_row]; 
    
    //Removed
    //U.fill();

    //Ask C++ guru if this is ok
    BASKER_MATRIX        *Bp;
    Int                  bbcol = k;
    //if(sep_flg == BASKER_FALSE)
    if(l == 0 )
      {
        Bp = &(AVM[U_col][U_row]);
	bbcol = Bp->scol;
      }
    else
      {
	Bp = &(thread_array[kid].C);
	//printf("Using temp matrix, kid: %d\n", kid);
	//Bp->print();
      }
    BASKER_MATRIX       B = *Bp;

    //B.print();

    //BASKER_MATRIX_VIEW   &B = AV[U_col][U_row];
    //B.init_perm(&gperm);  //provide perm vector
    //B.init_offset(k,0);  //col,offset

    #ifdef BASKER_2DL
    INT_1DARRAY ws  = LL[X_col][X_row].iws;
    Int    ws_size  = LL[X_col][X_row].iws_size;
    ENTRY_1DARRAY X = LL[X_col][X_row].ews;
    #else
    INT_1DARRAY  ws = thread_array[kid].iws;
    Int     ws_size = thread_array[kid].iws_size; 
    ENTRY_1DARRAY X = thread_array[X_col].ews;
    #endif
    
    Int brow = U.srow;
    Int bcol = U.scol;

    Int *color     = &(ws[0]);
    Int *pattern   = &(color[ws_size]);
    
    Int j, t, xnnz;
    Int top = ws_size;

    Int unnz = 0 ; 
    if(k != U.scol)
      {unnz = U.col_ptr[k-U.scol];}

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
    ASSERT(top == ws_size);
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

    //Bgood(remove)
    //Int i = B.offset;
    //Int mi = B.m_offset;
    //for(i = B.offset; i < mi; i++)
    for(Int i = B.col_ptr[k-bbcol]; i < B.col_ptr[k-bbcol+1]; i++)
      {
	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid>=0)
	printf("kid: %d index: %d %d\n", 
	       kid, i, B.row_idx(i));
	#endif

	//Bgood(remove)
	//Note: add a break here so we don't have to go over all
	//if(B.good(i) == L.max_idx)
	//if(B.good(i) == BASKER_MAX_IDX)
	//  {continue;} //Need a way to make this a break
	j = B.row_idx(i);
	
	#ifdef BASKER_2DL
	X[j-brow] = B.val(i);
	#else
	X[j] = B.val(i);
	#endif

	#ifdef BASKER_DEBUG_NFACTOR_COL
	if(kid>=0)
        printf("kid: %d i: %d  val: %g  top: %d \n", 
	       kid,i, B.val(i), top);
	#ifdef BASKER_2DL
	if(kid>=0)
        printf("kid: %d Nx in Ak %d %g color = %d \n",
               kid, j, X[j-brow],  ws[j-brow] );
	#else
	if(kid>=0)
        printf("kid: %d Nx in Ak %d %g color = %d \n",
               kid, j, X[j],  ws[0 + j] );
	#endif
        #endif

	#ifdef BASKER_2DL
	if(color[j-brow] == 0)
	#else
	if(color[j] == 0)
	#endif
	  {
	    #ifdef BASKER_INC_LVL
	    t_local_reach_selective(kid, l, l, j, &top);
	    #else
	    //t_local_reach(kid, l, l, j, &top); //Note: comeback
	    #endif
	  }//if not colored
      }//end over each nnz in column
    xnnz = ws_size - top;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid>=0)
    printf("xnnz: %d ws_size: %d top: %d , kid: %d\n",
	   xnnz, ws_size, top, kid);
    #endif


    
    //Count ops to show imbalance
    #ifdef BASKER_COUNT_OPS
    thread_array[kid].ops_counts[0][l] += xnnz;
    #endif


    //WE SHOUD DO A UNNZ COUNT
     //count number of nnz
     for(Int i=top; i < ws_size; i++)
       {
	 j = pattern[i];
	 t = gperm[j];
	 //if(t == L.max_idx)
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


     //#ifdef BASKER_DEBUG_NFACTOR_COL
     if(unnz+ucnt-1 > uunnz)
       {printf("kid: %d col: %d need to realloc, unnz: %d ucnt: %d uunnz: %d U_col: %d U_row: %d \n", kid, k, unnz, ucnt, uunnz, U_col, U_row);}
     //#endif


     #ifdef BASKER_INC_LVL
     //Note we have not made a t_back_solve_atomic_selective
     t_back_solve_selective(kid, l,l, k, top, xnnz);
     
     #else //BASKER_INC_LVL


     #ifdef BASKER_ATOMIC_2
     t_back_solve(kid, l,l, k, top, xnnz);
     #else
     t_back_solve_atomic(kid, team_leader,
		         lvl, l,k , top,
		         xnnz);
     #endif
     #endif //BASKER_INC_LVL

     //move over nnz to U 
     for(Int i = top; i < ws_size; i++)
       {
	 j = pattern[i];
	 t = gperm[j];
	 
	 #ifdef BASKER_DEBUG_NFACTOR_COL
	 #ifdef BASKER_2DL
	 if(kid>=0)
	   printf("considering j: %d t:%d val: %e, kid: %d \n",
		  j, t, X[j-brow], kid);
	 #else
	 if(kid>=0)
	   printf("considering j: %d t:%d val: %e, kid: %d \n",
		  j, t, X[j], kid);
	 #endif
         #endif

	 #ifdef BASKER_2DL
	 if(X[j-brow] !=0)
	 #else
	 if(X[j] != 0)
	 #endif
	  {

	    //SHOULD:: REMOVE, check if this fine
	    //if( (t != L.max_idx) && (t >= L.scol) && 
	    //(t<(L.scol+L.ncol)))
	    //if(t!=L.max_idx)
	    if(t != BASKER_MAX_IDX)
              {
                #ifdef BASKER_DEBUG_NFACTOR_COL
		if(kid>=0)
                printf("kid: %d adding x[%d] to U\n", kid, j); 
                #endif

                U.row_idx[unnz] = gperm[j];
		#ifdef BASKER_2DL
		U.val[unnz] = X[j-brow];
		#else
                U.val[unnz] = X[j];
		#endif
                unnz++;
		#ifdef BASKER_2DL
		X[j-brow] = 0;
		#else
                X[j] = 0;
		#endif
		
              }//if in U
            else
              {
		//ASSERT if 0==1
		#ifdef BASKER_2DL
		printf("----Error--- kid: %d extra L[%d]=%f \n",
		       kid, j, X[j-brow]);
                #endif
              }//lower
	  }//end not 0
       }//over all x

     U.col_ptr[k+1-bcol] = unnz;

     #ifdef BASKER_2DL
  
     
     #ifdef BASKER_DEBUG_NFACTOR_COL
     printf("JB TEST kid: %d lvl: %d l: %d L_col: %d size: %d Lrow: %d Xrow: %d\n",
	    kid, lvl, l, L_col, LL_size[X_col], L_row, X_row);
     #endif


     #ifndef BASKER_MULTIPLE_UPPER

     //----------------------Update offdiag-----------------//
     X_row++;
     L_row++;
     for(; X_row < LL_size[X_col]; X_row++, L_row++)
       {

	 //printf("xrow: %d \n", X_row);
	 BASKER_BOOL A_option = BASKER_FALSE;
	 //if((kid == team_leader) && blk_row == 1)
	 //  {A_option = BASKER_TRUE;}
	    
	 //back_solve 
	 #ifdef BASKER_INC_LVL
	 t_back_solve_offdiag_selective(kid,
			      L_col, L_row,
			      X_col, X_row,
			      k, col_idx_offset,
			      U.val,
			      U.row_idx,
		     U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			      U.col_ptr[k-bcol],
			      A_option);

	 #else

	  printf("t_bsolve_d test, kid: %d xsize: %d\n",
	 	kid, 
	 	U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol]);
	 
	 t_back_solve_offdiag(kid,
			      L_col, L_row,
			      X_col, X_row,
			      k, col_idx_offset,
			      U.val,
			      U.row_idx,
		     U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			      U.col_ptr[k-bcol],
			      A_option);
	 #endif

       }//end for over all offdiag
     #endif
     
     #endif

     //Bgood(removed)
     //B.flip_base();
     
     return 0;
  }//end t_upper_col_factor()

  //uses glb idx for local blks
  //This will do the off diagonal updates now
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag_old
  (
   Int kid,
   Int lvl,
   Int l,
   Int k
   )
  {

    // printf("t_upper_col_factor_offdiag called kid: %d \n", 
    //	   kid);

   
    //Note: We can remove this and pass it
    //Might be faster
    Int my_leader = kid ;
    if(l>0)
      {
	my_leader = find_leader(kid, l-1);
      }
    //Note: We can remove this and pass it
    //Might be faster
    //hard coding for 2 right now
    int lteam_size = pow(2, l);
    
    #ifdef BASKER_2DL
    Int L_col = S[l][my_leader];
    Int L_row = 0;
    Int U_col = S[lvl][kid];
    Int U_row = (lvl==1)?(kid%2):S[l][kid]%LU_size[U_col];
    Int X_col = S[0][my_leader];
    Int X_row = l; //this will change for us 
    Int col_idx_offset = 0;
    BASKER_MATRIX        &U = LU[U_col][U_row];
    Int bcol = U.scol;
    #else
    BASKER_ASSERT(0==1, "t_upper_col_factor_offdiag, only work with with 2D layout");
    #endif

    #ifdef BASKER_2DL
    INT_1DARRAY ws  = LL[X_col][X_row].iws;
    Int    ws_size  = LL[X_col][X_row].iws_size;
    ENTRY_1DARRAY X = LL[X_col][X_row].ews;
    #else
    BASKER_ASSERT(0==1, "t_upper_col_factor_offdiag, only works with 2D layout");
    #endif
    

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("Upper_fact_offdiag, kid: %d leader: %d l: %d lvl: %d ltsize: %d works_size: %d X: %d %d L_row: %d %d \n",
	   kid, my_leader, l, lvl,lteam_size, LL_size[X_col], X_col, X_row, L_col, L_row);
    #endif

     //----------------------Update offdiag-----------------//
    //Want to assign in a round-robin fashion

    //X_row++;
    //L_row++;
    X_row += (kid-my_leader)+1;
    L_row += (kid-my_leader)+1;
    for(; X_row < LL_size[X_col]; X_row+=lteam_size, L_row+=lteam_size)
       {

	 #ifdef BASKER_DEBUG_NFACTOR_COL
	 printf("OFF-DIAG, kid: %d, l: %d  X: %d %d L: %d %d \n",
		kid, l,X_col, X_row, L_col, L_row);
	 #endif

	 BASKER_BOOL A_option = BASKER_FALSE;
	    
	 //back_solve 
	 t_back_solve_offdiag(kid,
			      L_col, L_row,
			      X_col, X_row,
			      k, col_idx_offset,
			      U.val,
			      U.row_idx,
		     U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			      U.col_ptr[k-bcol],
			      A_option);


       }//end for over all offdiag     

    return 0;
  }//end t_upper_col_factor_offdiag()


  //old uses global idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::t_lower_col_factor_old
  (Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k, 
   Entry &opivot
   )
  {
    //Get needed variables
    Int L_col = S[lvl][kid];
    Int L_row = 0;
    Int U_col = S[lvl][kid];
    Int U_row = LU_size[U_col]-1;
    #ifdef BASKER_2DL
    Int X_col = S[0][kid];
    Int col_idx_offset = 0;
    #else
    Int X_col = team_leader;
    #endif

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("LOWER_COL_FACTOR kid: %d \n", kid);
    printf("kid %d using L: %d %d  U: %d %d  X %d \n",
	   kid, L_col, L_row, U_col, U_row, X_col);
    #endif
    //end get needed variables

    BASKER_MATRIX        &L = LL[L_col][L_row];
    BASKER_MATRIX        &U = LU[U_col][U_row]; 
    //U.fill(); removed
    BASKER_MATRIX        &B = thread_array[kid].C;

    //BASKER_MATRIX_VIEW   &B = AV[U_col][U_row];
    // B.init_perm(&gperm);  //provide perm vector
    //B.init_offset(k,0);   //col,offset
    
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

    #ifdef BASKER_2DL
    INT_1DARRAY  ws = LL[X_col][l+1].iws;
    Int     ws_size = LL[X_col][l+1].iws_size;
    ENTRY_1DARRAY X = LL[X_col][l+1].ews;
    #else
    INT_1DARRAY  ws = thread_array[kid].iws;
    Int     ws_size = thread_array[kid].iws_size; 
    ENTRY_1DARRAY X = thread_array[X_col].ews;
    #endif

    Int brow     = U.srow;
    Int bcol     = U.scol;

    Int lval       = L.col_ptr[k-bcol];
    Int uval       = U.col_ptr[k-bcol];
    
    Int *color     = &(ws[0]);
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
    printf("\n----------------- K = %d --------------\n", k);
    #endif
    	  
    value = 0.0;
    pivot = 0.0;
    //maxindex = A.max_idx;
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

   for(i = B.col_ptr[0]; i < B.col_ptr[1]; i++)
     {
       #ifdef BASKER_DEBUG_NFACTOR_COL
       //Bgood(remove)
       //B.good(i);
       #endif
       
       
       j = B.row_idx(i);
    
       if(j > (U.srow+U.nrow))
	 {
	   //printf("j continue \n");
	   break;
	 }
  
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid>=0)
       printf("j: %d i: %d \n", j, i);
       #endif

       #ifdef BASKER_2DL
       X[j-brow] = B.val(i);
       #else
       X[j] = B.val(i);
       #endif
  
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(kid>=0)
       printf("i: %d j: %d  val: %g  top: %d \n", 
	      i, gperm[j], B.val(i), top);
       #ifdef BASKER_2DL
       if(kid>=0)
       printf("Nxk in Ak %d %g color = %d \n",
	      j, X[j-brow],  
	      ws[j-brow ] );
       #else
       if(kid>=0)
       printf("Nx in Ak %d %g color = %d \n",
	      j, X[j],  
	      ws[j ] );
       #endif
       #endif

       #ifdef BASKER_2DL
       if(color[j-brow] == 0)
       #else
       if(color[j] == 0)
       #endif
	 {
	   //printf("doing reach: %d \n", j);
	   #ifdef BASKER_INC_LVL
	   t_local_reach_selective(kid,lvl,l+1,j, &top);
	   #else
	   //t_local_reach(kid,lvl,l+1, j, &top);
	   #endif
	 }
     }//over each nnz in the column
   xnnz = ws_size - top;
   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid>=0)
     {
       printf("xnnz: %d ws_size: %d top: %d \n",
	      xnnz, ws_size, top);
     }
   #endif


   #ifdef BASKER_OPS_COUNT
   thread_array[kid].ops_counts[0][l] += xnnz;
   #endif
   
  
   t_back_solve(kid, lvl,l+1, k, top, xnnz); // note: l not lvl given

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

       //printf("considering: %d %f \n", 
       //      j, value);

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
   ucnt = ws_size - top - lcnt +1;
   

   //----------------------------Sym-----
   //SYM
   if(Options.no_pivot == BASKER_TRUE)
     {
       maxindex = k;
       pivot = X[k-brow];
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
  
   gperm[maxindex] = k;
   gpermi[k] = maxindex;
    
   
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
   


   //printf("kid: %d lnnz: %d llnz: %d \n", 
   //kid, lnnz, llnnz);
   L.row_idx[lnnz] = maxindex;
   L.val[lnnz] = (Entry) 1.0;
   lnnz++;

   Entry lastU = (Entry) 0.0;
  
   //For every nnz in LHS
   for( i = top; i < ws_size; i++)
     {
       j = pattern[i];
       t = gperm[j];
       
       #ifdef BASKER_DEBUG_NFACTOR_COL
       if(k>=0)
	 printf("j %d  t %d , kid: %d \n", j, t, kid);
       #endif
       
       //if not zero
       #ifdef BASKER_2DL
       if(X[j-brow] !=0)
       #else
       if(X[j] != 0)
       #endif
	 {
           #ifdef BASKER_DEBUG_NFACTOR_COL
	   #ifdef BASKER_2DL
	   if(kid>=0)
	     printf("found value: %f at %d, kid: %d \n",
		    X[j-brow], j, kid);
	   #else
	   if(kid>=0)
	     printf("found value: %f at %d, kid: %d \n",
		    X[j], j, kid);
	   #endif
           #endif
	   
	   //if(t != L.max_idx)
	   if(t != BASKER_MAX_IDX)
	     {
	       if(t < k)
		 {
		   #ifdef BASKER_DEBUG_NFACTOR_COL
		   #ifdef BASKER_2DL
		   if(kid>=0)
		   printf("U insert: %f at %d \n",
			  X[j-brow], gperm[j]);
		   #else
		   if(kid>=0)
		   printf("U insert: %f at %d \n",
			  X[j], gperm[j]);
		   #endif
		   #endif
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
               #ifdef BASKER_DEBUG_NFACTOR_COL
	       #ifdef BASKER_2DL
	       if(kid>=0)
	       printf("inserting %f at %d into %d \n", 
		      X[j-brow]/pivot, j, lnnz );
	       #else
	       if(kid>=0)
	       printf("inserting %f at %d into %d \n", 
		      X[j]/pivot, j, lnnz );
	       #endif
               #endif
	       L.row_idx[lnnz] = j;
	       #ifdef BASKER_2DL
	       L.val[lnnz] = X[j-brow]/pivot;
	       #else
	       L.val[lnnz] = X[j]/pivot;
	       #endif
	       lnnz++;
	     }
	 }//end if() not zero             
       //move inside if not zero..... extra set ops not needed
  
       #ifdef BASKER_2DL
       X[j-brow] = 0;
       #else
       X[j] = 0;
       #endif
     }//if(x[j-brow] != 0)
   
   //Fill in last element of U
   U.row_idx[unnz] = k;
   U.val[unnz] = lastU;
   unnz++;
   
   xnnz = 0;
   top = ws_size;
   
   #ifdef BASKER_DEBUG_NFACTOR_COL
   printf("setting col: %d %d %d %d\n",
               k-bcol, cu_ltop, k+1-bcol, lnnz);
   #endif
   L.col_ptr[k-bcol] = cu_ltop;
   L.col_ptr[k+1-bcol] = lnnz;
   cu_ltop = lnnz;
   
   U.col_ptr[k-bcol] = cu_utop;
   U.col_ptr[k+1-bcol] = unnz;
   cu_utop = unnz;

   #ifdef BASKER_DEBUG_NFACTOR_COL
   if(kid>=0)
   printf("col_fact k: %d Unnz: %d   Lnnz: %d \n",
	  k, unnz, lnnz);
   #endif
  
   
   #ifndef BASKER_MULTIPLE_LOWER
   #ifdef BASKER_2DL
   //----Update offdiag (this can be done in parallel l>0)---//
   //-Need to fix
   #ifdef BASKER_DEBUG_NFACTOR_COL
   printf("k: %d -----TTTEST(blk_row-- %d %d \n",
     kid, lvl+1, LL_size[L_col]);
   printf("k: %d -----TTTEST(x_row)--- %d %d \n",
       kid, l+2, LL_size[X_col]);
   #endif
   for(Int blk_row = L_row+1, x_row = l+2; blk_row < LL_size[L_col]; blk_row++, x_row++)
     { 
       t_back_solve_offdiag(kid,
			    L_col, blk_row,
			    X_col, x_row,
			    k, col_idx_offset,
			    U.val, U.row_idx,
		  U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			    U.col_ptr[k-bcol],
			    BASKER_TRUE);
       t_move_offdiag_L(kid, 
			L_col, blk_row,
			X_col, x_row,
			k, pivot);
     }//end for over all offdiag blks
   #endif
   #endif

   #ifdef BASKER_DEBUG_NFACTOR_DEBUG
   print_factor(L,U);
   #endif

  //Bgood(remove)
  //B.flip_base();
   
   return 0;
  }//end t_lower_col_fact()

  //glb idx in local blk
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag_old
  (
   Int kid,
   Int lvl,
   Int l,
   Int k,
   Entry pivot
   )
  {
   
    #ifndef BASKER_MULTIPLE_LOWER
    BASKER_ASSERT(0==1,"BASKER_MULTIPLE_LOWER ERROR");
    return 0;
    #endif


    Int leader_id = find_leader(kid, l);
    Int lteam_size = pow(2,l+1);
    Int L_col = S[lvl][leader_id];
    Int L_row = 0;
    Int U_col = S[lvl][leader_id];
    Int U_row = LU_size[U_col]-1;
    Int X_col = S[0][leader_id];
    Int X_row = l+1;
    Int col_idx_offset = 0;
   
    BASKER_MATRIX        &L = LL[L_col][L_row];
    BASKER_MATRIX        &U = LU[U_col][U_row]; //U.fill();
    
    INT_1DARRAY  ws = LL[X_col][X_row].iws;
    Int     ws_size = LL[X_col][X_row].iws_size;
    ENTRY_1DARRAY X = LL[X_col][X_row].ews;

    Int brow     = U.srow;
    Int bcol     = U.scol;

    pivot = U.tpivot;

    //printf("lower_off, kid: %d leader_id: %d lsize: %d X: %d %d U %d %d L: %d %d \n",
    //	   kid, leader_id, lteam_size, 
    //	   X_col, X_row, U_col, U_row, L_col, L_row);


    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; L_row < LL_size[L_col]; X_row+=(lteam_size), L_row+=(lteam_size))
     { 


       // printf("OFF_DIAG_LOWER. kid: %d U: %d %d L: %d %d X: %d %d pivot: %f \n", kid, U_col, U_row, L_col, L_row, X_col, X_row, pivot);

       
       t_back_solve_offdiag(leader_id,
			    L_col, L_row,
			    X_col, X_row,
			    k, col_idx_offset,
			    U.val, U.row_idx,
		  U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
			    U.col_ptr[k-bcol],
			    BASKER_TRUE);
       t_move_offdiag_L(leader_id, 
			L_col, L_row,
			X_col, X_row,
			k, pivot);
       
     }//end for over all offdiag blks

    return 0;
  }//end t_lower_col_factor_offdiag()


  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_copy_update_matrix_old
  (Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k)
  {
    Int       leader_idx = S[0][kid];
    BASKER_MATRIX   &C = thread_array[kid].C;  
    Int nnz = 0;
    //COME BACK HERE

    //team_leader = find_leader(kid, l);
    //----------------Add A------------------
    //Over each blk
    Int last_blk = l+2;   
    if(lvl ==(l+1))
      {
	last_blk = LL_size[leader_idx];
      }
    // for(Int bl = l+1; bl < last_blk; bl++)
      {
	Int bl = l+1;
	Int A_col = S[lvl][kid];
	Int A_row = (lvl==1)?(2):S[bl][kid]%(LU_size[A_col]);
	Int CM_idx = kid;

	//Bgood(remove)
	//BASKER_MATRIX_VIEW &B = AV[A_col][A_row];
	//B.init_perm(&gperm);
	//B.init_offset(k, 0); //could be made faster
	BASKER_MATRIX  *Bp;
	if(A_row != (LU_size[A_col]-1))
	  {
	    //printf("upper picked, kid: %d \n", kid);
	    //printf("up: %d %d kid: %d \n",
	    //	   A_col, A_row, kid);
	    Bp = &(AVM[A_col][A_row]);
	  }
	else
	  {
	    //printf("lower picked, kid: %d\n", kid);
	    Bp = &(ALM[A_col][0]);
	  }

	BASKER_MATRIX   &B  = *Bp;


	//printf("ADDING UPDATES TO B\n");
	//B.info();
	//B.print();

	team_leader = find_leader(kid, l);
	//ENTRY_1DARRAY   X = LL[team_leader][bl].ews;
	ENTRY_1DARRAY   X = LL[leader_idx][bl].ews;
	//INT_1DARRAY    ws = LL[team_leader][bl].iws;
	INT_1DARRAY    ws = LL[leader_idx][bl].iws;
	//Int brow = LL[team_leader][bl].srow;
	//Int nrow = LL[team_leader][bl].nrow;
	Int brow = LL[leader_idx][bl].srow;
	Int nrow = LL[leader_idx][bl].nrow;
	//Int p_size = LL[team_leader][bl].p_size;
	Int p_size = LL[leader_idx][bl].p_size;
	//Int ws_size = LL[team_leader][bl].iws_size;
	Int ws_size = LL[leader_idx][bl].iws_size;
	Int *color = &(ws[0]);
	Int *pattern = &(color[ws_size]);


        #ifdef BASKER_DEBUG_NFACTOR_COL
	printf("copy, kid: %d bl: %d  A: %d %d \n", 
	       kid, bl, A_col, A_row);
        #endif
      
	
	//Bgood(remove)
	//Note: note where called
	//for(Int i = B.offset; i < B.m_offset; i++)
	Int bbcol = B.scol;
	//Int 
	for(Int i = B.col_ptr[k-bbcol]; i < B.col_ptr[k-bbcol+1]; i++)
	  {
	    Int B_row = B.row_idx(i);
	    Int j = gperm[B_row];

	    //Note: Do we need this anymore?
	    //This is being used, might want to check out why
	    //HERE!!!! (Check out G2)
	    //if(B_row < brow)
	    //  {
	    //	
	    //	printf("continue at : %d %d %d, kid: %d \n",
	    //	       B_row, B.srow, brow, kid);
	    //	continue;
	    //}
	    
	    //Note: Do we need this anymore?
	    //if(B_row > (brow+nrow))
	    //{
	    // printf("breaking at %d %d \n", B_row, j);
	    // break;
	    //}
	    
            #ifdef BASKER_DEBUG_NFACTOR_COL
	    printf("Scanning_2 A: %d %d lvl %d l: %d bl:%d brow: % d %d K: %d \n",
		 B_row, j, lvl, l, bl, brow, B.srow, kid);
	    printf("Adding Aval: %f to xval: %f \n", X[B_row-brow], B.val(i));
             #endif
	 
	    X[B_row-brow] += B.val(i);
	    if(color[B_row-brow] == 0)
	      {
		color[B_row-brow] = 1;
		pattern[p_size++] = B_row;
	      }

	  }//end for over all nnz
      }//end over all blks
        
    //-------------move into C 
    //(Right now have to do dense but like to do sparse)


    last_blk = LL_size[leader_idx];
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    for(Int bl=l+1; bl<last_blk; bl++)
      {
	Int A_col = S[lvl][kid];
	Int A_row = (lvl==1)?(2):S[bl][kid]%(LU_size[A_col]);
	Int CM_idx = kid;
	//ENTRY_1DARRAY   X = LL[team_leader][bl].ews;
	ENTRY_1DARRAY   X = LL[leader_idx][bl].ews;
	//INT_1DARRAY    ws = LL[team_leader][bl].iws;
	INT_1DARRAY    ws = LL[leader_idx][bl].iws;
	//Int        ws_size =LL[team_leader][bl].ews_size;
	Int        ws_size =LL[leader_idx][bl].ews_size;
	//Int brow = LL[team_leader][bl].srow;
	Int brow = LL[leader_idx][bl].srow;
	//Int nrow = LL[team_leader][bl].nrow;
	Int nrow = LL[leader_idx][bl].nrow;
	//Int p_size = LL[team_leader][bl].p_size;
	Int p_size = LL[leader_idx][bl].p_size;

	//For recounting patterns in dense blk
	//Need better sparse update
	Int p_count  =0 ; 
	
	Int *color = &(ws[0]);
	Int *pattern = &(color[ws_size]);

        #ifdef BASKER_DEBUG_NFACTOR_COL
	printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
	       kid, A_col, A_row, team_leader, bl,p_size);
        #endif
      
	//over all dim(S)
	for(Int jj=brow; jj < (brow+nrow); jj++)
	//for(Int jj = 0 ; jj < p_size; jj++)
	  {
	    //Int j = pattern[jj];
	    Int j = jj;
	    #ifdef BASKER_DEBUG_NFACTOR_COL
	    printf("considering: %d %d %f, kid: %d\n",
	       j,brow,X[j-brow], kid);
	    #endif

	    if(X[j-brow] != 0)
	      {
		if(bl == l+1)
		  {
		    #ifdef BASKER_DEBUG_NFACTOR_COL
		    printf("moving X[%d] %f kid: %d \n",
		     j, X[j-brow], kid);
		    #endif
		    C.row_idx[nnz] = j;
		    C.val[nnz] = X[j-brow];
		    nnz++;
		    X[j-brow] = 0;
		    color[j-brow]  = 0;
		  }
		else
		  {
		    #ifdef BASKER_DEBUG_NFACTOR_COL
		    printf("counting [%d] %f kid: %d \n",
		       j, X[j-brow], kid);
		    #endif
		    pattern[p_count++] = j;
		    color[j-brow] = 1;
		  }
	      }//if not empyt
	  }//over all dim(S)
	
	if(bl == l+1)
	  {
	    
	    #ifdef BASKER_DEBUG_NFACTOR_COL
	    printf("SETTING move_over set 0, L: %d %d kid: %d \n",
	       leader_idx, bl, kid);
	    #endif
	    LL[leader_idx][bl].p_size = 0;
	    p_count =0;
	  }
	else
	  {
	    #ifdef BASKER_DEBUG_NFACTOR_COL
	    printf("SETTING Re-pop pattern: %d %d size: %d \n",
	       leader_idx, bl, p_count);
	    #endif
	    LL[leader_idx][bl].p_size = p_count;
	  }

      }//over all blks

    //printf("kid: %d col_ptr: %d \n",
    //   kid, nnz);

    C.col_ptr[0] = 0;
    C.col_ptr[1] = nnz;

    //printf("Done with move, kid: %d found nnz: %d \n",
    //   kid, nnz);
    
    //C.info();
    //C.print();

    
    //set pointer 
    //Bgood(remove)
    /*
    for(Int bl = l+1; bl < l+2; bl++)
      {
	Int A_col = S[lvl][kid];
	Int A_row = (lvl==1)?(2):S[bl][kid]%(LU_size[A_col]);
	Int CM_idx = kid;
	BASKER_MATRIX_VIEW &B = AV[A_col][A_row];
	
	B.flip_base(&(thread_array[kid].C));
	B.k_offset = k;
    
	if(kid == 0)
	  {
	    //	    printf("\n\n----------------TEST, kid: %d-----------\n\n", kid);
	    // B.base->print();
	  }
      }
    */


    return 0;
  }//end t_copy_update_matrix()


}//end namespace BaskerNS

#endif //end def BASKER_NFACTOR_COL_INC_HPP
