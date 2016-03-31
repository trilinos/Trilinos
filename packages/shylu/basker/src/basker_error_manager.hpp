#ifndef BASKER_ERROR_MANAGER
#define BASKER_ERROR_MANAGER

/*Basker Includes*/
#include "basker_types.hpp"
#include "basker_util.hpp"
#include "basker_structs.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"


/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif


/*System Includes*/
#include <iostream>
#include <string>

namespace BaskerNS
{

  //===========DOMAIN ERROR HANDLE==========//
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::nfactor_domain_error
  (
   INT_1DARRAY threads_start
   )
  {
    Int nthread_remalloc = 0;
    for(Int ti = 0; ti < num_threads; ti++)
      {

	//Note: jdb we can make this into a switch
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOERROR)
	  {
	    threads_start(ti) = BASKER_MAX_IDX;
	    continue;
	  }//end if NOERROR

	if(thread_array(ti).error_type ==
	   BASKER_ERROR_SINGULAR)
	  {
	    printf("ERROR THREAD: %d DOMBLK SINGULAR: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if SINGULAR
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOMALLOC)
	  {
	    printf("ERROR THREADS: %d DOMBLK NOMALLOC: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if NOMALLOC
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_REMALLOC)
	  {
	    
	    BASKER_ASSERT(thread_array(ti).error_blk >= 0,
			  "nfactor_dom_error error_blk");

	    printf("ERROR THREADS: %d DOMBLK MALLOC: %d %d \n",
		   ti,
		   thread_array(ti).error_blk, 
		   thread_array(ti).error_subblk);
	    

	    //If on diagonal, want to compare L and U
	    Int resize_L = BASKER_MAX_IDX;   
	    Int resize_U = BASKER_MAX_IDX;
	    if(thread_array(ti).error_subblk != BASKER_MAX_IDX)
	      {
		
		BASKER_ASSERT(thread_array(ti).error_info >0,
			      "L) newsize not big enough");
		resize_L = thread_array(ti).error_info;
		printf("L(%d)  resize: %d \n", ti,  resize_L);


		//if L is already bigger and U, 
		//We will want re size U as, well
		if(thread_array(ti).error_subblk == 0)
		  {
		    Int blkcol = thread_array(ti).error_blk;
		    Int blkUrow = LU_size(blkcol)-1;
		    if(LL(blkcol)(0).nnz >=
		       LU(blkcol)(blkUrow).nnz)
		      {
			printf("resize U too \n");
			resize_U = thread_array(ti).error_info;
		      }
		  }//if - a domain
	      }
	    //We don't care about the other way since,
	    //L is already checked before U.
	    printf("now to check U\n");
	    if(thread_array(ti).error_subblk == -1)
	      {  
		resize_U = thread_array(ti).error_info;
		printf("resize U: %d \n", resize_U);
	      }

	    printf("before resize \n");
	    //Resize L
	    if(resize_L > BASKER_MAX_IDX)
	      {
		BASKER_MATRIX &L =
   LL(thread_array(ti).error_blk)(thread_array(ti).error_subblk);
		REALLOC_INT_1DARRAY(L.row_idx,
				    L.nnz,
				    resize_L);
		REALLOC_ENTRY_1DARRAY(L.val,
				      L.nnz,
				      resize_L);
		if(Options.incomplete == BASKER_TRUE)
		  {
		    REALLOC_INT_1DARRAY(L.inc_lvl,
					L.nnz,
					resize_L);

		  }
		L.nnz = resize_L;
	      }

	    //Resize U
	    if(resize_U > BASKER_MAX_IDX)
	      {
		BASKER_MATRIX &U = 
		  LU(thread_array(ti).error_blk)(0);
		REALLOC_INT_1DARRAY(U.row_idx,
			       U.nnz,
			       resize_U);
		REALLOC_ENTRY_1DARRAY(U.val,
				      U.nnz,
				      resize_U);
		U.nnz = resize_U;
  
	      }

	    //clean up workspace
	    //if(LL(thread_array(ti).error_blk)(0).w_fill == 
	    //  BASKER_TRUE)
	      {
		//Clear workspace, whole column
		for(Int sb = 0; 
		    sb < LL_size(thread_array(ti).error_blk);
		    sb++)
		  {
		    BASKER_MATRIX &SL = 
		      LL(thread_array(ti).error_blk)(sb);
		    for(Int i = 0; 
			i < SL.iws_size*SL.iws_mult; 
			++i)
		      {
			SL.iws(i) = (Int) 0;
		      }
		    for(Int i = 0;
			i < SL.ews_size*SL.ews_mult; 
			++i)
		      {
			SL.ews(i) = (Entry) 0;
		      }
		    if(sb == 0)
		      {
			printf("CLEARING PERM\n");
			//Clear perm
			for(Int i = SL.srow; 
			    i < SL.srow+SL.nrow; ++i)
			  {
			    gperm(i) = BASKER_MAX_IDX;
			  }

			//Clear incomplete ws
			if(Options.incomplete == BASKER_TRUE)
			  {
			    for(Int i = SL.srow;
				i < SL.srow+SL.nrow; ++i)
			      {
				INC_LVL_TEMP(i) = BASKER_MAX_IDX;
			      }
			  }
		      }
		  }//for - sb (subblks)
	      }//if ws is filled
	     
	    threads_start(ti) = thread_array(ti).error_blk;
	

	    //Reset 
	    thread_array(ti).error_type = BASKER_ERROR_NOERROR;
	    thread_array(ti).error_blk  = BASKER_MAX_IDX;
	    thread_array(ti).error_info = BASKER_MAX_IDX;

	    nthread_remalloc++;
    
	  }//if REMALLOC

      }//for all threads
    
    if(nthread_remalloc == 0)
      {
	return BASKER_SUCCESS;
      }
    else
      {
	return nthread_remalloc;
      }

    //Should never be here
    BASKER_ASSERT(0==1, "nfactor_diag_error, should never");
    return BASKER_SUCCESS;
  }//end nfactor_domain_error

 

  //========SEP ERROR HANDLE===========//
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::nfactor_sep_error
  (
   INT_1DARRAY thread_start
   )
  {
    Int nthread_remalloc = 0;
    for(Int ti = 0; ti < num_threads; ti++)
      {

	//Note: jdb we can make this into a switch
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOERROR)
	  {
	    thread_start(ti) = BASKER_MAX_IDX;
	    continue;
	  }//end if NOERROR

	if(thread_array(ti).error_type ==
	   BASKER_ERROR_SINGULAR)
	  {
	    printf("ERROR THREAD: %d DOMBLK SINGULAR: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if SINGULAR
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOMALLOC)
	  {
	    printf("ERROR THREADS: %d DOMBLK NOMALLOC: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if NOMALLOC
	
	
	//Find lvl in sep error happend
	Int error_sep_lvl = BASKER_MAX_IDX;
	for(Int l = 1; l < tree.nlvls; l++)
	  {
	    if(thread_array(ti).error_blk == 
	       S(l)(ti))
	      {
		error_sep_lvl = l;
		break;
	      }
	  }


	if(thread_array(ti).error_type ==
	   BASKER_ERROR_REMALLOC)
	  {

	    BASKER_ASSERT(thread_array(ti).error_blk > 0,
			  "nfactor_SEP_error error_blk");

	    printf("ERROR THREADS: %d SEPBLK MALLOC: %d %d \n",
		   ti,
		   thread_array(ti).error_blk, 
		   thread_array(ti).error_subblk);

	    printf("ERROR SEPLVL: %d \n",
		   error_sep_lvl);
	    
	    //If on diagonal, want to compare L and U
	    Int resize_L = BASKER_MAX_IDX;   
	    Int resize_U = BASKER_MAX_IDX;
	    if(thread_array(ti).error_subblk <= -1)
	      {
		resize_L = thread_array(ti).error_info;    
		printf("L size: %d \n", resize_L);
	      }
	    //We don't care about the other way since,
	    //L is already checked before U.
	    if(thread_array(ti).error_subblk > -1)
	      {
		resize_U = thread_array(ti).error_info;
		printf("U size: %d \n", resize_U);
	      }

	    //Resize L
	    if(resize_L > BASKER_MAX_IDX)
	      {
		//printf("Resizing L\n");
		const Int tsb = (-1*thread_array(ti).error_subblk)-1;
		BASKER_MATRIX &L =
   LL(thread_array(ti).error_blk)(tsb);
		REALLOC_INT_1DARRAY(L.row_idx,
				    L.nnz,
				    resize_L);
		REALLOC_ENTRY_1DARRAY(L.val,
				      L.nnz,
				  resize_L);
		L.nnz = resize_L;
	      }

	    //Resize U
	    if(resize_U > BASKER_MAX_IDX)
	      {
		//printf("Resizing U\n");
		const Int tsb = thread_array(ti).error_subblk;
		BASKER_MATRIX &U = 
		  LU(thread_array(ti).error_blk)(tsb);
		REALLOC_INT_1DARRAY(U.row_idx,
			       U.nnz,
			       resize_U);
		REALLOC_ENTRY_1DARRAY(U.val,
				      U.nnz,
				      resize_U);
		//printf("done resizing normal\n");
		if(Options.incomplete == BASKER_TRUE)
		  {
		    REALLOC_INT_1DARRAY(U.inc_lvl,
					U.nnz,
					resize_U);
		  }
		//printf("done resizing incomplete");
		U.nnz = resize_U;
	      }
	    
	    //printf("done resize matrix\n");
	    //clean up workspace
	    //No nice way to do this since multiple threads
	    //Though this could be done in parallel in the future
	    for(Int p = 0; p < num_threads; p++)
	      {
		Int blk = S(0)(p);
		//printf("clear blk: %d \n", blk);
		//if(LL(blk)(0).w_fill == BASKER_TRUE)
		  {
		    //Clear workspace, whole column
		    for(Int sb = 0; 
			sb < LL_size(blk);
			sb++)
		      {
			BASKER_MATRIX &SL =  LL(blk)(sb);
			for(Int i = 0; 
			    i < SL.iws_size*SL.iws_mult; 
			    ++i)
			  {
			    SL.iws(i) = (Int) 0;
			  }
			for(Int i = 0;
			    i < SL.ews_size*SL.ews_mult; 
			    ++i)
			  {
			    SL.ews(i) = (Entry) 0;
			  }
		      }//for - sb (subblks)
		  }//if ws is filled
	      }//for-other all threads
	     

	    //Clear perm
	      for(Int p = 0; p < num_threads; p++)
	      {
		Int blk = S(error_sep_lvl)(p);
		//printf("perm clear blk: %d \n", blk);
		//if(LL(blk)(0).w_fill == BASKER_TRUE)
		{
		  BASKER_MATRIX &TM = LL(blk)(0);
		  for(Int i = TM.scol; i < TM.scol+TM.ncol; i++)
		    {
		      gperm(i) = BASKER_MAX_IDX;
		    }
		  
		}//if ws is filled
	      }//for-other all threads
	  

	    //printf("done resize workspace\n");
	    //Note, will have to clear the perm in all sep blk in that level
	    //Clear permuation
	    BASKER_MATRIX &SL = 
	      LL(thread_array(ti).error_blk)(0);
	    
	    for(Int i = SL.srow; i < (SL.srow+SL.nrow);
		i++)
	      {
		gperm(i) = BASKER_MAX_IDX;
	      }//for--to clear perm
	
	
	    thread_start(ti) = thread_array(ti).error_blk;
	
	    //Reset 
	    thread_array(ti).error_type = BASKER_ERROR_NOERROR;
	    thread_array(ti).error_blk  = BASKER_MAX_IDX;
	    thread_array(ti).error_info = BASKER_MAX_IDX;
	    
	    for(Int i = 0; i < num_threads; i++)
	      {
		basker_barrier.ExitSet(i,BASKER_FALSE);
	      }
	    
	    nthread_remalloc++;
    
	  }//if REMALLOC

	//Reset Inc vector 
	if(Options.inc_lvl == BASKER_TRUE)
	  {
	    for(Int i = 0; i < INC_LVL_TEMP.dimension_0(); i++)
	      {
		INC_LVL_TEMP(i) = BASKER_MAX_IDX;
	      }

	  }


      }//for all threads
    
    if(nthread_remalloc == 0)
      {
	return BASKER_SUCCESS;
      }
    else
      {
	return nthread_remalloc;
      }

    //Should never be here
    BASKER_ASSERT(0==1, "nfactor_diag_error, should never");
    return BASKER_SUCCESS;

  }//end nfactor_sep_error


  //========BTF ERROR HANDLE==============//
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::nfactor_diag_error
  (
   INT_1DARRAY threads_start
   )
  {
    Int nthread_remalloc = 0;
    for(Int ti = 0; ti < num_threads; ti++)
      {

	//Note: jdb we can make this into a switch
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOERROR)
	  {
	    threads_start(ti) = BASKER_MAX_IDX;
	    continue;
	  }//end if NOERROR

	if(thread_array(ti).error_type ==
	   BASKER_ERROR_SINGULAR)
	  {
	    printf("ERROR THREAD: %d DIAGBLK SINGULAR: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if SINGULAR
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOMALLOC)
	  {
	    printf("ERROR THREADS: %d DIAGBLK NOMALLOC: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if NOMALLOC
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_REMALLOC)
	  {

	    BASKER_ASSERT(thread_array(ti).error_blk > 0,
			  "nfactor_diag_error error_blk");

	    printf("ERROR THREADS: %d DIAGBLK MALLOC: %d \n",
		   ti,
		   thread_array(ti).error_blk);
	    
	    	    //Clean the workspace
	    printf("test: %d %d \n",
		   thread_array(ti).iws_size*thread_array(ti).iws_mult,
		   thread_array(ti).ews_size*thread_array(ti).ews_mult);

	    for(Int i = 0; 
		i < thread_array(ti).iws_size*thread_array(ti).iws_mult;
		i++)
	      {
		thread_array(ti).iws(i) = (Int) 0;
	      }
	    for(Int i = 0;
		i < thread_array(ti).ews_size*thread_array(ti).ews_mult;
		i++)
	      {
		thread_array(ti).ews(i) = (Entry) 0;
	      }

	   

	    //Resize L
	    BASKER_MATRIX &L = LBTF(thread_array(ti).error_blk);
	    REALLOC_INT_1DARRAY(L.row_idx,
			       L.nnz,
			       thread_array(ti).error_info);
	    REALLOC_ENTRY_1DARRAY(L.val,
				 L.nnz,
				 thread_array(ti).error_info);
	    L.nnz = thread_array(ti).error_info;
	    for(Int i = 0; i < L.ncol; i++)
	      {
		L.col_ptr(i) = 0;
	      }

	    for(Int i = L.srow; i < (L.srow+L.nrow); i++)
	      {
		gperm(i) = BASKER_MAX_IDX;
	      }

	    //Resize U
	    BASKER_MATRIX &U = UBTF(thread_array(ti).error_blk);
	    REALLOC_INT_1DARRAY(U.row_idx,
			       U.nnz,
			       thread_array(ti).error_info);
	    REALLOC_ENTRY_1DARRAY(U.val,
				 U.nnz,
				 thread_array(ti).error_info);
	    U.nnz = thread_array(ti).error_info;
	    for(Int i = 0; i < U.ncol; i++)
	      {
		U.col_ptr(i) = 0;
	      }

	    
	    printf("Setting thread start(%d) %d \n",
		   ti, thread_array(ti).error_blk);

	    threads_start(ti) = thread_array(ti).error_blk;
	

	    //Reset 
	    thread_array(ti).error_type = BASKER_ERROR_NOERROR;
	    thread_array(ti).error_blk  = BASKER_MAX_IDX;
	    thread_array(ti).error_info = BASKER_MAX_IDX;

	    nthread_remalloc++;
    
	  }//if REMALLOC

      }//for all threads
    
    if(nthread_remalloc == 0)
      {
	return BASKER_SUCCESS;
      }
    else
      {
	return nthread_remalloc;
      }

    //Should never be here
    BASKER_ASSERT(0==1, "nfactor_diag_error, should never");
    return BASKER_SUCCESS;
  }//end nfactor_diag_error




}//end namespace BaskerNS

#endif //END BASER_ERROR_MANAGER
