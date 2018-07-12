#ifndef SHYLUBASKER_ERROR_MANAGER
#define SHYLUBASKER_ERROR_MANAGER

/*Basker Includes*/
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"


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
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR THREAD: " << ti 
            << " DOMBLK SINGULAR: " << thread_array(ti).error_blk
            << std::endl;
        }
        return BASKER_ERROR;
      }//end if SINGULAR

      if(thread_array(ti).error_type ==
          BASKER_ERROR_NOMALLOC)
      {
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR THREAD: " << ti 
            << " DOMBLK NOMALLOC : " << thread_array(ti).error_blk
            << std::endl;
        }
        return BASKER_ERROR;
      }//end if NOMALLOC

      if(thread_array(ti).error_type ==
          BASKER_ERROR_REMALLOC)
      {

        BASKER_ASSERT(thread_array(ti).error_blk >= 0,
            "nfactor_dom_error error_blk");
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR THREAD: " << ti 
            << " DOMBLK MALLOC : " << thread_array(ti).error_blk
            << " " << thread_array(ti).error_subblk
            << std::endl;
        }

        //If on diagonal, want to compare L and U
        Int resize_L = BASKER_MAX_IDX;   
        Int resize_U = BASKER_MAX_IDX;
        if(thread_array(ti).error_subblk != BASKER_MAX_IDX)
        {

          BASKER_ASSERT(thread_array(ti).error_info >0, "L) newsize not big enough");
          resize_L = thread_array(ti).error_info;
          if(Options.verbose == BASKER_TRUE)
          {
            std::cout << "L(" << ti << ")  resize: " << resize_L << std::endl;
          }

          //if L is already bigger and U, 
          //We will want re size U as, well
          if(thread_array(ti).error_subblk == 0)
          {
            Int blkcol = thread_array(ti).error_blk;
            Int blkUrow = LU_size(blkcol)-1;
            if(LL(blkcol)(0).nnz >=
                LU(blkcol)(blkUrow).nnz)
            {
              resize_U = thread_array(ti).error_info;
            }
          }//if - a domain
        }
        //We don't care about the other way since,
        //L is already checked before U.
        if(thread_array(ti).error_subblk == -1)
        {  
          resize_U = thread_array(ti).error_info;
          if(Options.verbose == BASKER_TRUE)
          {
            std::cout << "resize U: " << resize_U << std::endl;
          }
        }

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
          L.clear_pend();
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
          //Still need to clear pend
          BASKER_MATRIX &L = 
            LL(thread_array(ti).error_blk)(0);
          L.clear_pend();

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
      return BASKER_SUCCESS;
    else
      return nthread_remalloc;

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
      if(thread_array(ti).error_type == BASKER_ERROR_NOERROR)
      {
        thread_start(ti) = BASKER_MAX_IDX;
        continue;
      }//end if NOERROR

      if(thread_array(ti).error_type == BASKER_ERROR_SINGULAR)
      {
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR THREAD: " << ti 
            << " DOMBLK SINGULAR: " << thread_array(ti).error_blk
            << std::endl;
        }
        return BASKER_ERROR;
      }//end if SINGULAR

      if(thread_array(ti).error_type == BASKER_ERROR_NOMALLOC)
      {
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR THREADS: " << ti 
            << " DOMBLK NOMALLOC: " << thread_array(ti).error_blk
            << std::endl;
        }
        return BASKER_ERROR;
      }//end if NOMALLOC


      //Find lvl in sep error happend
      Int error_sep_lvl = BASKER_MAX_IDX;
      for(Int l = 1; l < tree.nlvls+1; l++)
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

        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR_THREADS: " << ti
            << " SEPBLK MALLOC: " << thread_array(ti).error_blk 
            << " " << thread_array(ti).error_subblk
            << std::endl;

          std::cout << "ERROR SEPLVL: " << error_sep_lvl << std::endl;
        }

        //If on diagonal, want to compare L and U
        Int resize_L = BASKER_MAX_IDX;   
        Int resize_U = BASKER_MAX_IDX;
        if(thread_array(ti).error_subblk <= -1)
        {
          resize_L = thread_array(ti).error_info;    
          if(Options.verbose == BASKER_TRUE)
          {
            std::cout << "L size: " << resize_L << std::endl;
          }
        }
        //We don't care about the other way since,
        //L is already checked before U.
        if(thread_array(ti).error_subblk > -1)
        {
          resize_U = thread_array(ti).error_info;
          if(Options.verbose == BASKER_TRUE)
          {
            std::cout << "U size: " << resize_U << std::endl;
          }
        }

        //Resize L
        if(resize_L > BASKER_MAX_IDX)
        {
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
          L.clear_pend();
        }

        //Resize U
        if(resize_U > BASKER_MAX_IDX)
        {
          const Int tsb = thread_array(ti).error_subblk;
          BASKER_MATRIX &U = 
            LU(thread_array(ti).error_blk)(tsb);
          REALLOC_INT_1DARRAY(U.row_idx,
              U.nnz,
              resize_U);
          REALLOC_ENTRY_1DARRAY(U.val,
              U.nnz,
              resize_U);
          if(Options.incomplete == BASKER_TRUE)
          {
            REALLOC_INT_1DARRAY(U.inc_lvl,
                U.nnz,
                resize_U);
          }
          U.nnz = resize_U;
        }

        //clean up workspace
        //No nice way to do this since multiple threads
        //Though this could be done in parallel in the future
        for(Int p = 0; p < num_threads; p++)
        {
          Int blk = S(0)(p);
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
          //if(LL(blk)(0).w_fill == BASKER_TRUE)
          {
            BASKER_MATRIX &TM = LL(blk)(0);
            for(Int i = TM.scol; i < TM.scol+TM.ncol; i++)
            {
              gperm(i) = BASKER_MAX_IDX;
            }

          }//if ws is filled
        }//for-other all threads


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
        //for(Int i = 0; i < INC_LVL_TEMP.extent(0); i++) //NDE - warning: comparison s and us
        for(Int i = 0; i < static_cast<Int>( INC_LVL_TEMP.extent(0) ); i++) 
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
    BASKER_ASSERT(0==1, "nfactor_sep_error, should never");
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
        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR_THREAD: " << ti
            << " DIAGBLK SINGULAR " << thread_array(ti).error_blk
            << std::endl;
        }
        return BASKER_ERROR;
      }//end if SINGULAR

      if(thread_array(ti).error_type ==
          BASKER_ERROR_NOMALLOC)
      {
        std::cout << "ERROR_THREADS: " << ti
          << " DIAGBLK NOMALLOC " << thread_array(ti).error_blk
          << std::endl;
        return BASKER_ERROR;
      }//end if NOMALLOC

      if(thread_array(ti).error_type ==
          BASKER_ERROR_REMALLOC)
      {

        BASKER_ASSERT(thread_array(ti).error_blk > 0,
            "nfactor_diag_error error_blk");

        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "ERROR_THREADS: " << ti
            << " DIAGBLK MALLOC " << thread_array(ti).error_blk
            << std::endl;

          //Clean the workspace
          std::cout << "test: "
            << thread_array(ti).iws_size*thread_array(ti).iws_mult
            << " " << thread_array(ti).ews_size*thread_array(ti).ews_mult
            << std::endl;
        }

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
        L.clear_pend();
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


        if(Options.verbose == BASKER_TRUE)
        {
          std::cout << "Setting thread start( " << ti
            << ")  " << thread_array(ti).error_blk
            << std::endl;
        }

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
