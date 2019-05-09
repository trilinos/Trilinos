#ifndef SHYLUBASKER_SFACTOR_INC_HPP
#define SHYLUBASKER_SFACTOR_INC_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_util.hpp"

//#include "shylubasker_sfactor.hpp"

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::sfactor_inc()
  {

    //printf("-======= TEST CALL Sfactor_inc====== \n");

    if(btf_tabs_offset!=0)
    {
      sfactor_nd_estimate();
    }

    if(Options.btf == BASKER_TRUE)
    {
      if(btf_nblks > 1)
      {
        BASKER_ASSERT(0==1, "Should not get here");
      }
    }

    //Add init stuff
    
    //if(btf_tabs_offset != 0)
     {
     #ifdef BASKER_KOKKOS
     kokkos_sfactor_init_factor<Int,Entry,Exe_Space> iF(this);
     Kokkos::parallel_for(TeamPolicy(num_threads,1), iF);
     Kokkos::fence();
     #else
     #endif     
     }

     //if(btf_tabs_offset != 0)
     {
     //printf("before allow workspace\n");
     //Allocate workspace
     #ifdef BASKER_KOKKOS
     typedef Kokkos::TeamPolicy<Exe_Space>      TeamPolicy;
     kokkos_sfactor_init_workspace<Int,Entry,Exe_Space> iWS(this);
     Kokkos::parallel_for(TeamPolicy(num_threads,1), iWS);
     Kokkos::fence();
     #else
     
     #endif
     }
     
     
     BASKER_ASSERT(A.nrow > 0, "Sfactor A.nrow");
     MALLOC_INT_1DARRAY(gperm, A.nrow);
     init_value(gperm,A.nrow, BASKER_MAX_IDX);
     MALLOC_INT_1DARRAY(gpermi,A.nrow);
     init_value(gpermi, A.nrow, BASKER_MAX_IDX);
    
    
     //Incomplete Factor Setup
     if(Options.incomplete == BASKER_TRUE)
     {
     Int lvl_nnz = 1.2*global_nnz;
     MALLOC_INT_1DARRAY(INC_LVL_ARRAY, lvl_nnz);
     init_value(INC_LVL_ARRAY, lvl_nnz, BASKER_MAX_IDX);
     MALLOC_INT_1DARRAY(INC_LVL_ARRAY_CNT, A.nrow);
     init_value(INC_LVL_ARRAY_CNT, A.nrow,(Int)0);
     MALLOC_INT_1DARRAY(INC_LVL_TEMP, A.nrow);
     init_value(INC_LVL_TEMP, A.nrow, BASKER_MAX_IDX);
     }

     return 0;

  }//end sfactor_inc()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_estimate()
  {
    //We estmate only by using a fraction of the number of nnz
    
    for(Int p=0; p < num_threads; ++p)
    {
      Int blk = S(0)(p);
      sfactor_nd_dom_estimate(ALM(blk)(0), 
          LL(blk)(0), 
          LU(blk)(LU_size(blk)-1));

      for(Int l=0; l < tree.nlvls; l++)
      {
        Int U_col = S(l+1)(p);

        Int my_row_leader = find_leader(p,l);
        Int my_new_row   = 
          blk - S(0)(my_row_leader);

        Int U_row = (l==0)?(p%2):S(0)(p)%LU_size(U_col);
        if((blk > 14) &&
            (blk > LU_size(U_col)) &&
            (l!=0))
        {
          Int tm = (blk+1)/16;
          U_row = ((blk+1) -tm*16)%LU_size(U_col);
        }

        //JDB TEST PASSED
        U_row = my_new_row;

        sfactor_nd_upper_estimate(AVM(U_col)(U_row),
            LU(U_col)(U_row));

        sfactor_nd_lower_estimate(ALM(blk)(l+1),
            LL(blk)(l+1));

      } // end for l

      //printf("============SFACTOR INC SEP=======\n");
      for(Int lvl = 0; lvl < tree.nlvls; lvl++)
      {
        Int p = pow(tree.nparts, tree.nlvls-lvl-1);

        for(Int pp=0; pp < p; pp++)
        {
          Int ppp = pp*pow(tree.nparts, lvl+1);
          Int U_col = S(lvl+1)(ppp);
          Int U_row = 0;

          sfactor_nd_sep_estimate(ALM(U_col)(U_row),
              LL(U_col)(U_row),
              LU(U_col)(LU_size(U_col)-1));

          Int innerblk = U_col;
          for(Int l = lvl+1; l < tree.nlvls; l++)
          {
            U_col = S(l+1)(ppp);

            Int my_row_leader = find_leader(ppp,l);
            Int my_new_row = 
              S(lvl+1)(ppp) - S(0)(my_row_leader);

            U_row = S(lvl+1)(ppp)%LU_size(U_col);	 
            if((S(lvl+1)(ppp) > 14) &&
                (S(lvl+1)(ppp) > LU_size(U_col)) 
              )
            {
              Int tm = (S(lvl+1)(ppp)+1)/16;
              U_row = ((S(lvl+1)(ppp)+1) - 
                  (tm*16))%LU_size(U_col);
            }

            //JDB TEST PASS
            U_row = my_new_row;

            sfactor_nd_sep_upper_estimate(AVM(U_col)(U_row),
                LU(U_col)(U_row));

            sfactor_nd_sep_lower_estimate(
                ALM(innerblk)(l-lvl),
                LL(innerblk)(l-lvl));

          }//for - l
        }//for -p
      }//for - lvl
    }//for-over all threads
  }//end sfactor_nd_estimate


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::sfactor_nd_dom_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &LM,
   BASKER_MATRIX &UM
  )
  {
    //Get number of nnz in top/lower halfs
    Int nnz_top = 0; 
    Int nnz_low = 0;
    for(Int k = 0; k < M.ncol; k++)
    {
      for(Int j = M.col_ptr(k); j < M.col_ptr(k+1); ++j)
      {
        if(M.row_idx(j) >= k)
        {
          nnz_top++;
        }
        if(M.row_idx(j) <= k)
        {
          nnz_low++;
        }
      }//for -j = col_ptr
    }//for - k...ncol
 
    //NDE: compiler error in debug mode; no op* for int and complex<double>
    LM.nnz = nnz_low *
      (((double)BASKER_FILL_LESTIMATE+(double)Options.user_fill)* 
       (double)(Options.inc_lvl+1));
    global_nnz += LM.nnz;
    UM.nnz = nnz_top *
      (((double)BASKER_FILL_UESTIMATE+(double)Options.user_fill)*
       (double)(Options.inc_lvl+1));
    global_nnz += UM.nnz;
    
  }//end sfactor_nd_dom_estimate


  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_lower_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &LM
  )
  {
    LM.nnz = M.nnz *
      (((double)BASKER_FILL_LLOWERESTIMATE+(double)Options.user_fill)*
       (double)(Options.inc_lvl+1));
    
    global_nnz += LM.nnz;

  }//end sfactor_nd_lower_estimate()


  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_upper_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &UM
  )
  {
    UM.nnz = M.nnz * 
      (((double)BASKER_FILL_UUPPERESTIMATE+(double)Options.user_fill)*
       (double)(Options.inc_lvl+1));
    
    global_nnz += UM.nnz;

  }//end sfactor_nd_upper_estimate()

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_sep_upper_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &UM
  )
  {
    UM.nnz = M.ncol*M.nrow;
    global_nnz += UM.nnz;
  }//end sfactor_nd_sep_upper_estimate


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_sep_lower_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &LM
  )
  {
    LM.nnz = M.ncol*M.nrow;
    global_nnz += LM.nnz;
  }//end sfactor_nd_sep_lower_estimate


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::sfactor_nd_sep_estimate
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &LM,
   BASKER_MATRIX &UM
  )
  {
    //Get number of nnz in top/lower halfs
    Int nnz_top = 0; 
    Int nnz_low = 0;
    for(Int k = 0; k < M.ncol; k++)
    {
      for(Int j = M.col_ptr(k); j < M.col_ptr(k+1); ++j)
      {
        if(M.row_idx(j) >= k)
        {
          nnz_top++;
        }
        if(M.row_idx(j) <= k)
        {
          nnz_low++;
        }
      }//for -j = col_ptr
    }//for - k...ncol

    LM.nnz = nnz_top * 
      (((double)BASKER_FILL_LSEPESTIMATE+Options.user_fill)*
       (Options.inc_lvl+1));

    global_nnz += LM.nnz;

    UM.nnz = nnz_low * 
      (((double)BASKER_FILL_USEPESTIMATE+Options.user_fill)*
       (Options.inc_lvl+1));

    LM.nnz = M.nrow*M.nrow;
    UM.nnz = M.nrow*M.nrow;
    global_nnz += UM.nnz;

  }//end sfactor_nd_sep_estimate

}//end namespace BaskerNS

#endif
