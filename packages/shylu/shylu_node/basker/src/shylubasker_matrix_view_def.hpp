#ifndef SHYLUBASKER_MATRIX_VIEW_DEF_HPP
#define SHYLUBASKER_MATRIX_VIEW_DEF_HPP

/*Basker Includes*/
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_types.hpp"
//#include "shylubasker_decl.hpp"

/*System Includes*/
#include <iostream>
#include <stdio.h>

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#else
#include <omp.h>
#endif

using namespace std;

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrixView<Int,Entry,Exe_Space>::BaskerMatrixView()
  {
    m_offset = 0;
    k_offset=0;
    offset = 0;
    roll_back = BASKER_FALSE;
  }//end BaskerMatrixView()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrixView<Int,Entry,Exe_Space>::BaskerMatrixView
  (
   BASKER_MATRIX *_base,
	 Int _sr, Int _m, Int _sc, Int _n
  )
  {
    init(_base, _sr, _m, _sc, _n);
  }//end BaskerMatrixView()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::init
  (
   BASKER_MATRIX *_base, Int _sr, Int _m, Int _sc, Int _n
  )
  {
    //printf("Matrix View Init.  sr: %d  m: %d  sc: %d  n: %d \n",
    //	   _sr, _m, _sc, _n);

    base = _base;
    scol = _sc;
    srow = _sr;
    ncol = _n;
    nrow = _m;
    k_offset = 0;
    //m_offset = _base->col_ptr[_base->ncol+1];
    m_offset = _base->nnz;
    offset = 0;
    perm = BASKER_FALSE;
    roll_back = BASKER_FALSE;
  }//end init()


  /*
  //To be finished
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::convert_to_2D(BASKER_MATRIX &M)
  {
    //Asummption is that M is a empty BASKER MATRIX
    //Will convert the current view into a standalone matrix

    INT_1DARRAY col_counts;
    MALLOC_1DARRAY(col_counts, ncol);
    init_value(col_counts, ncol, 0);

    Int t_nnz = 0;
    for(Int i = scol;  i < scol+ncol; i++)
    {
      for(Int i = base->col_ptr[i - k_offset];
          i < base->col_ptr[i+1-k_offset];
          i++)
      {
        j = base->row_idx[i];
        if(j > srow+nrow)
        {
          break;
        }
        if(j >= srow)
        {
          col_count[i-scol] = col_count[i-scol]+1;;
        }
      }//over all row
    }//over all col

    //sum up 
  }
  */


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::init_perm(INT_1DARRAY *pinv)
  {
    perm = BASKER_TRUE;
    lpinv = pinv;
  }//end init_perm()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::init_offset(Int k, Int prev)
  {
    if(!perm)
    {
      offset = 0;

      if(prev == 0)
      {
        prev = base->col_ptr[k - k_offset];
      }

      m_offset = base->col_ptr[k+1 -k_offset];

      while((prev < base->col_ptr[k+1 - k_offset])&&(base->row_idx[prev] < srow))
      {
        prev++;
      }

      offset = prev;
    }
    else
    {
      offset = base->col_ptr[k-k_offset];
      m_offset = base->col_ptr[k+1 -k_offset];
      // can't use a linear ordering of them because they 
      //may be anywhere!!!!!
    }
  }//end init_offset()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int BaskerMatrixView<Int,Entry,Exe_Space>::good(Int i)
  {
    if(!perm)
    {
      //printf("not permed %d \n", i);
      //printf("good: i: %d m_offset: %d row_idx: %d srow: %d nrow: %d \n", i, m_offset, base->row_idx[i], srow, nrow);
      if(
          (base->row_idx[i] < (srow+nrow)) &&
          (base->row_idx[i] >= srow))
        return 0;
      else
        return BASKER_MAX_IDX;
    }
    else
    {
      //printf("good: j: %d t: %d srow: %d\n",
      //     base->row_idx[i], 
      //   (*lpinv)[base->row_idx[i]], 
      //	 srow);
      if((i < m_offset) && 
          ((*lpinv)[base->row_idx[i]] < (srow+nrow)) &&
          ((*lpinv)[base->row_idx[i]] >= srow))
      {
        return 0;
      }
      else
      {
        return BASKER_MAX_IDX;
      }
    }
  }//end good()


  //Note: function defunct
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int BaskerMatrixView<Int,Entry,Exe_Space>::good_extend(Int i)
  {
    if(i < m_offset)
      return 0;
    else
      return BASKER_MAX_IDX;
  }//end good_extend
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int BaskerMatrixView<Int,Entry,Exe_Space>::col_ptr(Int k)
  {
    return base->col_ptr[k-k_offset];
  }//end col_ptr()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int BaskerMatrixView<Int,Entry,Exe_Space>::row_idx(Int i)
  {
    return base->row_idx[i];
  }//end row_idx()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Entry BaskerMatrixView<Int,Entry,Exe_Space>::val(Int i)
  {
    return base->val[i];
  }//end val()


  //This needs to be made faster or not done in code
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int BaskerMatrixView<Int,Entry,Exe_Space>::nnz()
  {
    Int _nnz = 0;

    for(Int k=scol; k<scol+ncol; k++)
    {
      for(Int j = col_ptr(k); j < col_ptr(k+1); j++)
      {
        if(good(j) == 0)
        {
          _nnz++;
        }
      }
    }//over all columns

    return _nnz;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::flip_base()
  {
    if(roll_back == BASKER_TRUE)
    {
      #ifdef BASKER_DEBUG
      printf("Flipping Basker for old one\n");
      #endif
      BASKER_MATRIX *temp = base;
      base = base_backup;
      base_backup = temp;
      offset = 0;
      k_offset = 0;
    }
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::flip_base(
                                            BASKER_MATRIX *_base)
  {
    #ifdef BASKER_DEBUG
    printf("Fliping Base for new one\n");
    #endif
    roll_back = BASKER_TRUE;
    base_backup = base;
    base = _base;
    offset = 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::info()
  {
    print();
  }//end info()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrixView<Int,Entry,Exe_Space>::print()
  {
      printf("\n scol: %d   srow: %d  \n ncol: %d nrow: %d \n offset: %d k_offset: %d \n", 
             scol, srow, ncol, nrow, offset, k_offset);
  }//end print()

}//end namespace Basker

#endif //end of ifndef basker_matrix_view_def
