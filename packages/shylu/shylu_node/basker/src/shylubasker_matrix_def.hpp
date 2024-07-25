// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_MATRIX_DEF_HPP
#define SHYLUBASKER_MATRIX_DEF_HPP

#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif
#include "Teuchos_ScalarTraits.hpp"


#include <iostream>
//using namespace std;

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry,Exe_Space>::BaskerMatrix()
  {
    label   = "default called";
    ncol    = 0;
    nrow    = 0;
    nnz     = 0;
    mnnz    = 0;
    v_fill  = BASKER_FALSE;
    tpivot  = 0;
    #ifdef BASKER_2DL
    p_size  = 0;
    w_fill  = BASKER_FALSE;
    #endif
    inc_lvl_flg = BASKER_FALSE;
    anorm = -1.0;
    gnorm = -1.0;
    //printf("matrix init\n");
  }//end BaskerMatrix()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry,Exe_Space>::BaskerMatrix(string _label)
  {
    label   = _label;
    ncol    = 0;
    nrow    = 0;
    nnz     = 0;
    v_fill  = BASKER_FALSE;
    tpivot  = 0;
    #ifdef BASKER_2DL
    p_size  = 0;
    w_fill  = BASKER_FALSE;
    #endif
    inc_lvl_flg = BASKER_FALSE;
    anorm = -1.0;
    gnorm = -1.0;
  }//end BaskerMatrix(label)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry, Exe_Space>::BaskerMatrix
  (
   Int _m, Int _n, Int _nnz,
   Int *colptr, Int *rowind, Entry *nzval
  )
  {
    tpivot = 0;
    anorm = -1.0;
    gnorm = -1.0;
    #ifdef BASKER_2DL
    p_size = 0;
    w_fill = BASKER_FALSE;
    #endif

    v_fill = BASKER_FALSE;
    init_matrix("matrix", _m,_n,_nnz, colptr, rowind, nzval);
  }//end BaskerMatrix(int, int, int, int*, int*, Entry*)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry, Exe_Space>::BaskerMatrix
  (
   string _label, Int _m, Int _n, Int _nnz,
   Int *colptr, Int *rowind, Entry *nzval
  )
  {
    tpivot = 0;
    anorm = -1.0;
    gnorm = -1.0;
    #ifdef BASKER_2DL
    p_size = 0;
    w_fill = BASKER_FALSE;
    #endif

    v_fill = BASKER_FALSE;
    init_matrix(_label, _m,_n,_nnz, colptr, rowind, nzval);
  }//end BaskerMatrix(int, int, int, int*, int*, Entry*)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry,Exe_Space>::~BaskerMatrix()
  {
    Finalize();
  }//end ~BaskerMatrix

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::Finalize()
  {
    if(v_fill == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(col_ptr);
      FREE_INT_1DARRAY(row_idx);
      FREE_ENTRY_1DARRAY(val);
      v_fill = BASKER_FALSE;
    }

    if(w_fill == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(iws);
      FREE_ENTRY_1DARRAY(ews);
      w_fill = BASKER_FALSE;
    }

    if(inc_lvl_flg == BASKER_TRUE)
    {
      FREE_INT_1DARRAY(inc_lvl);
      inc_lvl_flg = BASKER_FALSE;
    }
  }//end finalize()

  template <class Int, class Entry, class Exe_Space>
  void BaskerMatrix<Int,Entry,Exe_Space>::set_shape
  (
   Int _sr, Int _m, Int _sc, Int _n
  )
  {
    srow = _sr;
    nrow = _m;
    scol = _sc;
    ncol = _n;
    nnz  = 0;
    mnnz = 0;
  }//end set_shape()

  template <class Int, class Entry, class Exe_Space>
  void BaskerMatrix<Int, Entry, Exe_Space>::copy_vec
  (
   Int* ptr, Int size, INT_1DARRAY a
  )
  {
    //Does Kokkos have a copy function
    for(Int i=0; i < size; i++)
    {
      a(i) = ptr[i];
    }

  }//end copy_vec(Int*, 1d)

  template <class Int, class Entry, class Exe_Space>
  void BaskerMatrix<Int, Entry, Exe_Space>::copy_vec
  (
   Entry* ptr, Int size, ENTRY_1DARRAY a
  )
  {
    //Kokkos::deep_copy(a.data, ptr);
    for(Int i=0; i < size; i++)
    {
      a(i) = ptr[i];
    }
    //return 0;
  }//edn copy_vec(Entry*, 1d);

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int, Entry, Exe_Space>::init_col()
  {
    //printf( " init_col(n=%d)\n",ncol );
    BASKER_ASSERT(ncol >= 0, "INIT_COL, ncol > 0");
    MALLOC_INT_1DARRAY(col_ptr, ncol+1);
    MALLOC_INT_1DARRAY(col_idx, ncol+1);
    for(Int i = 0; i < ncol+1; ++i)
    {
      col_ptr(i) = (Int) BASKER_MAX_IDX;
      col_idx(i) = (Int) BASKER_MAX_IDX;
    }
  }//end init_col()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::clean_col()
  {
    //printf( " clean_col(ncol = %d)\n",ncol );
    for(Int i = 0; i < ncol+1; ++i)
    {
      col_ptr(i) = (Int) BASKER_MAX_IDX;
      col_idx(i) = (Int) BASKER_MAX_IDX;
    }
    nnz = 0;
    mnnz = 0;
  }//end clean_col()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int, Entry, Exe_Space>::init_vectors
  (
   Int _m, Int _n, Int _nnz
  )
  {
    nrow = _m;
    ncol = _n;
    nnz  = _nnz;
    mnnz = _nnz;

    if(ncol >= 0)
    {
      BASKER_ASSERT((ncol+1)>0, "matrix init_vector ncol");
      MALLOC_INT_1DARRAY(col_ptr,ncol+1);
    }
    if(nnz > 0)
    {
      BASKER_ASSERT(nnz > 0, "matrix init_vector nnz");
      MALLOC_INT_1DARRAY(row_idx,nnz);
      MALLOC_ENTRY_1DARRAY(val,nnz);
#ifdef BASKER_INC_LVL
      MALLOC_INT_1DARRAY(inc_lvl, nnz);
#endif
    }
    else if(nnz==0)
    {
      BASKER_ASSERT((nnz+1)>0, "nnz+1 init_vector ");
      MALLOC_INT_1DARRAY(row_idx, nnz+1);
      row_idx(0) = (Int) 0;
      MALLOC_ENTRY_1DARRAY(val, nnz+1);
      val(0) = (Entry) 0;
#ifdef BASKER_INC_LVL
      MALLOC_INT_1DARRAY(inc_lvl, nnz+1);
#endif
    }

  }//end init_vectors(Int, Int, Int)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int, Entry, Exe_Space>::init_matrix
  (
   string _label, Int _m, Int _n, 
   Int _nnz
  )
  {
    label = _label;
    init_vectors(_m, _n, _nnz);
  }//end init_matrix(string, Int, Int, Int)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::init_matrix
  (
   string _label, Int _m, Int _n, Int _nnz,
   Int *_col_ptr, Int *_row_idx, Entry *_val
  )
  {
    label = _label;
    init_vectors(_m, _n, _nnz);
    copy_vec(_col_ptr,  _n+1,  col_ptr);
    copy_vec(_row_idx,  _nnz,  row_idx);
    copy_vec(_val,      _nnz,  val);
  }//end init_matrix(string with ptrs)


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int, Entry, Exe_Space>::init_matrix
  (
   string _label,
   Int _sr, Int _m,
   Int _sc, Int _n,
   Int _nnz
  )
  {
    srow = _sr;
    erow = srow + _m;
    scol = _sc;
    ecol = scol + _n;

    init_vectors(_m, _n, _nnz);
  }//end init_matrix(Int, Int, Int, Int, Int)//srow and scol


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int BaskerMatrix<Int,Entry,Exe_Space>::copy_values
  (
   Int _sr, Int _m, 
   Int _sc, Int _n, 
   Int _nnz,
   Int *_col_ptr, Int *_row_idx, Entry *_val
  )
  {
    if(nrow!=_m)
    {
      std::cout << "Error: Changed m size" << std::endl;
      std::cout << "old m: " << nrow << " new m: " << _m << std::endl;
      return BASKER_ERROR;
    }
    if(ncol==_n)
    {
      copy_vec(_col_ptr, _n+1,  col_ptr);
    }
    else
    {
      std::cout << "Error: Changed n size" << std::endl;
      std::cout << "old n: " << ncol << " new n: " << _n << std::endl;
      return BASKER_ERROR;
    }
    if(nnz == _nnz)
    {
      copy_vec(_row_idx, _nnz, row_idx);
      copy_vec(_val,_nnz,     val);
    }
    else
    {
      std::cout << "Error: Changed number of entries " << std::endl;
      std::cout << "old nnz: " << nnz << " new nnz: " << _nnz << std::endl;
      return BASKER_ERROR;
    }

    srow = _sr;
    erow = srow + _m;
    scol = _sc;
    ecol = scol + _n;

    return 0;
  }//end copy_values()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int BaskerMatrix<Int,Entry,Exe_Space>::copy_values
  (
   Int _m, Int _n, Int _nnz,
   Int *_col_ptr, Int *_row_idx, Entry *_val
  )
  {
    if(nrow!=_m)
    {
      std::cout << "Error: Changed m size" << std::endl;
      std::cout << "old m: " << nrow << " new m: " << _m << std::endl;
      return BASKER_ERROR;
    }
    if(ncol==_n)
    {
      copy_vec(_col_ptr, _n+1,  col_ptr);
    }
    else
    {
      std::cout << "Error: Changed n size" << std::endl;
      std::cout << "old n: " << ncol << " new n: " << _n << std::endl;
      return BASKER_ERROR;
    }
    if(nnz == _nnz)
    {
      copy_vec(_row_idx, _nnz, row_idx);
      copy_vec(_val,_nnz,     val);
    }
    else
    {
      std::cout << "Error: Changed number of entries " << std::endl;
      std::cout << "old nnz: " << nnz << " new nnz: " << _nnz << std::endl;
      return BASKER_ERROR;
    }

    return 0;
  }//end copy_values()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::malloc_perm(Int n)
  {
    BASKER_ASSERT(n > 0, "matrix malloc_perm");
    MALLOC_INT_1DARRAY(lpinv,n);    
 
    //Fix later.  //NDE determine what the issue is...
    init_perm();

  }//end init_perm(Int)

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry, Exe_Space>::init_perm()
  {
    for(Int i =0 ; i < nrow; i++)
    {
      lpinv(i) = BASKER_MAX_IDX;
    }
  }//end init_perm

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::malloc_union_bit()
  {
    BASKER_ASSERT(nrow > 0, "matrix_malloc_union_bit");
    MALLOC_BOOL_1DARRAY(union_bit, nrow);
  }//end malloc_union_bit
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::init_union_bit()
  {
    for(Int i =0 ; i <nrow; i++)
    {
      union_bit[i] = BASKER_FALSE;
    }
  }//end init_union_bit

  template <class Int, class Entry, class Exe_Space>
  void BaskerMatrix<Int,Entry,Exe_Space>::init_pend()
  {
    //printf( " > init_pend(ncol = %d)\n",ncol );
    if(ncol > 0)
    {
      BASKER_ASSERT((ncol+1)>0, "matrix init_pend")
      MALLOC_INT_1DARRAY(pend,ncol+1);
      for(Int i =0; i < ncol+1; ++i)
      {
        pend(i) = BASKER_MAX_IDX;
      }
    }
  }//end init_pend()

  template <class Int, class Entry, class Exe_Space>
  void BaskerMatrix<Int,Entry,Exe_Space>::clear_pend()
  {
    if(ncol > 0)
    {
      for(Int i = 0 ; i < ncol+1; ++i)
      {
        pend(i) = BASKER_MAX_IDX;
      }
    }
  }// end clear_pend()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int BaskerMatrix<Int,Entry,Exe_Space>::fill()
  {
    if(v_fill == BASKER_TRUE)
      return -1;

    #if 1
    Kokkos::deep_copy(col_ptr, 0);
    Kokkos::deep_copy(row_idx, 0);

    const Entry zero(0.0);
    Kokkos::deep_copy(val, zero);
    #else
    for(Int i = 0; i < ncol+1; i++)
    {
      col_ptr(i) = 0;
    }

    for(Int i = 0; i < nnz; i++)
    {
      row_idx(i) = 0;
    }

    for(Int i = 0; i < nnz; i++)
    {
      val(i) = 0;
    }
    #endif

    //#ifdef BASKER_INC_LVL
    if(inc_lvl_flg == BASKER_TRUE)
    {
      for(Int i = 0; i < nnz; i++)
      {
        inc_lvl(i) = 0;
      }
    }
    //#endif

    v_fill = BASKER_TRUE;
    return 0;
  }
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::convert2D
  (
   BASKER_MATRIX &M,
   BASKER_BOOL   alloc,
   Int kid,
   bool keep_zeros
  )
  {
    using STS = Teuchos::ScalarTraits<Entry>;
    using Mag = typename STS::magnitudeType;

    if(nnz == 0)
    {
      for(Int i = 0; i < ncol+1; i++)
      {
        col_ptr(i) = 0;
      }

      MALLOC_INT_1DARRAY(row_idx, 1);
      row_idx(0) = (Int) 0;
      MALLOC_ENTRY_1DARRAY(val, 1);
      val(0) = (Entry) 0;
      return;
    }

    //info();
    //We could check some flag ??
    //We assume a pre-scan has already happened
    if(alloc == BASKER_TRUE)
    {
      //printf("ALLOC\n");
      if(nnz > 0)
      {
        BASKER_ASSERT(nnz > 0, "matrix row nnz 2");
        MALLOC_INT_1DARRAY(row_idx, nnz);
      }
      else if(nnz == 0)
      {
        BASKER_ASSERT((nnz+1)>0, "matrix row nnz 3");
        MALLOC_INT_1DARRAY(row_idx, nnz+1);
      }
    }

    if(alloc == BASKER_TRUE)
    {
      if(nnz > 0)
      {
        BASKER_ASSERT(nnz > 0, "matrix nnz 4");
        MALLOC_ENTRY_1DARRAY(val, nnz);
      }
      else if(nnz == 0)
      {
        BASKER_ASSERT((nnz+1) > 0, "matrix nnz 5");
        MALLOC_ENTRY_1DARRAY(val, nnz+1);
      }
    }

    const Entry zero(0.0);

    anorm = abs(zero);
    Int temp_count = 0;
    for(Int k = scol; k < scol+ncol; ++k)
    {
      //note col_ptr[k-scol] contains the starting index (by find_2D_convert)
      // i.e., start of the diagonal block for U, or start of the first off-diagonal block for L
      if(col_ptr(k-scol) == BASKER_MAX_IDX)
      {
        col_ptr(k-scol) = temp_count;
        continue;
      }

      Int skipped_count = 0;
      Int kept_count = 0;
      Mag anorm_k (0.0);
      for(Int i = col_ptr(k-scol); i < M.col_ptr(k+1); i++)
      {
        Int j = M.row_idx(i);
        if(j >= srow+nrow)
        {
          if (!keep_zeros && (kept_count == 0 && skipped_count > 0)) {
            // if all were zero, then add the last entry to avoid empty column.
            row_idx(temp_count) = M.row_idx(i-1)-srow;
            val(temp_count) = M.val(i-1);
            anorm_k += abs(M.val(i-1));
            temp_count++;
          }
          break;
        }

        if(j < srow)
        {
          std::cout << std::endl 
                    << "BaskerMatrix::convert2D(kid = " << kid 
                    << "): j " << j 
                    << " srow: " << srow 
                    << " scol: " << scol 
                    << " nrow: " << nrow 
                    << " ncol: " << ncol 
                    << " with k: " << k
                    << " idx: " << i
                    << std::endl << std::endl;
          char error_msg[100];
          sprintf(error_msg, " ERROR: j is less than srow (j=%d, srow=%d, kid=%d)",(int)j,(int)srow,(int)kid);
          BASKER_ASSERT(0 == 1, error_msg);
        }
        if (keep_zeros || M.val(i) != zero)
        {
          row_idx(temp_count) = j-srow;
          val(temp_count) = M.val(i);

          anorm_k += abs(M.val(i));
          temp_count++;
          kept_count ++;
        } else {
          skipped_count ++;
        }
      }
      anorm = (anorm > anorm_k ? anorm : anorm_k);
      col_ptr(k-scol) = temp_count;
    }

    //NO!!1
    //Slide over the col counts
    for(Int i = ncol; i > 0; i--)
    {
      col_ptr(i) = col_ptr(i-1);
    }
    col_ptr(0) = (Int) 0;

  }//end convert2d(Matrix)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::init_inc_lvl()
  {
    MALLOC_INT_1DARRAY(inc_lvl, nnz+1);
    inc_lvl_flg = BASKER_TRUE;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::info()
  {
    std::cout << "\n Matrix Information: scol: " << scol
              << " ncol: " << ncol
              << " srow: " << srow 
              << " nrow: " << nrow
              << " nnz: " << nnz
              << std::endl;
  }//end info()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::print()
  {
    std::cout << "\n Matrix Print: \n Label: " << label.c_str()
              << " sr: " << srow
              << " sc: " << scol
              << "\n ncol: " << ncol
              << " nnz: " << nnz
              << std::endl;

    std::cout << " col_ptr: ";
    for(Int i=0; i < ncol+1; i++)
    {
      std::cout << col_ptr(i) << " , ";
    }
    std::cout << "\n row_idx: ";
    for(Int i=0; i < nnz; i++)
    {
      std::cout << row_idx(i) << " , ";
    }
    std::cout << "\n val: ";
    for(Int i=0; i < nnz; i++)
    {
      std::cout << val(i) << " , ";
    }

    std::cout << "\n END PRINT " << std::endl;

  }//end print()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::print_matrix(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (nrow > 0 && ncol > 0) {
      for(Int j = 0; j < ncol; j++) {
        for(Int k = col_ptr[j]; k < col_ptr[j+1]; k++) {
          fprintf(fp,"%d %d %.16e\n", (int)row_idx[k], (int)j, val[k]);
        }
      }
    }
    fclose(fp);
  }
}//end namespace basker

#endif //end BASKER_MATRIX_DEF_HPP
