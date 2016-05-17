#ifndef BASKER_MATRIX_DEF_HPP
#define BASKER_MATRIX_DEF_HPP

#include "basker_matrix_decl.hpp"
#include "basker_types.hpp"
#include "basker_util.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif


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
  }//end BaskerMatrix(label)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry, Exe_Space>::BaskerMatrix
  (Int _m, Int _n, Int _nnz,
   Int *col_ptr, Int *row_idx, Entry *val
   )
  {
    tpivot = 0;
    #ifdef BASKER_2DL
    p_size = 0;
    w_fill = BASKER_FALSE;
    #endif

    v_fill = BASKER_FALSE;
    init_matrix("matrix", _m,_n,_nnz, col_ptr, row_idx, val);
  }//end BaskerMatrix(int, int, int, int*, int*, Entry*)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  BaskerMatrix<Int,Entry, Exe_Space>::BaskerMatrix
  (
   string label, Int _m, Int _n, Int _nnz,
   Int *col_ptr, Int *row_idx, Entry *val
   )
  {
    tpivot = 0;
    #ifdef BASKER_2DL
    p_size = 0;
    w_fill = BASKER_FALSE;
    #endif

    v_fill = BASKER_FALSE;
    init_matrix(label, _m,_n,_nnz, col_ptr, row_idx, val);
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
    BASKER_ASSERT(ncol >= 0, "INIT_COL, ncol > 0");
    MALLOC_INT_1DARRAY(col_ptr, ncol+1);
    for(Int i = 0; i < ncol+1; ++i)
      {
	col_ptr(i) = (Int) BASKER_MAX_IDX;
      }

  }//end init_col()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::clean_col()
  {
    for(Int i = 0; i < ncol+1; ++i)
      {
	col_ptr(i) = (Int) BASKER_MAX_IDX;
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
    copy_vec(_col_ptr, _n+1,  col_ptr);
    copy_vec(_row_idx, _nnz, row_idx);
    copy_vec(_val,_nnz,     val);
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
   Int *_col_ptr, Int *_row_idx, Entry *_val)
  {
    if(nrow!=_m)
      {
	printf("Error: Changed m size\n");
	printf("old m: %d new m: %d \n", 
	       nrow , _m);
	return BASKER_ERROR;
      }
    if(ncol==_n)
      {
	copy_vec(_col_ptr, _n+1,  col_ptr);
      }
    else
      {
	printf("ERROR: Changed n size \n");
	printf("old n: %d new n: %d \n",
	       ncol, _n);
	return BASKER_ERROR;
      }
    if(nnz == _nnz)
      {
	copy_vec(_row_idx, _nnz, row_idx);
	copy_vec(_val,_nnz,     val);
      }
    else
      {
	printf("ERROR: Changed number of entries \n");
	printf("old nnz: %d new nnz: %d \n",
	       nnz, _nnz);
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
  int BaskerMatrix<Int,Entry,Exe_Space>::copy_values(Int _m, Int _n, Int _nnz,
				   Int *_col_ptr, Int *_row_idx, Entry *_val)
  {

    if(nrow!=_m)
      {
	printf("Error: Changed m size\n");
	printf("old m: %d new m: %d \n",
	       nrow, _m);
	return BASKER_ERROR;
      }
    if(ncol==_n)
      {
	copy_vec(_col_ptr, _n+1,  col_ptr);
      }
    else
      {
	printf("ERROR: Changed n size \n");
	printf("old n: %d new n: %d \n",
	       ncol, _n);
	return BASKER_ERROR;
      }
    if(nnz == _nnz)
      {
	copy_vec(_row_idx, _nnz, row_idx);
	copy_vec(_val,_nnz,     val);
      }
    else
      {
	printf("ERROR: Changed number of entries \n");
	printf("old nnz: %d new nnz: %d \n",
	       nnz, _nnz);
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
 
    //Fix later.
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
   Int kid
   )
  {
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
	else if(nnz ==0)
	  {
	    BASKER_ASSERT((nnz+1)>0, "matrix row nnz 3");
	    MALLOC_INT_1DARRAY(row_idx, nnz+1);
	  }
      }
    //init_value(row_idx, nnz, (Int) 0);
    //printf("clear row: %d \n", nnz);
    for(Int i = 0; i < nnz; ++i)
      {
	//printf("clear row_idx(%d) \n", i);
	row_idx(i) = 0;
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
    //init_value(val, nnz, (Entry) 0);
    for(Int i = 0; i < nnz; ++i)
      {
	val(i) = 0;
      }
    
    Int temp_count = 0;
    
    for(Int k = scol; k < scol+ncol; ++k)
      {
	//note col_ptr[k-scol] contains the starting index
	if(col_ptr(k-scol) == BASKER_MAX_IDX)
	  {
	    col_ptr(k-scol) = temp_count;
	    //printf("continue called, k: %d  \n", k);
	    continue;
	  }
  
	for(Int i = col_ptr(k-scol); i < M.col_ptr(k+1); i++)
	  {
	    Int j = M.row_idx(i);
	    //printf("i: %d j:%d \n", i,j);
	    if(j >= srow+nrow)
	      {
		break;
	      }
	    //printf("writing row_dix: %d i: %d  val: %d nnz: %d srow: %d nrow: %d \n",
	    //	   temp_count, i, j, nnz, 
	    //	   srow, nrow);
	    //BASKER_ASSERT(temp_count < nnz, "2DConvert, too many values");

	    //row_idx[temp_count] = j;

	    if(j-srow <0)
	      {
		printf("kid: %d j: %ld srow: %ld k: %d idx: %d   \n",
		       kid, j, srow, k, i);
		BASKER_ASSERT(0==1, "j-srow NO");
	      }


	    row_idx(temp_count) = j-srow;
	    val(temp_count) = M.val(i);
	    temp_count++;
	  }
	col_ptr(k-scol) = temp_count;
      }
   
    //col_ptr[0] = 0;
    
    //NO!!1
    //Slide over the col counts
    for(Int i = ncol; i > 0; i--)
    {
      col_ptr(i) = col_ptr(i-1);
     }
    col_ptr(0) = (Int) 0;
    
    //info();
    //print();
    
  }//end convert2d(Matrix)

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::init_inc_lvl()
  {

    MALLOC_INT_1DARRAY(inc_lvl, nnz+1);
    //for(Int i = 0; i < nnz+1; i++)
    //  {
    //	inc_lvl(i) = 0; 
    // }
    inc_lvl_flg = BASKER_TRUE;

  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::info()
  {

    printf("\n Matrix Information: scol: %d ncol: %d srow: %d nrow: %d nnz: %d \n \n", 
	   scol, ncol, srow, nrow, nnz);

  }//end info()
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void BaskerMatrix<Int,Entry,Exe_Space>::print()
  {
    printf("\n Matrix Print: \n Label: %s sr: %d sc: %d \n", 
	   label.c_str(), srow, scol);
    printf("ncol: %d nnz: %d \n", 
	   ncol, nnz);
    printf(" col_ptr: ");
    for(Int i=0; i < ncol+1; i++)
      {
        printf(" %d , ", col_ptr(i));
      }
    printf("\n row_idx: ");
    for(Int i=0; i < nnz; i++)
      {
        printf(" %d , ", row_idx(i));
      }
    printf("\n val: ");
    for(Int i=0; i < nnz; i++)
      {
        printf(" %e , " , val(i));
      }

    printf("\n END PRINT \n");

  }//end print()

}//end namespace basker


#ifndef KOKKOS_USING_EXPERIMENTAL_VIEW

namespace Kokkos
{
namespace Impl
{

//----------------------------Specialization for BaskerMatrix----------------//
  template <class ExecSpace, class Int, class Entry>
  struct ViewDefaultConstruct<ExecSpace, BaskerNS::BaskerMatrix<Int,Entry,ExecSpace>, true>
  {
    typedef BaskerNS::BaskerMatrix<Int,Entry,ExecSpace> type;
    type * const m_ptr;
 
    
    KOKKOS_FORCEINLINE_FUNCTION
    void operator()(const typename ExecSpace::size_type& i) const
    {
      new(m_ptr+i) type();
    }

    ViewDefaultConstruct( type * pointer, size_t capacity )
      :m_ptr( pointer )
    {

      Kokkos::RangePolicy<ExecSpace> range(0, capacity);
      parallel_for(range, *this);
      ExecSpace::fence();

    }//end constructor
   
  };//end ViewDefault Construct BaskerMatrix

}//end namespace imp
}//end namepace kokkos

#endif //end ifndef Kokkos_Using_Experimetnal_View

#endif //end BASKER_MATRIX_DEF_HPP
