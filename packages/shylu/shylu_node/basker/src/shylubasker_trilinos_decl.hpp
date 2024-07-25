// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_TRILINOS_DECL_HPP
#define SHYLUBASKER_TRILINOS_DECL_HPP

#include "shylubasker_decl.hpp"

namespace BaskerNS {
  // This handles two different types for ints for example with Tpetra CRS Matrix
  template <class Int, class Entry, class Exe_Space>
  class BaskerTrilinosInterface : public BaskerNS::Basker<Int, Entry, Exe_Space>
  {
    friend class BaskerNS::Basker<Int,Entry,Exe_Space>;

    // Member type for copying ptr values of type non-matching with indices ptr to force matching types
    Int* matching_type_col_ptr;

  public:

    BASKER_INLINE
    BaskerTrilinosInterface() : BaskerNS::Basker<Int, Entry, Exe_Space>()
    { matching_type_col_ptr = nullptr; }

    // No polymorphic behavior so a virtual destructor is not necessary
    BASKER_INLINE
    ~BaskerTrilinosInterface() {
      if ( matching_type_col_ptr != nullptr )
        delete[] matching_type_col_ptr ;
    }

    template < typename U >
    BASKER_INLINE
    int Symbolic
    (Int nrow, Int ncol, typename std::enable_if< !std::is_same<Int,U>::value, Int>::type nnz, U *col_ptr, Int *row_idx, Entry *val, bool _crs_transpose_needed = false)
    {
      // NDE: Allocate a new array for the non-matching type; copy into that then pass that along
      //Int* matching_type_col_ptr = new Int[ncol+1];
      matching_type_col_ptr = new Int[ncol+1];

    #ifdef KOKKOS_ENABLE_OPENMP
    #pragma omp parallel for
    #endif
      for (Int i = 0; i < ncol+1; ++i)
      {
        matching_type_col_ptr[i] = col_ptr[i];
      }


      int return_value = 
        BaskerNS::Basker<Int, Entry, Exe_Space>::Symbolic
        (
         nrow,
         ncol,
         nnz, 
         matching_type_col_ptr, 
         row_idx, 
         val,
         _crs_transpose_needed
        );

      //delete [] matching_type_col_ptr;

      return return_value;
    }

    template < typename U >
    BASKER_INLINE
    int Symbolic
    (Int nrow, Int ncol, typename std::enable_if< std::is_same<Int,U>::value, Int>::type nnz, U *col_ptr, Int *row_idx, Entry *val, bool _crs_transpose_needed = false)
    {
      // NDE: Allocate a new array for the non-matching type; copy into that then pass that along

      int return_value = 
        BaskerNS::Basker<Int, Entry, Exe_Space>::Symbolic
        (
         nrow,
         ncol,
         nnz, 
         col_ptr, 
         row_idx, 
         val,
         _crs_transpose_needed
        );

      return return_value;
    }


    template < typename U >
    BASKER_INLINE
    int Factor
    (Int nrow, Int ncol, typename std::enable_if< !std::is_same<Int,U>::value, Int>::type nnz, U *col_ptr, Int *row_idx, Entry *val)
    {
      //Int* matching_type_col_ptr = new Int[ncol+1];

    #ifdef KOKKOS_ENABLE_OPENMP
    #pragma omp parallel for
    #endif
      for (Int i = 0; i < ncol+1; ++i)
      {
        matching_type_col_ptr[i] = col_ptr[i];
      }

      int return_value = 
        BaskerNS::Basker<Int,Entry,Exe_Space>::Factor
        (
         nrow, 
         ncol,
         nnz, 
         matching_type_col_ptr, 
         row_idx, 
         val
        ); 

      //delete [] matching_type_col_ptr;

      return return_value;
    }

    template < typename U >
    BASKER_INLINE
    int Factor
    (Int nrow, Int ncol, typename std::enable_if< std::is_same<Int,U>::value, Int>::type nnz, U *col_ptr, Int *row_idx, Entry *val)
    {

      int return_value = 
        BaskerNS::Basker<Int,Entry,Exe_Space>::Factor
        (
         nrow, 
         ncol,
         nnz, 
         col_ptr, 
         row_idx, 
         val
        ); 

      return return_value;
    }

  }; //end BaskerTrilinosInterface class

} //end namespace BaskerNS
#endif
