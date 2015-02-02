#pragma once
#ifndef __CRS_MATRIX_HELPER_IMPL_HPP__
#define __CRS_MATRIX_HELPER_IMPL_HPP__

/// \file crs_matrix_helper_impl.hpp
/// \brief This file includes utility functions to convert between flat and hierarchical matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include <Kokkos_Core.hpp> 

#include "util.hpp"

namespace Example { 

  using namespace std;

  template<typename CrsFlatBase,
           typename CrsHierBase>
  KOKKOS_INLINE_FUNCTION
  int
  CrsMatrixHelper::flat2hier(CrsFlatBase &flat, 
                             CrsHierBase &hier) {
    typedef typename CrsHierBase::ordinal_type           ordinal_type;
    typedef typename CrsHierBase::size_type              size_type;
    typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    
    hier.createInternalArrays(flat._m, flat._n, flat._nnz);
    
    size_type nnz = 0;
    for (ordinal_type i=0;i<flat._m;++i) {
      ordinal_type jsize = flat._ap[i+1] - flat._ap[i];
      
      hier._ap[i] = nnz;
      ordinal_type_array_ptr ci = flat.ColsInRow(i);
      for (ordinal_type j=0;j<jsize;++j,++nnz) {
        if (i >= ci[j]) {
          hier._aj[nnz] = ci[j];
          hier._ax[nnz].setView(&flat,     i, 1,
                                /**/   ci[j], 1);
        }
      }
    }
    
    hier._ap[flat._m] = nnz;
    hier._nnz = nnz;
    
    return 0;
  } 

  template<typename CrsFlatBase,
           typename CrsHierBase>
  KOKKOS_INLINE_FUNCTION
  int
  CrsMatrixHelper::flat2hier(int uplo,
                             CrsFlatBase &flat, 
                             CrsHierBase &hier,
                             const typename CrsHierBase::ordinal_type       nblks,
                             const typename CrsHierBase::ordinal_type_array range,
                             const typename CrsHierBase::ordinal_type_array tree) {
    switch(uplo) {
    case Uplo::Upper: return flat2hier_upper(flat, hier, nblks, range, tree); 
    case Uplo::Lower: return flat2hier_lower(flat, hier, nblks, range, tree); 
    }
    return -1;
  }
  
  template<typename CrsFlatBase,
           typename CrsHierBase>
  KOKKOS_INLINE_FUNCTION
  int
  CrsMatrixHelper::flat2hier_upper(CrsFlatBase &flat, 
                                   CrsHierBase &hier,
                                   const typename CrsHierBase::ordinal_type       nblks,
                                   const typename CrsHierBase::ordinal_type_array range,
                                   const typename CrsHierBase::ordinal_type_array tree) {
    typedef typename CrsHierBase::ordinal_type            ordinal_type;
    typedef typename CrsHierBase::size_type               size_type;
    
    typedef typename CrsHierBase::ordinal_type_array     ordinal_type_array;
    typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    typedef typename CrsHierBase::value_type_array_ptr   value_type_array_ptr;
    
    size_type nnz = 0;
    
    // count nnz and nnz in rows for the upper triangular hier matrix
    for (ordinal_type i=0;i<nblks;++i) 
      for (ordinal_type j=i;j != -1;++nnz,j=tree[j]) ;
    
    // create upper triangular block matrix
    hier.createInternalArrays(nblks, nblks, nnz);    

    nnz = 0;
    for (ordinal_type i=0;i<nblks;++i) {
      hier._ap[i] = nnz;
      for (ordinal_type j=i;j != -1;++nnz,j=tree[j]) {
        hier._aj[nnz] = j;
        hier._ax[nnz].setView(&flat, range[i], (range[i+1] - range[i]),
                              /**/   range[j], (range[j+1] - range[j]));
      }
    }
    
    hier._ap[nblks] = nnz;
    hier._nnz = nnz;
  }

  template<typename CrsFlatBase,
           typename CrsHierBase>
  KOKKOS_INLINE_FUNCTION
  int
  CrsMatrixHelper::flat2hier_lower(CrsFlatBase &flat, 
                                   CrsHierBase &hier,
                                   const typename CrsHierBase::ordinal_type       nblks,
                                   const typename CrsHierBase::ordinal_type_array range,
                                   const typename CrsHierBase::ordinal_type_array tree) {
    typedef typename CrsHierBase::ordinal_type           ordinal_type;
    typedef typename CrsHierBase::size_type              size_type;
    
    typedef typename CrsHierBase::ordinal_type_array     ordinal_type_array;
    typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    typedef typename CrsHierBase::value_type_array_ptr   value_type_array_ptr;
    
    ordinal_type_array tmp = ordinal_type_array("flat2hier:tmp", nblks+1);
    size_type nnz = 0;
    
    // count nnz and nnz in rows for lower triangular matrix
    for (ordinal_type i=0;i<nblks;++i)         
      for (ordinal_type j=i;j != -1;++nnz) { 
        ++tmp[j];
        j = tree[j];
      }
    
    // count nnz and nnz in rows for lower triangular matrix      
    hier.createInternalArrays(nblks, nblks, nnz);
    for (ordinal_type i=1;i<(nblks+1);++i) 
      hier._ap[i] = hier._ap[i-1] + tmp[i-1];
    
    for (ordinal_type i=0;i<(nblks+1);++i) 
      tmp[i] = hier._ap[i];
    
    for (ordinal_type i=0;i<nblks;++i) 
      for (ordinal_type j=i;j != -1;j=tree[j]) {
        hier._aj[tmp[j]] = i;
        hier._ax[tmp[j]].setView(&flat, range[j], (range[j+1] - range[j]),
                                 /**/   range[i], (range[i+1] - range[i]));
        ++tmp[j];
      }
    
    return 0;
  } 

}


#endif

