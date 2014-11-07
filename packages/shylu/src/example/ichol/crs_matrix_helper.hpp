#pragma once
#ifndef __CRS_MATRIX_HELPER_HPP__
#define __CRS_MATRIX_HELPER_HPP__

/// \file crs_matrix_helper.hpp
/// \brief This file includes utility functions to convert between flat and hierarchical matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include <Kokkos_Core.hpp> 

#include "util.hpp"
#include "coo.hpp"

namespace Example { 

  using namespace std;


  class CrsMatrixHelper {
  public:
    template<typename CrsFlatBase,
             typename CrsHierBase>
    KOKKOS_INLINE_FUNCTION
    static int
    flat2hier(CrsFlatBase &flat, 
              CrsHierBase &hier) {
      typedef typename CrsHierBase::ordinal_type            ordinal_type;
      typedef typename CrsHierBase::size_type               size_type;
      typedef typename CrsHierBase::ordinal_type_array_view ordinal_type_array_view;

      hier.createInternalArrays(flat._m, flat._n, flat._nnz);

      size_type nnz = 0;
      for (ordinal_type i=0;i<flat._m;++i) {
        ordinal_type jsize = flat._ap[i+1] - flat._ap[i];

        hier._ap[i] = nnz;
        ordinal_type_array_view ci = flat.ColsInRow(i);
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
    static int
    flat2hier(CrsFlatBase &flat, 
              CrsHierBase &hier,
              const typename CrsHierBase::ordinal_type       nblks,
              const typename CrsHierBase::ordinal_type_array range,
              const typename CrsHierBase::ordinal_type_array tree) {
      typedef typename CrsHierBase::ordinal_type            ordinal_type;
      typedef typename CrsHierBase::size_type               size_type;

      typedef typename CrsHierBase::ordinal_type_array      ordinal_type_array;
      typedef typename CrsHierBase::ordinal_type_array_view ordinal_type_array_view;
      typedef typename CrsHierBase::value_type_array_view   value_type_array_view;
      
      ordinal_type_array tmp = ordinal_type_array("flat2hier:tmp", nblks+1);
      size_type nnz = 0;

      // count nnz and nnz in rows for lower triangular matrix
      for (ordinal_type i=0;i<nblks;++i) {
        ordinal_type j = i;
        for ( ;j != -1;++nnz) { 
          ++tmp[j];
          j = tree[j];
        }
      }

      // count nnz and nnz in rows for lower triangular matrix      
      hier.createInternalArrays(nblks, nblks, nnz);
      for (ordinal_type i=1;i<(nblks+1);++i) 
        hier._ap[i] = hier._ap[i-1] + tmp[i-1];

      for (ordinal_type i=0;i<(nblks+1);++i) 
        tmp[i] = hier._ap[i];

      for (ordinal_type i=0;i<nblks;++i) {
        ordinal_type j = i;
        for ( ;j != -1; ) {
          hier._aj[tmp[j]] = i;
          hier._ax[tmp[j]].setView(&flat, range[j], (range[j+1] - range[j]),
                                   /**/   range[i], (range[i+1] - range[i]));
          ++tmp[j];
          j = tree[j];
        }
      }

      nnz = 0;
      for (ordinal_type i=0;i<hier._m;++i) {
        ordinal_type jsize = hier._ap[i+1] - hier._ap[i];
        ordinal_type_array_view ci = hier.ColsInRow(i);
        for (ordinal_type j=0;j<jsize;++j,++nnz) {
          cout << hier._aj[nnz] << endl;
        }
      }
      
      // // create upper triangular block matrix
      // CrsHierBase up("flat2hier::up", nblks, nblks, nnz);

      // nnz = 0;
      // for (ordinal_type i=0;i<nblks;++i) {
      //   up._ap[i] = nnz;
      //   ordinal_type j = i;
      //   for ( ;j != -1; ++nnz) {
      //     up._aj[nnz] = j;
      //     up._ax[nnz].setView(&flat, range[i], (range[i+1] - range[i]),
      //                          /**/   range[j], (range[j+1] - range[j]));
      //     j = tree[j];
      //   }
      // }
      
      // up._ap[nblks] = nnz;
      // up._nnz = nnz;

      return 0;
    } 

  };

}

#endif
