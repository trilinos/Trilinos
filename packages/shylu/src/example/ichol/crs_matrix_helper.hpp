#pragma once
#ifndef __CRS_MATRIX_HELPER_HPP__
#define __CRS_MATRIX_HELPER_HPP__

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
    flat2hier(CrsFlatBase &flat, CrsHierBase &hier) {
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
          hier._aj[nnz] = ci[j];
          hier._ax[nnz].setView(&flat,     i, 1,
                                /**/   ci[j], 1);
        }
      }

      hier._ap[flat._m] = nnz;
      hier._nnz = nnz;

      return 0;
    } 
  };

}

#endif
