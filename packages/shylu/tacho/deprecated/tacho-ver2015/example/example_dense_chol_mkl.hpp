#pragma once
#ifndef __EXAMPLE_DENSE_CHOL_MKL_HPP__
#define __EXAMPLE_DENSE_CHOL_MKL_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"
#include "Teuchos_ScalarTraits.hpp"

#include "util.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "chol.hpp"
#include "dense_flop.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleDenseCholMKL(const OrdinalType mmin,
                          const OrdinalType mmax,
                          const OrdinalType minc,
                          const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;
    
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "DenseCholMKL:: test matrices "
         <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << endl;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      DenseMatrixBaseType TT("TT", m, m), AA("AA", m, m);

      // random T matrix
      for (ordinal_type j=0;j<AA.NumCols();++j) {
        for (ordinal_type i=0;i<AA.NumRows();++i)
          TT.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        TT.Value(j,j) = abs(TT.Value(j,j));
      }

      {
        Teuchos::BLAS<ordinal_type,value_type> blas;

        // should be square
        const ordinal_type nn = AA.NumRows();
        const ordinal_type kk = TT.NumRows();

        blas.HERK(Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                  nn, kk,
                  1.0,
                  TT.ValuePtr(), TT.ColStride(),
                  0.0,
                  AA.ValuePtr(), AA.ColStride());
      }

      const double flop = get_flop_chol<value_type>(m);
      {
        timer.reset();
        const ordinal_type mm = AA.NumRows();

        int ierr = 0;
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        lapack.POTRF('U',
                     mm,
                     AA.ValuePtr(), AA.ColStride(),
                     &ierr);
        t = timer.seconds();
        if (ierr)
          ERROR(">> Fail to perform MKL Cholesky");
        cout << "DenseCholMKL:: Performance = " << (flop/t/1.0e9) << " [GFLOPs]" << endl;
      }
    }

    return r_val;
  }
}

#endif
