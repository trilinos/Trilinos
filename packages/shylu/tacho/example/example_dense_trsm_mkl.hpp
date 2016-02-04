#pragma once
#ifndef __EXAMPLE_DENSE_TRSM_MKL_HPP__
#define __EXAMPLE_DENSE_TRSM_MKL_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"  
#include "Teuchos_ScalarTraits.hpp" 

#include "util.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"
#include "dense_matrix_helper.hpp"

#include "task_view.hpp"
#include "task_factory.hpp"

#include "trsm.hpp"
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
  int exampleDenseTrsmMKL(const OrdinalType mmin,
                          const OrdinalType mmax,
                          const OrdinalType minc,
                          const OrdinalType k,
                          const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "DenseGemmMKL:: test matrices "
         <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << " , k = "<< k << endl;

    ostringstream os;
    os.precision(3);
    os << scientific;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      os.str("");

      DenseMatrixBaseType AA("AA", m, m), BB("BB", m, k), BC("BC", m, k);
      
      // setup upper triangular
      for (ordinal_type j=0;j<AA.NumCols();++j) {
        AA.Value(j,j) = 10.0;
        for (ordinal_type i=0;i<j;++i)
          AA.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
      }

      // setup one and right hand side is going to be overwritten by the product of AB
      for (ordinal_type j=0;j<BB.NumCols();++j)
        for (ordinal_type i=0;i<BB.NumRows();++i)
          BB.Value(i,j) = 1.0;

      Teuchos::BLAS<ordinal_type,value_type> blas;

      blas.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
                m, k, m,
                1.0,
                AA.ValuePtr(), AA.ColStride(),
                BB.ValuePtr(), BB.ColStride(),
                0.0,
                BC.ValuePtr(), BC.ColStride());
      BB.copy(BC);

      const double flop = get_flop_trsm_upper<value_type>(m, k);

      os << "DenseTrsmMKL:: m = " << m << " k = " << k;
      {
        timer.reset();
        Teuchos::BLAS<ordinal_type,value_type> blas;

        const ordinal_type mm = AA.NumRows();
        const ordinal_type nn = BB.NumCols();

        blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                  Teuchos::NON_UNIT_DIAG,
                  mm, nn,
                  1.0,
                  AA.ValuePtr(), AA.ColStride(),
                  BB.ValuePtr(), BB.ColStride());
        t = timer.seconds();
        os << ":: MKL Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }
      cout << os.str() << endl;
    }

    return r_val;
  }
}

#endif
