#ifndef __TACHO_EXAMPLE_DENSE_GEMM_MKL_HPP__
#define __TACHO_EXAMPLE_DENSE_GEMM_MKL_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "Tacho_Util.hpp"

#include "Tacho_DenseMatrixBase.hpp"
#include "Tacho_DenseFlopCount.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif

namespace Tacho {

  template<typename HostSpaceType>
  int exampleDenseGemmMKL(const ordinal_type mmin,
                          const ordinal_type mmax,
                          const ordinal_type minc,
                          const ordinal_type k,
                          const bool verbose) {
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> DenseMatrixBaseType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;

    double t = 0.0;

    std::cout << "DenseGemmMKL:: test matrices "
              <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << " , k = "<< k << std::endl;

    std::ostringstream os;
    os.precision(3);
    os << std::scientific;
        
    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      os.str("");
      
      DenseMatrixBaseType AA("AA", mmax, k), BB("BB", k, mmax), CC("CC", m, m);

      for (ordinal_type j=0;j<AA.NumCols();++j)
        for (ordinal_type i=0;i<AA.NumRows();++i)
          AA.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;

      for (ordinal_type j=0;j<BB.NumCols();++j)
        for (ordinal_type i=0;i<BB.NumRows();++i)
          BB.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;

      for (ordinal_type j=0;j<CC.NumCols();++j)
        for (ordinal_type i=0;i<CC.NumRows();++i)
          CC.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;

      const double flop = DenseFlopCount<value_type>::Gemm(m, m, k);

      os << "DenseGemmMKL:: m = " << m << " n = " << m << " k = " << k;
      {
        timer.reset();
        Teuchos::BLAS<ordinal_type,value_type> blas;

        const ordinal_type mm = CC.NumRows();
        const ordinal_type nn = CC.NumCols();
        const ordinal_type kk = BB.NumRows();

        blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                  mm, nn, kk,
                  1.0,
                  AA.ValuePtr(), AA.ColStride(),
                  BB.ValuePtr(), BB.ColStride(),
                  1.0,
                  CC.ValuePtr(), CC.ColStride());
        t = timer.seconds();
        os << ":: MKL Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }
      std::cout << os.str() << std::endl;
    }

    return r_val;
  }
}

#endif
