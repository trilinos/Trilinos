#ifndef __TACHO_EXAMPLE_PARDISO_CHOL_HPP__
#define __TACHO_EXAMPLE_PARDISO_CHOL_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"
#include "Tacho_DenseMatrixBase.hpp"
#include "Tacho_DenseMatrixTools.hpp"
#include "Tacho_MatrixMarket.hpp"
#include "Tacho_CrsData.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#include "Tacho_ExamplePardiso.hpp"

namespace Tacho {

  template<typename DeviceSpaceType>
  KOKKOS_INLINE_FUNCTION
  int examplePardisoChol(const string file_input,
                         const int  nrhs,
                         const int  nthreads,
                         const bool skip_factorize,
                         const bool skip_solve,
                         const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
    std::cout << std::endl;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   CrsMatrixBaseHostType;
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> DenseMatrixBaseHostType;

    int r_val = 0;

    // mkl nthreads setting 
    mkl_set_dynamic(0);
    mkl_set_num_threads(nthreads);

    Kokkos::Impl::Timer timer;
    double t = 0.0;
    Pardiso pardiso;

    std::cout << "PardisoChol:: init" << std::endl;
    {
      timer.reset();
      r_val = pardiso.init<value_type,AlgoChol::ExternalPardiso>();
      t = timer.seconds();
      
      if (r_val) {
        std::cout << "PardisoChol:: Pardiso init error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      }
    }
    std::cout << "PardisoChol:: init ::time = " << t << std::endl;
    
    std::cout << "PardisoChol:: import input file = " << file_input << std::endl;
    CrsMatrixBaseHostType AA("AA");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }

      const auto extension = file_input.substr(file_input.find_last_of(".") + 1);
      if        (extension == "mtx") {
        std::cout << "PardisoChol:: Input matrix is MatrixMarket format" << std::endl;
        MatrixMarket::read(AA, in);
      } else if (extension == "crs") {
        std::cout << "PardisoChol:: Input matrix is CRS data format" << std::endl;
        CrsData::read(AA, in);
        CrsMatrixTools::sortColumnsPerRow(AA);
      }
      
      // somehow pardiso does not like symmetric matrix
      // CrsMatrixBaseHostType AA_tmp("AA_tmp");
      // AA_tmp.createConfTo(AA);
      // CrsMatrixTools::copy(AA_tmp,
      //                      Uplo::Lower, 0, 
      //                      AA);
      // AA = AA_tmp;
    }
    t = timer.seconds();

    // row pointer array
    typename CrsMatrixBaseHostType::size_type_array row_ptr;
    {
      const auto m = AA.NumRows();
      row_ptr = typename CrsMatrixBaseHostType::size_type_array("PardisoChol::rowptr", m);
      const auto row_ptr_begin = AA.RowPtrBegin();
      const auto row_ptr_end   = AA.RowPtrEnd();
      for (auto i=0;i<m;++i) 
        row_ptr(i) = row_ptr_begin(i);
      row_ptr(m) = row_ptr_end(m-1);
    }
    std::cout << "PardisoChol:: import input file::time = " << t << std::endl;
    
    DenseMatrixBaseHostType BB("BB",  AA.NumRows(), nrhs), XX, RR;
    XX.createConfTo(BB);
    RR.createConfTo(XX);

    DenseMatrixBaseHostType PP("PP",  AA.NumRows(), 1);

    {
      const auto m = AA.NumRows();
      srand(time(NULL));
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) 
          XX.Value(i, rhs) = ((value_type)rand()/(RAND_MAX));
        
        // matvec
        HostSpaceType::execution_space::fence();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                             [&](const ordinal_type i) {
                               const auto nnz  = AA.NumNonZerosInRow(i);
                               const auto cols = AA.ColsInRow(i);
                               const auto vals = AA.ValuesInRow(i);

                               value_type tmp = 0;
                               for (ordinal_type j=0;j<nnz;++j)
                                 tmp += vals(j)*XX.Value(cols(j), rhs);
                               BB.Value(i, rhs) = tmp;
                             } );
        HostSpaceType::execution_space::fence();
      }
      DenseMatrixTools::copy(RR, XX); // keep solution on RR
    }

    pardiso.setProblem(AA.NumRows(),
                       (double*)AA.Values().data(),
                       (int*)row_ptr.data(),
                       (int*)AA.Cols().data(),
                       (int*)PP.Values().data(),
                       BB.NumCols(),
                       (double*)BB.Values().data(),
                       (double*)XX.Values().data());
    
    std::cout << "PardisoChol:: analyze matrix" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::Analyze);
      t = timer.seconds();
      
      if (r_val) {
        std::cout << "PardisoChol:: Pardiso analyze error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::Analyze) << std::endl;
      }
    }
    std::cout << "PardisoChol:: analyze matrix::time = " << t << std::endl;

    if (!skip_factorize) {
      std::cout << "PardisoChol:: factorize matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Factorize);
        t = timer.seconds();
        
        if (r_val) {
          std::cout << "PardisoChol:: Pardiso factorize error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Factorize) << std::endl;
        }
      }
      std::cout << "PardisoChol:: factorize matrix::time = " << t << std::endl;
    }

    if (!skip_factorize && !skip_solve) {
      std::cout << "PardisoChol:: solve matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Solve);
        t = timer.seconds();
        
        if (r_val) {
          std::cout << "PardisoChol:: Pardiso solve error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Solve) << std::endl;
        }
      }
      std::cout << "PardisoChol:: solve matrix::time = " << t << std::endl;
    }
    
    {
      double error = 0, norm = 0;
      const auto m = AA.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) {
          {
            const auto val = Util::abs(XX.Value(i, rhs) - RR.Value(i, rhs));
            error += val*val;
          }
          {
            const auto val = Util::abs(RR.Value(i, rhs));
            norm  += val*val;
          }
        }
      }
      std::cout << "PardisoChol:: error = " << error << " , norm = " << norm << std::endl;
    }

    std::cout << "PardisoChol:: release all" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::ReleaseAll);
      t = timer.seconds();

      if (r_val) {
        std::cout << "PardisoChol:: release error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::ReleaseAll) << std::endl;
      }
    }
    std::cout << "PardisoChol:: release all::time = " << t << std::endl;
    
    return r_val;
  }

}

#endif
#endif
