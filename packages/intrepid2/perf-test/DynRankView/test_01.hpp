// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Performance test comparing dynrankview overhead
\author Created by Kyungjoo Kim.
*/

#include "Intrepid2_config.h"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Intrepid2 {

  namespace Test {

    namespace Serial {
      
      // compute determinant for rank 2 array
      template<typename inMatViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      typename inMatViewType::value_type
      det( const inMatViewType inMat ) {
        const auto dim = inMat.extent(0);
        return ( dim == 3 ? ( inMat(0,0) * inMat(1,1) * inMat(2,2) +
                              inMat(1,0) * inMat(2,1) * inMat(0,2) +
                              inMat(2,0) * inMat(0,1) * inMat(1,2) -
                              inMat(2,0) * inMat(1,1) * inMat(0,2) -
                              inMat(0,0) * inMat(2,1) * inMat(1,2) -
                              inMat(1,0) * inMat(0,1) * inMat(2,2) ) :
                 dim == 2 ? ( inMat(0,0) * inMat(1,1) -
                              inMat(0,1) * inMat(1,0) ) :
                 /**/       ( inMat(0,0) ) );
      }

    }

    template<typename detArrayViewType,
             typename inMatViewType>
    struct F_det {
      detArrayViewType _detArray;
      inMatViewType    _inMats;
      
      KOKKOS_INLINE_FUNCTION
      F_det( detArrayViewType   detArray_,
             inMatViewType      inMats_ )
        : _detArray(detArray_), _inMats(inMats_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        // comapring left and right stride for serial
        // 15 % difference, right stride access is better (of course)

        // right
        const auto stride = _inMats.extent(1);
        const auto 
          i = iter / stride, 
          j = iter % stride;

        // left
        // const auto stride = _inMats.extent(0);
        // const auto 
        //   i = iter % stride, 
        //   j = iter / stride;

        // current implementation subview overhead exists
        // {
        //   auto mat = Kokkos::subview(_inMats, i, j, Kokkos::ALL(), Kokkos::ALL());
        //   _detArray(i, j) = Serial::det(mat);
        // }

        // direct access (base is 8 rank)
        {
          const auto dim = _inMats.extent(2);
          const auto val = ( dim == 3 ? ( _inMats(i,j,0,0) * _inMats(i,j,1,1) * _inMats(i,j,2,2) +
                                          _inMats(i,j,1,0) * _inMats(i,j,2,1) * _inMats(i,j,0,2) +
                                          _inMats(i,j,2,0) * _inMats(i,j,0,1) * _inMats(i,j,1,2) -
                                          _inMats(i,j,2,0) * _inMats(i,j,1,1) * _inMats(i,j,0,2) -
                                          _inMats(i,j,0,0) * _inMats(i,j,2,1) * _inMats(i,j,1,2) -
                                          _inMats(i,j,1,0) * _inMats(i,j,0,1) * _inMats(i,j,2,2) ) :
                             dim == 2 ? ( _inMats(i,j,0,0) * _inMats(i,j,1,1) -
                                          _inMats(i,j,0,1) * _inMats(i,j,1,0) ) :
                             /**/       ( _inMats(i,j,0,0) ) );
          
          _detArray(i, j) = val;
        }

      }
    };

    template<typename ValueType, typename DeviceSpaceType>
    int DynRankView_PerfTest01(const ordinal_type nworkset,
                               const ordinal_type C,
                               const ordinal_type P,
                               const ordinal_type D,
                               const bool verbose) {

      Teuchos::RCP<std::ostream> verboseStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        verboseStream = Teuchos::rcp(&std::cout, false);
      else
        verboseStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

      *verboseStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*verboseStream, false);
      *verboseStream << "HostSpace::    ";   HostSpaceType().print_configuration(*verboseStream, false);

      std::cout
        << "===============================================================================\n" 
        << " Performance Test for measuring DynRankView Overhead \n"
        << " # of workset = " << nworkset << "\n" 
        << " Test Array Structure (C,P,D) = " << C << ", " << P << ", " << D << "\n"
        << "===============================================================================\n";

      Kokkos::Timer timer;
      double t_dynrankview[20] = {}, t_view[20] = {};
      int errorFlag = 0, itest = 0;

      *verboseStream
        << "\n"
        << "===============================================================================\n"
        << "TEST 1: RealSpaceTools det \n"
        << "        input (C,P,D,D), output (C,P) \n"
        << "===============================================================================\n";
      
      try {
        *verboseStream << " -> Testing View \n";

        { // Kokkos View
          typedef Kokkos::View<ValueType****,DeviceSpaceType> inViewType;
          typedef Kokkos::View<ValueType**,  DeviceSpaceType> outViewType;
          //typedef Kokkos::View<ValueType********,DeviceSpaceType> inViewType;
          //typedef Kokkos::View<ValueType********,  DeviceSpaceType> outViewType;
          typedef F_det<outViewType,inViewType> FunctorType;

          inViewType  in ("inView",  C,P,D,D);
          outViewType out("outView", C,P);          

          const auto loopSize = C*P;
          Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

          Kokkos::deep_copy(in, 1.0);                          

          DeviceSpaceType().fence();
          timer.reset();

          for (ordinal_type i=0;i<nworkset;++i) 
            Kokkos::parallel_for( policy, FunctorType(out, in) );

          DeviceSpaceType().fence();
          t_view[itest] = timer.seconds();
        }

        *verboseStream << " -> Testing DynRankView \n";

        { // Kokkos DynRankView
          typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> ViewType;
          typedef F_det<ViewType,ViewType> FunctorType;

          ViewType in("inDynRankView", C,P,D,D), out("outDynRankView", C,P);          
          const auto loopSize = C*P;
          Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

          Kokkos::deep_copy(in, 1.0);                          


          DeviceSpaceType().fence();
          timer.reset();

          for (ordinal_type i=0;i<nworkset;++i) 
            Kokkos::parallel_for( policy, FunctorType(out, in) );

          DeviceSpaceType().fence();
          t_dynrankview[itest] = timer.seconds();
        }
      } catch (std::exception err) {
        *verboseStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *verboseStream << err.what() << '\n';
        *verboseStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      ++itest;

      std::cout 
        << "===============================================================================\n";

      for (auto i=0;i<itest;++i) {
        const auto ratio = (t_view[i]/t_dynrankview[i]);
        std::cout 
          << "TEST " << i 
          << ":  t_view = " << t_view[i] 
          << ",  t_dynrankview = " << t_dynrankview[i] 
          << ", ratio (v/drv) = " << ratio << "\n";

        const auto tol = 0.95;
        if (ratio < tol)  {
          std::cout << "Performance test failed as view/dynrankview is below " << tol << std::endl;
          errorFlag = -1;
        }
        
      }
      std::cout 
        << "===============================================================================\n";
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  } // end of namespace TEST
} // end of namespace Intrepid2
