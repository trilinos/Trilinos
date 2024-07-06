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
      template<typename OutputViewType,
               typename inputViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      void
      clone( /**/  OutputViewType output,
             const inputViewType  input ) {
        const ordinal_type iend = output.extent(0);
        const ordinal_type jend = output.extent(1);

        for (ordinal_type i=0;i<iend;++i)
          for (ordinal_type j=0;j<jend;++j)
            output(i, j) = input(i, j);
      }
      
    }

    template<typename outputFieldViewType,
             typename inputFieldViewType,
             int innerRank>
    struct F_clone {
      outputFieldViewType _outputFields;
      inputFieldViewType _inputFields;
      
      KOKKOS_INLINE_FUNCTION
      F_clone( outputFieldViewType outputFields_,
               inputFieldViewType inputFields_ )
        : _outputFields(outputFields_), _inputFields(inputFields_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        // right
        const auto stride = _outputFields.extent(1);
        const auto 
          i = iter / stride, 
          j = iter % stride;

        if (innerRank == -1) {
          // current implementation subview overhead exists
          auto out = Kokkos::subview(_outputFields, i, j, Kokkos::ALL(), Kokkos::ALL());
          auto in  = Kokkos::subview(_inputFields,     j, Kokkos::ALL(), Kokkos::ALL());
          Serial::clone( out, in );
        } else {
          // rank specific direct access to containers
          switch (innerRank) {
          case 0: {
            _outputFields(i, j) = _inputFields(j);
            break;
          }
          case 1: {
            const ordinal_type kend = _outputFields.extent(2);
            for (ordinal_type k=0;k<kend;++k)
              _outputFields(i, j, k) = _inputFields(j, k);
            break;
          }
          case 2: {
            const ordinal_type kend = _outputFields.extent(2);
            const ordinal_type lend = _outputFields.extent(3);
            for (ordinal_type k=0;k<kend;++k)
              for (ordinal_type l=0;l<lend;++l)
                _outputFields(i, j, k, l) = _inputFields(j, k, l);
            break;
          }
          }
        }
      }
    };
    
    template<typename ValueType, typename DeviceSpaceType>
    int DynRankView_PerfTest02(const ordinal_type nworkset,
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
      *verboseStream << "\n";

      std::cout
        << "===============================================================================\n" 
        << " Performance Test for measuring subview Overhead \n"
        << " # of workset = " << nworkset << "\n" 
        << " Test Array Structure (C,P,D) = " << C << ", " << P << ", " << D << "\n"
        << "===============================================================================\n";

      Kokkos::Timer timer;
      double t_without_subview[20] = {}, t_with_subview[20] = {};
      int errorFlag = 0, itest = 0;

      *verboseStream
        << "\n"
        << "===============================================================================\n"
        << "ArrayTools clone \n"
        << "===============================================================================\n";
      
      try {
        typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> ViewType;
        const auto loopSize = C*P;
        Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

        *verboseStream << "TEST " << itest << ":  input (P), output (C,P) \n" << "\n";        
        {
          ViewType in("inDynRankView", P), out("outDynRankView", C,P);          
          Kokkos::deep_copy(in, 1.0);                          
          {
            *verboseStream << " -> with subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,-1> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_with_subview[itest] = timer.seconds();
          }
          {
            *verboseStream << " -> without subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,0> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_without_subview[itest] = timer.seconds();
          }
        }
        ++itest;
        *verboseStream << "TEST " << itest << ":  input (P,D), output (C,P,D) \n" << "\n";        
        {
          ViewType in("inDynRankView", P,D), out("outDynRankView", C,P,D);          
          Kokkos::deep_copy(in, 1.0);                          
          {
            *verboseStream << " -> with subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,-1> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_with_subview[itest] = timer.seconds();
          }
          {
            *verboseStream << " -> without subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,1> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_without_subview[itest] = timer.seconds();
          }
        }
        ++itest;

        *verboseStream << "TEST " << itest << ":  input (P,D,D), output (C,P,D,D) \n" << "\n";        
        {
          ViewType in("inDynRankView", P,D,D), out("outDynRankView", C,P,D,D);          
          Kokkos::deep_copy(in, 1.0);                          
          {
            *verboseStream << " -> with subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,-1> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_with_subview[itest] = timer.seconds();
          }
          {
            *verboseStream << " -> without subview \n";
            DeviceSpaceType().fence();
            timer.reset();
            
            typedef F_clone<ViewType,ViewType,2> FunctorType;
            for (ordinal_type i=0;i<nworkset;++i) {
              Kokkos::parallel_for( policy, FunctorType(out, in) );
            }
            DeviceSpaceType().fence();
            t_without_subview[itest] = timer.seconds();
          }
        }
        ++itest;

      } catch (std::exception err) {
        *verboseStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *verboseStream << err.what() << '\n';
        *verboseStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      std::cout 
        << "===============================================================================\n";

      for (auto i=0;i<itest;++i) {
        std::cout 
          << "TEST " << i 
          << ": t_direct = " << std::setw(8) << t_without_subview[i] 
          << ", t_subview= " << std::setw(8) << t_with_subview[i] 
          << ", ratio(d/s) = " << std::setw(8) << (t_without_subview[i]/t_with_subview[i]) << "\n";
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
