// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#include "Intrepid2_CellTools_Serial.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  

    template<typename OutputViewType,
             typename inputViewType,
             typename worksetViewType>
    struct F_mapToPhysicalFrame {
      OutputViewType _output;
      inputViewType _input;
      worksetViewType _workset;

      KOKKOS_INLINE_FUNCTION
      F_mapToPhysicalFrame(OutputViewType output_,
                           inputViewType input_,
                           worksetViewType workset_) 
        : _output(output_), 
          _input(input_), 
          _workset(workset_) {}

      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type cl) const {
        const ordinal_type P = _output.extent(1);
        const ordinal_type N = _workset.extent(1);

        double buf[10]; Kokkos::View<double*,Kokkos::Impl::ActiveExecutionMemorySpace> val(&buf[0], N); // N
        auto nodes = Kokkos::subview(_workset, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D

        for (ordinal_type i=0;i<P;++i) {
          auto in  = Kokkos::subview(_input,      i, Kokkos::ALL()); // D
          auto out = Kokkos::subview(_output, cl, i, Kokkos::ALL()); // D
          Impl::Basis_HGRAD_QUAD_C1_FEM::Serial<OPERATOR_VALUE>::getValues(val, in);
          Impl::CellTools::Serial::mapToPhysicalFrame(out, val, nodes);
        }
      }
    };

    template<typename OutputViewType,
             typename inputViewType,
             typename worksetViewType>      
    struct F_mapToReferenceFrame {
      OutputViewType _output;
      inputViewType _input;
      worksetViewType _workset;


      KOKKOS_INLINE_FUNCTION
      F_mapToReferenceFrame(OutputViewType output_,
                            inputViewType input_,
                            worksetViewType workset_) 
        : _output(output_), 
          _input(input_), 
          _workset(workset_) {}
      
      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type cl) const {
        const ordinal_type P = _output.extent(1);

        auto nodes = Kokkos::subview(_workset, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D
        for (ordinal_type i=0;i<P;++i) {
          auto in  = Kokkos::subview(_input,  cl, i, Kokkos::ALL()); // D
          auto out = Kokkos::subview(_output, cl, i, Kokkos::ALL()); // D
          
          Impl::CellTools::Serial
            ::mapToReferenceFrame<Impl::Basis_HGRAD_QUAD_C1_FEM>(out, in, nodes);
        }
        
      }
    };
        
    template<typename ValueType, typename DeviceType>
    int CellTools_Test06(const bool verbose) {

      using ExecSpaceType = typename DeviceType::execution_space;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space HostSpaceType ;

      typedef typename ExecSpaceType::array_layout DeviceArrayLayout;
        

      *outStream << "DeviceSpace::  ";   ExecSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
      
      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) Serial interface tests                                               |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";
  
      const ValueType tol = tolerence()*100.0;

      int errorFlag = 0;
      
      try {
        
        *outStream
          << "\n"
          << "===============================================================================\n" 
          << "| Test 1: quad c1 element test :                                              |\n" 
          << "===============================================================================\n\n";

        {
          const ordinal_type C = 3, P = 5, N = 4, D = 2;
          
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,HostSpaceType> pts_on_ref_host("pts_on_ref_host", P, D);
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,HostSpaceType> workset_host("workset_host", C, N, D);
           
          pts_on_ref_host(0, 0) =  0.0;           pts_on_ref_host(0, 1) =  0.0; 
          pts_on_ref_host(1, 0) =  0.3;           pts_on_ref_host(1, 1) =  0.2; 
          pts_on_ref_host(2, 0) =  0.2;           pts_on_ref_host(2, 1) =  0.4; 
          pts_on_ref_host(3, 0) = -0.1;           pts_on_ref_host(3, 1) =  0.1; 
          pts_on_ref_host(4, 0) = -0.2;           pts_on_ref_host(4, 1) =  0.3; 
          /*
          workset_host(0, 0, 0) = -1.0;           workset_host(0, 0, 1) = -1.0; 
          workset_host(0, 1, 0) =  1.0;           workset_host(0, 1, 1) = -1.0; 
          workset_host(0, 2, 0) =  1.0;           workset_host(0, 2, 1) =  1.0; 
          workset_host(0, 3, 0) = -1.0;           workset_host(0, 3, 1) =  1.0; 

          workset_host(1, 0, 0) = -1.0;           workset_host(1, 0, 1) = -1.0; 
          workset_host(1, 1, 0) =  1.0;           workset_host(1, 1, 1) = -1.0; 
          workset_host(1, 2, 0) =  1.0;           workset_host(1, 2, 1) =  1.0; 
          workset_host(1, 3, 0) = -1.0;           workset_host(1, 3, 1) =  1.0; 

          workset_host(2, 0, 0) = -1.0;           workset_host(2, 0, 1) = -1.0; 
          workset_host(2, 1, 0) =  1.0;           workset_host(2, 1, 1) = -1.0; 
          workset_host(2, 2, 0) =  1.0;           workset_host(2, 2, 1) =  1.0; 
          workset_host(2, 3, 0) = -1.0;           workset_host(2, 3, 1) =  1.0; 
          */
          workset_host(0, 0, 0) =  0.0;           workset_host(0, 0, 1) =  0.0; 
          workset_host(0, 1, 0) =  2.0;           workset_host(0, 1, 1) =  0.0; 
          workset_host(0, 2, 0) =  2.0;           workset_host(0, 2, 1) =  2.0; 
          workset_host(0, 3, 0) =  0.0;           workset_host(0, 3, 1) =  2.0; 

          workset_host(1, 0, 0) = -3.0;           workset_host(1, 0, 1) = -4.0; 
          workset_host(1, 1, 0) =  2.0;           workset_host(1, 1, 1) =  0.0; 
          workset_host(1, 2, 0) =  1.0;           workset_host(1, 2, 1) =  2.0; 
          workset_host(1, 3, 0) = -1.0;           workset_host(1, 3, 1) =  1.0; 

          workset_host(2, 0, 0) = -0.5;           workset_host(2, 0, 1) = -0.5; 
          workset_host(2, 1, 0) =  2.0;           workset_host(2, 1, 1) = -0.5; 
          workset_host(2, 2, 0) =  3.0;           workset_host(2, 2, 1) =  2.0; 
          workset_host(2, 3, 0) =  1.0;           workset_host(2, 3, 1) =  2.0; 

          auto pts_on_ref = Kokkos::create_mirror_view(typename DeviceType::memory_space(), pts_on_ref_host);
          Kokkos::deep_copy(pts_on_ref, pts_on_ref_host);

          auto workset = Kokkos::create_mirror_view(typename DeviceType::memory_space(), workset_host);
          Kokkos::deep_copy(workset, workset_host);

          Kokkos::RangePolicy<ExecSpaceType> policy(0, C);
          
          ///
          /// mapToPhysicalFrame 
          ///

          // ** compute via front interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> a_pts_on_phy("a_pts_on_phy", C, P, D);
          {
            CellTools<DeviceType>
              ::mapToPhysicalFrame(a_pts_on_phy, pts_on_ref, workset, 
                                   shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >()));
          }
          auto a_pts_on_phy_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a_pts_on_phy);
          Kokkos::deep_copy(a_pts_on_phy_host, a_pts_on_phy);
          
          // ** compute via impl interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> b_pts_on_phy("b_pts_on_phy", C, P, D);


          {
            typedef F_mapToPhysicalFrame<decltype(b_pts_on_phy),
                                         decltype(pts_on_ref),
                                         decltype(workset)> FunctorType;
            Kokkos::parallel_for(policy, FunctorType(b_pts_on_phy, pts_on_ref, workset));
          }
          // I love lambda...
          // Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type cl) {
          //     double buf[N]; Kokkos::View<double*,Kokkos::Impl::ActiveExecutionMemorySpace> val(&buf[0], N); // N
          //     auto nodes = Kokkos::subview(workset, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D
          //     for (ordinal_type i=0;i<P;++i) {
          //       auto pt_on_ref = Kokkos::subview(  pts_on_ref,     i, Kokkos::ALL()); // D
          //       auto pt_on_phy = Kokkos::subview(b_pts_on_phy, cl, i, Kokkos::ALL()); // D
          //       Impl::Basis_HGRAD_QUAD_C1_FEM::Serial<OPERATOR_VALUE>::getValues(val, pt_on_ref);
          //       Impl::CellTools::Serial::mapToPhysicalFrame(pt_on_phy, val, nodes);
          //     }
          //   });
          auto b_pts_on_phy_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b_pts_on_phy);
          Kokkos::deep_copy(b_pts_on_phy_host, b_pts_on_phy);

          // ** compare
          {
            for (ordinal_type cl=0;cl<C;++cl)
              for (ordinal_type i=0;i<P;++i) 
                for (ordinal_type j=0;j<D;++j) {
                  const double diff = std::abs(a_pts_on_phy_host(cl, i, j) - b_pts_on_phy_host(cl, i, j));
                  if (diff > tol) {
                    *outStream << "Error : mapToPhysicalFrame at (" 
                               << cl << "," << i << "," << j 
                               << ") with diff = " << diff << "\n"; 
                    errorFlag++;
                  }
                }
          }

          ///
          /// mapToReferenceFrame 
          ///
          auto pts_on_phy = a_pts_on_phy;

          // ** compute via front interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> a_pts_on_ref("a_pts_on_ref", C, P, D);
          {
            CellTools<DeviceType>
              ::mapToReferenceFrame(a_pts_on_ref, pts_on_phy, workset, 
                                   shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >()));
          }
          auto a_pts_on_ref_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a_pts_on_ref);
          Kokkos::deep_copy(a_pts_on_ref_host, a_pts_on_ref);

          // ** compare
          {
            for (ordinal_type cl=0;cl<C;++cl)
              for (ordinal_type i=0;i<P;++i) 
                for (ordinal_type j=0;j<D;++j) {
                  const double diff = std::abs(pts_on_ref_host(i, j) - a_pts_on_ref_host(cl, i, j));
                  if (diff > tol) {
                    *outStream << "Error : mapToReferenceFrame (front version) at (" 
                               << cl << "," << i << "," << j 
                               << ") with diff = " << diff << "\n"; 
                    errorFlag++;
                  }
                }
          }
                    
          // ** compute via impl interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> b_pts_on_ref("b_pts_on_ref", C, P, D);

          {
            typedef F_mapToReferenceFrame<decltype(b_pts_on_ref),
                                          decltype(pts_on_phy),
                                          decltype(workset)> FunctorType;
            Kokkos::parallel_for(policy, FunctorType(b_pts_on_ref, pts_on_phy, workset));
          }
          // Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type cl) {
          //     auto nodes = Kokkos::subview(workset, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D
          //     for (ordinal_type i=0;i<P;++i) {
          //       auto pt_on_phy = Kokkos::subview(  pts_on_phy, cl, i, Kokkos::ALL()); // D
          //       auto pt_on_ref = Kokkos::subview(b_pts_on_ref, cl, i, Kokkos::ALL()); // D

          //       Impl::CellTools::Serial
          //         ::mapToReferenceFrame<Impl::Basis_HGRAD_QUAD_C1_FEM>(pt_on_ref, pt_on_phy, nodes);
          //     }
          //   });
          auto b_pts_on_ref_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b_pts_on_ref);
          Kokkos::deep_copy(b_pts_on_ref_host, b_pts_on_ref);

          // ** compare          
          {
            for (ordinal_type cl=0;cl<C;++cl)
              for (ordinal_type i=0;i<P;++i) 
                for (ordinal_type j=0;j<D;++j) {
                  const double diff = std::abs(pts_on_ref_host(i, j) - b_pts_on_ref_host(cl, i, j));
                  if (diff > tol) {
                    *outStream << "Error : mapToReferenceFrame (impl version) at (" 
                               << cl << "," << i << "," << j 
                               << ") with diff = " << diff << "\n"; 
                    errorFlag++;
                  }
                }
          }

          
        }
      } catch (std::logic_error &err) {
        //============================================================================================//
        // Wrap up test: check if the test broke down unexpectedly due to an exception                //
        //============================================================================================//
        *outStream << err.what() << "\n";
        errorFlag = -1000;
      }
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  }
}
    

