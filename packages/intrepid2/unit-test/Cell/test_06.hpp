// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

        double buf[10]; Kokkos::View<double*,Kokkos::AnonymousSpace> val(&buf[0], N); // N
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

      typedef typename ExecSpaceType::array_layout DeviceArrayLayout;
        
      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) Serial interface tests                                               |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";
  
      const ValueType tol = tolerence<ValueType>()*100.0;

      int errorFlag = 0;
      
      try {
        
        {
          *outStream
            << "\n"
            << "===============================================================================\n"
            << "| Test 1: reference quad c1 element mapped into 2D physical space:            |\n"
            << "===============================================================================\n\n";

          const ordinal_type C = 3, P = 5, N = 4, D = 2;
          
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,Kokkos::HostSpace> pts_on_ref_host("pts_on_ref_host", P, D);
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,Kokkos::HostSpace> workset_host("workset_host", C, N, D);
           
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
          auto a_pts_on_phy_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), a_pts_on_phy);
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
          //     double buf[N]; Kokkos::View<double*,Kokkos::AnonymousSpace> val(&buf[0], N); // N
          //     auto nodes = Kokkos::subview(workset, cl, Kokkos::ALL(), Kokkos::ALL()); // N,D
          //     for (ordinal_type i=0;i<P;++i) {
          //       auto pt_on_ref = Kokkos::subview(  pts_on_ref,     i, Kokkos::ALL()); // D
          //       auto pt_on_phy = Kokkos::subview(b_pts_on_phy, cl, i, Kokkos::ALL()); // D
          //       Impl::Basis_HGRAD_QUAD_C1_FEM::Serial<OPERATOR_VALUE>::getValues(val, pt_on_ref);
          //       Impl::CellTools::Serial::mapToPhysicalFrame(pt_on_phy, val, nodes);
          //     }
          //   });
          auto b_pts_on_phy_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), b_pts_on_phy);
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
          auto a_pts_on_ref_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), a_pts_on_ref);
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
          auto b_pts_on_ref_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), b_pts_on_ref);
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


          *outStream
            << "\n"
            << "===============================================================================\n"
            << "| Test 2: reference quad c1 element mapped into 3D physical space:            |\n"
            << "===============================================================================\n\n";


          Kokkos::DynRankView<ValueType,DeviceArrayLayout,Kokkos::HostSpace> workset3d_host("workset3d_host", C, N, D+1);

          //set x,y components
          for(size_t i=0; i<workset_host.extent(0);++i)
            for(size_t j=0; j<workset_host.extent(1);++j)
            for(int d=0; d<D;++d)
              workset3d_host(i,j,d) = workset_host(i,j,d);

          //set z component
          workset3d_host(0,0,2) = -1.0;
          workset3d_host(0,1,2) =  1.0;
          workset3d_host(0,2,2) =  4.0;
          workset3d_host(0,3,2) =  2.0;

          workset3d_host(1,0,2) =  1.0;
          workset3d_host(1,1,2) =  0.0;
          workset3d_host(1,2,2) =  0.0;
          workset3d_host(1,3,2) =  0.0;

          workset3d_host(2,0,2) =  0.0;
          workset3d_host(2,1,2) =  0.0;
          workset3d_host(2,2,2) =  5.0;
          workset3d_host(2,3,2) =  0.0;

          auto workset3d = Kokkos::create_mirror_view(typename DeviceType::memory_space(), workset3d_host);
          Kokkos::deep_copy(workset3d, workset3d_host);

          ///
          /// mapToPhysicalFrame
          ///

          // ** compute via impl interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> b_pts3d_on_phy("b_pts_on_phy", C, P, D+1);


          {
            typedef F_mapToPhysicalFrame<decltype(b_pts3d_on_phy),
                                         decltype(pts_on_ref),
                                         decltype(workset3d)> FunctorType;
            Kokkos::parallel_for(policy, FunctorType(b_pts3d_on_phy, pts_on_ref, workset3d));
          }

          ///
          /// mapToReferenceFrame
          ///

          // ** compute via impl interface
          Kokkos::DynRankView<ValueType,DeviceArrayLayout,DeviceType> b_pts3d_on_ref("b_pts3d_on_ref", C, P, D);

          {
            typedef F_mapToReferenceFrame<decltype(b_pts3d_on_ref),
                                          decltype(b_pts3d_on_phy),
                                          decltype(workset3d)> FunctorType;
            Kokkos::parallel_for(policy, FunctorType(b_pts3d_on_ref, b_pts3d_on_phy, workset3d));
          }

          auto b_pts3d_on_ref_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), b_pts3d_on_ref);
          Kokkos::deep_copy(b_pts3d_on_ref_host, b_pts3d_on_ref);

          // ** compare
          {
            for (ordinal_type cl=0;cl<C;++cl)
              for (ordinal_type i=0;i<P;++i)
                for (ordinal_type j=0;j<D;++j) {
                  const double diff = std::abs(pts_on_ref_host(i, j) - b_pts3d_on_ref_host(cl, i, j));
                  if (diff > tol) {
                    *outStream << "Error (3d physical space) : mapToReferenceFrame (impl version) at ("
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
    

