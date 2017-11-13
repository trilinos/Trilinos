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

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::G_TET_COMP12_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, J. Ostien and Kyungjoo Kim.
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_TET_COMP12_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include <random>

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename ValueType, typename DeviceSpaceType>
    int HGRAD_TET_COMP12_FEM_Test01(const bool verbose) {
      
      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing
      
      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

      *outStream                                                        
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|             Unit Test (Basis_HGRAD_TET_COMP12_FEM)                          |\n"
        << "|                                                                             |\n"
        << "|     1) Evaluation of Basis Function Values                                  |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
        << "|                      Jake Ostien   (jtostie@sandia.gov),                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_TET_COMP12_FEM<DeviceSpaceType,outputValueType,pointValueType> tetBasis;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: correctness of basis function values                                |\n"
        << "===============================================================================\n";
      
      // output precision
      outStream -> precision(20);
      
      // VALUE: Each row gives the 10 correct basis set values at an evaluation point
      const ValueType nodalBasisValues[] = {
        // first 4 vertices
        1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        // second 6 vertices
        0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 1.0
      };

      const ValueType pointBasisValues[] = {
        // pt 0 {0.25, 0.25, 0.25}
        0.0, 0.0, 0.0, 0.0,  1./6., 1./6., 1./6., 1./6., 1./6., 1./6.,
        // pt 1 {0.5, 1/6, 1/6}
        0.0, 0.0, 0.0, 0.0,  1./3., 1./3., 0.0, 0.0, 1./3., 0.0,
        // pt 2 {1/6, 0.5, 0.1/6}
        0.0, 0.0, 0.0, 0.0,  0.0, 1./3., 1./3., 0.0, 0.0, 1./3.,
        // pt 3 {1/6, 1/6, 0.5}
        0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1./3., 1./3., 1./3.,
        // pt 4 {1/6, 1/6, 1/6}
        0.0, 0.0, 0.0, 0.0,  1./3., 0.0, 1./3., 1./3., 0.0, 0.0,
        // pt 5
        0.170820393249936908922752100619382870632, 0.0, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0,
        // pt 6
        0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.0,
        // pt 7
        0.0, 0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0, 0.276393202250021030359082633126872376456,
        // pt 8
        0.0, 0.0, 0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456,
      };

      // GRAD and D1: each row gives the 3x10 correct values of the gradients of the 10 basis functions
      const ValueType pointBasisGrads[] = {
        // point 0
        -1./4.,   -1./4.,   -1./4.,
        1./4.,      0.0,      0.0, 
        0.0,    1./4.,      0.0,   
        0.0,      0.0,    1./4.,
        0.0,   -3./4.,   -3./4.,
        3./4.,    3./4.,      0.0,
        -3./4.,      0.0,   -3./4.,
        -3./4.,   -3./4.,      0.0,
        3./4.,      0.0,    3./4.,
        0.0,    3./4.,    3./4.,

        // point 1
        -1./24.,  -1./24.,  -1./24.,
        7./8.,      0.0,      0.0,
        0.0,   1./24.,      0.0,
        0.0,      0.0,   1./24.,
        -35./36., -19./12., -19./12.,
        11./18.,  19./12.,      0.0,
        -17./36.,      0.0,   -1./3.,
        -17./36.,   -1./3.,      0.0,
        11./18.,      0.0,  19./12.,
        -5./36.,    1./3.,    1./3.,

        // point 2
        -1./24.,  -1./24.,  -1./24.,
        1./24.,      0.0,      0.0,
        0.0,    7./8.,      0.0,
        0.0,      0.0,   1./24.,
        0.0, -17./36.,   -1./3.,
        19./12.,  11./18.,      0.0,
        -19./12., -35./36., -19./12.,
        -1./3., -17./36.,      0.0,
        1./3.,  -5./36.,    1./3.,
        0.0,  11./18.,  19./12.,

        // point 3
        -1./24.,  -1./24.,  -1./24.,
        1./24.,      0.0,      0.0,
        0.0,   1./24.,      0.0,
        0.0,      0.0,    7./8.,
        0.0,   -1./3., -17./36.,
        1./3.,    1./3.,  -5./36.,
        -1./3.,      0.0, -17./36.,
        -19./12., -19./12., -35./36.,
        19./12.,      0.0,  11./18.,
        0.0,  19./12.,  11./18.,

        // point 4
        -7./8.,   -7./8.,   -7./8.,
        1./24.,      0.0,      0.0,
        0.0,   1./24.,      0.0,
        0.0,      0.0,   1./24.,
        35./36., -11./18., -11./18.,
        17./36.,  17./36.,   5./36.,
        -11./18.,  35./36., -11./18.,
        -11./18., -11./18.,  35./36.,
        17./36.,   5./36.,  17./36.,
        5./36.,  17./36.,  17./36.,

        // point 5
        -1.088525491562421136153440125774228588290, -1.088525491562421136153440125774228588290, -1.088525491562421136153440125774228588290,
        -0.029508497187473712051146708591409529430, 0.0, 0.0,
        0.0, -0.029508497187473712051146708591409529430, 0.0,
        0.0, 0.0, -0.029508497187473712051146708591409529430,
        1.30437298687487732290535130675991113734, -0.563661001875017525299235527605726980380, -0.563661001875017525299235527605726980380,
        0.377322003750035050598471055211453960760, 0.377322003750035050598471055211453960760, 0.186338998124982474700764472394273019620,
        -0.563661001875017525299235527605726980380, 1.30437298687487732290535130675991113734, -0.563661001875017525299235527605726980380,
        -0.563661001875017525299235527605726980380, -0.563661001875017525299235527605726980380, 1.30437298687487732290535130675991113734,
        0.377322003750035050598471055211453960760, 0.186338998124982474700764472394273019620, 0.377322003750035050598471055211453960760,
        0.186338998124982474700764472394273019620, 0.377322003750035050598471055211453960760, 0.377322003750035050598471055211453960760,

        // point 6
        0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430,
        1.088525491562421136153440125774228588290, 0.0, 0.0,
        0.0, -0.029508497187473712051146708591409529430, 0.0,
        0.0, 0.0, -0.029508497187473712051146708591409529430,
        -1.30437298687487732290535130675991113734, -1.868033988749894848204586834365638117720, -1.868033988749894848204586834365638117720,
        0.563661001875017525299235527605726980380, 1.868033988749894848204586834365638117720, 0.0,
        -0.377322003750035050598471055211453960760, 0.0, -0.190983005625052575897706582817180941140,
        -0.377322003750035050598471055211453960760, -0.190983005625052575897706582817180941140, 0.0,
        0.563661001875017525299235527605726980380, 0.0, 1.868033988749894848204586834365638117720,
        -0.186338998124982474700764472394273019620, 0.190983005625052575897706582817180941140, 0.19098300562505257589770658281718094114,

        // point 7
        0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430,
        -0.029508497187473712051146708591409529430, 0.0, 0.0,          
        0.0, 1.088525491562421136153440125774228588290, 0.0,           
        0.0, 0.0, -0.029508497187473712051146708591409529430,          
        0.0, -0.377322003750035050598471055211453960760, -0.190983005625052575897706582817180941140,
        1.868033988749894848204586834365638117720, 0.563661001875017525299235527605726980380, 0.0,
        -1.868033988749894848204586834365638117720, -1.30437298687487732290535130675991113734, -1.868033988749894848204586834365638117720,
        -0.190983005625052575897706582817180941140, -0.377322003750035050598471055211453960760, 0.0,
        0.190983005625052575897706582817180941140, -0.186338998124982474700764472394273019620, 0.190983005625052575897706582817180941140,
        0.0, 0.563661001875017525299235527605726980380, 1.868033988749894848204586834365638117720,

        // point 8
        0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430,
        -0.029508497187473712051146708591409529430, 0.0, 0.0,
        0.0, -0.029508497187473712051146708591409529430, 0.0,
        0.0, 0.0, 1.088525491562421136153440125774228588290,
        0.0, -0.190983005625052575897706582817180941140, -0.377322003750035050598471055211453960760,
        0.190983005625052575897706582817180941140, 0.190983005625052575897706582817180941140, -0.186338998124982474700764472394273019620,
        -0.190983005625052575897706582817180941140, 0.0, -0.377322003750035050598471055211453960760,
        -1.868033988749894848204586834365638117720, -1.868033988749894848204586834365638117720, -1.30437298687487732290535130675991113734,
        1.868033988749894848204586834365638117720, 0.0, 0.563661001875017525299235527605726980380,
        0.0, 1.868033988749894848204586834365638117720, 0.563661001875017525299235527605726980380,
      };
      
      try{
        DynRankViewHost ConstructWithLabel(tetNodesHost, 10, 3);
        
        tetNodesHost(0,0) = 0.0;  tetNodesHost(0,1) = 0.0;  tetNodesHost(0,2) = 0.0;
        tetNodesHost(1,0) = 1.0;  tetNodesHost(1,1) = 0.0;  tetNodesHost(1,2) = 0.0;
        tetNodesHost(2,0) = 0.0;  tetNodesHost(2,1) = 1.0;  tetNodesHost(2,2) = 0.0;
        tetNodesHost(3,0) = 0.0;  tetNodesHost(3,1) = 0.0;  tetNodesHost(3,2) = 1.0;
        tetNodesHost(4,0) = 0.5;  tetNodesHost(4,1) = 0.0;  tetNodesHost(4,2) = 0.0;
        tetNodesHost(5,0) = 0.5;  tetNodesHost(5,1) = 0.5;  tetNodesHost(5,2) = 0.0;
        tetNodesHost(6,0) = 0.0;  tetNodesHost(6,1) = 0.5;  tetNodesHost(6,2) = 0.0;
        tetNodesHost(7,0) = 0.0;  tetNodesHost(7,1) = 0.0;  tetNodesHost(7,2) = 0.5;
        tetNodesHost(8,0) = 0.5;  tetNodesHost(8,1) = 0.0;  tetNodesHost(8,2) = 0.5;
        tetNodesHost(9,0) = 0.0;  tetNodesHost(9,1) = 0.5;  tetNodesHost(9,2) = 0.5;
        
        auto tetNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), tetNodesHost);
        Kokkos::deep_copy(tetNodes, tetNodesHost);

        DynRankViewHost ConstructWithLabel(tetPointsHost, 9, 3);
        
        // from the 5 point integration
        tetPointsHost(0,0) = 0.25;     tetPointsHost(0,1) = 0.25;     tetPointsHost(0,2) = 0.25;
        tetPointsHost(1,0) = 0.5;      tetPointsHost(1,1) = (1./6.);  tetPointsHost(1,2) = (1./6.);
        tetPointsHost(2,0) = (1./6.);  tetPointsHost(2,1) = 0.5;      tetPointsHost(2,2) = (1./6.);
        tetPointsHost(3,0) = (1./6.);  tetPointsHost(3,1) = (1./6.);  tetPointsHost(3,2) = 0.5;
        tetPointsHost(4,0) = (1./6.);  tetPointsHost(4,1) = (1./6.);  tetPointsHost(4,2) = (1./6.);
        
        // from the 4 point integration
        tetPointsHost(5,0) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(5,1) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(5,2) = 0.1381966011250105151795413165634361882280;
        
        tetPointsHost(6,0) = 0.5854101966249684544613760503096914353161;
        tetPointsHost(6,1) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(6,2) = 0.1381966011250105151795413165634361882280;
        
        tetPointsHost(7,0) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(7,1) = 0.5854101966249684544613760503096914353161;
        tetPointsHost(7,2) = 0.1381966011250105151795413165634361882280;
        
        tetPointsHost(8,0) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(8,1) = 0.1381966011250105151795413165634361882280;
        tetPointsHost(8,2) = 0.5854101966249684544613760503096914353161;

        auto tetPoints = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), tetPointsHost);
        Kokkos::deep_copy(tetPoints, tetPointsHost);

        // Dimensions for the output arrays:
        const ordinal_type numFields = tetBasis.getCardinality();
        const ordinal_type numNodes  = tetNodes.dimension(0);
        const ordinal_type spaceDim  = tetBasis.getBaseCellTopology().getDimension();
    
        // Check VALUE of basis functions at nodes: resize vals to rank-2 container:\n";
        {
          *outStream << " check VALUE of basis functions at nodes\n";
          DynRankView vals = DynRankView("vals", numFields, numNodes);
          tetBasis.getValues(vals, tetNodes, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);

          for (ordinal_type i=0;i<numFields;++i) {
            for (ordinal_type j=0;j<numNodes;++j) {
              const ordinal_type l =  i + j * numFields;
              if (std::abs(vals_host(i,j) - nodalBasisValues[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << vals_host(i,j)
                           << " but reference value: " << nodalBasisValues[l] << "\n";
              }
            }
          }
        }

        const ordinal_type numPoints = tetPoints.dimension(0);

        // Check VALUE of basis functions at points: resize vals to rank-2 container:\n";
        {
          *outStream << " check VALUE of basis functions at points\n";
          DynRankView vals = DynRankView("vals", numFields, numPoints);
          tetBasis.getValues(vals, tetPoints, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);

          for (ordinal_type i=0;i<numFields;++i) {
            for (ordinal_type j=0;j<numPoints;++j) {
              const ordinal_type l =  i + j * numFields;
              if (std::abs(vals_host(i,j) - pointBasisValues[l]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << vals_host(i,j)
                           << " but reference value: " << pointBasisValues[l] << "\n";
              }
            }
          }
        }
        
        // Check VALUE of basis functions at random points: resize vals to rank-2 container:\n";
        {
          *outStream << " check VALUE of basis functions at random points\n";
          const ordinal_type numRandomPoints = 16384;

          DynRankViewHost tetRandomPointsHost = DynRankViewHost("tetRandomPointsHost", numRandomPoints, 3);
          {
            ordinal_type point = 0;

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0, 1);
            while (point < numRandomPoints) {
              const ValueType r = dis(gen), s = dis(gen), t = dis(gen);
              if (r + s + t > 1.0) {
                // do nothing
              } else {
                tetRandomPointsHost(point, 0) = r;
                tetRandomPointsHost(point, 1) = s;
                tetRandomPointsHost(point, 2) = t;
                ++point;
              }
            }
          }
          
          auto tetRandomPoints = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), tetRandomPointsHost);
          Kokkos::deep_copy(tetRandomPoints, tetRandomPointsHost);

          DynRankView vals = DynRankView("vals", numFields, numRandomPoints);
        
          tetBasis.getValues(vals, tetRandomPoints, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
        
          for (ordinal_type j=0;j<numRandomPoints;++j) {
            ValueType sum = 0.0;
            for (ordinal_type i=0;i<numFields;++i)
              sum += vals_host(i,j);

            if (std::abs(sum - 1.0) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              
              // Just indicate something bad happened
              *outStream << " Composite tet basis functions";
              *outStream << " are not summing to 1.0\n";
              *outStream << " sum : " << sum << "\n";
            }
          }
        }
        
    
        // Check GRAD of basis functions at points: resize vals to rank-3 container:\n";
        {
          DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
          tetBasis.getValues(vals, tetPoints, OPERATOR_GRAD);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i=0;i<numFields;++i) {
            for (ordinal_type j=0;j<numPoints;++j) {
              for (ordinal_type k=0;k<spaceDim;++k) {
                const ordinal_type l = k + i * spaceDim + j * spaceDim * numFields;
                if (std::abs(vals_host(i,j,k) - pointBasisGrads[l]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  
                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed grad component: " << vals_host(i,j,k)
                             << " but reference grad component: " << pointBasisGrads[l] << "\n";
                }
              }
            }
          }
        }
      } catch (std::logic_error err) {
        *outStream << err.what() << "\n\n";
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
