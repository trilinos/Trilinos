// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_02.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_LINE_C1_FEM class.
    \author Created by Kyungjoo Kim, Mauro Perego
 */


#include "Intrepid2_config.h"
#include "Kokkos_Random.hpp"
#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HGRAD_LINE_C1_FEM_Test02(const bool verbose) {

      //! Setup test output stream.
      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "HGRAD_LINE_C1_FEM, Test 2", {}
      );

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| Testing Team-level Implemntation of getValues                               |\n"
      << "===============================================================================\n";

      using DeviceSpaceType = typename DeviceType::execution_space;
      Kokkos::print_configuration(std::cout, false);

      int errorFlag = 0;

      try { 
        using BasisType = Basis_HGRAD_LINE_C1_FEM<DeviceType,OutValueType,PointValueType>;
        auto basisPtr = Teuchos::rcp(new BasisType());
        
        // problem setup 
        //   let's say we want to evaluate 1000 points in parallel. output values are stored in outputValuesA and B.
        //   A is compuated via serial interface and B is computed with top-level interface.
        const int ncells = 5, npts = 10, ndim = 1;
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesA, ncells, basisPtr->getCardinality(), npts);
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesB, basisPtr->getCardinality(), npts);

        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputGradsA, ncells, basisPtr->getCardinality(), npts, ndim);
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputGradsB, basisPtr->getCardinality(), npts, ndim);
        Kokkos::DynRankView<PointValueType,DeviceType> ConstructWithLabelPointView(point, 1);
        
        using ScalarType = typename ScalarTraits<PointValueType>::scalar_type;
        
        Kokkos::View<ScalarType*,DeviceType> inputPointsViewToUseRandom("inputPoints", npts*ndim*get_dimension_scalar(point));
        auto vcprop = Kokkos::common_view_alloc_prop(point);
        Kokkos::DynRankView<PointValueType,DeviceType> inputPoints (Kokkos::view_wrap(inputPointsViewToUseRandom.data(), vcprop),  npts, ndim);
        
        // random values between (0,1)
        Kokkos::Random_XorShift64_Pool<DeviceType> random(13718);
        Kokkos::fill_random(inputPointsViewToUseRandom, random, 1.0);
        

        *outStream << "Computing values and gradients for " << ncells << " cells and " << npts << " points using team-level getValues function" <<std::endl;

        { // evaluation using parallel loop over cell
          auto basisPtr_device = copy_virtual_class_to_device<DeviceType,BasisType>(*basisPtr);
          auto basisRawPtr_device = basisPtr_device.get();

          int scratch_space_level =1;
          const int vectorSize = getVectorSizeForHierarchicalParallelism<PointValueType>();
          Kokkos::TeamPolicy<DeviceSpaceType> teamPolicy(ncells, Kokkos::AUTO,vectorSize);

          { //compute values
            auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                auto valsACell = Kokkos::subview(outputValuesA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                basisRawPtr_device->getValues(valsACell, inputPoints, OPERATOR_VALUE, team_member, team_member.team_scratch(scratch_space_level));
            };              
            
            //Get the required size of the scratch space per team and per thread.
            int perThreadSpaceSize(0), perTeamSpaceSize(0);
            basisPtr->getScratchSpaceSize(perTeamSpaceSize,perThreadSpaceSize,inputPoints, OPERATOR_VALUE);
            teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerTeam(perTeamSpaceSize), Kokkos::PerThread(perThreadSpaceSize));

            Kokkos::parallel_for (teamPolicy,functor);
          }

          { //compute gradients
            auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                auto gradsACell = Kokkos::subview(outputGradsA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                basisRawPtr_device->getValues(gradsACell, inputPoints, OPERATOR_GRAD, team_member, team_member.team_scratch(scratch_space_level));
            };              
            
            //Get the required size of the scratch space per team and per thread.
            int perThreadSpaceSize(0), perTeamSpaceSize(0);
            basisPtr->getScratchSpaceSize(perTeamSpaceSize,perThreadSpaceSize,inputPoints, OPERATOR_GRAD);
            teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerTeam(perTeamSpaceSize), Kokkos::PerThread(perThreadSpaceSize));

            Kokkos::parallel_for (teamPolicy,functor);
          }
        }

        *outStream << "Computing values and gradients for " << npts << " points using high-level getValues function" <<std::endl;


        // evaluation using high level interface (on one cell)
        basisPtr->getValues(outputValuesB, inputPoints, OPERATOR_VALUE);
        basisPtr->getValues(outputGradsB, inputPoints, OPERATOR_GRAD);
        
        *outStream << "Comparing values and gradients on host" <<std::endl;
        { // compare values
          const auto outputValuesA_Host = Kokkos::create_mirror_view(outputValuesA); Kokkos::deep_copy(outputValuesA_Host, outputValuesA);
          const auto outputValuesB_Host = Kokkos::create_mirror_view(outputValuesB); Kokkos::deep_copy(outputValuesB_Host, outputValuesB);
          
          OutValueType diff = 0; 
          auto tol = epsilon<double>();
          for (size_t ic=0;ic<outputValuesA_Host.extent(0);++ic)
            for (size_t i=0;i<outputValuesA_Host.extent(1);++i)
              for (size_t j=0;j<outputValuesA_Host.extent(2);++j) {
                diff = std::abs(outputValuesB_Host(i,j) - outputValuesA_Host(ic,i,j));
                if (diff > tol) {
                  ++errorFlag;
                  std::cout << ", ic: " << ic << ", i: " << i << ", j: " << j 
                            << ", val A: " << outputValuesA_Host(ic,i,j) 
                            << ", val B: " << outputValuesB_Host(i,j) 
                            << ", |diff|: " << diff
                            << ", tol: " << tol
                            << std::endl;
                }
              }
        }

        { 
          // compare grads
          const auto outputGradsA_Host = Kokkos::create_mirror_view(outputGradsA); Kokkos::deep_copy(outputGradsA_Host, outputGradsA);
          const auto outputGradsB_Host = Kokkos::create_mirror_view(outputGradsB); Kokkos::deep_copy(outputGradsB_Host, outputGradsB);
          
          OutValueType diff = 0;
          auto tol = epsilon<double>();
          for (size_t ic=0;ic<outputGradsA_Host.extent(0);++ic)
            for (size_t i=0;i<outputGradsA_Host.extent(1);++i)
              for (size_t j=0;j<outputGradsA_Host.extent(2);++j) {
                diff = 0;
                for (int d=0;d<ndim;++d)
                  diff += std::abs(outputGradsB_Host(i,j,d) - outputGradsA_Host(ic,i,j,d));
                if (diff > tol) {
                  ++errorFlag;
                  std::cout << ", ic: " << ic << ", i: " << i << ", j: " << j 
                            << ", grads A: [" << outputGradsA_Host(ic,i,j,0) << "]"
                            << ", grads B: [" << outputGradsB_Host(i,j,0)  <<"]"
                            << ", |diff|: " << diff
                            << ", tol: " << tol
                            << std::endl;
                }
              }
        }
      } catch (std::exception &err) {
        std::cout << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        std::cout << err.what() << '\n';
        std::cout << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      return errorFlag;
    }
  }
}
