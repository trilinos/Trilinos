// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_02.hpp
    \brief  Unit tests for the Intrepid2::HDIV_HEX_I1_FEM class.
    \author Created by Kyungjoo Kim, Mauro Perego
 */


#include "Intrepid2_config.h"
#include "Kokkos_Random.hpp"
#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    // This test evaluates the basis functions at a set of points on a batch of cells using the team-level getValues,
    // and compares the results with those obtained using the classic getValues function.
    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HDIV_HEX_I1_FEM_Test02(const bool verbose) {

      //! Setup test output stream.
      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "HDIV_HEX_I1_FEM, Test 2", {}
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
        using BasisType = Basis_HDIV_HEX_I1_FEM<DeviceType,OutValueType,PointValueType>;
        auto basisPtr = Teuchos::rcp(new BasisType());
        
        
        const int ncells = 5, npts = 10, ndim = 3;
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesA, ncells, basisPtr->getCardinality(), npts, ndim);
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesB, basisPtr->getCardinality(), npts, ndim);

        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputDivergencesA, ncells, basisPtr->getCardinality(), npts);
        Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputDivergencesB, basisPtr->getCardinality(), npts);

        Kokkos::DynRankView<PointValueType,DeviceType> ConstructWithLabelPointView(inputPoints, npts, ndim);

        { //randomly initialize inputPoints including derivatives for fad types
          auto inputPointsViewToUseRandom = as_scalar_1d_view(inputPoints);          
          Kokkos::Random_XorShift64_Pool<DeviceType> random(20260504); // random values between (0,1)
          Kokkos::fill_random(inputPointsViewToUseRandom, random, 0.0, 1.0);
        } 
        
        *outStream << "Computing values and divergences for " << ncells << " cells and " << npts << " points using team-level getValues function" <<std::endl;

        { // evaluation using parallel loop over cell
          auto basisPtr_device = copy_virtual_class_to_device<DeviceType,BasisType>(*basisPtr);
          auto basisRawPtr_device = basisPtr_device.get();

          int scratch_space_level =1;
          const int vectorSize = getVectorSizeForHierarchicalParallelism<PointValueType>();
          Kokkos::TeamPolicy<DeviceSpaceType> teamPolicy(ncells, Kokkos::AUTO,vectorSize);

          { //compute values
            auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                auto valsACell = Kokkos::subview(outputValuesA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                basisRawPtr_device->getValues(valsACell, inputPoints, OPERATOR_VALUE, team_member, scratch_space_level);
            };              
            
            //Get the required size of the scratch space per team and per thread.
            int perThreadSpaceSize(0);
            basisPtr->getScratchSpaceSize(perThreadSpaceSize,inputPoints, OPERATOR_VALUE);
            teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerThread(perThreadSpaceSize));

            Kokkos::parallel_for (teamPolicy,functor);
          }

          { //compute divergences
            auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                auto divergencesACell = Kokkos::subview(outputDivergencesA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                basisRawPtr_device->getValues(divergencesACell, inputPoints, OPERATOR_DIV, team_member, scratch_space_level);
            };              
            
            //Get the required size of the scratch space per team and per thread.
            int perThreadSpaceSize(0);
            basisPtr->getScratchSpaceSize(perThreadSpaceSize,inputPoints, OPERATOR_DIV);
            teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerThread(perThreadSpaceSize));

            Kokkos::parallel_for (teamPolicy,functor);
          }
        }

        *outStream << "Computing values and divergences for " << npts << " points using high-level getValues function" <<std::endl;


        // evaluation using high level interface (on one cell)
        basisPtr->getValues(outputValuesB, inputPoints, OPERATOR_VALUE);
        basisPtr->getValues(outputDivergencesB, inputPoints, OPERATOR_DIV);
        
        *outStream << "Comparing values and divergences on host" <<std::endl;
        { // compare values
          const auto outputValuesA_Host = Kokkos::create_mirror_view(outputValuesA); Kokkos::deep_copy(outputValuesA_Host, outputValuesA);
          const auto outputValuesB_Host = Kokkos::create_mirror_view(outputValuesB); Kokkos::deep_copy(outputValuesB_Host, outputValuesB);
          
          const auto tol = 100.0 * epsilon<double>();
          for (size_t ic=0;ic<outputValuesA_Host.extent(0);++ic)
            for (size_t i=0;i<outputValuesA_Host.extent(1);++i)
              for (size_t j=0;j<outputValuesA_Host.extent(2);++j) {
                auto maxBNorm = computeMaxNorm(outputValuesB_Host(i,j,0));
                auto diffNorm = computeMaxNorm(outputValuesB_Host(i,j,0) - outputValuesA_Host(ic,i,j,0));
                for (int d=1;d<ndim;++d) {
                  maxBNorm = std::max(maxBNorm, computeMaxNorm(outputValuesB_Host(i,j,d)));
                  diffNorm = std::max(diffNorm, computeMaxNorm(outputValuesB_Host(i,j,d)- outputValuesA_Host(ic,i,j,d)));
                }
                const auto diffRelNorm = diffNorm/std::max(1.0, maxBNorm);
                if (diffRelNorm > tol) {
                  ++errorFlag;
                  std::cout << ", ic: " << ic << ", i: " << i << ", j: " << j 
                            << ", val A: [" << outputValuesA_Host(ic,i,j,0) << ", " << outputValuesA_Host(ic,i,j,1) << ", " << outputValuesA_Host(ic,i,j,2) << "]"
                            << ", val B: [" << outputValuesB_Host(i,j,0) << ", " << outputValuesB_Host(i,j,1) << ", " << outputValuesB_Host(i,j,2) << "]"
                            << ", |rel diff|: " << diffRelNorm
                            << ", tol: " << tol
                            << std::endl;
                }
              }
        }

        { 
          // compare divergences
          const auto outputDivergencesA_Host = Kokkos::create_mirror_view(outputDivergencesA); Kokkos::deep_copy(outputDivergencesA_Host, outputDivergencesA);
          const auto outputDivergencesB_Host = Kokkos::create_mirror_view(outputDivergencesB); Kokkos::deep_copy(outputDivergencesB_Host, outputDivergencesB);
          
          const auto tol = 100.0 * epsilon<double>();
          for (size_t ic=0;ic<outputDivergencesA_Host.extent(0);++ic)
            for (size_t i=0;i<outputDivergencesA_Host.extent(1);++i)
              for (size_t j=0;j<outputDivergencesA_Host.extent(2);++j) {
                const auto valA = outputDivergencesA_Host(ic,i,j);
                const auto valB = outputDivergencesB_Host(i,j);
                const auto maxBNorm = computeMaxNorm(valB);
                const auto diffRelNorm = computeMaxNorm(valB - valA)/std::max(1.0, maxBNorm);
                if (diffRelNorm > tol) {
                  ++errorFlag;
                  std::cout << ", ic: " << ic << ", i: " << i << ", j: " << j 
                            << ", divergence A: " << valA
                            << ", divergence B: " << valB
                            << ", |rel diff|: " << diffRelNorm
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
