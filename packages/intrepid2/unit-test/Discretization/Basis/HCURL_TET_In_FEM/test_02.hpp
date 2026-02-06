// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_02.hpp
    \brief  Unit tests for the Intrepid2::HCURL_TET_In_FEM class.
    \author Created by Kyungjoo Kim, Mauro Perego
 */


#include "Intrepid2_config.h"
#include "Kokkos_Random.hpp"
#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "packages/intrepid2/unit-test/Discretization/Basis/Setup.hpp"

namespace Intrepid2 {

  namespace Test {

    // This test evaluates the basis functions at a set of points on a batch of cells using the team-level getValues,
    // and compares the results with those obtained using the classic getValues function.
    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HCURL_TET_In_FEM_Test02(const bool verbose) {

      //! Setup test output stream.
      Teuchos::RCP<std::ostream> outStream = setup_output_stream<DeviceType>(
        verbose, "HCURL_TET_In_FEM, Test 2", {}
      );

      *outStream
      << "\n"
      << "===============================================================================\n"
      << "| Testing Team-level Implemntation of getValues                               |\n"
      << "===============================================================================\n";

      using DeviceSpaceType = typename DeviceType::execution_space;
      Kokkos::print_configuration(std::cout, false);

      int errorFlag = 0;
      constexpr int maxOrder = 7;

      try { 
        for (int order=1;order <= maxOrder;++order) {
          using BasisType = Basis_HCURL_TET_In_FEM<DeviceType,OutValueType,PointValueType>;
          auto basisPtr = Teuchos::rcp(new BasisType(order));
          
          const int ncells = 5, npts = 10, ndim = 3;
          Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesA, ncells, basisPtr->getCardinality(), npts, ndim);
          Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputValuesB, basisPtr->getCardinality(), npts, ndim);

          Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputCurlsA, ncells, basisPtr->getCardinality(), npts, ndim);
          Kokkos::DynRankView<OutValueType,DeviceType> ConstructWithLabelOutView(outputCurlsB, basisPtr->getCardinality(), npts, ndim);

          Kokkos::DynRankView<PointValueType,DeviceType> ConstructWithLabelPointView(inputPoints, npts, ndim);

          using ScalarType = typename ScalarTraits<PointValueType>::scalar_type;
          Kokkos::View<ScalarType**,DeviceType> inputPointsViewToUseRandom("inputPoints", npts, ndim);

          // random values between (0,1)
          Kokkos::Random_XorShift64_Pool<DeviceType> random(20251125);
          Kokkos::fill_random(inputPointsViewToUseRandom, random, 0.0, 1.0);

          auto policy = Kokkos::MDRangePolicy<DeviceSpaceType,Kokkos::Rank<2>>({0,0},{npts,ndim});
          Kokkos::parallel_for("initialize view", policy, KOKKOS_LAMBDA (const int &i, const int &j) {inputPoints(i,j) = inputPointsViewToUseRandom(i,j);});
          

          *outStream << "Order: " << order << ": Computing values and curls for " << ncells << " cells and " << npts << " points using team-level getValues function" <<std::endl;

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
              
              // avoid using a team size larger than needed, to reduce allocated scrach space memory
              ordinal_type team_size = teamPolicy.team_size_recommended(functor, Kokkos::ParallelForTag());
              *outStream << "Max Recommended team size: " << team_size << ", Requested team size: " << npts <<std::endl;
              team_size = std::min(team_size, npts);
              teamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>(ncells, team_size,vectorSize);
              
              //Get the required size of the scratch space per team and per thread.
              int perThreadSpaceSize(0), perTeamSpaceSize(0);
              basisPtr->getScratchSpaceSize(perTeamSpaceSize,perThreadSpaceSize,inputPoints, OPERATOR_VALUE);
              teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerTeam(perTeamSpaceSize), Kokkos::PerThread(perThreadSpaceSize));

              Kokkos::parallel_for (teamPolicy,functor);
            }

            { //compute curls
              auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                  auto curlsACell = Kokkos::subview(outputCurlsA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                  basisRawPtr_device->getValues(curlsACell, inputPoints, OPERATOR_CURL, team_member, team_member.team_scratch(scratch_space_level));
              };      

              // avoid using a team size larger than needed, to reduce allocated scrach space memory
              ordinal_type team_size = teamPolicy.team_size_recommended(functor, Kokkos::ParallelForTag());
              *outStream << "Max Recommended team size: " << team_size << ", Requested team size: " << npts <<std::endl;
              team_size = std::min(team_size, npts);
              teamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>(ncells, team_size,vectorSize);        
              
              //Get the required size of the scratch space per team and per thread.
              int perThreadSpaceSize(0), perTeamSpaceSize(0);
              basisPtr->getScratchSpaceSize(perTeamSpaceSize,perThreadSpaceSize,inputPoints, OPERATOR_CURL);
              teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerTeam(perTeamSpaceSize), Kokkos::PerThread(perThreadSpaceSize));

              Kokkos::parallel_for (teamPolicy,functor);
            }
          }

          *outStream << "Order: " << order << ": Computing values and curls for " << npts << " points using high-level getValues function" <<std::endl;


          // evaluation using high level interface (on one cell)
          basisPtr->getValues(outputValuesB, inputPoints, OPERATOR_VALUE);
          basisPtr->getValues(outputCurlsB, inputPoints, OPERATOR_CURL);
          
          *outStream << "Order: " << order << ": Comparing values and curls on host" <<std::endl;
          { // compare values
            const auto outputValuesA_Host = Kokkos::create_mirror_view(outputValuesA); Kokkos::deep_copy(outputValuesA_Host, outputValuesA);
            const auto outputValuesB_Host = Kokkos::create_mirror_view(outputValuesB); Kokkos::deep_copy(outputValuesB_Host, outputValuesB);
            
            OutValueType diff = 0; 
            const auto tol = 100.0 * epsilon<double>();  
            for (size_t ic=0;ic<outputValuesA_Host.extent(0);++ic)
              for (size_t i=0;i<outputValuesA_Host.extent(1);++i)
                for (size_t j=0;j<outputValuesA_Host.extent(2);++j) {
                  diff = 0;
                  OutValueType maxMagnitude = 0;
                  for (int d=0;d<ndim;++d) {
                    diff += std::abs(outputValuesB_Host(i,j,d) - outputValuesA_Host(ic,i,j,d));
                    maxMagnitude = std::max(maxMagnitude, std::max(std::abs(outputValuesA_Host(ic,i,j,d)), std::abs(outputValuesB_Host(i,j,d))));
                  }
                  if (diff > tol * std::max(1.0, maxMagnitude)) {
                    ++errorFlag;
                    std::cout << " order: " << order
                              << ", ic: " << ic << ", i: " << i << ", j: " << j 
                              << ", val A: [" << outputValuesA_Host(ic,i,j,0) << ", " << outputValuesA_Host(ic,i,j,1) << "]"
                              << ", val B: [" << outputValuesB_Host(i,j,0) << ", " << outputValuesB_Host(i,j,1) << "]"
                              << ", |diff|: " << diff
                              << ", tol: " << tol
                              << std::endl;
                  }
                }
          }

          { 
            // compare curls
            const auto outputCurlsA_Host = Kokkos::create_mirror_view(outputCurlsA); Kokkos::deep_copy(outputCurlsA_Host, outputCurlsA);
            const auto outputCurlsB_Host = Kokkos::create_mirror_view(outputCurlsB); Kokkos::deep_copy(outputCurlsB_Host, outputCurlsB);
            
            OutValueType diff = 0;
            //Note, the PR intel 2021 serial build shows substantially higher errors (possibly due to operation rearrangements).
            auto tol = 1.0e6*epsilon<double>(); 
            for (size_t ic=0;ic<outputCurlsA_Host.extent(0);++ic)
              for (size_t i=0;i<outputCurlsA_Host.extent(1);++i)
                for (size_t j=0;j<outputCurlsA_Host.extent(2);++j) {
                  diff = 0;
                  OutValueType maxMagnitude = 0;
                  for (int d=0;d<ndim;++d) {
                    diff += std::abs(outputCurlsB_Host(i,j,d) - outputCurlsA_Host(ic,i,j,d));
                    maxMagnitude = std::max(maxMagnitude, std::max(std::abs(outputCurlsA_Host(ic,i,j,d)), std::abs(outputCurlsB_Host(i,j,d))));
                  }
                  if (diff > tol * std::max(1.0, maxMagnitude)) {
                    ++errorFlag;
                    std::cout << " order: " << order
                              << ", ic: " << ic << ", i: " << i << ", j: " << j 
                              << ", curls A: [" << outputCurlsA_Host(ic,i,j,0)<< ", " << outputCurlsA_Host(ic,i,j,1) << ", " << outputCurlsA_Host(ic,i,j,2) << "]"
                              << ", curls B: [" << outputCurlsB_Host(i,j,0) << ", " << outputCurlsB_Host(i,j,1)<< ", " << outputCurlsB_Host(i,j,2) << "]"
                              << ", |diff|: " << diff
                              << ", tol: " << tol
                              << std::endl;
                  }
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
