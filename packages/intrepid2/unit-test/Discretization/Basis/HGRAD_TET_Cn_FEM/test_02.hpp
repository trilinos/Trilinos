// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.hpp
    \brief  Unit tests for the Intrepid2::HGRAD_TET_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
 */


#include "Intrepid2_config.h"
#include "Kokkos_Random.hpp"
#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

namespace Intrepid2 {

  namespace Test {

    // This code provides an example to use serial interface of high order elements
    template<typename OutValueType, typename PointValueType, typename DeviceType>
    int HGRAD_TET_Cn_FEM_Test02(const bool verbose) {
      using DeviceSpaceType = typename DeviceType::execution_space;
      Kokkos::print_configuration(std::cout, false);

      int errorFlag = 0;

      try { 
        for (int order=1;order<Parameters::MaxOrder;++order) {
          using TetBasisType = Basis_HGRAD_TET_Cn_FEM<DeviceType,OutValueType,PointValueType>;
          auto tetBasisPtr = Teuchos::rcp(new TetBasisType(order));
          
          // problem setup 
          //   let's say we want to evaluate 1000 points in parallel. output values are stored in outputValuesA and B.
          //   A is compuated via serial interface and B is computed with top-level interface.
          const int ncells = 20, npts = 50, ndim = 3;
          Kokkos::DynRankView<OutValueType,DeviceType> outputValuesA("outputValuesA", ncells, tetBasisPtr->getCardinality(), npts);
          Kokkos::DynRankView<OutValueType,DeviceType> outputValuesB("outputValuesB", tetBasisPtr->getCardinality(), npts);
          
          Kokkos::View<PointValueType**,DeviceType> inputPointsViewToUseRandom("inputPoints", npts, ndim);
          Kokkos::DynRankView<PointValueType,DeviceType> inputPoints (inputPointsViewToUseRandom.data(),  npts, ndim);
          
          // random values between (-1,1) x (-1,1)
          Kokkos::Random_XorShift64_Pool<DeviceType> random(13718);
          Kokkos::fill_random(inputPointsViewToUseRandom, random, 1.0);
          

          { // evaluation using parallel loop over cell
            auto tetBasisPtr_device = copy_virtual_class_to_device<DeviceType,TetBasisType>(*tetBasisPtr);
            auto tetBasisRawPtr_device = tetBasisPtr_device.get();

            int scratch_space_level =1;
            auto functor = KOKKOS_LAMBDA (typename Kokkos::TeamPolicy<DeviceSpaceType>::member_type team_member) {
                auto valsACell = Kokkos::subview(outputValuesA, team_member.league_rank(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
                tetBasisRawPtr_device->getValues(valsACell, inputPoints, OPERATOR_VALUE, team_member, team_member.team_scratch(scratch_space_level));
              };

            const int vectorSize = getVectorSizeForHierarchicalParallelism<PointValueType>();
            Kokkos::TeamPolicy<DeviceSpaceType> teamPolicy(ncells, Kokkos::AUTO,vectorSize);
            //Get the required size of the scratch space per team and per thread.
            int perThreadSpaceSize(0), perTeamSpaceSize(0);
            tetBasisPtr->getScratchSpaceSize(perTeamSpaceSize,perThreadSpaceSize,inputPoints);
            teamPolicy.set_scratch_size(scratch_space_level, Kokkos::PerTeam(perTeamSpaceSize), Kokkos::PerThread(perThreadSpaceSize));

            Kokkos::parallel_for (teamPolicy,functor);
          }


          // evaluation using high level interface
          tetBasisPtr->getValues(outputValuesB, inputPoints, OPERATOR_VALUE);
          
          // compare 
          const auto outputValuesA_Host = Kokkos::create_mirror_view(outputValuesA); Kokkos::deep_copy(outputValuesA_Host, outputValuesA);
          const auto outputValuesB_Host = Kokkos::create_mirror_view(outputValuesB); Kokkos::deep_copy(outputValuesB_Host, outputValuesB);
          
          double sum = 0, diff = 0;
          for (size_t ic=0;ic<outputValuesA_Host.extent(0);++ic)
            for (size_t i=0;i<outputValuesA_Host.extent(1);++i)
              for (size_t j=0;j<outputValuesA_Host.extent(2);++j) {
                sum += std::abs(outputValuesB_Host(i,j));
                diff += std::abs(outputValuesB_Host(i,j) - outputValuesA_Host(ic,i,j));
                if (verbose) {
                  std::cout << " order = " << order
                            << " i = " << i << " j = " << j 
                            << " val A = " << outputValuesA_Host(ic,i,j) 
                            << " val B = " << outputValuesB_Host(i,j) 
                            << " diff  = " << (outputValuesA_Host(ic,i,j) - outputValuesB_Host(i,j)) 
                            << std::endl;
                }
              }
          if (diff/sum > 1.0e-9) {
            errorFlag = -1;
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
