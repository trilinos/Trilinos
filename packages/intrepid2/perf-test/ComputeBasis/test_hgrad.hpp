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

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "test_util.hpp"
#include "test_functors.hpp"

namespace Intrepid2 {
  
  namespace Test {

    template<typename ValueType, typename DeviceSpaceType>
    int ComputeBasis_HGRAD(const ordinal_type nworkset,
                           const ordinal_type C,
                           const ordinal_type order,
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

      using BasisTypeHost = Basis_HGRAD_HEX_C1_FEM<HostSpaceType,ValueType,ValueType>;
      using ImplBasisType = Impl::Basis_HGRAD_HEX_C1_FEM;
      using range_type = Kokkos::pair<ordinal_type,ordinal_type>;

      typedef typename DeviceSpaceType::array_layout DeviceArrayLayout;
      typedef typename HostSpaceType::array_layout HostArrayLayout;

      constexpr size_t LLC_CAPACITY = 32*1024*1024;
      Intrepid2::Test::Flush<LLC_CAPACITY,DeviceSpaceType> flush;
      
      Kokkos::Timer timer;
      double t_horizontal = 0, t_vertical = 0;
      int errorFlag = 0;

      BasisTypeHost hostBasis;
      const auto cellTopo = hostBasis.getBaseCellTopology();

      auto cubature = DefaultCubatureFactory::create<DeviceSpaceType,ValueType,ValueType>(cellTopo, order);

      const ordinal_type 
        numCells = C, 
        numVerts = cellTopo.getVertexCount(),
        numDofs = hostBasis.getCardinality(),
        numPoints = cubature->getNumPoints(), 
        spaceDim = cubature->getDimension();

      Kokkos::DynRankView<ValueType,DeviceArrayLayout,HostSpaceType> 
	dofCoordsHost("dofCoordsHost", numDofs, spaceDim);
      hostBasis.getDofCoords(dofCoordsHost);
      const auto refNodesHost = Kokkos::subview(dofCoordsHost, range_type(0, numVerts), Kokkos::ALL());
      
      // pertub nodes
      Kokkos::DynRankView<ValueType,DeviceArrayLayout,HostSpaceType> 
	worksetCellsHost("worksetCellsHost", numCells, numVerts, spaceDim);
      for (ordinal_type cell=0;cell<numCells;++cell) {
        for (ordinal_type i=0;i<numVerts;++i)
          for (ordinal_type j=0;j<spaceDim;++j) {
            ValueType val = (rand()/(RAND_MAX + 1.0))*0.2 -0.1;
            worksetCellsHost(cell, i, j) = refNodesHost(i, j) + val; 
          }
      }

      auto worksetCells = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), worksetCellsHost);
      Kokkos::deep_copy(worksetCells, worksetCellsHost);
      
      Kokkos::DynRankView<ValueType,DeviceSpaceType> 
	refPoints("refPoints", numPoints, spaceDim), 
	refWeights("refWeights", numPoints);
      cubature->getCubature(refPoints, refWeights);

      std::cout
        << "===============================================================================\n" 
        << " Performance Test evaluating ComputeBasis \n"
        << " # of workset = " << nworkset << "\n" 
        << " Test Array Structure (C,F,P,D) = " << numCells << ", " << numDofs << ", " << numPoints << ", " << spaceDim << "\n"
        << "===============================================================================\n";

      *verboseStream
        << "\n"
        << "===============================================================================\n"
        << "TEST 1: evaluateFields flat version\n"
        << "===============================================================================\n";
      
      try {
        Kokkos::DynRankView<ValueType,DeviceSpaceType> 
          refBasisValues("refBasisValues", numDofs, numPoints),
          refBasisGrads ("refBasisGrads",  numDofs, numPoints, spaceDim);
        
        auto space = typename DeviceSpaceType::execution_space();
        
        ImplBasisType::getValues<DeviceSpaceType>(space, refBasisValues, refPoints, OPERATOR_VALUE);
        ImplBasisType::getValues<DeviceSpaceType>(space, refBasisGrads,  refPoints, OPERATOR_GRAD);

	std::cout << " Ref completed\n";
        
        const ordinal_type ibegin = -3;
        // testing sequential appraoch
        {          
          Kokkos::DynRankView<ValueType,DeviceSpaceType> 
            jacobian           ("jacobian",            numCells,          numPoints, spaceDim, spaceDim),
            jacobianInv        ("jacobianInv",         numCells,          numPoints, spaceDim, spaceDim),
            jacobianDet        ("jacobianDet",         numCells,          numPoints),
            cellMeasure        ("cellMeasure",         numCells,          numPoints),
            phyBasisValues     ("phyBasisValues",      numCells, numDofs, numPoints),
            phyBasisGrads      ("phyBasisGrads",       numCells, numDofs, numPoints, spaceDim),
            weightedBasisValues("weightedBasisValues", numCells, numDofs, numPoints),
            weightedBasisGrads ("weightedBasisGrads",  numCells, numDofs, numPoints, spaceDim);
          
          typedef CellTools<DeviceSpaceType> cts;
          typedef FunctionSpaceTools<DeviceSpaceType> fts;

          for (ordinal_type iwork=ibegin;iwork<nworkset;++iwork) {
            
            flush.run();

            DeviceSpaceType().fence();
            timer.reset();

            cts::setJacobian(jacobian, refPoints, worksetCells, cellTopo);
            cts::setJacobianInv(jacobianInv, jacobian);
            cts::setJacobianDet(jacobianDet, jacobian);
            fts::computeCellMeasure(cellMeasure, jacobianDet, refWeights);
            
            fts::HGRADtransformVALUE(phyBasisValues, refBasisValues);
            fts::multiplyMeasure(weightedBasisValues, cellMeasure, phyBasisValues);
            
            fts::HGRADtransformGRAD(phyBasisGrads, jacobianInv, refBasisGrads);
            fts::multiplyMeasure(weightedBasisGrads, cellMeasure, phyBasisGrads);

            DeviceSpaceType().fence();
            t_horizontal += (iwork >= 0)*timer.seconds();
          }
        }

        // testing vertical approach
        {          
          Kokkos::DynRankView<ValueType,DeviceSpaceType> 
            weightedBasisValues("weightedBasisValues", numCells, numDofs, numPoints),
            weightedBasisGrads ("weightedBasisGrads",  numCells, numDofs, numPoints, spaceDim);

          typedef F_hgrad_eval<ValueType,ValueType,DeviceSpaceType> FunctorType;

          using range_policy_type = Kokkos::MDRangePolicy
            < DeviceSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
          range_policy_type policy( {        0,         0 },
                                    { numCells, numPoints } );

          FunctorType functor(weightedBasisValues,
                              weightedBasisGrads,
                              refBasisGrads,
                              worksetCells,
                              refWeights,
                              refBasisValues,
                              refBasisGrads);
          
          for (ordinal_type iwork=ibegin;iwork<nworkset;++iwork) {
            flush.run();

            DeviceSpaceType().fence();
            timer.reset();
            
            Kokkos::parallel_for(policy, functor);

            DeviceSpaceType().fence();
            t_vertical += (iwork >= 0)*timer.seconds();
          }

        }

      } catch (std::exception err) {
        *verboseStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *verboseStream << err.what() << '\n';
        *verboseStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }
      
      std::cout 
        << "TEST HGRAD " 
        << ": t_horizontal = " << (t_horizontal/nworkset) 
        << ", t_vertical = " << (t_vertical/nworkset)
        << std::endl;
      
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
