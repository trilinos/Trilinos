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
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_DefaultCubatureFactory.hpp"

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
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceType>
    int CellTools_Test04(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream                                                       
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) Mapping to and from reference cells with base and extended topologies|\n"
        << "|        using default initial guesses when computing the inverse F^{-1}      |\n"
        << "|     2) Repeat all tests from 1) using user-defined initial guess for F^{-1} |\n"
        << "|     3) Exception testing                                                    |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n"
        << "|                      Kara Peterson(kjpeter@sandia.gov), or                  |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov)                       |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      using ct = CellTools<DeviceType>;
      using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

      const ValueType tol = tolerence<ValueType>()*100.0;

      int errorFlag  = 0;

      // Collect all supported cell topologies
      std::vector<shards::CellTopology> standardBaseTopologies;
      shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);
      const auto topoSize = standardBaseTopologies.size();

      for (auto testCase=0;testCase<2;++testCase) {
        try {
          switch (testCase) {
          case 0: {
            *outStream
              << "\n"
              << "===============================================================================\n" 
              << "| Test 1: computing F(x) and F^{-1}(x) using default initial guesses.         |\n" 
              << "===============================================================================\n\n";
            break;
          }
          case 1: {
            *outStream
              << "\n"
              << "===============================================================================\n" 
              << "| Test 2: computing F(x) and F^{-1}(x) using user-defined initial guess.      |\n" 
              << "===============================================================================\n\n";
            break;
          }
          } 

          /*
           *  Test summary:
           *
           *    A reference point set is mapped to physical frame and then back to reference frame.
           *    Test passes if the final set of points matches the first set of points. The cell workset
           *    is generated by perturbing randomly the worksetCell of a reference cell with the specified 
           *    cell topology. 
           *
           */
          DefaultCubatureFactory cubFactory;   
          
          // Initialize testing env.
          const auto numCells = 10;
          const auto testAccuracy = 4;

          // Loop over cell topologies, make cell workset for each one by perturbing the worksetCell & test methods
          for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];
            
            if (!ct::hasReferenceCell(cell))
              continue;
      
            // 1.   Define a single reference point set using cubature factory with order 4 cubature
            const auto cellCubature = cubFactory.create<DeviceType,ValueType,ValueType>(cell, testAccuracy);
            const auto cubDim = cellCubature->getDimension();
            const auto cubNumPoints = cellCubature->getNumPoints();
            
            DynRankView ConstructWithLabel(cubPoints,  cubNumPoints, cubDim );
            DynRankView ConstructWithLabel(cubWeights, cubNumPoints );
            cellCubature->getCubature(cubPoints, cubWeights);

            // 2.   Define a cell workset by perturbing the worksetCell of the reference cell with the specified topology
            // 2.1  Resize dimensions of the rank-3 (C,N,D) cell workset array for the current topology
            const auto numNodes = cell.getNodeCount();
            const auto cellDim  = cell.getDimension();

            DynRankView ConstructWithLabel(worksetCell, numCells, numNodes, cellDim);
            
            // 2.2  Copy worksetCell of the reference cell with the same topology to temp rank-2 (N,D) array
            DynRankView ConstructWithLabel(refCellNodes, numNodes, cellDim);
            ct::getReferenceSubcellNodes(refCellNodes, cellDim, 0, cell);
            
            // create mirror host views
            auto hRefCellNodes = Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace(), refCellNodes);

            auto hWorksetCell = Kokkos::create_mirror_view(worksetCell);

            // 2.3  Create randomly perturbed version of the reference cell and save in the cell workset array
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) {
              // Move vertices +/-0.125 along their axes. Gives nondegenerate cells for base and extended topologies 
              for (size_type nodeOrd=0;nodeOrd<numNodes;++nodeOrd) 
                for(size_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<double>::random()/16.0;
                  hWorksetCell(cellOrd, nodeOrd, d) = hRefCellNodes(nodeOrd, d) + delta;
                } 
            }
            Kokkos::deep_copy(worksetCell,hWorksetCell);
            
            /* 
             * 3.1 Test 1: single point set to single physical cell: map ref. point set in rank-2 (P,D) array
             *      to a physical point set in rank-2 (P,D) array for a specified cell ordinal. Use the cub.
             *      points array for this test. Resize physPoints and controlPoints to rank-2 (P,D) arrays.
             */
            DynRankView ConstructWithLabel(physPoints,    numCells, cubNumPoints, cubDim );
            DynRankView ConstructWithLabel(controlPoints, numCells, cubNumPoints, cubDim );
            
            *outStream 
              << " Mapping a set of " << cubNumPoints << " points to one cell in a workset of " << numCells << " " 
              << cell.getName() << " cells. \n";

            // Forward map:: requires cell ordinal; input (P,D) -> output (C,P,D)
            ct::mapToPhysicalFrame(physPoints, cubPoints, worksetCell, cell);

            // Inverse map: requires cell ordinal; input (C,P,D) -> output (C,P,D)
            switch (testCase) {
            case 0: {
              ct::mapToReferenceFrame(controlPoints, physPoints, worksetCell, cell);
              break;
            }
            case 1: {
              DynRankView ConstructWithLabel(initGuess, numCells, cubNumPoints, cubDim );
              // cubPoints become the initial guess
              RealSpaceTools<DeviceType>::clone(initGuess, cubPoints);
              ct::mapToReferenceFrameInitGuess(controlPoints, initGuess, physPoints, worksetCell, cell);
              break;
            }
            }
            
            // copy to mirror host views
            auto hControlPoints = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), controlPoints);
            auto hCubPoints = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cubPoints);

            // Points in controlPoints should match the originals in cubPoints up to a tolerance
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) 
              for (auto pt=0;pt<cubNumPoints;++pt) 
                for (size_type d=0;d<cellDim;++d) 
                  if( std::abs( hControlPoints(cellOrd, pt, d) - hCubPoints(pt, d) ) > tol ) {
                    errorFlag++;
                    *outStream
                      << std::setw(70) << "^^^^----FAILURE!" << "\n"
                      << " Mapping a single point set to a single physical cell in a workset failed for: \n"
                      << "                    Cell Topology = " << cell.getName() << "\n"
                      << " Physical cell ordinal in workset = " << cellOrd << "\n"
                      << "          Reference point ordinal = " << std::setprecision(12) << pt << "\n"
                      << "    At reference point coordinate = " << std::setprecision(12) << d << "\n"
                      << "                   Original value = " << hCubPoints(pt, d) << "\n"
                      << "                     F^{-1}F(P_d) = " << hControlPoints(cellOrd, pt, d) <<"\n";
                  }
          }
        } catch (std::logic_error &err) {
          *outStream << err.what() << "\n";
          errorFlag = -1000;
        }
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);

      return errorFlag;
    }
  } // end test
} // end intrepid2
