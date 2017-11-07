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
    catch (std::logic_error err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceSpaceType>
    int CellTools_Test04(const bool verbose) {

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

      typedef CellTools<DeviceSpaceType> ct;
      typedef RealSpaceTools<DeviceSpaceType> rst;
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;

      const ValueType tol = tolerence()*100.0;

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

          DynRankView ConstructWithLabel(cubPointsMax,  Parameters::MaxIntegrationPoints, Parameters::MaxDimension);
          DynRankView ConstructWithLabel(cubWeightsMax, Parameters::MaxIntegrationPoints);
          // for arbitrary high order elements, this should be Parameters::MaxIntegrationPoints 
          const auto maxNumNodes = 64; 
          DynRankView ConstructWithLabel(refCellNodesMax,          maxNumNodes, Parameters::MaxDimension);
          DynRankView ConstructWithLabel(worksetCellMax, numCells, maxNumNodes, Parameters::MaxDimension);

          DynRankView ConstructWithLabel(physPointsMax,    numCells, Parameters::MaxIntegrationPoints, Parameters::MaxDimension);
          DynRankView ConstructWithLabel(controlPointsMax, numCells, Parameters::MaxIntegrationPoints, Parameters::MaxDimension);

          DynRankView ConstructWithLabel(initGuessMax, numCells, Parameters::MaxIntegrationPoints, Parameters::MaxDimension);

          // Loop over cell topologies, make cell workset for each one by perturbing the worksetCell & test methods
          for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];
            
            if (!ct::hasReferenceCell(cell))
              continue;
      
            // 1.   Define a single reference point set using cubature factory with order 4 cubature
            const auto cellCubature = cubFactory.create<DeviceSpaceType>(cell, testAccuracy); 
            const auto cubDim = cellCubature->getDimension();
            const auto cubNumPoints = cellCubature->getNumPoints();
            
            typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
            range_type cubPointRange(0, cubNumPoints), cubDimRange(0, cubDim);            
            auto cubPoints  = Kokkos::subdynrankview( cubPointsMax,  cubPointRange, cubDimRange );
            auto cubWeights = Kokkos::subdynrankview( cubWeightsMax, cubPointRange );

            cellCubature->getCubature(cubPoints, cubWeights);

            // 2.   Define a cell workset by perturbing the worksetCell of the reference cell with the specified topology
            // 2.1  Resize dimensions of the rank-3 (C,N,D) cell workset array for the current topology
            const auto numNodes = cell.getNodeCount();
            const auto cellDim  = cell.getDimension();
            
            range_type nodeRange(0, numNodes), cellDimRange(0, cellDim);            
            auto worksetCell = Kokkos::subdynrankview(worksetCellMax, Kokkos::ALL(), nodeRange, cellDimRange);
            
            // 2.2  Copy worksetCell of the reference cell with the same topology to temp rank-2 (N,D) array
            auto refCellNodes = Kokkos::subdynrankview(refCellNodesMax, nodeRange, cellDimRange);
            ct::getReferenceSubcellNodes(refCellNodes, cellDim, 0, cell);
            
            // 2.3  Create randomly perturbed version of the reference cell and save in the cell workset array
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) {
              // Move vertices +/-0.125 along their axes. Gives nondegenerate cells for base and extended topologies 
              for (size_type nodeOrd=0;nodeOrd<numNodes;++nodeOrd) 
                for(size_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<double>::random()/16.0;
                  worksetCell(cellOrd, nodeOrd, d) = refCellNodes(nodeOrd, d) + delta;
                } 
            }
            
            /* 
             * 3.1 Test 1: single point set to single physical cell: map ref. point set in rank-2 (P,D) array
             *      to a physical point set in rank-2 (P,D) array for a specified cell ordinal. Use the cub.
             *      points array for this test. Resize physPoints and controlPoints to rank-2 (P,D) arrays.
             */
            auto physPoints    = Kokkos::subdynrankview( physPointsMax,    Kokkos::ALL(), cubPointRange, cubDimRange );
            auto controlPoints = Kokkos::subdynrankview( controlPointsMax, Kokkos::ALL(), cubPointRange, cubDimRange );
            
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
              auto initGuess = Kokkos::subdynrankview( initGuessMax, Kokkos::ALL(), cubPointRange, cubDimRange );
              // cubPoints become the initial guess
              rst::clone(initGuess, cubPoints);
              ct::mapToReferenceFrameInitGuess(controlPoints, initGuess, physPoints, worksetCell, cell);
              break;
            }
            }
            
            // Points in controlPoints should match the originals in cubPoints up to a tolerance
            for (auto cellOrd=0;cellOrd<numCells;++cellOrd) 
              for (auto pt=0;pt<cubNumPoints;++pt) 
                for (size_type d=0;d<cellDim;++d) 
                  if( std::abs( controlPoints(cellOrd, pt, d) - cubPoints(pt, d) ) > tol ) {
                    errorFlag++;
                    *outStream
                      << std::setw(70) << "^^^^----FAILURE!" << "\n"
                      << " Mapping a single point set to a single physical cell in a workset failed for: \n"
                      << "                    Cell Topology = " << cell.getName() << "\n"
                      << " Physical cell ordinal in workset = " << cellOrd << "\n"
                      << "          Reference point ordinal = " << std::setprecision(12) << pt << "\n"
                      << "    At reference point coordinate = " << std::setprecision(12) << d << "\n"
                      << "                   Original value = " << cubPoints(pt, d) << "\n"
                      << "                     F^{-1}F(P_d) = " << controlPoints(cellOrd, pt, d) <<"\n";
                  }
          }
        } catch (std::logic_error err) {
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
