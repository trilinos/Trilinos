// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename DeviceType>
    int OrientationToolsQuadMesh(const bool verbose) {

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
        << "|                 Unit Test (OrientationTools, Quad, getOrientation)          |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      //const double tol = tolerence();

      typedef OrientationTools<DeviceType> ots;
      try {

        {
          *outStream << "\n -- Testing Quadrilateral \n\n";

          // 
          // 9 12 13 16
          // 4  3 11 15
          // 5  2  8 14
          // 1  6  7 10
          ordinal_type refMesh[9][4] = { { 1, 6, 2, 5 },
                                         { 6, 7, 8, 2 },
                                         { 7,10,14, 8 },
                                         { 5, 2, 3, 4 },
                                         { 2, 8,11, 3 },
                                         { 8,14,15,11 },
                                         { 4, 3,12, 9 },
                                         { 3,11,13,12 },
                                         {11,15,16,13 } };
          
          const ordinal_type refOrts[9][4] = { { 0,1,0,1 },
                                               { 0,0,1,0 },
                                               { 0,0,1,1 },
                                               { 1,0,0,0 },
                                               { 0,0,1,1 },
                                               { 0,0,1,1 },
                                               { 1,0,1,1 },
                                               { 0,0,1,1 },
                                               { 0,0,1,1 } };
          
          const ordinal_type numCells = 9, numVerts = 4, numEdges = 4;
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

          // view to store orientation 
          Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);

          // view to import refMesh from host
          Kokkos::DynRankView<ordinal_type,Kokkos::LayoutRight,Kokkos::HostSpace> 
            elemNodesHost(&refMesh[0][0], numCells, numVerts);
          auto elemNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), elemNodesHost);
          Kokkos::deep_copy(elemNodes, elemNodesHost);
          
          ots::getOrientation(elemOrts, elemNodes, cellTopo);

          // for comparison move data to host
          auto elemOrtsHost = Kokkos::create_mirror_view( elemOrts);
          
          for (auto cell=0;cell<numCells;++cell) {
            // decode edge orientations
            ordinal_type orts[4];
            elemOrtsHost(cell).getEdgeOrientation(orts, numEdges);

            int flag = 0;
            std::stringstream s1, s2;
            for (auto edgeId=0;edgeId<numEdges;++edgeId) {
              s1 << orts[edgeId] << "  ";
              s2 << refOrts[cell][edgeId] << "  ";
              flag += (orts[edgeId] != refOrts[cell][edgeId]);              
            }
            *outStream << " cell id = " << cell 
                       << "  computed edge ort = " << s1.str() 
                       << " :: ref edge ort = " << s2.str() 
                       << " \n";
            if (flag) {
              errorFlag += flag;
              *outStream << "                    ^^^^^^^^^^^^^^ FAILURE\n ";
            }
          }
        }

      } catch (std::exception &err) {
        std::cout << " Exeption\n";
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED = " << errorFlag << "\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  }
}
    
