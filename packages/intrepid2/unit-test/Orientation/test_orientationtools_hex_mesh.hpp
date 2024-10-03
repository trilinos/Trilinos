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
    int OrientationToolsHexMesh(const bool verbose) {

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
        << "|                 Unit Test (OrientationTools, Hex, getOrientation)           |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      //const double tol = tolerence();

      typedef OrientationTools<DeviceType> ots;
      try {

        {
          *outStream << "\n -- Testing Hexahedron \n\n";

          ordinal_type refMesh[27][8] = { {   8 ,   33 ,   45 ,   10 ,   11 ,   37 ,   57 ,   15  },
                                          {  33 ,   34 ,   46 ,   45 ,   37 ,   38 ,   58 ,   57  },
                                          {  34 ,   21 ,   25 ,   46 ,   38 ,   23 ,   32 ,   58  },
                                          {  10 ,   45 ,   43 ,    9 ,   15 ,   57 ,   59 ,   16  },
                                          {  45 ,   46 ,   44 ,   43 ,   57 ,   58 ,   60 ,   59  },
                                          {  46 ,   25 ,   26 ,   44 ,   58 ,   32 ,   31 ,   60  },
                                          {   9 ,   43 ,   41 ,    5 ,   16 ,   59 ,   51 ,    7  },
                                          {  43 ,   44 ,   42 ,   41 ,   59 ,   60 ,   52 ,   51  },
                                          {  44 ,   26 ,   24 ,   42 ,   60 ,   31 ,   27 ,   52  },
                                          {  11 ,   37 ,   57 ,   15 ,   12 ,   39 ,   61 ,   13  },
                                          {  37 ,   38 ,   58 ,   57 ,   39 ,   40 ,   62 ,   61  },
                                          {  38 ,   23 ,   32 ,   58 ,   40 ,   22 ,   30 ,   62  },
                                          {  15 ,   57 ,   59 ,   16 ,   13 ,   61 ,   63 ,   14  },
                                          {  57 ,   58 ,   60 ,   59 ,   61 ,   62 ,   64 ,   63  },
                                          {  58 ,   32 ,   31 ,   60 ,   62 ,   30 ,   29 ,   64  },
                                          {  16 ,   59 ,   51 ,    7 ,   14 ,   63 ,   49 ,    6  },
                                          {  59 ,   60 ,   52 ,   51 ,   63 ,   64 ,   50 ,   49  },
                                          {  60 ,   31 ,   27 ,   52 ,   64 ,   29 ,   28 ,   50  },
                                          {  12 ,   39 ,   61 ,   13 ,    1 ,   35 ,   53 ,    3  },
                                          {  39 ,   40 ,   62 ,   61 ,   35 ,   36 ,   54 ,   53  },
                                          {  40 ,   22 ,   30 ,   62 ,   36 ,   18 ,   20 ,   54  },
                                          {  13 ,   61 ,   63 ,   14 ,    3 ,   53 ,   55 ,    4  },
                                          {  61 ,   62 ,   64 ,   63 ,   53 ,   54 ,   56 ,   55  },
                                          {  62 ,   30 ,   29 ,   64 ,   54 ,   20 ,   19 ,   56  },
                                          {  14 ,   63 ,   49 ,    6 ,    4 ,   55 ,   47 ,    2  },
                                          {  63 ,   64 ,   50 ,   49 ,   55 ,   56 ,   48 ,   47  },
                                          {  64 ,   29 ,   28 ,   50 ,   56 ,   19 ,   17 ,   48  } };

          //const ordinal_type refOrts[27][6];

          const ordinal_type numCells = 27, numVerts = 8, numEdges = 12, numFaces = 6;
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );

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
            ordinal_type ortEdges[12], ortFaces[6];
            elemOrtsHost(cell).getEdgeOrientation(ortEdges, numEdges);
            elemOrtsHost(cell).getFaceOrientation(ortFaces, numFaces);

            //int flag = 0;
            std::stringstream s1; //, s2;
            s1 << "edge ";
            for (auto edgeId=0;edgeId<numEdges;++edgeId) {
              s1 << ortEdges[edgeId] << "  ";
              //s2 << refOrts[cell][edgeId] << "  ";
              //flag += (orts[edgeId] != refOrts[cell][edgeId]);
            }
            s1 << "face ";
            for (auto faceId=0;faceId<numFaces;++faceId) {
              s1 << ortFaces[faceId] << "  ";
              //s2 << refOrts[cell][edgeId] << "  ";
              //flag += (orts[edgeId] != refOrts[cell][edgeId]);
            }

            *outStream << " cell id = " << std::setw(3) << cell
                       << " orts = " << s1.str()
              //<< " :: ref = " << s2.str()
                       << " \n";
            // if (flag) {
            //   errorFlag += flag;
            //   *outStream << "                    ^^^^^^^^^^^^^^ FAILURE\n ";
            // }

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

