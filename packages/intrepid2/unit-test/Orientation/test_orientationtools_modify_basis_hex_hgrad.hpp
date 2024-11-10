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
    int OrientationToolsModifyBasis_HEX_HGRAD(const bool verbose) {

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
        << "|           Unit Test (OrientationTools, Hex, getModifiedHgradBasis)          |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      const double tol = 100*tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

      typedef OrientationTools<DeviceType> ots;
      try {
        constexpr ordinal_type order = 3;
        if(order <= maxOrder)
        {
          *outStream << "\n -- Testing Hexahedron \n\n";



          Basis_HGRAD_HEX_Cn_FEM<DeviceType> cellBasis(order);
          const auto cellTopo = cellBasis.getBaseCellTopology();
          const ordinal_type ndofBasis = cellBasis.getCardinality();

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
          
          const ordinal_type numCells = 27, numVerts = 8, numEdges = 12, numFaces = 6;

          // view to import refMesh from host          
          Kokkos::DynRankView<ordinal_type,Kokkos::LayoutRight,Kokkos::HostSpace> 
            elemNodesHost(&refMesh[0][0], numCells, numVerts);
          auto elemNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), elemNodesHost);
          Kokkos::deep_copy(elemNodes, elemNodesHost);
          
          // compute orientations for cells (one time computation)
          Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
          ots::getOrientation(elemOrts, elemNodes, cellTopo);

          auto elemOrtsHost = Kokkos::create_mirror_view( elemOrts);
          Kokkos::deep_copy(elemOrtsHost, elemOrts);

          // cell specific modified basis 
          Kokkos::DynRankView<double,DeviceType> outValues("outValues", numCells, ndofBasis);
          Kokkos::DynRankView<double,DeviceType> refValues("refValues", numCells, ndofBasis);

          auto refValuesHost = Kokkos::create_mirror_view( refValues);
          for (auto cell=0;cell<numCells;++cell)           
            for (auto bf=0;bf<ndofBasis;++bf) 
              refValuesHost(cell, bf) = (bf+1);
          Kokkos::deep_copy(refValues, refValuesHost);

          // modify refValues accounting for orientations
          ots::modifyBasisByOrientation(outValues,
                                        refValues,
                                        elemOrts,
                                        &cellBasis);

          auto outValuesHost = Kokkos::create_mirror_view( outValues);
          Kokkos::deep_copy(outValuesHost, outValues);

          for (auto cell=0;cell<numCells;++cell) {
            *outStream << "\n cell = " << cell << "\n";

            int flag = 0 ;
            std::stringstream s1, s2;
            
            s1 << " :: vert = " ;
            s2 << " :: vert = " ;
            for (auto vertId=0;vertId<numVerts;++vertId) {
              const auto ord = cellBasis.getDofOrdinal(0, vertId, 0);
              s1 << std::setw(4) << refValuesHost(cell, ord);              
              s2 << std::setw(4) << outValuesHost(cell, ord); 
              flag += (std::abs(outValuesHost(cell, ord) - refValuesHost(cell, ord)) > tol);
            }
            *outStream << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n\n";
            s1.str("");s2.str("");

            const ordinal_type reverse[numEdges][2] = { { 0,1 }, // bottom
                                                        { 0,1 },
                                                        { 1,0 },
                                                        { 1,0 },
                                                        { 0,1 }, // top
                                                        { 0,1 },
                                                        { 1,0 },
                                                        { 1,0 },
                                                        { 0,1 }, // vertical ones
                                                        { 0,1 },
                                                        { 0,1 },
                                                        { 0,1 } };
            ordinal_type orts[numEdges];
            elemOrtsHost(cell).getEdgeOrientation(orts, numEdges);
            
            s1 << " :: edge(000000000000) = " ;
            s2 << " :: edge(" 
               << orts[0] << orts[1] << orts[2] << orts[3] 
               << orts[4] << orts[5] << orts[6] << orts[7] 
               << orts[8] << orts[9] << orts[10] << orts[11] 
               << ") = ";
            for (auto edgeId=0;edgeId<numEdges;++edgeId) {
              const auto ndof = cellBasis.getDofTag(cellBasis.getDofOrdinal(1, edgeId, 0))(3);
              for (auto i=0;i<ndof;++i) {
                const ordinal_type ii = (reverse[edgeId][orts[edgeId]] ? ndof-i-1 : i);

                const auto refOrd = cellBasis.getDofOrdinal(1, edgeId, ii);
                const auto outOrd = cellBasis.getDofOrdinal(1, edgeId, i);
                s1 << std::setw(4) << refValuesHost(cell, outOrd);              
                s2 << std::setw(4) << outValuesHost(cell, outOrd);              

                flag += (std::abs(outValuesHost(cell, outOrd) - refValuesHost(cell, refOrd)) > tol);
              }

              s1 << " // ";
              s2 << " // ";
            }
            *outStream << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n\n";
            s1.str("");s2.str("");

            s1 << " :: face(000000) = " ;
            s2 << " :: face(" << orts[0] << orts[1] << orts[2] << orts[3] << orts[4] << orts[5] << ") = ";
            for (auto faceId=0;faceId<numFaces;++faceId) {
              const auto ndof = cellBasis.getDofTag(cellBasis.getDofOrdinal(2, faceId, 0))(3);
              for (auto i=0;i<ndof;++i) {
                //const auto refOrd = cellBasis.getDofOrdinal(2, faceId, i);
                const auto outOrd = cellBasis.getDofOrdinal(2, faceId, i);
                s1 << std::setw(4) << refValuesHost(cell, outOrd);              
                s2 << std::setw(4) << outValuesHost(cell, outOrd);              

                // i dont know how to verify this consistency
              }

              s1 << " // ";
              s2 << " // ";
            }
            *outStream << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n\n";
            s1.str("");s2.str("");

            s1 << " :: intr = " ;
            s2 << " :: intr = " ;
            {
              const auto ndof = cellBasis.getDofTag(cellBasis.getDofOrdinal(3, 0, 0))(3);
              for (auto i=0;i<ndof;++i) {
                const auto ord = cellBasis.getDofOrdinal(3, 0, i);
                s1 << std::setw(4) << refValuesHost(cell, ord);              
                s2 << std::setw(4) << outValuesHost(cell, ord);              
                flag += (std::abs(outValuesHost(cell, ord) - refValuesHost(cell, ord)) > tol);
              }
            }
            *outStream << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n\n";
            s1.str("");s2.str("");

            if (flag) {
              *outStream << "                      ^^^^^^^^^^^^ FAILURE\n";
              errorFlag += flag;
            }
          }
          ots::clearCoeffMatrix();
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
    
