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
    int OrientationToolsModifyBasis_QUAD_HGRAD(const bool verbose) {

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
        << "|           Unit Test (OrientationTools, Quad, getModifiedHgradBasis)          |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      const double tol = tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

      typedef OrientationTools<DeviceType> ots;
      try {
        constexpr ordinal_type order = 4;
        if(order <= maxOrder)
        {
          *outStream << "\n -- Testing Quadrilateral \n\n";



          Basis_HGRAD_QUAD_Cn_FEM<DeviceType> cellBasis(order);
          const auto cellTopo = cellBasis.getBaseCellTopology();
          const ordinal_type ndofBasis = cellBasis.getCardinality();

          
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
          const ordinal_type numCells = 9, numVerts = 4, numEdges = 4;

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

            const ordinal_type reverse[numEdges][2] = { { 0,1 },
                                                        { 0,1 },
                                                        { 1,0 },
                                                        { 1,0 } };
            ordinal_type orts[numEdges];
            elemOrtsHost(cell).getEdgeOrientation(orts, numEdges);

            s1 << " :: edge(0000) = " ;
            s2 << " :: edge(" << orts[0] << orts[1] << orts[2] << orts[3] << ") = ";
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

            s1 << " :: intr = " ;
            s2 << " :: intr = " ;
            {
              const auto ndof = cellBasis.getDofTag(cellBasis.getDofOrdinal(2, 0, 0))(3);
              for (auto i=0;i<ndof;++i) {
                const auto ord = cellBasis.getDofOrdinal(2, 0, i);
                s1 << std::setw(4) << refValuesHost(cell, ord);              
                s2 << std::setw(4) << outValuesHost(cell, ord);              
                flag += (std::abs(outValuesHost(cell, ord) - refValuesHost(cell, ord)) > tol);
              }
            }

            *outStream << "\n cell = " << cell << "\n"
                       << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n";
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
    
