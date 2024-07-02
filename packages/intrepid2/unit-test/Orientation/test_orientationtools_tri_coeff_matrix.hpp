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

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
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

    // test
    template<typename DeviceType>
    int OrientationToolsTriCoeffMatrix(const bool verbose) {

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
        << "|                 Unit Test (OrientationTools, initialize)                    |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      const double tol = tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;

      typedef OrientationTools<DeviceType> ots;
      try {
        const ordinal_type testOrderBegin = 1, testOrderEnd = std::min(4, maxOrder);
        for (ordinal_type testOrder=testOrderBegin;testOrder<=testOrderEnd;++testOrder) {
          *outStream << "\n -- Testing order " << testOrder << "\n"
                     << "===============================================================================\n"
                     << "\n";
          {
            // reference permutation is either inorder or inverse order
            const ordinal_type refCoeff[3][2] = { { 0, 1 },
                                                  { 0, 1 },
                                                  { 1, 0 } };

            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
            const ordinal_type order = testOrder;

            {
              *outStream << "\n -- Testing Triangle HGRAD \n\n";
              Basis_HGRAD_TRI_Cn_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                const ordinal_type lineDim = 1, numEdge = 3, numOrt = 2;
                for (auto edgeId=0;edgeId<numEdge;++edgeId) {
                  if (cellBasis.getDofCount(lineDim, edgeId)) {
                    const auto ordEdge = cellBasis.getDofOrdinal(lineDim, edgeId, 0);
                    const auto ndofEdge = cellBasis.getDofTag(ordEdge)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofEdge);
                    for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
                      *outStream << "\n edgeId = " << edgeId << " edgeOrt = " << edgeOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, edgeId, edgeOrt, range, range);
                      for (auto i=0;i<ndofEdge;++i) {
                        for (auto j=0;j<ndofEdge;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";

                        switch (refCoeff[edgeId][edgeOrt]) {
                        case 0: errorFlag += (std::abs(mat(i,i) - 1.0) > tol); break;
                        case 1: errorFlag += (std::abs(mat(i,ndofEdge-i-1) - 1.0) > tol); break;
                        }
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  edge " << edgeId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            {
              *outStream << "\n -- Testing Triangle HCURL \n\n";
              Basis_HCURL_TRI_In_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                if (testOrder == 1) {
                  const ordinal_type ndofBasis = cellBasis.getCardinality();
                  Kokkos::DynRankView<double,Kokkos::HostSpace>
                    outputValues("output", ndofBasis, 3, 2),
                    inputPoints("input", 3, 2);
                  
                  inputPoints(0, 0) =  0.25; inputPoints(0, 1) = 0;
                  inputPoints(1, 0) =  0.50; inputPoints(1, 1) = 0;
                  inputPoints(2, 0) =  0.75; inputPoints(2, 1) = 0;
                  
                  cellBasis.getValues(outputValues, inputPoints);

                  *outStream << "1st edge 3 point eval (0.25, 0.50, 0.75)\n";
                  for (int i=0;i<ndofBasis;++i) {                  
                    *outStream << std::setw(3) << i;
                    for (int j=0;j<3;++j) 
                      *outStream << " :: " << outputValues(i,j,0) << " " << outputValues(i,j,1) << " , ";
                    *outStream << "\n";
                  }
                }
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                const ordinal_type lineDim = 1, numEdge = 3, numOrt = 2;
                for (auto edgeId=0;edgeId<numEdge;++edgeId) {
                  if (cellBasis.getDofCount(lineDim, edgeId)) {
                    const auto ordEdge = cellBasis.getDofOrdinal(lineDim, edgeId, 0);
                    const auto ndofEdge = cellBasis.getDofTag(ordEdge)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofEdge);
                    for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
                      *outStream << "\n edgeId = " << edgeId << " edgeOrt = " << edgeOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, edgeId, edgeOrt, range, range);
                      for (auto i=0;i<ndofEdge;++i) {
                        for (auto j=0;j<ndofEdge;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";

                        // this is not anymore simple permutation;
                        // to verify orientation matrix, it needs to compare with ooutput file
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  edge " << edgeId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            {
              *outStream << "\n -- Testing Triangle HDIV \n\n";
              Basis_HDIV_TRI_In_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                const ordinal_type lineDim = 1, numEdge = 3, numOrt = 2;
                for (auto edgeId=0;edgeId<numEdge;++edgeId) {
                  if (cellBasis.getDofCount(lineDim, edgeId)) {
                    const auto ordEdge = cellBasis.getDofOrdinal(lineDim, edgeId, 0);
                    const auto ndofEdge = cellBasis.getDofTag(ordEdge)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofEdge);
                    for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
                      *outStream << "\n edgeId = " << edgeId << " edgeOrt = " << edgeOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, edgeId, edgeOrt, range, range);
                      for (auto i=0;i<ndofEdge;++i) {
                        for (auto j=0;j<ndofEdge;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";

                        // this is not anymore simple permutation;
                        // to verify orientation matrix, it needs to compare with ooutput file
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  edge " << edgeId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            ots::clearCoeffMatrix();
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

