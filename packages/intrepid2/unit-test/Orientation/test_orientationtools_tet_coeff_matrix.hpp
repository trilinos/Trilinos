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

    template<typename DeviceType>
    int OrientationToolsTetCoeffMatrix(const bool verbose) {

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
      //const double tol = tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;

      typedef OrientationTools<DeviceType> ots;
      try {
        const ordinal_type testOrderBegin = 1, testOrderEnd = std::min(4, maxOrder);
        for (ordinal_type testOrder=testOrderBegin;testOrder<=testOrderEnd;++testOrder) {
          *outStream << "\n -- Testing order " << testOrder << "\n"
                     << "===============================================================================\n"
                     << "\n";
          {
            // I do not know how I can provide reference permutation for tet or any 3D elements
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
            const ordinal_type order = testOrder;

            {
              *outStream << "\n -- Testing Tetrahedral HGRAD \n\n";
              Basis_HGRAD_TET_Cn_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                // check face only
                const ordinal_type faceDim = 2, numEdge = 6, numFace = 4, numOrt = 6;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, numEdge + faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            ots::clearCoeffMatrix();
          }

          {
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
            const ordinal_type order = testOrder;

            {
              *outStream << "\n -- Testing Tetrahedral HCURL \n\n";
              Basis_HCURL_TET_In_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                // check face only
                const ordinal_type faceDim = 2, numEdge = 6, numFace = 4, numOrt = 6;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, numEdge + faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            ots::clearCoeffMatrix();
          }

          {
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
            const ordinal_type order = testOrder;

            {
              *outStream << "\n -- Testing Tetrahedral HDIV \n\n";
              Basis_HDIV_TET_In_FEM<DeviceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view( matData);
                Kokkos::deep_copy(matDataHost, matData);

                // check face only
                const ordinal_type faceDim = 2, numFace = 4, numOrt = 6;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);

                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";
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

