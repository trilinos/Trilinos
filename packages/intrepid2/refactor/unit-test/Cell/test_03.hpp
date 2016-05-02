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

//#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureTensor.hpp"

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
    int CellTools_Test03(const bool verbose) {

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
        << "|     1) Edge parametrizations                                                |\n"
        << "|     2) Face parametrizations                                                |\n"
        << "|     3) Edge tangents                                                        |\n"
        << "|     4) Face tangents and normals                                            |\n"
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

      const ValueType tol = Parameters::Tolerence*100.0;

      int errorFlag = 0;

      try {

        /***********************************************************************************************
         *
         * Common for test 3 and 4: edge tangents and face normals for standard cells with base topo
         *
         **********************************************************************************************/

        // Allocate storage and extract all standard cells with base topologies
        std::vector<shards::CellTopology> standardBaseTopologies;
        shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);
        const auto topoSize = standardBaseTopologies.size();

        // Define topologies for the edge and face parametrization domains. (faces are Tri or Quad)
        const auto paramEdge     = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
        const auto paramTriFace  = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
        const auto paramQuadFace = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

        *outStream
          << "\n"
          << "===============================================================================\n"
          << "| Test 4: face/side normals for stand. 3D cells with base topologies:         |\n"
          << "===============================================================================\n\n";

        // This test loops over standard 3D cells with base topologies, creates a set of nodes and tests normals
        {
          // Define cubature on the edge parametrization domain:
          const auto testAccuracy = 6;
          CubatureDirectLineGauss<DeviceSpaceType> edgeCubature(testAccuracy);
          //CubatureDirectTriangleGauss<DeviceSpaceType> triFaceCubature(testAccuracy);
          CubatureTensor<DeviceSpaceType> quadFaceCubature( edgeCubature, edgeCubature );

          const auto faceCubDim = 2;
          //const auto numTriFaceCubPoints  = triFaceCubature.getNumPoints();
          const auto numQuadFaceCubPoints = quadFaceCubature.getNumPoints();

          // Allocate storage for cubature points and weights on face parameter domain and fill with points:
          //DynRankView ConstructWithLabel(paramTriFacePoints,   numTriFaceCubPoints,  faceCubDim);
          //DynRankView ConstructWithLabel(paramTriFaceWeights,  numTriFaceCubPoints);
          DynRankView ConstructWithLabel(paramQuadFacePoints,  numQuadFaceCubPoints, faceCubDim);
          DynRankView ConstructWithLabel(paramQuadFaceWeights, numQuadFaceCubPoints);

          //triFaceCubature.getCubature(paramTriFacePoints, paramTriFaceWeights);
          quadFaceCubature.getCubature(paramQuadFacePoints, paramQuadFaceWeights);

          // Loop over admissible topologies
          for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];

            // skip cells not supported
            if (!ct::hasReferenceCell(cell))
              continue;

            // Exclude 2D and Pyramid<5> cells
            if ( cell.getDimension() == 3 && cell.getKey() != shards::Pyramid<5>::key ) {
              const auto cellDim  = cell.getDimension();
              const auto nCount   = cell.getNodeCount();
              const auto vCount   = cell.getVertexCount();

              DynRankView ConstructWithLabel(refCellVertices, nCount, cellDim);
              ct::getReferenceSubcellVertices(refCellVertices, cellDim, 0, cell);

              *outStream << " Testing face/side normals for cell topology " <<  cell.getName() <<"\n";

              // Array for physical cell vertices ( must have rank 3 for setJacobians)
              DynRankView ConstructWithLabel(physCellVertices, 1, vCount, cellDim);

              // Randomize reference cell vertices by moving them up to +/- (1/8) units along their
              // coordinate axis. Guaranteed to be non-degenerate for standard cells with base topology
              for (size_type v=0;v<vCount;++v)
                for (size_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<double>::random()/8.0;
                  physCellVertices(0, v, d) = refCellVertices(v, d) + delta;
                }

              // Allocate storage for cub. points on a ref. face; Jacobians, phys. face normals and
              // benchmark normals.
              //DynRankView ConstructWithLabel(refTriFacePoints, numTriFaceCubPoints, cellDim);
              //DynRankView ConstructWithLabel(triFacePointsJacobians, 1, numTriFaceCubPoints, cellDim, cellDim);
              //DynRankView ConstructWithLabel(triFacePointNormals, 1, numTriFaceCubPoints, cellDim);
              //DynRankView ConstructWithLabel(triSidePointNormals, 1, numTriFaceCubPoints, cellDim);

              DynRankView ConstructWithLabel(refQuadFacePoints,          numQuadFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(quadFacePointsJacobians, 1, numQuadFaceCubPoints, cellDim, cellDim);
              DynRankView ConstructWithLabel(quadFacePointNormals,    1, numQuadFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(quadSidePointNormals,    1, numQuadFaceCubPoints, cellDim);

              // Loop over faces:
              for (size_type faceOrd=0;faceOrd<cell.getSideCount();++faceOrd) {

                // This test presently includes only Triangle<3> and Quadrilateral<4> faces. Once we support
                // cells with extended topologies we will add their faces as well.
                switch ( cell.getCellTopologyData(2, faceOrd)->key ) {
                  // case shards::Triangle<3>::key: {
                  //   // Compute face normals using CellTools
                  //   ct::mapToReferenceSubcell(refTriFacePoints, paramTriFacePoints, 2, faceOrd, cell);
                  //   ct::setJacobian(triFacePointsJacobians, refTriFacePoints, physCellVertices, cell);
                  //   ct::getPhysicalFaceNormals(triFacePointNormals, triFacePointsJacobians, faceOrd, cell);
                  //   ct::getPhysicalSideNormals(triSidePointNormals, triFacePointsJacobians, faceOrd, cell);

                  //   /*
                  //    * Compute face normals using direct linear parametrization of the face: the map from
                  //    * standard 2-simplex to physical Triangle<3> face in 3D is
                  //    * F(x,y) = V0 + (V1-V0)x + (V2-V0)*y
                  //    * Face normal is vector product Tx X Ty where Tx = (V1-V0); Ty = (V2-V0)
                  //    */
                  //   const auto v0ord = cell.getNodeMap(2, faceOrd, 0);
                  //   const auto v1ord = cell.getNodeMap(2, faceOrd, 1);
                  //   const auto v2ord = cell.getNodeMap(2, faceOrd, 2);

                  //   // Loop over face points: redundant for affine faces, but CellTools gives one vector
                  //   // per point so need to check all points anyways.
                  //   for (auto pt=0;pt<numTriFaceCubPoints;++pt) {
                  //     DynRankView ConstructWithLabel(tanX, 3);
                  //     DynRankView ConstructWithLabel(tanY, 3);
                  //     DynRankView ConstructWithLabel(faceNormal, 3);
                  //     for (auto d=0;d<cellDim;++d) {
                  //       tanX(d) = physCellVertices(0, v1ord, d) - physCellVertices(0, v0ord, d);
                  //       tanY(d) = physCellVertices(0, v2ord, d) - physCellVertices(0, v0ord, d);
                  //     }

                  //     rst::vecprod(faceNormal, tanX, tanY);

                  //     // Compare direct normal with d-component of the face/side normal by CellTools
                  //     for (auto d=0;d<cellDim;++d) {

                  //       // face normal method
                  //       if ( std::abs(faceNormal(d) - triFacePointNormals(0, pt, d)) > tol ){
                  //         errorFlag++;
                  //         *outStream
                  //           << std::setw(70) << "^^^^----FAILURE!" << "\n"
                  //           << " Face normal computation by CellTools failed for: \n"
                  //           << "       Cell Topology = " << cell.getName() << "\n"
                  //           << "       Face Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                  //           << "        Face ordinal = " << faceOrd << "\n"
                  //           << "   Face point number = " << pt << "\n"
                  //           << "   Normal coordinate = " << d  << "\n"
                  //           << "     CellTools value = " <<  triFacePointNormals(0, pt, d)
                  //           << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
                  //       }
                  //       //side normal method
                  //       if( std::abs(faceNormal(d) - triSidePointNormals(0, pt, d)) > tol ){
                  //         errorFlag++;
                  //         *outStream
                  //           << std::setw(70) << "^^^^----FAILURE!" << "\n"
                  //           << " Side normal computation by CellTools failed for: \n"
                  //           << "       Cell Topology = " << cell.getName() << "\n"
                  //           << "       Side Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                  //           << "        Side ordinal = " << faceOrd << "\n"
                  //           << "   Side point number = " << pt << "\n"
                  //           << "   Normal coordinate = " << d  << "\n"
                  //           << "     CellTools value = " <<  triSidePointNormals(0, pt, d)
                  //           << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
                  //       }
                  //     }
                  //   }
                  //   break;
                  // }
                case shards::Quadrilateral<4>::key: {
                  // Compute face normals using CellTools
                  ct::mapToReferenceSubcell(refQuadFacePoints, paramQuadFacePoints, 2, faceOrd, cell);
                  ct::setJacobian(quadFacePointsJacobians, refQuadFacePoints, physCellVertices, cell);
                  ct::getPhysicalFaceNormals(quadFacePointNormals, quadFacePointsJacobians, faceOrd, cell);
                  ct::getPhysicalSideNormals(quadSidePointNormals, quadFacePointsJacobians, faceOrd, cell);
                  /*
                   * Compute face normals using direct bilinear parametrization of the face: the map from
                   * [-1,1]^2 to physical Quadrilateral<4> face in 3D is
                   * F(x,y) = ((V0+V1+V2+V3) + (-V0+V1+V2-V3)*X + (-V0-V1+V2+V3)*Y + (V0-V1+V2-V3)*X*Y)/4
                   * Face normal is vector product Tx X Ty where
                   *          Tx = ((-V0+V1+V2-V3) + (V0-V1+V2-V3)*Y)/4
                   *          Ty = ((-V0-V1+V2+V3) + (V0-V1+V2-V3)*X)/4
                   */
                  const auto v0ord = cell.getNodeMap(2, faceOrd, 0);
                  const auto v1ord = cell.getNodeMap(2, faceOrd, 1);
                  const auto v2ord = cell.getNodeMap(2, faceOrd, 2);
                  const auto v3ord = cell.getNodeMap(2, faceOrd, 3);

                  // Loop over face points (redundant for affine faces, but needed for later when we handle non-affine ones)
                  for (auto pt=0;pt<numQuadFaceCubPoints;++pt) {
                    DynRankView ConstructWithLabel(tanX, 3);
                    DynRankView ConstructWithLabel(tanY, 3);
                    DynRankView ConstructWithLabel(faceNormal, 3);

                    for (size_type d=0;d<cellDim;++d) {
                      tanX(d) = ( physCellVertices(0, v0ord, d)*(-1.0 + paramQuadFacePoints(pt,1) )  +
                                  physCellVertices(0, v1ord, d)*( 1.0 - paramQuadFacePoints(pt,1) ) +
                                  physCellVertices(0, v2ord, d)*( 1.0 + paramQuadFacePoints(pt,1) ) +
                                  physCellVertices(0, v3ord, d)*(-1.0 - paramQuadFacePoints(pt,1) ) )/4.0;

                      tanY(d) = ( physCellVertices(0, v0ord, d)*(-1.0 + paramQuadFacePoints(pt,0) ) +
                                  physCellVertices(0, v1ord, d)*(-1.0 - paramQuadFacePoints(pt,0) ) +
                                  physCellVertices(0, v2ord, d)*( 1.0 + paramQuadFacePoints(pt,0) ) +
                                  physCellVertices(0, v3ord, d)*( 1.0 - paramQuadFacePoints(pt,0) ) )/4.0;
                    }

                    rst::vecprod(faceNormal, tanX, tanY);

                    // Compare direct normal with d-component of the face/side normal by CellTools
                    for (size_type d=0;d<cellDim;++d) {

                      // face normal method
                      if( std::abs(faceNormal(d) - quadFacePointNormals(0, pt, d)) > tol ) {
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Face normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Face Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Face ordinal = " << faceOrd << "\n"
                          << "   Face point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  quadFacePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
                      }
                      //side normal method
                      if ( std::abs(faceNormal(d) - quadSidePointNormals(0, pt, d)) > tol ) {
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Side normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Side Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Side ordinal = " << faceOrd << "\n"
                          << "   Side point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  quadSidePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
                      }
                    }
                  }
                  break;
                }
                default: {
                  errorFlag++;
                  *outStream << " Face normals test failure: face topology not supported \n\n";
                }
                }
              }
            }
          }
        }
      } catch (std::logic_error err) {
        //============================================================================================//
        // Wrap up test: check if the test broke down unexpectedly due to an exception                //
        //============================================================================================//
        *outStream << err.what() << "\n";
        errorFlag = -1000;
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);

      return errorFlag;
    }
  }
}


