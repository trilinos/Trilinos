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

//#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"
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
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

    template<typename ValueType, typename DeviceType>
    int CellTools_Test03(const bool verbose) {
      using ExecSpaceType = typename DeviceType::execution_space;

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

      using HostDeviceType = Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>;
      using ct = CellTools<DeviceType>;
      using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;
      using DynRankViewHost = Kokkos::DynRankView<ValueType,Kokkos::HostSpace>;

      const ValueType tol = tolerence<ValueType>()*100.0;

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
        const ordinal_type topoSize = standardBaseTopologies.size();

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
          const ordinal_type testAccuracy = 6;
          CubatureDirectLineGauss<DeviceType,ValueType,ValueType> edgeCubature(testAccuracy);
          CubatureDirectTriDefault<DeviceType,ValueType,ValueType> triFaceCubature(testAccuracy);
          CubatureTensor<DeviceType,ValueType,ValueType> quadFaceCubature( edgeCubature, edgeCubature );

          const ordinal_type faceCubDim = 2;
          const ordinal_type numTriFaceCubPoints  = triFaceCubature.getNumPoints();
          const ordinal_type numQuadFaceCubPoints = quadFaceCubature.getNumPoints();

          // Allocate storage for cubature points and weights on face parameter domain and fill with points:
          DynRankView ConstructWithLabel(paramTriFacePoints,   numTriFaceCubPoints,  faceCubDim);
          DynRankView ConstructWithLabel(paramTriFaceWeights,  numTriFaceCubPoints);
          DynRankView ConstructWithLabel(paramQuadFacePoints,  numQuadFaceCubPoints, faceCubDim);
          DynRankView ConstructWithLabel(paramQuadFaceWeights, numQuadFaceCubPoints);

          triFaceCubature.getCubature(paramTriFacePoints, paramTriFaceWeights);
          quadFaceCubature.getCubature(paramQuadFacePoints, paramQuadFaceWeights);
          // create mirror host view
          auto hParamQuadFacePoints = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), paramQuadFacePoints);

          // Loop over admissible topologies
          for (ordinal_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];

            // skip cells not supported
            if (!ct::hasReferenceCell(cell))
              continue;

            // Exclude 2D and Pyramid<5> cells
            if ( cell.getDimension() == 3 ) { // && cell.getKey() != shards::Pyramid<5>::key ) {
              *outStream << " Testing face/side normals for cell topology " <<  cell.getName() <<"\n";

              const ordinal_type cellDim  = cell.getDimension();
              const ordinal_type nCount   = cell.getNodeCount();
              const ordinal_type vCount   = cell.getVertexCount();

              DynRankView ConstructWithLabel(refCellVertices, nCount, cellDim);
              ct::getReferenceSubcellVertices(refCellVertices, cellDim, 0, cell);

              // create mirror host view
              auto hRefCellVertices = Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace(), refCellVertices);

              // Array for physical cell vertices ( must have rank 3 for setJacobians)
              DynRankView ConstructWithLabel(physCellVertices, 1, vCount, cellDim);
              auto hPhysCellVertices = Kokkos::create_mirror_view(physCellVertices);

              // Randomize reference cell vertices by moving them up to +/- (1/8) units along their
              // coordinate axis. Guaranteed to be non-degenerate for standard cells with base topology
              for (ordinal_type v=0;v<vCount;++v)
                for (ordinal_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<double>::random()/8.0;
                  hPhysCellVertices(0, v, d) = hRefCellVertices(v, d) + delta;
                }
              Kokkos::deep_copy(physCellVertices,hPhysCellVertices);

              // Allocate storage for cub. points on a ref. face; Jacobians, phys. face normals and
              // benchmark normals.
              DynRankView ConstructWithLabel(refTriFacePoints, numTriFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(triFacePointsJacobians, 1, numTriFaceCubPoints, cellDim, cellDim);
              DynRankView ConstructWithLabel(triFacePointNormals,    1, numTriFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(triSidePointNormals,    1, numTriFaceCubPoints, cellDim);
              auto hTriFacePointNormals = Kokkos::create_mirror_view(triFacePointNormals);
              auto hTriSidePointNormals = Kokkos::create_mirror_view(triSidePointNormals);

              DynRankView ConstructWithLabel(refQuadFacePoints,          numQuadFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(quadFacePointsJacobians, 1, numQuadFaceCubPoints, cellDim, cellDim);
              DynRankView ConstructWithLabel(quadFacePointNormals,    1, numQuadFaceCubPoints, cellDim);
              DynRankView ConstructWithLabel(quadSidePointNormals,    1, numQuadFaceCubPoints, cellDim);
              auto hQuadFacePointNormals = Kokkos::create_mirror_view(quadFacePointNormals);
              auto hQuadSidePointNormals = Kokkos::create_mirror_view(quadSidePointNormals);
              
              // Loop over faces:
              for (size_type faceOrd=0;faceOrd<cell.getSideCount();++faceOrd) {
                
                // This test presently includes only Triangle<3> and Quadrilateral<4> faces. Once we support
                // cells with extended topologies we will add their faces as well.
                switch ( cell.getCellTopologyData(2, faceOrd)->key ) {
                case shards::Triangle<3>::key: {
                  // Compute face normals using CellTools
                  ct::mapToReferenceSubcell(refTriFacePoints, paramTriFacePoints, 2, faceOrd, cell);
                  ct::setJacobian(triFacePointsJacobians, refTriFacePoints, physCellVertices, cell);
                  ct::getPhysicalFaceNormals(triFacePointNormals, triFacePointsJacobians, faceOrd, cell);
                  ct::getPhysicalSideNormals(triSidePointNormals, triFacePointsJacobians, faceOrd, cell);
                  Kokkos::deep_copy(hTriFacePointNormals,triFacePointNormals);
                  Kokkos::deep_copy(hTriSidePointNormals,triSidePointNormals);
                  /*
                   * Compute face normals using direct linear parametrization of the face: the map from
                   * standard 2-simplex to physical Triangle<3> face in 3D is
                   * F(x,y) = V0 + (V1-V0)x + (V2-V0)*y
                   * Face normal is vector product Tx X Ty where Tx = (V1-V0); Ty = (V2-V0)
                   */
                  const auto v0ord = cell.getNodeMap(2, faceOrd, 0);
                  const auto v1ord = cell.getNodeMap(2, faceOrd, 1);
                  const auto v2ord = cell.getNodeMap(2, faceOrd, 2);
                  
                  // Loop over face points: redundant for affine faces, but CellTools gives one vector
                  // per point so need to check all points anyways.
                  for (ordinal_type pt=0;pt<numTriFaceCubPoints;++pt) {
                    DynRankViewHost ConstructWithLabel(hTanX, 3);
                    DynRankViewHost ConstructWithLabel(hTanY, 3);
                    DynRankViewHost ConstructWithLabel(hFaceNormal, 3);
                    for (ordinal_type d=0;d<cellDim;++d) {
                      hTanX(d) = hPhysCellVertices(0, v1ord, d) - hPhysCellVertices(0, v0ord, d);
                      hTanY(d) = hPhysCellVertices(0, v2ord, d) - hPhysCellVertices(0, v0ord, d);
                    }
                    
                    RealSpaceTools<HostDeviceType>::vecprod(hFaceNormal, hTanX, hTanY);

                    ExecSpaceType().fence();

                    // Compare direct normal with d-component of the face/side normal by CellTools
                    for (ordinal_type d=0;d<cellDim;++d) {
                      
                      // face normal method
                      if ( std::abs(hFaceNormal(d) - hTriFacePointNormals(0, pt, d)) > tol ){
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Face normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Face Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Face ordinal = " << faceOrd << "\n"
                          << "   Face point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  hTriFacePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  hFaceNormal(d) << "\n\n";
                      }
                      //side normal method
                      if( std::abs(hFaceNormal(d) - hTriSidePointNormals(0, pt, d)) > tol ){
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Side normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Side Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Side ordinal = " << faceOrd << "\n"
                          << "   Side point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  hTriSidePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  hFaceNormal(d) << "\n\n";
                      }
                    }
                  }
                  break;
                }
                case shards::Quadrilateral<4>::key: {
                  // Compute face normals using CellTools
                  ct::mapToReferenceSubcell(refQuadFacePoints, paramQuadFacePoints, 2, faceOrd, cell);
                  ct::setJacobian(quadFacePointsJacobians, refQuadFacePoints, physCellVertices, cell);
                  ct::getPhysicalFaceNormals(quadFacePointNormals, quadFacePointsJacobians, faceOrd, cell);
                  ct::getPhysicalSideNormals(quadSidePointNormals, quadFacePointsJacobians, faceOrd, cell);
                  Kokkos::deep_copy(hQuadFacePointNormals,quadFacePointNormals);
                  Kokkos::deep_copy(hQuadSidePointNormals,quadSidePointNormals);
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
                  for (ordinal_type pt=0;pt<numQuadFaceCubPoints;++pt) {
                    DynRankViewHost ConstructWithLabel(hTanX, 3);
                    DynRankViewHost ConstructWithLabel(hTanY, 3);
                    DynRankViewHost ConstructWithLabel(hFaceNormal, 3);

                    for (ordinal_type d=0;d<cellDim;++d) {
                      hTanX(d) = ( hPhysCellVertices(0, v0ord, d)*(-1.0 + hParamQuadFacePoints(pt,1) )  +
                                  hPhysCellVertices(0, v1ord, d)*( 1.0 - hParamQuadFacePoints(pt,1) ) +
                                  hPhysCellVertices(0, v2ord, d)*( 1.0 + hParamQuadFacePoints(pt,1) ) +
                                  hPhysCellVertices(0, v3ord, d)*(-1.0 - hParamQuadFacePoints(pt,1) ) )/4.0;

                      hTanY(d) = ( hPhysCellVertices(0, v0ord, d)*(-1.0 + hParamQuadFacePoints(pt,0) ) +
                                  hPhysCellVertices(0, v1ord, d)*(-1.0 - hParamQuadFacePoints(pt,0) ) +
                                  hPhysCellVertices(0, v2ord, d)*( 1.0 + hParamQuadFacePoints(pt,0) ) +
                                  hPhysCellVertices(0, v3ord, d)*( 1.0 - hParamQuadFacePoints(pt,0) ) )/4.0;
                    }

                    RealSpaceTools<HostDeviceType>::vecprod(hFaceNormal, hTanX, hTanY);

                    ExecSpaceType().fence();

                    // Compare direct normal with d-component of the face/side normal by CellTools
                    for (ordinal_type d=0;d<cellDim;++d) {

                      // face normal method
                      if( std::abs(hFaceNormal(d) - hQuadFacePointNormals(0, pt, d)) > tol ) {
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Face normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Face Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Face ordinal = " << faceOrd << "\n"
                          << "   Face point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  hQuadFacePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  hFaceNormal(d) << "\n\n";
                      }
                      //side normal method
                      if ( std::abs(hFaceNormal(d) - hQuadSidePointNormals(0, pt, d)) > tol ) {
                        errorFlag++;
                        *outStream
                          << std::setw(70) << "^^^^----FAILURE!" << "\n"
                          << " Side normal computation by CellTools failed for: \n"
                          << "       Cell Topology = " << cell.getName() << "\n"
                          << "       Side Topology = " << cell.getCellTopologyData(2, faceOrd) -> name << "\n"
                          << "        Side ordinal = " << faceOrd << "\n"
                          << "   Side point number = " << pt << "\n"
                          << "   Normal coordinate = " << d  << "\n"
                          << "     CellTools value = " <<  hQuadSidePointNormals(0, pt, d)
                          << "     Benchmark value = " <<  hFaceNormal(d) << "\n\n";
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
      } catch (std::logic_error &err) {
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


