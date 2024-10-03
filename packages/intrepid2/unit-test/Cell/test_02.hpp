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
    int CellTools_Test02(const bool verbose) {

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
  
      using ct = CellTools<DeviceType>;
      using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

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
        const auto topoSize = standardBaseTopologies.size();
        
        // Define topologies for the edge and face parametrization domains. (faces are Tri or Quad)
        const auto paramEdge     = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
        const auto paramTriFace  = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
        const auto paramQuadFace = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
        
        *outStream
          << "\n"
          << "===============================================================================\n" 
          << "| Test 3: edge tangents/normals for stand. cells with base topologies:        |\n" 
          << "===============================================================================\n\n";

        {
          // Define cubature on the edge parametrization domain:
          const auto testAccuracy = 6;
          CubatureDirectLineGauss<DeviceType,ValueType,ValueType> edgeCubature(testAccuracy);

          const auto cubDim       = edgeCubature.getDimension();
          const auto numCubPoints = edgeCubature.getNumPoints();
          
          // Allocate storage for cubature points and weights on edge parameter domain and fill with points:
          DynRankView ConstructWithLabel(paramEdgePoints,  numCubPoints, cubDim);
          DynRankView ConstructWithLabel(paramEdgeWeights, numCubPoints);

          edgeCubature.getCubature(paramEdgePoints, paramEdgeWeights);
          
          // Loop over admissible topologies 
          for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
            const auto cell = standardBaseTopologies[topoOrd];

            // skip cells not supported
            if (!ct::hasReferenceCell(cell))
              continue;
            
            // Exclude 0D (node), 1D (Line) 
            if ( cell.getDimension() >= 2 ) { //&& cell.getKey() != shards::Pyramid<5>::key ) { 
              const auto cellDim  = cell.getDimension();
              const auto nCount   = cell.getNodeCount();
              const auto vCount   = cell.getVertexCount();

              DynRankView ConstructWithLabel(refCellVertices, nCount, cellDim);
              ct::getReferenceSubcellVertices(refCellVertices, cellDim, 0, cell);
              
              *outStream << " Testing edge tangents (and normals for cellDim = 2) for " 
                         <<  cell.getName() << " cellDim = " << cellDim <<"\n";
                            
              // Array for physical cell vertices ( must have rank 3 for setJacobians)
              DynRankView ConstructWithLabel(physCellVertices, 1, vCount, cellDim);
              
              // create mirror host views
              auto hRefCellVertices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), refCellVertices);

              auto hPhysCellVertices = Kokkos::create_mirror_view(physCellVertices);

              // Randomize reference cell vertices by moving them up to +/- (1/8) units along their
              // coordinate axis. Guaranteed to be non-degenerate for standard cells with base topology 
              for (size_type v=0;v<vCount;++v) 
                for (size_type d=0;d<cellDim;++d) {
                  const auto delta = Teuchos::ScalarTraits<ValueType>::random()/8.0;
                  hPhysCellVertices(0, v, d) = hRefCellVertices(v, d) + delta;
                }

              Kokkos::deep_copy(physCellVertices,hPhysCellVertices);

        
              // Allocate storage for cub. points on a ref. edge; Jacobians, phys. edge tangents/normals
              DynRankView ConstructWithLabel(refEdgePoints,          numCubPoints, cellDim);        

              // here, 1 means that the container includes a single cell
              DynRankView ConstructWithLabel(edgePointsJacobians, 1, numCubPoints, cellDim, cellDim);
              DynRankView ConstructWithLabel(edgePointTangents,   1, numCubPoints, cellDim);
              DynRankView ConstructWithLabel(edgePointNormals,    1, numCubPoints, cellDim);

              // create mirror host views
              auto hEdgePointTangents = Kokkos::create_mirror_view(edgePointTangents);
              auto hEdgePointNormals = Kokkos::create_mirror_view(edgePointNormals);
              
              // Loop over edges:
              for (size_type edgeOrd=0;edgeOrd<cell.getEdgeCount();++edgeOrd) {
                /* 
                 * Compute tangents on the specified physical edge using CellTools:
                 *    1. Map points from edge parametrization domain to ref. edge with specified ordinal
                 *    2. Compute parent cell Jacobians at ref. edge points
                 *    3. Compute physical edge tangents
                 */
                ct::mapToReferenceSubcell(refEdgePoints, paramEdgePoints, 1, edgeOrd, cell);
                ct::setJacobian(edgePointsJacobians, refEdgePoints, physCellVertices, cell);
                ct::getPhysicalEdgeTangents(edgePointTangents, edgePointsJacobians, edgeOrd, cell); 
                Kokkos::deep_copy( hEdgePointTangents, edgePointTangents);

                // Compute side normals (here only for sides of dim 2):
                if(cellDim == 2) {
                  ct::getPhysicalSideNormals(edgePointNormals, edgePointsJacobians, edgeOrd, cell);
                  Kokkos::deep_copy( hEdgePointNormals, edgePointNormals);
                }

                /*
                 * Compute tangents directly using parametrization of phys. edge and compare with CellTools tangents.
                 *    1. Get edge vertices
                 *    2. For affine edges tangent coordinates are given by F'(t) = (V1-V0)/2
                 *       (for now we only test affine edges, but later we will test edges for cells 
                 *        with extended topologies.)
                 */
                const auto v0ord = cell.getNodeMap(1, edgeOrd, 0);
                const auto v1ord = cell.getNodeMap(1, edgeOrd, 1);
                
                for (auto pt=0;pt<numCubPoints;++pt) {
                  
                  // Temp storage for directly computed edge tangents
                  Kokkos::DynRankView<ValueType,Kokkos::HostSpace> ConstructWithLabel(hEdgeBenchmarkTangents, 3);
                  
                  for (size_type d=0;d<cellDim;++d) {
                    hEdgeBenchmarkTangents(d) = (hPhysCellVertices(0, v1ord, d) - hPhysCellVertices(0, v0ord, d))/2.0;
                    
                    // Compare with d-component of edge tangent by CellTools
                    if ( std::abs(hEdgeBenchmarkTangents(d) - hEdgePointTangents(0, pt, d)) > tol ) {
                      errorFlag++;
                      *outStream
                        << std::setw(70) << "^^^^----FAILURE!" << "\n"
                        << " Edge tangent computation by CellTools failed for: \n"
                        << "       Cell Topology = " << cell.getName() << "\n"
                        << "        Edge ordinal = " << edgeOrd << "\n"
                        << "   Edge point number = " << pt << "\n"
                        << "  Tangent coordinate = " << d << "\n"
                        << "     CellTools value = " <<  hEdgePointTangents(0, pt, d) << "\n"
                        << "     Benchmark value = " <<  hEdgeBenchmarkTangents(d) << "\n\n";
                    }
                  } // for d
                  
                  // Test side normals for 2D cells only: edge normal has coordinates (t1, -t0)
                  if (cellDim == 2) {
                    if ( std::abs(hEdgeBenchmarkTangents(1) - hEdgePointNormals(0, pt, 0)) > tol ) {
                      errorFlag++;
                      *outStream
                        << std::setw(70) << "^^^^----FAILURE!" << "\n"
                        << " Edge Normal computation by CellTools failed for: \n"
                        << "       Cell Topology = " << cell.getName() << "\n"
                        << "        Edge ordinal = " << edgeOrd << "\n"
                        << "   Edge point number = " << pt << "\n"
                        << "   Normal coordinate = " << 0 << "\n"
                        << "     CellTools value = " <<  hEdgePointNormals(0, pt, 0) << "\n"
                        << "     Benchmark value = " <<  hEdgeBenchmarkTangents(1) << "\n\n";
                    }
                    if ( std::abs(hEdgeBenchmarkTangents(0) + hEdgePointNormals(0, pt, 1)) > tol ) {
                      errorFlag++;
                      *outStream
                        << std::setw(70) << "^^^^----FAILURE!" << "\n"
                        << " Edge Normal computation by CellTools failed for: \n"
                        << "       Cell Topology = " << cell.getName() << "\n"
                        << "        Edge ordinal = " << edgeOrd << "\n"
                        << "   Edge point number = " << pt << "\n"
                        << "   Normal coordinate = " << 1  << "\n"
                        << "     CellTools value = " <<  hEdgePointNormals(0, pt, 1) << "\n"
                        << "     Benchmark value = " << -hEdgeBenchmarkTangents(0) << "\n\n";
                    }
                  } // edge normals            
                } // for pt
              }// for edgeOrd
            }// if admissible cell
          }// topOrd
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
    

