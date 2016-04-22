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

        
    template<typename ValueType, typename DeviceSpaceType>
    int CellTools_Test02(const bool verbose) {
      typedef ValueType value_type;

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
      typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      int errorFlag = 0;

      
      // Vertices of the parametrization domain for 1-subcells: standard 1-cube [-1,1]
      DynRankView ConstructWithLabel(cube_1, 2, 1);
      cube_1(0,0) = -1.0; 
      cube_1(1,0) = 1.0;
      
  
      // Vertices of the parametrization domain for triangular faces: the standard 2-simplex
      DynRankView ConstructWithLabel(simplex_2, 3, 2);
      simplex_2(0, 0) = 0.0;   simplex_2(0, 1) = 0.0;
      simplex_2(1, 0) = 1.0;   simplex_2(1, 1) = 0.0;
      simplex_2(2, 0) = 0.0;   simplex_2(2, 1) = 1.0;
      
      
      // Vertices of the parametrization domain for quadrilateral faces: the standard 2-cube
      DynRankView ConstructWithLabel(cube_2, 4, 2);
      cube_2(0, 0) =  -1.0;    cube_2(0, 1) =  -1.0;
      cube_2(1, 0) =   1.0;    cube_2(1, 1) =  -1.0;
      cube_2(2, 0) =   1.0;    cube_2(2, 1) =   1.0;
      cube_2(3, 0) =  -1.0;    cube_2(3, 1) =   1.0;
      
  
      try {
        // // Pull all available topologies from Shards
        // std::vector<shards::CellTopology> allTopologies;
        // shards::getTopologies(allTopologies);
        
        // const auto topoSize = allTopologies.size();

        // /***********************************************************************************************
        //  *
        //  * Common for test 3 and 4: edge tangents and face normals for standard cells with base topo
        //  *
        //  **********************************************************************************************/
        
        // // Allocate storage and extract all standard cells with base topologies
        // std::vector<shards::CellTopology> standardBaseTopologies;    
        // shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);
        
        // // Define topologies for the edge and face parametrization domains. (faces are Tri or Quad)
        // CellTopology paramEdge    (shards::getCellTopologyData<shards::Line<2> >() );
        // CellTopology paramTriFace (shards::getCellTopologyData<shards::Triangle<3> >() );
        // CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
        
        // // Define CubatureFactory:
        // DefaultCubatureFactory<DeviceSpaceType>  cubFactory;   
        
        // *outStream 
        //   << "\n"
        //   << "===============================================================================\n" 
        //   << "| Test 4: face/side normals for stand. 3D cells with base topologies:         |\n" 
        //   << "===============================================================================\n\n";

        // // This test loops over standard 3D cells with base topologies, creates a set of nodes and tests normals 
        // {
        //   // Define cubature on the edge parametrization domain:
        //   const auto triFaceCubature  = cubFactory.create(paramTriFace, 6); 
        //   const auto quadFaceCubature = cubFactory.create(paramQuadFace, 6); 
          
        //   const auto faceCubDim           = triFaceCubature->getDimension();
        //   const auto numTriFaceCubPoints  = triFaceCubature->getNumPoints();
        //   const auto numQuadFaceCubPoints = quadFaceCubature->getNumPoints();    
          
        //   // Allocate storage for cubature points and weights on face parameter domain and fill with points:
        //   DynRankView ConstructWithLabel(paramTriFacePoints, numTriFaceCubPoints, faceCubDim);
        //   DynRankView ConstructWithLabel(paramTriFaceWeights, numTriFaceCubPoints);
        //   DynRankView ConstructWithLabel(paramQuadFacePoints, numQuadFaceCubPoints, faceCubDim);
        //   DynRankView ConstructWithLabel(paramQuadFaceWeights, numQuadFaceCubPoints);
          
        //   triFaceCubature->getCubature(paramTriFacePoints, paramTriFaceWeights);
        //   quadFaceCubature->getCubature(paramQuadFacePoints, paramQuadFaceWeights);
    
        //   // Loop over admissible topologies 
        //   for(cti = standardBaseTopologies.begin(); cti !=standardBaseTopologies.end(); ++cti){
            
        //     // Exclude 2D and Pyramid<5> cells
        //     if( ( (*cti).getDimension() == 3) && ( (*cti).getKey() != shards::Pyramid<5>::key) ){ 
              
        //       int cellDim = (*cti).getDimension();
        //       int vCount  = (*cti).getVertexCount();
        //       DynRankView ConstructWithLabel refCellVertices(vCount, cellDim);
        //       CellTools::getReferenceSubcellVertices(refCellVertices, cellDim, 0, (*cti) );
              
        //       *outStream << " Testing face/side normals for cell topology " <<  (*cti).getName() <<"\n";
              
        //       // Array for physical cell vertices ( must have rank 3 for setJacobians)
        //       DynRankView ConstructWithLabel physCellVertices(1, vCount, cellDim);
              
        //       // Randomize reference cell vertices by moving them up to +/- (1/8) units along their
        //       // coordinate axis. Guaranteed to be non-degenerate for standard cells with base topology 
        //       for(int v = 0; v < vCount; v++){
        //         for(int d = 0; d < cellDim; d++){
        //           double delta = Teuchos::ScalarTraits<double>::random()/8.0;
        //           physCellVertices(0, v, d) = refCellVertices(v, d) + delta;
        //         } //for d
        //       }// for v     
              
        //       // Allocate storage for cub. points on a ref. face; Jacobians, phys. face normals and 
        //       // benchmark normals.
        //       DynRankView ConstructWithLabel refTriFacePoints(numTriFaceCubPoints, cellDim);        
        //       DynRankView ConstructWithLabel refQuadFacePoints(numQuadFaceCubPoints, cellDim);        
        //       DynRankView ConstructWithLabel triFacePointsJacobians(1, numTriFaceCubPoints, cellDim, cellDim);
        //       DynRankView ConstructWithLabel quadFacePointsJacobians(1, numQuadFaceCubPoints, cellDim, cellDim);
        //       DynRankView ConstructWithLabel triFacePointNormals(1, numTriFaceCubPoints, cellDim);
        //       DynRankView ConstructWithLabel triSidePointNormals(1, numTriFaceCubPoints, cellDim);
        //       DynRankView ConstructWithLabel quadFacePointNormals(1, numQuadFaceCubPoints, cellDim);
        //       DynRankView ConstructWithLabel quadSidePointNormals(1, numQuadFaceCubPoints, cellDim);
              
              
        //       // Loop over faces:
        //       for(int faceOrd = 0; faceOrd < (int)(*cti).getSideCount(); faceOrd++){
                
        //         // This test presently includes only Triangle<3> and Quadrilateral<4> faces. Once we support
        //         // cells with extended topologies we will add their faces as well.
        //         switch( (*cti).getCellTopologyData(2, faceOrd) -> key ) {
                  
        //         case shards::Triangle<3>::key: 
        //           {
        //             // Compute face normals using CellTools
        //             CellTools::mapToReferenceSubcell(refTriFacePoints, paramTriFacePoints, 2, faceOrd, (*cti) );
        //             CellTools::setJacobian(triFacePointsJacobians, refTriFacePoints, physCellVertices, (*cti) );
        //             CellTools::getPhysicalFaceNormals(triFacePointNormals, triFacePointsJacobians, faceOrd, (*cti));               
        //             CellTools::getPhysicalSideNormals(triSidePointNormals, triFacePointsJacobians, faceOrd, (*cti));               
        //             /* 
        //              * Compute face normals using direct linear parametrization of the face: the map from
        //              * standard 2-simplex to physical Triangle<3> face in 3D is 
        //              * F(x,y) = V0 + (V1-V0)x + (V2-V0)*y 
        //              * Face normal is vector product Tx X Ty where Tx = (V1-V0); Ty = (V2-V0)
        //              */
        //             int v0ord = (*cti).getNodeMap(2, faceOrd, 0);
        //             int v1ord = (*cti).getNodeMap(2, faceOrd, 1);
        //             int v2ord = (*cti).getNodeMap(2, faceOrd, 2);
                    
        //             // Loop over face points: redundant for affine faces, but CellTools gives one vector 
        //             // per point so need to check all points anyways.
        //             for(int pt = 0; pt < numTriFaceCubPoints; pt++){
        //               DynRankView ConstructWithLabel tanX(3), tanY(3), faceNormal(3);
        //               for(int d = 0; d < cellDim; d++){
        //                 tanX(d) = (physCellVertices(0, v1ord, d) - physCellVertices(0, v0ord, d));
        //                 tanY(d) = (physCellVertices(0, v2ord, d) - physCellVertices(0, v0ord, d));
        //               }// for d
                      
        //               RealSpaceTools<double>::vecprod(faceNormal, tanX, tanY); 
                      
        //               // Compare direct normal with d-component of the face/side normal by CellTools
        //               for(int d = 0; d < cellDim; d++){
                        
        //                 // face normal method
        //                 if( abs(faceNormal(d) - triFacePointNormals(0, pt, d)) > INTREPID2_THRESHOLD ){
        //                   errorFlag++;
        //                   *outStream
        //                     << std::setw(70) << "^^^^----FAILURE!" << "\n"
        //                     << " Face normal computation by CellTools failed for: \n"
        //                     << "       Cell Topology = " << (*cti).getName() << "\n"
        //                     << "       Face Topology = " << (*cti).getCellTopologyData(2, faceOrd) -> name << "\n"
        //                     << "        Face ordinal = " << faceOrd << "\n"
        //                     << "   Face point number = " << pt << "\n"
        //                     << "   Normal coordinate = " << d  << "\n"
        //                     << "     CellTools value = " <<  triFacePointNormals(0, pt, d)
        //                     << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
        //                 }
        //                 //side normal method
        //                 if( abs(faceNormal(d) - triSidePointNormals(0, pt, d)) > INTREPID2_THRESHOLD ){
        //                   errorFlag++;
        //                   *outStream
        //                     << std::setw(70) << "^^^^----FAILURE!" << "\n"
        //                     << " Side normal computation by CellTools failed for: \n"
        //                     << "       Cell Topology = " << (*cti).getName() << "\n"
        //                     << "       Side Topology = " << (*cti).getCellTopologyData(2, faceOrd) -> name << "\n"
        //                     << "        Side ordinal = " << faceOrd << "\n"
        //                     << "   Side point number = " << pt << "\n"
        //                     << "   Normal coordinate = " << d  << "\n"
        //                     << "     CellTools value = " <<  triSidePointNormals(0, pt, d)
        //                     << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
        //                 }
        //               } // for d
        //             } // for pt
        //           }
        //           break;
                  
        //         case shards::Quadrilateral<4>::key:
        //           {
        //             // Compute face normals using CellTools
        //             CellTools::mapToReferenceSubcell(refQuadFacePoints, paramQuadFacePoints, 2, faceOrd, (*cti) );
        //             CellTools::setJacobian(quadFacePointsJacobians, refQuadFacePoints, physCellVertices, (*cti) );
        //             CellTools::getPhysicalFaceNormals(quadFacePointNormals, quadFacePointsJacobians, faceOrd, (*cti));               
        //             CellTools::getPhysicalSideNormals(quadSidePointNormals, quadFacePointsJacobians, faceOrd, (*cti)); 
        //             /*
        //              * Compute face normals using direct bilinear parametrization of the face: the map from
        //              * [-1,1]^2 to physical Quadrilateral<4> face in 3D is 
        //              * F(x,y) = ((V0+V1+V2+V3) + (-V0+V1+V2-V3)*X + (-V0-V1+V2+V3)*Y + (V0-V1+V2-V3)*X*Y)/4 
        //              * Face normal is vector product Tx X Ty where
        //              *          Tx = ((-V0+V1+V2-V3) + (V0-V1+V2-V3)*Y)/4
        //              *          Ty = ((-V0-V1+V2+V3) + (V0-V1+V2-V3)*X)/4
        //              */
        //             int v0ord = (*cti).getNodeMap(2, faceOrd, 0);
        //             int v1ord = (*cti).getNodeMap(2, faceOrd, 1);
        //             int v2ord = (*cti).getNodeMap(2, faceOrd, 2);
        //             int v3ord = (*cti).getNodeMap(2, faceOrd, 3);
                    
        //             // Loop over face points (redundant for affine faces, but needed for later when we handle non-affine ones)
        //             for(int pt = 0; pt < numTriFaceCubPoints; pt++){
        //               DynRankView ConstructWithLabel tanX(3), tanY(3), faceNormal(3);
        //               for(int d = 0; d < cellDim; d++){
        //                 tanX(d) = (physCellVertices(0, v0ord, d)*(-1.0 + paramQuadFacePoints(pt,1) )  +
        //                            physCellVertices(0, v1ord, d)*( 1.0 - paramQuadFacePoints(pt,1) ) + 
        //                            physCellVertices(0, v2ord, d)*( 1.0 + paramQuadFacePoints(pt,1) ) + 
        //                            physCellVertices(0, v3ord, d)*(-1.0 - paramQuadFacePoints(pt,1) ) )/4.0;
                        
        //                 tanY(d) = (physCellVertices(0, v0ord, d)*(-1.0 + paramQuadFacePoints(pt,0) ) +
        //                            physCellVertices(0, v1ord, d)*(-1.0 - paramQuadFacePoints(pt,0) ) + 
        //                            physCellVertices(0, v2ord, d)*( 1.0 + paramQuadFacePoints(pt,0) ) + 
        //                            physCellVertices(0, v3ord, d)*( 1.0 - paramQuadFacePoints(pt,0) ) )/4.0;
        //               }// for d
                      
        //               RealSpaceTools<double>::vecprod(faceNormal, tanX, tanY); 
        //               // Compare direct normal with d-component of the face/side normal by CellTools
        //               for(int d = 0; d < cellDim; d++){
                        
        //                 // face normal method
        //                 if( abs(faceNormal(d) - quadFacePointNormals(0, pt, d)) > INTREPID2_THRESHOLD ){
        //                   errorFlag++;
        //                   *outStream
        //                     << std::setw(70) << "^^^^----FAILURE!" << "\n"
        //                     << " Face normal computation by CellTools failed for: \n"
        //                     << "       Cell Topology = " << (*cti).getName() << "\n"
        //                     << "       Face Topology = " << (*cti).getCellTopologyData(2, faceOrd) -> name << "\n"
        //                     << "        Face ordinal = " << faceOrd << "\n"
        //                     << "   Face point number = " << pt << "\n"
        //                     << "   Normal coordinate = " << d  << "\n"
        //                     << "     CellTools value = " <<  quadFacePointNormals(0, pt, d)
        //                     << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
        //                 }
        //                 //side normal method
        //                 if( abs(faceNormal(d) - quadSidePointNormals(0, pt, d)) > INTREPID2_THRESHOLD ){
        //                   errorFlag++;
        //                   *outStream
        //                     << std::setw(70) << "^^^^----FAILURE!" << "\n"
        //                     << " Side normal computation by CellTools failed for: \n"
        //                     << "       Cell Topology = " << (*cti).getName() << "\n"
        //                     << "       Side Topology = " << (*cti).getCellTopologyData(2, faceOrd) -> name << "\n"
        //                     << "        Side ordinal = " << faceOrd << "\n"
        //                     << "   Side point number = " << pt << "\n"
        //                     << "   Normal coordinate = " << d  << "\n"
        //                     << "     CellTools value = " <<  quadSidePointNormals(0, pt, d)
        //                     << "     Benchmark value = " <<  faceNormal(d) << "\n\n";
        //                 }
        //               } // for d
        //             }// for pt
        //           }// case Quad
        //           break;
        //         default:
        //           errorFlag++;
        //           *outStream << " Face normals test failure: face topology not supported \n\n";
        //         } // switch 
        //       }// for faceOrd
        //     }// if admissible
        //   }// for cti
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
    

