// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Example of the CellTools class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson
*/

#include "Intrepid_CellTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Shards_CellTopology.hpp"

#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;
using namespace shards;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  typedef CellTools<double>       CellTools;
  typedef shards::CellTopology    CellTopology;
  
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the CellTools class                        |\n" \
  << "|                                                                             |\n" \
  << "|  1) Reference edge parametrizations                                         |\n" \
  << "|  2) Reference face parametrizations                                         |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n" \
  << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| Summary:                                                                    |\n"\
  << "| Reference edge parametrizations map [-1,1] to the edges of reference cells. |\n"\
  << "| They are used to define, e.g., integration points on the edges of 2D and 3D |\n"\
  << "| reference cells. Edge parametrizations for special 2D  cells such as Beam   |\n"\
  << "| and ShellLine, are also supported.                                          |\n"\
  << "===============================================================================\n";
 
  /* 
    Specification of integration points on 1-subcells (edges) of reference cells. Edges are 
    parametrized by [-1,1] and integration points on an edge are defined by mapping integration
    points from the parametrization domain [-1,1] to a specific edge on the reference cell.
   
   1. Common tasks: definition of integration points in the edge parametrization domain [-1,1]
      These steps are independent of parent cell topology:
   
      a. Instantiate a CubatureFactory object to create cubatures (needed for face maps too)
      b. Define parametrization domain for the edges as having Line<2> cell topology. This is 
         required by the CubatureFactory in order to select cubature points and weights from 
         the reference line [-1,1]
      c. Use CubatureFactory to select cubature of the desired degree for the Line<2> topology
      d. Allocate containers for the cubature points and weights.
   
   2. Parent cell topology specific tasks
   
      a. Select the parent cell topology
      b. Allocate containers for the images of the integration points on [-1,1] on the edges
      c. Apply the edge parametrization map to the pointss in [-1,1]
   */
  
  // Step 1.a (Define CubatureFactory)
  DefaultCubatureFactory<double>  cubFactory;   
  
  
  // Step 1.b (Define the topology of the parametrization domain)
  CellTopology edgeParam(shards::getCellTopologyData<shards::Line<2> >() );
  
  
  // Step 1.c (selects Gauss rule of order 6 on [-1,1]) 
  Teuchos::RCP<Cubature<double> > edgeParamCubature = cubFactory.create(edgeParam, 6); 
  
  
  // Step 1.d (allocate storage for cubature points)
  int cubDim       = edgeParamCubature -> getDimension();
  int numCubPoints = edgeParamCubature -> getNumPoints();
  
  FieldContainer<double> edgeParamCubPts(numCubPoints, cubDim);
  FieldContainer<double> edgeParamCubWts(numCubPoints);
  edgeParamCubature -> getCubature(edgeParamCubPts, edgeParamCubWts);
    
  
  
  std::cout \
  << "===============================================================================\n"\
  << "| EXAMPLE 1.1                                                                 |\n"
  << "| Edge parametrizations for standard 2D cells: Triangle                       |\n"\
  << "===============================================================================\n";
    
  // Step 2.a (select reference cell topology)
  CellTopology triangle_3(getCellTopologyData<Triangle<3> >() );
  
  
  // Step 2.b (allocate storage for points on an edge of the reference cell)
  FieldContainer<double> triEdgePoints(numCubPoints, triangle_3.getDimension() );
    
  
  // Step 2.c (same points are mapped to all edges, can also map different set to each edge) 
  for(int edgeOrd = 0; edgeOrd < (int)triangle_3.getEdgeCount(); edgeOrd++){
    
    CellTools::mapToReferenceSubcell(triEdgePoints, 
                                     edgeParamCubPts,
                                     1,
                                     edgeOrd,
                                     triangle_3);
    
    // Optional: print the vertices of the reference edge
    CellTools::printSubcellVertices(1, edgeOrd, triangle_3);
        
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t Parameter point " 
      << std::setw(12) << std::right << edgeParamCubPts(pt, 0) << std::setw(10) << "  -->  " << "(" 
      << std::setw(10) << std::right << triEdgePoints(pt, 0) << " , " 
      << std::setw(10) << std::right << triEdgePoints(pt, 1) << ")\n";
    }    
    std::cout << "\n";
  }
  
  
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1.2                                                                 |\n"
    << "| Edge parametrizations for standard 2D cells: Quadrilateral                  |\n"\
    << "===============================================================================\n";
  
  // Step 2.a (select reference cell topology)
  CellTopology quad_4(getCellTopologyData<Quadrilateral<4> >() );
  
  
  // Step 2.b (allocate storage for points on an edge of the reference cell)
  FieldContainer<double> quadEdgePoints(numCubPoints, quad_4.getDimension() );
  
  
  // Step 2.c (same points are mapped to all edges, can also map different set to each edge) 
  for(int edgeOrd = 0; edgeOrd < (int)quad_4.getEdgeCount(); edgeOrd++){
    
    CellTools::mapToReferenceSubcell(quadEdgePoints, 
                                     edgeParamCubPts,
                                     1,
                                     edgeOrd,
                                     quad_4);
    
    // Optional: print the vertices of the reference edge
    CellTools::printSubcellVertices(1, edgeOrd, quad_4);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t Parameter point " 
      << std::setw(12) << std::right << edgeParamCubPts(pt, 0) << std::setw(10) << "  -->  " << "(" 
      << std::setw(10) << std::right << quadEdgePoints(pt, 0) << " , " 
      << std::setw(10) << std::right << quadEdgePoints(pt, 1) << ")\n";
    }    
    std::cout << "\n";
  }
  
  
  /* 
    Specification of integration points on 2-subcells (faces) of reference cells. Reference cells
    can have triangular, quadrilateral or a mixture of triangular and quadrilateral faces. Thus,
    parametrization domain of a face depends on that face's topology and is either the standard 
    2-simplex {(0,0), (1,0), (0,1)} for triangular faces or the standard 2-cube [-1,1]^2 for 
    quadrilateral faces. 
      
   1. Common tasks: definition of integration points in the standard 2-simplex and the standard
      2-cube. These steps are independent of parent cell topology:
   
      a. Instantiate a CubatureFactory object to create cubatures (already done!)
      b. Define parametrization domain for traingular faces as having Triangle<3> topology and for
         quadrilateral faces - as having Quadrilateral<4> topology. This is required by the 
         CubatureFactory in order to select cubature points and weights from the appropriate
         face parametrization domain. 
      c. Use CubatureFactory to select cubature of the desired degree for Triangle<3> and 
         Quadrilateral<4> topologies
      d. Allocate containers for the cubature points and weights on the parametrization domains.
   
   2. Parent cell topology specific tasks
   
      a. Select the parent cell topology
      b. Allocate containers for the images of the integration points from the parametrization
         domains on the reference faces
      c. Apply the face parametrization map to the points in the parametrization domain
   */
  
  // Step 1.b (Define the topology of the parametrization domain)
  CellTopology triFaceParam(shards::getCellTopologyData<shards::Triangle<3> >() );
  CellTopology quadFaceParam(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  
  // Step 1.c (selects Gauss rule of order 3 on [-1,1]^2 and a rule of order 3 on Triangle) 
  Teuchos::RCP<Cubature<double> > triFaceParamCubature = cubFactory.create(triFaceParam, 3); 
  Teuchos::RCP<Cubature<double> > quadFaceParamCubature = cubFactory.create(quadFaceParam, 3); 
  
  
  // Step 1.d - Triangle faces (allocate storage for cubature points)
  int triFaceCubDim    = triFaceParamCubature -> getDimension();
  int triFaceNumCubPts = triFaceParamCubature -> getNumPoints();
  
  FieldContainer<double> triFaceParamCubPts(triFaceNumCubPts, triFaceCubDim);
  FieldContainer<double> triFaceParamCubWts(triFaceNumCubPts);
  triFaceParamCubature -> getCubature(triFaceParamCubPts, triFaceParamCubWts);
  
  
  // Step 1.d - Quadrilateral faces (allocate storage for cubature points)
  int quadFaceCubDim    = quadFaceParamCubature -> getDimension();
  int quadFaceNumCubPts = quadFaceParamCubature -> getNumPoints();
  
  FieldContainer<double> quadFaceParamCubPts(quadFaceNumCubPts, quadFaceCubDim);
  FieldContainer<double> quadFaceParamCubWts(quadFaceNumCubPts);
  quadFaceParamCubature -> getCubature(quadFaceParamCubPts, quadFaceParamCubWts);
  
  
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 2.1                                                                 |\n"
    << "| Face parametrizations for standard 3D cells: Tetrahedron                    |\n"\
    << "===============================================================================\n";
  
  // Step 2.a (select reference cell topology)
  CellTopology tet_4(getCellTopologyData<Tetrahedron<4> >() );
  
  
  // Step 2.b (allocate storage for points on a face of the reference cell)
  FieldContainer<double> tetFacePoints(triFaceNumCubPts, tet_4.getDimension() );
  
  
  // Step 2.c (same points are mapped to all faces, can also map different set to each face) 
  for(int faceOrd = 0; faceOrd < (int)tet_4.getSideCount(); faceOrd++){
    
    CellTools::mapToReferenceSubcell(tetFacePoints, 
                                     triFaceParamCubPts,
                                     2,
                                     faceOrd,
                                     tet_4);
    
    // Optional: print the vertices of the reference face 
    CellTools::printSubcellVertices(2, faceOrd, tet_4);
    
    for(int pt = 0; pt < triFaceNumCubPts; pt++){
      std::cout << "\t Parameter point (" 
      << std::setw(10) << std::right <<  triFaceParamCubPts(pt, 0) << " , "
      << std::setw(10) << std::right <<  triFaceParamCubPts(pt, 1) << ")  " 
      << std::setw(10) << " -->  " << "(" 
      << std::setw(10) << std::right << tetFacePoints(pt, 0) << " , " 
      << std::setw(10) << std::right << tetFacePoints(pt, 1) << " , " 
      << std::setw(10) << std::right << tetFacePoints(pt, 2) << ")\n";
    }    
    std::cout << "\n";
  }
  
  
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 2.2                                                                 |\n"
    << "| Face parametrizations for standard 3D cells: Wedge                          |\n"\
    << "| Example of a reference cell that has two different kinds of faces           |\n"\
    << "===============================================================================\n";
  
  
  // Step 2.a (select reference cell topology)
  CellTopology wedge_6(getCellTopologyData<Wedge<6> >() );
  
  
  // Step 2.b (allocate storage for points on a face of the reference cell)
  //          Wedge<6> has Triangle<3> and Quadrilateral<4> faces. Two different arrays are needed
  //          to store the points on each face because different types integration rules are used
  //          on their respective parametrization domains and numbers of points defined by these
  //          rules do not necessarily match.
  FieldContainer<double> wedgeTriFacePoints(triFaceNumCubPts, wedge_6.getDimension() );
  FieldContainer<double> wedgeQuadFacePoints(quadFaceNumCubPts, wedge_6.getDimension() );

  
  // Step 2.c (for Wedge<6> need to distinguish Triangle<3> and Quadrilateral<4> faces) 
  for(int faceOrd = 0; faceOrd < (int)wedge_6.getSideCount(); faceOrd++){
    
    // Optional: print the vertices of the reference face 
    CellTools::printSubcellVertices(2, faceOrd, wedge_6);
    
    if( wedge_6.getKey(2, faceOrd) == shards::Triangle<3>::key ){
      CellTools::mapToReferenceSubcell(wedgeTriFacePoints, 
                                       triFaceParamCubPts,
                                       2,
                                       faceOrd,
                                       wedge_6);
      
      for(int pt = 0; pt < triFaceNumCubPts; pt++){
        std::cout << "\t Parameter point (" 
        << std::setw(10) << std::right <<  triFaceParamCubPts(pt, 0) << " , "
        << std::setw(10) << std::right <<  triFaceParamCubPts(pt, 1) << ")  "
        << std::setw(10) << "   -->    " << "(" 
        << std::setw(10) << std::right << wedgeTriFacePoints(pt, 0) << " , " 
        << std::setw(10) << std::right << wedgeTriFacePoints(pt, 1) << " , " 
        << std::setw(10) << std::right << wedgeTriFacePoints(pt, 2) << ")\n";
      }    
      std::cout << "\n";
    }
    else if(wedge_6.getKey(2, faceOrd) == shards::Quadrilateral<4>::key) {
      CellTools::mapToReferenceSubcell(wedgeQuadFacePoints, 
                                       quadFaceParamCubPts,
                                       2,
                                       faceOrd,
                                       wedge_6);
      
      for(int pt = 0; pt < quadFaceNumCubPts; pt++){
        std::cout << "\t Parameter point (" 
        << std::setw(10) << std::right <<  quadFaceParamCubPts(pt, 0) << " , "
        << std::setw(10) << std::right <<  quadFaceParamCubPts(pt, 1) << ")  "
        << std::setw(10) << "   -->    " << "(" 
        << std::setw(10) << std::right << wedgeQuadFacePoints(pt, 0) << " , " 
        << std::setw(10) << std::right << wedgeQuadFacePoints(pt, 1) << " , " 
        << std::setw(10) << std::right << wedgeQuadFacePoints(pt, 2) << ")\n";
      }    
      std::cout << "\n";
    }
    else {
      std::cout << " Invalid face encountered \n"; 
    }
  }
  
  return 0;
}
