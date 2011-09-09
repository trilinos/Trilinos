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
    \author Created by P. Bochev, D. Ridzal and K. Peterson, March 20, 2009
*/
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  typedef CellTools<double>       CellTools;
  typedef shards::CellTopology    CellTopology;

 std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the CellTools class                        |\n" \
  << "|                                                                             |\n" \
  << "|  1) Definition of integration points on edge and face worksets              |\n" \
  << "|  2) Computation of face normals and edge tangents on face and edge worksets |\n" \
  << "|  2) Computation of side normals on face and edge worksets                   |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n" \
  << "|                      Kara Peterson (kjpeter@sandia.gov)                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "|  EXAMPLE 1: Definition of integration points on a Triangle edge workset     |\n"\
  << "===============================================================================\n";
 
 /*
  *  1. Common tasks for getting integration points on edge worksets
  */
  
  // Step 1.a: Define CubatureFactory
  DefaultCubatureFactory<double>  cubFactory;   
  
  // Step 1.b: Define topology of the edge parametrization domain [-1,1]
  CellTopology paramEdge(shards::getCellTopologyData<shards::Line<2> >() );
  
  // Step 1.c: Define storage for cubature points and weights on [-1,1]
  FieldContainer<double> paramEdgeWeights;
  FieldContainer<double> paramEdgePoints;
  
  // Step 1.d: Define storage for cubature points on a reference edge
  FieldContainer<double> refEdgePoints;
  
  // Step 1.f: Define storage for cubature points on workset edges
  FieldContainer<double> worksetEdgePoints;


  
  /*
   *  2. Define edge workset comprising of 4 edges corresponding to reference edge 2 on Triangle<3>
   */
  
  // Step 2.a: Specify cell topology of the parent cell
  CellTopology triangle_3( shards::getCellTopologyData<shards::Triangle<3> >() );
  
  // Step 2.b: Specify the vertices of 4 parent Triangle<3> cells
  int worksetSize    = 4;
  int pCellNodeCount = triangle_3.getVertexCount();
  int pCellDim       = triangle_3.getDimension();
  
  FieldContainer<double> triNodes(worksetSize, pCellNodeCount, pCellDim);
  // Triangle 0
  triNodes(0, 0, 0) = 0.5;   triNodes(0, 0, 1) = 2.0;
  triNodes(0, 1, 0) = 0.0;   triNodes(0, 1, 1) = 1.0;
  triNodes(0, 2, 0) = 1.0;   triNodes(0, 2, 1) = 1.0;
  // Triangle 1
  triNodes(1, 0, 0) = 1.0;   triNodes(1, 0, 1) = 1.0;               
  triNodes(1, 1, 0) = 1.0;   triNodes(1, 1, 1) = 0.0;
  triNodes(1, 2, 0) = 0.0;   triNodes(1, 2, 1) = 0.0;
  // Triangle 2
  triNodes(2, 0, 0) = 0.0;   triNodes(2, 0, 1) = 0.0;
  triNodes(2, 1, 0) = 1.0;   triNodes(2, 1, 1) = 1.0;
  triNodes(2, 2, 0) = 0.0;   triNodes(2, 2, 1) = 1.0;
  // Triangle 3
  triNodes(3, 0, 0) = 1.0;   triNodes(3, 0, 1) = 1.0;
  triNodes(3, 1, 0) = 2.5;   triNodes(3, 1, 1) = 1.5;
  triNodes(3, 2, 0) = 0.5;   triNodes(3, 2, 1) = 2.0;
  
  // Step 2.c: Specify dimension and ordinal (relative to ref. cell) of the subcells in the workset
  int subcellDim = 1;
  int subcellOrd = 2;  
  
  
  
  /*
   *  3. Obtain the desired Gauss rule on the edge parametrization domain [-1,1]
   */
  
  // Step 3.a: selects Gauss rule of order 6 on [-1,1]
  Teuchos::RCP<Cubature<double> > triEdgeCubature = cubFactory.create(paramEdge, 6); 
  
  // Step 3.b allocate storage for cubature points on [-1,1]
  int cubDim       = triEdgeCubature -> getDimension();
  int numCubPoints = triEdgeCubature -> getNumPoints();
  
  // Arrays must be properly sized for the specified set of Gauss points
  paramEdgePoints.resize(numCubPoints, cubDim);
  paramEdgeWeights.resize(numCubPoints);
  triEdgeCubature -> getCubature(paramEdgePoints, paramEdgeWeights);

  
  
  /*
   *  4. Map Gauss points from [-1,1] to every edge in the workset
   */

  // Step 4.a: Allocate storage for integration points on the reference edge
  refEdgePoints.resize(numCubPoints, pCellDim);
 
  // Step 4.b: Allocate storage for integration points on the edge workset
  worksetEdgePoints.resize(worksetSize, numCubPoints, pCellDim);
  
  // Step 4.c: Map Gauss points to reference edge: paramEdgePoints -> refEdgePoints
  CellTools::mapToReferenceSubcell(refEdgePoints,
                                   paramEdgePoints,
                                   subcellDim,                      
                                   subcellOrd,
                                   triangle_3);
  
  // Step 4.d: Map Gauss points from ref. edge to edge workset: refEdgePoints -> worksetEdgePoints
  CellTools::mapToPhysicalFrame(worksetEdgePoints,
                                refEdgePoints,
                                triNodes,
                                triangle_3);
  
  
  
  /*
   *  5. Print Gauss points on the edges of the workset
   */
  for(int pCell = 0; pCell < worksetSize; pCell++){
      
    CellTools::printWorksetSubcell(triNodes, triangle_3, pCell, subcellDim, subcellOrd);
     
      for(int pt = 0; pt < numCubPoints; pt++){
        std::cout << "\t 1D Gauss point " 
        << std::setw(12) << std::right << paramEdgePoints(pt, 0) << std::setw(10) << "  -->  " << "(" 
        << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
        << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 1) << ")\n";
      }    
      std::cout << "\n";      
  }//pCell
  
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "|  EXAMPLE 2: Definition of integration points on a Quadrilateral edge workset|\n"\
    << "===============================================================================\n";
  
  /*
   *  2. Define edge workset comprising of 2 edges corresponding to reference edge 2 on Quadrilateral<4>
   *     (Can reuse Step 1, Example 1 because edge parametrization domain is always [-1,1])
   */
  
  // Step 2.a: Specify cell topology of the parent cell
  CellTopology quadrilateral_4( shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  // Step 2.b: Specify the vertices of 2 parent Quadrilateral<4> cells
  worksetSize    = 2;
  pCellNodeCount = quadrilateral_4.getVertexCount();
  pCellDim       = quadrilateral_4.getDimension();
  
  FieldContainer<double> quadNodes(worksetSize, pCellNodeCount, pCellDim);
  // Quadrilateral 0
  quadNodes(0, 0, 0) = 1.00;   quadNodes(0, 0, 1) = 1.00;               
  quadNodes(0, 1, 0) = 2.00;   quadNodes(0, 1, 1) = 0.75;
  quadNodes(0, 2, 0) = 1.75;   quadNodes(0, 2, 1) = 2.00;  
  quadNodes(0, 3, 0) = 1.25;   quadNodes(0, 3, 1) = 2.00; 
  // Quadrilateral 1
  quadNodes(1, 0, 0) = 2.00;   quadNodes(1, 0, 1) = 0.75;               
  quadNodes(1, 1, 0) = 3.00;   quadNodes(1, 1, 1) = 1.25;
  quadNodes(1, 2, 0) = 2.75;   quadNodes(1, 2, 1) = 2.25;
  quadNodes(1, 3, 0) = 1.75;   quadNodes(1, 3, 1) = 2.00;
  
  // Step 2.c: Specify dimension and ordinal (relative to ref. cell) of the subcells in the workset
  subcellDim = 1;
  subcellOrd = 2;  
  
  /*
   *  3. Obtain the desired Gauss rule on the edge parametrization domain [-1,1]
   */
  
  // Step 3.a: selects Gauss rule of order 4 on [-1,1]
  Teuchos::RCP<Cubature<double> > quadEdgeCubature = cubFactory.create(paramEdge, 6); 
  
  // Step 3.b: allocate storage for cubature points 
  cubDim       = quadEdgeCubature -> getDimension();
  numCubPoints = quadEdgeCubature -> getNumPoints();
  
  // Arrays must be properly sized for the specified set of Gauss points
  paramEdgePoints.resize(numCubPoints, cubDim);
  paramEdgeWeights.resize(numCubPoints);
  quadEdgeCubature -> getCubature(paramEdgePoints, paramEdgeWeights);
  
  /*
   *  4. Map Gauss points from [-1,1] to every edge in the workset
   */
  
  // Step 4.a: Allocate storage for integration points on the reference edge
  refEdgePoints.resize(numCubPoints, pCellDim);
  
  // Step 4.b: Allocate storage for integration points on the edge workset
  worksetEdgePoints.resize(worksetSize, numCubPoints, pCellDim);
  
  // Step 4.c: Map Gauss points to reference edge: paramEdgePoints -> refEdgePoints
  CellTools::mapToReferenceSubcell(refEdgePoints,
                                   paramEdgePoints,
                                   subcellDim,                      
                                   subcellOrd,
                                   quadrilateral_4);
  
  // Step 4.d: Map Gauss points from ref. edge to edge workset: refEdgePoints -> worksetEdgePoints
  CellTools::mapToPhysicalFrame(worksetEdgePoints,
                                refEdgePoints,
                                quadNodes,
                                quadrilateral_4);
  
  /*
   *  5. Print Gauss points on the edges of the workset
   */
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(quadNodes, quadrilateral_4, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 1D Gauss point " 
      << std::setw(12) << std::right << paramEdgePoints(pt, 0) << std::setw(10) << "  -->  " << "(" 
      << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
      << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 1) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 3: Edge tangents at Gauss points on a Quadrilateral edge workset    |\n"\
    << "===============================================================================\n";
  
  /*  This task requires Gauss points on edge parametrization domain [-1,1] and on the 
   *  reference edge whose ordinal matches the edge workset ordinal. This repeats the first few
   *  steps from Example 2:
   *  
   *  1. Define cubature factory and topology for edge parametrization domain [-1,1];
   *     (Can reuse Step 1, Example 1 because edge parametrization domain is always [-1,1])
   *  2. Define an edge workset;
   *  3. Obtain the desired Gauss rule on the edge parametrization domain [-1,1];
   *  4. Map Gauss points from [-1,1] to reference edge 
   *     NOTE: this example only demonstrates computation of the edge tangents and so, Gauss 
   *     points on the edge workset are not needed. Thus we skip mapping Gauss points from 
   *     reference edge to edge workset
   *
   *  5. Compute Jacobians at Gauss points on reference edge for all cells in the workset:
   */
  
  // Step 5.a: Define and allocate storage for workset Jacobians
  FieldContainer<double> worksetJacobians(worksetSize, numCubPoints, pCellDim, pCellDim);
  
  // Step 5.b: Compute Jacobians at Gauss pts. on reference edge for all parent cells:
  CellTools::setJacobian(worksetJacobians,
                         refEdgePoints,
                         quadNodes,
                         quadrilateral_4);
  /*
   * 6. Get the (non-normalized) edge tangents for the edge workset:
   */
  // Step 6.a: Allocate storage for edge tangents
  FieldContainer<double> edgeTangents(worksetSize, numCubPoints, pCellDim);
  
  // Step 6.b: Compute the edge tangents:
  CellTools::getPhysicalEdgeTangents(edgeTangents,
                                     worksetJacobians,
                                     subcellOrd,
                                     quadrilateral_4); 
  
  // Step 6.c: Print edge tangents at Gauss points on workset edges (these Gauss points were computed in Example 2)
  std::cout 
    << "Edge tangents computed by CellTools::getPhysicalEdgeTangents.\n"
    << "Local edge ordinal = " << subcellOrd <<"\n";
  
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(quadNodes, quadrilateral_4, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 2D Gauss point: (" 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 1) << ")  " 
      << std::setw(8) << " edge tangent:  " << "(" 
      << std::setw(8) << std::right << edgeTangents(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << edgeTangents(pCell, pt, 1) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell

  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 4: Side normals at Gauss points on a Quadrilateral side workset     |\n"\
    << "===============================================================================\n";

  /* For this task we reuse the edge workset from Example 3 as a side workset and compute
   * normals to the sides in that set. This task repeats the first 5 steps from Example 3
   * The only difference is that in Step 6 we call CellTools::getPhysicalSideNormals
   */
  /*
   * 6. Get the (non-normalized) side normals for the side (edge) workset:
   */
  // Step 6.a: Allocate storage for side normals
  FieldContainer<double> sideNormals(worksetSize, numCubPoints, pCellDim);
  
  // Step 6.b: Compute the side normals:
  CellTools::getPhysicalSideNormals(sideNormals,
                                    worksetJacobians,
                                    subcellOrd,
                                    quadrilateral_4); 
  
  // Step 6.c: Print side normals at Gauss points on workset sides (these Gauss points were computed in Example 2)
  std::cout 
    << "Side normals computed by CellTools::getPhysicalSideNormals.\n"
    << "Local side (edge) ordinal = " << subcellOrd <<"\n";
  
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(quadNodes, quadrilateral_4, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 2D Gauss point: (" 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 1) << ")  " 
      << std::setw(8) << " side normal:  " << "(" 
      << std::setw(8) << std::right << sideNormals(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << sideNormals(pCell, pt, 1) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
    
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 5: Definition of integration points on a Hexahedron edge workset    |\n"\
    << "===============================================================================\n";
  
  /*
   *  2. Define edge workset comprising of 1 edge corresponding to reference edge 10 on Hexahedron<8>
   *     (Can reuse Step 1, Example 1 because edge parametrization domain is always [-1,1])
   */
  
  // Step 2.a: Specify cell topology of the parent cell
  CellTopology hexahedron_8( shards::getCellTopologyData<shards::Hexahedron<8> >() );
  
  // Step 2.b: Specify the vertices of the parent Hexahedron<8> cell
  worksetSize    = 1;
  pCellNodeCount = hexahedron_8.getVertexCount();
  pCellDim       = hexahedron_8.getDimension();

  FieldContainer<double> hexNodes(worksetSize, pCellNodeCount, pCellDim);
  // bottom face vertices
  hexNodes(0, 0, 0) = 0.00;   hexNodes(0, 0, 1) = 0.00,   hexNodes(0, 0, 2) = 0.00;          
  hexNodes(0, 1, 0) = 1.00;   hexNodes(0, 1, 1) = 0.00,   hexNodes(0, 1, 2) = 0.00;
  hexNodes(0, 2, 0) = 1.00;   hexNodes(0, 2, 1) = 1.00,   hexNodes(0, 2, 2) = 0.00;
  hexNodes(0, 3, 0) = 0.00;   hexNodes(0, 3, 1) = 1.00,   hexNodes(0, 3, 2) = 0.00;
  // top face vertices
  hexNodes(0, 4, 0) = 0.00;   hexNodes(0, 4, 1) = 0.00,   hexNodes(0, 4, 2) = 1.00;          
  hexNodes(0, 5, 0) = 1.00;   hexNodes(0, 5, 1) = 0.00,   hexNodes(0, 5, 2) = 1.00;
  hexNodes(0, 6, 0) = 1.00;   hexNodes(0, 6, 1) = 1.00,   hexNodes(0, 6, 2) = 0.75;
  hexNodes(0, 7, 0) = 0.00;   hexNodes(0, 7, 1) = 1.00,   hexNodes(0, 7, 2) = 1.00;
  
  // An alternative hex obtained by intersection of the unit cube [0,1]^3 and the plane
  // z = 1 - 1/4x - 1/4y. The top face (local ordinal 5) of the resulting hex lies in this plane
  // and has the same normal vector parallel to (1/4,1/4,1). This workset allows to test face normals
  FieldContainer<double> hexNodesAlt(worksetSize, pCellNodeCount, pCellDim);
  // bottom face vertices
  hexNodesAlt(0, 0, 0) = 0.00;   hexNodesAlt(0, 0, 1) = 0.00,   hexNodesAlt(0, 0, 2) = 0.00;          
  hexNodesAlt(0, 1, 0) = 1.00;   hexNodesAlt(0, 1, 1) = 0.00,   hexNodesAlt(0, 1, 2) = 0.00;
  hexNodesAlt(0, 2, 0) = 1.00;   hexNodesAlt(0, 2, 1) = 1.00,   hexNodesAlt(0, 2, 2) = 0.00;
  hexNodesAlt(0, 3, 0) = 0.00;   hexNodesAlt(0, 3, 1) = 1.00,   hexNodesAlt(0, 3, 2) = 0.00;
  // top face vertices
  hexNodesAlt(0, 4, 0) = 0.00;   hexNodesAlt(0, 4, 1) = 0.00,   hexNodesAlt(0, 4, 2) = 1.00;          
  hexNodesAlt(0, 5, 0) = 1.00;   hexNodesAlt(0, 5, 1) = 0.00,   hexNodesAlt(0, 5, 2) = 0.75;
  hexNodesAlt(0, 6, 0) = 1.00;   hexNodesAlt(0, 6, 1) = 1.00,   hexNodesAlt(0, 6, 2) = 0.50;
  hexNodesAlt(0, 7, 0) = 0.00;   hexNodesAlt(0, 7, 1) = 1.00,   hexNodesAlt(0, 7, 2) = 0.75;  
  
  // Step 2.c: Specify the edge ordinal, relative to the reference cell, of the edge workset
  subcellDim = 1;
  subcellOrd = 5;  
  
  /*
   *  3. Obtain the desired Gauss rule on the edge parametrization domain [-1,1]
   */
  
  // Step 3.a: selects Gauss rule of order 4 on [-1,1]
  Teuchos::RCP<Cubature<double> > hexEdgeCubature = cubFactory.create(paramEdge, 4); 
  
  // Step 3.b allocate storage for cubature points
  cubDim       = hexEdgeCubature -> getDimension();
  numCubPoints = hexEdgeCubature -> getNumPoints();
  
  // Arrays must be properly sized for the specified set of Gauss points
  paramEdgePoints.resize(numCubPoints, cubDim);
  paramEdgeWeights.resize(numCubPoints);
  hexEdgeCubature -> getCubature(paramEdgePoints, paramEdgeWeights);
  
  /*
   *  4. Map Gauss points from [-1,1] to every edge in the workset
   */
  
  // Step 4.a: Allocate storage for integration points on the reference edge
  refEdgePoints.resize(numCubPoints, pCellDim);
  
  // Step 4.b: Allocate storage for integration points on the edge workset
  worksetEdgePoints.resize(worksetSize, numCubPoints, pCellDim);
  
  // Step 4.c: Map Gauss points to reference edge: paramEdgePoints -> refEdgePoints
  CellTools::mapToReferenceSubcell(refEdgePoints,
                                   paramEdgePoints,
                                   subcellDim,                      
                                   subcellOrd,
                                   hexahedron_8);
  
  // Step 4.d: Map Gauss points from ref. edge to edge workset: refEdgePoints -> worksetEdgePoints
  CellTools::mapToPhysicalFrame(worksetEdgePoints,
                                refEdgePoints,
                                hexNodes,
                                hexahedron_8);
  
  /*
   *  5. Print Gauss points on the edges of the workset
   */
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 1D Gauss point " 
      << std::setw(12) << std::right << paramEdgePoints(pt, 0) << std::setw(10) << "  -->  " << "(" 
      << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
      << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 1) << ", " 
      << std::setw(10) << std::right << worksetEdgePoints(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 6: Edge tangents at Gauss points on a Hexahedron edge workset       |\n"\
    << "===============================================================================\n";
  
  /*  This task requires Gauss points on edge parametrization domain [-1,1] and on the 
   *  reference edge whose ordinal matches the edge workset ordinal. This repeats the first few
   *  steps from Example 5:
   *  
   *  1. Define cubature factory and topology for edge parametrization domain [-1,1];
   *     (Can reuse Step 1, Example 1 because edge parametrization domain is always [-1,1])
   *  2. Define an edge workset;
   *  3. Obtain the desired Gauss rule on the edge parametrization domain [-1,1];
   *  4. Map Gauss points from [-1,1] to reference edge 
   *     NOTE: this example only demonstrates computation of the edge tangents and so, Gauss 
   *     points on the edge workset are not needed. Thus we skip mapping Gauss points from 
   *     reference edge to edge workset
   *
   *  5. Compute Jacobians at Gauss points on reference edge for all cells in the workset:
   */
  
  // Step 5.a: Define and allocate storage for workset Jacobians
  worksetJacobians.resize(worksetSize, numCubPoints, pCellDim, pCellDim);
  
  // Step 5.b: Compute Jacobians at Gauss pts. on reference edge for all parent cells:
  CellTools::setJacobian(worksetJacobians,
                         refEdgePoints,
                         hexNodes,
                         hexahedron_8);
  
  /*
   * 6. Get the (non-normalized) edge tangents for the edge workset:
   */
  // Step 6.a: Allocate storage for edge tangents
  edgeTangents.resize(worksetSize, numCubPoints, pCellDim);
  
  // Step 6.b: Compute the edge tangents:
  CellTools::getPhysicalEdgeTangents(edgeTangents,
                                     worksetJacobians,
                                     subcellOrd,
                                     hexahedron_8); 
  
  // Step 6.c: Print edge tangents at Gauss points on workset edges (these Gauss points were computed in Example 5)
  std::cout 
    << "Edge tangents computed by CellTools::getPhysicalEdgeTangents.\n"
    << "Local edge ordinal = " << subcellOrd <<"\n";
  
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 3D Gauss point: (" 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetEdgePoints(pCell, pt, 2) << ")  " 
      << std::setw(8) << " edge tangent:  " << "(" 
      << std::setw(8) << std::right << edgeTangents(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << edgeTangents(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << edgeTangents(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
   
 std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 7: Definition of integration points on a Hexahedron face workset    |\n"\
    << "===============================================================================\n";
  /*
   *  1. Common tasks for getting integration points on face worksets: 
   *  1.a: can reuse the cubature factory, but need parametrization domain for the faces
   */
  
  // Step 1.b: Define topology of the face parametrization domain as [-1,1]x[-1,1]
  CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  // Step 1.c: Define storage for cubature points and weights on [-1,1]x[-1,1]
  FieldContainer<double> paramFaceWeights;
  FieldContainer<double> paramFacePoints;
  
  // Step 1.d: Define storage for cubature points on a reference face
  FieldContainer<double> refFacePoints;
  
  // Step 1.f: Define storage for cubature points on workset faces
  FieldContainer<double> worksetFacePoints;
    
  /*
   *  2. Define face workset comprising of 1 face corresponding to reference face 5 on Hexahedron<8>
   *  2.a: Reuse the parent cell topology from Example 3
   *  2.b: Reuse the vertices from Example 3
   */
  
  // Step 2.c: Specify dimension and ordinal (relative to ref. cell) of the subcells in the workset
  subcellDim = 2;
  subcellOrd = 5;  
  
  /*
   *  3. Obtain the desired Gauss rule on the face parametrization domain [-1,1]x[-1,1]
   */
  
  // Step 3.a: selects Gauss rule of order 3 on [-1,1]x[-1,1]
  Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactory.create(paramQuadFace, 3); 
  
  // Step 3.b allocate storage for cubature points on [-1,1]x[-1,1]
  cubDim       = hexFaceCubature -> getDimension();
  numCubPoints = hexFaceCubature -> getNumPoints();
  
  // Arrays must be properly sized for the specified set of Gauss points
  paramFacePoints.resize(numCubPoints, cubDim);
  paramFaceWeights.resize(numCubPoints);
  hexFaceCubature -> getCubature(paramFacePoints, paramFaceWeights);
  
  /*
   *  4. Map Gauss points from [-1,1]x[-1,1] to every face in the workset
   */
    
  // Step 4.a: Allocate storage for integration points on the reference face
  refFacePoints.resize(numCubPoints, pCellDim);

  // Step 4.b: Allocate storage for integration points on the face in the workset
  worksetFacePoints.resize(worksetSize, numCubPoints, pCellDim);

  // Step 4.c: Map Gauss points to reference face: paramFacePoints -> refFacePoints
  CellTools::mapToReferenceSubcell(refFacePoints,
                                   paramFacePoints,
                                   subcellDim,                      
                                   subcellOrd,
                                   hexahedron_8);

  // Step 4.d: Map Gauss points from ref. face to face workset: refFacePoints -> worksetFacePoints
  CellTools::mapToPhysicalFrame(worksetFacePoints,
                                refFacePoints,
                                hexNodesAlt,
                                hexahedron_8);
  
  
  /*
   *  5. Print Gauss points on the faces of the workset
   */
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodesAlt, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 2D Gauss point (" 
      << std::setw(8) << std::right <<  paramFacePoints(pt, 0) << ", "
      << std::setw(8) << std::right <<  paramFacePoints(pt, 1) << ")  " 
      << std::setw(8) << " -->  " << "(" 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 8: Face normals at Gauss points on a Hexahedron face workset        |\n"\
    << "===============================================================================\n";
  
  /*  This task requires Gauss points on face parametrization domain [-1,1]x[-1,1] and on the 
   *  reference face whose ordinal matches the face workset ordinal. This repeats the first few
   *  steps from Example 5:
   *  
   *  1. Define cubature factory and topology for face parametrization domain [-1,1]x[-1,1];
   *  2. Define a face workset;
   *  3. Obtain the desired Gauss rule on the face parametrization domain [-1,1]x[-1,1];
   *  4. Map Gauss points from [-1,1]x[-1,1] to reference face
   *     NOTE: this example only demonstrates computation of the face normals and so, Gaus 
   *     points on the face workset are not needed. Thus we skip mapping Gauss points from 
   *     reference face to face workset
   *
   *  5. Compute Jacobians at Gauss points on reference face for all cells in the workset:
   */
  
  // Step 5.a: Define and allocate storage for workset Jacobians (reuse FC from example 5)
  worksetJacobians.resize(worksetSize, numCubPoints, pCellDim, pCellDim);
  
  // Step 5.b: Compute Jacobians at Gauss pts. on reference face for all parent cells:
  CellTools::setJacobian(worksetJacobians,
                         refFacePoints,
                         hexNodesAlt,
                         hexahedron_8);
  
  
  /*
   * 6. Get the (non-normalized) face normals for the face workset directly
   */
  // Step 6.a: Allocate storage for face normals
  FieldContainer<double> faceNormals(worksetSize, numCubPoints, pCellDim);
  
  // Step 6.b: Compute the face normals
  CellTools::getPhysicalFaceNormals(faceNormals,
                                    worksetJacobians,
                                    subcellOrd,
                                    hexahedron_8);
  
  // Step 6.c: Print face normals at Gauss points on workset faces (these Gauss points were computed in Example 5)
  std::cout 
    << "Face normals computed by CellTools::getPhysicalFaceNormals\n"
    << "Local face ordinal = " << subcellOrd <<"\n";

  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodesAlt, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 3D Gauss point: (" 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 2) << ")  " 
      << std::setw(8) << " outer normal:  " << "(" 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
    
  /*
   * 7. Get the (non-normalized) face normals for the face workset using the face tangents. This may
   *    be useful if, for whatever reason,  face tangents are needed independently 
   */
  // Step 7.a: Allocate storage for face tangents
  FieldContainer<double> uFaceTan(worksetSize, numCubPoints, pCellDim);
  FieldContainer<double> vFaceTan(worksetSize, numCubPoints, pCellDim);
  
  // Step 7.b: Compute face tangents
  CellTools::getPhysicalFaceTangents(uFaceTan,
                                     vFaceTan,
                                     worksetJacobians,
                                     subcellOrd,
                                     hexahedron_8);
  
  // Step 7.c: Face outer normals (relative to parent cell) are uTan x vTan:
  RealSpaceTools<double>::vecprod(faceNormals, uFaceTan, vFaceTan);
  
  // Step 7.d: Print face normals at Gauss points on workset faces (these points were computed in Example 7)
  std::cout << "Face normals computed by CellTools::getPhysicalFaceTangents followed by vecprod\n";
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 3D Gauss point: (" 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 2) << ")  " 
      << std::setw(8) << " outer normal:  " << "(" 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << faceNormals(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 8: Side normals at Gauss points on a Hexahedron side workset        |\n"\
    << "===============================================================================\n";
  
  /* For this task we reuse the edge workset from Example 7 as a side workset and compute
   * normals to the sides in that set. This task repeats the first 5 steps from Example 3
   * The only difference is that in Step 6 we call CellTools::getPhysicalSideNormals
   */
  
  /*
   * 6. Get the (non-normalized) side normals for the side (face) workset:
   */
  // Step 6.a: Allocate storage for side normals
  sideNormals.resize(worksetSize, numCubPoints, pCellDim);
  
  // Step 6.b: Compute the side normals:
  CellTools::getPhysicalSideNormals(sideNormals,
                                    worksetJacobians,
                                    subcellOrd,
                                    hexahedron_8); 
  
  // Step 7.d: Print side normals at Gauss points on workset sides (these points were computed in Example 7)
  std::cout << "Side normals computed by CellTools::getPhysicalSideNormals\n";
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 3D Gauss point: (" 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetFacePoints(pCell, pt, 2) << ")  " 
      << std::setw(8) << " side normal:  " << "(" 
      << std::setw(8) << std::right << sideNormals(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << sideNormals(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << sideNormals(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  return 0;
}
