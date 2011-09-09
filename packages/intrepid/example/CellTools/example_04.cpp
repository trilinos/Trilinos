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

/** \brief  Evaluation of a 3D vector field in physical coordinate frame
    \param  v1, v2, v3    [out] - vector mfield evaluated at the argument point
    \param  x, y, z       [in]  - argument, a point in 3D Euclidean space
*/
void vField(double& v1, double& v2, double& v3, 
            const double& x, const double& y, const double& z);

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  typedef CellTools<double>       CellTools;
  typedef RealSpaceTools<double>  RealSpaceTools;
  typedef shards::CellTopology    CellTopology;

 std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the CellTools class                        |\n" \
  << "|                                                                             |\n" \
  << "|  1) Computation of face flux, for a given vector field, on a face workset   |\n" \
  << "|  2) Computation of edge circulation, for a given vector field, on a face    |\n" \
  << "|     workset.                                                                |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n" \
  << "|                      Kara Peterson (kjpeter@sandia.gov)                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "|  EXAMPLE 1: Computation of face flux on a face workset                      |\n"\
  << "===============================================================================\n";

  
  /**  Given a vector field u(x,y,z) and a face workset we want to compute the flux of u on every
    *  face in this workset. A face workset is a set of faces that are images of the same reference
    *  face. It is defined by the following items:
    *    1. cell topology of a parent cell
    *    2. a set of nodes in physical frame defining the parenct cells in the workset
    *    3. subcell dimension and ordinal, relative to the reference cell in 1)
    *
    *  Given a face workset, the key steps to accomplish the task, , are as follows:
    *    1. Obtain cubature points on workset faces, i.e., in physical frame;
    *    2. Obtain (non-normalized) face normals at cubature points on workset faces
    *    3. Evaluate the vector field u(x,y,z) at cubature points on workset faces
    *    4. Compute dot product of u(x,y,z) and the face normals, times the cubature weights
    */
  
  /*************************************************************************************************
    *
    *  Step 0. Face workset comprising of 1 face of a Hexahedron<8> cell
    *
    ************************************************************************************************/
  
  //   Step 0.a: Specify cell topology of the parent cell
  CellTopology hexahedron_8( shards::getCellTopologyData<shards::Hexahedron<8> >() );
  
  //   Step 0.b: Specify the vertices of the parent Hexahedron<8> cell
  int worksetSize    = 2;
  int pCellNodeCount = hexahedron_8.getVertexCount();
  int pCellDim       = hexahedron_8.getDimension();
  
  FieldContainer<double> hexNodes(worksetSize, pCellNodeCount, pCellDim);
  // cell 0 bottom face vertices:
  hexNodes(0, 0, 0) = 0.00;   hexNodes(0, 0, 1) = 0.00,   hexNodes(0, 0, 2) = 0.00;          
  hexNodes(0, 1, 0) = 1.00;   hexNodes(0, 1, 1) = 0.00,   hexNodes(0, 1, 2) = 0.00;
  hexNodes(0, 2, 0) = 1.00;   hexNodes(0, 2, 1) = 1.00,   hexNodes(0, 2, 2) = 0.00;
  hexNodes(0, 3, 0) = 0.00;   hexNodes(0, 3, 1) = 1.00,   hexNodes(0, 3, 2) = 0.00;
  // cell 0 top face vertices
  hexNodes(0, 4, 0) = 0.00;   hexNodes(0, 4, 1) = 0.00,   hexNodes(0, 4, 2) = 1.00;          
  hexNodes(0, 5, 0) = 1.00;   hexNodes(0, 5, 1) = 0.00,   hexNodes(0, 5, 2) = 1.00;
  hexNodes(0, 6, 0) = 1.00;   hexNodes(0, 6, 1) = 1.00,   hexNodes(0, 6, 2) = 1.00;
  hexNodes(0, 7, 0) = 0.00;   hexNodes(0, 7, 1) = 1.00,   hexNodes(0, 7, 2) = 1.00;
  
  // cell 1 bottom face vertices:
  hexNodes(1, 0, 0) = 0.00;   hexNodes(1, 0, 1) = 0.00,   hexNodes(1, 0, 2) = 0.00;          
  hexNodes(1, 1, 0) = 1.00;   hexNodes(1, 1, 1) = 0.00,   hexNodes(1, 1, 2) = 0.00;
  hexNodes(1, 2, 0) = 1.00;   hexNodes(1, 2, 1) = 1.00,   hexNodes(1, 2, 2) = 0.00;
  hexNodes(1, 3, 0) = 0.00;   hexNodes(1, 3, 1) = 1.00,   hexNodes(1, 3, 2) = 0.00;
  // cell 1 top face vertices
  hexNodes(1, 4, 0) = 0.00;   hexNodes(1, 4, 1) = 0.00,   hexNodes(1, 4, 2) = 1.00;          
  hexNodes(1, 5, 0) = 1.00;   hexNodes(1, 5, 1) = 0.00,   hexNodes(1, 5, 2) = 1.00;
  hexNodes(1, 6, 0) = 1.00;   hexNodes(1, 6, 1) = 1.00,   hexNodes(1, 6, 2) = 0.75;
  hexNodes(1, 7, 0) = 0.00;   hexNodes(1, 7, 1) = 1.00,   hexNodes(1, 7, 2) = 1.00;
  
  
  //   Step 0.c: Specify the face ordinal, relative to the reference cell, of the face workset
  int subcellDim = 2;
  int subcellOrd = 1;  
  
  
  
  /*************************************************************************************************
    *
    *  Step 1:    Obtain Gauss points on workset faces for Hexahedron<8> topology
    *       1.1   Define cubature factory, face parametrization domain and arrays for cubature points
    *       1.2   Select Gauss rule on D = [-1,1]x[-1,1] 
    *       1.3   Map Gauss points from D to reference face and workset faces
    *
    ************************************************************************************************/
  
  //   Step 1.1.a: Define CubatureFactory
  DefaultCubatureFactory<double>  cubFactory;   
  
  //   Step 1.1.b: Define topology of the face parametrization domain as [-1,1]x[-1,1]
  CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  //   Step 1.1.c: Define storage for cubature points and weights on [-1,1]x[-1,1]
  FieldContainer<double> paramGaussWeights;
  FieldContainer<double> paramGaussPoints;
  
  //   Step 1.1.d: Define storage for cubature points on a reference face
  FieldContainer<double> refGaussPoints;
  
  //   Step 1.1.f: Define storage for cubature points on workset faces
  FieldContainer<double> worksetGaussPoints;

  //----------------
  
  //   Step 1.2.a: selects Gauss rule of order 3 on [-1,1]x[-1,1]
  Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactory.create(paramQuadFace, 3); 
  
  //   Step 1.2.b allocate storage for cubature points on [-1,1]x[-1,1]
  int cubDim       = hexFaceCubature -> getDimension();
  int numCubPoints = hexFaceCubature -> getNumPoints();
  
  // Arrays must be properly sized for the specified set of Gauss points
  paramGaussPoints.resize(numCubPoints, cubDim);
  paramGaussWeights.resize(numCubPoints);
  hexFaceCubature -> getCubature(paramGaussPoints, paramGaussWeights);
  
  //----------------
  
  //   Step 1.3.a: Allocate storage for Gauss points on the reference face
  refGaussPoints.resize(numCubPoints, pCellDim);

  //   Step 1.3.b: Allocate storage for Gauss points on the face in the workset
  worksetGaussPoints.resize(worksetSize, numCubPoints, pCellDim);

  //   Step 1.3.c: Map Gauss points to reference face: paramGaussPoints -> refGaussPoints
  CellTools::mapToReferenceSubcell(refGaussPoints,
                                   paramGaussPoints,
                                   subcellDim,                      
                                   subcellOrd,
                                   hexahedron_8);

  //   Step 1.3.d: Map Gauss points from ref. face to face workset: refGaussPoints -> worksetGaussPoints
  CellTools::mapToPhysicalFrame(worksetGaussPoints,
                                refGaussPoints,
                                hexNodes,
                                hexahedron_8);
  
  
  
  /*************************************************************************************************
    *
    *  Step 2.   Obtain (non-normalized) face normals at cubature points on workset faces
    *       2.1  Compute parent cell Jacobians at Gauss points on workset faces
    *       2.2  Compute face tangents on workset faces and their vector product
    *
    ************************************************************************************************/
  
  //   Step 2.1.a: Define and allocate storage for workset Jacobians
  FieldContainer<double> worksetJacobians(worksetSize, numCubPoints, pCellDim, pCellDim);
  
  //   Step 2.1.b: Compute Jacobians at Gauss pts. on reference face for all parent cells:
  CellTools::setJacobian(worksetJacobians,
                         refGaussPoints,
                         hexNodes,
                         hexahedron_8);
  
  //----------------
  
  //   Step 2.2.a: Allocate storage for face tangents and face normals
  FieldContainer<double> worksetFaceTu(worksetSize, numCubPoints, pCellDim);
  FieldContainer<double> worksetFaceTv(worksetSize, numCubPoints, pCellDim);
  FieldContainer<double> worksetFaceN(worksetSize, numCubPoints, pCellDim);
  
  //   Step 2.2.b: Compute face tangents
  CellTools::getPhysicalFaceTangents(worksetFaceTu,
                                     worksetFaceTv,
                                     worksetJacobians,
                                     subcellOrd,
                                     hexahedron_8);
  
  //   Step 2.2.c: Face outer normals (relative to parent cell) are uTan x vTan:
  RealSpaceTools::vecprod(worksetFaceN, worksetFaceTu, worksetFaceTv);
  
  
  
  /*************************************************************************************************
    *
    *  Step 3.   Evaluate the vector field u(x,y,z) at cubature points on workset faces
    *
    ************************************************************************************************/
  
  //   Step 3.a:  Allocate storage for vector field values at Gauss points on workset faces
  FieldContainer<double> worksetVFieldVals(worksetSize, numCubPoints, pCellDim);
  
  //   Step 3.b:  Compute vector field at Gauss points: here we take u(x,y,z) = (x,y,z)
  for(int pCellOrd = 0; pCellOrd < worksetSize; pCellOrd++){
    for(int ptOrd = 0; ptOrd < numCubPoints; ptOrd++){
      
      double x = worksetGaussPoints(pCellOrd, ptOrd, 0);
      double y = worksetGaussPoints(pCellOrd, ptOrd, 1);
      double z = worksetGaussPoints(pCellOrd, ptOrd, 2);

      vField(worksetVFieldVals(pCellOrd, ptOrd, 0), 
             worksetVFieldVals(pCellOrd, ptOrd, 1), 
             worksetVFieldVals(pCellOrd, ptOrd, 2),  x, y, z);
      
    }// pt
  }//cell
  

  /*************************************************************************************************
    *
    *  Step 4.   Compute dot product of u(x,y,z) and the face normals, times the cubature weights
    *
    ************************************************************************************************/
  
  // Allocate storage for dot product of vector field and face normals at Gauss points
  FieldContainer<double> worksetFieldDotNormal(worksetSize, numCubPoints);
  
  // Compute the dot product
  RealSpaceTools::dot(worksetFieldDotNormal, worksetVFieldVals, worksetFaceN);
  
  // Allocate storage for face fluxes on the workset
  FieldContainer<double> worksetFluxes(worksetSize);
  
  //----------------
  
  // Integration loop (temporary)
  for(int pCellOrd = 0; pCellOrd < worksetSize; pCellOrd++){
    worksetFluxes(pCellOrd) = 0.0;
    for(int pt = 0; pt < numCubPoints; pt++ ){
      worksetFluxes(pCellOrd) += worksetFieldDotNormal(pCellOrd, pt)*paramGaussWeights(pt);
    }// pt
  }//cell
  
  std::cout << "Face fluxes on workset faces : \n\n";
  for(int pCellOrd = 0; pCellOrd < worksetSize; pCellOrd++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCellOrd, subcellDim, subcellOrd);
    std::cout << " Flux = " << worksetFluxes(pCellOrd) << "\n\n";
    
  }
  
  
  
  /*************************************************************************************************
    *
    *  Optional: print Gauss points and face normals at Gauss points
    *
    ************************************************************************************************/

  // Print Gauss points on [-1,1]x[-1,1] and their images on workset faces
  std::cout \
    << "===============================================================================\n" \
    << "|                        Gauss points on workset faces:                       |\n" \
    << "===============================================================================\n";
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 2D Gauss point (" 
      << std::setw(8) << std::right <<  paramGaussPoints(pt, 0) << ", "
      << std::setw(8) << std::right <<  paramGaussPoints(pt, 1) << ")  " 
      << std::setw(8) << " -->  " << "(" 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n\n";      
  }//pCell
  
  
  // Print face normals at Gauss points on workset faces
  std::cout \
    << "===============================================================================\n" \
    << "|          Face normals (non-unit) at Gauss points on workset faces:          |\n" \
    << "===============================================================================\n";
  for(int pCell = 0; pCell < worksetSize; pCell++){
    
    CellTools::printWorksetSubcell(hexNodes, hexahedron_8, pCell, subcellDim, subcellOrd);
    
    for(int pt = 0; pt < numCubPoints; pt++){
      std::cout << "\t 3D Gauss point: (" 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetGaussPoints(pCell, pt, 2) << ")" 
      << std::setw(8) << " out. normal:  " << "(" 
      << std::setw(8) << std::right << worksetFaceN(pCell, pt, 0) << ", " 
      << std::setw(8) << std::right << worksetFaceN(pCell, pt, 1) << ", " 
      << std::setw(8) << std::right << worksetFaceN(pCell, pt, 2) << ")\n";
    }    
    std::cout << "\n";      
  }//pCell
  
  return 0;
}

/*************************************************************************************************
 *
 *  Definition of the vector field function
 *
 ************************************************************************************************/


void vField(double& v1, double& v2, double& v3, const double& x, const double& y, const double& z)
{
  v1 = x;
  v2 = y;
  v3 = z;
}


