// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
\brief  Example of the Multicell class.
\author Created by P. Bochev and D. Ridzal, December 20, 2007
*/
#include "Intrepid_MultiCell.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Teuchos_Array.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the MultiCell class                        |\n" \
  << "|                                                                             |\n" \
  << "|  1) Working with cell types that have reference cells: setting the atlas    |\n" \
  << "|  2) Mapping points to and from the reference cell                           |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "|  EXAMPLE 1: setting the atlas of a MultiCell                                |\n"\
  << "|    We are creating  MultiCell objects whose generating cell types have      |\n"\
  << "|    reference cells. The atlas contains charts to map reference cell to      |\n"\
  << "|    each one of the physical cells that comprise the MultiCell object        |\n"\
  << "===============================================================================\n";
 
  // Define vertex data for 4 CELL_TRI cells
  double triNodes[] = {
    // triangle 0
    0.5, 2.0,
    0.0, 1.0,
    1.0, 1.0,
    // triangle 1
    1.0, 1.0,               
    1.0, 0.0,
    0.0, 0.0,
    // triangle 2
    0.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    // triangle 3
    1.0, 1.0,
    2.5, 1.5,
    0.5, 2.0
  };
  
  /*
   
   Need to fix point assignment operation shapePoints_[0][0] = Point<double> (data) 
   because shapePoints have the wrong space dimensions and operator = throws an exception
   
  // Define an array of 4 shape point sets to make some triangles curvilinear
  Teuchos::Array< ShapePoints<double> > triShapePoints(4);
  
  // The first TRI cell will be affine, so we don't provide any extra shape points:
  triShapePoints[0].shapePoints_.resize(0);
  triShapePoints[0].chartDegree_ = 1;
  
  // The second TRI cell will have chart of degree 2 which requires 3 edge points. Resize the
  // shapePoints_ for this cell to make room for edge points
  triShapePoints[1].shapePoints_.resize(1);
  triShapePoints[1].chartDegree_ = 1;
  
  //The third TRI cell will have chart of degree 3 which requires 6 edge points. . Resize the
  // shapePoints_ for this cell to make room for edge points
  triShapePoints[2].shapePoints_.resize(1);
  triShapePoints[2].chartDegree_ = 2;

  // The las TRI cell will be affine, so we don't provide any extra shape points:
  triShapePoints[0].shapePoints_.resize(0);
  triShapePoints[0].chartDegree_ = 1;
  
  
  // For 2nd TRI resize shapePoints_[0] (where edge points are stored) to 3
  triShapePoints[1].shapePoints_[0].resize(3);    
  triShapePoints[1].shapePoints_[0][0] = Point<double> (1.25, 0.5,FRAME_PHYSICAL);
  triShapePoints[1].shapePoints_[0][1] = Point<double> (0.5, -0.25,FRAME_PHYSICAL);
  triShapePoints[1].shapePoints_[0][2] = Point<double> (0.5,  0.5, FRAME_PHYSICAL);
  
  // For 3rd TRI resize shapePoints_[0] (where edge points are stored) to 6
  triShapePoints[2].shapePoints_[0].resize(6);    
  
  // 2 points for edge 0
  triShapePoints[2].shapePoints_[0][0] = Point<double> (0.3,  0.3 ,FRAME_PHYSICAL);
  triShapePoints[2].shapePoints_[0][1] = Point<double> (0.6,  0.6 ,FRAME_PHYSICAL);
  
  // 2 Points for edge 1
 // triShapePoints[2].shapePoints_[0][2] = Point<double> ;
  //triShapePoints[2].shapePoints_[0][3] = Point<double> ;
   
   */

  
  // Define vertex data for 2 CELL_QUAD cells
  double quadNodes[] = {
    // 1st QUAD
    1.00, 1.00,               
    2.00, 0.75,
    1.75, 2.00,  
    1.25, 2.00, 
    // 2ND QUAD
    2.00, 0.75,               
    3.00, 1.25,
    2.75, 2.25,
    1.75, 2.00};
  
  // Define vertex data for a single CELL_HEX cell
  double hexNodes[] = {
    // bottom face vertices
    1.00, 1.00,  0.00,          
    2.00, 0.75, -0.25,
    1.75, 2.00,  0.00,
    1.25, 2.00,  0.25,
    // top face vertices
    1.25, 0.75,  0.75,          
    1.75, 1.00,  1.00,
    2.00, 2.00,  1.25,
    1.00, 2.00,  1.00};
  
  // MultiCell holding 4 traingle cells
  MultiCell<double> triMcell(4,CELL_TRI, triNodes);
  cout << triMcell << endl;
  
  // MultiCell holding 2 quadrilateral cells
  MultiCell<double> quadMcell(2,CELL_QUAD,quadNodes);
  cout << quadMcell << endl;
  
  // Multicell holding a single hexahedron
  MultiCell<double> hexMcell(1,CELL_HEX,hexNodes);
  cout << hexMcell << endl;
  
  // Now set the atlas for each MultiCell
  cout << "Setting the atlas for each MultiCell...\n";
  triMcell.setAtlas();
  quadMcell.setAtlas();
  hexMcell.setAtlas();
  
  cout << "\tAtlas status of  triMcell is: " << triMcell.getAtlasStatusName() << "\n";
  cout << "\tAtlas status of quadMcell is: " << triMcell.getAtlasStatusName() << "\n";
  cout << "\tAtlas status of  hexMcell is: " << triMcell.getAtlasStatusName() << "\n";
  
  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: Using getJacobian to compute Jacobians at cubature point sets   |\n"\
    << "===============================================================================\n";
  
  // Test getJacobian method: define cubature factory and the necesssary point/weight arrays
  DefaultCubatureFactory<double> CFactory;
  int numPts;
  Teuchos::Array<Point<double> > cubPts;
  Teuchos::Array<double> cubWeights;
  int cellId=0;
  int subcellDim = 2;
  int subcellId  = 0;
  
  // Define cubature rule of degree 2 for a TRI cell:
  Teuchos::RCP<Cubature<double> > triCub =  CFactory.create(CELL_TRI, 2);
  triCub -> getCubature(numPts, cubPts,cubWeights);
  
  // Get Jacobians for subcDim = 2 (the tri), with cellId = 0
  triMcell.initializeMeasures(subcellDim,
                              subcellId,cubPts,
                              cubWeights);
  cout << " DF, DF^{-T} and det(DF) at all cubature points in the TRI cell with cellId = "<< cellId << "\n";
  for(int cp = 0; cp < numPts; cp++ ){
    cout << " At cubature point " << cp << " ---------------> \n\n";
    cout << "\t      DF = " << triMcell.getJacobian(cellId,subcellDim,subcellId)[cp];
    cout << "\t det(DF) = " << triMcell.getMeasure(cellId,subcellDim,subcellId)[cp]  <<"\n\n";
    cout << "\t DF^{-T} = " << triMcell.getJacobianTInv(cellId,subcellDim,subcellId)[cp]  <<"\n";
  }
    
  // Define cubature of degree 2 for a QUAD cell:
  Teuchos::RCP<Cubature<double> > quadCub =  CFactory.create(CELL_QUAD, 2);
  quadCub -> getCubature(numPts, cubPts,cubWeights);
  
  // Get Jacobians for subcDim = 2 (the QUAD), with cellId = 0
  quadMcell.initializeMeasures(subcellDim,subcellId,cubPts,cubWeights);
  cout << " Jacobians at all cubature points in the QUAD cell with cellId = "<< cellId << "\n";
  for(int cp = 0; cp < numPts; cp++ ){
    cout << " At cubature point " << cp << " ---------------> \n\n";
    cout << "\t      DF = " << quadMcell.getJacobian(cellId,subcellDim,subcellId)[cp];
    cout << "\t det(DF) = " << quadMcell.getMeasure(cellId,subcellDim,subcellId)[cp]  <<"\n\n";
    cout << "\t DF^{-T} = " << quadMcell.getJacobianTInv(cellId,subcellDim,subcellId)[cp]  <<"\n";
    cout << "\t DF*DF^{-t} = " <<
      quadMcell.getJacobian(cellId,subcellDim,subcellId)[cp]*\
      quadMcell.getJacobianTInv(cellId,subcellDim,subcellId)[cp] << "\n";
   }
  
  
  // This is how to access individual elements in the Jacobians
  int i = 0; int j = 1; cellId = 0;
  cout << "CellId = " << cellId << " Cubature Point Id = " << 0
    << "  Jacobian["<<i<<"]["<<j<<"] = "
    << quadMcell.getJacobian(cellId,subcellDim,subcellId)[0](i,j) <<"\n";


  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: Using charts with triMcell                                       |\n"\
    << "===============================================================================\n";
  
  // Using charts with the triMcell: we will map the vertices of the reference triangle
  Point<double> ref_tri_v0(0.0,0.0,FRAME_REFERENCE), tri_v0(2);
  Point<double> ref_tri_v1(1.0,0.0,FRAME_REFERENCE), tri_v1(2);
  Point<double> ref_tri_v2(0.0,1.0,FRAME_REFERENCE), tri_v2(2);
  
  // Mapping reference triangle vertices to the first cell in triMcell
  cout << "\nMapping reference triangle vertices to the first cell in triMcell:\n";
  tri_v0 = triMcell.mapToPhysicalCell(0,ref_tri_v0);
  tri_v1 = triMcell.mapToPhysicalCell(0,ref_tri_v1);
  tri_v2 = triMcell.mapToPhysicalCell(0,ref_tri_v2);
  //
  cout << "\t" << ref_tri_v0 << " maps to " << tri_v0 << endl;
  cout << "\t" << ref_tri_v1 << " maps to " << tri_v1 << endl;
  cout << "\t" << ref_tri_v2 << " maps to " << tri_v2 << endl;

  // Inverting the map: physical vertices to reference vertices
  cout << "\nInverting the map: physical vertices to reference vertices:\n";
  ref_tri_v0 = triMcell.mapToReferenceCell(0,tri_v0);
  ref_tri_v1 = triMcell.mapToReferenceCell(0,tri_v1);
  ref_tri_v2 = triMcell.mapToReferenceCell(0,tri_v2);
  //
  cout << "\t" << tri_v0 << " maps to " << ref_tri_v0 << endl;
  cout << "\t" << tri_v1 << " maps to " << ref_tri_v1 << endl;
  cout << "\t" << tri_v2 << " maps to " << ref_tri_v2 << endl;
  
  // Mapping reference triangle vertices to the third cell in triMcell
  cout << "\nMapping reference triangle vertices to the third cell in triMcell:\n";
  tri_v0 = triMcell.mapToPhysicalCell(2,ref_tri_v0);
  tri_v1 = triMcell.mapToPhysicalCell(2,ref_tri_v1);
  tri_v2 = triMcell.mapToPhysicalCell(2,ref_tri_v2);
  //
  cout << "\t" << ref_tri_v0 << " maps to " << tri_v0 << endl;
  cout << "\t" << ref_tri_v1 << " maps to " << tri_v1 << endl;
  cout << "\t" << ref_tri_v2 << " maps to " << tri_v2 << endl;
  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 3: Using charts with quadMcell                                       |\n"\
    << "===============================================================================\n";
  
  // Using charts with quadMcell: we will map the center of the reference quad cell
  Point<double> ref_quad_center(0.0,0.0,FRAME_REFERENCE), quad_center(2);
  
  // Mapping the center of the reference quad to the first quad in quadMcell and back:
  cout << "\nMapping the center of the reference quad to the first quad in quadMcell and back:\n";
  quad_center = quadMcell.mapToPhysicalCell(0,ref_quad_center);
  cout << "\t " << ref_quad_center << " maps to " << quad_center << endl;
  //
  ref_quad_center = quadMcell.mapToReferenceCell(0,quad_center);
  cout << "\t " << quad_center << " maps to " << ref_quad_center << endl;
  
  // Mapping the center of the reference quad to the second quad in quadMcell and back:
  cout << "\nMapping the center of the reference quad to the second quad in quadMcell and back:\n";
  quad_center = quadMcell.mapToPhysicalCell(1,ref_quad_center);
  cout << "\t " << ref_quad_center << " maps to " << quad_center << endl;
  //
  ref_quad_center = quadMcell.mapToReferenceCell(1,quad_center);
  cout << "\t " << quad_center << " maps to " << ref_quad_center << endl;

  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: Using charts with hexMcell                                       |\n"\
    << "===============================================================================\n";
  
  /*
  
  // Using charts with the hexMcell: we will map the vertices of the reference hex
  Point<double> ref_hex_v0(-1.0,-1.0,-1.0,FRAME_REFERENCE), hex_v0(3);
  Point<double> ref_hex_v1( 1.0,-1.0,-1.0,FRAME_REFERENCE), hex_v1(3);
  Point<double> ref_hex_v2( 1.0, 1.0,-1.0,FRAME_REFERENCE), hex_v2(3);
  Point<double> ref_hex_v3(-1.0, 1.0,-1.0,FRAME_REFERENCE), hex_v3(3);

  Point<double> ref_hex_v4(-1.0,-1.0, 1.0,FRAME_REFERENCE), hex_v4(3);
  Point<double> ref_hex_v5( 1.0,-1.0, 1.0,FRAME_REFERENCE), hex_v5(3);
  Point<double> ref_hex_v6( 1.0, 1.0, 1.0,FRAME_REFERENCE), hex_v6(3);
  Point<double> ref_hex_v7(-1.0, 1.0, 1.0,FRAME_REFERENCE), hex_v7(3);
  
  // Mapping reference hex vertices to the only cell in hexMcell
  cout << "\nMapping reference hex vertices to the third cell in triMcell:\n";
  hex_v0 = hexMcell.mapToPhysicalCell(0,ref_hex_v0);
  hex_v1 = hexMcell.mapToPhysicalCell(0,ref_hex_v1);
  hex_v2 = hexMcell.mapToPhysicalCell(0,ref_hex_v2);
  hex_v3 = hexMcell.mapToPhysicalCell(0,ref_hex_v3);

  hex_v4 = hexMcell.mapToPhysicalCell(0,ref_hex_v4);
  hex_v5 = hexMcell.mapToPhysicalCell(0,ref_hex_v5);
  hex_v6 = hexMcell.mapToPhysicalCell(0,ref_hex_v6);
  hex_v7 = hexMcell.mapToPhysicalCell(0,ref_hex_v7);
  
  cout << "\t" << ref_hex_v0 << " maps to " << hex_v0 << endl;
  cout << "\t" << ref_hex_v1 << " maps to " << hex_v1 << endl;
  cout << "\t" << ref_hex_v2 << " maps to " << hex_v2 << endl;
  cout << "\t" << ref_hex_v3 << " maps to " << hex_v3 << endl;
  
  cout << "\t" << ref_hex_v4 << " maps to " << hex_v4 << endl;
  cout << "\t" << ref_hex_v5 << " maps to " << hex_v5 << endl;
  cout << "\t" << ref_hex_v6 << " maps to " << hex_v6 << endl;
  cout << "\t" << ref_hex_v7 << " maps to " << hex_v7 << endl;
  
  // Mapping physical hex vertices back to the reference cell!
  cout << "\nMapping physical hex vertices back to the reference cell:\n";
  ref_hex_v0 = hexMcell.mapToReferenceCell(0,hex_v0);
  ref_hex_v1 = hexMcell.mapToReferenceCell(0,hex_v1);
  ref_hex_v2 = hexMcell.mapToReferenceCell(0,hex_v2);
  ref_hex_v3 = hexMcell.mapToReferenceCell(0,hex_v3);
  
  ref_hex_v4 = hexMcell.mapToReferenceCell(0,hex_v4);
  ref_hex_v5 = hexMcell.mapToReferenceCell(0,hex_v5);
  ref_hex_v6 = hexMcell.mapToReferenceCell(0,hex_v6);
  ref_hex_v7 = hexMcell.mapToReferenceCell(0,hex_v7);
  
  cout << "\t " << hex_v0 <<" maps to " << ref_hex_v0 << endl;
  cout << "\t " << hex_v1 <<" maps to " << ref_hex_v1 << endl;
  cout << "\t " << hex_v2 <<" maps to " << ref_hex_v2 << endl;
  cout << "\t " << hex_v3 <<" maps to " << ref_hex_v3 << endl;

  cout << "\t " << hex_v4 <<" maps to " << ref_hex_v4 << endl;
  cout << "\t " << hex_v5 <<" maps to " << ref_hex_v5 << endl;
  cout << "\t " << hex_v6 <<" maps to " << ref_hex_v6 << endl;
  cout << "\t " << hex_v7 <<" maps to " << ref_hex_v7 << endl;

   */
  
  return 0;
}
