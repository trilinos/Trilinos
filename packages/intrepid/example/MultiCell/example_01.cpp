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
\brief  Example of the MultiCell class.
\author Created by P. Bochev and D. Ridzal
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Intrepid_Cell.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the MultiCell class                        |\n" \
  << "|                                                                             |\n" \
  << "|     1) Creating MultiCells and accessing their data                         |\n" \
  << "|     2) Testing points for inclusion in reference and physical cells         |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: class MultiCell in 2D: CELL_TRI generating cell                  |\n"\
  << "===============================================================================\n";
   
  // Vertex data for 3 CELL_TRI cells in (x,y) interleaved format.
   double triNodes[] = 
    {0.0, 0.0,                      // nodes of the first triangle
     1.0, 0.0,
     1.0, 1.0,
     0.0, 0.0,                      // nodes of the second triangle
     1.0, 1.0,
     0.0, 1.0,
     0.0, 1.0,                      // nodes of the third triangle
     1.5, 1.0,
     0.5, 2.0};
   
   // Edge sign data for 3 cells in flat format
   short triEdgeSigns[] = 
     {1, 1, -1,                     // edge signs for the edges of the first triangle
     -1, -1, 1,                     // edge signs for the edges of the second triangle
      1, 1, -1};                    // edge signs for the edges of the third triangle
    
   // Create multicell consisting of 3 CELL_TRI cells with vertices specified in triNodes:
   MultiCell<double> triMcell(3,CELL_TRI,triNodes);   

   // Set the edge signs
   triMcell.setEdgeSigns(triEdgeSigns);                  
   
   // Display the newly created MultiCell
   cout << triMcell << endl;     
   
   
   // Uisng MultiCell constructor with ArrayType
   FieldContainer<double> triNodesFC(3, 3, 2);
   triNodesFC.setValues(triNodes, 18);
   
   MultiCell<double> triMcellFC(CELL_TRI, triNodesFC);
   cout << triMcellFC << endl;
   
   

   cout << "Testing multicell interface for the generating cell type...\n\n";
   cout << "\t # of cells in MCell    = " << triMcell.getMyNumCells() << "\n"; 
   cout << "\t generating cell type   = " << triMcell.getMyCellType() << "\n";
   cout << "\t generating cell name   = " << triMcell.getMyCellName() << "\n";
   cout << "\t generating cell dim    = " << triMcell.getMyCellDim() <<"\n";
   cout << "\t # of 0-subcells        = " << triMcell.getMyCellNumSubcells(0) << "\n";
   cout << "\t # of 1-subcells        = " << triMcell.getMyCellNumSubcells(1) << "\n";
   cout << "\t # of 2-subcells        = " << triMcell.getMyCellNumSubcells(2) << "\n";
   cout << "\t # of 3-subcells        = " << triMcell.getMyCellNumSubcells(3) << "\n";
   cout << "\t 1-subcell with index 0 = " << triMcell.getCellName(triMcell.getMySubcellType(1,0)) <<"\n";
   cout << "\t 1-subcell with index 1 = " << triMcell.getCellName(triMcell.getMySubcellType(1,1)) <<"\n";
   cout << "\t 1-subcell with index 2 = " << triMcell.getCellName(triMcell.getMySubcellType(1,2)) <<"\n";
   cout << "\t 2-subcell with index 0 = " << triMcell.getCellName(triMcell.getMySubcellType(2,0)) <<"\n\n";
   
   // Space for the node connectivities of subcells
   Teuchos::Array<int> subcellNodeConn;               
   triMcell.getMySubcellVertexIDs(subcellNodeConn,      // output - contains list of local node IDs
                                1,                    // dimension of the subcell whose nodes we want
                                0);                   // local order (relative to cell template) of the subcell
   
   cout << "Node connectivity of the 0th 1-subcell -> { ";
   for(unsigned int i=0; i<subcellNodeConn.size(); i++) cout << subcellNodeConn[i]<<" ";
   cout << "}\n\n";
   
   // Using the overloaded [] operator to access the vertices of the cells in the MultiCell
   cout << "Using the overloaded [] operator to access the vertices of the cell with cell ID = 1...\n";
   Point<double> vertex_0 = triMcell[1][0];            // vertex 0 of cell with cellID = 1
   Point<double> vertex_1 = triMcell[1][1];            // vertex 1 of cell with cellID = 1
   Point<double> vertex_2 = triMcell[1][2];            // vertex 2 of cell with cellID = 1
   cout << "\t triMcell[1][0] =" << vertex_0 << "\n";   // << is overloaded for Point objects
   cout << "\t triMcell[1][1] =" << vertex_1 << "\n";
   cout << "\t triMcell[1][2] =" << vertex_2 << "\n\n";
   
   // Using getCellVertex(i,j) to access the vertices of the cells in the MultiCell 
   cout << "Testing that triMcell[i][j] = triMcell.getCellVertex(i,j) for the cell with cell ID = 1...\n";
   Point<double> vertex_0a = triMcell.getCellVertex(1,0);    
   Point<double> vertex_1a = triMcell.getCellVertex(1,1);    
   Point<double> vertex_2a = triMcell.getCellVertex(1,2);    
   cout << "\t triMcell.getCellVertex(1,0) =" << vertex_0a << "\n";   
   cout << "\t triMcell.getCellVertex(1,1) =" << vertex_1a << "\n";   
   cout << "\t triMcell.getCellVertex(1,2) =" << vertex_2a << "\n\n";   

   // Using overloaded [] operator to access coordinates of the vertices stored as Point objects:
   cout << "Testing [] operator for vertex_0 \n";
   cout << "\t vertex_0[0] = "<< vertex_0[0] <<"\n";
   cout << "\t vertex_0[1] = "<< vertex_0[1] <<"\n\n";
   
   cout << "Testing [] operator for vertex_1 \n";
   cout << "\t vertex_1[0] = "<< vertex_1[0] <<"\n";
   cout << "\t vertex_1[1] = "<< vertex_1[1] <<"\n\n";
   
   cout \
   << "===============================================================================\n"\
   << "| EXAMPLE 2: class MultiCell in 3D: CELL_TRIPRISM generating cell             |\n"\
   << "===============================================================================\n";
   
   // Define an array to store cell coordinates in an interleaved format. Cell type is TRIPRISM.
   double prismnodes[] = 
    {0.0, 0.0, 0.0,         // nodes of the first prism
     1.0, 0.0, 0.0,
     0.5, 0.5, 0.0,
     0.0, 0.0, 1.0,
     1.0, 0.0, 1.0,
     0.5, 0.5, 1.0, 
     0.0, 0.0, 1.0,         // nodes of the second prism
     1.0, 0.0, 1.0,         // = nodes of the first prism offset by 1 in the z-axes
     0.5, 0.5, 1.0,         // i.e., two prisms stacked on top of each other
     0.0, 0.0, 2.0,
     1.0, 0.0, 2.0,
     0.5, 0.5, 2.0};
   
   // Define an array with the edge signs for the two prisms
   short prismEdgeSigns[] = 
    { 1,  1,  1, 1, 1, 1,  1,  1,  1,             // signs for the edges of prism #1
     -1, -1, -1, 1, 1, 1, -1, -1, -1};            // signs for the edges of prism #2
   
   // Define an array with the face signs for the two prisms
   short prismFaceSigns[] = 
    { 1,  1,  1,  1,  1,                          // signs for the faces of prism #1
     -1, -1, -1, -1, -1};                         // signs for the faces of prism #2
   
   // Define a multicell consisting of two prisms
   MultiCell<double> prismMcell(2, CELL_TRIPRISM, prismnodes);      
   
   // Set edge and face signs:
   prismMcell.setEdgeSigns(prismEdgeSigns);
   prismMcell.setFaceSigns(prismFaceSigns);
   cout << prismMcell << endl;
   
   cout << "Testing multicell interface for the generating cell type...\n\n";
   cout << "\t # of cells in MCell    = " << prismMcell.getMyNumCells() << "\n"; 
   cout << "\t generating cell type   = " << prismMcell.getMyCellType() << "\n";
   cout << "\t generating cell name   = " << prismMcell.getMyCellName() << "\n";
   cout << "\t ambient dimension      = " << prismMcell.getMyCellDim() <<"\n";
   cout << "\t # of 0-subcells        = " << prismMcell.getMyCellNumSubcells(0) << "\n";
   cout << "\t # of 1-subcells        = " << prismMcell.getMyCellNumSubcells(1) << "\n";
   cout << "\t # of 2-subcells        = " << prismMcell.getMyCellNumSubcells(2) << "\n";
   cout << "\t # of 3-subcells        = " << prismMcell.getMyCellNumSubcells(3) << "\n";
   cout << "\t 2-subcell with index 0 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,0)) <<"\n";
   cout << "\t 2-subcell with index 1 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,1)) <<"\n";
   cout << "\t 2-subcell with index 2 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,2)) <<"\n";
   cout << "\t 2-subcell with index 3 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,3)) <<"\n";
   cout << "\t 2-subcell with index 4 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,4)) <<"\n";
   cout << "\t 3-subcell with index 0 = " << prismMcell.getCellName(prismMcell.getMySubcellType(3,0)) <<"\n\n";
   
   // Accessing local node IDs (node lists) of the subcells
   prismMcell.getMySubcellVertexIDs(subcellNodeConn,    // output - contains list of local node IDs
                                  2,                  // dimension of the subcell whose nodes we want
                                  3);                  // local order (relative to cell template) of the subcell
   
   cout << "Node connectivity of the 2-subcell with ID = 3 -> { ";
   for(unsigned int i=0; i<subcellNodeConn.size(); i++) cout << subcellNodeConn[i]<<" ";
   cout << "}\n\n";
   
   // Using overloaded [] to access vertex coordinates
   cout << "Accessing vertex coordinates of cell with cellID = 1 ...\n";
   for(int i=0; i<prismMcell.getMyCellNumSubcells(0); i++){
     cout << "prismMcell[1]["<< i <<"] = " << prismMcell[1][i] << "\n";
   }
   
   cout << "\n" \
     << "===============================================================================\n"\
     << "| EXAMPLE 3: Using inclusion test methods in MultiCell                             |\n"\
     << "===============================================================================\n";
   
   //The static member function inReferenceCell can be used without a MultiCell instantiation
   cout << "\nUsing static member inReferenceCell to check if a Point is inside a reference cell: " << endl;
   
   // Points below illustrate using constructors that take 1,2 or 3 point coordinates 
   // and the user sets the frame kind explicitely to override the default (FRAME_PHYSICAL).
   
   // 1D point that is close to the right endpoint of the reference edge cell
   Point<double> p_in_edge(1.0-INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 2D point that is close to the top right corner of the reference quad cell
   Point<double> p_in_quad(1.0,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 2D point that is close to the midpoint of the slanted edge of the reference tri cell
   Point<double> p_in_tri(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 3D point that is close to 0th vertex of the reference hex cell
   Point<double> p_in_hex(1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 3D point that is close to the slanted face of the reference tet cell
   Point<double> p_in_tet(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,2.0*INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 3D point close to the top face of the reference prism 
   Point<double> p_in_prism(0.5,0.25,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
   
   // 3D point close to the top of the reference pyramid
   Point<double> p_in_pyramid(-INTREPID_EPSILON,INTREPID_EPSILON,(1.0-INTREPID_EPSILON),FRAME_REFERENCE);
   
   // Check if the points are in their respective reference cells
   bool in_edge    = MultiCell<double>::inReferenceCell(CELL_EDGE, p_in_edge);   
   bool in_tri     = MultiCell<double>::inReferenceCell(CELL_TRI, p_in_tri);
   bool in_quad    = MultiCell<double>::inReferenceCell(CELL_QUAD, p_in_quad);
   bool in_tet     = MultiCell<double>::inReferenceCell(CELL_TET, p_in_tet);
   bool in_hex     = MultiCell<double>::inReferenceCell(CELL_HEX, p_in_hex);
   bool in_prism   = MultiCell<double>::inReferenceCell(CELL_TRIPRISM, p_in_prism);
   bool in_pyramid = MultiCell<double>::inReferenceCell(CELL_PYRAMID, p_in_pyramid);
   //
   if(in_edge) {
     cout <<  p_in_edge << " is inside reference edge " << endl;
   }
   if(in_tri) {
     cout << p_in_tri << " is inside reference triangle " << endl;
   }
   if(in_quad) {
     cout << p_in_quad << " is inside reference quad " << endl;
   }
   if(in_tet) {
     cout << p_in_tet << " is inside reference tet " << endl;
   }
   if(in_hex) {
     cout << p_in_hex << " is inside reference hex " << endl;
   }
   if(in_prism) {
     cout << p_in_prism << " is inside reference prism " << endl;
   }
   if(in_pyramid) {
     cout << p_in_pyramid << " is inside reference pyramid " << endl;
   }
   
   // Now make 1,2 and 3D points with very small coefficients, but larger than threshold
   double small = 2.0*INTREPID_THRESHOLD;
   Point<double> p_eps_1D(small,FRAME_REFERENCE);
   Point<double> p_eps_2D(small,small,FRAME_REFERENCE);
   Point<double> p_eps_3D(small,small,small,FRAME_REFERENCE);
   
   // Add these points to the good reference points above:
   cout << "\nAdding small perturbations to these points..." << endl;
   p_in_edge    += p_eps_1D;
   p_in_tri     += p_eps_2D;
   p_in_quad    += p_eps_2D;
   p_in_tet     += p_eps_3D;
   p_in_hex     += p_eps_3D;
   p_in_prism   += p_eps_3D;
   p_in_pyramid += p_eps_3D;
   
   // Now check again if the points are in their respective reference cells.
   cout << "\nChecking if the perturbed Points belongs to reference cell: " << endl;
   in_edge    = MultiCell<double>::inReferenceCell(CELL_EDGE, p_in_edge);
   in_tri     = MultiCell<double>::inReferenceCell(CELL_TRI, p_in_tri);
   in_quad    = MultiCell<double>::inReferenceCell(CELL_QUAD, p_in_quad);
   in_tet     = MultiCell<double>::inReferenceCell(CELL_TET, p_in_tet);
   in_hex     = MultiCell<double>::inReferenceCell(CELL_HEX, p_in_hex);
   in_prism   = MultiCell<double>::inReferenceCell(CELL_TRIPRISM, p_in_prism);
   in_pyramid = MultiCell<double>::inReferenceCell(CELL_PYRAMID, p_in_pyramid);
   
   //
   if(! in_edge) {
     cout <<  p_in_edge << " is NOT inside reference edge " << endl;
   }
   if(! in_tri) {
     cout << p_in_tri << " is NOT inside reference triangle " << endl;
   }
   if(! in_quad) {
     cout << p_in_quad << " is NOT inside reference quad " << endl;
   }
   if(! in_tet) {
     cout << p_in_tet << " is NOT inside reference tet " << endl;
   }
   if(! in_hex) {
     cout << p_in_hex << " is NOT inside reference hex " << endl;
   }
   if(! in_prism) {
     cout << p_in_prism << " is NOT inside reference prism " << endl;
   }
   if(! in_pyramid) {
     cout << p_in_pyramid << " is NOT inside reference pyramid " << endl;
   }
   cout << "You may need at least 16 digits to see the added perturbation!" << endl;


   cout << "\n" \
     << "===============================================================================\n"\
     << "| EXAMPLE 4: Specialized methods of the Cell class                            |\n"\
     << "===============================================================================\n";
   
   // Construct a single CELL_TRI cell with vertices from triNodes (only the first set is used)
   Cell<double> triCell(CELL_TRI,triNodes);                 
   
   // Display the newly created Cell
   cout << triCell << endl;         
   
   //triCell.getMyEdgeSigns();
   cout << "Get vertex 2:\n" << triCell.getCellVertex(2) << "\n\n";
   triCell.getCellVertices();
   triCell.setAtlas();
   triCell.setChart();
   Point<double> refPt(0.4,0.4,FRAME_REFERENCE);
   cout << "Transformations based on: " << refPt << "\n";
   cout << " Jacobian:\n" << triCell.jacobian(refPt) << endl;
   cout << " Map to physical point:" << triCell.mapToPhysicalCell(refPt) << endl;
   cout << " Map back to reference point:" << triCell.mapToReferenceCell(triCell.mapToPhysicalCell(refPt)) << endl;
   cout << " Is physical point inside physical cell?";
   if (triCell.inPhysicalCell(triCell.mapToPhysicalCell(refPt)) == 0)
     cout << " YES.\n";
  else
    cout << " NO.\n";
 
  // Using the overloaded [] operator to access the vertices of the Cell
  cout << "Using the overloaded [] operator to access the vertices of the cell \n";
  vertex_0 = triCell[0][0];            // vertex 0 of the cell
  vertex_1 = triCell[0][1];            // vertex 1 of the cell
  vertex_2 = triCell[0][2];            // vertex 2 of the cell
  cout << "\t triCell[0][0] =" << vertex_0 << "\n";   // << is overloaded for Point objects
  cout << "\t triCell[0][1] =" << vertex_1 << "\n";
  cout << "\t triCell[0][2] =" << vertex_2 << "\n\n";
  
  // Using getCellVertex(i) to access the vertices of the cell:
  cout << "Testing that triCell[0][j] = triCell.getCellVertex(i) \n";
  vertex_0a = triCell.getCellVertex(0);    
  vertex_1a = triCell.getCellVertex(1);    
  vertex_2a = triCell.getCellVertex(2);    
  cout << "\t triCell.getCellVertex(0) =" << vertex_0a << "\n";   
  cout << "\t triCell.getCellVertex(1) =" << vertex_1a << "\n";   
  cout << "\t triCell.getCellVertex(2) =" << vertex_2a << "\n\n";   
  
  
  return 0;
}
