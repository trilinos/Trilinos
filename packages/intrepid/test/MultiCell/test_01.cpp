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
#include "Intrepid_MultiCell.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the MultiCell class                        |\n" \
  << "|              Creating MultiCells and accessing their data                   |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: class MultiCell in 2D                                            |\n"\
  << "===============================================================================\n";
   
  // Define an array to store cell coordinates in an interleaved format. Cell type is triangle.
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
   
   // Suppose that edge signs are also needed. Define an array for the edge signs
   short triEdgeSigns[] = 
     {1, 1, -1,                     // edge signs for the edges of the first triangle
     -1, -1, 1,                     // edge signs for the edges of the second triangle
      1, 1, -1};                    // edge signs for the edges of the third triangle
   
   // Invoke a ctor that takes an array of subcell signs and the dimension of the subcell
   MultiCell<double> triMcell(
      3,                            // number of cells (triangles) in the multicell instance
      2,                            // ambient dimension
      CELL_TRI,                     // type of cells forming the multicell
      triNodes,                     // array with interleaved node coordinates
      triEdgeSigns,                 // array with edge signs
      1);                           // dimension of the subcells for which sign data is provided
   
   // Display the newly created MultiCell
   cout << triMcell << endl;         
   
   cout << "Testing multicell interface for the generating cell type...\n\n";
   
   cout << "\t type                   = " << triMcell.getMyCellType() << "\n";
   cout << "\t name                   = " << triMcell.getMyCellName() << "\n";
   cout << "\t ambient dimension      = " << triMcell.getMyAmbientDim() <<"\n";
   cout << "\t topological dimension  = " << triMcell.getMyTopologicalDim() << "\n";
   cout << "\t # of nodes             = " << triMcell.getMyNumNodes() << "\n"; 
   cout << "\t # of 0-subcells        = " << triMcell.getMyNumSubcells(0) << "\n";
   cout << "\t # of 1-subcells        = " << triMcell.getMyNumSubcells(1) << "\n";
   cout << "\t # of 2-subcells        = " << triMcell.getMyNumSubcells(2) << "\n";
   cout << "\t # of 3-subcells        = " << triMcell.getMyNumSubcells(3) << "\n";
   cout << "\t 1-subcell with index 0 = " << triMcell.getCellName(triMcell.getMySubcellType(1,0)) <<"\n";
   cout << "\t 1-subcell with index 1 = " << triMcell.getCellName(triMcell.getMySubcellType(1,1)) <<"\n";
   cout << "\t 1-subcell with index 2 = " << triMcell.getCellName(triMcell.getMySubcellType(1,2)) <<"\n";
   
   cout << "\t 2-subcell with index 0 = " << triMcell.getCellName(triMcell.getMySubcellType(2,0)) <<"\n\n";
   
   // Space for the node connectivities of subcells
   Teuchos::Array<int> subcellNodeConn;               
   triMcell.getMySubcellNodeIDs(1,                    // dimension of the subcell whose nodes we want
                                0,                    // local order (relative to cell template) of the subcell
                                subcellNodeConn);     // output - contains list of local node IDs
   
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
   
   // Using getVertex(i,j) to access the vertices of the cells in the MultiCell 
   cout << "Testing that triMcell[i][j] = triMcell.getVertex(i,j) for the cell with cell ID = 1...\n";
   Point<double> vertex_0a = triMcell.getVertex(1,0);    
   Point<double> vertex_1a = triMcell.getVertex(1,1);    
   Point<double> vertex_2a = triMcell.getVertex(1,2);    
   cout << "\t triMcell.getVertex(1,0) =" << vertex_0a << "\n";   
   cout << "\t triMcell.getVertex(1,1) =" << vertex_1a << "\n";   
   cout << "\t triMcell.getVertex(1,2) =" << vertex_2a << "\n\n";   

   // Using overloaded [] operator to access coordinates of the vertices stored as Point objects:
   cout << "Testing [] operator for vertex_0 \n";
   cout << "\t vertex_0[0] = "<< vertex_0[0] <<"\n";
   cout << "\t vertex_0[1] = "<< vertex_0[1] <<"\n\n";
   
   cout << "Testing [] operator for vertex_1 \n";
   cout << "\t vertex_1[0] = "<< vertex_1[0] <<"\n";
   cout << "\t vertex_1[1] = "<< vertex_1[1] <<"\n\n";
   
   cout \
   << "===============================================================================\n"\
   << "| EXAMPLE 2: class MultiCell in 3D                                            |\n"\
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
   
   // Use the ctor that takes both edge AND face sign data
   MultiCell<double> prismMcell(2,                // number of cells (prisms) in the multicell 
                                3,                // ambient dimension
                                CELL_TRIPRISM,    // generating cell type
                                prismnodes,       // list of node coordinates
                                prismEdgeSigns,   // list of edge signs
                                prismFaceSigns);  // list of face signs
   cout << prismMcell << endl;
   
   cout << "Testing multicell interface for the generating cell type...\n\n";
   cout << "\t type                   = " << prismMcell.getMyCellType() << "\n";
   cout << "\t name                   = " << prismMcell.getMyCellName() << "\n";
   cout << "\t ambient dimension      = " << prismMcell.getMyAmbientDim() <<"\n";
   cout << "\t topological dimension  = " << prismMcell.getMyTopologicalDim() << "\n";
   cout << "\t # of nodes             = " << prismMcell.getMyNumNodes() << "\n"; 
   cout << "\t # of 0-subcells        = " << prismMcell.getMyNumSubcells(0) << "\n";
   cout << "\t # of 1-subcells        = " << prismMcell.getMyNumSubcells(1) << "\n";
   cout << "\t # of 2-subcells        = " << prismMcell.getMyNumSubcells(2) << "\n";
   cout << "\t # of 3-subcells        = " << prismMcell.getMyNumSubcells(3) << "\n";
   cout << "\t 2-subcell with index 0 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,0)) <<"\n";
   cout << "\t 2-subcell with index 1 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,1)) <<"\n";
   cout << "\t 2-subcell with index 2 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,2)) <<"\n";
   cout << "\t 2-subcell with index 3 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,3)) <<"\n";
   cout << "\t 2-subcell with index 4 = " << prismMcell.getCellName(prismMcell.getMySubcellType(2,4)) <<"\n";
   cout << "\t 3-subcell with index 0 = " << prismMcell.getCellName(prismMcell.getMySubcellType(3,0)) <<"\n\n";
   
   // Accessing local node IDs (node lists) of the subcells
   prismMcell.getMySubcellNodeIDs(2,                  // dimension of the subcell whose nodes we want
                                  3,                  // local order (relative to cell template) of the subcell
                                  subcellNodeConn);   // output - contains list of local node IDs
   
   cout << "Node connectivity of the 2-subcell with ID = 3 -> { ";
   for(unsigned int i=0; i<subcellNodeConn.size(); i++) cout << subcellNodeConn[i]<<" ";
   cout << "}\n\n";
   
   // Using overloaded [] to access vertex coordinates
   cout << "Accessing vertex coordinates of cell with cellID = 1 ...\n";
   for(int i=0; i<prismMcell.getMyNumSubcells(0); i++){
     cout << "prismMcell[1]["<< i <<"] = " << prismMcell[1][i] << "\n";
   }
  
  return 0;
}
