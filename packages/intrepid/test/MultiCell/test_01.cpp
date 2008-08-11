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
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                              Unit Test MultiCell                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Test Accessing MultiCell data                                        |\n" \
    << "|     2) Test mapping to and from a physical cell                             |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;

  try{
    
    *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 1: MultiCell interface:      CELL_TRI generating cell                  |\n"\
    << "===============================================================================\n\n";
    
    // Vertex data for 3 CELL_TRI cells in (x,y) interleaved format
    double triNodes[] = 
      { 0.0, 0.0,                      // nodes of the first triangle
        1.0, 0.0,
        1.0, 1.0,
        0.0, 0.0,                      // nodes of the second triangle
        1.0, 1.0,
        0.0, 1.0,
        0.0, 1.0,                      // nodes of the third triangle
        1.5, 1.0,
        0.5, 2.0};
    
    // Edge & "face" sign and tag data for 3 cells in flat format
    short triEdgeSigns[] = {  1, 1,-1,  -1,-1, 1,  1, 1,-1}; 
    short triFaceSigns[] = { -1,-1, 1,   1, 1,-1, -1,-1, 1}; 
    short triEdgeTags[]  = {  0, 0, 1,   0, 1, 0,  0, 0, 0}; 
    short triFaceTags[]  = {  1, 1, 0,   1, 0, 1,  1, 1, 1}; 

    // Multicell consisting of 3 CELL_TRI cells with vertices specified in triNodes:
    MultiCell<double> triMcell(3,CELL_TRI,triNodes);    

    
    
    //============================================================================================//
    //  Testing getMySubcellVertexIDs                                                               //
    //============================================================================================//
    Teuchos::Array<int> subcellNodeIDs;
    int correctSubcellNodeID = -1;
    for(int cellID = 0; cellID < triMcell.getMyNumCells(); cellID++) {
      for(int subcellDim = 0; subcellDim < triMcell.getMyCellDim(); subcellDim++) {
        for(int subcellID = 0; subcellID < triMcell.getMyCellNumSubcells(subcellDim); subcellID++){
          
          triMcell.getMySubcellVertexIDs(subcellNodeIDs,subcellDim,subcellID);
          
          ECell subcellType = triMcell.getMySubcellType(subcellDim,subcellID);
          int numSubcellNodes = MultiCell<double>::getCellNumSubcells(subcellType, 0);
          
          for(int subcellNode = 0; subcellNode < numSubcellNodes; subcellNode++){
            
            // Hardcode correct subcell node IDs for a triangle generating cell
            if(subcellDim == 0) {
              
              // 0-subcells are the nodes, they have 1 node whose ID equals the subcellID
              correctSubcellNodeID = subcellID;
            }
            else if(subcellDim == 1) {
              
              // 1-subcells of CELL_TRI are the edges. They have 2 nodes.
              if(subcellNode == 0) {
                if(subcellID == 0) {
                  correctSubcellNodeID = 0;
                }
                else if(subcellID == 1) {
                  correctSubcellNodeID = 1;
                }
                else if(subcellID == 2) {
                  correctSubcellNodeID = 2;
                }
              }
              else if(subcellNode == 1){
                if(subcellID == 0) {
                  correctSubcellNodeID = 1;
                }
                else if(subcellID == 1) {
                  correctSubcellNodeID = 2;
                }
                else if(subcellID == 2) {
                  correctSubcellNodeID = 0;
                }
              }
            }
            
            // 2-subcell is the triangle itself
            else if(subcellDim == 2) {
              correctSubcellNodeID = subcellNode;
            }

            if(correctSubcellNodeID != subcellNodeIDs[subcellNode] ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Subcell Node ID retrived by getMySubcellVertexIDs = " << subcellNodeIDs[subcellNode] << "\n";
              *outStream << "                         Correct subcell node ID = " << correctSubcellNodeID << "\n";
              
            }
          }
        }
      }
    }
    
    //============================================================================================//
    // Test setting and accessing edge/face signs/tags.                                           //
    //============================================================================================//
    
    // For 2D generating cell edges and faces are 1-subcells and tests can be done together
    triMcell.setEdgeSigns(triEdgeSigns);  
    triMcell.setFaceSigns(triFaceSigns);  
    triMcell.setEdgeTags(triEdgeTags);
    triMcell.setFaceTags(triFaceTags);
    
    for(int cellID = 0; cellID < triMcell.getMyNumCells(); cellID++) {
      for(int edgeID = 0; edgeID < triMcell.getMyCellNumSubcells(1); edgeID++) {
        
        // Compute the flat index of the edge/face sign for the given cellID and edgeID
        int flatIndex = triMcell.getMyCellNumSubcells(1)*cellID + edgeID;
        
        // Check edge signs
        if(triEdgeSigns[flatIndex] != triMcell.getCellEdgeSigns(cellID)[edgeID]) {
          errorFlag++ ;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Edge sign value in data array = " << triEdgeSigns[flatIndex] << "\n";
          *outStream << " Edge sign value in MultiCell  = " << triMcell.getCellEdgeSigns(cellID)[edgeID] << "\n";
        }
        
        // Check edge tags
        if(triEdgeTags[flatIndex] != triMcell.getCellEdgeTags(cellID)[edgeID]) {
          errorFlag++ ;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Edge tag value in data array = " << triEdgeTags[flatIndex] << "\n";
          *outStream << " Edge tag value in MultiCell  = " << triMcell.getCellEdgeTags(cellID)[edgeID] << "\n";
        }
        
        // Check face signs
        if(triFaceSigns[flatIndex] != triMcell.getCellFaceSigns(cellID)[edgeID]) {
          errorFlag++ ;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Face sign value in data array = " << triFaceSigns[flatIndex] << "\n";
          *outStream << " Face sign value in MultiCell  = " << triMcell.getCellFaceSigns(cellID)[edgeID] << "\n";
        }
        
        // Check face tags
        if(triFaceTags[flatIndex] != triMcell.getCellFaceTags(cellID)[edgeID]) {
          errorFlag++ ;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Face tag value in data array = " << triFaceTags[flatIndex] << "\n";
          *outStream << " Face tag value in MultiCell  = " << triMcell.getCellFaceTags(cellID)[edgeID] << "\n";
        }
      }
    }
    
    //============================================================================================//
    //  Test getMyCellType and getCellType                                                        //
    //============================================================================================//
    if(triMcell.getMyCellType() != CELL_TRI ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Cell type by getMyCellType = " << MultiCell<double>::getCellName(triMcell.getMyCellType()) << "\n";
      *outStream << "          Correct cell type = " << MultiCell<double>::getCellName(CELL_TRI) << "\n";
    }
    if(MultiCell<double>::getCellType("Triangle") != CELL_TRI ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Cell type by getCellType   = " << MultiCell<double>::getCellName(MultiCell<double>::getCellType("Triangle")) << "\n";
      *outStream << "          Correct cell type = " << MultiCell<double>::getCellName(CELL_TRI) << "\n";
    }

    //============================================================================================//
    //  Test getMyCellName and getCellName                                                        //
    //============================================================================================//
    if(triMcell.getMyCellName() != "Triangle" ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Cell name by getMyCellName = " << triMcell.getMyCellName() << "\n";
      *outStream << "          Correct cell name = Triangle \n";
    }
    if(MultiCell<double>::getCellName(triMcell.getMyCellType()) != "Triangle" ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Cell name by getCellName   = " << MultiCell<double>::getCellName(triMcell.getMyCellType()) << "\n";
      *outStream << "          Correct cell name = Triangle \n";
    }
    
    //============================================================================================//
    //  Test getMyCellDim and getCellDim                                                          //
    //============================================================================================//
    if(triMcell.getMyCellDim() != 2 ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Generating cell dimension by getMyCellDim = " << triMcell.getMyCellDim() << "\n";
      *outStream << "         Correct generating cell dimension = 2 \n";
    }
    if(MultiCell<double>::getCellDim(triMcell.getMyCellType()) != 2 ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Generating cell dimension by getCellDim   = " << MultiCell<double>::getCellDim(triMcell.getMyCellType()) << "\n";
      *outStream << "         Correct generating cell dimension = 2 \n";
    }
    
    //============================================================================================//
    //  Test getMyCellNumSubcells and getCellNumSubcells                                                  //
    //============================================================================================//
    for(int subcellDim = 0; subcellDim < triMcell.getMyCellDim(); subcellDim++){
      int numSubcells=-1;
      if(subcellDim ==0 ){
        numSubcells = 3 ;
      }
      else if(subcellDim == 1) {
        numSubcells = 3;
      }
      else if(subcellDim == 2) {
        numSubcells = 1;
      }
      
      if(triMcell.getMyCellNumSubcells(subcellDim) != numSubcells ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " # of " << subcellDim << "-subcells by getMyCellNumSubcells = " << triMcell.getMyCellNumSubcells(subcellDim) << "\n";
        *outStream << " # of " << subcellDim << "-subcells should be           = " << numSubcells << "\n";
      }
      if(MultiCell<double>::getCellNumSubcells(CELL_TRI,subcellDim) != numSubcells ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " # of " << subcellDim << "-subcells by getCellNumSubcells   = " << MultiCell<double>::getCellNumSubcells(CELL_TRI,subcellDim) << "\n";
        *outStream << " # of " << subcellDim << "-subcells should be           = " << numSubcells << "\n";
      }
    }
    
    //============================================================================================//
    //  Test getMySubcellType and getSubcellType                                                  //
    //============================================================================================//
    for(int subcellDim = 0; subcellDim <= triMcell.getMyCellDim(); subcellDim++) {
      
      // Loop over subcell ID
      for(int subcellID = 0; subcellID < triMcell.getMyCellNumSubcells(subcellDim); subcellID++) {
        
        ECell subcellType = CELL_NODE;
        if(subcellDim == 0 ){
          subcellType = CELL_NODE;
        }
        else if(subcellDim == 1) {
          subcellType = CELL_EDGE;
        }
        else if(subcellDim == 2) {
          subcellType = CELL_TRI;
        }
        
        if(triMcell.getMySubcellType(subcellDim,subcellID) != subcellType ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Subcell type by getMySubcellType = " << MultiCell<double>::getCellName(triMcell.getMySubcellType(subcellDim,subcellID)) << "\n";
          *outStream << "           Subcell type should be = " << MultiCell<double>::getCellName(subcellType) << "\n";
        }
        if(MultiCell<double>::getSubcellType(CELL_TRI,subcellDim,subcellID) != subcellType ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Subcell type by getSubcellType   = " << MultiCell<double>::getCellName(MultiCell<double>::getSubcellType(CELL_TRI,subcellDim,subcellID)) << "\n";
          *outStream << "           Subcell type should be = " << MultiCell<double>::getCellName(subcellType) << "\n";
        }
      }
    }
    
    //============================================================================================//
    //  Test getCellVertex, getCellVertices and overloaded []                                     //
    //============================================================================================//
    Teuchos::Array< Point<double> > cellVertices;
    int numCells     = triMcell.getMyNumCells();
    int cellNumNodes = triMcell.getMyCellNumSubcells(0);
    int cellDim      = triMcell.getMyCellDim();
    
    for(int cellID = 0; cellID < numCells; cellID++) {
      for(int vertexID = 0; vertexID < cellNumNodes; vertexID++) {
        
        Point<double> vertex = triMcell.getCellVertex(cellID,vertexID);
        for(int dim = 0; dim < cellDim; dim++ ) {
          
          // Compute index of the dim-th coordinate of the node with vertexID in the cell with cellID
          int flatIndex = cellID*(cellNumNodes*cellDim) + vertexID*cellDim + dim;
          
          if(triMcell.getCellVertex(cellID,vertexID)[dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by getCellVertex  = " << triMcell.getCellVertex(cellID,vertexID)[dim]<< "\n";
            *outStream << "          Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }
          cellVertices = triMcell.getCellVertices(cellID);
          if(cellVertices[vertexID][dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by getCellVertices  = " << cellVertices[vertexID][dim]<< "\n";
            *outStream << "            Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }
          if(triMcell[cellID][vertexID][dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by overloaded []  = " << triMcell[cellID][vertexID][dim]<< "\n";
            *outStream << "          Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }          
        }
      }
    }
    
    // Same test but for a MultiCell defined by a ctor that takes FieldContainer of vertices:
    FieldContainer<double> triNodesFC(3, 3, 2);
    triNodesFC.setValues(triNodes, 18);
    MultiCell<double> triMcellFC(CELL_TRI, triNodesFC);   
    
    numCells     = triMcellFC.getMyNumCells();
    cellNumNodes = triMcellFC.getMyCellNumSubcells(0);
    cellDim      = triMcellFC.getMyCellDim();
    
    for(int cellID = 0; cellID < numCells; cellID++) {
      for(int vertexID = 0; vertexID < cellNumNodes; vertexID++) {
        
        Point<double> vertex = triMcellFC.getCellVertex(cellID,vertexID);
        for(int dim = 0; dim < cellDim; dim++ ) {
          
          // Compute index of the dim-th coordinate of the node with vertexID in the cell with cellID
          int flatIndex = cellID*(cellNumNodes*cellDim) + vertexID*cellDim + dim;
          
          if(triMcellFC.getCellVertex(cellID, vertexID)[dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by getCellVertex  = " << triMcellFC.getCellVertex(cellID,vertexID)[dim]<< "\n";
            *outStream << "          Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }
          cellVertices = triMcellFC.getCellVertices(cellID);
          if(cellVertices[vertexID][dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by getCellVertices  = " << cellVertices[vertexID][dim]<< "\n";
            *outStream << "            Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }
          if(triMcellFC[cellID][vertexID][dim] != triNodes[flatIndex]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Node coordinate by overloaded []  = " << triMcellFC[cellID][vertexID][dim]<< "\n";
            *outStream << "          Correct node coordinate  = " << triNodes[flatIndex] << "\n";
          }          
        }
      }
    }
    
    
    
    //============================================================================================//
    //  Test mapToPhysicalCell and mapToReferenceCell (implicit test of setAtlas() and jacobian)  //
    //============================================================================================//
    
    // The test makes random points in the reference CELL_TRI, maps them to a physical cell and back.
    triMcell.setAtlas();
    Point<double> randomRefPoint(2,FRAME_REFERENCE);
    
    // set seed
    std::srand( (long)(&randomRefPoint) );
    for(int numPts = 0; numPts < 500; numPts++){
      // Use when time() becomes platform-independent.
      // std::srand( std::time(NULL)*numPts );			
      
      // Two random numbers between 0 and 1
      double rx = ((double)std::rand())/RAND_MAX;	
      double ry = ((double)std::rand())/RAND_MAX;	
      
      // Make a random point inside the reference CELL_TRI
      if( rx + ry <= 1.0 ) {
        randomRefPoint = Point<double>::Point(rx,ry,FRAME_REFERENCE); 
      }
      else {
        double sum = rx + ry;
        rx = rx/sum;
        ry = ry/sum;
      }
      
      //Loop over the cells in the MultiCell
      for(int cellID = 0; cellID < triMcell.getMyNumCells(); cellID++) {
        Point<double> randomPhysPoint = triMcell.mapToPhysicalCell(cellID,randomRefPoint);
        Point<double> tempPoint       = triMcell.mapToReferenceCell(cellID,randomPhysPoint);
        
        if( randomRefPoint.distance(tempPoint) > INTREPID_TOL ) {
          errorFlag++; 
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Original reference point " << randomRefPoint << "\n";
          *outStream << " Result after mapToPhysicalCell followed by mapToReferenceCell: " << tempPoint << "\n";
          *outStream << " Distance between the two points " << randomRefPoint.distance(tempPoint) 
            << "exceeds the default INTREPID_TOL value of " << INTREPID_TOL << "\n"; 
        }
      }
    }
    
    //============================================================================================//
    //  Test inReferenceCell (a static member function)                                           //
    //============================================================================================//
    
    // Create reference points near the boundaries of the reference cells
    Point<double> p_in_edge(1.0-INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_quad(1.0,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_tri(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_hex(1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_tet(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,2.0*INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_prism(0.5,0.25,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
    Point<double> p_in_pyramid(-INTREPID_EPSILON,INTREPID_EPSILON,(1.0-INTREPID_EPSILON),FRAME_REFERENCE);
    
    // Check if the points are in their respective reference cells
    bool in_edge    = MultiCell<double>::inReferenceCell(CELL_EDGE, p_in_edge);   
    bool in_tri     = MultiCell<double>::inReferenceCell(CELL_TRI, p_in_tri);
    bool in_quad    = MultiCell<double>::inReferenceCell(CELL_QUAD, p_in_quad);
    bool in_tet     = MultiCell<double>::inReferenceCell(CELL_TET, p_in_tet);
    bool in_hex     = MultiCell<double>::inReferenceCell(CELL_HEX, p_in_hex);
    bool in_prism   = MultiCell<double>::inReferenceCell(CELL_TRIPRISM, p_in_prism);
    bool in_pyramid = MultiCell<double>::inReferenceCell(CELL_PYRAMID, p_in_pyramid);
    
    if( !in_edge ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_edge << " incorrectly classified as not in the closed reference CELL_EDGE \n";
    }
    if( !in_tri ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_tri << " incorrectly classified as not in the closed reference CELL_TRI \n";
    }
    if( !in_quad ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_quad << " incorrectly classified as not in the closed reference CELL_QUAD \n";
    }
    if( !in_tet ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_tet << " incorrectly classified as not in the closed reference CELL_TET \n";
    }
    if( !in_hex ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_hex << " incorrectly classified as not in the closed reference CELL_HEX \n";
    }
    if( !in_prism ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_prism << " incorrectly classified as not in the closed reference CELL_TRIPRISM \n";
    }
    if( !in_pyramid ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_pyramid << " incorrectly classified as not in the closed reference CELL_PYRAMID \n";
    }
    
    // Now make 1,2 and 3D points with very small coefficients, but larger than threshold
    double small = 2.0*INTREPID_THRESHOLD;
    Point<double> p_eps_1D(small,FRAME_REFERENCE);
    Point<double> p_eps_2D(small,small,FRAME_REFERENCE);
    Point<double> p_eps_3D(small,small,small,FRAME_REFERENCE);
    
    // Add these points to the good reference points above:
    p_in_edge    += p_eps_1D;
    p_in_tri     += p_eps_2D;
    p_in_quad    += p_eps_2D;
    p_in_tet     += p_eps_3D;
    p_in_hex     += p_eps_3D;
    p_in_prism   += p_eps_3D;
    p_in_pyramid += p_eps_3D;
    
    // Now check again if the points are in their respective reference cells.
    in_edge    = MultiCell<double>::inReferenceCell(CELL_EDGE, p_in_edge);
    in_tri     = MultiCell<double>::inReferenceCell(CELL_TRI, p_in_tri);
    in_quad    = MultiCell<double>::inReferenceCell(CELL_QUAD, p_in_quad);
    in_tet     = MultiCell<double>::inReferenceCell(CELL_TET, p_in_tet);
    in_hex     = MultiCell<double>::inReferenceCell(CELL_HEX, p_in_hex);
    in_prism   = MultiCell<double>::inReferenceCell(CELL_TRIPRISM, p_in_prism);
    in_pyramid = MultiCell<double>::inReferenceCell(CELL_PYRAMID, p_in_pyramid);
    
    if( in_edge ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_edge << " incorrectly classified as in the closed reference CELL_EDGE \n";
    }
    if( in_tri ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_tri << " incorrectly classified as in the closed reference CELL_TRI \n";
    }
    if( in_quad ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_quad << " incorrectly classified as in the closed reference CELL_QUAD \n";
    }
    if( in_tet ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_tet << " incorrectly classified as in the closed reference CELL_TET \n";
    }
    if( in_hex ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_hex << " incorrectly classified as in the closed reference CELL_HEX \n";
    }
    if( in_prism ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_prism << " incorrectly classified as in the closed reference CELL_TRIPRISM \n";
    }
    if( in_pyramid ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << p_in_pyramid << " incorrectly classified as in the closed reference CELL_PYRAMID \n";
    }

  }// try
  
  //============================================================================================//
  // Wrap up test: check if the test broke down unexpectedly due to an exception                //
  //============================================================================================//
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
