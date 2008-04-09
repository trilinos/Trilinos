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


/** \file   Intrepid_MultiCellDef.hpp
\brief  Definition file for the Intrepid::MultiCell class.
\author Created by P. Bochev, D. Ridzal, and D. Day.
*/

#include "Intrepid_CellTemplates.hpp"

namespace Intrepid {
  
  template<class Scalar>
  const char* MultiCell<Scalar>::cellNames_[] = {
    "Node",
    "Edge",
    "Triangle",
    "Quadrilateral",
    "Tetrahedron",
    "Hexahedron",
    "Pyramid",
    "Pentagon",
    "Hexagon",
    "Heptagon",
    "Octagon",
    "Nonagon",
    "Decagon",
    "TriPrism",
    "PentaPrism",
    "HexaPrism",
    "Max_Canonical_Cell",
    "Poly0",
    "Poly1",
    "Poly2",
    "Poly3",
    "Poly4",
    "Poly5",
    "Poly6",
    "Poly7",
    "Poly8",
    "Poly9",
    "Max_Cell_Types"
  };
  
  
  
  template<class Scalar>
  MultiCell<Scalar>::MultiCell(const int      numCells,
                               const ECell    generatingCellType,
                               const Scalar*  vertices) : 
  numCells_(numCells), myCellType_(generatingCellType) {
   

#ifdef HAVE_INTREPID_DEBUG
    
    // Number of cells has to be positive
    TEST_FOR_EXCEPTION( !(0 < numCells), 
                        std::invalid_argument, 
                        ">>> ERROR (MultiCell): Number of cells cannot be 0 or negative.")
    
    // Admissible generating cell types start with CELL_EDGE and end with CELL_POLY9
    TEST_FOR_EXCEPTION( !( (CELL_EDGE <= generatingCellType) && (generatingCellType < CELL_MAX) ),
                        std::invalid_argument,
                        ">>> ERROR (MultiCell): Invalid generating cell type");  
#endif
    
    // Find out how many nodes does the generating cell have
    int numNodesPerCell = this -> getMyCellNumNodes();
    int ambientDim      = this -> getMyCellDim();
    
    // Array of points sized to hold as many points as there are vertices in the generating cell
    Point<Scalar> tempPoint(ambientDim);
    Teuchos::Array< Point<Scalar> > tempPointArray(numNodesPerCell, tempPoint);
    
    // Resize data member to hold as many such arrays as there are cells in the MultiCell
    vertices_.assign(numCells_, tempPointArray);
    
    // fill in the coordinate info
    for (int i=0; i < numCells_; i++) {
      for (int j=0; j < numNodesPerCell; j++) {
        
        // do some pointer arithmetic to extract vertex coordinates (they are interleaved)
        Point<Scalar> workPoint(vertices + (i*numNodesPerCell+j)*ambientDim, ambientDim);
        vertices_[i][j] = workPoint;
      }
    }
  }
  
  
  
  template<class Scalar>
  void MultiCell<Scalar>::setEdgeSigns(const short*  edgeSigns) {
    
    int numEdgesPerCell = 0;
    int ambientDim      = this -> getMyCellDim();

    // Edge signs are not meaningful unless dimension of the generating cell is 2 or 3
    if(ambientDim == 2 || ambientDim == 3) {
      numEdgesPerCell = this -> getMyNumSubcells(1);
    }
    else {
      TEST_FOR_EXCEPTION( (ambientDim == 1), 
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): Edge signs cannot be defined for 1-dimensional cells");  
    }
    
    // Allocate space for the edge signs
    Teuchos::Array<short> tempEdgeArray(numEdgesPerCell);
    edgeSigns_.assign(numCells_, tempEdgeArray);
    
    // Fill in the edge signs
    for (int i=0; i<numCells_; i++) {
      edgeSigns_[i].assign(edgeSigns+i*numEdgesPerCell, edgeSigns+(i+1)*numEdgesPerCell);
    }
  }
  
  
  
  template<class Scalar>
  void  MultiCell<Scalar>::setFaceSigns(const short*   faceSigns) {
    
    int numFacesPerCell = 0;
    int ambientDim      = this -> getMyCellDim();

    if(ambientDim == 3) {
      
      // in 3D cell "faces" are 2 dimensional subcells of the cell (true faces)
      numFacesPerCell = this -> getMyNumSubcells(2);
    }
    else {
      if(ambientDim == 2) {
        
        // In 2D "faces" are 1-dimensional subcells (edges)
        numFacesPerCell = this -> getMyNumSubcells(1);
      }
      else {
        
        // Face signs are not meaningful unless dimension of the generating cell is 2 or 3
        TEST_FOR_EXCEPTION( (ambientDim == 1), 
                            std::invalid_argument,
                            ">>> ERROR (MultiCell): Face signs cannot be defined for 1-dimensional cells"); 
      }
    }
    
    // Allocate a zero array of 'short' arrays for face signs:
    Teuchos::Array<short> tempFaceArray(numFacesPerCell);
    faceSigns_.assign(numCells_, tempFaceArray);
    
    // Fill in the face signs
    for (int i=0; i<numCells_; i++) {
      faceSigns_[i].assign(faceSigns+i*numFacesPerCell, faceSigns+(i+1)*numFacesPerCell);
    }
  }
  
  
  
  template<class Scalar>
    void MultiCell<Scalar>::setEdgeTags(const short*  edgeTags) {
      
      int numEdgesPerCell = 0;
      int ambientDim      = this -> getMyCellDim();
      
      // Edge tags are not meaningful unless dimension of the generating cell is 2 or 3
      if(ambientDim == 2 || ambientDim == 3) {
        numEdgesPerCell = this -> getMyNumSubcells(1);
      }
      else {
        TEST_FOR_EXCEPTION( (ambientDim == 1), 
                            std::invalid_argument,
                            ">>> ERROR (MultiCell): Edge tags cannot be defined for 1-dimensional cells");  
      }
      
      // Allocate space for the edge tags
      Teuchos::Array<short> tempEdgeArray(numEdgesPerCell);
      edgeTags_.assign(numCells_, tempEdgeArray);
      
      // Fill in the edge tags
      for (int i=0; i<numCells_; i++) {
        edgeTags_[i].assign(edgeTags+i*numEdgesPerCell, edgeTags+(i+1)*numEdgesPerCell);
      }
    }
  
  
  
  template<class Scalar>
    void  MultiCell<Scalar>::setFaceTags(const short*   faceTags) {
      
      int numFacesPerCell = 0;
      int ambientDim      = this -> getMyCellDim();
      
      if(ambientDim == 3) {
        
        // in 3D cell "faces" are 2 dimensional subcells of the cell (true faces)
        numFacesPerCell = this -> getMyNumSubcells(2);
      }
      else {
        if(ambientDim == 2) {
          
          // In 2D "faces" are 1-dimensional subcells (edges)
          numFacesPerCell = this -> getMyNumSubcells(1);
        }
        else {
          
          // Face tags are not meaningful unless dimension of the generating cell is 2 or 3
          TEST_FOR_EXCEPTION( (ambientDim == 1), 
                              std::invalid_argument,
                              ">>> ERROR (MultiCell): Face tags cannot be defined for 1-dimensional cells"); 
        }
      }
      
      // Allocate a zero array of 'short' arrays for face tags:
      Teuchos::Array<short> tempFaceArray(numFacesPerCell);
      faceTags_.assign(numCells_, tempFaceArray);
      
      // Fill in the face signs
      for (int i=0; i<numCells_; i++) {
        faceTags_[i].assign(faceTags+i*numFacesPerCell, faceTags+(i+1)*numFacesPerCell);
      }
    }
  

  
  template<class Scalar>
  void MultiCell<Scalar>::setConnMapCustom(const ECell           customCellType,
                                           const ConnMapTemplate customCellTemplate[]) {
      
      // Verify input: can only write to a custom template
#ifdef INTREPID_HAVE_DEBUG
    TEST_FOR_EXCEPTION( !( (customCellType > CELL_CANONICAL_MAX) && (customCellType < CELL_MAX) ),
                        std::invalid_argument,
                        ">>> ERROR (MultiCell): Invalid custom cell type");
#endif
                      
      // Copy connectivity map of the new cell type into the custom connectivity map storage
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][0] = customCellTemplate[0];
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][1] = customCellTemplate[1];
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][2] = customCellTemplate[2];
  }
  
  
  
  template<class Scalar>
  void MultiCell<Scalar>::setChart(const int cellID) {
    
    // A chart can be set only into an existing atlas.
    TEST_FOR_EXCEPTION( (atlas_.size() == 0),
                        std::invalid_argument,
                        ">>> ERROR (MultiCell): A chart can be set only into an existing atlas. ");

    // Verify that cellID is valid
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                        std::invalid_argument,
                        ">>> ERROR (MultiCell): Invalid cellID value.");
#endif   
    int ambientDim = this -> getMyCellDim();
    
    // Set the default chart of degree 1 for admissible generating cell types.
    switch(myCellType_) {
      
      case CELL_EDGE: {
        Scalar v0 = this -> getCellVertex(cellID,0)[0];
        Scalar v1 = this -> getCellVertex(cellID,1)[0];
        
        atlas_[cellID].refCellType_    = myCellType_;
        atlas_[cellID].mapping_[0][0]  = (v1 - v0)/2.0;
        atlas_[cellID].mapping_[0][1]  = (v1 + v0)/2.0;
        atlas_[cellID].mappingType_    = MAPPING_AFFINE;
        break;
      }
        
      case CELL_TRI:
        atlas_[cellID].refCellType_     = myCellType_;        
        for(int dim=0; dim < ambientDim; dim++){
          Scalar v0 = this -> getCellVertex(cellID,0)[dim];
          Scalar v1 = this -> getCellVertex(cellID,1)[dim];
          Scalar v2 = this -> getCellVertex(cellID,2)[dim];
          
          atlas_[cellID].mapping_[dim][0] = (-v0 + v1);                     // x coefficient
          atlas_[cellID].mapping_[dim][1] = (-v0 + v2);                     // y coefficient
          atlas_[cellID].mapping_[dim][2] =   v0;                           // free term
        }
          atlas_[cellID].mappingType_ = MAPPING_AFFINE;
        break;
        
      case CELL_QUAD:
        atlas_[cellID].refCellType_     = myCellType_;                
        for(int dim=0; dim < ambientDim; dim++){
          Scalar v0 = this -> getCellVertex(cellID,0)[dim];
          Scalar v1 = this -> getCellVertex(cellID,1)[dim];
          Scalar v2 = this -> getCellVertex(cellID,2)[dim];
          Scalar v3 = this -> getCellVertex(cellID,3)[dim];
          
          atlas_[cellID].mapping_[dim][0] = ( v0 - v1 + v2 - v3)/4.0;            // xy coefficient
          atlas_[cellID].mapping_[dim][1] = (-v0 + v1 + v2 - v3)/4.0;            // x  coefficient
          atlas_[cellID].mapping_[dim][2] = (-v0 - v1 + v2 + v3)/4.0;            // y  coefficient
          atlas_[cellID].mapping_[dim][3] = ( v0 + v1 + v2 + v3)/4.0;            // free term
        }
          atlas_[cellID].mappingType_ = MAPPING_NON_AFFINE;
        break;
        
      case CELL_TET:
        atlas_[cellID].refCellType_ = myCellType_;        
        for(int dim=0; dim < ambientDim; dim++){
          Scalar v0 = this -> getCellVertex(cellID,0)[dim];
          Scalar v1 = this -> getCellVertex(cellID,1)[dim];
          Scalar v2 = this -> getCellVertex(cellID,2)[dim];
          Scalar v3 = this -> getCellVertex(cellID,3)[dim];
          
          atlas_[cellID].mapping_[dim][0] = (-v0 + v1);                          // x coefficient
          atlas_[cellID].mapping_[dim][1] = (-v0 + v2);                          // y coefficient
          atlas_[cellID].mapping_[dim][2] = (-v0 + v3);                          // z coefficient
          atlas_[cellID].mapping_[dim][3] = (-v0 + v3);                          // free term
        }
          atlas_[cellID].mappingType_ = MAPPING_AFFINE;
        break;
        
      case CELL_HEX:
        atlas_[cellID].refCellType_ = myCellType_;        
        for(int dim=0; dim < ambientDim; dim++){
          Scalar v0 = this -> getCellVertex(cellID,0)[dim];
          Scalar v1 = this -> getCellVertex(cellID,1)[dim];
          Scalar v2 = this -> getCellVertex(cellID,2)[dim];
          Scalar v3 = this -> getCellVertex(cellID,3)[dim];
          
          Scalar v4 = this -> getCellVertex(cellID,4)[dim];
          Scalar v5 = this -> getCellVertex(cellID,5)[dim];
          Scalar v6 = this -> getCellVertex(cellID,6)[dim];
          Scalar v7 = this -> getCellVertex(cellID,7)[dim];
          
          atlas_[cellID].mapping_[dim][0] = (-v0 + v1 - v2 + v3 + v4 - v5 + v6 - v7)/8.0;   // xyz
          atlas_[cellID].mapping_[dim][1] = ( v0 - v1 + v2 - v3 + v4 - v5 + v6 - v7)/8.0;   // xy
          atlas_[cellID].mapping_[dim][2] = ( v0 - v1 - v2 + v3 - v4 + v5 + v6 - v7)/8.0;   // xz
          atlas_[cellID].mapping_[dim][3] = ( v0 + v1 - v2 - v3 - v4 - v5 + v6 + v7)/8.0;   // yz
          atlas_[cellID].mapping_[dim][4] = (-v0 + v1 + v2 - v3 - v4 + v5 + v6 - v7)/8.0;   // x
          atlas_[cellID].mapping_[dim][5] = (-v0 - v1 + v2 + v3 - v4 - v5 + v6 + v7)/8.0;   // y
          atlas_[cellID].mapping_[dim][6] = (-v0 - v1 - v2 - v3 + v4 + v5 + v6 + v7)/8.0;   // z
          atlas_[cellID].mapping_[dim][7] = ( v0 + v1 + v2 + v3 + v4 + v5 + v6 + v7)/8.0;   // free term
        }
          atlas_[cellID].mappingType_ = MAPPING_NON_AFFINE;
        break;
        
      case CELL_TRIPRISM:
      case CELL_PYRAMID:
        TEST_FOR_EXCEPTION( ( (myCellType_ == CELL_TRIPRISM) || (myCellType_ == CELL_PYRAMID) ),
                            std::invalid_argument,
                            ">>> ERROR (MultiCell): Default chart not yet implemented for this cell type. ");
        break;
        
      default:
        TEST_FOR_EXCEPTION( !( (myCellType_ == CELL_EDGE) ||
                               (myCellType_ == CELL_TRI)  ||
                               (myCellType_ == CELL_QUAD) ||
                               (myCellType_ == CELL_TET)  ||
                               (myCellType_ == CELL_HEX)  ||
                               (myCellType_ == CELL_TRIPRISM)  ||
                               (myCellType_ == CELL_PYRAMID) ),
                            std::invalid_argument,
                            ">>> ERROR (MultiCell): Default chart not available for this cell type. ");
        break;
    } // switch    
  }
  
  
  template<class Scalar>
  void MultiCell<Scalar>::setChart(const int cellID, const ShapePoints<Scalar>& shapePoints) {
    
    // A chart can be set only into an existing atlas.
    TEST_FOR_EXCEPTION( (atlas_.size() == 0),
                        std::invalid_argument,
                        ">>> ERROR (MultiCell): A chart can be set only into an existing atlas. ");
    
    // Verify that cellID is valid
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): Invalid cellID value.");
#endif   
    int ambientDim = this -> getMyCellDim();
                       
    // If size of shapePoints = 0 there are no additional shape points provided, set default chart.
    if( shapePoints.shapePoints_.size() == 0 ) {
      this -> setChart(cellID);      
    } 
    
    // If size of shapePoints != 0 there are additional shape points provided, set higher deg. chart.
    else {
      
      // Temporary: call default chart method:
      this -> setChart(cellID);
      
      // Chart of degree 2
      
      // Chart of degree 3
    }
  }
    
    
  
  template<class Scalar>
  void MultiCell<Scalar>::setAtlas() {    
      
    // Resize atlas to match the number of cells in the MultiCell
    atlas_.resize(numCells_);
      
    // Compute default cell charts of degree 1 and store them in the atlas. 
    for(int cellID = 0; cellID < numCells_; cellID++){
      this -> setChart(cellID);
    }
  }
  
  
  
  template<class Scalar>
    void MultiCell<Scalar>::setAtlas(const Teuchos::Array< ShapePoints<Scalar> >& shapePoints) {
      
      // Resize atlas to match the number of cells
      atlas_.resize(numCells_);
      
      // Compute cell charts using the provided shapePoints set and store in the atlas.
      for(int cellID = 0; cellID < numCells_; cellID++){
        this -> setChart(cellID, shapePoints[cellID]);
      }
    }
  
  


template<class Scalar>
Matrix<Scalar> MultiCell<Scalar>::jacobian(const int cellID, 
                                           const Point<Scalar>& refPoint,
                                           double threshold) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (atlas_.size() == 0),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): jacobian requires an existing atlas. ");  
  TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid cellID value.");
  TEST_FOR_EXCEPTION( ( refPoint.getFrameKind() != FRAME_REFERENCE ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point has to be in FRAME_REFERENCE coordinate frame.");
  TEST_FOR_EXCEPTION( !(inReferenceCell(myCellType_, refPoint, threshold) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point is not inside its reference cell.");
#endif
  
  int ambientDim = this -> getMyCellDim();
  
  // Temp storage for the Matrix coefficients by row
  Scalar DF[9];                                   
  switch(myCellType_) {    
    
    // the EDGE chart is always affine
    case CELL_EDGE:  
      DF[0] = atlas_[cellID].mapping_[0][0];                    
      break;
      
      // TRI and TET charts are affine and can be computed by the same loop
    case CELL_TRI:
    case CELL_TET:
      for(int row = 0; row < ambientDim; row++){
        for(int col = 0; col < ambientDim; col++){
          DF[col + row*ambientDim] = atlas_[cellID].mapping_[row][col]; 
        }
      }
      break;
      
      // For QUAD and HEX rows contain grad of row-th coordinate function evaluated at refPoint
    case CELL_QUAD:                                     
      for(int row = 0; row < ambientDim; row++){
        DF[0 + row*ambientDim] = \
        atlas_[cellID].mapping_[row][0]*refPoint[1] + \
        atlas_[cellID].mapping_[row][1];
        DF[1 + row*ambientDim] = \
          atlas_[cellID].mapping_[row][0]*refPoint[0] + \
          atlas_[cellID].mapping_[row][2];
      }
      break;
      
    case CELL_HEX:
      for(int row = 0; row < ambientDim; row++){
        DF[0 + row*ambientDim] = \
        atlas_[cellID].mapping_[row][0]*refPoint[1]*refPoint[2] + \
        atlas_[cellID].mapping_[row][1]*refPoint[1] + \
        atlas_[cellID].mapping_[row][2]*refPoint[2] + \
        atlas_[cellID].mapping_[row][4];
        //
        DF[1 + row*ambientDim] = \
          atlas_[cellID].mapping_[row][0]*refPoint[0]*refPoint[2] + \
          atlas_[cellID].mapping_[row][1]*refPoint[0] + \
          atlas_[cellID].mapping_[row][3]*refPoint[2] + \
          atlas_[cellID].mapping_[row][5];
        //
        DF[2 + row*ambientDim] = \
          atlas_[cellID].mapping_[row][0]*refPoint[0]*refPoint[1] + \
          atlas_[cellID].mapping_[row][2]*refPoint[0] + \
          atlas_[cellID].mapping_[row][3]*refPoint[1] + \
          atlas_[cellID].mapping_[row][6];
      }
      break;
      
    case CELL_TRIPRISM:
    case CELL_PYRAMID:
      TEST_FOR_EXCEPTION( ( (myCellType_ == CELL_TRIPRISM) || (myCellType_ == CELL_PYRAMID) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): Jacobian not yet implemented for this cell type. ");
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (myCellType_ == CELL_EDGE) ||
                             (myCellType_ == CELL_TRI)  ||
                             (myCellType_ == CELL_QUAD) ||
                             (myCellType_ == CELL_TET)  ||
                             (myCellType_ == CELL_HEX)  ||
                             (myCellType_ == CELL_TRIPRISM)  ||
                             (myCellType_ == CELL_PYRAMID) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): Jacobian not available for this cell type. ");
      break;
  }
  Matrix<Scalar> DF_temp(DF,ambientDim);
  return DF_temp;
}


  
template<class Scalar>
Point<Scalar> MultiCell<Scalar>::mapToPhysicalCell(const int cellID, 
                                                   const Point<Scalar>&  refPoint, 
                                                   const double threshold) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid cellID value.");
  TEST_FOR_EXCEPTION( ( refPoint.getFrameKind() != FRAME_REFERENCE ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point has to be in FRAME_REFERENCE coordinate frame.");
  TEST_FOR_EXCEPTION( !(inReferenceCell(myCellType_, refPoint, threshold) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point is not inside its reference cell.");
#endif
  
  int ambientDim = getMyCellDim();
  
  // Temp array for the image of refPoint
  Scalar physCoords[3];                                    
  switch(myCellType_){
    
    case CELL_EDGE:
      physCoords[0] = \
      atlas_[cellID].mapping_[0][0]*refPoint[0] +\
      atlas_[cellID].mapping_[0][1];
      break;
      
    case CELL_TRI:
      for(int dim = 0;dim < ambientDim; dim++){
        physCoords[dim] = \
        atlas_[cellID].mapping_[dim][0]*refPoint[0] + \
        atlas_[cellID].mapping_[dim][1]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][2];
      }
      break;
      
    case CELL_QUAD:
      for(int dim = 0; dim < ambientDim; dim++){
        physCoords[dim] = \
        atlas_[cellID].mapping_[dim][0]*refPoint[0]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][1]*refPoint[0] + \
        atlas_[cellID].mapping_[dim][2]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][3];
      }
      break;
      
    case CELL_TET:
      for(int dim = 0; dim < ambientDim; dim++){
        physCoords[dim] = \
        atlas_[cellID].mapping_[dim][0]*refPoint[0] + \
        atlas_[cellID].mapping_[dim][1]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][2]*refPoint[2] + \
        atlas_[cellID].mapping_[dim][3];
      }
      break;
      
    case CELL_HEX:
      for(int dim = 0; dim < ambientDim; dim++){
        physCoords[dim] = \
        atlas_[cellID].mapping_[dim][0]*refPoint[0]*refPoint[1]*refPoint[2]+\
        atlas_[cellID].mapping_[dim][1]*refPoint[0]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][2]*refPoint[0]*refPoint[2] + \
        atlas_[cellID].mapping_[dim][3]*refPoint[1]*refPoint[2] + \
        atlas_[cellID].mapping_[dim][4]*refPoint[0] + \
        atlas_[cellID].mapping_[dim][5]*refPoint[1] + \
        atlas_[cellID].mapping_[dim][6]*refPoint[2] + \
        atlas_[cellID].mapping_[dim][7]; 
      }
      break;
      
    case CELL_TRIPRISM:
    case CELL_PYRAMID:
      TEST_FOR_EXCEPTION( ( (myCellType_ == CELL_TRIPRISM) || (myCellType_ == CELL_PYRAMID) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): map to physical cell not yet implemented for this cell type. ");
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (myCellType_ == CELL_EDGE) ||
                             (myCellType_ == CELL_TRI)  ||
                             (myCellType_ == CELL_QUAD) ||
                             (myCellType_ == CELL_TET)  ||
                             (myCellType_ == CELL_HEX)  ||
                             (myCellType_ == CELL_TRIPRISM)  ||
                             (myCellType_ == CELL_PYRAMID) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): map to physical cell not available for this cell type. ");
      break;
  }
  
  // The return point is in PHYSICAL space, set its type accordingly:
  Point<Scalar> physPoint(physCoords,ambientDim,FRAME_PHYSICAL);
  return physPoint;
}

  

template<class Scalar>
Point<Scalar> MultiCell<Scalar>::mapToReferenceCell(const int cellID, 
                                                    const Point<Scalar>& physPoint) const 
{
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid cellID value.");
  TEST_FOR_EXCEPTION( ( physPoint.getFrameKind() != FRAME_PHYSICAL ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point has to be in FRAME_PHYSICAL coordinate frame.");
#endif
  
  // This method cannot work if <var>atlas_<var> is undefined, i.e., its size is 0
  TEST_FOR_EXCEPTION( (atlas_.size() == 0),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): The atlas of this MultiCell has not been defined.");
  
  // Initialize old (xOld=0) and new (refPoint) Newton iterates. Must be FRAME_REFERENCE type
  int ambientDim = getMyCellDim();
  Point<Scalar> xOld(ambientDim,FRAME_REFERENCE);         
  Point<Scalar> refPoint = xOld;                                         
  
  // Newton method to compute the inverse of the mapping between reference and physicall HEX
  for(int iter = 0; iter < INTREPID_MAX_NEWTON; ++iter)	{	
    
    // First iterates may fail the inclusion tests using the tighter INTREPID_THRESHOLD value.
    // Define a dynamic threshold value that decreases as the iteration count grows:
    double threshold = (1.0 + 10.0 * exp(-(double)iter))*INTREPID_TOL;	
    Matrix<Scalar> jacobian_inv(ambientDim);
    
    // Compute Jacobian matrix at old iterate and invert in place. Use dynamic threshold!
//    try{
      jacobian_inv = this -> jacobian(cellID, xOld, threshold);	
      jacobian_inv.invert();
//    }
//    catch (std::invalid_argument err) {
      
      // If after INTREPID_MAX_NEWTON - 10 iterations jacobian still keeps throwing excetions
      // because xOld is not in the reference cell, throw an exception here.
//#ifdef HAVE_INTREPID_DEBUG
//      TEST_FOR_EXCEPTION( (iter > INTREPID_MAX_NEWTON - 10), 
//                          std::invalid_argument,
//                          ">>> ERROR (MultiCell): mapToReferenceCell failed to converge.");
//#endif
//    }
    
    // The Newton step. Use dynamic threshold in mapToPhysicalCell to avoid warnings
    refPoint = xOld + jacobian_inv*(physPoint - this -> mapToPhysicalCell(cellID,xOld,threshold));
    
    // Compute Euclidean distance (error) between old and new iterates: |xOld - refPoint|
    Scalar error = refPoint.distance(xOld);	
    
    // Check convergence and whether the max number of iterations has been exceeded 
    if (error < INTREPID_TOL) {		          
      break;
    } 
    else if ( iter > INTREPID_MAX_NEWTON - 5) {
#ifdef HAVE_INTREPID_DEBUG
      std::cout << " MultiCell::mapToReferenceCell warning: for " << physPoint
      << " method failed to converge to desired tolerance within " << INTREPID_MAX_NEWTON
      << " iterations\n";
#endif
      break;
    }
    
    // initialize next Newton step
    xOld = refPoint;	
  }
  
  // In DEBUG mode warn if the computed Point is not inside its reference cell.
  // Uses the "looser" INTREPID_TOL instead of the default INTREPID_THRESHOLD.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !(inReferenceCell(myCellType_, refPoint,INTREPID_TOL) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Computed reference point is outside its reference cell.");
#endif
  return refPoint;
}



template<class Scalar>
bool MultiCell<Scalar>::inReferenceCell(const ECell cellType, 
                                            const Point<Scalar>& refPoint,
                                            const double threshold) {
  
  // Verify input: Point has to be in FRAME_REFERENCE
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( refPoint.getFrameKind() != FRAME_REFERENCE ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point has to be in FRAME_REFERENCE coordinate frame.");
  TEST_FOR_EXCEPTION( ( refPoint.getDim() != MultiCell<Scalar>::getCellDim(cellType) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point dimension does not match reference cell dimension. ");
  
  
#endif
  bool testResult = true;
  
  // Use <var>threshold<var> to set reference element bounds for the tests
  double minus_one = -1.0 - threshold;
  double plus_one  =  1.0 + threshold;
  double minus_zero = - threshold;
  
  switch(cellType) {
    case CELL_EDGE:
      if( !(minus_one <= refPoint[0] && refPoint[0] <= plus_one)) \
        testResult = false;
      break;
    case CELL_TRI:{
      Scalar distance = std::max( std::max( -refPoint[0], -refPoint[1] ), refPoint[0] + refPoint[1] - 1 );
      if( distance > threshold ) testResult = false;
      break;
    }
    case CELL_QUAD:
      if(!((minus_one <= refPoint[0] && refPoint[0] <= plus_one) && \
           (minus_one <= refPoint[1] && refPoint[1] <= plus_one)))  \
             testResult = false;   
      break;
    case CELL_TET:{
      Scalar distance = std::max(  std::max(-refPoint[0],-refPoint[1]), std::max(-refPoint[2],refPoint[0] + refPoint[1] + refPoint[2] - 1)  );
      if( distance > threshold ) testResult = false;
      break;
    }
    case CELL_HEX:
      if(!((minus_one <= refPoint[0] && refPoint[0] <= plus_one) && \
           (minus_one <= refPoint[1] && refPoint[1] <= plus_one) && \
           (minus_one <= refPoint[2] && refPoint[2] <= plus_one)))  \
             testResult = false;
      break;
    case CELL_TRIPRISM: {
      // The base of the reference prism is the same as the reference triangle. Thus, 
      // the point is inside we apply the tri test to the first two 
      // coordinates and test the third coordinate to be between 
      Scalar distance = std::max( std::max( -refPoint[0], -refPoint[1] ), refPoint[0] + refPoint[1] - 1 );
      if( distance > threshold  || \
          refPoint[2] < minus_zero || refPoint[2] > plus_one) \
            testResult = false;
      break;
    }
    case CELL_PYRAMID:{
      // The base of the reference pyramid is the same as the reference quad cell.
      // Thus, a horizontal plane through a point P(x,y,z) inside the pyramid intersects
      // it at a quad that equals the base quad scaled by (1-z). Therefore, to check if
      // P(x,y,z) is inside the pyramid we check if (x,y) is inside [-1+z,1-z]^2
      Scalar left  = minus_one + refPoint[2];
      Scalar right = plus_one  - refPoint[2];
      if(!((left <= refPoint[0] && refPoint[0] <= right) && \
           (left <= refPoint[1] && refPoint[1] <= right)))  \
             testResult = false;  
      break;
    }
    default:
      TEST_FOR_EXCEPTION( !( (cellType == CELL_EDGE) ||
                             (cellType == CELL_TRI)  ||
                             (cellType == CELL_QUAD) ||
                             (cellType == CELL_TET)  ||
                             (cellType == CELL_HEX)  ||
                             (cellType == CELL_TRI)  ||
                             (cellType == CELL_TRIPRISM)  ||
                             (cellType == CELL_PYRAMID) ),
                          std::invalid_argument,
                          ">>> ERROR (MultiCell): Invalid cell type. ");
  }
  return testResult;
}



template<class Scalar>
bool MultiCell<Scalar>::inPhysicalCell(const int cellID, 
                                           const Point<Scalar>& physPoint,
                                           const double threshold) const 
{
  // Verify input: Point has to be in FRAME_PHYSICAL
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( ( 0 <= cellID) && (cellID < numCells_ ) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid cellID value.");
  TEST_FOR_EXCEPTION( ( physPoint.getFrameKind() != FRAME_PHYSICAL ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point has to be in FRAME_PHYSICAL coordinate frame.");
  
  TEST_FOR_EXCEPTION( ( physPoint.getDim() != MultiCell<Scalar>::getMyCellDim() ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Point dimension does not match generating cell dimension. ");
#endif
  
  // Map the point back to the reference cell and invoke the inclusion test for reference cells
  return inReferenceCell(myCellType_, mapToReferenceCell(cellID, physPoint), threshold);
}


      
  template<class Scalar>
  void MultiCell<Scalar>::printMyInfo(std::ostream & out) const {
      out.setf(std::ios_base::scientific, std::ios_base::floatfield);
      out.precision(6);
      
      EStatus atlasStatus = STATUS_UNDEFINED;
      EStatus edgeSignStatus = STATUS_UNDEFINED;
      EStatus faceSignStatus = STATUS_UNDEFINED;
      EStatus edgeTagsStatus = STATUS_UNDEFINED;
      EStatus faceTagsStatus = STATUS_UNDEFINED;
      
      if(atlas_.size() > 0)     {atlasStatus = STATUS_DEFINED; }
      if(edgeSigns_.size() > 0) {edgeSignStatus = STATUS_DEFINED;}
      if(faceSigns_.size() > 0) {faceSignStatus = STATUS_DEFINED;}
      if(edgeTags_.size() > 0)  {edgeTagsStatus = STATUS_DEFINED;}
      if(faceTags_.size() > 0) {faceTagsStatus = STATUS_DEFINED;}
      
       
      out  << "\n============================= MultiCell info: =============================\n\n"; 
      out << std::left;
      out << std::setw(30) << "\t Generating cell type:" << this -> getMyCellName() << "\n";
      out << std::setw(30) << "\t Atlas status:"         << EStatusToString(atlasStatus)    <<"\n";
      out << std::setw(30) << "\t Edge signs status:"    << EStatusToString(edgeSignStatus) <<"\n";
      out << std::setw(30) << "\t Face signs status:"    << EStatusToString(faceSignStatus) <<"\n";
      out << std::setw(30) << "\t Edge tags status: "    << EStatusToString(edgeTagsStatus) <<"\n";
      out << std::setw(30) << "\t Face tags status: "    << EStatusToString(faceTagsStatus) <<"\n";
      out << std::setw(30) << "\t Number of cells:"      << numCells_ << "\n\n";
      out << "Cell vertices:\n\n";
      
      // Print the vertices of all cells in the MultiCell
      int numNodesPerCell = this -> getMyCellNumNodes();
      for (int i=0; i < numCells_; i++) {
        out << std::setw(16) << "CELL ID = " << i << "\n";
        for (int j=0; j < numNodesPerCell; j++) {
          out << vertices_[i][j] << "\n";
        }
      }
      
      // Print edge connectivities
      out << "Edge template:\n";
      int numEdgesPerCell = this -> getMyNumSubcells(1);
      for (int i=0; i < numEdgesPerCell; i++) {
        Teuchos::Array<int> tempNodes;
        this -> getMySubcellNodeIDs(tempNodes, 1, i);
        out << "  " << std::setw(3) << i << " -> {" <<  tempNodes[0] << ", " << tempNodes[1] << "}" << "\n";
      }
      
      // Print edge signs only if they are defined!
      if(edgeSignStatus ==STATUS_DEFINED) {
        out << "Edge signs:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numEdgesPerCell; j++) {
            out << std::setw(5) << edgeSigns_[i][j];
          }
          out << "\n";
        }
      }
      
      // Print edge tags only if they are defined!
      if(edgeTagsStatus ==STATUS_DEFINED) {
        out << "Edge tags:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numEdgesPerCell; j++) {
            out << std::setw(5) << edgeTags_[i][j];
          }
          out << "\n";
        }
      }
      
      // Print face connectivities
      out << "Face template:\n";
      int numFacesPerCell = this -> getMyNumSubcells(2);
      for (int i=0; i < numFacesPerCell; i++) {
        Teuchos::Array<int> tempNodes;
        ECell tempType = this -> getMySubcellType(2, i);
        this -> getMySubcellNodeIDs(tempNodes, 2, i);
        out << "    " << i << " -> " << std::setw(13) << getCellName(tempType) << ": ";
        for (int j=0; j < this -> getCellNumNodes(tempType); j++) {
          out << std::setw(4) << tempNodes[j];
        }
        out << "\n";
      }
      
      // Print face signs only if they are defined
      if(faceSignStatus == STATUS_DEFINED) {
        out << "Face signs:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numFacesPerCell; j++) {
            out << std::setw(5) << faceSigns_[i][j];
          }
          out << "\n";
        }
      }
      
      // Print face tags only if they are defined
      if(faceTagsStatus == STATUS_DEFINED) {
        out << "Face tags:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numFacesPerCell; j++) {
            out << std::setw(5) << faceTags_[i][j];
          }
          out << "\n";
        }
      }
      out  << "\n=========================== End MultiCell info: ===========================\n"; 
  } // end printMyInfo
  
  
  template<class Scalar>
  std::ostream& operator << (std::ostream& os, const MultiCell<Scalar>& base) {
      base.printMyInfo(os);
      return os;
  }
  
} // namespace Intrepid
