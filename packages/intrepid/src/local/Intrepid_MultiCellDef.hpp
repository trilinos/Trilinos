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
    // initializations
    numCells_(numCells),
    myCellType_(generatingCellType) 
  {
      // Ambient dimension is always equal to the topological dimension of the generating cell
      ambientDim_ = getTopologicalDim(generatingCellType);
      
      // By default, the ctor never computes the charts of the cells
      atlasStatus_ = STATUS_UNDEFINED;
      
      // This constructor does not accept edge/face signs data and their status is undefined
      edgeSignsStatus_ = STATUS_UNDEFINED;                             
      faceSignsStatus_ = STATUS_UNDEFINED;                             
      
      int numNodesPerCell = this -> getMyNumNodes();
      
      // allocate a chunk of memory first, creating a zero array of Point arrays
      Point<Scalar> tempPoint(ambientDim_);
      Teuchos::Array< Point<Scalar> > tempPointArray(numNodesPerCell, tempPoint);
      vertices_.assign(numCells_, tempPointArray);
      
      // fill in the coordinate info
      for (int i=0; i<numCells_; i++) {
        for (int j=0; j<numNodesPerCell; j++) {
          
          // do some pointer arithmetic
          Point<Scalar> workPoint(vertices + (i*numNodesPerCell+j)*ambientDim_, ambientDim_);
          vertices_[i][j] = workPoint;
        }
      }
      
      // allocate a zero array of 'short' arrays for the edge signs
      int numEdgesPerCell = this -> getMyNumSubcells(1);
      Teuchos::Array<short> tempEdgeArray(numEdgesPerCell);
      edgeSigns_.assign(numCells_, tempEdgeArray);
      
      // allocate a zero array of 'short' arrays for the face signs
      int numFacesPerCell = this -> getMyNumSubcells(2);
      Teuchos::Array<short> tempFaceArray(numFacesPerCell);
      faceSigns_.assign(numCells_, tempFaceArray);
  }
  
  
  template<class Scalar>
  MultiCell<Scalar>::MultiCell(const int      numCells,
                               const ECell    generatingCellType,
                               const Scalar*  vertices,
                               const short*   subcellSigns,
                               const int      subcellDim) :
    // intializations
    numCells_(numCells),
    myCellType_(generatingCellType) 
  {
      // Ambient dimension is always equal to the topological dimension of the generating cell
      ambientDim_ = getTopologicalDim(generatingCellType);

      // By default, the ctor never computes the charts of the cells
      atlasStatus_ = STATUS_UNDEFINED;
      
      // This constructor accepts edge OR face signs, depending on what subcellDim is:
      switch(subcellDim){ 
        case 1: {          
          // Allocate space and fill edge orientations
          int numEdgesPerCell = this -> getMyNumSubcells(1);
          Teuchos::Array<short> tempEdgeArray(numEdgesPerCell);
          edgeSigns_.assign(numCells_, tempEdgeArray);
          
          // fill in the orientation info
          for (int i=0; i<numCells_; i++) {
            edgeSigns_[i].assign(subcellSigns+i*numEdgesPerCell, subcellSigns+(i+1)*numEdgesPerCell);
          }
          edgeSignsStatus_ = STATUS_DEFINED; 
          faceSignsStatus_ = STATUS_UNDEFINED; 
          break;
        }
        case 2: {
          // Allocate space and fill face orientations
          int numFacesPerCell = this -> getMyNumSubcells(2);
          Teuchos::Array<short> tempFaceArray(numFacesPerCell);
          faceSigns_.assign(numCells_, tempFaceArray);
          
          // fill in the orientation info
          for (int i=0; i<numCells_; i++) {
            faceSigns_[i].assign(subcellSigns+i*numFacesPerCell, subcellSigns+(i+1)*numFacesPerCell);
          }
          faceSignsStatus_ = STATUS_DEFINED; 
          edgeSignsStatus_ = STATUS_UNDEFINED; 
          break;
        }
        default:
          std::cerr << "MultiCell ctor error: invalid subcell dimension for setting edge or face signs";
          exit(1);
      }
      
      // allocate a chunk of memory first, creating a zero vector of Point vectors
      int numNodesPerCell = this -> getNumNodes(myCellType_);
      Point<Scalar> tempPoint(ambientDim_);
      Teuchos::Array< Point<Scalar> > tempPointArray(numNodesPerCell, tempPoint);
      vertices_.assign(numCells_, tempPointArray);
      
      // fill in the coordinate info
      for (int i=0; i<numCells_; i++) {
        for (int j=0; j<numNodesPerCell; j++) {
          
          // do some pointer arithmetic
          Point<Scalar> workPoint(vertices+(i*numNodesPerCell+j)*ambientDim_, ambientDim_);
          vertices_[i][j] = workPoint;
        }
      }      
  }
  
  
  template<class Scalar>
  MultiCell<Scalar>::MultiCell(const int      numCells,
                               const ECell    generatingCellType,
                               const Scalar*  vertices,
                               const short*   edgeSigns,
                               const short*   faceSigns) :
    // intializations
    numCells_(numCells),
    myCellType_(generatingCellType) 
  {
      // Ambient dimension is always equal to the topological dimension of the generating cell
      ambientDim_ = getTopologicalDim(generatingCellType);

      // By default, the ctor never computes the charts of the cells
      atlasStatus_ = STATUS_UNDEFINED;
      
      // allocate a chunk of memory first, creating a zero vector of Point vectors
      int numNodesPerCell = this -> getNumNodes(myCellType_);
      Point<Scalar> tempPoint(ambientDim_);
      Teuchos::Array< Point<Scalar> > tempPointArray(numNodesPerCell, tempPoint);
      vertices_.assign(numCells_, tempPointArray);
      
      // fill in the coordinate info
      for (int i=0; i<numCells_; i++) {
        for (int j=0; j<numNodesPerCell; j++) {
          
          // do some pointer arithmetic
          Point<Scalar> workPoint(vertices+(i*numNodesPerCell+j)*ambientDim_, ambientDim_);
          vertices_[i][j] = workPoint;
        }
      }
      
      // This constructor accepts edge and face signs: Allocate space and fill edge orientations
      int numEdgesPerCell = this -> getMyNumSubcells(1);
      Teuchos::Array<short> tempEdgeArray(numEdgesPerCell);
      edgeSigns_.assign(numCells_, tempEdgeArray);
      
      // fill in the orientation info
      for (int i=0; i<numCells_; i++) {
        edgeSigns_[i].assign(edgeSigns+i*numEdgesPerCell, edgeSigns+(i+1)*numEdgesPerCell);
      }
      edgeSignsStatus_ = STATUS_DEFINED;                             
      
      
      // Allocate space and fill face orientations
      int numFacesPerCell = this -> getMyNumSubcells(2);
      Teuchos::Array<short> tempFaceArray(numFacesPerCell);
      faceSigns_.assign(numCells_, tempFaceArray);
      
      // fill in the orientation info
      for (int i=0; i<numCells_; i++) {
        faceSigns_[i].assign(faceSigns+i*numFacesPerCell, faceSigns+(i+1)*numFacesPerCell);
      }
      faceSignsStatus_ = STATUS_DEFINED;                             
  }
  

  template<class Scalar>
  void MultiCell<Scalar>::setConnMapCustom(const ECell           customCellType,
                                           const ConnMapTemplate customCellTemplate[]) {
      
      // Can only write to a custom template
      assert((customCellType > CELL_CANONICAL_MAX) && (customCellType < CELL_MAX)); 
      
      // define connectivity map of the new cell type
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][0] = customCellTemplate[0];
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][1] = customCellTemplate[1];
      connMapCustom_[customCellType - CELL_CANONICAL_MAX - 1][2] = customCellTemplate[2];
  }

  
  template<class Scalar>
  EFailCode MultiCell<Scalar>::insideReferenceCell(const ECell cellType, 
                                                   const Point<Scalar>& refPoint,
                                                   const double threshold) {
      
      // First make sure the Point is in FRAME_REFERENCE
      if(refPoint.getFrameKind() != FRAME_REFERENCE) {
        std::cerr<<"MultiCell::insideReferenceCell error: \n"
        "\t Point frame is not of FRAME_REFERENCE kind\n";
        exit(1);
      }
      
      // Then, check if Point and cell dimensions coincide
      if(refPoint.getDim() != MultiCell<Scalar>::getTopologicalDim(cellType) ) {
        std::cerr<<"MultiCell::insideReferenceCell error: \n"
        "\t Point dimension is not the same as the cell dimension\n";
        exit(1);
      }
      EFailCode code = FAIL_CODE_SUCCESS;
      
      // Use <var>threshold<var> to set reference element bounds for the tests
      double minus_one = -1.0 - threshold;
      double plus_one  =  1.0 + threshold;
      double minus_zero = - threshold;
      
      switch(cellType) {
        case CELL_EDGE:
          if( !(minus_one <= refPoint[0] && refPoint[0] <= plus_one)) \
            code = FAIL_CODE_NOT_IN_REF_CELL;
          break;
        case CELL_TRI:{
          Scalar distance = std::max( std::max( -refPoint[0], -refPoint[1] ), refPoint[0] + refPoint[1] - 1 );
          if( distance > threshold ) code = FAIL_CODE_NOT_IN_REF_CELL;
          break;
        }
        case CELL_QUAD:
          if(!((minus_one <= refPoint[0] && refPoint[0] <= plus_one) && \
               (minus_one <= refPoint[1] && refPoint[1] <= plus_one)))  \
                 code = FAIL_CODE_NOT_IN_REF_CELL;   
          break;
        case CELL_TET:{
          Scalar distance = std::max(  std::max(-refPoint[0],-refPoint[1]), std::max(-refPoint[2],refPoint[0] + refPoint[1] + refPoint[2] - 1)  );
          if( distance > threshold ) code = FAIL_CODE_NOT_IN_REF_CELL;
          break;
        }
        case CELL_HEX:
          if(!((minus_one <= refPoint[0] && refPoint[0] <= plus_one) && \
               (minus_one <= refPoint[1] && refPoint[1] <= plus_one) && \
               (minus_one <= refPoint[2] && refPoint[2] <= plus_one)))  \
                 code = FAIL_CODE_NOT_IN_REF_CELL;
          break;
        case CELL_TRIPRISM: {
          // The base of the reference prism is the same as the reference triangle. Thus,
          // to check if the point is inside we apply the tri test to the first two 
          // coordinates and test the third coordinate to be between 
          Scalar distance = std::max( std::max( -refPoint[0], -refPoint[1] ), refPoint[0] + refPoint[1] - 1 );
          if( distance > threshold  || \
              refPoint[2] < minus_zero || refPoint[2] > plus_one) \
                code = FAIL_CODE_NOT_IN_REF_CELL;
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
                 code = FAIL_CODE_NOT_IN_REF_CELL;  
          break;
        }
        default:
          std::cerr<< "MultiCell::inclusion error: invalid cell type \n";
          exit(1);
      }
      
      return code;
  }
  

  template<class Scalar>
  void MultiCell<Scalar>::setAtlas() {
      
      // First check if <var>atlas_<var> was already set
      if(atlasStatus_ == STATUS_DEFINED) return; 
      
      // make sure there's enough room to hold all charts to avoid resizing of the array
      atlas_.resize(numCells_);
      
      // Loop over all cells in the MultiCell, compute their charts and store them in the <var>atlas_<var>
      for(int cellID = 0; cellID < numCells_; cellID++){
        switch(myCellType_) {
          case CELL_EDGE: {
            Scalar v0 = this -> getVertex(cellID,0)[0];
            Scalar v1 = this -> getVertex(cellID,1)[0];
            
            atlas_[cellID].refCellType_    = myCellType_;
            atlas_[cellID].mapping_[0][0]  = (v1 - v0)/2.0;
            atlas_[cellID].mapping_[0][1]  = (v1 + v0)/2.0;
            atlas_[cellID].mappingType_    = MAPPING_AFFINE;
            break;
          }
          case CELL_TRI:
            atlas_[cellID].refCellType_     = myCellType_;        
            for(int dim=0; dim < ambientDim_; dim++){
              Scalar v0 = this -> getVertex(cellID,0)[dim];
              Scalar v1 = this -> getVertex(cellID,1)[dim];
              Scalar v2 = this -> getVertex(cellID,2)[dim];
              
              atlas_[cellID].mapping_[dim][0] = (-v0 + v1);                     // x coefficient
              atlas_[cellID].mapping_[dim][1] = (-v0 + v2);                     // y coefficient
              atlas_[cellID].mapping_[dim][2] =   v0;                           // free term
            }
            atlas_[cellID].mappingType_ = MAPPING_AFFINE;
            break;
          case CELL_QUAD:
            atlas_[cellID].refCellType_     = myCellType_;                
            for(int dim=0; dim < ambientDim_; dim++){
              Scalar v0 = this -> getVertex(cellID,0)[dim];
              Scalar v1 = this -> getVertex(cellID,1)[dim];
              Scalar v2 = this -> getVertex(cellID,2)[dim];
              Scalar v3 = this -> getVertex(cellID,3)[dim];
              
              atlas_[cellID].mapping_[dim][0] = ( v0 - v1 + v2 - v3)/4.0;            // xy coefficient
              atlas_[cellID].mapping_[dim][1] = (-v0 + v1 + v2 - v3)/4.0;            // x  coefficient
              atlas_[cellID].mapping_[dim][2] = (-v0 - v1 + v2 + v3)/4.0;            // y  coefficient
              atlas_[cellID].mapping_[dim][3] = ( v0 + v1 + v2 + v3)/4.0;            // free term
            }
            atlas_[cellID].mappingType_ = MAPPING_NON_AFFINE;
            break;
          case CELL_TET:
            atlas_[cellID].refCellType_ = myCellType_;        
            for(int dim=0; dim < ambientDim_; dim++){
              Scalar v0 = this -> getVertex(cellID,0)[dim];
              Scalar v1 = this -> getVertex(cellID,1)[dim];
              Scalar v2 = this -> getVertex(cellID,2)[dim];
              Scalar v3 = this -> getVertex(cellID,3)[dim];
              
              atlas_[cellID].mapping_[dim][0] = (-v0 + v1);                          // x coefficient
              atlas_[cellID].mapping_[dim][1] = (-v0 + v2);                          // y coefficient
              atlas_[cellID].mapping_[dim][2] = (-v0 + v3);                          // z coefficient
              atlas_[cellID].mapping_[dim][3] = (-v0 + v3);                          // free term
            }
            atlas_[cellID].mappingType_ = MAPPING_AFFINE;
            break;
          case CELL_HEX:
            atlas_[cellID].refCellType_ = myCellType_;        
            for(int dim=0; dim < ambientDim_; dim++){
              Scalar v0 = this -> getVertex(cellID,0)[dim];
              Scalar v1 = this -> getVertex(cellID,1)[dim];
              Scalar v2 = this -> getVertex(cellID,2)[dim];
              Scalar v3 = this -> getVertex(cellID,3)[dim];
              
              Scalar v4 = this -> getVertex(cellID,4)[dim];
              Scalar v5 = this -> getVertex(cellID,5)[dim];
              Scalar v6 = this -> getVertex(cellID,6)[dim];
              Scalar v7 = this -> getVertex(cellID,7)[dim];
              
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
            std::cout << " Pullback method not implemented...\n";
            exit(1);
            break;
          default:
            std::cerr<<"MultiCell::getAtlas error: atlas_ not available for this cell type \n";
            exit(1);
            break;
        }
      }
      atlasStatus_ = STATUS_DEFINED;
  }
  
  template<class Scalar>
  const ChartTemplate<Scalar>& MultiCell<Scalar>::getChart(const int cellID) const{
      if(atlasStatus_ == STATUS_DEFINED) {
        return atlas_[cellID]; 
      }
      else {
        std::cerr<<"MultiCell::getChart error: \n"
        "\t cell Atlas has not been defined\n";
        exit(1);
      }
  }

  
  template<class Scalar>
  Matrix<Scalar> MultiCell<Scalar>::jacobian(const int cellID, 
                                             const Point<Scalar>& refPoint,
                                             double threshold) const{

      // Check if the PointType is in the FRAME_REFERENCE frame
      if(refPoint.getFrameKind() == FRAME_PHYSICAL){
        std::cerr << "MultiCell::Jacobian warning: \n"
        "\t Admissible frame type for refPoint is FRAME_REFERENCE, but refPoint has FRAME_PHYSICAL /n";
      }
      
      // Check if refPoint is indeed inside its reference cell - disable for efficiency.
      if(insideReferenceCell(myCellType_, refPoint, threshold) != FAIL_CODE_SUCCESS) {
        std::cerr << "MultiCell::Jacobian warning: \n"
        "\t refPoint has FRAME_REFERENCE type but is not inside its reference cell \n";
      }
      
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
          for(int row = 0; row < ambientDim_; row++){
            for(int col = 0; col < ambientDim_; col++){
              DF[col + row*ambientDim_] = atlas_[cellID].mapping_[row][col]; 
            }
          }
          break;
          
        // For QUAD and HEX rows contain grad of row-th coordinate function evaluated at refPoint
        case CELL_QUAD:                                     
          for(int row = 0; row < ambientDim_; row++){
            DF[0 + row*ambientDim_] = \
              atlas_[cellID].mapping_[row][0]*refPoint[1] + \
              atlas_[cellID].mapping_[row][1];
            DF[1 + row*ambientDim_] = \
              atlas_[cellID].mapping_[row][0]*refPoint[0] + \
              atlas_[cellID].mapping_[row][2];
          }
          break;
        case CELL_HEX:
          for(int row = 0; row < ambientDim_; row++){
            DF[0 + row*ambientDim_] = \
              atlas_[cellID].mapping_[row][0]*refPoint[1]*refPoint[2] + \
              atlas_[cellID].mapping_[row][1]*refPoint[1] + \
              atlas_[cellID].mapping_[row][2]*refPoint[2] + \
              atlas_[cellID].mapping_[row][4];
            //
            DF[1 + row*ambientDim_] = \
              atlas_[cellID].mapping_[row][0]*refPoint[0]*refPoint[2] + \
              atlas_[cellID].mapping_[row][1]*refPoint[0] + \
              atlas_[cellID].mapping_[row][3]*refPoint[2] + \
              atlas_[cellID].mapping_[row][5];
            //
            DF[2 + row*ambientDim_] = \
              atlas_[cellID].mapping_[row][0]*refPoint[0]*refPoint[1] + \
              atlas_[cellID].mapping_[row][2]*refPoint[0] + \
              atlas_[cellID].mapping_[row][3]*refPoint[1] + \
              atlas_[cellID].mapping_[row][6];
          }
          break;
        default:
          std::cerr << "MultiCell::Jacobian error: not implemented\n";
          exit(1);
      }
      Matrix<Scalar> DF_temp(DF,ambientDim_);
      return DF_temp;
  }
  
  
  template<class Scalar>
  Point<Scalar> MultiCell<Scalar>::mapToPhysicalCell(const int cellID, 
                                                       const Point<Scalar>&  refPoint,
                                                       const double threshold) const {

      // Check if refPoint has the correct frame kind - disable for efficiency
      if(refPoint.getFrameKind() != FRAME_REFERENCE){
        std::cerr <<  " Multicell::mapToPhysicalCell warning: \n"
        "\t Admissible Point type is REFERENCE but argument has PHYSICAL type.\n";
      }
      
      // Check if refPoint is inside its reference cell - disable for efficiency.
      if(insideReferenceCell(myCellType_, refPoint,threshold) != FAIL_CODE_SUCCESS) {
        std::cerr << "MultiCell::mapToPhysicalCell warning: \n"
        "\t refPoint has FRAME_REFERENCE type but is not inside its reference cell \n";
      }
      
      // Temp array for the image of refPoint
      Scalar physCoords[3];                                    
      switch(myCellType_){
        case CELL_EDGE:
          physCoords[0] = \
            atlas_[cellID].mapping_[0][0]*refPoint[0] +\
            atlas_[cellID].mapping_[0][1];
          break;
        case CELL_TRI:
          for(int dim = 0;dim < ambientDim_; dim++){
            physCoords[dim] = \
              atlas_[cellID].mapping_[dim][0]*refPoint[0] + \
              atlas_[cellID].mapping_[dim][1]*refPoint[1] + \
              atlas_[cellID].mapping_[dim][2];
          }
          break;
        case CELL_QUAD:
          for(int dim = 0; dim < ambientDim_; dim++){
            physCoords[dim] = \
              atlas_[cellID].mapping_[dim][0]*refPoint[0]*refPoint[1] + \
              atlas_[cellID].mapping_[dim][1]*refPoint[0] + \
              atlas_[cellID].mapping_[dim][2]*refPoint[1] + \
              atlas_[cellID].mapping_[dim][3];
          }
          break;
        case CELL_TET:
          for(int dim = 0; dim < ambientDim_; dim++){
            physCoords[dim] = \
              atlas_[cellID].mapping_[dim][0]*refPoint[0] + \
              atlas_[cellID].mapping_[dim][1]*refPoint[1] + \
              atlas_[cellID].mapping_[dim][2]*refPoint[2] + \
              atlas_[cellID].mapping_[dim][3];
          }
          break;
        case CELL_HEX:
          for(int dim = 0; dim < ambientDim_; dim++){
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
          std::cout<<"Not Implemented \n";
          exit(1);
          break;
        default:
          std::cerr<<"MultiCell::mapToPhysicallCell error: unexpected cell type\n";
          exit(1);
      }
      
      // The return point is in PHYSICAL space, set its type accordingly:
      Point<Scalar> physPoint(physCoords,ambientDim_,FRAME_PHYSICAL);
      return physPoint;
  }
  
  
  template<class Scalar>
  Point<Scalar> MultiCell<Scalar>::mapToReferenceCell(const int cellID,const Point<Scalar>& physPoint) const {
      
      // Warn if the coordinate frame of physPoint argument is not FRAME_PHYSICAL:
      if(physPoint.getFrameKind() != FRAME_PHYSICAL){
        std::cerr <<  " MultiCell::mapToReferenceCell warning: \n"
        "\t Admissible Point frame is FRAME_PHYSICAL, argument has FRAME_REFERENCE frame.\n";
      }
      
      // This method cannot work if <var>atlas_<var> is undefined. Check its status.
      if(atlasStatus_ == STATUS_UNDEFINED) {
        std::cerr << " MultiCell::mapToReferenceCell error: \n"
        "\t The atlas of this MultiCell is undefined.\n";
        exit(1);
      }
      
      // Initialize old (xOld=0) and new (refPoint) Newton iterates. Must be FRAME_REFERENCE type
      Point<Scalar> xOld(ambientDim_,FRAME_REFERENCE);         
      Point<Scalar> refPoint = xOld;                                         
      
      // Newton method to compute the inverse of the mapping between reference and physicall HEX
      for(int iter = 0; iter < INTREPID_MAX_NEWTON; ++iter)	{	
        
        // First iterates may fail the inclusion tests using the tighter INTREPID_THRESHOLD value.
        // Define a dynamic threshold value that decreases as the itretaion count grows:
        double threshold = std::max( (1.0 + 0.5 * exp(-(double)iter)), INTREPID_THRESHOLD);	
        
        // Compute Jacobian matrix at old iterate and invert in place. Use dynamic threshold!
        Matrix<Scalar> jacobian_inv = this -> jacobian(cellID, xOld, threshold);	
        jacobian_inv.invert();
        
        // The Newton step. Use dynamic threshold in mapToPhysicalCell to avoid warnings
        refPoint = xOld + jacobian_inv*(physPoint - this -> mapToPhysicalCell(cellID,xOld,threshold));
        
        // Compute Euclidean distance (error) between old and new iterates: |xOld - refPoint|
        Scalar error = refPoint.distance(xOld);	
        
        // Check convergence and whether the max number of iterations has been exceeded 
        if (error < INTREPID_TOL) {		          
          break;
        } 
        else if ( iter > INTREPID_MAX_NEWTON - 5) {
          std::cerr << " MultiCell::mapToReferenceCell warning: for " << physPoint
          << " method failed to converge to desired tolerance within " << INTREPID_MAX_NEWTON
          << " iterations\n";
          break;
        }
          
          // initialize next Newton step
          xOld = refPoint;	
      }
      
      // Warn if the computed Point is not inside its reference cell. Disable for efficiency.
      // Uses the "looser" INTREPID_TOL instead of the default INTREPID_THRESHOLD.
      if( insideReferenceCell(myCellType_, refPoint,INTREPID_TOL) == FAIL_CODE_ERROR) {
        std::cerr<< " MultiCell::mapToReferenceCell warning: \n"
        " \t Computed reference point is outside its reference cell\n";
      }
      return refPoint;
  }


  template<class Scalar>
  EFailCode MultiCell<Scalar>::insidePhysicalCell(const int cellID,
                                                  const Point<Scalar>& physPoint,
                                                  const double threshold) const {
      return insideReferenceCell(myCellType_, mapToReferenceCell(cellID, physPoint), threshold);
  }

      
  template<class Scalar>
  void MultiCell<Scalar>::printMyInfo(std::ostream & out) const {
      out.setf(std::ios_base::scientific, std::ios_base::floatfield);
      out.precision(6);
      
      out  << "\n============================= MultiCell info: =============================\n\n"; 
      out << std::left;
      out << std::setw(30) << "\t Generating cell type:" << this -> getMyCellName() << "\n";
      out << std::setw(30) << "\t Atlas status:"         << EStatusToString(atlasStatus_)    <<"\n";
      out << std::setw(30) << "\t Edge signs status:"    << EStatusToString(edgeSignsStatus_) <<"\n";
      out << std::setw(30) << "\t Face signs status:"    << EStatusToString(faceSignsStatus_) <<"\n";
      out << std::setw(30) << "\t Number of cells:"      << numCells_ << "\n\n";
      out << "Cell vertices:\n\n";
      
      // Print the vertices of all cells in the MultiCell
      int numNodesPerCell = this -> getMyNumNodes();
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
        this -> getMySubcellNodeIDs(1, i, tempNodes);
        out << "  " << std::setw(3) << i << " -> {" <<  tempNodes[0] << ", " << tempNodes[1] << "}" << "\n";
      }
      
      // Print edge signs only if they are defined!
      if(edgeSignsStatus_ ==STATUS_DEFINED) {
        out << "Edge signs:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numEdgesPerCell; j++) {
            out << std::setw(5) << edgeSigns_[i][j];
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
        this -> getMySubcellNodeIDs(2, i, tempNodes);
        out << "    " << i << " -> " << std::setw(13) << getCellName(tempType) << ": ";
        for (int j=0; j < this -> getNumNodes(tempType); j++) {
          out << std::setw(4) << tempNodes[j];
        }
        out << "\n";
      }
      
      // Print face signs only if they are defined
      if(faceSignsStatus_ == STATUS_DEFINED) {
        out << "Face signs:\n";
        for (int i=0; i < numCells_; i++) {
          for (int j=0; j < numFacesPerCell; j++) {
            out << std::setw(5) << faceSigns_[i][j];
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
