#ifndef INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_QUAD_C1_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on QUAD cells.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::Basis_HGRAD_QUAD_C1_FEM()
  {
    this -> basisCardinality_  = 4;
    this -> basisDegree_       = 1;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
    
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = { 0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1};
  
  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                             const ArrayScalar &  inputPoints,
                                                             const EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
  
  // Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);  
  
  // Temporaries: (x,y) coordinates of the evaluation point
  Scalar x = 0.0;                                    
  Scalar y = 0.0;                                    
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = (1.0 - x)*(1.0 - y)/4.0;
        outputValues(1, i0) = (1.0 + x)*(1.0 - y)/4.0;
        outputValues(2, i0) = (1.0 + x)*(1.0 + y)/4.0;
        outputValues(3, i0) = (1.0 - x)*(1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = -(1.0 - y)/4.0;
        outputValues(0, i0, 1) = -(1.0 - x)/4.0;
        
        outputValues(1, i0, 0) =  (1.0 - y)/4.0;
        outputValues(1, i0, 1) = -(1.0 + x)/4.0;
        
        outputValues(2, i0, 0) =  (1.0 + y)/4.0;
        outputValues(2, i0, 1) =  (1.0 + x)/4.0;
 
        outputValues(3, i0, 0) = -(1.0 + y)/4.0;
        outputValues(3, i0, 1) =  (1.0 - x)/4.0;
      }
      break;
      
    case OPERATOR_CURL:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = -(1.0 - x)/4.0;
        outputValues(0, i0, 1) =  (1.0 - y)/4.0;
        
        outputValues(1, i0, 0) = -(1.0 + x)/4.0;
        outputValues(1, i0, 1) = -(1.0 - y)/4.0;
        
        outputValues(2, i0, 0) =  (1.0 + x)/4.0;
        outputValues(2, i0, 1) = -(1.0 + y)/4.0;
        
        outputValues(3, i0, 0) =  (1.0 - x)/4.0;
        outputValues(3, i0, 1) =  (1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): DIV is invalid operator for rank-0 (scalar) functions in 2D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality=3) 
        outputValues(0, i0, 0) =  0.0;
        outputValues(0, i0, 1) =  0.25;
        outputValues(0, i0, 2) =  0.0;

        outputValues(1, i0, 0) =  0.0;
        outputValues(1, i0, 1) = -0.25;
        outputValues(1, i0, 2) =  0.0;
        
        outputValues(2, i0, 0) =  0.0;
        outputValues(2, i0, 1) =  0.25;
        outputValues(2, i0, 2) =  0.0;
        
        outputValues(3, i0, 0) =  0.0;
        outputValues(3, i0, 1) = -0.25;
        outputValues(3, i0, 2) =  0.0;
      }
      break;
      
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid::getDkCardinality(operatorType, 
                                                       this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
      }
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): FEM Basis calling an FVD member function");
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID_DEBUG
  // Verify rank of output array.
  TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) rank = 2 required for DofCoords array");
  // Verify 0th dimension of output array.
  TEST_FOR_EXCEPTION( !( DofCoords.dimension(0) == this -> basisCardinality_ ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
  // Verify 1st dimension of output array.
  TEST_FOR_EXCEPTION( !( DofCoords.dimension(1) == (int)(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

  DofCoords(0,0) = -1.0;   DofCoords(0,1) = -1.0;
  DofCoords(1,0) =  1.0;   DofCoords(1,1) = -1.0;
  DofCoords(2,0) =  1.0;   DofCoords(2,1) =  1.0;
  DofCoords(3,0) = -1.0;   DofCoords(3,1) =  1.0;
}

}// namespace Intrepid
#endif
