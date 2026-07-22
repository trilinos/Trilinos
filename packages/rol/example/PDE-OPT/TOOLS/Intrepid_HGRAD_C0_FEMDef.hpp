// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __INTREPID_HGRAD_C0_FEM_DEF_HPP__
#define __INTREPID_HGRAD_C0_FEM_DEF_HPP__

// Define a piecewise constant basis function for quads.

namespace Intrepid{
  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_C0_FEM<Scalar,ArrayScalar>::Basis_HGRAD_C0_FEM()
  {
    this -> basisCardinality_ = 1;
    this -> basisDegree_ = 0;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_ = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_ = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_ = false;
  }
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_C0_FEM<Scalar,ArrayScalar>::initializeTags(){
    int tagSize = 4;   // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;  // position in the tag, counting from 0, of the subcell dim
    int posScOrd = 1;  // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;  // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
    int tags[4] = {2, 0, 0, 1};
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
void Basis_HGRAD_C0_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                       const ArrayScalar &  inputPoints,
                                                       const EOperator      operatorType) const {

  #ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
  #endif
  //Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);
  switch (operatorType){
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++){
            outputValues(0,i0) = 1.0;
      }
      break;

    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = 0.0;
        outputValues(0, i0, 1) = 0.0;
      }
      break;

    case OPERATOR_CURL:
    case OPERATOR_DIV:
    case OPERATOR_D2:
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
    TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                      ">>> ERROR (Basis_HGRAD_C0_FEM): Invalid operator type");
  }
}

template<class Scalar, class ArrayScalar>
void Basis_HGRAD_C0_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_C0_FEM): FEM Basis calling an FVD member function");
}


template<class Scalar, class ArrayScalar>
void Basis_HGRAD_C0_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const{
  //One degree of freedom in the center of the cell.
  DofCoords(0,0) = 0.0; DofCoords(0,1) = 0.0; // specific to quads
}


}
#endif
