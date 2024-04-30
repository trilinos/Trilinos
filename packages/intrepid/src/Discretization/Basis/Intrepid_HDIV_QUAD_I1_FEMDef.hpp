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

/** \file   Intrepid_HDIV_QUAD_I1_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(div) functions on HEX cells.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
*/

namespace Intrepid {

template<class Scalar, class ArrayScalar>
Basis_HDIV_QUAD_I1_FEM<Scalar,ArrayScalar>::Basis_HDIV_QUAD_I1_FEM()
  {
    this -> basisCardinality_  = 4;
    this -> basisDegree_       = 1;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
template<class Scalar, class ArrayScalar>
void Basis_HDIV_QUAD_I1_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = {
                  1, 0, 0, 1,
                  1, 1, 0, 1,
                  1, 2, 0, 1,
                  1, 3, 0, 1
                };
  
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
void Basis_HDIV_QUAD_I1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                           const ArrayScalar &  inputPoints,
                                                           const EOperator      operatorType) const {
  
// Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HDIV_Args<Scalar, ArrayScalar>(outputValues,
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
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = 0.0;
        outputValues(0, i0, 1) = (y - 1.0)/4.0;

        outputValues(1, i0, 0) = (1.0 + x)/4.0;
        outputValues(1, i0, 1) = 0.0; 

        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = (1.0 + y)/4.0;

        outputValues(3, i0, 0) = (x - 1.0)/4.0;
        outputValues(3, i0, 1) = 0.0;
      }
      break;

    case OPERATOR_DIV:
      // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
      for (int i0 = 0; i0 < dim0; i0++) {
         outputValues(0, i0) = 0.25;
         outputValues(1, i0) = 0.25;
         outputValues(2, i0) = 0.25;
         outputValues(3, i0) = 0.25;
       }
      break;

    case OPERATOR_CURL:
       TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HDIV_QUAD_I1_FEM): CURL is invalid operator for HDIV Basis Functions");
      break;
      
    case OPERATOR_GRAD:
       TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_GRAD), std::invalid_argument,
                          ">>> ERROR (Basis_HDIV_QUAD_I1_FEM): GRAD is invalid operator for HDIV Basis Functions");
      break;

    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType == OPERATOR_D1)    ||
                            (operatorType == OPERATOR_D2)    ||
                            (operatorType == OPERATOR_D3)    ||
                            (operatorType == OPERATOR_D4)    ||
                            (operatorType == OPERATOR_D5)    ||
                            (operatorType == OPERATOR_D6)    ||
                            (operatorType == OPERATOR_D7)    ||
                            (operatorType == OPERATOR_D8)    ||
                            (operatorType == OPERATOR_D9)    ||
                            (operatorType == OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (Basis_HDIV_QUAD_I1_FEM): Invalid operator type");
      break;
      
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  &&
                            (operatorType != OPERATOR_CURL)  &&
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    &&
                            (operatorType != OPERATOR_D2)    &&
                            (operatorType != OPERATOR_D3)    &&
                            (operatorType != OPERATOR_D4)    &&
                            (operatorType != OPERATOR_D5)    &&
                            (operatorType != OPERATOR_D6)    &&
                            (operatorType != OPERATOR_D7)    &&
                            (operatorType != OPERATOR_D8)    &&
                            (operatorType != OPERATOR_D9)    &&
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (Basis_HDIV_QUAD_I1_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HDIV_QUAD_I1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                           const ArrayScalar &    inputPoints,
                                                           const ArrayScalar &    cellVertices,
                                                           const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HDIV_QUAD_I1_FEM): FEM Basis calling an FVD member function");
}

}// namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

