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

/** \file   Intrepid_BasisDef.hpp
    \brief  Implementation file for the abstract base class Intrepid::Basis.
    \author Created by P. Bochev and D. Ridzal.
*/


template<class Scalar, class ArrayScalar>
int Basis<Scalar, ArrayScalar>::getDofOrdinal(const int subcDim,
                                              const int subcOrd,
                                              const int subcDofOrd) {
  if (!basisTagsAreSet_) {
    initializeTags();
    basisTagsAreSet_ = true;
  }
  // Use .at() for bounds checking
  int dofOrdinal = tagToOrdinal_.at(subcDim).at(subcOrd).at(subcDofOrd);
  TEST_FOR_EXCEPTION( (dofOrdinal == -1), std::invalid_argument, 
                      ">>> ERROR (Basis): Invalid DoF tag");
  return dofOrdinal;
}


template<class Scalar,class ArrayScalar>
const std::vector<std::vector<std::vector<int> > > & Basis<Scalar, ArrayScalar>::getDofOrdinalData( ) 
{
  if (!basisTagsAreSet_) {
    initializeTags();
    basisTagsAreSet_ = true;
  }
  return tagToOrdinal_;
}


template<class Scalar, class ArrayScalar>
const std::vector<int>&  Basis<Scalar, ArrayScalar>::getDofTag(int dofOrd) {
  if (!basisTagsAreSet_) {
    initializeTags();
    basisTagsAreSet_ = true;
  }
  // Use .at() for bounds checking
  return ordinalToTag_.at(dofOrd);
}


template<class Scalar, class ArrayScalar>
const std::vector<std::vector<int> > & Basis<Scalar, ArrayScalar>::getAllDofTags() {
  if (!basisTagsAreSet_) {
    initializeTags();
    basisTagsAreSet_ = true;
  }
  return ordinalToTag_;
}



template<class Scalar, class ArrayScalar>
inline int Basis<Scalar, ArrayScalar>::getCardinality() const {
  return basisCardinality_;   
}


template<class Scalar, class ArrayScalar>
inline EBasis Basis<Scalar, ArrayScalar>::getBasisType() const {
  return basisType_;
}


template<class Scalar, class ArrayScalar>
inline const shards::CellTopology Basis<Scalar, ArrayScalar>::getBaseCellTopology() const {
  return basisCellTopology_;
}


template<class Scalar, class ArrayScalar>
inline int Basis<Scalar,ArrayScalar>::getDegree() const {
  return basisDegree_;
}


template<class Scalar, class ArrayScalar>
inline ECoordinates Basis<Scalar, ArrayScalar>::getCoordinateSystem() const {
  return basisCoordinates_;
}
  

//--------------------------------------------------------------------------------------------//
//                                                                                            //            
//                            Helper functions of the Basis class                             //
//                                                                                            //            
//--------------------------------------------------------------------------------------------//

template<class Scalar, class ArrayScalar>
void getValues_HGRAD_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard){
  
  int spaceDim = cellTopo.getDimension();
  
  // Verify inputPoints array
  TEST_FOR_EXCEPTION( !(inputPoints.rank() == 2), std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 2 required for inputPoints array");
  
  TEST_FOR_EXCEPTION(  (inputPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::getValues_HGRAD_Args): dim 0 (number of points) > 0 required for inputPoints array");

  TEST_FOR_EXCEPTION( !(inputPoints.dimension(1) == spaceDim), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");

  
  // Verify that all inputPoints are in the reference cell
  /*
   TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
                       ">>> ERROR: (Intrepid::getValues_HGRAD_Args) One or more points are outside the " 
                       << cellTopo <<" reference cell");
   */
  
  
  // Verify that operatorType is admissible for HGRAD fields
  TEST_FOR_EXCEPTION( ( (spaceDim == 2) && (operatorType == OPERATOR_DIV) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) DIV is invalid operator for rank-0 (scalar) fields in 2D."); 
  
  TEST_FOR_EXCEPTION( ( (spaceDim == 3) && ( (operatorType == OPERATOR_DIV) ||
                                             (operatorType == OPERATOR_CURL) ) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) DIV and CURL are invalid operators for rank-0 (scalar) fields in 3D."); 
  
  
  // Check rank of outputValues (all operators are admissible in 1D) and its dim 2 when operator is
  // GRAD, CURL (only in 2D), or Dk.
  
  if(spaceDim == 1) {
    switch(operatorType){
      case OPERATOR_VALUE:
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 2), std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_DIV:
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
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 3), std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 3 required for outputValues in 1D when operator = GRAD, CURL, DIV, or Dk.");
        
        TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == 1 ),
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 2 of outputValues must equal 1 when operator = GRAD, CURL, DIV, or Dk.");
        break;
      default:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid::getValues_HGRAD_Args) Invalid operator");
    }
  }
  else if(spaceDim > 1) {
    switch(operatorType){
      case OPERATOR_VALUE:
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 2), std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_D1:
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 3), std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");
        
        TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == spaceDim ),
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 2 of outputValues must equal cell dimension when operator = GRAD, CURL (in 2D), or D1.");
        break;
      case OPERATOR_D2:
      case OPERATOR_D3:
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10:
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 3), std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");
        
        TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == Intrepid::getDkCardinality(operatorType, spaceDim) ),
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 2 of outputValues must equal cardinality of the Dk multiset.");
        break;
      default:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid::getValues_HGRAD_Args) Invalid operator");
    }
  }

  
  // Verify dim 0 and dim 1 of outputValues
  TEST_FOR_EXCEPTION( !(outputValues.dimension(1) == inputPoints.dimension(0) ), 
                      std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
  
  TEST_FOR_EXCEPTION( !(outputValues.dimension(0) == basisCard ),
                      std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HGRAD_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
}



template<class Scalar, class ArrayScalar>
void getValues_HCURL_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard){
  
  int spaceDim = cellTopo.getDimension();
  
  // Verify that cell is 2D or 3D (this is redundant for default bases where we use correct cells,
  //  but will force any user-defined basis for HCURL spaces to use 2D or 3D cells
  TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) cell dimension = 2 or 3 required for HCURL basis"); 

  
  // Verify inputPoints array
  TEST_FOR_EXCEPTION( !(inputPoints.rank() == 2), std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) rank = 2 required for inputPoints array"); 
  TEST_FOR_EXCEPTION(  (inputPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::getValues_HCURL_Args): dim 0 (number of points) > 0 required for inputPoints array");

  TEST_FOR_EXCEPTION( !(inputPoints.dimension(1) == spaceDim), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");
  
  // Verify that all inputPoints are in the reference cell
  /*
  TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) One or more points are outside the " 
                      << cellTopo <<" reference cell");
   */
  
  // Verify that operatorType is admissible for HCURL fields
  TEST_FOR_EXCEPTION( !( (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_CURL) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) operator = VALUE or CURL required for HCURL fields."); 
  
  
  // Check rank of outputValues: for VALUE should always be rank-3 array with (F,P,D) layout 
  switch(operatorType) {
    
    case OPERATOR_VALUE:
      TEST_FOR_EXCEPTION( !(outputValues.rank() == 3), std::invalid_argument,
                          ">>> ERROR: (Intrepid::getValues_HCURL_Args) rank = 3 required for outputValues when operator is VALUE");
      TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == spaceDim ),
                          std::invalid_argument,
                          ">>> ERROR: (Intrepid::getValues_HCURL_Args) dim 2 of outputValues must equal cell dimension when operator is VALUE.");
      break;
      
    case OPERATOR_CURL:
      
      // in 3D we need an (F,P,D) container because CURL gives a vector field:
      if(spaceDim == 3) {
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 3 ) ,
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HCURL_Args) rank = 3 required for outputValues in 3D when operator is CURL");
        TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == spaceDim),
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HCURL_Args) dim 2 of outputValues must equal cell dimension in 3D when operator is CURL.");
      }
      // In 2D we need an (F,P) container because CURL gives a scalar field
      else if(spaceDim == 2) {
        TEST_FOR_EXCEPTION( !(outputValues.rank() == 2 ) ,
                            std::invalid_argument,
                            ">>> ERROR: (Intrepid::getValues_HCURL_Args) rank = 2 required for outputValues in 2D when operator is CURL");
      }
      break;
      
    default:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid::getValues_HCURL_Args) Invalid operator");
  }
  
  
  // Verify dim 0 and dim 1 of outputValues
  TEST_FOR_EXCEPTION( !(outputValues.dimension(1) == inputPoints.dimension(0) ), 
                      std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
  
  TEST_FOR_EXCEPTION( !(outputValues.dimension(0) == basisCard ),
                      std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HCURL_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");

}



template<class Scalar, class ArrayScalar>
void getValues_HDIV_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard){
  
  int spaceDim = cellTopo.getDimension();
  
  // Verify inputPoints array
  TEST_FOR_EXCEPTION( !(inputPoints.rank() == 2), std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) rank = 2 required for inputPoints array"); 
  TEST_FOR_EXCEPTION(  (inputPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::getValues_HDIV_Args): dim 0 (number of points) > 0 required for inputPoints array");

  TEST_FOR_EXCEPTION( !(inputPoints.dimension(1) == spaceDim), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");
  
  // Verify that all inputPoints are in the reference cell
  /*
  TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) One or more points are outside the " 
                      << cellTopo <<" reference cell");
   */
  
  // Verify that operatorType is admissible for HDIV fields
  TEST_FOR_EXCEPTION( !( (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_DIV) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) operator = VALUE or DIV required for HDIV fields."); 
  
  
  // Check rank of outputValues 
  switch(operatorType) {
    case OPERATOR_VALUE:
      TEST_FOR_EXCEPTION( !(outputValues.rank() == 3), std::invalid_argument,
                          ">>> ERROR: (Intrepid::getValues_HDIV_Args) rank = 3 required for outputValues when operator is VALUE.");
      
      TEST_FOR_EXCEPTION( !(outputValues.dimension(2) == spaceDim ),
                          std::invalid_argument,
                          ">>> ERROR: (Intrepid::getValues_HDIV_Args) dim 2 of outputValues must equal cell dimension for operator VALUE.");
      break;
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( !(outputValues.rank() == 2), std::invalid_argument,
                          ">>> ERROR: (Intrepid::getValues_HDIV_Args) rank = 2 required for outputValues when operator is DIV.");
      break;

    default:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid::getValues_HDIV_Args) Invalid operator");
  }

  
  // Verify dim 0 and dim 1 of outputValues
  TEST_FOR_EXCEPTION( !(outputValues.dimension(1) == inputPoints.dimension(0) ), 
                      std::invalid_argument, 
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
  
  TEST_FOR_EXCEPTION( !(outputValues.dimension(0) == basisCard ),
                      std::invalid_argument,
                      ">>> ERROR: (Intrepid::getValues_HDIV_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
}

// Pure virtual destructor (gives warnings if not included).
// Following "Effective C++: 3rd Ed." item 7 the implementation
// is included in the definition file.
template<class ArrayScalar>
DofCoordsInterface<ArrayScalar>::~DofCoordsInterface() {}
