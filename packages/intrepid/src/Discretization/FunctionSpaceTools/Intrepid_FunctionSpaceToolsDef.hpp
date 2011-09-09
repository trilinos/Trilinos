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

/** \file   Intrepid_FunctionSpaceToolsDef.hpp
    \brief  Definition file for the Intrepid::FunctionSpaceTools class.
    \author Created by P. Bochev and D. Ridzal.
*/


namespace Intrepid {

template<class Scalar, class ArrayTypeOut, class ArrayTypeIn>
void FunctionSpaceTools::HGRADtransformVALUE(ArrayTypeOut       & outVals,
                                             const ArrayTypeIn  & inVals) {

  ArrayTools::cloneFields<Scalar>(outVals, inVals);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
void FunctionSpaceTools::HGRADtransformGRAD(ArrayTypeOut       & outVals,
                                            const ArrayTypeJac & jacobianInverse,
                                            const ArrayTypeIn  & inVals,
                                            const char           transpose) {

  ArrayTools::matvecProductDataField<Scalar>(outVals, jacobianInverse, inVals, transpose);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
void FunctionSpaceTools::HCURLtransformVALUE(ArrayTypeOut        & outVals,
                                             const ArrayTypeJac  & jacobianInverse,
                                             const ArrayTypeIn   & inVals,
                                             const char            transpose) {

  ArrayTools::matvecProductDataField<Scalar>(outVals, jacobianInverse, inVals, transpose);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HCURLtransformCURL(ArrayTypeOut        & outVals,
                                            const ArrayTypeJac  & jacobian,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeIn   & inVals,
                                            const char            transpose) {

  ArrayTools::matvecProductDataField<Scalar>(outVals, jacobian, inVals, transpose);
  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, jacobianDet, outVals, true);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HDIVtransformVALUE(ArrayTypeOut        & outVals,
                                            const ArrayTypeJac  & jacobian,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeIn   & inVals,
                                            const char            transpose) {

  ArrayTools::matvecProductDataField<Scalar>(outVals, jacobian, inVals, transpose);
  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, jacobianDet, outVals, true);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HDIVtransformDIV(ArrayTypeOut        & outVals,
                                          const ArrayTypeDet  & jacobianDet,
                                          const ArrayTypeIn   & inVals) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, jacobianDet, inVals, true);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HVOLtransformVALUE(ArrayTypeOut        & outVals,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeIn   & inVals) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, jacobianDet, inVals, true);

}


template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight>
void FunctionSpaceTools::integrate(ArrayOut            & outputValues,
                                   const ArrayInLeft   & leftValues,
                                   const ArrayInRight  & rightValues,
                                   const ECompEngine     compEngine,
                                   const bool            sumInto) {
  int outRank = outputValues.rank();

  switch (outRank) {
    case 1: 
      dataIntegral<Scalar>(outputValues, leftValues, rightValues, compEngine, sumInto);
    break;  
    case 2: 
      functionalIntegral<Scalar>(outputValues, leftValues, rightValues, compEngine, sumInto);
    break;  
    case 3: 
      operatorIntegral<Scalar>(outputValues, leftValues, rightValues, compEngine, sumInto);
    break;
    default:
      TEST_FOR_EXCEPTION( ((outRank != 1) && (outRank != 2) && (outRank != 3)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::integrate): Output container must have rank 1, 2 or 3.");
  }

} // integrate


template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void FunctionSpaceTools::operatorIntegral(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {
  int lRank = leftFields.rank();

  switch (lRank) {
    case 3: 
      ArrayTools::contractFieldFieldScalar<Scalar>(outputFields, leftFields, rightFields, compEngine, sumInto);
    break;  
    case 4: 
      ArrayTools::contractFieldFieldVector<Scalar>(outputFields, leftFields, rightFields, compEngine, sumInto);
    break;  
    case 5: 
      ArrayTools::contractFieldFieldTensor<Scalar>(outputFields, leftFields, rightFields, compEngine, sumInto);
    break;
    default:
      TEST_FOR_EXCEPTION( ((lRank != 3) && (lRank != 4) && (lRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::operatorIntegral): Left fields input container must have rank 3, 4 or 5.");
  }

} // operatorIntegral


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::functionalIntegral(ArrayOutFields &       outputFields,
                                            const ArrayInData &    inputData,
                                            const ArrayInFields &  inputFields,
                                            const ECompEngine      compEngine,
                                            const bool             sumInto) {
  int dRank = inputData.rank();

  switch (dRank) {
    case 2: 
      ArrayTools::contractDataFieldScalar<Scalar>(outputFields, inputData, inputFields, compEngine, sumInto);
    break;  
    case 3: 
      ArrayTools::contractDataFieldVector<Scalar>(outputFields, inputData, inputFields, compEngine, sumInto);
    break;  
    case 4: 
      ArrayTools::contractDataFieldTensor<Scalar>(outputFields, inputData, inputFields, compEngine, sumInto);
    break;
    default:
      TEST_FOR_EXCEPTION( ((dRank != 2) && (dRank != 3) && (dRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::functionalIntegral): Data input container must have rank 2, 3 or 4.");
  }

} // functionalIntegral


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::dataIntegral(ArrayOutData &            outputData,
                                      const ArrayInDataLeft &   inputDataLeft,
                                      const ArrayInDataRight &  inputDataRight,
                                      const ECompEngine         compEngine,
                                      const bool                sumInto) {
  int lRank = inputDataLeft.rank();

  switch (lRank) {
    case 2: 
      ArrayTools::contractDataDataScalar<Scalar>(outputData, inputDataLeft, inputDataRight, compEngine, sumInto);
    break;  
    case 3: 
      ArrayTools::contractDataDataVector<Scalar>(outputData, inputDataLeft, inputDataRight, compEngine, sumInto);
    break;  
    case 4: 
      ArrayTools::contractDataDataTensor<Scalar>(outputData, inputDataLeft, inputDataRight, compEngine, sumInto);
    break;
    default:
      TEST_FOR_EXCEPTION( ((lRank != 2) && (lRank != 3) && (lRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::dataIntegral): Left data input container must have rank 2, 3 or 4.");
  }

} // dataIntegral



template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
inline void FunctionSpaceTools::computeCellMeasure(ArrayOut             & outVals,
                                                   const ArrayDet       & inDet,
                                                   const ArrayWeights   & inWeights) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inDet.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Input determinants container must have rank 2.");
#endif

  ArrayTools::scalarMultiplyDataData<Scalar>(outVals, inDet, inWeights);
  // must use absolute value of inDet, so flip sign where needed
  for (int cell=0; cell<outVals.dimension(0); cell++) {
    if (inDet(cell,0) < 0.0) {
      for (int point=0; point<outVals.dimension(1); point++) {
        outVals(cell, point) *= -1.0;
      }
    }
  }

} // computeCellMeasure



template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
void FunctionSpaceTools::computeFaceMeasure(ArrayOut                   & outVals,
                                            const ArrayJac             & inJac,
                                            const ArrayWeights         & inWeights,
                                            const int                    whichFace,
                                            const shards::CellTopology & parentCell) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Input Jacobian container must have rank 4.");
#endif

  // temporary storage for face normals
  FieldContainer<Scalar> faceNormals(inJac.dimension(0), inJac.dimension(1), inJac.dimension(2));

  // compute normals
  CellTools<Scalar>::getPhysicalFaceNormals(faceNormals, inJac, whichFace, parentCell);

  // compute lenghts of normals
  RealSpaceTools<Scalar>::vectorNorm(outVals, faceNormals, NORM_TWO);

  // multiply with weights
  ArrayTools::scalarMultiplyDataData<Scalar>(outVals, outVals, inWeights);

}



template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
void FunctionSpaceTools::computeEdgeMeasure(ArrayOut                   & outVals,
                                            const ArrayJac             & inJac,
                                            const ArrayWeights         & inWeights,
                                            const int                    whichEdge,
                                            const shards::CellTopology & parentCell) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Input Jacobian container must have rank 4.");
#endif

  // temporary storage for edge tangents
  FieldContainer<Scalar> edgeTangents(inJac.dimension(0), inJac.dimension(1), inJac.dimension(2));

  // compute normals
  CellTools<Scalar>::getPhysicalEdgeTangents(edgeTangents, inJac, whichEdge, parentCell);

  // compute lenghts of tangents
  RealSpaceTools<Scalar>::vectorNorm(outVals, edgeTangents, NORM_TWO);

  // multiply with weights
  ArrayTools::scalarMultiplyDataData<Scalar>(outVals, outVals, inWeights);

}



template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
void FunctionSpaceTools::multiplyMeasure(ArrayTypeOut             & outVals,
                                         const ArrayTypeMeasure   & inMeasure,
                                         const ArrayTypeIn        & inVals) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, inMeasure, inVals);

} // multiplyMeasure


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::scalarMultiplyDataField(ArrayOutFields &     outputFields,
                                                 ArrayInData &        inputData,
                                                 ArrayInFields &      inputFields,
                                                 const bool           reciprocal) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outputFields, inputData, inputFields, reciprocal);

} // scalarMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::scalarMultiplyDataData(ArrayOutData &           outputData,
                                                ArrayInDataLeft &        inputDataLeft,
                                                ArrayInDataRight &       inputDataRight,
                                                const bool               reciprocal) {

  ArrayTools::scalarMultiplyDataData<Scalar>(outputData, inputDataLeft, inputDataRight, reciprocal);

} // scalarMultiplyDataData


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::dotMultiplyDataField(ArrayOutFields &       outputFields,
                                              const ArrayInData &    inputData,
                                              const ArrayInFields &  inputFields) {

  ArrayTools::dotMultiplyDataField<Scalar>(outputFields, inputData, inputFields);

} // dotMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::dotMultiplyDataData(ArrayOutData &            outputData,
                                             const ArrayInDataLeft &   inputDataLeft,
                                             const ArrayInDataRight &  inputDataRight) {

  ArrayTools::dotMultiplyDataData<Scalar>(outputData, inputDataLeft, inputDataRight);

} // dotMultiplyDataData


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::vectorMultiplyDataField(ArrayOutFields &       outputFields,
                                                 const ArrayInData &    inputData,
                                                 const ArrayInFields &  inputFields) {

  int outRank = outputFields.rank();

  switch (outRank) {
    case 3:
    case 4:
      ArrayTools::crossProductDataField<Scalar>(outputFields, inputData, inputFields);
      break;
    case 5:
      ArrayTools::outerProductDataField<Scalar>(outputFields, inputData, inputFields);
      break;
    default:
      TEST_FOR_EXCEPTION( ((outRank != 3) && (outRank != 4) && (outRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataField): Output container must have rank 3, 4 or 5.");
  }

} // vectorMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::vectorMultiplyDataData(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight) {

  int outRank = outputData.rank();

  switch (outRank) {
    case 2:
    case 3:
      ArrayTools::crossProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight);
      break;
    case 4:
      ArrayTools::outerProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight);
      break;
    default:
      TEST_FOR_EXCEPTION( ((outRank != 2) && (outRank != 3) && (outRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataData): Output container must have rank 2, 3 or 4.");
  }

} // vectorMultiplyDataData


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::tensorMultiplyDataField(ArrayOutFields &       outputFields,
                                                 const ArrayInData &    inputData,
                                                 const ArrayInFields &  inputFields,
                                                 const char             transpose) {

  int outRank = outputFields.rank();

  switch (outRank) {
    case 4:
      ArrayTools::matvecProductDataField<Scalar>(outputFields, inputData, inputFields, transpose);
      break;
    case 5:
      ArrayTools::matmatProductDataField<Scalar>(outputFields, inputData, inputFields, transpose);
      break;
    default:
      TEST_FOR_EXCEPTION( ((outRank != 4) && (outRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
  }

} // tensorMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::tensorMultiplyDataData(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight,
                                                const char                transpose) {

  int outRank = outputData.rank();

  switch (outRank) {
    case 3:
      ArrayTools::matvecProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
      break;
    case 4:
      ArrayTools::matmatProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
      break;
    default:
      TEST_FOR_EXCEPTION( ((outRank != 3) && (outRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataData): Output container must have rank 3 or 4.");
  }

} // tensorMultiplyDataData


template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
void FunctionSpaceTools::applyLeftFieldSigns(ArrayTypeInOut        & inoutOperator,
                                             const ArrayTypeSign   & fieldSigns) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inoutOperator.rank() != 3), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input operator container must have rank 3.");
  TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input field signs container must have rank 2.");
  TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
  TEST_FOR_EXCEPTION( (inoutOperator.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): First dimensions (number of left fields) of the operator and field signs containers must agree!");
#endif

  for (int cell=0; cell<inoutOperator.dimension(0); cell++) {
    for (int lbf=0; lbf<inoutOperator.dimension(1); lbf++) {
      for (int rbf=0; rbf<inoutOperator.dimension(2); rbf++) {
        inoutOperator(cell, lbf, rbf) *= fieldSigns(cell, lbf);
      }
    }
  }

} // applyLeftFieldSigns


template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
void FunctionSpaceTools::applyRightFieldSigns(ArrayTypeInOut        & inoutOperator,
                                              const ArrayTypeSign   & fieldSigns) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inoutOperator.rank() != 3), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input operator container must have rank 3.");
  TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input field signs container must have rank 2.");
  TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
  TEST_FOR_EXCEPTION( (inoutOperator.dimension(2) != fieldSigns.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Second dimension of the operator container and first dimension of the field signs container (number of right fields) must agree!");
#endif

  for (int cell=0; cell<inoutOperator.dimension(0); cell++) {
    for (int lbf=0; lbf<inoutOperator.dimension(1); lbf++) {
      for (int rbf=0; rbf<inoutOperator.dimension(2); rbf++) {
        inoutOperator(cell, lbf, rbf) *= fieldSigns(cell, rbf);
      }
    }
  }

} // applyRightFieldSigns


template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
void FunctionSpaceTools::applyFieldSigns(ArrayTypeInOut        & inoutFunction,
                                         const ArrayTypeSign   & fieldSigns) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ((inoutFunction.rank() < 2) || (inoutFunction.rank() > 5)), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input function container must have rank 2, 3, 4, or 5.");
  TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
  TEST_FOR_EXCEPTION( (inoutFunction.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the function and field signs containers must agree!");
  TEST_FOR_EXCEPTION( (inoutFunction.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): First dimensions (number of fields) of the function and field signs containers must agree!");
#endif

  int numCells  = inoutFunction.dimension(0);
  int numFields = inoutFunction.dimension(1);
  int fRank     = inoutFunction.rank();

  switch (fRank) {
    case 2: {
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          inoutFunction(cell, bf) *= fieldSigns(cell, bf);
        }
      }
    }
    break;
  
    case 3: {
      int numPoints = inoutFunction.dimension(2);
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          for (int pt=0; pt<numPoints; pt++) {
            inoutFunction(cell, bf, pt) *= fieldSigns(cell, bf);
          }
        }
      }
    }
    break;
  
    case 4: {
      int numPoints = inoutFunction.dimension(2);
      int spaceDim1 = inoutFunction.dimension(3);
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          for (int pt=0; pt<numPoints; pt++) {
            for (int d1=0; d1<spaceDim1; d1++) {
             inoutFunction(cell, bf, pt, d1) *= fieldSigns(cell, bf);
            }
          }
        }
      }
    }
    break;
  
    case 5: {
      int numPoints = inoutFunction.dimension(2);
      int spaceDim1 = inoutFunction.dimension(3);
      int spaceDim2 = inoutFunction.dimension(4);
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          for (int pt=0; pt<numPoints; pt++) {
            for (int d1=0; d1<spaceDim1; d1++) {
              for (int d2=0; d2<spaceDim2; d2++) {
                inoutFunction(cell, bf, pt, d1, d2) *= fieldSigns(cell, bf);
              }
            }
          }
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( !( (fRank == 2) || (fRank == 3) || (fRank == 4) || (fRank == 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Method defined only for rank-2, 3, 4, or 5 input function containers.");
  
  }  // end switch fRank

} // applyFieldSigns


template<class Scalar, class ArrayOutPointVals, class ArrayInCoeffs, class ArrayInFields>
void FunctionSpaceTools::evaluate(ArrayOutPointVals     & outPointVals,
                                  const ArrayInCoeffs   & inCoeffs,
                                  const ArrayInFields   & inFields) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ((inFields.rank() < 3) || (inFields.rank() > 5)), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Input fields container must have rank 3, 4, or 5.");
  TEST_FOR_EXCEPTION( (inCoeffs.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Input coefficient container must have rank 2.");
  TEST_FOR_EXCEPTION( (outPointVals.rank() != inFields.rank()-1), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Output values container must have rank one less than the rank of the input fields container.");
  TEST_FOR_EXCEPTION( (inCoeffs.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the coefficient and fields input containers must agree!");
  TEST_FOR_EXCEPTION( (inCoeffs.dimension(1) != inFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): First dimensions (number of fields) of the coefficient and fields input containers must agree!");
  TEST_FOR_EXCEPTION( (outPointVals.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the input fields container and the output values container must agree!");
  for (int i=1; i<outPointVals.rank(); i++) {
    std::string errmsg  = ">>> ERROR (FunctionSpaceTools::evaluate): Dimensions ";
    errmsg += (char)(48+i);
    errmsg += " and ";
    errmsg += (char)(48+i+1);
    errmsg += " of the output values and input fields containers must agree!";
    TEST_FOR_EXCEPTION( (outPointVals.dimension(i) != inFields.dimension(i+1)), std::invalid_argument, errmsg );
  }
#endif

  int numCells  = inFields.dimension(0);
  int numFields = inFields.dimension(1);
  int numPoints = inFields.dimension(2);
  int fRank     = inFields.rank();

  switch (fRank) {
    case 3: {
      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int bf=0; bf<numFields; bf++) {
            outPointVals(cell, pt) += inCoeffs(cell, bf) * inFields(cell, bf, pt);
          }
        }
      }
    }
    break;
  
    case 4: {
      int spaceDim1 = inFields.dimension(3);
      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d1=0; d1<spaceDim1; d1++) {
            for (int bf=0; bf<numFields; bf++) {
              outPointVals(cell, pt, d1) += inCoeffs(cell, bf) * inFields(cell, bf, pt, d1);
            }
          }
        }
      }
    }
    break;
  
    case 5: {
      int spaceDim1 = inFields.dimension(3);
      int spaceDim2 = inFields.dimension(4);
      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d1=0; d1<spaceDim1; d1++) {
            for (int d2=0; d2<spaceDim2; d2++) {
              for (int bf=0; bf<numFields; bf++) {
                outPointVals(cell, pt, d1, d2) += inCoeffs(cell, bf) * inFields(cell, bf, pt, d1, d2);
              }
            }
          }
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( !( (fRank == 3) || (fRank == 4) || (fRank == 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::evaluate): Method defined only for rank-3, 4, or 5 input fields containers.");
  
  }  // end switch fRank

} // evaluate

} // end namespace Intrepid
