// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
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
inline void FunctionSpaceTools::computeMeasure(ArrayOut             & outVals,
                                               const ArrayDet       & inDet,
                                               const ArrayWeights   & inWeights) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inDet.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): Input determinants container must have rank 2.");
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

} // computeMeasure



template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
void FunctionSpaceTools::multiplyMeasure(ArrayTypeOut             & outVals,
                                         const ArrayTypeMeasure   & inMeasure,
                                         const ArrayTypeIn        & inVals) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outVals, inMeasure, inVals);

} // multiplyMeasure


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::scalarMultiplyDataField(ArrayOutFields &     outputFields,
                                                 const ArrayInData &  inputData,
                                                 ArrayInFields &      inputFields,
                                                 const bool           reciprocal) {

  ArrayTools::scalarMultiplyDataField<Scalar>(outputFields, inputData, inputFields, reciprocal);

} // scalarMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::scalarMultiplyDataData(ArrayOutData &           outputData,
                                                const ArrayInDataLeft &  inputDataLeft,
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
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of integration domains) of the operator and field signs containers must agree!");
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
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of integration domains) of the operator and field signs containers must agree!");
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
void FunctionSpaceTools::applyFieldSigns(ArrayTypeInOut        & inoutFunctional,
                                         const ArrayTypeSign   & fieldSigns) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inoutFunctional.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input functional container must have rank 2.");
  TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
  TEST_FOR_EXCEPTION( (inoutFunctional.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the operator and field signs containers must agree!");
  TEST_FOR_EXCEPTION( (inoutFunctional.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyFieldSigns): First dimensions (number of fields) of the functional and field signs containers must agree!");
#endif

  for (int cell=0; cell<inoutFunctional.dimension(0); cell++) {
    for (int bf=0; bf<inoutFunctional.dimension(1); bf++) {
      inoutFunctional(cell, bf) *= fieldSigns(cell, bf);
    }
  }

} // applyFieldSigns


} // end namespace Intrepid
