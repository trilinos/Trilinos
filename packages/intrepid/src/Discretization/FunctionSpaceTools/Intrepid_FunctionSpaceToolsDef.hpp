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

  ArrayTools::cloneValues<Scalar>(outVals, inVals);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeIn>
void FunctionSpaceTools::HGRADtransformGRAD(ArrayTypeOut       & outVals,
                                            const ArrayTypeJac & jacobianInverse,
                                            const ArrayTypeIn  & inVals,
                                            const char           transpose) {

  ArrayTools::multiplyTensorData<Scalar>(outVals, jacobianInverse, inVals, transpose);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeSign, class ArrayTypeIn>
void FunctionSpaceTools::HCURLtransformVALUE(ArrayTypeOut        & outVals,
                                             const ArrayTypeJac  & jacobianInverse,
                                             const ArrayTypeSign & fieldSigns,
                                             const ArrayTypeIn   & inVals,
                                             const char            transpose) {

  ArrayTools::multiplyTensorData<Scalar>(outVals, jacobianInverse, inVals, transpose);
  ArrayTools::scaleValues<Scalar>(outVals, fieldSigns);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
void FunctionSpaceTools::HCURLtransformCURL(ArrayTypeOut        & outVals,
                                            const ArrayTypeJac  & jacobian,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeSign & fieldSigns,
                                            const ArrayTypeIn   & inVals,
                                            const char            transpose) {

  ArrayTools::multiplyTensorData<Scalar>(outVals, jacobian, inVals, transpose);
  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, outVals);
  ArrayTools::scaleValues<Scalar>(outVals, fieldSigns);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeJac, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
void FunctionSpaceTools::HDIVtransformVALUE(ArrayTypeOut        & outVals,
                                            const ArrayTypeJac  & jacobian,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeSign & fieldSigns,
                                            const ArrayTypeIn   & inVals,
                                            const char            transpose) {

  ArrayTools::multiplyTensorData<Scalar>(outVals, jacobian, inVals, transpose);
  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, outVals);
  ArrayTools::scaleValues<Scalar>(outVals, fieldSigns);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeSign, class ArrayTypeIn>
void FunctionSpaceTools::HDIVtransformDIV(ArrayTypeOut        & outVals,
                                          const ArrayTypeDet  & jacobianDet,
                                          const ArrayTypeSign & fieldSigns,
                                          const ArrayTypeIn   & inVals) {

  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, inVals);
  ArrayTools::scaleValues<Scalar>(outVals, fieldSigns);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HVOLtransformVALUE(ArrayTypeOut        & outVals,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeIn   & inVals) {

  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, inVals);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
void FunctionSpaceTools::multiplyData(ArrayTypeOut             & outVals,
                                      const ArrayTypeData      & inData,
                                      const ArrayTypeIn        & inVals,
                                      const char               transpose) {

  int dataRank = inData.rank();

  switch (dataRank) {
    case 2:
      ArrayTools::multiplyScalarData<Scalar>(outVals, inData, inVals);
      break;
    case 3:
      ArrayTools::multiplyVectorData<Scalar>(outVals, inData, inVals);
      break;
    case 4:
      ArrayTools::multiplyTensorData<Scalar>(outVals, inData, inVals, transpose);
      break;
    default:
      TEST_FOR_EXCEPTION( ((dataRank != 2) && (dataRank != 3) && (dataRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::multiplyData): Data input container must have rank 2, 3 or 4.");
  }

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeMeasure, class ArrayTypeIn>
void FunctionSpaceTools::multiplyMeasure(ArrayTypeOut             & outVals,
                                         const ArrayTypeMeasure   & inMeasure,
                                         const ArrayTypeIn        & inVals) {

  ArrayTools::multiplyScalarData<Scalar>(outVals, inMeasure, inVals);

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
}


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
}


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
}


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
}


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

} // end namespace Intrepid
