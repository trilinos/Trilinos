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

  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, outVals);
  ArrayTools::scaleValues<Scalar>(outVals, fieldSigns);

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeDet, class ArrayTypeIn>
void FunctionSpaceTools::HVOLtransformVALUE(ArrayTypeOut        & outVals,
                                            const ArrayTypeDet  & jacobianDet,
                                            const ArrayTypeIn   & inVals) {

  ArrayTools::divideByScalarData<Scalar>(outVals, jacobianDet, outVals);

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


template<class Scalar, class ArrayTypeOut, class ArrayTypeInLeft, class ArrayTypeInRight>
void FunctionSpaceTools::integrate(ArrayTypeOut            & outputValues,
                                   const ArrayTypeInLeft   & leftValues,
                                   const ArrayTypeInRight  & rightValues,
                                   const ECompEngine         compEngine) {

  int lRank = leftValues.rank();

  if (lRank == rightValues.rank()) {
    switch (lRank) {
      case 3: 
        ArrayTools::contractScalar<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;  
      case 4: 
        ArrayTools::contractVector<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;  
      case 5: 
        ArrayTools::contractTensor<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;
      default:
        TEST_FOR_EXCEPTION( ((lRank != 3) && (lRank != 4) && (lRank != 5)), std::invalid_argument,
                            ">>> ERROR (FunctionSpaceTools::integrate): Left input container must have rank 3, 4 or 5.");
    }
  }
  else {
    switch (lRank) {
      case 3: 
        ArrayTools::contractScalarData<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;  
      case 4: 
        ArrayTools::contractVectorData<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;  
      case 5: 
        ArrayTools::contractTensorData<Scalar>(outputValues, leftValues, rightValues, compEngine);
      break;
      default:
        TEST_FOR_EXCEPTION( ((lRank != 3) && (lRank != 4) && (lRank != 5)), std::invalid_argument,
                            ">>> ERROR (FunctionSpaceTools::integrate): Left input container must have rank 3, 4 or 5.");
    }
  }

}


template<class Scalar, class ArrayTypeOut, class ArrayTypeWeights, class ArrayTypeDet>
void FunctionSpaceTools::computeMeasure(ArrayTypeOut             & outVals,
                                        const ArrayTypeWeights   & inWeights,
                                        const ArrayTypeDet       & inDet) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inWeights.rank() != 1), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): Input weights container must have rank 1.");
  TEST_FOR_EXCEPTION( (inDet.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): Input determinants container must have rank 2.");
  TEST_FOR_EXCEPTION( (outVals.rank() != inDet.rank()), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): Weighted measures and determinants arrays must have the same rank.");
  TEST_FOR_EXCEPTION( ( (inDet.dimension(1) != inWeights.dimension(0)) && (inDet.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): First dimension of the determinants container and zeroth dimension of the weights container (number of integration points) must agree or first determinants dimension must be 1!");
  TEST_FOR_EXCEPTION( (outVals.dimension(0) != inDet.dimension(0)), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): Zeroth dimensions of weighted measures and determinants arrays (number of integration domains) must agree.");
  TEST_FOR_EXCEPTION( (outVals.dimension(1) != inWeights.dimension(0)), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::computeMeasure): First dimension of weighted measures array and zeroth dimension of weights array (number of integration points) must agree.");
#endif

  int numCells         = inDet.dimension(0);
  int numWeightPoints  = inWeights.dimension(0);
  int numDetPoints     = inDet.dimension(1);

  if (numDetPoints > 1) {
    for(int cl = 0; cl < numCells; cl++) {
      for(int pt = 0; pt < numWeightPoints; pt++) {
        outVals(cl, pt) = std::abs(inDet(cl, pt))*inWeights(pt);
      } // P-loop
    } // C-loop
  }
  else {
    for(int cl = 0; cl < numCells; cl++) {
      for(int pt = 0; pt < numWeightPoints; pt++) {
        outVals(cl, pt) = std::abs(inDet(cl, 0))*inWeights(pt);
      } // P-loop
    } // C-loop
  }

} // computeMeasure



} // end namespace Intrepid
