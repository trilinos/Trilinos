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
template<class Scalar>
void FunctionSpaceTools::integrate(Intrepid::FieldContainer<Scalar>            & outputValues,
                                   const Intrepid::FieldContainer<Scalar>   & leftValues,
                                   const Intrepid::FieldContainer<Scalar>  & rightValues,
                                   const ECompEngine           compEngine,
                                   const bool            sumInto) {
int outRank = getrank(outputValues);
int lRank = getrank(leftValues);
switch (outRank) {
 case 1:{
      switch (lRank){
    case 2:{
 #ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells      = leftValues.dimension(0);
  int numPoints     = leftValues.dimension(1);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            tmpVal += leftValues(cl, qp)*rightValues(cl, qp);
          } // P-loop
          outputValues(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            tmpVal += leftValues(cl, qp)*rightValues(cl, qp);
          } // P-loop
          outputValues(cl) = tmpVal;
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      int incr = 1;              // increment
      if (sumInto) {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) += myblas.DOT(numPoints, &leftValues[cl*numPoints], incr, &rightValues[cl*numPoints], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) = myblas.DOT(numPoints, &leftValues[cl*numPoints], incr, &rightValues[cl*numPoints], incr);
        }
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataScalar): Computational engine not defined!");
  } // switch(compEngine)  
    }
    break;
    case 3:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numPoints       = leftValues.dimension(1);
  int dimVec          = leftValues.dimension(2);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += leftValues(cl, qp, iVec)*rightValues(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputValues(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += leftValues(cl, qp, iVec)*rightValues(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputValues(cl) = tmpVal;
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      int skip = numPoints*dimVec;  // size of the left data chunk per cell
      int incr = 1;                 // increment
      if (sumInto) {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) += myblas.DOT(skip, &leftValues[cl*skip], incr, &rightValues[cl*skip], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) = myblas.DOT(skip, &leftValues[cl*skip], incr, &rightValues[cl*skip], incr);
        }
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataVector): Computational engine not defined!");
  } // switch(compEngine)


}

    break;
    case 4:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(3) != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numPoints       = leftValues.dimension(1);
  int dim1Tensor      = leftValues.dimension(2);
  int dim2Tensor      = leftValues.dimension(3);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) { 
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += leftValues(cl, qp, iTens1, iTens2)*rightValues(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputValues(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += leftValues(cl, qp, iTens1, iTens2)*rightValues(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputValues(cl) = tmpVal;
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      int skip = numPoints*dim1Tensor*dim2Tensor;  // size of the left data chunk per cell
      int incr = 1;                                // increment
      if (sumInto) {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) += myblas.DOT(skip, &leftValues[cl*skip], incr, &rightValues[cl*skip], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputValues(cl) = myblas.DOT(skip, &leftValues[cl*skip], incr, &rightValues[cl*skip], incr);
        }
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataTensor): Computational engine not defined!");
  } // switch(compEngine)

}
    break;
    default:
#ifdef HAVE_INTREPID_DEBUG

      TEUCHOS_TEST_FOR_EXCEPTION( ((lRank != 2) && (lRank != 3) && (lRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::dataIntegral): Left data input container must have rank 2, 3 or 4.");

#endif
    break;

}

}
 break;
 case 2:{
  switch (lRank) {
    case 2:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(0) != leftValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (rightValues.dimension(2) != leftValues.dimension(1)) && (leftValues.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
#endif

  // get sizes
  int numCells       = rightValues.dimension(0);
  int numFields      = rightValues.dimension(1);
  int numPoints      = rightValues.dimension(2);
  int numDataPoints  = leftValues.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += rightValues(cl, lbf, qp)*leftValues(cl, qp);
              } // P-loop
              outputValues(cl, lbf) += tmpVal;
            } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += rightValues(cl, lbf, qp)*leftValues(cl, 0);
              } // P-loop
              outputValues(cl, lbf) += tmpVal;
            } // F-loop
          } // C-loop
        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += rightValues(cl, lbf, qp)*leftValues(cl, qp);
              } // P-loop
              outputValues(cl, lbf) = tmpVal;
            } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += rightValues(cl, lbf, qp)*leftValues(cl, 0);
              } // P-loop
              outputValues(cl, lbf) = tmpVal;
            } // F-loop
          } // C-loop
        } // numDataPoints
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x f) and B(p x 1), respectively,
              even though the indexing is ((C),F,P) and ((C),P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = 1
       N        #cols(A)                              = number of input fields
       K        #cols(B^T)                            = number of integration points * size of data
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints;
      int skipL    = numFields*numPoints;       // size of the left data chunk per cell
      int skipR    = numPoints;                 // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &leftValues[cl*skipR], numData,
                    &rightValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputFields[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputFields[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldScalar): Computational engine not defined!");
  } // switch(compEngine)
    }
    break;
    case 3:{
    #ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(0) != leftValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (rightValues.dimension(2) != leftValues.dimension(1)) && (leftValues.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(3) != leftValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
#endif

  // get sizes
  int numCells       = rightValues.dimension(0);
  int numFields      = rightValues.dimension(1);
  int numPoints      = rightValues.dimension(2);
  int dimVec         = rightValues.dimension(3);
  int numDataPoints  = leftValues.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += rightValues(cl, lbf, qp, iVec)*leftValues(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputValues(cl, lbf) += tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += rightValues(cl, lbf, qp, iVec)*leftValues(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputValues(cl, lbf) += tmpVal;
              } // F-loop
          } // C-loop
        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += rightValues(cl, lbf, qp, iVec)*leftValues(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputValues(cl, lbf) = tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += rightValues(cl, lbf, qp, iVec)*leftValues(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputValues(cl, lbf) = tmpVal;
              } // F-loop
          } // C-loop
        } // numDataPoints
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x f) and B(p x 1), respectively,
              even though the indexing is ((C),F,P) and ((C),P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = 1
       N        #cols(A)                              = number of input fields
       K        #cols(B^T)                            = number of integration points * size of data
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints*dimVec;
      int skipL    = numFields*numData;         // size of the left data chunk per cell
      int skipR    = numData;                   // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &leftValues[cl*skipR], numData,
                    &rightValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputFields[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputFields[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldVector): Computational engine not defined!");
  } // switch(compEngine)
    }
    break;
    case 4:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(0) != leftValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (rightValues.dimension(2) != leftValues.dimension(1)) && (leftValues.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(3) != leftValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.dimension(4) != leftValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
#endif

  // get sizes
  int numCells       = rightValues.dimension(0);
  int numFields      = rightValues.dimension(1);
  int numPoints      = rightValues.dimension(2);
  int dim1Tens       = rightValues.dimension(3);
  int dim2Tens       = rightValues.dimension(4);
  int numDataPoints  = leftValues.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (int iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += rightValues(cl, lbf, qp, iTens1, iTens2)*leftValues(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputValues(cl, lbf) += tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += rightValues(cl, lbf, qp, iTens1, iTens2)*leftValues(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputValues(cl, lbf) += tmpVal;
              } // F-loop
          } // C-loop
        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (int iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += rightValues(cl, lbf, qp, iTens1, iTens2)*leftValues(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputValues(cl, lbf) = tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += rightValues(cl, lbf, qp, iTens1, iTens2)*leftValues(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputValues(cl, lbf) = tmpVal;
              } // F-loop
          } // C-loop
        } // numDataPoints
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x f) and B(p x 1), respectively,
              even though the indexing is ((C),F,P) and ((C),P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = 1
       N        #cols(A)                              = number of input fields
       K        #cols(B^T)                            = number of integration points * size of data
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints*dim1Tens*dim2Tens;
      int skipL    = numFields*numData;         // size of the left data chunk per cell
      int skipR    = numData;                   // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &leftValues[cl*skipR], numData,
                    &rightValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputFields[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputFields[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldTensor): Computational engine not defined!");
  } // switch(compEngine)    
}
    break;
  default:
#ifdef HAVE_INTREPID_DEBUG

TEUCHOS_TEST_FOR_EXCEPTION( ((lRank != 2) && (lRank != 3) && (lRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::functionalIntegral): Data input container must have rank 2, 3 or 4.");

#endif
    break;
}

}
 break;
 case 3:{
  switch (lRank) {
    case 3:{ 
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numLeftFields   = leftValues.dimension(1);
  int numRightFields  = rightValues.dimension(1);
  int numPoints       = leftValues.dimension(2);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += leftValues(cl, lbf, qp)*rightValues(cl, rbf, qp);
              } // P-loop
              outputValues(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += leftValues(cl, lbf, qp)*rightValues(cl, rbf, qp);
              } // P-loop
              outputValues(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x l) and B(p x r), respectively,
              even though the indexing is ((C),L,P) and ((C),R,P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = number of right fields
       N        #cols(A)                              = number of left fields
       K        #cols(B^T)                            = number of integration points
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points 
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int skipL    = numLeftFields*numPoints;       // size of the left data chunk per cell
      int skipR    = numRightFields*numPoints;      // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numRightFields, numLeftFields, numPoints,
                    alpha, &rightValues[cl*skipR], numPoints,
                    &leftValues[cl*skipL], numPoints,
                    beta, &outputValues[cl*skipOp], numRightFields);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numLeftFields, numRightFields, numPoints,
                    alpha, &leftFields[cl*skipL], numPoints,
                    &rightFields[cl*skipR], numPoints,
                    beta, &outputFields[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");
  } // switch(compEngine)

} 
   break;
    case 4:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(3) != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numLeftFields   = leftValues.dimension(1);
  int numRightFields  = rightValues.dimension(1);
  int numPoints       = leftValues.dimension(2);
  int dimVec          = leftValues.dimension(3);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += leftValues(cl, lbf, qp, iVec)*rightValues(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputValues(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += leftValues(cl, lbf, qp, iVec)*rightValues(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputValues(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x l) and B(p x r), respectively,
              even though the indexing is ((C),L,P) and ((C),R,P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = number of right fields
       N        #cols(A)                              = number of left fields
       K        #cols(B^T)                            = number of integration points * size of vector
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of vector
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of vector
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int numData  = numPoints*dimVec;              
      int skipL    = numLeftFields*numData;         // size of the left data chunk per cell
      int skipR    = numRightFields*numData;        // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numRightFields, numLeftFields, numData,
                    alpha, &rightValues[cl*skipR], numData,
                    &leftValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], numRightFields);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numLeftFields, numRightFields, numData,
                    alpha, &leftFields[cl*skipL], numData,
                    &rightFields[cl*skipR], numData,
                    beta, &outputFields[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldVector): Computational engine not defined!");
  } // switch(compEngine)
}
    break;
    case 5:{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (rightValues.rank() != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(3) != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftValues.dimension(4) != rightValues.dimension(4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numLeftFields   = leftValues.dimension(1);
  int numRightFields  = rightValues.dimension(1);
  int numPoints       = leftValues.dimension(2);
  int dim1Tensor      = leftValues.dimension(3);
  int dim2Tensor      = leftValues.dimension(4);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
                  for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                    tmpVal += leftValues(cl, lbf, qp, iTens1, iTens2)*rightValues(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputValues(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
                  for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                    tmpVal += leftValues(cl, lbf, qp, iTens1, iTens2)*rightValues(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputValues(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      /*
       GEMM parameters and their values.
       (Note: It is assumed that the result needs to be transposed into row-major format.
              Think of left and right input matrices as A(p x l) and B(p x r), respectively,
              even though the indexing is ((C),L,P) and ((C),R,P). Due to BLAS formatting
              assumptions, we are computing (A^T*B)^T = B^T*A.)
       TRANSA   TRANS
       TRANSB   NO_TRANS
       M        #rows(B^T)                            = number of right fields
       N        #cols(A)                              = number of left fields
       K        #cols(B^T)                            = number of integration points * size of tensor
       ALPHA    1.0
       A        right data for cell cl                = &rightFields[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of tensor
       B        left data for cell cl                 = &leftFields[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of tensor
       BETA     0.0
       C        result for cell cl                    = &outputFields[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int numData  = numPoints*dim1Tensor*dim2Tensor;              
      int skipL    = numLeftFields*numData;         // size of the left data chunk per cell
      int skipR    = numRightFields*numData;        // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);
      if (sumInto) {
        beta = 1.0;
      }

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numRightFields, numLeftFields, numData,
                    alpha, &rightValues[cl*skipR], numData,
                    &leftValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], numRightFields);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numLeftFields, numRightFields, numData,
                    alpha, &leftFields[cl*skipL], numData,
                    &rightFields[cl*skipR], numData,
                    beta, &outputFields[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldTensor): Computational engine not defined!");
  } // switch(compEngine)
}
    break;
    default:
#ifdef HAVE_INTREPID_DEBUG

      TEUCHOS_TEST_FOR_EXCEPTION( ((lRank != 3) && (lRank != 4) && (lRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::operatorIntegral): Left fields input container must have rank 3, 4 or 5.");

#endif
     break;
}



}
 break;
default:
#ifdef HAVE_INTREPID_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 1) && (outRank != 2) && (outRank != 3)), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::integrate): Output container must have rank 1, 2 or 3.");

#endif
break;
}
}
template<class Scalar, class ArrayOut, class ArrayInLeft, class ArrayInRight>
void FunctionSpaceTools::integrate(ArrayOut            & outputValues,
                                   const ArrayInLeft   & leftValues,
                                   const ArrayInRight  & rightValues,
                                   const ECompEngine           compEngine,
                                   const bool            sumInto) {
	 ArrayWrapper<Scalar,ArrayOut, Rank<ArrayOut >::value, false>outputValuesWrap(outputValues);
     ArrayWrapper<Scalar,ArrayInLeft, Rank<ArrayInLeft >::value, true>leftValuesWrap(leftValues);
	 ArrayWrapper<Scalar,ArrayInRight, Rank<ArrayInRight >::value, true>rightValuesWrap(rightValues);
	 int outRank = getrank(outputValues);
  switch (outRank) {
    case 1: 
      dataIntegral<Scalar>(outputValuesWrap, leftValuesWrap, rightValuesWrap, compEngine, sumInto);
    break;  
    case 2: 
      functionalIntegral<Scalar>(outputValuesWrap, leftValuesWrap, rightValuesWrap, compEngine, sumInto);
    break;  
    case 3: 
      operatorIntegral<Scalar>(outputValuesWrap, leftValuesWrap, rightValuesWrap, compEngine, sumInto);
    break;
  default:
#ifdef HAVE_INTREPID_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 1) && (outRank != 2) && (outRank != 3)), std::invalid_argument,
				">>> ERROR (FunctionSpaceTools::integrate): Output container must have rank 1, 2 or 3.");
    
#endif
   break;
  }									   

} // integrate
template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void FunctionSpaceTools::operatorIntegral(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {
  int lRank = getrank(leftFields);

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
#ifdef HAVE_INTREPID_DEBUG

      TEUCHOS_TEST_FOR_EXCEPTION( ((lRank != 3) && (lRank != 4) && (lRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::operatorIntegral): Left fields input container must have rank 3, 4 or 5.");

#endif
     break;
}

} // operatorIntegral


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void FunctionSpaceTools::functionalIntegral(ArrayOutFields &       outputFields,
                                            const ArrayInData &    inputData,
                                            const ArrayInFields &  inputFields,
                                            const ECompEngine           compEngine,
                                            const bool             sumInto) {
  int dRank = getrank(inputData);

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
#ifdef HAVE_INTREPID_DEBUG

TEUCHOS_TEST_FOR_EXCEPTION( ((dRank != 2) && (dRank != 3) && (dRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::functionalIntegral): Data input container must have rank 2, 3 or 4.");

#endif  
    break;
}

} // functionalIntegral

template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::dataIntegral(ArrayOutData &            outputData,
                                      const ArrayInDataLeft &   inputDataLeft,
                                      const ArrayInDataRight &  inputDataRight,
                                      const ECompEngine           compEngine,
                                      const bool                sumInto) {
  int lRank = getrank(inputDataLeft);

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
#ifdef HAVE_INTREPID_DEBUG

      TEUCHOS_TEST_FOR_EXCEPTION( ((lRank != 2) && (lRank != 3) && (lRank != 4)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::dataIntegral): Left data input container must have rank 2, 3 or 4.");

#endif
    break;
}

} // dataIntegral

template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
inline void FunctionSpaceTools::computeCellMeasure(ArrayOut             & outVals,
                                                   const ArrayDet       & inDet,
                                                   const ArrayWeights   & inWeights) {
#ifdef HAVE_INTREPID_DEBUG

  TEUCHOS_TEST_FOR_EXCEPTION( (inDet.rank() != 2), std::invalid_argument,
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

  TEUCHOS_TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
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

  TEUCHOS_TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
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

  int outRank = getrank(outputFields);

  switch (outRank) {
    case 3:
    case 4:
      ArrayTools::crossProductDataField<Scalar>(outputFields, inputData, inputFields);
      break;
    case 5:
      ArrayTools::outerProductDataField<Scalar>(outputFields, inputData, inputFields);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 3) && (outRank != 4) && (outRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataField): Output container must have rank 3, 4 or 5.");
  }

} // vectorMultiplyDataField


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void FunctionSpaceTools::vectorMultiplyDataData(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight) {

  int outRank = getrank(outputData);

  switch (outRank) {
    case 2:
    case 3:
      ArrayTools::crossProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight);
      break;
    case 4:
      ArrayTools::outerProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 2) && (outRank != 3) && (outRank != 4)), std::invalid_argument,
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
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 4) && (outRank != 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
  }

} // tensorMultiplyDataField

template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
struct FunctionSpaceTools::tensorMultiplyDataDataTempSpec<Scalar,  ArrayOutData,  ArrayInDataLeft,  ArrayInDataRight,-1>{
	tensorMultiplyDataDataTempSpec(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight,
                                                const char                transpose) {
	 int outRank = getrank(outputData);

  switch (outRank) {
    case 3:
      ArrayTools::matvecProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
      break;
    case 4:
      ArrayTools::matmatProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
      break;
	}
}
};
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
struct FunctionSpaceTools::tensorMultiplyDataDataTempSpec<Scalar,  ArrayOutData,  ArrayInDataLeft,  ArrayInDataRight,3>{
		tensorMultiplyDataDataTempSpec(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight,
                                                const char                transpose) {
	ArrayTools::matvecProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
	
}
};
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
struct FunctionSpaceTools::tensorMultiplyDataDataTempSpec<Scalar,  ArrayOutData,  ArrayInDataLeft,  ArrayInDataRight,4>{
		tensorMultiplyDataDataTempSpec(ArrayOutData &            outputData,
                                                const ArrayInDataLeft &   inputDataLeft,
                                                const ArrayInDataRight &  inputDataRight,
                                                const char                transpose) {
	 ArrayTools::matmatProductDataData<Scalar>(outputData, inputDataLeft, inputDataRight, transpose);
}

};
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  void FunctionSpaceTools::tensorMultiplyDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft &   inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight,
                                     const char                transpose){
		FunctionSpaceTools::tensorMultiplyDataDataTempSpec<Scalar,ArrayOutData,ArrayInDataLeft,ArrayInDataRight,Rank<ArrayOutData>::value>(outputData,inputDataLeft,inputDataRight,transpose);								 
										 
									 }
template<class Scalar, class ArrayTypeInOut, class ArrayTypeSign>
void FunctionSpaceTools::applyLeftFieldSigns(ArrayTypeInOut        & inoutOperator,
                                             const ArrayTypeSign   & fieldSigns) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.rank() != 3), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input operator container must have rank 3.");
  TEUCHOS_TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input field signs container must have rank 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
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
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inoutOperator) != 3), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input operator container must have rank 3.");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(fieldSigns) != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input field signs container must have rank 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(2) != fieldSigns.dimension(1) ), std::invalid_argument,
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

  TEUCHOS_TEST_FOR_EXCEPTION( ((inoutFunction.rank() < 2) || (inoutFunction.rank() > 5)), std::invalid_argument,
                              ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input function container must have rank 2, 3, 4, or 5.");
  TEUCHOS_TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                              ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutFunction.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                              ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the function and field signs containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inoutFunction.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
                              ">>> ERROR (FunctionSpaceTools::applyFieldSigns): First dimensions (number of fields) of the function and field signs containers must agree!");

#endif

 ArrayWrapper<Scalar,ArrayTypeInOut, Rank<ArrayTypeInOut >::value, false>inoutFunctionWrap(inoutFunction);
 ArrayWrapper<Scalar,ArrayTypeSign, Rank<ArrayTypeSign >::value, true>fieldSignsWrap(fieldSigns);


  int numCells  = inoutFunction.dimension(0);
  int numFields = inoutFunction.dimension(1);
  int fRank     = getrank(inoutFunction);

  switch (fRank) {
    case 2: {
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          inoutFunctionWrap(cell, bf) *= fieldSignsWrap(cell, bf);
        }
      }
    }
    break;
  
    case 3: {
      int numPoints = inoutFunction.dimension(2);
      for (int cell=0; cell<numCells; cell++) {
        for (int bf=0; bf<numFields; bf++) {
          for (int pt=0; pt<numPoints; pt++) {
            inoutFunctionWrap(cell, bf, pt) *= fieldSignsWrap(cell, bf);
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
             inoutFunctionWrap(cell, bf, pt, d1) *= fieldSignsWrap(cell, bf);
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
                inoutFunctionWrap(cell, bf, pt, d1, d2) *= fieldSignsWrap(cell, bf);
              }
            }
          }
        }
      }
    }
    break;

    default:
#ifdef HAVE_INTREPID_DEBUG

      TEUCHOS_TEST_FOR_EXCEPTION( !( (inoutFunction.rank() == 2) || (inoutFunction.rank() == 3) || (inoutFunction.rank() == 4) || (inoutFunction.rank() == 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Method defined only for rank-2, 3, 4, or 5 input function containers.");
#endif
    break;
  }  // end switch fRank

     

} // applyFieldSigns

template<class Scalar, class ArrayOutPointVals, class ArrayInCoeffs, class ArrayInFields>
void FunctionSpaceTools::evaluate(ArrayOutPointVals     & outPointVals,
                                  const ArrayInCoeffs   & inCoeffs,
                                  const ArrayInFields   & inFields) {

#ifdef HAVE_INTREPID_DEBUG

  TEUCHOS_TEST_FOR_EXCEPTION( ((inFields.rank() < 3) || (inFields.rank() > 5)), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Input fields container must have rank 3, 4, or 5.");
  TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Input coefficient container must have rank 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.rank() != inFields.rank()-1), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Output values container must have rank one less than the rank of the input fields container.");
  TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the coefficient and fields input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.dimension(1) != inFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): First dimensions (number of fields) of the coefficient and fields input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the input fields container and the output values container must agree!");
  for (int i=1; i<outPointVals.rank(); i++) {
    std::string errmsg  = ">>> ERROR (FunctionSpaceTools::evaluate): Dimensions ";
    errmsg += (char)(48+i);
    errmsg += " and ";
    errmsg += (char)(48+i+1);
    errmsg += " of the output values and input fields containers must agree!";
    TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.dimension(i) != inFields.dimension(i+1)), std::invalid_argument, errmsg );
  }

#endif
    ArrayWrapper<Scalar,ArrayOutPointVals, Rank<ArrayOutPointVals >::value, false>outPointValsWrap(outPointVals);
    ArrayWrapper<Scalar,ArrayInCoeffs, Rank<ArrayInCoeffs>::value, true>inCoeffsWrap(inCoeffs);
    ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>inFieldsWrap(inFields);

  int numCells  = inFields.dimension(0);
  int numFields = inFields.dimension(1);
  int numPoints = inFields.dimension(2);
    int fRank     = getrank(inFields);

  switch (fRank) {
    case 3: {
      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int bf=0; bf<numFields; bf++) {
            outPointValsWrap(cell, pt) += inCoeffsWrap(cell, bf) * inFieldsWrap(cell, bf, pt);
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
              outPointValsWrap(cell, pt, d1) += inCoeffsWrap(cell, bf) * inFieldsWrap(cell, bf, pt, d1);
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
                outPointValsWrap(cell, pt, d1, d2) += inCoeffsWrap(cell, bf) * inFieldsWrap(cell, bf, pt, d1, d2);
              }
            }
          }
        }
      }
    }
    break;

    default:
#ifdef HAVE_INTREPID_DEBUG

 TEUCHOS_TEST_FOR_EXCEPTION( !( (fRank == 3) || (fRank == 4) || (fRank == 5)), std::invalid_argument,
                          ">>> ERROR (FunctionSpaceTools::evaluate): Method defined only for rank-3, 4, or 5 input fields containers.");

#endif
 break;
  }  // end switch fRank

} // evaluate

} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

