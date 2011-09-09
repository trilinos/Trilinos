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

/** \file   Intrepid_ArrayToolsDefContractions.hpp
    \brief  Definition file for contraction (integration) operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {
  
template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldScalar(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {


#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftFields.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
  TEST_FOR_EXCEPTION( (rightFields.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftFields.dimension(0);
  int numLeftFields   = leftFields.dimension(1);
  int numRightFields  = rightFields.dimension(1);
  int numPoints       = leftFields.dimension(2);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += leftFields(cl, lbf, qp)*rightFields(cl, rbf, qp);
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
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
                tmpVal += leftFields(cl, lbf, qp)*rightFields(cl, rbf, qp);
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
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
                    alpha, &rightFields[cl*skipR], numPoints,
                    &leftFields[cl*skipL], numPoints,
                    beta, &outputFields[cl*skipOp], numRightFields);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");
  } // switch(compEngine)
} // end contractFieldFieldScalar


template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldVector(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftFields.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
  TEST_FOR_EXCEPTION( (rightFields.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(3) != rightFields.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftFields.dimension(0);
  int numLeftFields   = leftFields.dimension(1);
  int numRightFields  = rightFields.dimension(1);
  int numPoints       = leftFields.dimension(2);
  int dimVec          = leftFields.dimension(3);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numLeftFields; lbf++) {
            for (int rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += leftFields(cl, lbf, qp, iVec)*rightFields(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
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
                  tmpVal += leftFields(cl, lbf, qp, iVec)*rightFields(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
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
                    alpha, &rightFields[cl*skipR], numData,
                    &leftFields[cl*skipL], numData,
                    beta, &outputFields[cl*skipOp], numRightFields);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldVector): Computational engine not defined!");
  } // switch(compEngine)
} // end contractFieldFieldVector

    
template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldTensor(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftFields.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
  TEST_FOR_EXCEPTION( (rightFields.rank() != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(3) != rightFields.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftFields.dimension(4) != rightFields.dimension(4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftFields.dimension(0);
  int numLeftFields   = leftFields.dimension(1);
  int numRightFields  = rightFields.dimension(1);
  int numPoints       = leftFields.dimension(2);
  int dim1Tensor      = leftFields.dimension(3);
  int dim2Tensor      = leftFields.dimension(4);

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
                    tmpVal += leftFields(cl, lbf, qp, iTens1, iTens2)*rightFields(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
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
                    tmpVal += leftFields(cl, lbf, qp, iTens1, iTens2)*rightFields(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
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
                    alpha, &rightFields[cl*skipR], numData,
                    &leftFields[cl*skipL], numData,
                    beta, &outputFields[cl*skipOp], numRightFields);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractFieldFieldTensor): Computational engine not defined!");
  } // switch(compEngine)
} // end contractFieldFieldTensor
    

template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldScalar(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine      compEngine,
                                         const bool             sumInto) {


#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputFields.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
  TEST_FOR_EXCEPTION( (inputData.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
#endif

  // get sizes
  int numCells       = inputFields.dimension(0);
  int numFields      = inputFields.dimension(1);
  int numPoints      = inputFields.dimension(2);
  int numDataPoints  = inputData.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, qp);
              } // P-loop
              outputFields(cl, lbf) += tmpVal;
            } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, 0);
              } // P-loop
              outputFields(cl, lbf) += tmpVal;
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
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, qp);
              } // P-loop
              outputFields(cl, lbf) = tmpVal;
            } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, 0);
              } // P-loop
              outputFields(cl, lbf) = tmpVal;
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
                    alpha, &inputData[cl*skipR], numData,
                    &inputFields[cl*skipL], numData,
                    beta, &outputFields[cl*skipOp], 1);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldScalar): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataFieldScalar


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldVector(ArrayOutFields &      outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine      compEngine,
                                         const bool             sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputFields.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
  TEST_FOR_EXCEPTION( (inputData.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
#endif

  // get sizes
  int numCells       = inputFields.dimension(0);
  int numFields      = inputFields.dimension(1);
  int numPoints      = inputFields.dimension(2);
  int dimVec         = inputFields.dimension(3);
  int numDataPoints  = inputData.dimension(1);

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
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
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
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
              } // F-loop
          } // C-loop
        }
        else { // constant data
          for (int cl = 0; cl < numCells; cl++) {
              for (int lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (int qp = 0; qp < numPoints; qp++) {
                  for (int iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
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
                    alpha, &inputData[cl*skipR], numData,
                    &inputFields[cl*skipL], numData,
                    beta, &outputFields[cl*skipOp], 1);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldVector): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataFieldVector


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldTensor(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine      compEngine,
                                         const bool             sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputFields.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
  TEST_FOR_EXCEPTION( (inputData.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
  TEST_FOR_EXCEPTION( (outputFields.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
  TEST_FOR_EXCEPTION( (inputFields.dimension(4) != inputData.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
#endif

  // get sizes
  int numCells       = inputFields.dimension(0);
  int numFields      = inputFields.dimension(1);
  int numPoints      = inputFields.dimension(2);
  int dim1Tens       = inputFields.dimension(3);
  int dim2Tens       = inputFields.dimension(4);
  int numDataPoints  = inputData.dimension(1);

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
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
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
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
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
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
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
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
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
                    alpha, &inputData[cl*skipR], numData,
                    &inputFields[cl*skipL], numData,
                    beta, &outputFields[cl*skipOp], 1);
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
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataFieldTensor): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataFieldTensor


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataScalar(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine         compEngine,
                                        const bool                sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputDataLeft.rank()  != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputDataRight.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
  TEST_FOR_EXCEPTION( (outputData.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells      = inputDataLeft.dimension(0);
  int numPoints     = inputDataLeft.dimension(1);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            tmpVal += inputDataLeft(cl, qp)*inputDataRight(cl, qp);
          } // P-loop
          outputData(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            tmpVal += inputDataLeft(cl, qp)*inputDataRight(cl, qp);
          } // P-loop
          outputData(cl) = tmpVal;
        } // C-loop
      }
    }
    break;

    case COMP_BLAS: {
      int incr = 1;              // increment
      if (sumInto) {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputData(cl) += myblas.DOT(numPoints, &inputDataLeft[cl*numPoints], incr, &inputDataRight[cl*numPoints], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputData(cl) = myblas.DOT(numPoints, &inputDataLeft[cl*numPoints], incr, &inputDataRight[cl*numPoints], incr);
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataScalar): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataDataScalar


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataVector(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine         compEngine,
                                        const bool                sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputDataLeft.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
  TEST_FOR_EXCEPTION( (inputDataRight.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
  TEST_FOR_EXCEPTION( (outputData.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(2) != inputDataRight.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells        = inputDataLeft.dimension(0);
  int numPoints       = inputDataLeft.dimension(1);
  int dimVec          = inputDataLeft.dimension(2);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += inputDataLeft(cl, qp, iVec)*inputDataRight(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputData(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += inputDataLeft(cl, qp, iVec)*inputDataRight(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputData(cl) = tmpVal;
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
          outputData(cl) += myblas.DOT(skip, &inputDataLeft[cl*skip], incr, &inputDataRight[cl*skip], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputData(cl) = myblas.DOT(skip, &inputDataLeft[cl*skip], incr, &inputDataRight[cl*skip], incr);
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataVector): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataDataVector

    
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataTensor(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine         compEngine,
                                        const bool                sumInto) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputDataLeft.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
  TEST_FOR_EXCEPTION( (inputDataRight.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
  TEST_FOR_EXCEPTION( (outputData.rank() != 1 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(2) != inputDataRight.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (inputDataLeft.dimension(3) != inputDataRight.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif

  // get sizes
  int numCells        = inputDataLeft.dimension(0);
  int numPoints       = inputDataLeft.dimension(1);
  int dim1Tensor      = inputDataLeft.dimension(2);
  int dim2Tensor      = inputDataLeft.dimension(3);

  switch(compEngine) {
    case COMP_CPP: {
      if (sumInto) { 
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += inputDataLeft(cl, qp, iTens1, iTens2)*inputDataRight(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputData(cl) += tmpVal;
        } // C-loop
      }
      else {
        for (int cl = 0; cl < numCells; cl++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < numPoints; qp++) {
            for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (int iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += inputDataLeft(cl, qp, iTens1, iTens2)*inputDataRight(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputData(cl) = tmpVal;
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
          outputData(cl) += myblas.DOT(skip, &inputDataLeft[cl*skip], incr, &inputDataRight[cl*skip], incr);
        }
      }
      else {
        for (int cl=0; cl < numCells; cl++) {
          Teuchos::BLAS<int, Scalar> myblas;
          outputData(cl) = myblas.DOT(skip, &inputDataLeft[cl*skip], incr, &inputDataRight[cl*skip], incr);
        }
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractDataDataTensor): Computational engine not defined!");
  } // switch(compEngine)
} // end contractDataDataTensor


} // end namespace Intrepid
