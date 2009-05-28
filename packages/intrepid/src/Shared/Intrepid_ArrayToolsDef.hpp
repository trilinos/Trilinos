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

/** \file   Intrepid_ArrayToolsDef.hpp
    \brief  Definition file for utility class to provide array tools,
            such as tensor contractions, etc.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {
  
template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
void ArrayTools::contractScalar(ArrayTypeOut &         outputValues,
                                const ArrayTypeIn1 &   leftValues,
                                const ArrayTypeIn2 &   rightValues,
                                const ECompEngine      compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftValues.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Rank of the left input argument must equal 3");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Rank of right input argument must equal 3!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalar): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numLeftFields   = leftValues.dimension(1);
  int numRightFields  = rightValues.dimension(1);
  int numPoints       = leftValues.dimension(2);

  switch(compEngine) {
    case COMP_CPP: {
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points 
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int skipL    = numLeftFields*numPoints;       // size of the left data chunk per cell
      int skipR    = numRightFields*numPoints;      // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);

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
                    alpha, &leftValues[cl*skipL], numPoints,
                    &rightValues[cl*skipR], numPoints,
                    beta, &outputValues[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractScalar): Computational engine not defined!");
  } // switch(compEngine)
} // end contractScalar


template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
void ArrayTools::contractVector(ArrayTypeOut &         outputValues,
                                const ArrayTypeIn1 &   leftValues,
                                const ArrayTypeIn2 &   rightValues,
                                const ECompEngine      compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftValues.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Rank of the left input argument must equal 4");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Rank of right input argument must equal 4!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(3) != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVector): Second dimension of output container and first dimension of right input container must agree!");
#endif

  // get sizes
  int numCells        = leftValues.dimension(0);
  int numLeftFields   = leftValues.dimension(1);
  int numRightFields  = rightValues.dimension(1);
  int numPoints       = leftValues.dimension(2);
  int dimVec          = leftValues.dimension(3);

  switch(compEngine) {
    case COMP_CPP: {
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of vector
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of vector
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int numData  = numPoints*dimVec;              
      int skipL    = numLeftFields*numData;         // size of the left data chunk per cell
      int skipR    = numRightFields*numData;        // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);

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
                    alpha, &leftValues[cl*skipL], numData,
                    &rightValues[cl*skipR], numData,
                    beta, &outputValues[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractVector): Computational engine not defined!");
  } // switch(compEngine)
} // end contractVector

    
template<class Scalar, class ArrayTypeOut, class ArrayTypeIn1, class ArrayTypeIn2>
void ArrayTools::contractTensor(ArrayTypeOut &         outputValues,
                                const ArrayTypeIn1 &   leftValues,
                                const ArrayTypeIn2 &   rightValues,
                                const ECompEngine      compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (leftValues.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Rank of the left input argument must equal 5");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Rank of right input argument must equal 5!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Rank of output argument must equal 3!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(2) != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(3) != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (leftValues.dimension(4) != rightValues.dimension(4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != leftValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): First dimension of output container and first dimension of left input container must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(2) != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensor): Second dimension of output container and first dimension of right input container must agree!");
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
      for (int cl = 0; cl < numCells; cl++) {
        for (int lbf = 0; lbf < numLeftFields; lbf++) {
          for (int rbf = 0; rbf < numRightFields; rbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < numPoints; qp++) {
              for (int iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
                for (int iTens2 =0; iTens2 < dim2Tensor; iTens2++) {
                  tmpVal += leftValues(cl, lbf, qp, iTens1, iTens2)*rightValues(cl, rbf, qp, iTens1, iTens2);
                } // D2-loop
              } // D1-loop
            } // P-loop
            outputValues(cl, lbf, rbf) = tmpVal;
          } // R-loop
        } // L-loop
      } // C-loop
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of tensor
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of tensor
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = number of right fields
      */
      int numData  = numPoints*dim1Tensor*dim2Tensor;              
      int skipL    = numLeftFields*numData;         // size of the left data chunk per cell
      int skipR    = numRightFields*numData;        // size of the right data chunk per cell
      int skipOp   = numLeftFields*numRightFields;  // size of the output data chunk per cell
      Scalar alpha(1.0);                            // these are left unchanged by GEMM
      Scalar beta(0.0);

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
                    alpha, &leftValues[cl*skipL], numData,
                    &rightValues[cl*skipR], numData,
                    beta, &outputValues[cl*skipOp], numLeftFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractTensor): Computational engine not defined!");
  } // switch(compEngine)
} // end contractTensor
    

template<class Scalar, class ArrayTypeOut, class ArrayTypeIn, class ArrayTypeData>
void ArrayTools::contractScalarData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeIn &     inputValues,
                                    const ArrayTypeData &   inputData,
                                    const ECompEngine       compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputValues.rank()  != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Rank of the fields input argument must equal 3");
  TEST_FOR_EXCEPTION( (inputData.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Rank of the data input argument must equal 2!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != inputValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != inputValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractScalarData): First dimensions (number of fields) of output container and fields input container must agree!");
#endif

  // get sizes
  int numCells       = inputValues.dimension(0);
  int numFields      = inputValues.dimension(1);
  int numPoints      = inputValues.dimension(2);
  int numDataPoints  = inputData.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (numDataPoints != 1) { // nonconstant data
        for (int cl = 0; cl < numCells; cl++) {
          for (int lbf = 0; lbf < numFields; lbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < numPoints; qp++) {
              tmpVal += inputValues(cl, lbf, qp)*inputData(cl, qp);
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
              tmpVal += inputValues(cl, lbf, qp)*inputData(cl, 0);
            } // P-loop
            outputValues(cl, lbf) = tmpVal;
          } // F-loop
        } // C-loop
      } // numDataPoints
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints;
      int skipL    = numFields*numPoints;       // size of the left data chunk per cell
      int skipR    = numPoints;                 // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &inputData[cl*skipR], numData,
                    &inputValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputValues[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputValues[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractScalarData): Computational engine not defined!");
  } // switch(compEngine)
} // end contractScalarData


template<class Scalar, class ArrayTypeOut, class ArrayTypeIn, class ArrayTypeData>
void ArrayTools::contractVectorData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeIn &     inputValues,
                                    const ArrayTypeData &   inputData,
                                    const ECompEngine       compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputValues.rank()  != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Rank of the fields input argument must equal 4");
  TEST_FOR_EXCEPTION( (inputData.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Rank of the data input argument must equal 3!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != inputValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != inputValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractVectorData): First dimensions of output container and fields input container (number of fields) must agree!");
#endif

  // get sizes
  int numCells       = inputValues.dimension(0);
  int numFields      = inputValues.dimension(1);
  int numPoints      = inputValues.dimension(2);
  int dimVec         = inputValues.dimension(3);
  int numDataPoints  = inputData.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (numDataPoints != 1) { // nonconstant data
        for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += inputValues(cl, lbf, qp, iVec)*inputData(cl, qp, iVec);
                } //D-loop
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
                  tmpVal += inputValues(cl, lbf, qp, iVec)*inputData(cl, 0, iVec);
                } //D-loop
              } // P-loop
              outputValues(cl, lbf) = tmpVal;
            } // F-loop
        } // C-loop
      } // numDataPoints
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints*dimVec;
      int skipL    = numFields*numData;         // size of the left data chunk per cell
      int skipR    = numData;                   // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &inputData[cl*skipR], numData,
                    &inputValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputValues[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputValues[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractVectorData): Computational engine not defined!");
  } // switch(compEngine)
} // end contractVectorData


template<class Scalar, class ArrayTypeOut, class ArrayTypeIn, class ArrayTypeData>
void ArrayTools::contractTensorData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeIn &     inputValues,
                                    const ArrayTypeData &   inputData,
                                    const ECompEngine       compEngine) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputValues.rank()  != 5 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Rank of the fields input argument must equal 5");
  TEST_FOR_EXCEPTION( (inputData.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Rank of the data input argument must equal 4!");
  TEST_FOR_EXCEPTION( (outputValues.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Rank of output argument must equal 2!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
  TEST_FOR_EXCEPTION( (inputValues.dimension(4) != inputData.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(0) != inputValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEST_FOR_EXCEPTION( (outputValues.dimension(1) != inputValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractTensorData): First dimensions (number of fields) of output container and fields input container must agree!");
#endif

  // get sizes
  int numCells       = inputValues.dimension(0);
  int numFields      = inputValues.dimension(1);
  int numPoints      = inputValues.dimension(2);
  int dim1Tens       = inputValues.dimension(3);
  int dim2Tens       = inputValues.dimension(4);
  int numDataPoints  = inputData.dimension(1);

  ECompEngine myCompEngine = (numDataPoints == 1 ? COMP_CPP : compEngine);

  switch(myCompEngine) {
    case COMP_CPP: {
      if (numDataPoints != 1) { // nonconstant data
        for (int cl = 0; cl < numCells; cl++) {
            for (int lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (int qp = 0; qp < numPoints; qp++) {
                for (int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for (int iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                    tmpVal += inputValues(cl, lbf, qp, iTens1, iTens2)*inputData(cl, qp, iTens1, iTens2);
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
                  for (int iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                    tmpVal += inputValues(cl, lbf, qp, iTens1, iTens2)*inputData(cl, 0, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputValues(cl, lbf) = tmpVal;
            } // F-loop
        } // C-loop
      } // numDataPoints
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
       A        right data for cell cl                = &rightValues[cl*skipR]
       LDA      #rows(B)                              = number of integration points * size of data
       B        left data for cell cl                 = &leftValues[cl*skipL]
       LDB      #rows(A)                              = number of integration points * size of data
       BETA     0.0
       C        result for cell cl                    = &outputValues[cl*skipOp]
       LDC      #rows(C)                              = 1
      */
      int numData  = numPoints*dim1Tens*dim2Tens;
      int skipL    = numFields*numData;         // size of the left data chunk per cell
      int skipR    = numData;                   // size of the right data chunk per cell
      int skipOp   = numFields;                 // size of the output data chunk per cell
      Scalar alpha(1.0);                        // these are left unchanged by GEMM
      Scalar beta(0.0);

      for (int cl=0; cl < numCells; cl++) {
        /* Use this if data is used in row-major format */
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    1, numFields, numData,
                    alpha, &inputData[cl*skipR], numData,
                    &inputValues[cl*skipL], numData,
                    beta, &outputValues[cl*skipOp], 1);
        /* Use this if data is used in column-major format */
        /*
        myblas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS,
                    numFields, 1, numData,
                    alpha, &inputValues[cl*skipL], numData,
                    &inputData[cl*skipR], numData,
                    beta, &outputValues[cl*skipOp], numFields);
        */
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( ( ~isValidCompEngine(compEngine) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::contractTensorData): Computational engine not defined!");
  } // switch(compEngine)
} // end contractTensorData


template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
void ArrayTools::multiplyScalarData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeData &   inputData,
                                    ArrayTypeIn &           inputValues) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputData.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::multiplyScalarData): Input data container must have rank 2.");
  if (outputValues.rank() <= inputValues.rank()) {
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 3) || (inputValues.rank() > 5) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Input fields container must have rank 3, 4, or 5.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Input and output fields containers must have the same rank.");
    TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::multiplyScalarData): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 2) || (inputValues.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Input fields container must have rank 2, 3, or 4.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): The rank of the input fields container must be one less than the rank of the output fields container.");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
    TEST_FOR_EXCEPTION( ( inputData.dimension(0) != outputValues.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyScalarData): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::multiplyScalarData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif

  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int numDataPoints  = inputData.dimension(1);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim1Tens = outputValues.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputValues.dimension(4);
    }
  }
  
  if (outvalRank == invalRank) {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(cl, bf, pt)*inputData(cl, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(cl, bf, pt, iVec)*inputData(cl, pt);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(cl, bf, pt, iTens1, iTens2)*inputData(cl, pt);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyScalarData): This branch of the method is defined only for rank-3,4 or 5 containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(cl, bf, pt)*inputData(cl, 0);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(cl, bf, pt, iVec)*inputData(cl, 0);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(cl, bf, pt, iTens1, iTens2)*inputData(cl, 0);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyScalarData): This branch of the method is defined only for rank-3, 4 or 5 input containers.");

      } // invalRank
    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(bf, pt)*inputData(cl, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec)*inputData(cl, pt);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2)*inputData(cl, pt);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyScalarData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(bf, pt)*inputData(cl, 0);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec)*inputData(cl, 0);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2)*inputData(cl, 0);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyScalarData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");

      } // invalRank
    } // numDataPoints

  } // end if (outvalRank = invalRank)

}// multiplyScalarData


template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
void ArrayTools::divideByScalarData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeData &   inputData,
                                    ArrayTypeIn &           inputValues) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputData.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::divideByScalarData): Input data container must have rank 2.");
  if (outputValues.rank() <= inputValues.rank()) {
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 3) || (inputValues.rank() > 5) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Input fields container must have rank 3, 4, or 5.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Input and output fields containers must have the same rank.");
    TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::divideByScalarData): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 2) || (inputValues.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Input fields container must have rank 2, 3, or 4.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): The rank of the input fields container must be one less than the rank of the output fields container.");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
    TEST_FOR_EXCEPTION( ( inputData.dimension(0) != outputValues.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::divideByScalarData): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::divideByScalarData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif

  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int numDataPoints  = inputData.dimension(1);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim1Tens = outputValues.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputValues.dimension(4);
    }
  }
  
  if (outvalRank == invalRank) {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(cl, bf, pt)/inputData(cl, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(cl, bf, pt, iVec)/inputData(cl, pt);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(cl, bf, pt, iTens1, iTens2)/inputData(cl, pt);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::divideByScalarData): This branch of the method is defined only for rank-3,4 or 5 containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(cl, bf, pt)/inputData(cl, 0);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(cl, bf, pt, iVec)/inputData(cl, 0);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(cl, bf, pt, iTens1, iTens2)/inputData(cl, 0);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::divideByScalarData): This branch of the method is defined only for rank-3, 4 or 5 input containers.");

      } // invalRank
    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(bf, pt)/inputData(cl, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec)/inputData(cl, pt);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2)/inputData(cl, pt);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::divideByScalarData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputValues(cl, bf, pt) = inputValues(bf, pt)/inputData(cl, 0);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec)/inputData(cl, 0);
                } // D1-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                    outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2)/inputData(cl, 0);
                  }// D2-loop
                } // D1-loop
              }// F-loop
            } // P-loop
          }// C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::divideByScalarData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");

      } // invalRank
    } // numDataPoints

  } // end if (outvalRank = invalRank)

}// divideByScalarData


template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
void ArrayTools::multiplyVectorData(ArrayTypeOut &         outputValues,
                                    const ArrayTypeData &  inputData,
                                    const ArrayTypeIn &    inputValues) {

#ifdef HAVE_INTREPID_DEBUG
  if (outputValues.rank() < inputValues.rank()) {
  TEST_FOR_EXCEPTION( (inputData.rank() != 3), std::invalid_argument,
                      ">>> ERROR (ArrayTools::multiplyVectorData): Input data container must have rank 3.");
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 4) || (inputValues.rank() > 5) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Input fields container must have rank 4 or 5.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()-1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): The rank of output fields container must be one less than the rank of the input fields container.");
    TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    TEST_FOR_EXCEPTION( (inputValues.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Third dimension of the fields input container and second dimension of data input container (contracting spatial dimension) must agree!");
    for (int i=0; i<outputValues.rank(); i++) {
      int shift = 0;
      if (i > 2) {shift = 1;}
      std::string errmsg  = ">>> ERROR (ArrayTools::multiplyVectorData): Dimensions ";
      errmsg += (char)(48+i+shift);
      errmsg += " and ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i+shift) != outputValues.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEST_FOR_EXCEPTION( (inputData.rank() != 3), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Input data container must have rank 3.");
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 3) || (inputValues.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Input fields container must have rank 3 or 4.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): The rank of output fields container must equal the rank of the input fields container.");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): First dimensions (number of integration points) of fields input container and data input container must agree or first data dimension must be 1!");
    TEST_FOR_EXCEPTION( ( inputValues.dimension(1) != outputValues.dimension(2) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): First dimension of fields input container and second dimension of fields output container (number of integration points) must agree!");
    TEST_FOR_EXCEPTION( (inputValues.dimension(2) != inputData.dimension(2) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Second dimensions (contracting spatial dimension) of the fields and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( inputData.dimension(0) != outputValues.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Zeroth dimensions of fields output container and data input containers (numbers of integration domains) must agree!");
    TEST_FOR_EXCEPTION( ( inputValues.dimension(0) != outputValues.dimension(1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyVectorData): Zeroth dimension of fields input container and first dimension of fields output container (number of fields) must agree!");
    if (inputValues.rank() > 3) {
      TEST_FOR_EXCEPTION( (inputValues.dimension(3) != outputValues.dimension(3)), std::invalid_argument,
                           ">>> ERROR (ArrayTools::multiplyVectorData): Third dimensions of the fields input container and the fields output container (second spatial dimension) must agree.");
     }
  }
#endif

  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int numDataPoints  = inputData.dimension(1);
  int dim1Tens       = inputData.dimension(2);
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim2Tens = outputValues.dimension(3);
  }

  Scalar temp(0);

  if (invalRank == outvalRank + 1) {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                // All but the contracting dimension are the same as in the input fields container
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputData(cl, pt, iVec)*inputValues(cl, bf, pt, iVec);
                } // D1-loop
                outputValues(cl, bf, pt) = temp;
              }// P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  temp = 0;
                  // All but the contracting dimension are the same as in the input fields container
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputData(cl, pt, iTens1)*inputValues(cl, bf, pt, iTens1, iTens2);
                  } // D1-loop
                  outputValues(cl, bf, pt, iTens2) =  temp;
                } // D2-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyVectorData): This branch of the method is defined only for rank-4 or 5 input fields containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                // All but the contracting dimension are the same as in the input fields container
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputData(cl, 0, iVec)*inputValues(cl, bf, pt, iVec);
                } // D1-loop
                outputValues(cl, bf, pt) = temp;
              }// P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  temp = 0;
                  // All but the contracting dimension are the same as in the input fields container
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputData(cl, 0, iTens1)*inputValues(cl, bf, pt, iTens1, iTens2);
                  } // D1-loop
                  outputValues(cl, bf, pt, iTens2) =  temp;
                }// D2-loop
              } // P-loop
            } // F-loop
          }// C-loop
        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyVectorData): This branch of the method is defined only for rank-4 or 5 input fields containers.");
      }// invalRank

    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                // All but the contracting dimension are the same as in the input fields container
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputData(cl, pt, iVec)*inputValues(bf, pt, iVec);
                } // D1-loop
                outputValues(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  temp = 0;
                  // All but the contracting dimension are the same as in the input fields container
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputData(cl, pt, iTens1)*inputValues(bf, pt, iTens1, iTens2);
                  } // D1-loop
                  outputValues(cl, bf, pt, iTens2) =  temp;
                } // D2-loop
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyVectorData): This branch of the method is defined only for rank-3 or 4 input fields containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                // All but the contracting dimension are the same as in the input fields container
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputData(cl, 0, iVec)*inputValues(bf, pt, iVec);
                } // D1-loop
                outputValues(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  temp = 0;
                  // All but the contracting dimension are the same as in the input fields container
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputData(cl, 0, iTens1)*inputValues(bf, pt, iTens1, iTens2);
                  } // D1-loop
                  outputValues(cl, bf, pt, iTens2) =  temp;
                } // D2-loop
              } // P-loop
            } // F-loop
          }// C-loop
        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyVectorData): This branch of the method is defined only for rank-3 or 4 input fields containers.");
      }// invalRank

    } // numDataPoints

  } // end if (invalRank == outvalRank + 1)

}// multiplyVectorData


template<class Scalar, class ArrayTypeOut, class ArrayTypeData, class ArrayTypeIn>
void ArrayTools::multiplyTensorData(ArrayTypeOut &          outputValues,
                                    const ArrayTypeData &   inputData,
                                    const ArrayTypeIn &     inputValues, 
                                    const char              transpose) {

#ifdef HAVE_INTREPID_DEBUG
  if (outputValues.rank() <= inputValues.rank()) {
    TEST_FOR_EXCEPTION( ( (inputData.rank() < 2) || (inputData.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Input data container must have rank 2, 3 or 4.");
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 4) || (inputValues.rank() > 5) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Input fields container must have rank 4 or 5.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): The rank of output fields container must equal the rank of the input fields container.");
    TEST_FOR_EXCEPTION( ((transpose != 'n') && (transpose != 'N') && (transpose != 't') && (transpose != 'T')), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
    TEST_FOR_EXCEPTION( (inputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    if (inputData.rank() > 2) {
      TEST_FOR_EXCEPTION( (inputValues.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Third dimension of the fields input container and second dimension of data input container (contracting spatial dimension) must agree!");
    }
    if (inputData.rank() == 4) {
      TEST_FOR_EXCEPTION( (inputData.dimension(2) != inputData.dimension(3) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Second and third dimensions of the data input container must agree!");
    }
    if (inputValues.rank() == 5) {
      TEST_FOR_EXCEPTION( (inputValues.dimension(3) != inputValues.dimension(4) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Third and fourth dimensions of the fields input container must agree!");
    }
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::multiplyTensorData): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEST_FOR_EXCEPTION( ( (inputData.rank() < 2) || (inputData.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Input data container must have rank 2, 3 or 4.");
    TEST_FOR_EXCEPTION( ( (inputValues.rank() < 3) || (inputValues.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Input fields container must have rank 3 or 4.");
    TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): The rank of input fields container must be one less than the rank of the output fields container.");
    TEST_FOR_EXCEPTION( ((transpose != 'n') && (transpose != 'N') && (transpose != 't') && (transpose != 'T')), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
    TEST_FOR_EXCEPTION( (outputValues.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): Zeroth dimensions (number of integration domains) of the fields output and data input containers must agree!");
    TEST_FOR_EXCEPTION( ( (inputValues.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::multiplyTensorData): First dimensions of the fields input and the data input container (number of integration points) must agree or first data dimension must be 1!");
    if (inputData.rank() > 2) {
      TEST_FOR_EXCEPTION( (inputValues.dimension(2) != inputData.dimension(2) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Second dimensions (contracting spatial dimension) of the fields and data input containers must agree!");
    }
    if (inputData.rank() == 4) {
      TEST_FOR_EXCEPTION( (inputData.dimension(2) != inputData.dimension(3) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Second and third dimensions of the data input container must agree!");
    }
    if (inputValues.rank() == 4) {
      TEST_FOR_EXCEPTION( (inputValues.dimension(2) != inputValues.dimension(3) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::multiplyTensorData): Second and third dimensions of the fields input container must agree!");
    }
    for (int i=0; i<inputValues.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::multiplyTensorData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the input and output fields containers must agree!";
      TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif
           
  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int dataRank       = inputData.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int numDataPoints  = inputData.dimension(1);
  int dimTens        = outputValues.dimension(3);

  Scalar temp(0);

  if (invalRank == outvalRank) {

    if (numDataPoints != 1) { // nonconstant data  

      switch(invalRank) {

        case 4: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, pt)*inputValues(cl, bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, pt, iVec)*inputValues(cl, bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, pt, iTens, jTens)*inputValues(cl, bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, pt, jTens, iTens)*inputValues(cl, bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        } // case 4
        break;

        case 5: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, pt)*inputValues(cl, bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, pt, iTens)*inputValues(cl, bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, pt, iTens, kTens)*inputValues(cl, bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, pt, kTens, iTens)*inputValues(cl, bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyTensorData): This branch of the method is defined only for rank-4 or 5 input fields containers.");
      } // invalRank

    }
    else { // constant data

      switch(invalRank) {

        case 4: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, 0)*inputValues(cl, bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, 0, iVec)*inputValues(cl, bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, 0, iTens, jTens)*inputValues(cl, bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, 0, jTens, iTens)*inputValues(cl, bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        } // case 4
        break;

        case 5: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, 0)*inputValues(cl, bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, 0, iTens)*inputValues(cl, bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, 0, iTens, kTens)*inputValues(cl, bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, 0, kTens, iTens)*inputValues(cl, bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        }// case 5
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyTensorData): This branch of the method is defined only for rank-4 or 5 input fields containers.");
      } // invalRank

    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data  

      switch(invalRank) {

        case 3: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, pt)*inputValues(bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, pt, iVec)*inputValues(bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, pt, iTens, jTens)*inputValues(bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, pt, jTens, iTens)*inputValues(bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        } // case 3
        break;

        case 4: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, pt)*inputValues(bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, pt, iTens)*inputValues(bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, pt, iTens, kTens)*inputValues(bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, pt, kTens, iTens)*inputValues(bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyTensorData): This branch of the method is defined only for rank-3 or 4 input fields containers.");
      } // invalRank

    }
    else { // constant data

      switch(invalRank) {

        case 3: {

          switch(dataRank) {
            case 2: 
   for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, 0)*inputValues(bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iVec = 0; iVec < dimTens; iVec++) {
                      outputValues(cl, bf, pt, iVec) = inputData(cl, 0, iVec)*inputValues(bf, pt, iVec);
                    } // iVec-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, 0, iTens, jTens)*inputValues(bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        temp = 0;
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp += inputData(cl, 0, jTens, iTens)*inputValues(bf, pt, jTens);
                        } // jTens-loop
                        outputValues(cl, bf, pt, iTens) = temp;
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        } // case 3
        break;

        case 4: {

          switch(dataRank) {
            case 2: 
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, 0)*inputValues(bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 3:
              for(int cl = 0; cl < numCells; cl++) {
                for(int bf = 0; bf < numFields; bf++) {
                  for(int pt = 0; pt < numPoints; pt++) {
                    for( int iTens = 0; iTens < dimTens; iTens++) {
                      for( int jTens = 0; jTens < dimTens; jTens++) {
                        outputValues(cl, bf, pt, iTens, jTens) = inputData(cl, 0, iTens)*inputValues(bf, pt, iTens, jTens);
                      } // jTens-loop
                    } // iTens-loop
                  } // P-loop
                } // F-loop
              }// C-loop
            break;

            case 4:
              if ((transpose == 'n') || (transpose == 'N')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, 0, iTens, kTens)*inputValues(bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else if ((transpose == 't') || (transpose == 'T')) {
                for(int cl = 0; cl < numCells; cl++) {
                  for(int bf = 0; bf < numFields; bf++) {
                    for(int pt = 0; pt < numPoints; pt++) {
                      for( int iTens = 0; iTens < dimTens; iTens++) {
                        for( int jTens = 0; jTens < dimTens; jTens++) {
                          temp = 0;
                          for( int kTens = 0; kTens < dimTens; kTens++) {
                            temp += inputData(cl, 0, kTens, iTens)*inputValues(bf, pt, kTens, jTens);
                          } // kTens-loop
                          outputValues(cl, bf, pt, iTens, jTens) = temp;
                        } // jTens-loop
                      } // iTens-loop
                    } // P-loop
                  } // F-loop
                }// C-loop
              }
              else {
                TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): The transpose flag must be 'n', 'N', 't' or 'T'.");
              }
            break;

            default:
                TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::multiplyTensorData): This method is defined only for rank-2, 3 or 4 input data containers.");
          } // switch dataRank

        }// case 4
        break;

        default:
              TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::multiplyTensorData): This branch of the method is defined only for rank-3 or 4 input fields containers.");
      } // invalRank

    } // numDataPoints

  } // end if (invalRank == outvalRank)

} // multiplyTensorData


template<class Scalar, class ArrayTypeOut, class ArrayTypeIn>
void ArrayTools::cloneValues(ArrayTypeOut &          outputValues,
                                 const ArrayTypeIn &     inputValues) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (inputValues.rank() < 2) || (inputValues.rank() > 4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneValues): Input fields container must have rank 2, 3, or 4.");
  TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()+1), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneValues): The rank of the input fields container must be one less than the rank of the output fields container.");
  for (int i=0; i<inputValues.rank(); i++) {
    std::string errmsg  = ">>> ERROR (ArrayTools::cloneValues): Dimensions ";
    errmsg += (char)(48+i);
    errmsg += " and ";
    errmsg += (char)(48+i+1);
    errmsg += " of the input and output fields containers must agree!";
    TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i+1)), std::invalid_argument, errmsg );
  }
#endif

  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim1Tens = outputValues.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputValues.dimension(4);
    }
  }
  
  switch(invalRank) {
    case 2: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            outputValues(cl, bf, pt) = inputValues(bf, pt);
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 2
    break;

    case 3: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iVec = 0; iVec < dim1Tens; iVec++) {
              outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec);
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 3
    break;

    case 4: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
              for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 4
    break;

    default:
      TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneValues): This method is defined only for rank-2, 3 or 4 input containers.");
  }// invalRank

} // cloneValues


template<class Scalar, class ArrayTypeOut, class ArrayTypeFactors, class ArrayTypeIn>
void ArrayTools::cloneScaleValues(ArrayTypeOut &            outputValues,
                                  const ArrayTypeFactors &  inputFactors,
                                  const ArrayTypeIn &       inputValues) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputFactors.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleValues): The rank of the input factors container must be 2.");
  TEST_FOR_EXCEPTION( ( (inputValues.rank() < 2) || (inputValues.rank() > 4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleValues): Input fields container must have rank 2, 3, or 4.");
  TEST_FOR_EXCEPTION( (outputValues.rank() != inputValues.rank()+1), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleValues): The rank of the input fields container must be one less than the rank of the output fields container.");
  TEST_FOR_EXCEPTION( ( inputFactors.dimension(0) != outputValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleValues): Zeroth dimensions of input factors container and output fields container (numbers of integration domains) must agree!");
  TEST_FOR_EXCEPTION( ( inputFactors.dimension(1) != outputValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleValues): First dimensions of input factors container and output fields container (numbers of fields) must agree!");
  for (int i=0; i<inputValues.rank(); i++) {
    std::string errmsg  = ">>> ERROR (ArrayTools::cloneScaleValues): Dimensions ";
    errmsg += (char)(48+i);
    errmsg += " and ";
    errmsg += (char)(48+i+1);
    errmsg += " of the input and output fields containers must agree!";
    TEST_FOR_EXCEPTION( (inputValues.dimension(i) != outputValues.dimension(i+1)), std::invalid_argument, errmsg );
  }
#endif

  // get sizes
  int invalRank      = inputValues.rank();
  int outvalRank     = outputValues.rank();
  int numCells       = outputValues.dimension(0);
  int numFields      = outputValues.dimension(1);
  int numPoints      = outputValues.dimension(2);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim1Tens = outputValues.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputValues.dimension(4);
    }
  }
  
  switch(invalRank) {
    case 2: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            outputValues(cl, bf, pt) = inputValues(bf, pt) * inputFactors(cl, bf);
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 2
    break;

    case 3: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iVec = 0; iVec < dim1Tens; iVec++) {
              outputValues(cl, bf, pt, iVec) = inputValues(bf, pt, iVec) * inputFactors(cl, bf);
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 3
    break;

    case 4: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
              for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                outputValues(cl, bf, pt, iTens1, iTens2) = inputValues(bf, pt, iTens1, iTens2) * inputFactors(cl, bf);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 4
    break;

    default:
      TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneScaleValues): This method is defined only for rank-2, 3 or 4 input containers.");
  }// invalRank

} // cloneScaleValues


template<class Scalar, class ArrayTypeInOut, class ArrayTypeFactors>
void ArrayTools::scaleValues(ArrayTypeInOut &          inoutValues,
                             const ArrayTypeFactors &  inputFactors) {

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputFactors.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleValues): The rank of the input factors container must be 2.");
  TEST_FOR_EXCEPTION( ( (inoutValues.rank() < 3) || (inoutValues.rank() > 5) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleValues): Input/output fields container must have rank 3, 4, or 5.");
  TEST_FOR_EXCEPTION( ( inputFactors.dimension(0) != inoutValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleValues): Zeroth dimensions of input factors container and input/output fields container (numbers of integration domains) must agree!");
  TEST_FOR_EXCEPTION( ( inputFactors.dimension(1) != inoutValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleValues): First dimensions (number of fields) of input factors and input/output fields containers must agree!");
#endif

  // get sizes
  int inoutRank      = inoutValues.rank();
  int numCells       = inoutValues.dimension(0);
  int numFields      = inoutValues.dimension(1);
  int numPoints      = inoutValues.dimension(2);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (inoutRank > 3) {
    dim1Tens = inoutValues.dimension(3);
    if (inoutRank > 4) {
      dim2Tens = inoutValues.dimension(4);
    }
  }
  
  switch(inoutRank) {
    case 3: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            inoutValues(cl, bf, pt) = inoutValues(cl, bf, pt) * inputFactors(cl, bf);
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 2
    break;

    case 4: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iVec = 0; iVec < dim1Tens; iVec++) {
              inoutValues(cl, bf, pt, iVec) = inoutValues(cl, bf, pt, iVec) * inputFactors(cl, bf);
            } // D1-loop
          }// P-loop
        } // F-loop
      } // C-loop
    }// case 3
    break;

    case 5: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
              for( int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                inoutValues(cl, bf, pt, iTens1, iTens2) = inoutValues(cl, bf, pt, iTens1, iTens2) * inputFactors(cl, bf);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 4
    break;

    default:
      TEST_FOR_EXCEPTION( !( (inoutRank == 3) || (inoutRank == 4) || (inoutRank == 5) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneScaleValues): This method is defined only for rank-3, 4 or 5 input/output containers.");
  }// inoutRank

} // scaleValues


} // end namespace Intrepid
