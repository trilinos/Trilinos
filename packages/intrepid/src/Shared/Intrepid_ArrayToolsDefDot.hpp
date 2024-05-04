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

/** \file   Intrepid_ArrayToolsDefDot.hpp
    \brief  Definition file for dot-multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::dotMultiplyDataField(ArrayOutFields &       outputFields,
                                      const ArrayInData &    inputData,
                                      const ArrayInFields &  inputFields) {
#ifdef HAVE_INTREPID_DEBUG
  if (getrank(inputFields) > getrank(inputData)) {
    TEUCHOS_TEST_FOR_EXCEPTION( ((getrank(inputData) < 2) || (getrank(inputData) > 4)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFields) != getrank(inputData)+1), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Input fields container must have rank one larger than the rank of the input data container.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 3), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    for (size_t i=2; i<getrank(inputData); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataField): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the input data and fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputData.dimension(i) != inputFields.dimension(i+1)), std::invalid_argument, errmsg );
    }
    for (size_t i=0; i<getrank(outputFields); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataField): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(i) != outputFields.dimension(i)), std::invalid_argument, errmsg );
    }
  }else {
    TEUCHOS_TEST_FOR_EXCEPTION( ((getrank(inputData) < 2) || (getrank(inputData) > 4)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFields) != getrank(inputData)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): The rank of fields input container must equal the rank of data input container.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 3), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): First dimensions of the fields and data input containers (number of integration points) must agree or first data dimension must be 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != outputFields.dimension(1)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the fields input container and first dimension of the fields output container (number of fields) must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(1) != outputFields.dimension(2)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): First dimension of the fields input container and second dimension of the fields output container (number of integration points) must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputData.dimension(0)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions of the fields output and data input containers (number of integration domains) must agree!");
    for (size_t i=2; i<getrank(inputData); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataField): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " of the input data and fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputData.dimension(i) != inputFields.dimension(i)), std::invalid_argument, errmsg );
    }
  }
#endif


   ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);
   ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);
   ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value, true>inputFieldsWrap(inputFields);


  // get sizes
  int invalRank      = getrank(inputFields);
  int dataRank       = getrank(inputData);
  int numCells       = outputFields.dimension(0);
  int numFields      = outputFields.dimension(1);
  int numPoints      = outputFields.dimension(2);
  int numDataPoints  = inputData.dimension(1);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (dataRank > 2) {
    dim1Tens = inputData.dimension(2);
    if (dataRank > 3) {
      dim2Tens = inputData.dimension(3);
    }
  }

  Scalar temp(0);

  if (invalRank == dataRank + 1) {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputFieldsWrap(cl, bf, pt) = inputDataWrap(cl, pt)*inputFieldsWrap(cl, bf, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataWrap(cl, pt, iVec)*inputFieldsWrap(cl, bf, pt, iVec);
                } // D1-loop
                outputFieldsWrap(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataWrap(cl, pt, iTens1, iTens2)*inputFieldsWrap(cl, bf, pt, iTens1, iTens2);
                  } // D1-loop
                } // D2-loop
                outputFieldsWrap(cl, bf, pt) =  temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 5
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): This branch of the method is defined only for rank-3, 4 or 5 input fields containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputFieldsWrap(cl, bf, pt) = inputDataWrap(cl, 0)*inputFieldsWrap(cl, bf, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataWrap(cl, 0, iVec)*inputFieldsWrap(cl, bf, pt, iVec);
                } // D1-loop
                outputFieldsWrap(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        case 5: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataWrap(cl, 0, iTens1, iTens2)*inputFieldsWrap(cl, bf, pt, iTens1, iTens2);
                  } // D1-loop
                } // D2-loop
                outputFieldsWrap(cl, bf, pt) =  temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 5
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): This branch of the method is defined only for rank-3, 4 or 5 input fields containers.");
      }// invalRank

    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputFieldsWrap(cl, bf, pt) = inputDataWrap(cl, pt)*inputFieldsWrap(bf, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataWrap(cl, pt, iVec)*inputFieldsWrap(bf, pt, iVec);
                } // D1-loop
                outputFieldsWrap(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataWrap(cl, pt, iTens1, iTens2)*inputFieldsWrap(bf, pt, iTens1, iTens2);
                  } // D1-loop
                } // D2-loop
                outputFieldsWrap(cl, bf, pt) =  temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): This branch of the method is defined only for rank-2, 3 or 4 input fields containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                outputFieldsWrap(cl, bf, pt) = inputDataWrap(cl, 0)*inputFieldsWrap(bf, pt);
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataWrap(cl, 0, iVec)*inputFieldsWrap(bf, pt, iVec);
                } // D1-loop
                outputFieldsWrap(cl, bf, pt) = temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int bf = 0; bf < numFields; bf++) {
              for(int pt = 0; pt < numPoints; pt++) {
                temp = 0;
                for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                  for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataWrap(cl, 0, iTens1, iTens2)*inputFieldsWrap(bf, pt, iTens1, iTens2);
                  } // D1-loop
                } // D2-loop
                outputFieldsWrap(cl, bf, pt) =  temp;
              } // P-loop
            } // F-loop
          } // C-loop
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): This branch of the method is defined only for rank-2, 3 or 4 input fields containers.");
      }// invalRank

    } // numDataPoints

  } // end if (invalRank == dataRank + 1)

}// dotMultiplyDataField



template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::dotMultiplyDataData(ArrayOutData &            outputData,
                                     const ArrayInDataLeft  &  inputDataLeft,
                                     const ArrayInDataRight &  inputDataRight) {
#ifdef HAVE_INTREPID_DEBUG
  if (getrank(inputDataRight) >= getrank(inputDataLeft)) {
    TEUCHOS_TEST_FOR_EXCEPTION( ((getrank(inputDataLeft) < 2) || (getrank(inputDataLeft) > 4)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataRight) != getrank(inputDataLeft)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): The rank of the right data input container must equal the rank of the left data input container.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != 2), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): First dimensions of the left and right data input containers (number of integration points) must agree or first left data dimension must be 1!");
    for (size_t i=0; i<getrank(inputDataLeft); i++) {
      if (i != 1) {
	std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataData): Dimensions ";
        errmsg += (char)(48+i);
        errmsg += " of the left and right data input containers must agree!";
        TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(i) != inputDataRight.dimension(i)), std::invalid_argument, errmsg );
      }
    }
    for (size_t i=0; i<getrank(outputData); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " of the output and right input data containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( ((getrank(inputDataLeft) < 2) || (getrank(inputDataLeft) > 4)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataRight) != getrank(inputDataLeft)-1), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): Right data input container must have rank one less than the rank of left data input container.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != 2), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the right data input container and first dimension of left data input container (number of integration points) must agree or first left data dimension must be 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != outputData.dimension(1)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the right data input container and first dimension of output data container (number of integration points) must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != outputData.dimension(0)), std::invalid_argument,
				">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions of the left data input and data output containers (number of integration domains) must agree!");
    for (size_t i=1; i<getrank(inputDataRight); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::dotMultiplyDataData): Dimensions ";
      errmsg += (char)(48+i+1);
      errmsg += " and ";
      errmsg += (char)(48+i);
      errmsg += " of the left and right data input containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(i+1) != inputDataRight.dimension(i)), std::invalid_argument, errmsg );
    }
  }
#endif



   ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);
   ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);
   ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value, true>inputDataRightWrap(inputDataRight);


  // get sizes
  size_t rightDataRank  = getrank(inputDataRight);
  size_t leftDataRank   = getrank(inputDataLeft);
  int numCells       = outputData.dimension(0);
  int numPoints      = outputData.dimension(1);
  int numDataPoints  = inputDataLeft.dimension(1);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (leftDataRank > 2) {
    dim1Tens = inputDataLeft.dimension(2);
    if (leftDataRank > 3) {
      dim2Tens = inputDataLeft.dimension(3);
    }
  }

  Scalar temp(0);

  if (rightDataRank == leftDataRank) {

    if (numDataPoints != 1) { // nonconstant data

      switch(rightDataRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
                outputDataWrap(cl, pt) = inputDataLeftWrap(cl, pt)*inputDataRightWrap(cl, pt);
            } // P-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataLeftWrap(cl, pt, iVec)*inputDataRightWrap(cl, pt, iVec);
              } // D1-loop
              outputDataWrap(cl, pt) = temp;
            } // P-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataLeftWrap(cl, pt, iTens1, iTens2)*inputDataRightWrap(cl, pt, iTens1, iTens2);
                } // D1-loop
              } // D2-loop
              outputDataWrap(cl, pt) =  temp;
            } // P-loop
          } // C-loop
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (rightDataRank == 2) || (rightDataRank == 3) || (rightDataRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 right data input containers.");
      }// rightDataRank

    }
    else { //constant data

      switch(rightDataRank) {
        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
                outputDataWrap(cl, pt) = inputDataLeftWrap(cl, 0)*inputDataRightWrap(cl, pt);
            } // P-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataLeftWrap(cl, 0, iVec)*inputDataRightWrap(cl, pt, iVec);
              } // D1-loop
              outputDataWrap(cl, pt) = temp;
            } // P-loop
          } // C-loop
        }// case 3
        break;

        case 4: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataLeftWrap(cl, 0, iTens1, iTens2)*inputDataRightWrap(cl, pt, iTens1, iTens2);
                } // D1-loop
              } // D2-loop
              outputDataWrap(cl, pt) =  temp;
            } // P-loop
          } // C-loop
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (rightDataRank == 2) || (rightDataRank == 3) || (rightDataRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 right data input containers.");
      }// rightDataRank

    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(rightDataRank) {
        case 1: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
                outputDataWrap(cl, pt) = inputDataLeftWrap(cl, pt)*inputDataRightWrap(pt);
            } // P-loop
          } // C-loop
        }// case 1
        break;

        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataLeftWrap(cl, pt, iVec)*inputDataRightWrap(pt, iVec);
              } // D1-loop
              outputDataWrap(cl, pt) = temp;
            } // P-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataLeftWrap(cl, pt, iTens1, iTens2)*inputDataRightWrap(pt, iTens1, iTens2);
                } // D1-loop
              } // D2-loop
              outputDataWrap(cl, pt) =  temp;
            } // P-loop
          } // C-loop
        }// case 3
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (rightDataRank == 1) || (rightDataRank == 2) || (rightDataRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 right data input containers.");
      }// rightDataRank

    }
    else { //constant data

      switch(rightDataRank) {
        case 1: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
                outputDataWrap(cl, pt) = inputDataLeftWrap(cl, 0)*inputDataRightWrap(pt);
            } // P-loop
          } // C-loop
        }// case 1
        break;

        case 2: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for( int iVec = 0; iVec < dim1Tens; iVec++) {
                  temp += inputDataLeftWrap(cl, 0, iVec)*inputDataRightWrap(pt, iVec);
              } // D1-loop
              outputDataWrap(cl, pt) = temp;
            } // P-loop
          } // C-loop
        }// case 2
        break;

        case 3: {
          for(int cl = 0; cl < numCells; cl++) {
            for(int pt = 0; pt < numPoints; pt++) {
              temp = 0;
              for(int iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                for( int iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    temp += inputDataLeftWrap(cl, 0, iTens1, iTens2)*inputDataRightWrap(pt, iTens1, iTens2);
                } // D1-loop
              } // D2-loop
              outputDataWrap(cl, pt) =  temp;
            } // P-loop
          } // C-loop
        }// case 3
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (rightDataRank == 1) || (rightDataRank == 2) || (rightDataRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 right data input containers.");
      }// rightDataRank

    } // numDataPoints

  } // end if (rightDataRank == leftDataRank)

}// dotMultiplyDataData



} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

