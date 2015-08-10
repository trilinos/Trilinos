
namespace ArrayToolsKokkos{
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::scalarMultiplyDataData2(ArrayOutData &           outputData,
                                        ArrayInDataLeft &        inputDataLeft,
                                        ArrayInDataRight &       inputDataRight,
                                        const bool               reciprocal) {

#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
  if (outputData.rank() <= inputDataRight.rank()) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() < 2) || (inputDataRight.rank() > 4) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2, 3, or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != inputDataLeft.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first dimension of the left data input container must be 1!");
    for (int i=0; i<inputDataRight.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the right input and output data containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() < 1) || (inputDataRight.rank() > 3) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1, 2, or 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( inputDataLeft.dimension(0) != outputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
    for (int i=0; i<inputDataRight.rank(); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the right input and output data containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif

  // get sizes
  int invalRank      = inputDataRight.rank();
  int outvalRank     = outputData.rank();
  int numCells       = outputData.dimension(0);
  int numPoints      = outputData.dimension(1);
  int numDataPoints  = inputDataLeft.dimension(1);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }

  if (outvalRank == invalRank) {

    if (numDataPoints != 1) { // nonconstant data

          if (reciprocal) {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)/inputDataLeft(cl, pt);
              } // P-loop
            } // C-loop
          }
          else {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)*inputDataLeft(cl, pt);
              } // P-loop
            } // C-loop
          }
        }// case 2
      

     /*   default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 containers.");*/
      }// invalRank

    }
    else { // constant left data

          if (reciprocal) {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)/inputDataLeft(cl, 0);
              } // P-loop
            } // C-loop
          }
          else {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)*inputDataLeft(cl, 0);
              } // P-loop
            } // C-loop
          }
        }// case 2
      



      /*  default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");*/

      } // invalRank
    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

          if (reciprocal) {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)/inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
            } // C-loop
          }
          else {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)*inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
            } // C-loop
          }
        }// case 2
    


       // default:
       //       TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 1) || (invalRank == 2) || (invalRank == 3) ), std::invalid_argument,
       //                           ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 input containers.");
      }// invalRank

    }
    else { //constant data

if (reciprocal) {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                  for( int iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)/inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
            } // C-loop
          }
          else {
            for(int cl = 0; cl < numCells; cl++) {
              for(int pt = 0; pt < numPoints; pt++) {
                for( int iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)*inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
            } // C-loop
          }
        }// case 2
        

  

     /*   default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 1) || (invalRank == 2) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 input containers.");*/

      } // invalRank
    } // numDataPoints

  } // end if (outvalRank = invalRank)

} // scalarMultiplyDataData

} //ArrayToolsKokkos
