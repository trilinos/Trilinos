// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_ArrayToolsDefContractions.hpp
    \brief  Definition file for contraction (integration) operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
*/
struct foo {
  typedef Kokkos::DefaultExecutionSpace execution_space ;
};



namespace Intrepid2 {
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldScalar_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldScalar_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 
  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += leftFields(cl, lbf, qp)*rightFields(cl, rbf, qp);
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldScalar_no_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldScalar_no_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {

  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += leftFields(cl, lbf, qp)*rightFields(cl, rbf, qp);
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop

   }
};


template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldScalar(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,                                          
                                          const bool                  sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(leftFields)  != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(rightFields) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");

  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");
#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>leftFieldsWrap(leftFields);
  ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>rightFieldsWrap(rightFields);
  // get sizes
  index_type numCells        = leftFields.dimension(0);


      if (sumInto) {

Kokkos::parallel_for (numCells, contractFieldFieldScalar_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));
      }
      else {
  Kokkos::parallel_for (numCells, contractFieldFieldScalar_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));    
}
} // end contractFieldFieldScalar


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldVector_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldVector_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
  index_type dimVec          = leftFields.dimension(3);

          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                for (index_type iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += leftFields(cl, lbf, qp, iVec)*rightFields(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldVector_no_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldVector_no_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {

  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
  index_type dimVec          = leftFields.dimension(3);
          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                for (index_type iVec = 0; iVec < dimVec; iVec++) {
                  tmpVal += leftFields(cl, lbf, qp, iVec)*rightFields(cl, rbf, qp, iVec);
                } //D-loop
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop


   }
};


template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldVector(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(leftFields)  != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(rightFields) != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(3) != rightFields.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>leftFieldsWrap(leftFields);
  ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>rightFieldsWrap(rightFields);

  // get sizes
  index_type numCells        = leftFields.dimension(0);


      if (sumInto) {
Kokkos::parallel_for (numCells, contractFieldFieldVector_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));
      }
      else {
Kokkos::parallel_for (numCells, contractFieldFieldVector_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));      

}
} // end contractFieldFieldVector


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldTensor_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldTensor_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
  index_type dim1Tensor      = leftFields.dimension(3);
  index_type dim2Tensor      = leftFields.dimension(4);


          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                for (index_type iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
                  for (index_type iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                    tmpVal += leftFields(cl, lbf, qp, iTens1, iTens2)*rightFields(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputFields(cl, lbf, rbf) += tmpVal;
            } // R-loop
          } // L-loop


   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsLeftWrap,class ArrayInFieldsRightWrap,class ArrayOutFields>
struct contractFieldFieldTensor_no_suminto {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsLeftWrap leftFields;
  ArrayInFieldsRightWrap rightFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractFieldFieldTensor_no_suminto (ArrayOutFieldsWrap outputFields_, ArrayInFieldsLeftWrap leftFields_, ArrayInFieldsRightWrap rightFields_) :
    outputFields (outputFields_),leftFields (leftFields_),rightFields(rightFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numLeftFields   = leftFields.dimension(1);
  index_type numRightFields  = rightFields.dimension(1);
  index_type numPoints       = leftFields.dimension(2);
  index_type dim1Tensor      = leftFields.dimension(3);
  index_type dim2Tensor      = leftFields.dimension(4);

          for (index_type lbf = 0; lbf < numLeftFields; lbf++) {
            for (index_type rbf = 0; rbf < numRightFields; rbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                for (index_type iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
                  for (index_type iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                    tmpVal += leftFields(cl, lbf, qp, iTens1, iTens2)*rightFields(cl, rbf, qp, iTens1, iTens2);
                  } // D2-loop
                } // D1-loop
              } // P-loop
              outputFields(cl, lbf, rbf) = tmpVal;
            } // R-loop
          } // L-loop


   }
};

    
template<class Scalar, class ArrayOutFields, class ArrayInFieldsLeft, class ArrayInFieldsRight>
void ArrayTools::contractFieldFieldTensor(ArrayOutFields &            outputFields,
                                          const ArrayInFieldsLeft &   leftFields,
                                          const ArrayInFieldsRight &  rightFields,
                                          const ECompEngine           compEngine,
                                          const bool                  sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(leftFields)  != 5 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(rightFields) != 5 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(2) != rightFields.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(3) != rightFields.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (leftFields.dimension(4) != rightFields.dimension(4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != rightFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != leftFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(2) != rightFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>leftFieldsWrap(leftFields);
  ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>rightFieldsWrap(rightFields);

  // get sizes
  index_type numCells        = leftFields.dimension(0);

      if (sumInto) {
      
Kokkos::parallel_for (numCells, contractFieldFieldTensor_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));
}
      else {
      Kokkos::parallel_for (numCells, contractFieldFieldTensor_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInFieldsLeft, Rank<ArrayInFieldsLeft>::value, true>, ArrayWrapper<Scalar,ArrayInFieldsRight, Rank<ArrayInFieldsRight>::value, true>,  ArrayOutFields> (outputFieldsWrap,leftFieldsWrap, rightFieldsWrap));

}
} // end contractFieldFieldTensor


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldScalar_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldScalar_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
          for (index_type lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, qp);
              } // P-loop
              outputFields(cl, lbf) += tmpVal;
            } // F-loop


   }
};
    
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldScalar_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldScalar_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);

            for (index_type lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, 0);
              } // P-loop
              outputFields(cl, lbf) += tmpVal;
            } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldScalar_no_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldScalar_no_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);

            for (index_type lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, qp);
              } // P-loop
              outputFields(cl, lbf) = tmpVal;
            } // F-loop

   }
};
    
 template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldScalar_no_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldScalar_no_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
            for (index_type lbf = 0; lbf < numFields; lbf++) {
              Scalar tmpVal(0);
              for (index_type qp = 0; qp < numPoints; qp++) {
                tmpVal += inputFields(cl, lbf, qp)*inputData(cl, 0);
              } // P-loop
              outputFields(cl, lbf) = tmpVal;
            } // F-loop


   }
};



template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldScalar(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine           compEngine,
                                         const bool             sumInto) {



#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFields)  != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputData) != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>inputDataWrap(inputData);
  ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>inputFieldsWrap(inputFields);

  // get sizes
  index_type numCells       = inputFields.dimension(0);
  index_type numDataPoints  = inputData.dimension(1);

      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldScalar_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldScalar_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldScalar_no_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldScalar_no_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));
       } // numDataPoints
      }
    
} // end contractDataFieldScalar


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldVector_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldVector_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dimVec         = inputFields.dimension(3);


              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
              } // F-loop


   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldVector_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldVector_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dimVec         = inputFields.dimension(3);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
              } // F-loop




   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldVector_no_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldVector_no_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {

  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dimVec         = inputFields.dimension(3);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, qp, iVec);
                  } // D-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
              } // F-loop


   }
};

 template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldVector_no_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldVector_no_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
            
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dimVec         = inputFields.dimension(3);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iVec = 0; iVec < dimVec; iVec++) {
                    tmpVal += inputFields(cl, lbf, qp, iVec)*inputData(cl, 0, iVec);
                  } //D-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
              } // F-loop



   }
};



template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldVector(ArrayOutFields &      outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine      compEngine,
                                         const bool             sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFields)  != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputData) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");


#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>inputDataWrap(inputData);
  ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>inputFieldsWrap(inputFields);

  // get sizes
  index_type numCells       = inputFields.dimension(0);
  index_type numDataPoints  = inputData.dimension(1);
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldVector_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldVector_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldVector_no_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldVector_no_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));
       } // numDataPoints
      }



 } // end contractDataFieldVector


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldTensor_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldTensor_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dim1Tens       = inputFields.dimension(3);
  index_type dim2Tens       = inputFields.dimension(4);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (index_type iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldTensor_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldTensor_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dim1Tens       = inputFields.dimension(3);
  index_type dim2Tens       = inputFields.dimension(4);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) += tmpVal;
              } // F-loop



   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldTensor_no_suminto_nonconst {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldTensor_no_suminto_nonconst (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dim1Tens       = inputFields.dimension(3);
  index_type dim2Tens       = inputFields.dimension(4);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (index_type iTens2 =0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, qp, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
              } // F-loop


   }
};

 template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct contractDataFieldTensor_no_suminto_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  contractDataFieldTensor_no_suminto_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = inputFields.dimension(1);
  index_type numPoints      = inputFields.dimension(2);
  index_type dim1Tens       = inputFields.dimension(3);
  index_type dim2Tens       = inputFields.dimension(4);

              for (index_type lbf = 0; lbf < numFields; lbf++) {
                Scalar tmpVal(0);
                for (index_type qp = 0; qp < numPoints; qp++) {
                  for (index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                    for (index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      tmpVal += inputFields(cl, lbf, qp, iTens1, iTens2)*inputData(cl, 0, iTens1, iTens2);
                    } // D2-loop
                  } // D1-loop
                } // P-loop
                outputFields(cl, lbf) = tmpVal;
              } // F-loop



   }
};


template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::contractDataFieldTensor(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields,
                                         const ECompEngine           compEngine,
                                         const bool             sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFields)  != 5 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputData) != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(3) != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(4) != inputData.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(0) != inputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputFields.dimension(1) != inputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif
  ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields>::value, false>outputFieldsWrap(outputFields);
  ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>inputDataWrap(inputData);
  ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>inputFieldsWrap(inputFields);

  // get sizes
  index_type numCells       = inputFields.dimension(0);
  index_type numDataPoints  = inputData.dimension(1);
      if (sumInto) {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldTensor_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldTensor_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        } // numDataPoints
      }
      else {
        if (numDataPoints != 1) { // nonconstant data
Kokkos::parallel_for (numCells, contractDataFieldTensor_no_suminto_nonconst<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));

        }
        else { // constant data
Kokkos::parallel_for (numCells, contractDataFieldTensor_no_suminto_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData>::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields>::value, true>,  ArrayOutFields> (outputFieldsWrap,inputData, inputFields));
       } // numDataPoints
      }



} // end contractDataFieldTensor

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataScalar_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataScalar_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 
  index_type numPoints     = inputDataLeft.dimension(1);
          Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            tmpVal += inputDataLeft(cl, qp)*inputDataRight(cl, qp);
          } // P-loop
          outputData(cl) += tmpVal;
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataScalar_no_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataScalar_no_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 
  index_type numPoints     = inputDataLeft.dimension(1);
          Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            tmpVal += inputDataLeft(cl, qp)*inputDataRight(cl, qp);
          } // P-loop
          outputData(cl) = tmpVal;
   }
};


template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataScalar(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine           compEngine,
                                        const bool                sumInto) {


#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataLeft)  != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataRight) != 2 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != 1 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif
  ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>outputDataWrap(outputData);
  ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>inputDataLeftWrap(inputDataLeft);
  ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>inputDataRightWrap(inputDataRight);
 // get sizes
  index_type numCells      = inputDataLeft.dimension(0);

      if (sumInto) {
 Kokkos::parallel_for (numCells, contractDataDataScalar_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }
      else {
 Kokkos::parallel_for (numCells, contractDataDataScalar_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }
    

} // end contractDataDataScalar


template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataVector_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataVector_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 
  index_type numPoints       = inputDataLeft.dimension(1);
  index_type dimVec          = inputDataLeft.dimension(2);

          Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            for (index_type iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += inputDataLeft(cl, qp, iVec)*inputDataRight(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputData(cl) += tmpVal;
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataVector_no_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataVector_no_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
   index_type numPoints       = inputDataLeft.dimension(1);
  index_type dimVec          = inputDataLeft.dimension(2);
          Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            for (index_type iVec = 0; iVec < dimVec; iVec++) {
              tmpVal += inputDataLeft(cl, qp, iVec)*inputDataRight(cl, qp, iVec);
            } // D-loop
          } // P-loop
          outputData(cl) = tmpVal;
   }
};

template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataVector(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine           compEngine,
                                        const bool                sumInto) {



#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataLeft)  != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataRight) != 3 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != 1 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(2) != inputDataRight.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif

  ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>outputDataWrap(outputData);
  ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>inputDataLeftWrap(inputDataLeft);
  ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>inputDataRightWrap(inputDataRight);

  // get sizes
  index_type numCells        = inputDataLeft.dimension(0);


      if (sumInto) {
 Kokkos::parallel_for (numCells, contractDataDataVector_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }
      else {
 Kokkos::parallel_for (numCells, contractDataDataVector_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }
 } // end contractDataDataVector

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataTensor_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataTensor_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
   index_type numPoints       = inputDataLeft.dimension(1);
   index_type dim1Tensor      = inputDataLeft.dimension(2);
   index_type dim2Tensor      = inputDataLeft.dimension(3);
           Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            for (index_type iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (index_type iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += inputDataLeft(cl, qp, iTens1, iTens2)*inputDataRight(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputData(cl) += tmpVal;
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayOutData>
struct contractDataDataTensor_no_suminto {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayOutData>::execution_space execution_space;
  contractDataDataTensor_no_suminto (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_, ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
   index_type numPoints       = inputDataLeft.dimension(1);
   index_type dim1Tensor      = inputDataLeft.dimension(2);
   index_type dim2Tensor      = inputDataLeft.dimension(3);
          Scalar tmpVal(0);
          for (index_type qp = 0; qp < numPoints; qp++) {
            for (index_type iTens1 = 0; iTens1 < dim1Tensor; iTens1++) {
              for (index_type iTens2 = 0; iTens2 < dim2Tensor; iTens2++) {
                tmpVal += inputDataLeft(cl, qp, iTens1, iTens2)*inputDataRight(cl, qp, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
          outputData(cl) = tmpVal;
   }
};
    
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::contractDataDataTensor(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight,
                                        const ECompEngine           compEngine,
                                        const bool                sumInto) {



#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataLeft)  != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataRight) != 4 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != 1 ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(1) != inputDataRight.dimension(1) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(2) != inputDataRight.dimension(2) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.dimension(3) != inputDataRight.dimension(3) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( (outputData.dimension(0) != inputDataRight.dimension(0) ), std::invalid_argument,
			      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( !isValidCompEngine(compEngine), std::invalid_argument,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Computational engine not defined!");

#endif

  ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>outputDataWrap(outputData);
  ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>inputDataLeftWrap(inputDataLeft);
  ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>inputDataRightWrap(inputDataRight);

  // get sizes
  index_type numCells        = inputDataLeft.dimension(0);



         if (sumInto) {
 Kokkos::parallel_for (numCells, contractDataDataTensor_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }
      else {
 Kokkos::parallel_for (numCells, contractDataDataTensor_no_suminto<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData>::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft>::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight>::value, true>,  ArrayOutData> (outputDataWrap,inputDataLeftWrap, inputDataRightWrap));
      }

 } // end contractDataDataTensor  


}
