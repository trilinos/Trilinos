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

/** \file   Intrepid_ArrayToolsDefCloneScale.hpp
    \brief  Definition file for clone / scale operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid2 {

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct cloneFields_2 {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  cloneFields_2 (ArrayOutFieldsWrap outputFields_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputFields (inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
        for(index_type bf = 0; bf < numFields; bf++) {
          for(index_type pt = 0; pt < numPoints; pt++) {
            outputFields(cl, bf, pt) = inputFields(bf, pt);
          } // P-loop
        } // F-loop
 
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct cloneFields_3 {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  cloneFields_3 (ArrayOutFieldsWrap outputFields_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputFields (inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank     = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
   
    }
  }
       for(index_type bf = 0; bf < numFields; bf++) {
          for(index_type pt = 0; pt < numPoints; pt++) {
            for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
              outputFields(cl, bf, pt, iVec) = inputFields(bf, pt, iVec);
            } // D1-loop
          } // P-loop
        } // F-loop
 
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInFieldsWrap,class ArrayOutFields>
struct cloneFields_4 {
  ArrayOutFieldsWrap outputFields;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayOutFields>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  cloneFields_4 (ArrayOutFieldsWrap outputFields_, ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputFields (inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
index_type outvalRank     = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputFields.dimension(4);
    }
  }
         for(index_type bf = 0; bf < numFields; bf++) {
          for(index_type pt = 0; pt < numPoints; pt++) {
            for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
              for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(bf, pt, iTens1, iTens2);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
 
   }
};


	template<class Scalar, class ArrayOutFields, class ArrayInFields>
void ArrayTools::cloneFields(ArrayOutFields &       outputFields,
                             const ArrayInFields &  inputFields) {

#ifdef HAVE_INTREPID2_DEBUG
	  TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputFields) < 2) || (getrank(inputFields) > 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::cloneFields): Input fields container must have rank 2, 3, or 4.");
	  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != getrank(inputFields)+1), std::invalid_argument,
				      ">>> ERROR (ArrayTools::cloneFields): The rank of the input fields container must be one less than the rank of the output fields container.");
	  for (index_type i=0; i<getrank(inputFields); i++) {
	    std::string errmsg  = ">>> ERROR (ArrayTools::cloneFields): Dimensions ";
	    errmsg += (char)(48+i);
	    errmsg += " and ";
	    errmsg += (char)(48+i+1);
	    errmsg += " of the input and output fields containers must agree!";
	    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(i) != outputFields.dimension(i+1)), std::invalid_argument, errmsg );
	  }
#endif
    ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value,false>outputFieldsWrap(outputFields);
    ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>inputFieldsWrap(inputFields);


  // get sizes
  index_type invalRank      = getrank(inputFields);
  int numCells       = outputFields.dimension(0);


  switch(invalRank) {
    case 2: {
 Kokkos::parallel_for (numCells, cloneFields_2<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value,false>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>, ArrayOutFields > (outputFieldsWrap,inputFieldsWrap));
    }// case 2
    break;

    case 3: {
 Kokkos::parallel_for (numCells, cloneFields_3<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value,false>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>, ArrayOutFields > (outputFieldsWrap,inputFieldsWrap));
    }// case 3
    break;

    case 4: {
Kokkos::parallel_for (numCells, cloneFields_4<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value,false>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>, ArrayOutFields > (outputFieldsWrap,inputFieldsWrap));
    }// case 4
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneFields): This method is defined only for rank-2, 3 or 4 input containers.");
  }// invalRank

} // cloneFields


template<class Scalar, class ArrayOutFields, class ArrayInFactors, class ArrayInFields>
void ArrayTools::cloneScaleFields(ArrayOutFields &        outputFields,
                                  const ArrayInFactors &  inputFactors,
                                  const ArrayInFields &   inputFields) {

#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFactors) != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleFields): The rank of the input factors container must be 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputFields) < 2) || (getrank(inputFields) > 4) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleFields): Input fields container must have rank 2, 3, or 4.");
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != getrank(inputFields)+1), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleFields): The rank of the input fields container must be one less than the rank of the output fields container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( inputFactors.dimension(0) != outputFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleFields): Zeroth dimensions of input factors container and output fields container (numbers of integration domains) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( inputFactors.dimension(1) != outputFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::cloneScaleFields): First dimensions of input factors container and output fields container (numbers of fields) must agree!");
  for (index_type i=0; i<getrank(inputFields); i++) {
    std::string errmsg  = ">>> ERROR (ArrayTools::cloneScaleFields): Dimensions ";
    errmsg += (char)(48+i);
    errmsg += " and ";
    errmsg += (char)(48+i+1);
    errmsg += " of the input and output fields containers must agree!";
    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(i) != outputFields.dimension(i+1)), std::invalid_argument, errmsg );
  }
#endif
   ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value,false>outputFieldsWrap(outputFields);
   ArrayWrapper<Scalar,ArrayInFactors, Rank<ArrayInFactors >::value,true>inputFactorswrap(inputFactors);
   ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>inputFieldsWrap(inputFields);
   
  // get sizes
  index_type invalRank      = getrank(inputFields);
  index_type outvalRank     = getrank(outputFields);
  int numCells       = outputFields.dimension(0);
  int numFields      = outputFields.dimension(1);
  int numPoints      = outputFields.dimension(2);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
      dim2Tens = outputFields.dimension(4);
    }
  }

  switch(invalRank) {
    case 2: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            outputFieldsWrap(cl, bf, pt) = inputFieldsWrap(bf, pt) * inputFactorswrap(cl, bf);
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
              outputFieldsWrap(cl, bf, pt, iVec) = inputFieldsWrap(bf, pt, iVec) * inputFactorswrap(cl, bf);
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
                outputFieldsWrap(cl, bf, pt, iTens1, iTens2) = inputFieldsWrap(bf, pt, iTens1, iTens2) * inputFactorswrap(cl, bf);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 4
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneScaleFields): This method is defined only for rank-2, 3 or 4 input containers.");
  }// invalRank

} // cloneScaleFields


template<class Scalar, class ArrayInOutFields, class ArrayInFactors>
void ArrayTools::scaleFields(ArrayInOutFields &      inoutFields,
                             const ArrayInFactors &  inputFactors) {

#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputFactors) != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleFields): The rank of the input factors container must be 2.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inoutFields) < 3) || (getrank(inoutFields) > 5) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleFields): Input/output fields container must have rank 3, 4, or 5.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( inputFactors.dimension(0) != inoutFields.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleFields): Zeroth dimensions of input factors container and input/output fields container (numbers of integration domains) must agree!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( inputFactors.dimension(1) != inoutFields.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scaleFields): First dimensions (number of fields) of input factors and input/output fields containers must agree!");
#endif
   ArrayWrapper<Scalar,ArrayInOutFields, Rank<ArrayInOutFields >::value,false>inoutFieldsWrap(inoutFields);
   ArrayWrapper<Scalar,ArrayInFactors, Rank<ArrayInFactors >::value,true>inputFactorsWrap(inputFactors);
  // get sizes
  index_type inoutRank      = getrank(inoutFields);
  int numCells       = inoutFields.dimension(0);
  int numFields      = inoutFields.dimension(1);
  int numPoints      = inoutFields.dimension(2);
  int dim1Tens       = 0;
  int dim2Tens       = 0;
  if (inoutRank > 3) {
    dim1Tens = inoutFields.dimension(3);
    if (inoutRank > 4) {
      dim2Tens = inoutFields.dimension(4);
    }
  }

  switch(inoutRank) {
    case 3: {
      for(int cl = 0; cl < numCells; cl++) {
        for(int bf = 0; bf < numFields; bf++) {
          for(int pt = 0; pt < numPoints; pt++) {
            inoutFieldsWrap(cl, bf, pt) = inoutFieldsWrap(cl, bf, pt) * inputFactorsWrap(cl, bf);
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
              inoutFieldsWrap(cl, bf, pt, iVec) = inoutFieldsWrap(cl, bf, pt, iVec) * inputFactorsWrap(cl, bf);
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
                inoutFieldsWrap(cl, bf, pt, iTens1, iTens2) = inoutFieldsWrap(cl, bf, pt, iTens1, iTens2) * inputFactorsWrap(cl, bf);
              } // D2-loop
            } // D1-loop
          } // P-loop
        } // F-loop
      } // C-loop
    }// case 4
    break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( (inoutRank == 3) || (inoutRank == 4) || (inoutRank == 5) ), std::invalid_argument,
                          ">>> ERROR (ArrayTools::cloneScaleFields): This method is defined only for rank-3, 4 or 5 input/output containers.");
  }// inoutRank

} // scaleFields


} // end namespace Intrepid2
