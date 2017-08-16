#ifndef INTREPID2_ARRAYTOOLSDEFSCALAR_HPP
#define INTREPID2_ARRAYTOOLSDEFSCALAR_HPP
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

/** \file   Intrepid_ArrayToolsDefScalar.hpp
    \brief  Definition file for scalar multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
*/


namespace Intrepid2 {
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_3 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_3 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 index_type numFields      = outputFields.dimension(1);
 index_type numPoints      = outputFields.dimension(2);
            for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(cl, bf, pt)/inputData(cl, pt);
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_3 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_3 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 index_type numFields      = outputFields.dimension(1);
 index_type numPoints      = outputFields.dimension(2);
            for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(cl, bf, pt)*inputData(cl, pt);
                } // P-loop
              } // F-loop
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_4 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_4 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
      
    }
  }


            for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(cl, bf, pt, iVec)/inputData(cl, pt);
                  } // D1-loop
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_4 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_4 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
      
    }
  }


            for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(cl, bf, pt, iVec)*inputData(cl, pt);
                  } // D1-loop
                } // P-loop
              } // F-loop
   }
};


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_5 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_5 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(cl, bf, pt, iTens1, iTens2)/inputData(cl, pt);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_5 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_5 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(cl, bf, pt, iTens1, iTens2)*inputData(cl, pt);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_3_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_3_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(cl, bf, pt)/inputData(cl, 0);
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_3_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_3_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(cl, bf, pt)*inputData(cl, 0);
                } // P-loop
              } // F-loop
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_4_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_4_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
  
    }
  }
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(cl, bf, pt, iVec)/inputData(cl, 0);
                  } // D1-loop
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_4_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_4_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}
  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
  
    }
  }
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(cl, bf, pt, iVec)*inputData(cl, 0);
                  } // D1-loop
                } // P-loop
              } // F-loop
   }
};


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_recip_invalr_5_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_recip_invalr_5_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(cl, bf, pt, iTens1, iTens2)/inputData(cl, 0);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_eqrank_nrecip_invalr_5_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_eqrank_nrecip_invalr_5_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(cl, bf, pt, iTens1, iTens2)*inputData(cl, 0);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_2 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_2 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 index_type numFields      = outputFields.dimension(1);
 index_type numPoints      = outputFields.dimension(2);
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(bf, pt)/inputData(cl, pt);
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_2 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_nrecip_invalr_2 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}
  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
 index_type numFields      = outputFields.dimension(1);
 index_type numPoints      = outputFields.dimension(2);
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(bf, pt)*inputData(cl, pt);
                } // P-loop
              } // F-loop
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_3 {
   ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_3 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
  
    }
  }
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(bf, pt, iVec)/inputData(cl, pt);
                  } // D1-loop
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_3 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_nrecip_invalr_3 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
  
    }
  }
              for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(bf, pt, iVec)*inputData(cl, pt);
                  } // D1-loop
                } // P-loop
              } // F-loop
   }
};


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_4 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_4 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(bf, pt, iTens1, iTens2)/inputData(cl, pt);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_4 {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_nrecip_invalr_4 (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(bf, pt, iTens1, iTens2)*inputData(cl, pt);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop
   }
};
template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_2_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_2_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(bf, pt)/inputData(cl, 0);
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_2_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_nrecip_invalr_2_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  outputFields(cl, bf, pt) = inputFields(bf, pt)*inputData(cl, 0);
                } // P-loop
              } // F-loop
   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_3_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_3_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;

  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
   
    }
  }
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(bf, pt, iVec)/inputData(cl, 0);
                  } // D1-loop
                } // P-loop
              } // F-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_3_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_nrecip_invalr_3_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;

  index_type outvalRank = outputFields.rank();
  if (outvalRank > 3) {
    dim1Tens = outputFields.dimension(3);
    if (outvalRank > 4) {
   
    }
  }
             for(index_type bf = 0; bf < numFields; bf++) {
                for(index_type pt = 0; pt < numPoints; pt++) {
                  for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputFields(cl, bf, pt, iVec) = inputFields(bf, pt, iVec)*inputData(cl, 0);
                  } // D1-loop
                } // P-loop
              } // F-loop

   }
};


template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_recip_invalr_4_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  scalarMultiplyDataField_neqrank_recip_invalr_4_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(bf, pt, iTens1, iTens2)/inputData(cl, 0);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop

   }
};

template <class Scalar,class ArrayOutFieldsWrap,class ArrayInDataWrap,class ArrayInFieldsWrap,class ArrayInFields>
struct scalarMultiplyDataField_neqrank_nrecip_invalr_4_const {
  ArrayOutFieldsWrap outputFields;
  ArrayInDataWrap inputData;
  ArrayInFieldsWrap inputFields;
typedef typename conditional_eSpace<ArrayInFields>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataField_neqrank_nrecip_invalr_4_const (ArrayOutFieldsWrap outputFields_, ArrayInDataWrap inputData_,ArrayInFieldsWrap inputFields_) :
    outputFields (outputFields_),inputData (inputData_),inputFields(inputFields_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numFields      = outputFields.dimension(1);
  index_type numPoints      = outputFields.dimension(2);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  index_type outvalRank = outputFields.rank();
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
                      outputFields(cl, bf, pt, iTens1, iTens2) = inputFields(bf, pt, iTens1, iTens2)*inputData(cl, 0);
                    } // D2-loop
                  } // D1-loop
                } // F-loop
              } // P-loop
   }
};              

template<class Scalar, class ArrayOutFields, class ArrayInData, class ArrayInFields>
void ArrayTools::scalarMultiplyDataField(ArrayOutFields &     outputFields,
                                         const ArrayInData &  inputData,
                                         const ArrayInFields &      inputFields,
                                         const bool           reciprocal) {
#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputData) != 2), std::invalid_argument,
			      ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input data container must have rank 2.");
  if (getrank(outputFields) <= getrank(inputFields)) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputFields) < 3) || (getrank(inputFields) > 5) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 3, 4, or 5.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != getrank(inputFields)), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Input and output fields containers must have the same rank.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(0) != inputData.dimension(0) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(2) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    for (index_type i=0; i<getrank(inputFields); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataField): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the input and output fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(i) != outputFields.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputFields) < 2) || (getrank(inputFields) > 4) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 2, 3, or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputFields) != getrank(inputFields)+1), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): The rank of the input fields container must be one less than the rank of the output fields container.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputFields.dimension(1) != inputData.dimension(1)) && (inputData.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( inputData.dimension(0) != outputFields.dimension(0) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
    for (index_type i=0; i<getrank(inputFields); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataField): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the input and output fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputFields.dimension(i) != outputFields.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif
   ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);
   ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);
   ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>inputFieldsWrap(inputFields);

  // get sizes
  index_type invalRank      = getrank(inputFields);
  index_type outvalRank     = getrank(outputFields);
  index_type numCells       = outputFields.dimension(0);
  index_type numDataPoints  = inputData.dimension(1);
  if (outvalRank == invalRank) {

    if (numDataPoints != 1) { // nonconstant data
      switch(invalRank) {
        case 3: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
			  
			
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
 
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 4
        break;

        case 5: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_5<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_5<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 5
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): This branch of the method is defined only for rank-3,4 or 5 containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 3: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 4
        break;

        case 5: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_recip_invalr_5_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_eqrank_nrecip_invalr_5_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 5
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 3) || (invalRank == 4) || (invalRank == 5) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): This branch of the method is defined only for rank-3, 4 or 5 input containers.");

      } // invalRank
    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 2: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 2
        break;

        case 3: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): This branch of the method is defined only for rank-2, 3 or 4 input containers.");
      }// invalRank

    }
    else { //constant data

      switch(invalRank) {
        case 2: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 2
        break;

        case 3: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_recip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataField_neqrank_nrecip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>, ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>, ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>,  ArrayInFields> (outputFieldsWrap,inputDataWrap, inputFieldsWrap));
          }
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): This branch of the method is defined only for rank-2, 3 or 4 input containers.");

      } // invalRank
    } // numDataPoints

  } // end if (outvalRank = invalRank)

} // scalarMultiplyDataField

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_2 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_2 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numPoints      = outputData.dimension(1);
               for(index_type pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)/inputDataLeft(cl, pt);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_2 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_2 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  int numPoints      = outputData.dimension(1);
               for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)*inputDataLeft(cl, pt);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_3 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_3 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
 
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(cl, pt, iVec)/inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_3 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_3 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(cl, pt, iVec)*inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
   }
};


template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_4 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_4 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(cl, pt, iTens1, iTens2)/inputDataLeft(cl, pt);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_4 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_4 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(cl, pt, iTens1, iTens2)*inputDataLeft(cl, pt);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_2_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_2_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numPoints      = outputData.dimension(1);
               for(index_type pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)/inputDataLeft(cl, 0);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_2_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_2_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  int numPoints      = outputData.dimension(1);
               for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(cl, pt)*inputDataLeft(cl, 0);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_3_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_3_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
 
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                   outputData(cl, pt, iVec) = inputDataRight(cl, pt, iVec)/inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_3_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_3_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;

  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(cl, pt, iVec)*inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
   }
};


template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_recip_invalr_4_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_recip_invalr_4_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(cl, pt, iTens1, iTens2)/inputDataLeft(cl, 0);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_eqrank_nrecip_invalr_4_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_eqrank_nrecip_invalr_4_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                     outputData(cl, pt, iTens1, iTens2) = inputDataRight(cl, pt, iTens1, iTens2)*inputDataLeft(cl, 0);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_1 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_1 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {

  index_type numPoints      = outputData.dimension(1);
               for(index_type pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(pt)/inputDataLeft(cl, pt);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_1 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_1 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  int numPoints      = outputData.dimension(1);
               for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(pt)*inputDataLeft(cl, pt);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_2 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_2 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;

  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)/inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_2 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_2 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;

  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)*inputDataLeft(cl, pt);
                } // D1-loop
              } // P-loop
   }
};


template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_3 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_3 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(pt, iTens1, iTens2)/inputDataLeft(cl, pt);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_3 {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_3 (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(pt, iTens1, iTens2)*inputDataLeft(cl, pt);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};
template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_1_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_1_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type numPoints      = outputData.dimension(1);
               for(index_type pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(pt)/inputDataLeft(cl, 0);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_1_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_1_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  int numPoints      = outputData.dimension(1);
               for(int pt = 0; pt < numPoints; pt++) {
                  outputData(cl, pt) = inputDataRight(pt)*inputDataLeft(cl, 0);
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_2_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_2_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;

  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
 
    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                   outputData(cl, pt, iVec) = inputDataRight(pt, iVec)/inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_2_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_2_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {

    }
  }
              for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iVec = 0; iVec < dim1Tens; iVec++) {
                    outputData(cl, pt, iVec) = inputDataRight(pt, iVec)*inputDataLeft(cl, 0);
                } // D1-loop
              } // P-loop
   }
};


template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_recip_invalr_3_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_recip_invalr_3_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                      outputData(cl, pt, iTens1, iTens2) = inputDataRight(pt, iTens1, iTens2)/inputDataLeft(cl, 0);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};

template <class Scalar,class ArrayOutDataWrap,class ArrayInDataLeftWrap,class ArrayInDataRightWrap,class ArrayInDataLeft>
struct scalarMultiplyDataData_neqrank_nrecip_invalr_3_const {
  ArrayOutDataWrap outputData;
  ArrayInDataLeftWrap inputDataLeft;
  ArrayInDataRightWrap inputDataRight;
typedef typename conditional_eSpace<ArrayInDataLeft>::execution_space execution_space;
  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  scalarMultiplyDataData_neqrank_nrecip_invalr_3_const (ArrayOutDataWrap outputData_, ArrayInDataLeftWrap inputDataLeft_,ArrayInDataRightWrap inputDataRight_) :
    outputData (outputData_),inputDataLeft (inputDataLeft_),inputDataRight(inputDataRight_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const index_type cl) const {
  index_type outvalRank     = getrank(outputData);
  index_type numPoints      = outputData.dimension(1);
  index_type dim1Tens       = 0;
  index_type dim2Tens       = 0;
  if (outvalRank > 2) {
    dim1Tens = outputData.dimension(2);
    if (outvalRank > 3) {
      dim2Tens = outputData.dimension(3);
    }
  }
               for(index_type pt = 0; pt < numPoints; pt++) {
                for( index_type iTens1 = 0; iTens1 < dim1Tens; iTens1++) {
                  for( index_type iTens2 = 0; iTens2 < dim2Tens; iTens2++) {
                     outputData(cl, pt, iTens1, iTens2) = inputDataRight(pt, iTens1, iTens2)*inputDataLeft(cl, 0);
                  } // D2-loop
                } // D1-loop
              } // P-loop
   }
};
template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
void ArrayTools::scalarMultiplyDataData(ArrayOutData &           outputData,
                                        const ArrayInDataLeft &        inputDataLeft,
                                        const ArrayInDataRight &       inputDataRight,
                                        const bool               reciprocal) {

#ifdef HAVE_INTREPID2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inputDataLeft) != 2), std::invalid_argument,
			      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
  if (getrank(outputData) <= getrank(inputDataRight)) {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputDataRight) < 2) || (getrank(inputDataRight) > 4) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2, 3, or 4.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != getrank(inputDataRight)), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != inputDataLeft.dimension(0) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first dimension of the left data input container must be 1!");
    for (index_type i=0; i<getrank(inputDataRight); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimension ";
      errmsg += (char)(48+i);
      errmsg += " of the right input and output data containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( ( (getrank(inputDataRight) < 1) || (getrank(inputDataRight) > 3) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1, 2, or 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(outputData) != getrank(inputDataRight)+1), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
    TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
    TEUCHOS_TEST_FOR_EXCEPTION( ( inputDataLeft.dimension(0) != outputData.dimension(0) ), std::invalid_argument,
				">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
    for (index_type i=0; i<getrank(inputDataRight); i++) {
      std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the right input and output data containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i+1)), std::invalid_argument, errmsg );
    }
  }
#endif


   ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);
   ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);
   ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>inputDataRightWrap(inputDataRight);										


  // get sizes
  index_type invalRank      = getrank(inputDataRight);
  index_type outvalRank     = getrank(outputData);
  int numCells       = outputData.dimension(0);
  int numDataPoints  = inputDataLeft.dimension(1);

  if (outvalRank == invalRank) {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 2: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
        }// case 2
}
        break;

        case 3: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_4<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft>(outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 containers.");
      }// invalRank

    }
    else { // constant left data

      switch(invalRank) {
        case 2: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 2
        break;

        case 3: {
          if (reciprocal) {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 3
        break;

        case 4: {
          if (reciprocal) {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_recip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_eqrank_nrecip_invalr_4_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 4
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 2) || (invalRank == 3) || (invalRank == 4) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-2, 3 or 4 input containers.");

      } // invalRank
    } // numDataPoints

  }
  else {

    if (numDataPoints != 1) { // nonconstant data

      switch(invalRank) {
        case 1: {
if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_1<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_1<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 1
        break;

        case 2: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_2<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 2
        break;

        case 3: {
          if (reciprocal) {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
  Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_3<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 3
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 1) || (invalRank == 2) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 input containers.");
      }// invalRank
	
    }
    else { //constant data

      switch(invalRank) {
        case 1: {
          if (reciprocal) {
 Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_1_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_1_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 1
        break;

        case 2: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_2_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 2
        break;

        case 3: {
          if (reciprocal) {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_recip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
          else {
Kokkos::parallel_for (numCells, scalarMultiplyDataData_neqrank_nrecip_invalr_3_const<Scalar,ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>, ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>, ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value,true>, ArrayInDataLeft> (outputDataWrap,inputDataLeftWrap,inputDataRightWrap));
          }
        }// case 3
        break;

        default:
              TEUCHOS_TEST_FOR_EXCEPTION( !( (invalRank == 1) || (invalRank == 2) || (invalRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): This branch of the method is defined only for rank-1, 2 or 3 input containers.");

      } // invalRank
    } // numDataPoints

  } // end if (outvalRank = invalRank)

} // scalarMultiplyDataData
} // end namespace Intrepid2
#endif

