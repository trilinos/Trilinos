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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayToolsDefDot.hpp
    \brief  Definition file for dot-multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_DOT_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_DOT_HPP__

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::Internal::
  dotMultiply( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
               const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
               const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput, 
               const bool hasField ) {
    typedef leftInputValueType value_type;

    struct Functor {
      Kokkos::DynRankView<outputValueType,    outputProperties...>     _output;
      Kokkos::DynRankView<leftInputValueType, leftInputProperties...>  _leftInput;
      Kokkos::DynRankView<rightInputValueType,rightInputProperties...> _rightInput;
      const bool _hasField;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputValueType,    outputProperties...>     output_,
              Kokkos::DynRankView<leftInputValueType, leftInputProperties...>  leftInput_,
              Kokkos::DynRankView<rightInputValueType,rightInputProperties...> rightInput_,
              const bool hasField_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, bf, pt;
        size_type leftRank(_leftInput.rank()), rightRank(_rightInput.rank());

        if (_hasField) 
          Util::unrollIndex( cl, bf, pt, 
                       _output.dimension(0),
                       _output.dimension(1), 
                       iter );
        else          
          Util::unrollIndex( cl, pt,
                       _output.dimension(0),
                       iter);

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cl, bf, pt) :
                        /**/        Kokkos::subdynrankview(_output, cl,     pt));
        
        const auto left = (_leftInput.dimension(1) == 1) ? Kokkos::subdynrankview(_leftInput, cl, 0, Kokkos::ALL(), Kokkos::ALL()) :
                            /**/                           Kokkos::subdynrankview(_leftInput, cl, pt, Kokkos::ALL(), Kokkos::ALL());

        
        const auto right = (rightRank == leftRank + int(_hasField)) ?
                             ( _hasField ? Kokkos::subdynrankview(_rightInput, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                             /**/          Kokkos::subdynrankview(_rightInput, cl,     pt, Kokkos::ALL(), Kokkos::ALL())) :
                             ( _hasField ? Kokkos::subdynrankview(_rightInput,     bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                             /**/          Kokkos::subdynrankview(_rightInput,         pt, Kokkos::ALL(), Kokkos::ALL()));
        
        const size_type iend  = left.dimension(0);
        const size_type jend  = left.dimension(1);

        value_type tmp(0);
        for(size_type i = 0; i < iend; ++i)
          for(size_type j = 0; j < jend; ++j)
            tmp += left(i, j)*right(i, j);
        result() = tmp;
      }
    };

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                 /**/       output.dimension(0)*output.dimension(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(output, leftInput, rightInput, hasField) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  dotMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                        const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                        const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      if (inputFields.rank() > inputData.rank()) {
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.rank() < 2 || inputData.rank() > 4, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.rank() != (inputData.rank()+1), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input fields container must have rank one larger than the rank of the input data container.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 3, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(0) != inputData.dimension(0), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                                        inputData.dimension(1) != 1, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
        for (size_type i=2;i<inputData.rank();++i) {
          INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(i) != inputFields.dimension(i+1), dbgInfo,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i+1) of inputFields");
        }
        for (size_t i=0;i<outputFields.rank();++i) {
          INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(i) != outputFields.dimension(i), dbgInfo,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
        }
      } else {
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.rank() < 2 || inputData.rank() > 4, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.rank() != inputData.rank(), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): The rank of fields input container must equal the rank of data input container.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 3, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(1) != inputFields.dimension(1) &&
                                        inputData.dimension(1) != 1, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimensions of the fields and data input containers (number of integration points) must agree or first data dimension must be 1!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(0) != outputFields.dimension(1), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the fields input container and first dimension of the fields output container (number of fields) must agree!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(1) != outputFields.dimension(2), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimension of the fields input container and second dimension of the fields output container (number of integration points) must agree!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != inputData.dimension(0), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions of the fields output and data input containers (number of integration domains) must agree!");
        for (size_type i=2;i<inputData.rank();++i) {
          INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(i) != inputFields.dimension(i+1), dbgInfo,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i+1) of inputFields");
        }
      }

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::dotMultiply( outputFields,
                                                      inputData,
                                                      inputFields,
                                                      true );
  }



  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  dotMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                       const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      if (inputDataRight.rank() >= inputDataLeft.rank()) {
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.rank() != inputDataLeft.rank(), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): The rank of the right data input container must equal the rank of the left data input container.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.rank() != 2, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1) &&
                                        inputDataLeft.dimension(1) != 1, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first left data dimension must be 1!");
        for (size_type i=0;i<inputDataLeft.rank();++i) {
          if (i != 1) {
            INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(i) != inputDataRight.dimension(i), dbgInfo,
                                            ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i) does not match to the dimension (i) of inputDataRight");
          }
        }
        for (size_type i=0;i<outputData.rank();++i) {
          INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.dimension(i) != outputData.dimension(i), dbgInfo,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
        }
      } else {
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.rank() != (inputDataLeft()-1), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Right data input container must have rank one less than the rank of left data input container.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.rank() != 2,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(0) &&
                                        inputDataLeft.dimension(1) != 1, dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of left data input container (number of integration points) must agree or first left data dimension must be 1!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.dimension(0) != outputData.dimension(1), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of output data container (number of integration points) must agree!");
        INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(0) != outputData.dimension(0), dbgInfo,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimensions of the left data input and data output containers (number of integration domains) must agree!");
        for (size_type i=1;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(i+1) != inputDataRight.dimension(i), dbgInfo,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i+1) does not match to the dimension (i) of inputDataRight");
        }
      }

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::dotMultiply( outputData,
                                                      inputDataLeft,
                                                      inputDataRight,
                                                      false );
  }

} // end namespace Intrepid2
#endif
