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

/** \file   Intrepid2_ArrayToolsDefScalar.hpp
    \brief  Definition file for scalar multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_SCALAR_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_SCALAR_HPP__

namespace Intrepid2 {

    namespace FunctorArrayTools {
    template < typename outputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_scalarMultiply{
      outputViewType _output;
      leftInputViewType _leftInput;
      rightInputViewType _rightInput;
      const bool _hasField, _reciprocal;

      KOKKOS_INLINE_FUNCTION
      F_scalarMultiply(outputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_,
              const bool reciprocal_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_), 
          _hasField(hasField_), _reciprocal(reciprocal_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl = 0, bf = 0, pt = 0;
        size_type outputRank(_output.rank()), rightRank(_rightInput.rank());

        if (_hasField) 
          Util::unrollIndex( cl, bf, pt, 
                       _output.dimension(0),
                       _output.dimension(1),
                       iter );
        else          
          Util::unrollIndex( cl, pt,
                       _output.dimension(0),
                       iter );

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_output, cl,     pt, Kokkos::ALL(), Kokkos::ALL()));
        
        const auto right = ( outputRank == rightRank ? ( _hasField ? Kokkos::subdynrankview(_rightInput, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                                                         /**/        Kokkos::subdynrankview(_rightInput, cl,     pt, Kokkos::ALL(), Kokkos::ALL()) ) : 
                             /**/                      ( _hasField ? Kokkos::subdynrankview(_rightInput,     bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                                                         /**/        Kokkos::subdynrankview(_rightInput,         pt, Kokkos::ALL(), Kokkos::ALL()) ) );
                                                                
        const auto left = (_leftInput.dimension(1) == 1) ? Kokkos::subdynrankview(_leftInput, cl, 0) :
                            /**/                           Kokkos::subdynrankview(_leftInput, cl, pt);

        const size_type iend  = result.dimension(0);
        const size_type jend  = result.dimension(1);

        const auto val = left();
        if (_reciprocal)
          for(size_type i = 0; i < iend; ++i)
            for(size_type j = 0; j < jend; ++j)
              result(i, j) = right(i, j)/val;
        else
          for(size_type i = 0; i < iend; ++i)
            for(size_type j = 0; j < jend; ++j)
              result(i, j) = right(i, j)*val;
      }
    };
    } // namespace

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::
  scalarMultiply( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>     output,
                  const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>  leftInput,
                  const Kokkos::DynRankView<rightInputValueType,rightInputProperties...> rightInput,
                  const bool hasField,
                  const bool reciprocal ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>     outputViewType;
    typedef Kokkos::DynRankView<leftInputValueType, leftInputProperties...>  leftInputViewType;
    typedef Kokkos::DynRankView<rightInputValueType,rightInputProperties...> rightInputViewType;
    typedef FunctorArrayTools::F_scalarMultiply<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                 /**/       output.dimension(0)*output.dimension(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, reciprocal) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
  scalarMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool reciprocal ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 2, std::invalid_argument,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input data container must have rank 2.");

      if (outputFields.rank() <= inputFields.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 3, 4, or 5.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != inputFields.rank(), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input and output fields containers must have the same rank.");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(2) &&
                                  inputData.dimension(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree, or first data dimension must be 1!");
        for (size_t i=0;i<inputFields.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(i) != outputFields.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i) of outputFields");
        }
      }
      else {
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 2 || inputFields.rank() > 4, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 2, 3, or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != (inputFields.rank()+1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): The rank of the input fields container must be one less than the rank of the output fields container.");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(1) &&
                                  inputData.dimension(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(0) != outputFields.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
        for (size_type i=0;i<inputFields.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(i) != outputFields.dimension(i+1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
        }
      }
    }
#endif
    // body
    ArrayTools<ExecSpaceType>::Internal::scalarMultiply( outputFields,
                                                         inputData,
                                                         inputFields,
                                                         true,
                                                         reciprocal );
  }



  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<ExecSpaceType>::
  scalarMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool reciprocal ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() != 2, std::invalid_argument,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");

      if (outputData.rank() <= inputDataRight.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 4, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2, 3, or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != inputDataRight.rank(), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(0) != inputDataLeft.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1) &&
                                  inputDataLeft.dimension(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree, or first dimension of the left data input container must be 1!");
        for (size_type i=0;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(i) != outputData.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
        }
      } else {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 1 || inputDataRight.rank() > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1, 2, or 3.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != (inputDataRight.rank()+1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(0)  &&
                                  inputDataLeft.dimension(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != outputData.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
        for (size_type i=0;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(i) != outputData.dimension(i+1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i+1) of outputData");
        }
      }
    }
#endif
    // body
    ArrayTools<ExecSpaceType>::Internal::scalarMultiply( outputData,
                                                         inputDataLeft,
                                                         inputDataRight,
                                                         false,
                                                         reciprocal );
  }

} // end namespace Intrepid2

#endif
