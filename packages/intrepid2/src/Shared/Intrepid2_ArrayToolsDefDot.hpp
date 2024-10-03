// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayToolsDefDot.hpp
    \brief  Definition file for dot-multiply operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_DOT_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_DOT_HPP__

namespace Intrepid2 {

    namespace FunctorArrayTools {
    /**
       \brief Functor for dotMultiply see Intrepid2::ArrayTools for more
    */
    template < typename OutputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_dotMultiply {
      OutputViewType _output;
      leftInputViewType _leftInput;
      rightInputViewType _rightInput;
      const bool _hasField;
      typedef typename OutputViewType::value_type value_type;

      KOKKOS_INLINE_FUNCTION
      F_dotMultiply(OutputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl(0), bf(0), pt(0);
        size_type leftRank(_leftInput.rank()), rightRank(_rightInput.rank());

        if (_hasField) 
          unrollIndex( cl, bf, pt, 
                             _output.extent(0),
                             _output.extent(1), 
                             _output.extent(2), 
                             iter );
        else          
          unrollIndex( cl, pt,
                             _output.extent(0),
                             _output.extent(1),
                             iter);
        
        auto result = ( _hasField ? Kokkos::subview(_output, cl, bf, pt) :
                                    Kokkos::subview(_output, cl,     pt));
        
        const auto left = (_leftInput.extent(1) == 1) ? Kokkos::subview(_leftInput, cl, 0, Kokkos::ALL(), Kokkos::ALL()) :
                                                           Kokkos::subview(_leftInput, cl, pt, Kokkos::ALL(), Kokkos::ALL());

        
        const auto right = (rightRank == leftRank + ordinal_type(_hasField)) ?
                             ( _hasField ? Kokkos::subview(_rightInput, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                                           Kokkos::subview(_rightInput, cl,     pt, Kokkos::ALL(), Kokkos::ALL())) :
                             ( _hasField ? Kokkos::subview(_rightInput,     bf, pt, Kokkos::ALL(), Kokkos::ALL()) :
                                           Kokkos::subview(_rightInput,         pt, Kokkos::ALL(), Kokkos::ALL()));
        
        const ordinal_type iend  = left.extent(0);
        const ordinal_type jend  = left.extent(1);

        value_type tmp(0);
        for(ordinal_type i = 0; i < iend; ++i)
          for(ordinal_type j = 0; j < jend; ++j)
            tmp += left(i, j)*right(i, j);
        result() = tmp;
      }
    };
    } //namespace

  template<typename DeviceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<DeviceType>::Internal::
  dotMultiply(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
               const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
               const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput, 
               const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      OutputViewType;
    typedef Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_dotMultiply<OutputViewType, leftInputViewType, rightInputViewType> FunctorType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField) );
  }



  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
  dotMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                        const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                        const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      if (inputFields.rank() > inputData.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() < 2 || inputData.rank() > 4, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() != (inputData.rank()+1), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input fields container must have rank one larger than the rank of the input data container.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 3, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2) &&
                                        inputData.extent(1) != 1, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
        for (size_type i=2;i<inputData.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(i) != inputFields.extent(i+1), std::invalid_argument,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i+1) of inputFields");
        }
        for (size_type i=0;i<outputFields.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(i) != outputFields.extent(i), std::invalid_argument,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
        }
      } else {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() < 2 || inputData.rank() > 4, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() != inputData.rank(), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): The rank of fields input container must equal the rank of data input container.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 3, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(1) &&
                                        inputData.extent(1) != 1, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimensions of the fields and data input containers (number of integration points) must agree or first data dimension must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != outputFields.extent(1), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the fields input container and first dimension of the fields output container (number of fields) must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(1) != outputFields.extent(2), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimension of the fields input container and second dimension of the fields output container (number of integration points) must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions of the fields output and data input containers (number of integration domains) must agree!");
        for (size_type i=2;i<inputData.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(i) != inputFields.extent(i), std::invalid_argument,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i) of inputFields");
        }
      }
    }
#endif

    ArrayTools<DeviceType>::Internal::dotMultiply( outputFields,
                                                      inputData,
                                                      inputFields,
                                                      true );
  }



  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
  dotMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                       const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      if (inputDataRight.rank() >= inputDataLeft.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() != inputDataLeft.rank(), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): The rank of the right data input container must equal the rank of the left data input container.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 2, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1) &&
                                        inputDataLeft.extent(1) != 1, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first left data dimension must be 1!");
        for (size_type i=0;i<inputDataLeft.rank();++i) {
          if (i != 1) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(i) != inputDataRight.extent(i), std::invalid_argument,
                                            ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i) does not match to the dimension (i) of inputDataRight");
          }
        }
        for (size_type i=0;i<outputData.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(i) != outputData.extent(i), std::invalid_argument,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
        }
      } else {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() != (inputDataLeft.rank()-1), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Right data input container must have rank one less than the rank of left data input container.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 2, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(0) &&
                                        inputDataLeft.extent(1) != 1, std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of left data input container (number of integration points) must agree or first left data dimension must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(0) != outputData.extent(1), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of output data container (number of integration points) must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != outputData.extent(0), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimensions of the left data input and data output containers (number of integration domains) must agree!");
        for (size_type i=1;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(i+1) != inputDataRight.extent(i), std::invalid_argument,
                                          ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i+1) does not match to the dimension (i) of inputDataRight");
        }
      }
    }
#endif

    ArrayTools<DeviceType>::Internal::dotMultiply( outputData,
                                                      inputDataLeft,
                                                      inputDataRight,
                                                      false );
  }

} // end namespace Intrepid2
#endif
