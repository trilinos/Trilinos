// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayToolsDefTensor.hpp
    \brief  Definition file for tensor multiply operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__

namespace Intrepid2 {

    namespace FunctorArrayTools {
    /**
       \brief Functor for crossProduct see Intrepid2::ArrayTools for more
    */
    template < typename OutputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_crossProduct{
      OutputViewType _output;
      const leftInputViewType _leftInput;
      const rightInputViewType _rightInput;
      const bool _hasField, _isCrossProd3D;

      KOKKOS_INLINE_FUNCTION
      F_crossProduct(OutputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_,
              const bool isCrossProd3D_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_), _isCrossProd3D(isCrossProd3D_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cell, field(0), point;
        size_type rightRank = _rightInput.rank();

        if (_hasField) 
          unrollIndex( cell, field, point, 
                             _output.extent(0),
                             _output.extent(1),
                             _output.extent(2),
                             iter );
        else           
          unrollIndex( cell, point, 
                             _output.extent(0), 
                             _output.extent(1), 
                             iter );

        auto result = ( _hasField ? Kokkos::subview(_output, cell, field, point, Kokkos::ALL()) :
                                    Kokkos::subview(_output, cell,        point, Kokkos::ALL()));

        auto left   = Kokkos::subview(_leftInput, cell, point, Kokkos::ALL());

        auto right  = (rightRank == 2 + size_type(_hasField)) ?
                      ( _hasField ? Kokkos::subview(_rightInput,       field, point, Kokkos::ALL()) :
                                    Kokkos::subview(_rightInput,              point, Kokkos::ALL())) :
                      ( _hasField ? Kokkos::subview(_rightInput, cell, field, point, Kokkos::ALL()) :
                                    Kokkos::subview(_rightInput, cell,        point, Kokkos::ALL()));

        // branch prediction is possible
        if (_isCrossProd3D) {
          result(0) = left(1)*right(2) - left(2)*right(1);
          result(1) = left(2)*right(0) - left(0)*right(2);
          result(2) = left(0)*right(1) - left(1)*right(0);
        } else {
          result(0) = left(0)*right(1) - left(1)*right(0);
        }
      }
    };
    } //namespace

  template<typename DeviceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<DeviceType>::Internal::
  crossProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      OutputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_crossProduct<OutputViewType, leftInputViewType, rightInputViewType> FunctorType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    const bool isCrossProd3D = (leftInput.extent(2) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isCrossProd3D) );
  }



  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
  crossProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputData(C,P,D);
       *      (2) inputFields(C,F,P,D) or (F,P,D);
       *      (3) outputFields(C,F,P,D) in 3D, or (C,F,P) in 2D
       */
      // (1) inputData is (C, P, D) and 2 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputData must have rank 3");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(2) < 2 || inputData.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension(2) must be 2 or 3");

      // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputFields must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(inputFields.rank()-1) < 2 ||
                                inputFields.extent(inputFields.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

      // (3) outputFields is (C,F,P,D) in 3D and (C,F,P) in 2D => rank = inputData.extent(2) + 1
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != inputData.extent(2)+1, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): outputFields rank must match to inputData dimension(2)+1");
      /*
       *   Dimension cross-checks:
       *      (1) inputData    vs. inputFields
       *      (2) outputFields vs. inputData
       *      (3) outputFields vs. inputFields
       *
       *   Cross-check (1):
       */
      if (inputFields.rank() == 4) {
        // inputData(C,P,D) vs. inputFields(C,F,P,D): dimensions C, P, D must match
        const size_type f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
        }
      } else {
        // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
        const size_type f1[] = { 1, 2 }, f2[] = { 1, 2 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
        }
      }
      /*
       *  Cross-check (2):
       */
      if (inputData.extent(2) == 2) {
        //  in 2D: outputFields(C,F,P) vs. inputData(C,P,D): dimensions C,P must match
        // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
        const size_type f1[] = { 0, 2 }, f2[] = { 0, 1 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputData");
        }
      } else {
        // in 3D: outputFields(C,F,P,D) vs. inputData(C,P,D): dimensions C,P,D must match
        const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 1, 2 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputData");
        }
      }
      /*
       *  Cross-check (3):
       */
      if (inputData.extent(2) == 2) {
        // In 2D:
        if (inputFields.rank() == 4) {
          //  and rank-4 inputFields: outputFields(C,F,P) vs. inputFields(C,F,P,D): dimensions C,F,P must match
          const size_type f1[] = { 0, 1, 2 }, f2[] = { 0, 1, 2 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        } else {
          //  and rank-3 inputFields: outputFields(C,F,P) vs. inputFields(F,P,D): dimensions F,P must match
          const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      } else {
        // In 3D:
        if (inputFields.rank() == 4) {
          //  and rank-4 inputFields: outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C,F,P,D must match
          for (size_type i=0; i<outputFields.rank(); ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(i) != inputFields.extent(i), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        } else {
          // and rank-3 inputFields: outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F,P,D must match
          const size_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::crossProduct( outputFields,
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
  crossProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputDataLeft(C,P,D);
       *      (2) inputDataRight(C,P,D) or (P,D);
       *      (3) outputData(C,P,D) in 3D, or (C,P) in 2D
       */
      // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft must have rank 3");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) < 2 || inputDataLeft.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension(2) must be 2 or 3");

      // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(inputDataRight.rank()-1) < 2 ||
                                inputDataRight.extent(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight dimension (rank-1) must have rank 2 or 3" );

      // (3) outputData is (C,P,D) in 3D and (C,P) in 2D => rank = inputDataLeft.extent(2)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != inputDataLeft.extent(2), std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): outputData rank must match to inputDataLeft dimension(2)");
      /*
       *   Dimension cross-checks:
       *      (1) inputDataLeft vs. inputDataRight
       *      (2) outputData    vs. inputDataLeft
       *      (3) outputData    vs. inputDataRight
       *
       *   Cross-check (1):
       */
      if (inputDataRight.rank() == 3) {
        // inputDataLeft(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
        for (size_type i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(i) != inputDataRight.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
        }
      }
      // inputDataLeft(C, P,D) vs. inputDataRight(P,D): dimensions P, D must match
      else {
        const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
        }
      }
      /*
       *  Cross-check (2):
       */
      if (inputDataLeft.extent(2) == 2) {
        // in 2D: outputData(C,P) vs. inputDataLeft(C,P,D): dimensions C, P must match
        const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != outputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to outputData");
        }
      } else {
        // in 3D: outputData(C,P,D) vs. inputDataLeft(C,P,D): all dimensions C, P, D must match
        for (size_type i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(i) != outputData.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to outputData");
        }
      }
      /*
       *  Cross-check (3):
       */
      if (inputDataLeft.extent(2) == 2) {
        // In 2D:
        if (inputDataRight.rank() == 3) {
          //  and rank-3 inputDataRight: outputData(C,P) vs. inputDataRight(C,P,D): dimensions C,P must match
          const size_type f1[] = { 0, 1 }, f2[] = { 0, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
          }
        } else {
          //  and rank-2 inputDataRight: outputData(C,P) vs. inputDataRight(P,D): dimension P must match
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(1) != inputDataRight.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
        }
      } else {
        // In 3D:
        if (inputDataRight.rank() == 3) {
          //  and rank-3 inputDataRight: outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C,P,D must match
          for (size_type i=0; i<outputData.rank(); ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(i) != inputDataRight.extent(i), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
          }
        } else {
          //  and rank-2 inputDataRight: outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
          const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::crossProduct( outputData,
                                                   inputDataLeft,
                                                   inputDataRight,
                                                   false );
  }


    namespace FunctorArrayTools {
    /**
       \brief Functor for outerProduct see Intrepid2::ArrayTools for more
    */
    template < typename OutputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_outerProduct {
      OutputViewType _output;
      const leftInputViewType _leftInput;
      const rightInputViewType _rightInput;
      const bool _hasField;

      KOKKOS_INLINE_FUNCTION
      F_outerProduct(OutputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cell, field(0), point;
        size_type rightRank = _rightInput.rank();

        if (_hasField) 
          unrollIndex( cell, field, point, 
                             _output.extent(0),
                             _output.extent(1),
                             _output.extent(2),
                             iter );
        else           
          unrollIndex( cell, point, 
                             _output.extent(0), 
                             _output.extent(1), 
                             iter );
        
        auto result = ( _hasField ? Kokkos::subview(_output, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                                    Kokkos::subview(_output, cell,        point, Kokkos::ALL(), Kokkos::ALL()));

        auto left   = Kokkos::subview(_leftInput, cell, point, Kokkos::ALL());

        auto right  = (rightRank == 2 + size_type(_hasField)) ?
                      ( _hasField ? Kokkos::subview(_rightInput,       field, point, Kokkos::ALL()) :
                                    Kokkos::subview(_rightInput,              point, Kokkos::ALL())) :
                      ( _hasField ? Kokkos::subview(_rightInput, cell, field, point, Kokkos::ALL()) :
                                    Kokkos::subview(_rightInput, cell,        point, Kokkos::ALL()));

        const ordinal_type iend = result.extent(0);
        const ordinal_type jend = result.extent(1);
        for (ordinal_type i=0; i<iend; ++i)
          for (ordinal_type j=0; j<jend; ++j)
            result(i, j) = left(i)*right(j);
      }
    };
    }

  template<typename DeviceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<DeviceType>::Internal::
  outerProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      OutputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_outerProduct<OutputViewType, leftInputViewType, rightInputViewType> FunctorType;

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
  outerProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputData(C,P,D);
       *      (2) inputFields(C,F,P,D) or (F,P,D);
       *      (3) outputFields(C,F,P,D,D)
       */
      // (1) inputData is (C, P, D) and 2 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputData must have rank 3");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(2) < 2 || inputData.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension(2) must be 2 or 3");

      // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputFields must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(inputFields.rank()-1) < 2 ||
                                inputFields.extent(inputFields.rank()-1) < 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

      // (3) outputFields is (C,F,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputFields must have rank 5");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(3) < 2 ||
                                outputFields.extent(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension(3) must be 2 or 3");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(4) < 2 ||
                                outputFields.extent(4) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension(4) must be 2 or 3");

      /*
       *   Dimension cross-checks:
       *      (1) inputData    vs. inputFields
       *      (2) outputFields vs. inputData
       *      (3) outputFields vs. inputFields
       *
       *   Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P,D): dimensions C, P, D must match
       */
      {
        const size_type f1[] = { 0, 2, 3, 4 }, f2[] = { 0, 1, 2, 2 };
        for (size_type i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputData");
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 4) {
        // Cross-check (1): inputData(C,P,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        {
          const size_type f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D): dimensions C, F, P, D must match
        {
          const size_type f1[] = { 0, 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3, 3 };
          for (size_type i=0; i<5; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      } else {
        // Cross-check (1): inputData(C,P,D) vs. inputFields(F,P,D): dimensions  P, D must match
        {
          const size_type f1[] = { 1, 2 }, f2[] = { 1, 2 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        {
          const size_type f1[] = { 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 2 };
          for (size_type i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::outerProduct(  outputFields,
                                                                 inputData,
                                                                 inputFields,
                                                                 true );
  }


  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValuetype,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
  outerProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftValuetype, inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputDataLeft(C,P,D);
       *      (2) inputDataRight(C,P,D) or (P,D);
       *      (3) outputData(C,P,D,D)
       */
      // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft must have rank 3");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) < 2 || inputDataLeft.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension(2) must be 2 or 3");

      // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(inputDataRight.rank()-1) < 2 ||
                                inputDataRight.extent(inputDataRight.rank()-1) < 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight dimension (rank-1) must have rank 2 or 3" );

      // (3) outputData is (C,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputData must have rank 5");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(2) < 2 ||
                                outputData.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension(3) must be 2 or 3");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(3) < 2 ||
                                outputData.extent(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension(4) must be 2 or 3");

      /*
       *   Dimension cross-checks:
       *      (1) inputDataLeft vs. inputDataRight
       *      (2) outputData    vs. inputDataLeft
       *      (3) outputData    vs. inputDataRight
       *
       *   Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P,D): dimensions C, P, D must match
       */
      {
        const size_type f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
        for (size_type i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataLeft.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataLeft");
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 3) {
        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(C,P,D):  all dimensions  C, P, D must match
        for (size_type i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(i) != inputDataRight.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D): dimensions C, P, D must match
        {
          const size_type f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
          for (size_type i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataRight");
          }
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(P,D): dimensions  P, D must match
        {
          const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
          }
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D): dimensions P, D must match
        {
          const size_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 1 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataRight");
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::outerProduct( outputData,
                                                               inputDataLeft,
                                                               inputDataRight,
                                                               false );
  }


  namespace FunctorArrayTools {
  /**
     \brief Functor for matvecProduct; new version avoids both subviews and branching.  See Intrepid2::ArrayTools for more.
  */
  template < typename OutputViewType,
             typename leftInputViewType,
             typename rightInputViewType,
             ordinal_type leftInputRank,
             ordinal_type rightInputRank,
             bool hasField,
             bool isTranspose>
  struct F_matvecProduct {
    /**/  OutputViewType     _output;
    const leftInputViewType  _leftInput;
    const rightInputViewType _rightInput;
    
    const ordinal_type _iend;
    const ordinal_type _jend;
    
    using value_type = typename OutputViewType::value_type;
    
    KOKKOS_INLINE_FUNCTION
    F_matvecProduct(OutputViewType     output_,
                        leftInputViewType  leftInput_,
                        rightInputViewType rightInput_)
      : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
        _iend(output_.extent_int(rank(output_)-1)), _jend(rightInput_.extent_int(rightInputRank-1))
    {}
    
    //  ****** hasField == true cases ******
    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_type cl,
                    const ordinal_type bf,
                    const ordinal_type pt) const
    {
      apply_matvec_product(cl, bf, pt);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==4 && r==4 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      if (isTranspose) {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, j,i)*_rightInput(cl,bf,pt, j);
          _output(cl,bf,pt, i) = tmp;
        }
      } else {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, i,j)*_rightInput(cl,bf,pt, j);
          _output(cl,bf,pt, i) = tmp;
        }
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==4 && r==3 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      if (isTranspose) {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, j,i)*_rightInput(bf,pt, j);
          _output(cl,bf,pt, i) = tmp;
        }
      } else {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, i,j)*_rightInput(bf,pt, j);
          _output(cl,bf,pt, i) = tmp;
        }
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==3 && r==4 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      for (ordinal_type i=0;i<_iend;++i)
        _output(cl,bf,pt, i) = _leftInput(cl,lpt, i)*_rightInput(cl,bf,pt, i);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==3 && r==3 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      for (ordinal_type i=0;i<_iend;++i)
        _output(cl,bf,pt, i) = _leftInput(cl,lpt, i)*_rightInput(bf,pt, i);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==2 && r==4 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      const value_type & val = _leftInput(cl,lpt);
      for (ordinal_type i=0;i<_iend;++i) {
        _output(cl,bf,pt, i) = val*_rightInput(cl,bf,pt, i);
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==2 && r==3 && hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &bf,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      const value_type & val = _leftInput(cl,lpt);
      for (ordinal_type i=0;i<_iend;++i) {
        _output(cl,bf,pt, i) = val*_rightInput(bf,pt, i);
      }
    }
    
    //  ****** hasField == false cases ******
    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_type cl,
                    const ordinal_type pt) const
    {
      apply_matvec_product(cl, pt);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==4 && r==3 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      if (isTranspose) {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, j,i)*_rightInput(cl,pt, j);
          _output(cl,pt, i) = tmp;
        }
      } else {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, i,j)*_rightInput(cl,pt, j);
          _output(cl,pt, i) = tmp;
        }
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==4 && r==2 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      if (isTranspose) {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, j,i)*_rightInput(pt, j);
          _output(cl,pt, i) = tmp;
        }
      } else {
        for (ordinal_type i=0;i<_iend;++i) {
          value_type tmp(0);
          for (ordinal_type j=0;j<_jend;++j)
            tmp += _leftInput(cl,lpt, i,j)*_rightInput(pt, j);
          _output(cl,pt, i) = tmp;
        }
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==3 && r==3 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      for (ordinal_type i=0;i<_iend;++i)
        _output(cl,pt, i) = _leftInput(cl,lpt, i)*_rightInput(cl,pt, i);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==3 && r==2 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      for (ordinal_type i=0;i<_iend;++i)
        _output(cl,pt, i) = _leftInput(cl,lpt, i)*_rightInput(pt, i);
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==2 && r==3 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      const value_type & val = _leftInput(cl,lpt);
      for (ordinal_type i=0;i<_iend;++i) {
        _output(cl,pt, i) = val*_rightInput(cl,pt, i);
      }
    }
    
    template <ordinal_type l=leftInputRank, ordinal_type r=rightInputRank, bool hf=hasField>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if_t<l==2 && r==2 && !hf, void>
    apply_matvec_product(const ordinal_type   &cl,
                         const ordinal_type   &pt) const {
      const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
      const value_type & val = _leftInput(cl,lpt);
      for (ordinal_type i=0;i<_iend;++i) {
        _output(cl,pt, i) = val*_rightInput(pt, i);
      }
    }
  };
  } //namespace

  template<typename DeviceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<DeviceType>::Internal::
  matvecProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {
    
    using Output = Kokkos::DynRankView<outputValueType,    outputProperties...>;
    using Left   = const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>;
    using Right  = const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>;
    
    // FTNMAB: FunctorType with left rank N, right rank M, hasField = (A==T), isTranspose = (B==T)
    
    // hasField == true
    using FT44TT = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,4,  true,  true>;
    using FT43TT = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,3,  true,  true>;
    using FT44TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,4,  true, false>;
    using FT43TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,3,  true, false>;
    
    // for left rank 3, 2, isTranspose does not matter
    using FT34TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,3,4,  true, false>;
    using FT33TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,3,3,  true, false>;
    using FT24TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,2,4,  true, false>;
    using FT23TF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,2,3,  true, false>;
    
    // hasField == false
    using FT43FT = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,3, false,  true>;
    using FT42FT = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,2, false,  true>;
    using FT43FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,3, false, false>;
    using FT42FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,4,2, false, false>;
    
    // for left rank 3, 2, isTranspose does not matter
    using FT33FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,3,3, false, false>;
    using FT32FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,3,2, false, false>;
    using FT23FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,2,3, false, false>;
    using FT22FF = FunctorArrayTools::F_matvecProduct<Output,Left,Right,2,2, false, false>;

    const ordinal_type l = leftInput.rank();
    const ordinal_type r = rightInput.rank();
    
    using range_policy2 = Kokkos::MDRangePolicy< ExecSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    range_policy2 policy2( { 0, 0 }, { output.extent(0), output.extent(1) } );
    
    using range_policy3 = Kokkos::MDRangePolicy< ExecSpaceType, Kokkos::Rank<3>, Kokkos::IndexType<ordinal_type> >;
    range_policy3 policy3( { 0, 0, 0 }, { output.extent(0), output.extent(1), output.extent(2) } );
    
    // just to make the below a little easier to read, we pack l and r together:
    const ordinal_type lr = l * 10 + r;
    const auto &ov = output;
    const auto &lv = leftInput;
    const auto &rv = rightInput;
    
    if (hasField) // this means we want policy3
    {
      if (l == 4)
      {
        // isTranspose matters
        if (isTranspose)
        {
          switch (r)
          {
            case 4:  { FT44TT f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } break;
            default: { FT43TT f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); }
          }
        }
        else
        {
          switch (r)
          {
            case 4:  { FT44TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } break;
            default: { FT43TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); }
          }
        }
      }
      else // l == 3 or 2; isTranspose does not matter
      {
        switch (lr)
        {
          case 34: { FT34TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } break;
          case 33: { FT33TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } break;
          case 24: { FT24TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } break;
          default: { FT23TF f(ov,lv,rv); Kokkos::parallel_for( policy3, f ); } // 23
        }
      }
    }
    else // hasField is false; use policy2
    {
      if (l == 4)
      {
        // isTranspose matters
        if (isTranspose)
        {
          switch (r)
          {
            case 3:  { FT43FT f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } break;
            default: { FT42FT f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } // 2
          }
        }
        else
        {
          switch (r)
          {
            case 3:  { FT43FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } break;
            default: { FT42FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } // 2
          }
        }
      }
      else // l == 3 or 2; isTranspose does not matter
      {
        switch (lr)
        {
          case 33: { FT33FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } break;
          case 32: { FT32FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } break;
          case 23: { FT23FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } break;
          default: { FT22FF f(ov,lv,rv); Kokkos::parallel_for( policy2, f ); } // 22
        }
      }
    }
  }

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::matvecProductDataField(      Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields,
                          const char transpose ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
       *      (2) inputFields(C,F,P,D) or (F,P,D);
       *      (3) outputFields(C,F,P,D)
       */
      // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() < 2 || inputData.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): inputData must have rank 2,3 or 4" );
      if (inputData.rank() > 2) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(2) < 1 ||
                                  inputData.extent(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) must be 1,2 or 3");
      }
      if (inputData.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(3) < 1 ||
                                  inputData.extent(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(3) must be 1,2 or 3");
      }

      // (2) inputFields is (C, F, P, D) or (F, P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 ||
                                inputFields.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): inputFields must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( (inputFields.rank()-1) < 1 ||
                                (inputFields.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): inputFields dimension(rank-1) must be 1,2, or 3" );

      // (3) outputFields is (C,F,P,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): outputFields must have rank 4" );
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(3) < 1 ||
                                outputFields.extent(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(3) must be 1,2 or 3" );

      /*
       *   Dimension cross-checks:
       *      (1) inputData    vs. inputFields
       *      (2) outputFields vs. inputData
       *      (3) outputFields vs. inputFields
       *
       *   Cross-check (2): outputFields(C,F,P,D) vs. inputData(C,P), (C,P,D) or (C,P,D,D):
       *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
       *   data is specified (P>1). Do not check P dimensions with constant data, i.e., when P=1 in
       *   inputData(C,1,...)
       */
      if (inputData.extent(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(2) != inputData.extent(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(2) must match to inputData dimension(1)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(0) must match to inputData dimension(0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
        const size_type f1[] = { 0, 3 }, f2[] = { 0, 2 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }
      if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
        const size_type f1[] = { 0, 3, 3 }, f2[] = { 0, 2, 3 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 4) {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        if (inputData.extent(1) > 1) { // check P dimension if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(2) != inputData.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputFields dimension (2) does not match to inputData dimension (1)" );
        }
        if (inputData.rank() == 2) { // inputData(C,P) -> C match
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(0) != inputFields.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
        }
        if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
          const size_type f1[] = { 0, 2 }, f2[] = { 0, 3 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
          const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 3 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C, F, P, D must match
        for (size_type i=0; i<outputFields.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(i) != inputFields.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      } else {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D): dimensions  P, D must match
        if (inputData.extent(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(1) does not match to inputFields dimension (1)" );
        }
        if (inputData.rank() == 3) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(2) != inputFields.extent(2), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) does not match to inputFields dimension (2)" );
        }
        if (inputData.rank() == 4) {
          const size_type f1[] = { 2, 3 }, f2[] = { 2, 2 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        {
          const size_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputFields dimension" );
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::matvecProduct( outputFields,
                                                        inputData,
                                                        inputFields,
                                                        true,
                                                        transpose == 't' || transpose == 'T' );
  }



  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
  matvecProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>    outputData,
                         const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...> inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...>    inputDataRight,
                         const char transpose ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D1,D2); P=1 is admissible to allow multiply by const. left data
       *      (2) inputDataRight(C,P,D) or (P,D);
       *      (3) outputData(C,P,D)
       */
      // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D1,D2) and 1 <= D,D1,D2 <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() < 2 ||
                                inputDataLeft.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft must have rank 2,3 or 4" );

      if (inputDataLeft.rank() > 2) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) < 1 ||
                                  inputDataLeft.extent(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(2) must be 1, 2 or 3");
      }
      if (inputDataLeft.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(3) < 1 ||
                                  inputDataLeft.extent(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(3) must be 1, 2 or 3");
      }

      // (2) inputDataRight is (C, P, D) or (P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 ||
                                inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(inputDataRight.rank()-1) < 1 ||
                                inputDataRight.extent(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight dimension (rank-1) must be 1,2 or 3" );

      // (3) outputData is (C,P,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): outputData must have rank 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(2) < 1 ||
                                outputData.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(2) must be 1, 2 or 3");

      /*
       *   Dimension cross-checks:
       *      (1) inputDataLeft vs. inputDataRight
       *      (2) outputData    vs. inputDataLeft
       *      (3) outputData    vs. inputDataRight
       *
       *   Cross-check (2): outputData(C,P,D) vs. inputDataLeft(C,P), (C,P,D) or (C,P,D,D):
       *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
       *   data is specified (P>1). Do not check P dimensions with constant left data, i.e., when P=1 in
       *   inputDataLeft(C,1,...)
       */
      if (inputDataLeft.extent(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(1) != inputDataLeft.extent(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(1) muat match inputDataLeft dimension(1)" );
      }
      if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(0) != inputDataLeft.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(0) muat match inputDataLeft dimension(0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const size_type f1[] = { 0, 2 }, f2[] = { 0, 2 };
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataLeft.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D1,D2): check C and D
        size_type f1[] = { 0, 2}, f2[] = { 0, 2};
        if (transpose) f2[1] = 3;
        for (size_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataLeft.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 3) {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D):  dimensions  C, P, D must match
        if (inputDataLeft.extent(1) > 1) { // check P dimension if P>1 in inputDataLeft
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) muat match to inputDataRight dimension (1)" );
        }
        if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (0) muat match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
          const size_type f1[] = { 0, 2 }, f2[] = { 0, 2 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
          size_type f1[] = { 0, 3}, f2[] = { 0, 2};
          if (transpose) f1[1] = 2;
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P) vs. inputDataRight(C,P): all dimensions C, P must match
        for (size_type i=0; i<outputData.rank()-1; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(i) != inputDataRight.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match to inputDataRight dimension" );
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D): dimensions  P, D must match
        if (inputDataLeft.extent(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) does mot match to inputDataright dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check D
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) != inputDataRight.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (2) does mot match to inputDataright dimension (1)" );
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
          const size_type f1[] = { 2, 3 }, f2[] = { 1, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
        {
          const size_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match to inputDataRight dimension" );
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::matvecProduct( outputData,
                                                        inputDataLeft,
                                                        inputDataRight,
                                                        false,
                                                        transpose == 't' || transpose == 'T' );
  }


    namespace FunctorArrayTools {
    /**
       \brief Functor for matmatProduct see Intrepid2::ArrayTools for more
    */
    template < typename OutputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_matmatProduct{
      OutputViewType _output;
      leftInputViewType _leftInput;
      rightInputViewType _rightInput;
      const bool _hasField, _isTranspose;
      typedef typename leftInputViewType::value_type value_type; 

      KOKKOS_INLINE_FUNCTION
      F_matmatProduct(OutputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_,
              const bool isTranspose_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_), _isTranspose(isTranspose_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cell(0), field(0), point(0);
        size_type leftRank = _leftInput.rank();
        size_type rightRank = _rightInput.rank();
        
        if (_hasField) 
          unrollIndex( cell, field, point, 
                             _output.extent(0),
                             _output.extent(1),
                             _output.extent(2),
                             iter );
        else           
          unrollIndex( cell, point, 
                             _output.extent(0), 
                             _output.extent(1), 
                             iter );
        
        auto result = ( _hasField ? Kokkos::subview(_output, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                                    Kokkos::subview(_output, cell,        point, Kokkos::ALL(), Kokkos::ALL()));
        
        const auto lpoint = (_leftInput.extent(1) == 1 ? 0 : point);
        auto left   = ( leftRank == 4 ? Kokkos::subview(_leftInput, cell, lpoint, Kokkos::ALL(), Kokkos::ALL()) :
                        leftRank == 3 ? Kokkos::subview(_leftInput, cell, lpoint, Kokkos::ALL()) :
                                        Kokkos::subview(_leftInput, cell, lpoint) );
        
        //temporary 
        const bool hasCell = (_hasField ? rightRank == 5 : rightRank == 4); 
        auto right  = ( _hasField ? ( hasCell ? Kokkos::subview(_rightInput, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                                                Kokkos::subview(_rightInput,       field, point, Kokkos::ALL(), Kokkos::ALL())):
                                    ( hasCell ? Kokkos::subview(_rightInput, cell,        point, Kokkos::ALL(), Kokkos::ALL()) :
                                                Kokkos::subview(_rightInput,              point, Kokkos::ALL(), Kokkos::ALL())));
        
        const ordinal_type iend = result.extent(0);
        const ordinal_type jend = result.extent(1);

        switch (leftRank) {
        case 4: {
          if (_isTranspose) {
            const size_type kend = right.extent(0);
            for (ordinal_type i=0; i<iend; ++i)
              for (ordinal_type j=0; j<jend; ++j) {
                result(i, j) = value_type(0);
                for (size_type k=0; k<kend; ++k)
                  result(i, j) += left(k, i) * right(k, j);
              }
          } else {
            const size_type kend = right.extent(0);
            for (ordinal_type i=0; i<iend; ++i)
              for (ordinal_type j=0; j<jend; ++j) {
                result(i, j) = value_type(0);
                for (size_type k=0; k<kend; ++k)
                  result(i, j) += left(i, k) * right(k, j);
              }
          }
          break;
        }
        case 3: { //matrix is diagonal
          for (ordinal_type i=0; i<iend; ++i)
            for (ordinal_type j=0; j<jend; ++j)
              result(i, j) = left(i) * right(i, j);
          break;
        }
        case 2: { //matrix is a scaled identity
          for (ordinal_type i=0; i<iend; ++i) 
            for (ordinal_type j=0; j<jend; ++j) 
              result(i, j) = left() * right(i, j);
          break;
        }
        }
      }
    };
    } //namespace

  template<typename DeviceType>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<DeviceType>::Internal::
  matmatProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      OutputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_matmatProduct<OutputViewType, leftInputViewType, rightInputViewType> FunctorType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isTranspose) );
  }




  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
  matmatProductDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                          const char transpose ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
       *      (2) inputFields(C,F,P,D,D) or (F,P,D,D);
       *      (3) outputFields(C,F,P,D,D)
       */
      // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() < 2 ||
                                inputData.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputData must have rank 2,3 or 4" );
      if (inputData.rank() > 2) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(2) < 1 ||
                                  inputData.extent(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (2) must be 1,2 or 3" );
      }
      if (inputData.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(3) < 1 ||
                                  inputData.extent(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (3) must be 1,2 or 3" );
      }

      // (2) inputFields is (C,F,P,D,D) or (F,P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 4 ||
                                inputFields.rank() > 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields must have rank 4 or 5");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(inputFields.rank()-1) < 1 ||
                                inputFields.extent(inputFields.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-1) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(inputFields.rank()-2) < 1 ||
                                inputFields.extent(inputFields.rank()-2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-2) must be 1,2 or 3" );

      // (3) outputFields is (C,F,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields must have rank 5" );
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(3) < 1 ||
                                outputFields.extent(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (3) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(4) < 1 ||
                                outputFields.extent(4) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (4) must be 1,2 or 3" );

      /*
       *   Dimension cross-checks:
       *      (1) inputData    vs. inputFields
       *      (2) outputFields vs. inputData
       *      (3) outputFields vs. inputFields
       *
       *   Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P), (C,P,D) or (C,P,D,D):
       *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
       *   data is specified (P>1). Do not check P dimensions with constant data, i.e., when P=1 in
       *   inputData(C,1,...)
       */
      if (inputData.extent(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(2) != inputData.extent(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (2) does not match to inputData dimension (1)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (0) does not match to inputData dimension (0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
        const size_type f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 2 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }
      if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
        const size_type f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 3 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputData.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 5) {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D,D):  dimensions  C, P, D must match
        if (inputData.extent(1) > 1) { // check P dimension if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (2)" );
        }
        if (inputData.rank() == 2) { // inputData(C,P) -> C match
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(0) != inputFields.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
        }
        if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match

          const size_type f1[] = { 0, 2, 2 }, f2[] = { 0, 3, 4 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
          const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 4 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D,D): all dimensions C, F, P, D must match
        for (size_type i=0; i<outputFields.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(i) != inputFields.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      } else {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D,D): dimensions  P, D must match
        if (inputData.extent(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (1)" );
        }
        if (inputData.rank() == 3) {
          const size_type f1[] = { 2, 2 }, f2[] = { 2, 3 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {
          const size_type f1[] = { 2, 3 }, f2[] = { 2, 3 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D,D): dimensions F, P, D must match
        {
          const size_type f1[] = { 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3 };
          for (size_type i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(f1[i]) != inputFields.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputFields dimension" );
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::matmatProduct( outputFields,
                                                        inputData,
                                                        inputFields,
                                                        true,
                                                        transpose == 't' || transpose == 'T' );
  }



  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::matmatProductDataData(     Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                         const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                         const char transpose ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data
       *      (2) inputDataRight(C,P,D,D) or (P,D,D);
       *      (3) outputData(C,P,D,D)
       */
      // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() < 2 ||
                                inputDataLeft.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft must have rank 2,3 or 4" );
      if (inputDataLeft.rank() > 2) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) < 1 ||
                                  inputDataLeft.extent(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (2) must be 1,2 or 3" );
      }
      if (inputDataLeft.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(3) < 1 ||
                                  inputDataLeft.extent(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (3) must be 1,2 or 3" );
      }

      // (2) inputDataRight is (C,P,D,D) or (P,D,D) and 1 <= (D=dimension(rank-1),(rank-2)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 3 ||
                                inputDataRight.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(inputDataRight.rank()-1) < 1 ||
                                inputDataRight.extent(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight (rank-1) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(inputDataRight.rank()-2) < 1 ||
                                inputDataRight.extent(inputDataRight.rank()-2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight (rank-2) must be 1,2 or 3" );

      // (3) outputData is (C,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): outputData must have rank 4" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(2) < 1 ||
                                outputData.extent(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): outputData (2) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(3) < 1 ||
                                outputData.extent(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): outputData (3) must be 1,2 or 3" );

      /*
       *   Dimension cross-checks:
       *      (1) inputDataLeft vs. inputDataRight
       *      (2) outputData    vs. inputDataLeft
       *      (3) outputData    vs. inputDataRight
       *
       *   Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P), (C,P,D) or (C,P,D,D):
       *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
       *   data is specified (P>1). Do not check P dimensions with constant left data, i.e., when P=1 in
       *   inputDataLeft(C,1,...)
       */
      if (inputDataLeft.extent(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(1) != inputDataLeft.extent(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension (1) does not match to inputDataLeft dimension (1)" );
      }
      if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(0) != inputDataLeft.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension (0) does not match to inputDataLeft dimension (0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataLeft.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataLeft dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
        for (size_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataLeft.extent(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataLeft dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 4) {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D,D):  dimensions  C, P, D must match
        if (inputDataLeft.extent(1) > 1) { // check P dimension if P>1 in inputDataLeft
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (1) does not match to inputDataRight dimension (1)" );
        }
        if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (0) does not match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
          const size_type f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
          const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D,D): all dimensions C, P, D must match
        for (size_type i=0; i< outputData.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(i) !=  inputDataRight.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataRight dimension" );
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D,D): dimensions  P, D must match
        if (inputDataLeft.extent(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (1) does not match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check D
          const size_type f1[] = { 2, 2 }, f2[] = { 1, 2 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
          const size_type f1[] = { 2, 3 }, f2[] = { 1, 2 };
          for (size_type i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D,D): dimensions P, D must match
        {
          const size_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataRight dimension" );
          }
        }
      }
    }
#endif
    // body
    ArrayTools<DeviceType>::Internal::matmatProduct( outputData,
                                                        inputDataLeft,
                                                        inputDataRight,
                                                        false,
                                                        transpose == 't' || transpose == 'T' );
  }

} // end namespace Intrepid2
#endif
