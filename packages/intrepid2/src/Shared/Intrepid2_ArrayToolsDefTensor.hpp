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
    template < typename outputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_crossProduct{
      outputViewType _output;
      const leftInputViewType _leftInput;
      const rightInputViewType _rightInput;
      const bool _hasField, _isCrossProd3D;

      KOKKOS_INLINE_FUNCTION
      F_crossProduct(outputViewType output_,
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

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::
  crossProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_crossProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename rightInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    const bool isCrossProd3D = (leftInput.extent(2) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isCrossProd3D) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
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
    ArrayTools<ExecSpaceType>::Internal::crossProduct( outputFields,
                                                   inputData,
                                                   inputFields,
                                                   true );
  }


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<ExecSpaceType>::
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
    ArrayTools<ExecSpaceType>::Internal::crossProduct( outputData,
                                                   inputDataLeft,
                                                   inputDataRight,
                                                   false );
  }


    namespace FunctorArrayTools {
    /**
       \brief Functor for outerProduct see Intrepid2::ArrayTools for more
    */
    template < typename outputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_outerProduct {
      outputViewType _output;
      const leftInputViewType _leftInput;
      const rightInputViewType _rightInput;
      const bool _hasField;

      KOKKOS_INLINE_FUNCTION
      F_outerProduct(outputViewType output_,
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

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::
  outerProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_outerProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField) );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
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
    ArrayTools<ExecSpaceType>::Internal::outerProduct(  outputFields,
                                                                 inputData,
                                                                 inputFields,
                                                                 true );
  }


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValuetype,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<ExecSpaceType>::
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
    ArrayTools<ExecSpaceType>::Internal::outerProduct( outputData,
                                                               inputDataLeft,
                                                               inputDataRight,
                                                               false );
  }


  namespace FunctorArrayTools {
    /**
       \brief Functor for matvecProduct see Intrepid2::ArrayTools for more
    */
    template < typename outputViewType, 
               typename leftInputViewType, 
               typename rightInputViewType>
    struct F_matvecProduct {
      /**/  outputViewType     _output;
      const leftInputViewType  _leftInput;
      const rightInputViewType _rightInput;

      const bool _isTranspose;
      
      KOKKOS_INLINE_FUNCTION
      F_matvecProduct(outputViewType     output_,
                      leftInputViewType  leftInput_,
                      rightInputViewType rightInput_,
                      const bool isTranspose_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_), _isTranspose(isTranspose_) {}

      template<typename resultViewType,
               typename leftViewType,
               typename rightViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void 
      apply_matvec_product(      resultViewType &result, 
                           const leftViewType   &left,
                           const rightViewType  &right,
                           const bool isTranspose) {
        const ordinal_type iend = result.extent(0);
        const ordinal_type jend = right.extent(0);

        typedef typename resultViewType::value_type value_type; 

        switch (left.rank()) {
        case 2:
          if (isTranspose) {
            for (ordinal_type i=0;i<iend;++i) {
              value_type tmp(0);
              for (ordinal_type j=0;j<jend;++j)
                tmp += left(j, i)*right(j);
              result(i) = tmp;
            }
          } else {
            for (ordinal_type i=0;i<iend;++i) {
              value_type tmp(0);
              for (ordinal_type j=0;j<jend;++j)
                tmp += left(i, j)*right(j);
              result(i) = tmp;
            }
          }
          break;
        case 1: { //matrix is diagonal
          for (ordinal_type i=0;i<iend;++i)
            result(i) = left(i)*right(i);
          break;
        }
        case 0:  { //matrix is a scaled identity
          const value_type val = left();
          for (ordinal_type i=0;i<iend;++i) {
            result(i) = val*right(i);
          }
          break;
        }
        }
      }
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type pt) const {
        const auto rightRank = _rightInput.rank();
        const auto leftRank  = _leftInput.rank();
        
        auto result = Kokkos::subview(_output, cl, pt, Kokkos::ALL());
        
        const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
        const auto left = ( leftRank == 4 ? Kokkos::subview(_leftInput, cl, lpt, Kokkos::ALL(), Kokkos::ALL()) :
                            leftRank == 3 ? Kokkos::subview(_leftInput, cl, lpt, Kokkos::ALL()) :
                                            Kokkos::subview(_leftInput, cl, lpt));
        
        const auto right = ( rightRank == 2 ? Kokkos::subview(_rightInput,     pt, Kokkos::ALL()) :
                                              Kokkos::subview(_rightInput, cl, pt, Kokkos::ALL()) );
        apply_matvec_product( result, left, right, _isTranspose );
      }
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type bf,
                      const ordinal_type pt) const {
        const auto rightRank = _rightInput.rank();
        const auto leftRank  = _leftInput.rank();

        auto result = Kokkos::subview(_output, cl, bf, pt, Kokkos::ALL());

        const auto lpt  = (_leftInput.extent(1) == 1 ? size_type(0) : pt);
        const auto left = ( leftRank == 4 ? Kokkos::subview(_leftInput, cl, lpt, Kokkos::ALL(), Kokkos::ALL()) :
                            leftRank == 3 ? Kokkos::subview(_leftInput, cl, lpt, Kokkos::ALL()) :
                                            Kokkos::subview(_leftInput, cl, lpt));
        
        const auto right = ( rightRank == 3 ? Kokkos::subview(_rightInput,     bf, pt, Kokkos::ALL()) :
                                              Kokkos::subview(_rightInput, cl, bf, pt, Kokkos::ALL())); 
        
        apply_matvec_product( result, left, right, _isTranspose );
      }
    };
  } //namespace

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::
  matvecProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {

    typedef       Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_matvecProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    if (hasField) {
      using range_policy_type = Kokkos::Experimental::MDRangePolicy
        < ExecSpaceType, Kokkos::Experimental::Rank<3>, Kokkos::IndexType<ordinal_type> >;
      range_policy_type policy( { 0, 0, 0 },
                                { output.extent(0), output.extent(1), output.extent(2) } );
      Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, isTranspose) );
    } else {
      using range_policy_type = Kokkos::Experimental::MDRangePolicy
        < ExecSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
      range_policy_type policy( { 0, 0 },
                                { output.extent(0), output.extent(1) } );
      Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, isTranspose) );
    }
  }

  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::matvecProductDataField(      Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
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
    ArrayTools<ExecSpaceType>::Internal::matvecProduct( outputFields,
                                                        inputData,
                                                        inputFields,
                                                        true,
                                                        transpose == 't' || transpose == 'T' );
  }



  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<ExecSpaceType>::
  matvecProductDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>    outputData,
                         const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...> inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...>    inputDataRight,
                         const char transpose ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data
       *      (2) inputDataRight(C,P,D) or (P,D);
       *      (3) outputData(C,P,D)
       */
      // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required
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
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const size_type f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
        for (size_type i=0; i<3; ++i) {
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
          const size_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
          for (size_type i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(f1[i]) != inputDataRight.extent(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
        for (size_type i=0; i<outputData.rank(); ++i) {
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
    ArrayTools<ExecSpaceType>::Internal::matvecProduct( outputData,
                                                        inputDataLeft,
                                                        inputDataRight,
                                                        false,
                                                        transpose == 't' || transpose == 'T' );
  }


    namespace FunctorArrayTools {
    /**
       \brief Functor for matmatProduct see Intrepid2::ArrayTools for more
    */
    template < typename outputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_matmatProduct{
      outputViewType _output;
      leftInputViewType _leftInput;
      rightInputViewType _rightInput;
      const bool _hasField, _isTranspose;
      typedef typename leftInputViewType::value_type value_type; 

      KOKKOS_INLINE_FUNCTION
      F_matmatProduct(outputViewType output_,
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

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::
  matmatProduct(       Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_matmatProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.extent(0)*output.extent(1)*output.extent(2) :
                                            output.extent(0)*output.extent(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isTranspose) );
  }




  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
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
    ArrayTools<ExecSpaceType>::Internal::matmatProduct( outputFields,
                                                        inputData,
                                                        inputFields,
                                                        true,
                                                        transpose == 't' || transpose == 'T' );
  }



  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<ExecSpaceType>::matmatProductDataData(     Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
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
    ArrayTools<ExecSpaceType>::Internal::matmatProduct( outputData,
                                                        inputDataLeft,
                                                        inputDataRight,
                                                        false,
                                                        transpose == 't' || transpose == 'T' );
  }

} // end namespace Intrepid2
#endif
