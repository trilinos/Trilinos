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
    \brief  Definition file for tensor multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__

namespace Intrepid2 {

    namespace FunctorArrayTools {
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
        size_type cell, field, point;
        size_type rightRank = _rightInput.rank();

        if (_hasField) 
          Util::unrollIndex( cell, field, point, 
                       _output.dimension(0),
                       _output.dimension(1),
                       iter );
        else           
          Util::unrollIndex( cell, point, 
                       _output.dimension(0), 
                       iter );

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cell, field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_output, cell,        point, Kokkos::ALL()));

        auto left   = Kokkos::subdynrankview(_leftInput, cell, point, Kokkos::ALL());

        auto right  = (rightRank == 2 + size_type(_hasField)) ?
                      ( _hasField ? Kokkos::subdynrankview(_rightInput,       field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput,              point, Kokkos::ALL())) :
                      ( _hasField ? Kokkos::subdynrankview(_rightInput, cell, field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput, cell,        point, Kokkos::ALL()));

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
  crossProduct( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_crossProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename rightInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                            output.dimension(0)*output.dimension(1) );
    const bool isCrossProd3D = (leftInput.dimension(2) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isCrossProd3D) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
  crossProductDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
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
      INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(2) < 2 || inputData.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension(2) must be 2 or 3");

      // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputFields must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(inputFields.rank()-1) < 2 ||
                                inputFields.dimension(inputFields.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

      // (3) outputFields is (C,F,P,D) in 3D and (C,F,P) in 2D => rank = inputData.dimension(2) + 1
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != inputData.dimension(2)+1, std::invalid_argument,
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
        const size_t f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
        }
      } else {
        // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
        const size_t f1[] = { 1, 2 }, f2[] = { 1, 2 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
        }
      }
      /*
       *  Cross-check (2):
       */
      if (inputData.dimension(2) == 2) {
        //  in 2D: outputFields(C,F,P) vs. inputData(C,P,D): dimensions C,P must match
        // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
        const size_t f1[] = { 0, 2 }, f2[] = { 0, 1 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputData");
        }
      } else {
        // in 3D: outputFields(C,F,P,D) vs. inputData(C,P,D): dimensions C,P,D must match
        const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 1, 2 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputData");
        }
      }
      /*
       *  Cross-check (3):
       */
      if (inputData.dimension(2) == 2) {
        // In 2D:
        if (inputFields.rank() == 4) {
          //  and rank-4 inputFields: outputFields(C,F,P) vs. inputFields(C,F,P,D): dimensions C,F,P must match
          const size_t f1[] = { 0, 1, 2 }, f2[] = { 0, 1, 2 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        } else {
          //  and rank-3 inputFields: outputFields(C,F,P) vs. inputFields(F,P,D): dimensions F,P must match
          const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      } else {
        // In 3D:
        if (inputFields.rank() == 4) {
          //  and rank-4 inputFields: outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C,F,P,D must match
          for (size_t i=0; i<outputFields.rank(); ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(i) != inputFields.dimension(i), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
          }
        } else {
          // and rank-3 inputFields: outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F,P,D must match
          const size_t f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
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
  crossProductDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
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
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) < 2 || inputDataLeft.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension(2) must be 2 or 3");

      // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(inputDataRight.rank()-1) < 2 ||
                                inputDataRight.dimension(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight dimension (rank-1) must have rank 2 or 3" );

      // (3) outputData is (C,P,D) in 3D and (C,P) in 2D => rank = inputDataLeft.dimension(2)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != inputDataLeft.dimension(2), std::invalid_argument,
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
        for (size_t i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(i) != inputDataRight.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
        }
      }
      // inputDataLeft(C, P,D) vs. inputDataRight(P,D): dimensions P, D must match
      else {
        const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
        }
      }
      /*
       *  Cross-check (2):
       */
      if (inputDataLeft.dimension(2) == 2) {
        // in 2D: outputData(C,P) vs. inputDataLeft(C,P,D): dimensions C, P must match
        const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != outputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to outputData");
        }
      } else {
        // in 3D: outputData(C,P,D) vs. inputDataLeft(C,P,D): all dimensions C, P, D must match
        for (size_t i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(i) != outputData.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to outputData");
        }
      }
      /*
       *  Cross-check (3):
       */
      if (inputDataLeft.dimension(2) == 2) {
        // In 2D:
        if (inputDataRight.rank() == 3) {
          //  and rank-3 inputDataRight: outputData(C,P) vs. inputDataRight(C,P,D): dimensions C,P must match
          const size_t f1[] = { 0, 1 }, f2[] = { 0, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
          }
        } else {
          //  and rank-2 inputDataRight: outputData(C,P) vs. inputDataRight(P,D): dimension P must match
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(1) != inputDataRight.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
        }
      } else {
        // In 3D:
        if (inputDataRight.rank() == 3) {
          //  and rank-3 inputDataRight: outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C,P,D must match
          for (size_t i=0; i<outputData.rank(); ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(i) != inputDataRight.dimension(i), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
          }
        } else {
          //  and rank-2 inputDataRight: outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
          const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
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
        size_type cell, field, point;
        size_type rightRank = _rightInput.rank();

        if (_hasField) 
          Util::unrollIndex( cell, field, point, 
                       _output.dimension(0),
                       _output.dimension(1),
                       iter );
        else           
          Util::unrollIndex( cell, point, 
                       _output.dimension(0), 
                       iter );

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_output, cell,        point, Kokkos::ALL(), Kokkos::ALL()));

        auto left   = Kokkos::subdynrankview(_leftInput, cell, point, Kokkos::ALL());

        auto right  = (rightRank == 2 + size_type(_hasField)) ?
                      ( _hasField ? Kokkos::subdynrankview(_rightInput,       field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput,              point, Kokkos::ALL())) :
                      ( _hasField ? Kokkos::subdynrankview(_rightInput, cell, field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput, cell,        point, Kokkos::ALL()));

        const size_type iend = result.dimension(0);
        const size_type jend = result.dimension(1);
        for (size_type i=0; i<iend; ++i)
          for (size_type j=0; j<jend; ++j)
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
  outerProduct( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                const bool hasField ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_outerProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                 /**/       output.dimension(0)*output.dimension(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField) );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
  outerProductDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
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
      INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(2) < 2 || inputData.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension(2) must be 2 or 3");

      // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputFields must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(inputFields.rank()-1) < 2 ||
                                inputFields.dimension(inputFields.rank()-1) < 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

      // (3) outputFields is (C,F,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputFields must have rank 5");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(3) < 2 ||
                                outputFields.dimension(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension(3) must be 2 or 3");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(4) < 2 ||
                                outputFields.dimension(4) > 3, std::invalid_argument,
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
        const size_t f1[] = { 0, 2, 3, 4 }, f2[] = { 0, 1, 2, 2 };
        for (size_t i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputData");
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 4) {
        // Cross-check (1): inputData(C,P,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        {
          const size_t f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D): dimensions C, F, P, D must match
        {
          const size_t f1[] = { 0, 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3, 3 };
          for (size_t i=0; i<5; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputFields");
          }
        }
      } else {
        // Cross-check (1): inputData(C,P,D) vs. inputFields(F,P,D): dimensions  P, D must match
        {
          const size_t f1[] = { 1, 2 }, f2[] = { 1, 2 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        {
          const size_t f1[] = { 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 2 };
          for (size_t i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
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
  outerProductDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
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
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) < 2 || inputDataLeft.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension(2) must be 2 or 3");

      // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(inputDataRight.rank()-1) < 2 ||
                                inputDataRight.dimension(inputDataRight.rank()-1) < 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight dimension (rank-1) must have rank 2 or 3" );

      // (3) outputData is (C,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputData must have rank 5");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(2) < 2 ||
                                outputData.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension(3) must be 2 or 3");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(3) < 2 ||
                                outputData.dimension(3) > 3, std::invalid_argument,
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
        const size_t f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
        for (size_t i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataLeft");
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 3) {
        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(C,P,D):  all dimensions  C, P, D must match
        for (size_t i=0; i<inputDataLeft.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(i) != inputDataRight.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D): dimensions C, P, D must match
        {
          const size_t f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
          for (size_t i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataRight");
          }
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(P,D): dimensions  P, D must match
        {
          const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
          }
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D): dimensions P, D must match
        {
          const size_t f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 1 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
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
    template < typename outputViewType, typename leftInputViewType, typename rightInputViewType >
    struct F_matvecProduct {
      outputViewType _output;
      leftInputViewType _leftInput;
      rightInputViewType _rightInput;
      const bool _hasField, _isTranspose;
      typedef typename leftInputViewType::value_type value_type; 

      KOKKOS_INLINE_FUNCTION
      F_matvecProduct(outputViewType output_,
              leftInputViewType leftInput_,
              rightInputViewType rightInput_,
              const bool hasField_,
              const bool isTranspose_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_),
          _hasField(hasField_), _isTranspose(isTranspose_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cell, field, point;
        size_type rightRank = _rightInput.rank();
        size_type leftRank = _leftInput.rank();

        if (_hasField) 
          Util::unrollIndex( cell, field, point, 
                       _output.dimension(0),
                       _output.dimension(1),
                       iter );
        else           
          Util::unrollIndex( cell, point, 
                       _output.dimension(0), 
                       iter );

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cell, field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_output, cell,        point, Kokkos::ALL()));

        auto lpoint = (_leftInput.dimension(1) == 1) ? size_type(0) : point;
        auto left   = ( leftRank == 4 ? Kokkos::subdynrankview(_leftInput, cell, lpoint, Kokkos::ALL(), Kokkos::ALL()) :
                      ((leftRank == 3) ? Kokkos::subdynrankview(_leftInput, cell, lpoint, Kokkos::ALL()) :
                       Kokkos::subdynrankview(_leftInput, cell, lpoint)));

        auto right  = (rightRank == 2 + size_type(_hasField)) ?
                      ( _hasField ? Kokkos::subdynrankview(_rightInput,       field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput,              point, Kokkos::ALL())) :
                      ( _hasField ? Kokkos::subdynrankview(_rightInput, cell, field, point, Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput, cell,        point, Kokkos::ALL()));

        const size_type iend = result.dimension(0);
        const size_type jend = right.dimension(0);
        switch (leftRank) {
          case 4:
            if (!_isTranspose)
              for (size_type i=0; i<iend; ++i) {
                result(i) = value_type(0);
                for (size_type j=0; j<jend; ++j)
                  result(i) += left(i, j) * right(j);
            } else
              for (size_type i=0; i<iend; ++i) {
                result(i) = value_type(0);
                for (size_type j=0; j<jend; ++j)
                  result(i) += left(j, i) * right(j);
            }
            break;
          case 3:  //matrix is diagonal
            for (size_type i=0; i<iend; ++i)
                result(i) = left(i) * right(i);
            break;
          case 2:  //matrix is a scaled identity
            for (size_type i=0; i<iend; ++i) {
                result(i) = left() * right(i);
            }
            break;
        }
      }
    };
    } //namespace

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::matvecProduct( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_matvecProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                 /**/       output.dimension(0)*output.dimension(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isTranspose) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::matvecProductDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
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
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(2) < 1 ||
                                  inputData.dimension(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) must be 1,2 or 3");
      }
      if (inputData.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(3) < 1 ||
                                  inputData.dimension(3) > 3, std::invalid_argument,
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
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(3) < 1 ||
                                outputFields.dimension(3) > 3, std::invalid_argument,
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
      if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(2) != inputData.dimension(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(2) must match to inputData dimension(1)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(0) must match to inputData dimension(0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
        const size_t f1[] = { 0, 3 }, f2[] = { 0, 2 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }
      if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
        const size_t f1[] = { 0, 3, 3 }, f2[] = { 0, 2, 3 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 4) {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(2) != inputData.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputFields dimension (2) does not match to inputData dimension (1)" );
        }
        if (inputData.rank() == 2) { // inputData(C,P) -> C match
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(0) != inputFields.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
        }
        if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
          const size_t f1[] = { 0, 2 }, f2[] = { 0, 3 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
          const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 3 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C, F, P, D must match
        for (size_t i=0; i<outputFields.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(i) != inputFields.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      } else {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D): dimensions  P, D must match
        if (inputData.dimension(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(1) does not match to inputFields dimension (1)" );
        }
        if (inputData.rank() == 3) {
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(2) != inputFields.dimension(2), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) does not match to inputFields dimension (2)" );
        }
        if (inputData.rank() == 4) {
          const size_t f1[] = { 2, 3 }, f2[] = { 2, 2 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        {
          const size_t f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
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
  matvecProductDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>    outputData,
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
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) < 1 ||
                                  inputDataLeft.dimension(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(2) must be 1, 2 or 3");
      }
      if (inputDataLeft.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(3) < 1 ||
                                  inputDataLeft.dimension(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(3) must be 1, 2 or 3");
      }

      // (2) inputDataRight is (C, P, D) or (P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 ||
                                inputDataRight.rank() > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight must have rank 2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(inputDataRight.rank()-1) < 1 ||
                                inputDataRight.dimension(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight dimension (rank-1) must be 1,2 or 3" );

      // (3) outputData is (C,P,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): outputData must have rank 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(2) < 1 ||
                                outputData.dimension(2) > 3, std::invalid_argument,
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
      if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(1) != inputDataLeft.dimension(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(1) muat match inputDataLeft dimension(1)" );
      }
      if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(0) != inputDataLeft.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(0) muat match inputDataLeft dimension(0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const size_t f1[] = { 0, 2 }, f2[] = { 0, 2 };
        for (size_t i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const size_t f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 3) {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D):  dimensions  C, P, D must match
        if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) muat match to inputDataRight dimension (1)" );
        }
        if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (0) muat match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
          const size_t f1[] = { 0, 2 }, f2[] = { 0, 2 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
          const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
        for (size_t i=0; i<outputData.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(i) != inputDataRight.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match to inputDataRight dimension" );
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D): dimensions  P, D must match
        if (inputDataLeft.dimension(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) does mot match to inputDataright dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check D
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) != inputDataRight.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (2) does mot match to inputDataright dimension (1)" );
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
          const size_t f1[] = { 2, 3 }, f2[] = { 1, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
        {
          const size_t f1[] = { 1, 2 }, f2[] = { 0, 1 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
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
        size_type cell, field, point;
        size_type leftRank = _leftInput.rank();

        if (_hasField) 
          Util::unrollIndex( cell, field, point, 
                       _output.dimension(0),
                       _output.dimension(1),
                       iter );
        else           
          Util::unrollIndex( cell, point, 
                       _output.dimension(0), 
                       iter );

        auto result = ( _hasField ? Kokkos::subdynrankview(_output, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_output, cell,        point, Kokkos::ALL(), Kokkos::ALL()));

        auto lpoint = (_leftInput.dimension(1) == 1) ? size_type(0) : point;
        auto left   = ( leftRank == 4 ? Kokkos::subdynrankview(_leftInput, cell, lpoint, Kokkos::ALL(), Kokkos::ALL()) :
                      ((leftRank == 3) ? Kokkos::subdynrankview(_leftInput, cell, lpoint, Kokkos::ALL()) :
                      Kokkos::subdynrankview(_leftInput, cell, lpoint)));

        auto right  = ( _hasField ? Kokkos::subdynrankview(_rightInput, cell, field, point, Kokkos::ALL(), Kokkos::ALL()) :
                        /**/        Kokkos::subdynrankview(_rightInput, cell,        point, Kokkos::ALL(), Kokkos::ALL()));

        const size_type iend = result.dimension(0);
        const size_type jend = result.dimension(1);

        switch (leftRank) {
          case 4:
            if (_isTranspose) {
              const size_type kend = right.dimension(0);
              for (size_type i=0; i<iend; ++i)
                for (size_type j=0; j<jend; ++j) {
                  result(i, j) = value_type(0);
                  for (size_type k=0; k<kend; ++k)
                    result(i, j) += left(k, i) * right(k, j);
                }
            } else {
              const size_type kend = right.dimension(0);
              for (size_type i=0; i<iend; ++i)
                for (size_type j=0; j<jend; ++j) {
                  result(i, j) = value_type(0);
                  for (size_type k=0; k<kend; ++k)
                    result(i, j) += left(i, k) * right(k, j);
                }
            }
            break;
          case 3:  //matrix is diagonal
            for (size_type i=0; i<iend; ++i)
              for (size_type j=0; j<jend; ++j)
                result(i, j) = left(i) * right(i, j);
            break;
          case 2:  //matrix is a scaled identity
            for (size_type i=0; i<iend; ++i) {
              for (size_type j=0; j<jend; ++j)
                result(i, j) = left() * right(i, j);
            }
            break;
        }


      }
    };
    } //namespace

  template<typename SpT>
  template<typename outputValueType,     class ...outputProperties,
           typename leftInputValueType,  class ...leftInputProperties,
           typename rightInputValueType, class ...rightInputProperties>
  void
  ArrayTools<SpT>::Internal::matmatProduct( /**/  Kokkos::DynRankView<outputValueType,    outputProperties...>      output,
                 const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInput,
                 const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInput,
                 const bool hasField,
                 const bool isTranspose ) {

    typedef Kokkos::DynRankView<outputValueType,    outputProperties...>      outputViewType;
    typedef const Kokkos::DynRankView<leftInputValueType, leftInputProperties...>   leftInputViewType;
    typedef const Kokkos::DynRankView<rightInputValueType,rightInputProperties...>  rightInputViewType;
    typedef FunctorArrayTools::F_matmatProduct<outputViewType, leftInputViewType, rightInputViewType> FunctorType;
    typedef typename ExecSpace< typename leftInputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;

    const size_type loopSize = ( hasField ? output.dimension(0)*output.dimension(1)*output.dimension(2) :
                                 /**/       output.dimension(0)*output.dimension(1) );
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, leftInput, rightInput, hasField, isTranspose) );
  }




  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
  matmatProductDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
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
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(2) < 1 ||
                                  inputData.dimension(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (2) must be 1,2 or 3" );
      }
      if (inputData.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(3) < 1 ||
                                  inputData.dimension(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (3) must be 1,2 or 3" );
      }

      // (2) inputFields is (C,F,P,D,D) or (F,P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 4 ||
                                inputFields.rank() > 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields must have rank 4 or 5");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(inputFields.rank()-1) < 1 ||
                                inputFields.dimension(inputFields.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-1) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(inputFields.rank()-2) < 1 ||
                                inputFields.dimension(inputFields.rank()-2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-2) must be 1,2 or 3" );

      // (3) outputFields is (C,F,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 5, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields must have rank 5" );
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(3) < 1 ||
                                outputFields.dimension(3) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (3) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(4) < 1 ||
                                outputFields.dimension(4) > 3, std::invalid_argument,
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
      if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(2) != inputData.dimension(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (2) does not match to inputData dimension (1)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (0) does not match to inputData dimension (0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
        const size_t f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 2 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }
      if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
        const size_t f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 3 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputFields.rank() == 5) {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D,D):  dimensions  C, P, D must match
        if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(2), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (2)" );
        }
        if (inputData.rank() == 2) { // inputData(C,P) -> C match
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(0) != inputFields.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
        }
        if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match

          const size_t f1[] = { 0, 2, 2 }, f2[] = { 0, 3, 4 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
          const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 4 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D,D): all dimensions C, F, P, D must match
        for (size_t i=0; i<outputFields.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(i) != inputFields.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      } else {
        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D,D): dimensions  P, D must match
        if (inputData.dimension(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (1)" );
        }
        if (inputData.rank() == 3) {
          const size_t f1[] = { 2, 2 }, f2[] = { 2, 3 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }
        if (inputData.rank() == 4) {
          const size_t f1[] = { 2, 3 }, f2[] = { 2, 3 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
          }
        }

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D,D): dimensions F, P, D must match
        {
          const size_t f1[] = { 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3 };
          for (size_t i=0; i<4; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]), std::invalid_argument,
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
  ArrayTools<ExecSpaceType>::matmatProductDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
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
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) < 1 ||
                                  inputDataLeft.dimension(2) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (2) must be 1,2 or 3" );
      }
      if (inputDataLeft.rank() == 4) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(3) < 1 ||
                                  inputDataLeft.dimension(3) > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (3) must be 1,2 or 3" );
      }

      // (2) inputDataRight is (C,P,D,D) or (P,D,D) and 1 <= (D=dimension(rank-1),(rank-2)) <= 3 is required.
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 3 ||
                                inputDataRight.rank() > 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight must have rank 3 or 4" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(inputDataRight.rank()-1) < 1 ||
                                inputDataRight.dimension(inputDataRight.rank()-1) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight (rank-1) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.dimension(inputDataRight.rank()-2) < 1 ||
                                inputDataRight.dimension(inputDataRight.rank()-2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight (rank-2) must be 1,2 or 3" );

      // (3) outputData is (C,P,D,D)
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 4, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): outputData must have rank 4" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(2) < 1 ||
                                outputData.dimension(2) > 3, std::invalid_argument,
                                ">>> ERROR (ArrayTools::matmatProductDataData): outputData (2) must be 1,2 or 3" );
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(3) < 1 ||
                                outputData.dimension(3) > 3, std::invalid_argument,
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
      if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(1) != inputDataLeft.dimension(1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension (1) does not match to inputDataLeft dimension (1)" );
      }
      if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(0) != inputDataLeft.dimension(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension (0) does not match to inputDataLeft dimension (0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataLeft dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
        for (size_t i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataLeft dimension" );
        }
      }

      /*
       *   Cross-checks (1,3):
       */
      if (inputDataRight.rank() == 4) {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D,D):  dimensions  C, P, D must match
        if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (1) does not match to inputDataRight dimension (1)" );
        }
        if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (0) does not match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
          const size_t f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
          const size_t f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D,D): all dimensions C, P, D must match
        for (size_t i=0; i< outputData.rank(); ++i) {
          INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(i) !=  inputDataRight.dimension(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): outputData dimension does not match to inputDataRight dimension" );
        }
      } else {
        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D,D): dimensions  P, D must match
        if (inputDataLeft.dimension(1) > 1) { // check P if P>1 in inputData
          INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(0), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft dimension (1) does not match to inputDataRight dimension (0)" );
        }
        if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check D
          const size_t f1[] = { 2, 2 }, f2[] = { 1, 2 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
          const size_t f1[] = { 2, 3 }, f2[] = { 1, 2 };
          for (size_t i=0; i<2; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::matmatProductDataData): intputDataLeft dimension does not match to inputDataRight dimension" );
          }
        }
        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D,D): dimensions P, D must match
        {
          const size_t f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
          for (size_t i=0; i<3; ++i) {
            INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]), std::invalid_argument,
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
