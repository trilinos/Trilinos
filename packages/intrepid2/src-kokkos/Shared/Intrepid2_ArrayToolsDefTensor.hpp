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

  template<typename ExecSpaceType>
  template<class ...outputProperties,
           class ...leftInputProperties,
           class ...rightInputProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  crossProduct( /**/  Kokkos::DynRankView<outputProperties...>      output,
                const Kokkos::DynRankView<leftInputProperties...>   leftInput,
                const Kokkos::DynRankView<rightInputProperties...>  rightInput ) {
    
    struct Functor {
      /**/  Kokkos::DynRankView<outputProperties...>    _output;
      const Kokkos::DynRankView<leftInputProperties...> _leftInput;
      const Kokkos::DynRankView<leftInputProperties...> _rightInput;
      const bool _is_crossprod_3d;
      
      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputProperties...>    &output_,
              Kokkos::DynRankView<leftInputProperties...> &leftInput_,
              Kokkos::DynRankView<leftInputProperties...> &rightInput_,
              const bool is_crossprod_3d_);
      : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_), _is_crossprod_3d(is_crossprod_3d_) {}
      
      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        size_type i, j, k;
        Util::unrollIndex(i, j, k, _output, iter);
        
        auto result = Kokkos::subdynrankview(_output,     i, j, k, Kokkos::ALL());
        auto left   = Kokkos::subdynrankview(_leftInput,  i,    k, Kokkos::ALL());
    
        const size_type r = _rightInput.rank();
        auto right  = ( r == 3 ? Kokkos::subdynrankview(_rightInput,     j, k, Kokkos::ALL()) :
                        /**/     Kokkos::subdynrankview(_rightInput,  i, j, k, Kokkos::ALL()) );
        
        // branch prediction is possible
        if (_is_crossprod_3d) {
          result(0) = left(1)*right(2) - left(2)*right(1);
          result(1) = left(2)*right(0) - left(0)*right(2);
          result(2) = left(0)*right(1) - left(1)*right(0);
        } else {
          result(0) = left(0)*right(1) - left(1)*right(0);
        }
      }
    };
    const size_type loopSize = output.dimension(0)*output.dimension(1)*output.dimension(2);
    const bool is_crossprod_3d = (leftInput.dimension(2) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(output, leftInput, rightInput, is_crossprod_3d) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  crossProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P,D);
     *      (2) inputFields(C,F,P,D) or (F,P,D);
     *      (3) outputFields(C,F,P,D) in 3D, or (C,F,P) in 2D
     */
    // (1) inputData is (C, P, D) and 2 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 3
                              ">>> ERROR (ArrayTools::crossProductDataField): inputData must have rank 3");
    INTREPID2_TEST_FOR_ABORT( inputData.dimension(2) < 2 || inputData.dimension(2) > 3,
                              ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension(2) must be 2 or 3");

    // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 3 || inputFields.rank() > 4,
                              ">>> ERROR (ArrayTools::crossProductDataField): inputFields must have rank 3 or 4" );
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(inputFields.rank()-1) < 2 ||
                              inputFields.dimension(inputFields.rank()-1) < 3,
                              ">>> ERROR (ArrayTools::crossProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

    // (3) outputFields is (C,F,P,D) in 3D and (C,F,P) in 2D => rank = inputData.dimension(2) + 1
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != inputData.dimension(2)+1,
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
      const ordinal_type f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
      }
    } else {
      // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
      const ordinal_type f1[] = { 1, 2 }, f2[] = { 1, 2 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::crossProductDataField): inputData dimension does not match with inputFields");
      }
    }
    /*
     *  Cross-check (2):
     */
    if (inputData.dimension(2) == 2) {
      //  in 2D: outputFields(C,F,P) vs. inputData(C,P,D): dimensions C,P must match
      // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match
      const ordinal_type f1[] = { 0, 2 }, f2[] = { 0, 1 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputData");
      }
    } else {
      // in 3D: outputFields(C,F,P,D) vs. inputData(C,P,D): dimensions C,P,D must match
      const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 1, 2 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
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
        const ordinal_type f1[] = { 0, 1, 2 }, f2[] = { 0, 1, 2 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
        }
      } else {
        //  and rank-3 inputFields: outputFields(C,F,P) vs. inputFields(F,P,D): dimensions F,P must match
        const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
        }
      }
    } else {
      // In 3D:
      if (inputFields.rank() == 4) {
        //  and rank-4 inputFields: outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C,F,P,D must match
        for (ordinal_type i=0; i<outputFields.rank(); ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(i) != inputFields.dimension(i),
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
        }
      } else {
        // and rank-3 inputFields: outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F,P,D must match
        const ordinal_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::crossProductDataField): outputFields dimension does not match with inputFields");
        }
      }
    }
#endif

    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  crossProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::crossProductDataData):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P,D);
     *      (2) inputDataRight(C,P,D) or (P,D);
     *      (3) outputData(C,P,D) in 3D, or (C,P) in 2D
     */
    // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() != 3
                              ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft must have rank 3");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) < 2 || inputDataLeft.dimension(2) > 3,
                              ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension(2) must be 2 or 3");

    // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 2 || inputDataRight.rank() > 3,
                              ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight must have rank 2 or 3" );
    INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(inputDataRight.rank()-1) < 2 ||
                              inputDataRight.dimension(inputDataRight.rank()-1) < 3,
                              ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight dimension (rank-1) must have rank 2 or 3" );

    // (3) outputData is (C,P,D) in 3D and (C,P) in 2D => rank = inputDataLeft.dimension(2)
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != inputDataLeft.dimension(2),
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
      for (ordinal_type i=0; i<inputDataLeft.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(i) != inputDataRight.dimension(i),
                                  ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
      }
    }
    // inputDataLeft(C, P,D) vs. inputDataRight(P,D): dimensions P, D must match
    else {
      const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to inputDataRight");
      }
    }
    /*
     *  Cross-check (2):
     */
    if (inputDataLeft.dimension(2) == 2) {
      // in 2D: outputData(C,P) vs. inputDataLeft(C,P,D): dimensions C, P must match
      const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != outputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::crossProductDataData): inputDataLeft dimension does not match to outputData");
      }
    } else {
      // in 3D: outputData(C,P,D) vs. inputDataLeft(C,P,D): all dimensions C, P, D must match
      for (ordinal_type i=0; i<inputDataLeft.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(i) != outputData.dimension(i),
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
        const ordinal_type f1[] = { 0, 1 }, f2[] = { 0, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
        }
      } else {
        //  and rank-2 inputDataRight: outputData(C,P) vs. inputDataRight(P,D): dimension P must match
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(1) != inputDataRight.dimension(0),
                                  ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
      }
    } else {
      // In 3D:
      if (inputDataRight.rank() == 3) {
        //  and rank-3 inputDataRight: outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C,P,D must match
        for (ordinal_type i=0; i<outputData.rank(); ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(i) != inputDataRight.dimension(i),
                                    ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
        }
      } else {
        //  and rank-2 inputDataRight: outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
        const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::crossProductDataData): outputData dimension does not match to inputDataRight");
        }
      }
    }
#endif
    // body

  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  outerProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P,D);
     *      (2) inputFields(C,F,P,D) or (F,P,D);
     *      (3) outputFields(C,F,P,D,D)
     */
    // (1) inputData is (C, P, D) and 2 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 3
                              ">>> ERROR (ArrayTools::outerProductDataField): inputData must have rank 3");
    INTREPID2_TEST_FOR_ABORT( inputData.dimension(2) < 2 || inputData.dimension(2) > 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension(2) must be 2 or 3");

    // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 3 || inputFields.rank() > 4,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputFields must have rank 3 or 4" );
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(inputFields.rank()-1) < 2 ||
                              inputFields.dimension(inputFields.rank()-1) < 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputFields dimension (rank-1) must have rank 2 or 3" );

    // (3) outputFields is (C,F,P,D,D)
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 5,
                              ">>> ERROR (ArrayTools::outerProductDataField): outputFields must have rank 5");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(3) < 2 ||
                              outputFields.dimension(3) > 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension(3) must be 2 or 3");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(4) < 2 ||
                              outputFields.dimension(4) > 3,
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
      const ordinal_type f1[] = { 0, 2, 3, 4 }, f2[] = { 0, 1, 2, 2 };
      for (ordinal_type i=0; i<4; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputData");
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (inputFields.rank() == 4) {
      // Cross-check (1): inputData(C,P,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
      {
        const ordinal_type f1[] = { 0, 1, 2 }, f2[] = { 0, 2, 3 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
        }
      }

      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D): dimensions C, F, P, D must match
      {
        const ordinal_type f1[] = { 0, 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3, 3 };
        for (ordinal_type i=0; i<5; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputFields");
        }
      }
    } else {
      // Cross-check (1): inputData(C,P,D) vs. inputFields(F,P,D): dimensions  P, D must match
      {
        const ordinal_type f1[] = { 1, 2 }, f2[] = { 1, 2 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): inputData dimension does not match with inputFields");
        }
      }

      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D): dimensions F, P, D must match
      {
        const ordinal_type f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
        for (ordinal_type i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputFields dimension does not match with inputFields");
        }
      }
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  outerProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P,D);
     *      (2) inputDataRight(C,P,D) or (P,D);
     *      (3) outputData(C,P,D,D)
     */
    // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() != 3
                              ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft must have rank 3");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) < 2 || inputDataLeft.dimension(2) > 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension(2) must be 2 or 3");

    // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 2 || inputDataRight.rank() > 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight must have rank 2 or 3" );
    INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(inputDataRight.rank()-1) < 2 ||
                              inputDataRight.dimension(inputDataRight.rank()-1) < 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): inputDataRight dimension (rank-1) must have rank 2 or 3" );

    // (3) outputData is (C,P,D,D)
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 4,
                              ">>> ERROR (ArrayTools::outerProductDataField): outputData must have rank 5");
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(2) < 2 ||
                              outputData.dimension(2) > 3,
                              ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension(3) must be 2 or 3");
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(3) < 2 ||
                              outputData.dimension(3) > 3,
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
      const ordinal_type f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
      for (ordinal_type i=0; i<4; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataLeft");
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (getrank(inputDataRight) == 3) {
      // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(C,P,D):  all dimensions  C, P, D must match
      for (ordinal_type i=0; i<inputDataLeft.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(i) != inputDataRight.dimension(i),
                                  ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
      }

      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D): dimensions C, P, D must match
      {
        const ordinal_type f1[] = { 0, 1, 2, 3 }, f2[] = { 0, 1, 2, 2 };
        for (ordinal_type i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataRight");
        }
      }
    } else {
      // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(P,D): dimensions  P, D must match
      {
        const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): inputDataLeft dimension does not match with inputDataRight");
        }
      }

      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D): dimensions P, D must match
      {
        const ordinal_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 1 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::outerProductDataField): outputData dimension does not match with inputDataRight");
        }
      }
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>
  matvecProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<intputFieldProperties...> inputFields,
                          const char transpose ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
     *      (2) inputFields(C,F,P,D) or (F,P,D);
     *      (3) outputFields(C,F,P,D)
     */
    // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputData.rank() < 2 || inputData.rank() > 4,
                              ">>> ERROR (ArrayTools::matvecProductDataField): inputData must have rank 2,3 or 4" );
    if (inputData.rank() > 2) {
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(2) < 1 ||
                                inputData.dimension(2) > 3,
                                ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) must be 1,2 or 3");
    }
    if (inputData.rank() == 4) {
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(3) < 1 ||
                                inputData.dimension(3) > 3,
                                ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(3) must be 1,2 or 3");
    }

    // (2) inputFields is (C, F, P, D) or (F, P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 3 ||
                              inputFields.rank() > 4,
                              ">>> ERROR (ArrayTools::matvecProductDataField): inputFields must have rank 3 or 4" );
    INTREPID2_TEST_FOR_ABORT( inputFields.rank(inputfields.rank()-1) < 1 ||
                              inputFields.rank(inputfields.rank()-1) > 3,
                              ">>> ERROR (ArrayTools::matvecProductDataField): inputFields dimension(rank-1) must be 1,2, or 3" );

    // (3) outputFields is (C,F,P,D)
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 4,
                              ">>> ERROR (ArrayTools::matvecProductDataField): outputFields must have rank 4" );
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(3) < 1 ||
                              outputFields.dimension(3) > 3,
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
      INTREPID2_TEST_FOR_ABORT( outputFields.dimension(2) != inputData.dimension(1),
                                ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(2) must match to inputData dimension(1)" );
    }
    if (getrank(inputData) == 2) { // inputData(C,P) -> C match
      INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputData.dimension(0),
                                ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension(0) must match to inputData dimension(0)" );
    }
    if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
      const ordinal_type f1[] = { 0, 3 }, f2[] = { 0, 2 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
      }
    }
    if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
      const ordinal_type f1[] = { 0, 3, 3 }, f2[] = { 0, 2, 3 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputData dimension" );
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (inputFields.rank() == 4) {
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
      if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputFields.dimension(2) != inputData(1),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputFields dimension (2) does not match to inputData dimension (1)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(0) != inputFields.dimension(0),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
        const ordinal_type f1[] = { 0, 2 }, f2[] = { 0, 3 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }
      if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
        const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 3 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }

      // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C, F, P, D must match
      for (ordinal_type i=0; i<outputFields.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(i) != inputFields.dimension(i),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputFields dimension" );
      }
    } else {
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D): dimensions  P, D must match
      if (inputData.dimension(1) > 1) { // check P if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(1),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(1) does not match to inputFields dimension (1)" );
      }
      if (inputData.rank() == 3) {
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(2) != inputFields.dimension(2),
                                  ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension(2) does not match to inputFields dimension (2)" );
      }
      if (inputData.rank() == 4) {
        const ordinal_type f1[] = { 2, 3 }, f2[] = { 2, 2 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }

      // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F, P, D must match
      {
        const ordinal_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      }
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  matvecProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>    outputData,
                         const Kokkos::DynRankView<inputDataLeftProperties...> inputDataLeft,
                         const Kokkos::DynRankView<intputFieldProperties...>   inputDataRight,
                         const char transpose ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data
     *      (2) inputDataRight(C,P,D) or (P,D);
     *      (3) outputData(C,P,D)
     */
    // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() < 2 ||
                              inputDataLeft.rank() > 4,
                              ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft must have rank 2,3 or 4" );

    if (inputDataLeft.rank() > 2) {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) < 1 ||
                                inputDataLeft.dimension(2) > 3,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(2) must be 1, 2 or 3");
    }
    if (inputDataLeft.rank() == 4) {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(3) < 1 ||
                                inputDataLeft.dimension(3) > 3,
                                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension(3) must be 1, 2 or 3");
    }

    // (2) inputDataRight is (C, P, D) or (P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 2 ||
                              inputDataRight.rank() > 3,
                              ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight must have rank 2 or 3" );
    INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(inputDataRight.rank()-1) < 1 ||
                              inputDataRight.dimension(inputDataRight.rank()-1) > 3,
                              ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight dimension (rank-1) must be 1,2 or 3" );

    // (3) outputData is (C,P,D)
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 3,
                              ">>> ERROR (ArrayTools::matvecProductDataData): outputData must have rank 3" );
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(2) < 1 ||
                              outputData.dimension(2) > 3,
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
      INTREPID2_TEST_FOR_ABORT( outputData.dimension(1) != inputDataLeft.dimension(1),
                                ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(1) muat match inputDataLeft dimension(1)" );
    }
    if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
      INTREPID2_TEST_FOR_ABORT( outputData.dimension(0) != inputDataLeft.dimension(0),
                                ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension(0) muat match inputDataLeft dimension(0)" );
    }
    if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
      const ordinal_type f1[] = { 0, 2 }, f2[] = { 0, 2 };
      for (ordinal_type i=0; i<2; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
      }
    }
    if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
      const ordinal_type f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match inputDataLeft dimension" );
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (inputDataRight.rank() == 3) {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D):  dimensions  C, P, D must match
      if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) muat match to inputDataRight dimension (1)" );
      }
      if (getrank(inputDataLeft) == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (0) muat match to inputDataRight dimension (0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const ordinal_type f1[] = { 0, 2 }, f2[] = { 0, 2 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
        }
      }

      // Cross-check (3): outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
      for (ordinal_type i=0; i<outputData.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(i) != inputDataRight.dimension(i),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match to inputDataRight dimension" );
      }
    } else {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D): dimensions  P, D must match
      if (inputDataLeft.dimension(1) > 1) { // check P if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(0),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (1) does mot match to inputDataright dimension (0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check D
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) != inputDataRight.dimension(1),
                                  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension (2) does mot match to inputDataright dimension (1)" );
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
        const ordinal_type f1[] = { 2, 3 }, f2[] = { 1, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft dimension muat match to inputDataRight dimension" );
        }
      }

      // Cross-check (3): outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
      {
        const ordinal_type f1[] = { 1, 2 }, f2[] = { 0, 1 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matvecProductDataData): outputData dimension muat match to inputDataRight dimension" );
        }
      }
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldsProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  matmatProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<inputFieldProperties...>  inputFields,
                          const char transpose = 'N' ) {

#ifdef HAVE_INTREPID_DEBUG

    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
     *      (2) inputFields(C,F,P,D,D) or (F,P,D,D);
     *      (3) outputFields(C,F,P,D,D)
     */
    // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputData.rank() < 2 ||
                              inputData.rank() > 4,
                              ">>> ERROR (ArrayTools::matmatProductDataField): inputData must have rank 2,3 or 4" );
    if (inputData.rank() > 2) {
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(2) < 1 ||
                                inputData.dimension(2) > 3,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (2) must be 1,2 or 3" );
    }
    if (inputData.rank() == 4) {
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(3) < 1 ||
                                inputData.dimension(3) > 3,
                                ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (3) must be 1,2 or 3" );
    }

    // (2) inputFields is (C,F,P,D,D) or (F,P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 4 ||
                              inputFields.rank() > 5,
                              ">>> ERROR (ArrayTools::matmatProductDataField): inputFields must have rank 4 or 5");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(inputFields.rank()-1) < 1 ||
                              inputFields.dimension(inputFields.rank()-1) > 3,
                              ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-1) must be 1,2 or 3" );
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(inputFields.rank()-2) < 1 ||
                              inputFields.dimension(inputFields.rank()-2) > 3,
                              ">>> ERROR (ArrayTools::matmatProductDataField): inputFields dimension (rank-2) must be 1,2 or 3" );

    // (3) outputFields is (C,F,P,D,D)
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 5,
                              ">>> ERROR (ArrayTools::matmatProductDataField): outputFields must have rank 5" );
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(3) < 1 ||
                              outputFields.dimension(3) > 3,
                              ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (3) must be 1,2 or 3" );
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(4) < 1 ||
                              outputFields.dimension(4) > 3,
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
      INTREPID2_TEST_FOR_ABORT( outputFields.dimension(2) != inputData.dimension(1),
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (2) does not match to inputData dimension (1)" );
    }
    if (inputData.rank() == 2) { // inputData(C,P) -> C match
      INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputData.dimension(0),
                                ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension (0) does not match to inputData dimension (0)" );
    }
    if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match
      const ordinal_type f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 2 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
      }
    }
    if (inputData.rank() == 4) { // inputData(C,P,D,D) -> C, D, D match
      const ordinal_type f1[] = { 0, 3, 4 }, f2[] = { 0, 2, 3 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputData.dimension(f2[i]),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputData dimension" );
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (inputFields.rank() == 5) {
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D,D):  dimensions  C, P, D must match
      if (inputData.dimension(1) > 1) { // check P dimension if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (2)" );
      }
      if (inputData.rank() == 2) { // inputData(C,P) -> C match
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(0) != inputFields.dimension(0),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (0) does not match to inputFields dimension (0)" );
      }
      if (inputData.rank() == 3) { // inputData(C,P,D) -> C, D match

        const ordinal_type f1[] = { 0, 2, 2 }, f2[] = { 0, 3, 4 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }
      if (inputData.rank() == 4) {  // inputData(C,P,D,D) -> C, D, D match
        const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 3, 4 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }

      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D,D): all dimensions C, F, P, D must match
      for (ordinal_type i=0; i<outputFields.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( outputFields.dimension(i) != inputFields.dimension(i),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputFields dimension" );
      }
    } else {
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D,D): dimensions  P, D must match
      if (inputData.dimension(1) > 1) { // check P if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(1),
                                  ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension (1) does not match to inputFields dimension (1)" );
      }
      if (inputData.rank() == 3) {
        const ordinal_type f1[] = { 2, 2 }, f2[] = { 2, 3 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }
      if (inputData.rank() == 4) {
        const ordinal_type f1[] = { 2, 3 }, f2[] = { 2, 3 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputData.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matmatProductDataField): inputData dimension does not match to inputFields dimension" );
        }
      }

      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D,D): dimensions F, P, D must match
      {
        const ordinal_type f1[] = { 1, 2, 3, 4 }, f2[] = { 0, 1, 2, 3 };
        for (ordinal_type i=0; i<4; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputFields.dimension(f1[i]) != inputFields.dimension(f2[i]),
                                    ">>> ERROR (ArrayTools::matmatProductDataField): outputFields dimension does not match to inputFields dimension" );
        }
      }
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>
  matmatProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                         const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                         const char transpose  ) {

#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data
     *      (2) inputDataRight(C,P,D,D) or (P,D,D);
     *      (3) outputData(C,P,D,D)
     */
    // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() < 2 ||
                              inputDataLeft.rank() > 4,
                              "h:: inputDataLeft must have rank 2,3 or 4" );
    if (inputDataLeft.rank() > 2) {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) < 1 || 
                                inputDataLeft.dimension(2) > 3,
                                "h:: inputDataLeft dimension (2) must be 1,2 or 3" );
    }
    if (inputDataLeft.rank() == 4) {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(3) < 1 || 
                                inputDataLeft.dimension(3) > 3,
                                "h:: inputDataLeft dimension (3) must be 1,2 or 3" );
    }

    // (2) inputDataRight is (C,P,D,D) or (P,D,D) and 1 <= (D=dimension(rank-1),(rank-2)) <= 3 is required.
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 3 ||
                              inputDataRight.rank() > 4,
                              "h:: inputDataRight must have rank 3 or 4" );
    INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(inputDataRight.rank()-1) < 1 ||
                              inputDataRight.dimension(inputDataRight.rank()-1) > 3 ||
                              "h:: inputDataRight (rank-1) must be 1,2 or 3" );
    INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(inputDataRight.rank()-2) < 1 ||
                              inputDataRight.dimension(inputDataRight.rank()-2) > 3 ||
                              "h:: inputDataRight (rank-2) must be 1,2 or 3" );

    // (3) outputData is (C,P,D,D)
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 4,
                              "h:: outputData must have rank 4" );                              
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(2) < 1 ||
                              outputData.dimension(2) > 3 ||
                              "h:: outputData (2) must be 1,2 or 3" );
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(3) < 1 ||
                              outputData.dimension(3) > 3 ||
                              "h:: outputData (3) must be 1,2 or 3" );

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
      INTREPID2_TEST_FOR_ABORT( outputData.dimension(1) != inputDataLeft.dimension(1),
                                "h:: outputData dimension (1) does not match to inputDataLeft dimension (1)" );                                
    }
    if (inputDataLeft.rank() == 2) { // inputDataLeft(C,P): check C
      INTREPID2_TEST_FOR_ABORT( outputData.dimension(0) != inputDataLeft.dimension(0),
                                "h:: outputData dimension (0) does not match to inputDataLeft dimension (0)" );                                
    }
    if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
      const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 2 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]),
                                    "h:: outputData dimension does not match to inputDataLeft dimension" );
      }
    }
    if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
      const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
      for (ordinal_type i=0; i<3; ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataLeft.dimension(f2[i]),
                                  "h:: outputData dimension does not match to inputDataLeft dimension" );
      }
    }

    /*
     *   Cross-checks (1,3):
     */
    if (inputDataRight.rank() == 4) {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D,D):  dimensions  C, P, D must match
      if (inputDataLeft.dimension(1) > 1) { // check P dimension if P>1 in inputDataLeft
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1),
                                  "h:: inputDataLeft dimension (1) does not match to inputDataRight dimension (1)" );
      }
      if (getrank(inputDataLeft) == 2) { // inputDataLeft(C,P): check C
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0),
                                  "h:: inputDataLeft dimension (0) does not match to inputDataRight dimension (0)" );
      }
      if (inputDataLeft.rank() == 3) {  // inputDataLeft(C,P,D): check C and D
        const ordinal_type f1[] = { 0, 2, 2 }, f2[] = { 0, 2, 3 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    "h:: intputDataLeft dimension does not match to inputDataRight dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check C and D
        const ordinal_type f1[] = { 0, 2, 3 }, f2[] = { 0, 2, 3 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    "h:: intputDataLeft dimension does not match to inputDataRight dimension" );
        }
      }

      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D,D): all dimensions C, P, D must match
      for (ordinal_type i=0; i< outputData.rank(); ++i) {
        INTREPID2_TEST_FOR_ABORT( outputData.dimension(i) !=  inputDataRight.dimension(i),
                                  "h:: outputData dimension does not match to inputDataRight dimension" );
    } else {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D,D): dimensions  P, D must match
      if (inputDataLeft.dimension(1) > 1) { // check P if P>1 in inputData
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(0),
                                  "h:: inputDataLeft dimension (1) does not match to inputDataRight dimension (0)" );
      }
      if (getrank(inputDataLeft) == 3) {  // inputDataLeft(C,P,D): check D
        const ordinal_type f1[] = { 2, 2 }, f2[] = { 1, 2 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    "h:: intputDataLeft dimension does not match to inputDataRight dimension" );
        }
      }
      if (inputDataLeft.rank() == 4) {  // inputDataLeft(C,P,D,D): check D
        const ordinal_type f1[] = { 2, 3 }, f2[] = { 1, 2 };
        for (ordinal_type i=0; i<2; ++i) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    "h:: intputDataLeft dimension does not match to inputDataRight dimension" );
        }
      }
      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D,D): dimensions P, D must match
      {
        const ordinal_type f1[] = { 1, 2, 3 }, f2[] = { 0, 1, 2 };
        for (ordinal_type i=0; i<3; ++i) {
          INTREPID2_TEST_FOR_ABORT( outputData.dimension(f1[i]) != inputDataRight.dimension(f2[i]),
                                    "h:: outputData dimension does not match to inputDataRight dimension" );
        }
      }
    }
#endif
    //body
  }

} // end namespace Intrepid
#endif
