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

namespace Intrepid {

  template<typename ExecSpaceType>
  template<class ...outputProperties,
           class ...leftInputProperties,
           class ...rightInputProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  dotMultiply( /**/  Kokkos::DynRankView<outputProperties...>     output,
               const Kokkos::DynRankView<leftInputProperties...>  leftInput,
               const Kokkos::DynRankView<rightInputProperties...> rightInput ) {
    struct Functor {
      /**/  Kokkos::DynRankView<outputProperties...>     _output;
      const Kokkos::DynRankView<leftInputProperties...>  _leftInput;
      const Kokkos::DynRankView<rightInputProperties...> _rightInput;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputProperties...>     &output_,
              Kokkos::DynRankView<leftInputProperties...>  &leftInput_,
              Kokkos::DynRankView<rightInputProperties...> &rightInput_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl) {
        auto result = Kokkos::subdynrankview(_output, cl, Kokkos::ALL(), Kokkos::ALL());

        const auto left = Kokkos::subdynrankview(_leftInput, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        const auto right = ( _rightInput.rank() == (_leftInput.rank() +1) ?
                             Kokkos::subdynrankview(_rightInput, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) :
                             Kokkos::subdynrankview(_rightInput,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );

        const size_type numFields = result.dimension(0);
        const size_type numPoints = result.dimension(1);
        const size_type dim1Tens  = left.dimension(1);
        const size_type dim2Tens  = left.dimension(2);

        temp += inputDataWrap(cl, pt, iTens1, iTens2)*inputFieldsWrap(cl, bf, pt, iTens1, iTens2);
        outputFieldsWrap(cl, bf, pt) =  temp;

          for(size_type bf = 0; bf < numFields; ++bf)
            for(size_type pt = 0; pt < numPoints; ++pt) {
              value_type tmp(0);
              for(size_type iTens1 = 0; iTens1 < dim1Tens; ++iTens1)
                for(size_type iTens2 = 0; iTens2 < dim2Tens; ++iTens2)
                  temp += left(pt, iTens1, iTens2)*right(bf, pt, iTens1, iTens2);
              result(bf, pt) = tmp;
            }

      }
    };

    const size_type numCells = output.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(output, leftInput, rightInput) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  dotMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                        const Kokkos::DynRankView<inputDataProperties...>   inputData,
                        const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID_DEBUG
    if (inputFields.rank() > inputData.rank()) {
      INTREPID2_TEST_FOR_ABORT( inputData.rank() < 2 || inputData.rank() > 4,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_ABORT( inputFields.rank() != (inputData.rank()+1),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Input fields container must have rank one larger than the rank of the input data container.");
      INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 3,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
      INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != inputData.dimension(0),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                                inputData.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      for (size_type i=2;i<inputData.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(i) != inputFields.dimension(i+1),
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i+1) of inputFields");
      }
      for (size_t i=0;i<outputFields.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputFields.dimension(i) != outputFields.dimension(i),
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
      }
    } else {
      INTREPID2_TEST_FOR_ABORT( inputData.rank() < 2 || inputData.rank() > 4,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Input data container must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_ABORT( inputFields.rank() != inputData.rank(),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): The rank of fields input container must equal the rank of data input container.");
      INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 3,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Output fields container must have rank 3.");
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(1) &&
                                inputData.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimensions of the fields and data input containers (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != outputFields.dimension(1),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimension of the fields input container and first dimension of the fields output container (number of fields) must agree!");
      INTREPID2_TEST_FOR_ABORT( inputFields.dimension(1) != outputFields.dimension(2),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): First dimension of the fields input container and second dimension of the fields output container (number of integration points) must agree!");
      INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputData.dimension(0),
                                ">>> ERROR (ArrayTools::dotMultiplyDataField): Zeroth dimensions of the fields output and data input containers (number of integration domains) must agree!");
      for (size_type i=2;i<inputData.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputData.dimension(i) != inputFields.dimension(i+1),
                                  ">>> ERROR (ArrayTools::dotMultiplyDataField): inputData dimension (i) does not match to the dimension (i+1) of inputFields");
      }
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::dotMultiply( outputFields,
                                                      inputData,
                                                      inputFields );
  }


  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  dotMultiplyDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                       const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    if (inputDataRight.rank() >= inputDataLeft.rank()) {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() != inputDataLeft.rank(),
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): The rank of the right data input container must equal the rank of the left data input container.");
      INTREPID2_TEST_FOR_ABORT( outputData.rank() != 2,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1) &&
                                inputDataLeft.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first left data dimension must be 1!");
      for (size_type i=0;i<inputDataLeft.rank();++i) {
        if (i != 1) {
          INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(i) != inputDataRight.dimension(i),
                                    ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i) does not match to the dimension (i) of inputDataRight");
        }
      }
      for (size_type i=0;i<outputData.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(i) != outputData.dimension(i),
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
      }
    } else {
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() < 2 || inputDataLeft.rank() > 4,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Left data input container must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() != (inputDataLeft()-1),
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Right data input container must have rank one less than the rank of left data input container.");
      INTREPID2_TEST_FOR_ABORT( outputData.rank() != 2,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Data output container must have rank 2.");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(0) &&
                                inputDataLeft.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of left data input container (number of integration points) must agree or first left data dimension must be 1!");
      INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(0) != outputData.dimension(1),
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimension of the right data input container and first dimension of output data container (number of integration points) must agree!");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != outputData.dimension(0),
                                ">>> ERROR (ArrayTools::dotMultiplyDataData): Zeroth dimensions of the left data input and data output containers (number of integration domains) must agree!");
      for (size_type i=1;i<inputDataRight.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(i+1) != inputDataRight.dimension(i),
                                  ">>> ERROR (ArrayTools::dotMultiplyDataData): inputDataLeft dimension (i+1) does not match to the dimension (i) of inputDataRight");
      }
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::dotMultiply( outputData,
                                                      inputDataLeft,
                                                      inputDataRight );
  }

} // end namespace Intrepid2
