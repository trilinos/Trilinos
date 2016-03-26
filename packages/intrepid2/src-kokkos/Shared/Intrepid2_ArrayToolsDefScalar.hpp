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

  template<typename ExecSpaceType>
  template<class ...outputProperties,
           class ...leftInputProperties,
           class ...righInputtProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  scalarMultiply( /**/  Kokkos::DynRankView<outputProperties...>     output,
                  const Kokkos::DynRankView<leftInputProperties...>  leftInput,
                  const Kokkos::DynRankView<rightInputProperties...> rightInput,
                  const bool reciprocal = false ) {

    struct Functor {
      /**/  Kokkos::DynRankView<outputProperties...>     _output;
      const Kokkos::DynRankView<leftInputProperties...>  _leftInput;
      const Kokkos::DynRankView<rightInputProperties...> _rightInput;
      const bool _reciprocal;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputProperties...>     &output_,
              Kokkos::DynRankView<leftInputProperties...>  &leftInput_,
              Kokkos::DynRankView<rightInputProperties...> &rightInput_,
              const bool reciprocal_)
        : _output(output_), _leftInput(leftInput_), _rightInput(rightInput_), _reciprocal(reciprocal_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl) {
        auto result = Kokkos::subdynrankview(_output, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        const auto left = ( _leftInput.rank() == _output.rank() ?
                            Kokkos::subdynrankview(_leftInput, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) :
                            Kokkos::subdynrankview(_leftInput,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );

        const auto right = Kokkos::subdynrankview(_rightInput, cl, Kokkos::ALL());

        const size_type numFields = result.dimension(0);
        const size_type numPoints = result.dimension(1);
        const size_type dim1Tens  = result.dimension(2);
        const size_type dim2Tens  = result.dimension(3);

        if (_reciprocal)
          for(size_type bf = 0; bf < numFields; ++bf)
            for(size_type pt = 0; pt < numPoints; ++pt)
              for(size_type iTens1 = 0; iTens1 < dim1Tens; ++iTens1)
                for(size_type iTens2 = 0; iTens2 < dim2Tens; ++iTens2)
                  result(bf, pt, iTens1, iTens2) = left(bf, pt, iTens1, iTens2)/right(pt);
        else
          for(size_type bf = 0; bf < numFields; ++bf)
            for(size_type pt = 0; pt < numPoints; ++pt)
              for(size_type iTens1 = 0; iTens1 < dim1Tens; ++iTens1)
                for(size_type iTens2 = 0; iTens2 < dim2Tens; ++iTens2)
                  result(bf, pt, iTens1, iTens2) = left(bf, pt, iTens1, iTens2)*right(pt);
      }
    };

    const size_type numCells = output.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(output, leftInput, rightInput, reciprocal) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  scalarMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldProperties...>  inputFields,
                           const bool reciprocal = false ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 2,
                              ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input data container must have rank 2.");

    if (outputFields.rank() <= inputFields.rank()) {
      INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 3 || inputFields.rank() > 5,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 3, 4, or 5.");
      INTREPID2_TEST_FOR_ABORT( outputFields.rank() != inputFields.rank(),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input and output fields containers must have the same rank.");
      INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != inputData.dimension(0),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                                inputData.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree, or first data dimension must be 1!");
      for (size_t i=0;i<inputFields.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputFields.dimension(i) != outputFields.dimension(i),
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i) of outputFields");
      }
    }
    else {
      INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 2 || inputFields.rank() > 4,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 2, 3, or 4.");
      INTREPID2_TEST_FOR_ABORT( outputFields.rank() != (inputFields.rank()+1)
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): The rank of the input fields container must be one less than the rank of the output fields container.");
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(1) &&
                                inputData.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_ABORT( inputData.dimension(0) != outputFields.dimension(0),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
      for (size_type i=0;i<inputFields.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputFields.dimension(i) != outputFields.dimension(i+1),
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
      }
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::scalarMultiply( outputFields,
                                                         inputData,
                                                         inputFields,
                                                         reciprocal );
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  scalarMultiplyDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                          const bool reciprocal = false ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank() != 2,
                              ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");

    if (outputData.rank() <= inputDataRight.rank()) {
      INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 2 || inputDataRight.rank() > 4,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2, 3, or 4.");
      INTREPID2_TEST_FOR_ABORT( outputData.rank() != inputDataRight.rank(),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
      INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(0) != inputDataLeft.dimension(0),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1) &&
                                inputDataLeft.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree, or first dimension of the left data input container must be 1!");
      for (size_type i=0;i<inputDataRight.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(i) != outputData.dimension(i),
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
      }
    } else {
      INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() < 1 || inputDataRight.rank() > 3,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1, 2, or 3.");
      INTREPID2_TEST_FOR_ABORT( outputData.rank() != (inputDataRight.rank()+1),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(0)  &&
                                inputDataLeft.dimension(1) != 1,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
      INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != outputData.dimension(0),
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
      for (size_type i=0;i<inputDataRight.rank();++i) {
        INTREPID2_TEST_FOR_ABORT( inputDataRight.dimension(i) != outputData.dimension(i+1),
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i+1) of outputData");
      }
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::scalarMultiply( outputData,
                                                         inputDataLeft,
                                                         inputDataRight,
                                                         reciprocal );
  }

} // end namespace Intrepid2

#endif
