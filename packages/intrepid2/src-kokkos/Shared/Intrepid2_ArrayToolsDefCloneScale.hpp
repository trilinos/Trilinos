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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_ArrayToolsDefCloneScale.hpp
    \brief  Definition file for clone / scale operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void ArrayTools<ExecSpaceType>::
  cloneFields( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
               const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

    //static_assert
    // - two memory space of input and output fields should be matched

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( ( inputFields.rank() < 2 || inputFields.rank() > 4 ),
                              ">>> ERROR (ArrayTools::cloneFields): Input fields container must have rank 2, 3, or 4.");
    INTREPID2_TEST_FOR_ABORT( ( outputFields.rank() != (inputFields.rank()+1) ),
                              ">>> ERROR (ArrayTools::cloneFields): The rank of the input fields container must be one less than the rank of the output fields container.");
    for (ordinal_type i=0;i<inputFields.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( (inputFields.dimension(i) != outputFields.dimension(i+1)),
                                ">>> ERROR (ArrayTools::cloneFields): Dimensions of input and output fields containers do not match.");
    }
#endif

    struct Functor {
      Kokkos::DynRankView<outputFieldProperties...> _outputFields;
      Kokkos::DynRankView<intputFieldProperties...> _inputFields;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldProperties...> &outputFields_,
              Kokkos::DynRankView<intputFieldProperties...> &inputFields_)
        : _outputFields(outputFields_), _inputFields(inputFields_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) {
        Kokkos::deep_copy(Kokkos::subview(_outputFields, cl /* , Kokkos::All() */ ), _inputFields);
      }
    };

    const size_type numCells = outputFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(outputFields, inputFields) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputFactorProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void ArrayTools<ExecSpaceType>::
  cloneScaleFields( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                    const Kokkos::DynRankView<inputFactorProperties...> inputFactors,
                    const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

    // static assert
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputFactors.rank() != 2,
                              ">>> ERROR (ArrayTools::cloneScaleFields): The rank of the input factors container must be 2.");
    INTREPID2_TEST_FOR_ABORT( inputFields.rank() < 2 || inputFields.rank() > 4,
                              ">>> ERROR (ArrayTools::cloneScaleFields): Input fields container must have rank 2, 3, or 4.");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != (inputFields.rank()+1),
                              ">>> ERROR (ArrayTools::cloneScaleFields): The rank of the input fields container must be one less than the rank of the output fields container.");
    INTREPID2_TEST_FOR_ABORT( inputFactors.dimension(0) != outputFields.dimension(0),
                              ">>> ERROR (ArrayTools::cloneScaleFields): Zeroth dimensions of input factors container and output fields container (numbers of integration domains) must agree!");
    INTREPID2_TEST_FOR_ABORT( inputFactors.dimension(1) != outputFields.dimension(1),
                              ">>> ERROR (ArrayTools::cloneScaleFields): First dimensions of input factors container and output fields container (numbers of fields) must agree!");
    for (ordinal_type i=0;i<inputFields.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( (inputFields.dimension(i) != outputFields.dimension(i+1)),
                                ">>> ERROR (ArrayTools::cloneScaleFields): Dimensions of input and output fields containers do not match.");
    }
#endif

    struct Functor {
      Kokkos::DynRankView<outputFieldProperties...> _outputFields;
      Kokkos::DynRankView<inputFactorProperties...> _inputFactors;
      Kokkos::DynRankView<intputFieldProperties...> _inputFields;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldProperties...> &outputFields_,
              Kokkos::DynRankView<inputFactorProperties...> &inputFactors_,
              Kokkos::DynRankView<intputFieldProperties...> &inputFields_)
        : _outputFields(outputFields_), _inputFactors(inputFactors_), _inputFields(inputFields_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) {
        const ordinal_type numFields = outputFields.dimension(1);
        const ordinal_type numPoints = outputFields.dimension(2);
        const ordinal_type dim1Tens  = outputFields.dimension(3);
        const ordinal_type dim2Tens  = outputFields.dimension(4);

        for(ordinal_type bf = 0; bf < numFields; ++bf)
          for(ordinal_type pt = 0; pt < numPoints; ++pt)
            for(ordinal_type iTens1 = 0; iTens1 < dim1Tens; ++iTens1)
              for(ordinal_type iTens2 = 0; iTens2 < dim2Tens; ++iTens2)
                outputFieldsWrap(cl, bf, pt, iTens1, iTens2)
                  = inputFieldsWrap(bf, pt, iTens1, iTens2) * inputFactorswrap(cl, bf);
      }
    };

    const size_type numCells = outputFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(outputFields, inputFactors, inputFields) );
  }

  template<typename ExecSpaceType>
  template<class ...inoutFieldProperties,
           class ...inputFactorProperties>
  KOKKOS_INLINE_FUNCTION
  void ArrayTools<ExecSpaceType>::
  scaleFields( /**/  Kokkos::DynRankView<inoutFieldProperties...>  inoutFields,
               const Kokkos::DynRankView<inputFactorProperties...> inputFactors ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputFactors.rank() != 2,
                              ">>> ERROR (ArrayTools::scaleFields): The rank of the input factors container must be 2.");
    INTREPID2_TEST_FOR_ABORT( inoutFields.rank() < 3 || inoutFields.rank() > 5,
                              ">>> ERROR (ArrayTools::scaleFields): Input/output fields container must have rank 3, 4, or 5.");
    INTREPID2_TEST_FOR_ABORT( inputFactors.dimension(0) != inoutFields.dimension(0),
                              ">>> ERROR (ArrayTools::scaleFields): Zeroth dimensions of input factors container and input/output fields container (numbers of integration domains) must agree!");
    INTREPID2_TEST_FOR_ABORT( inputFactors.dimension(1) != inoutFields.dimension(1),
                              ">>> ERROR (ArrayTools::scaleFields): First dimensions (number of fields) of input factors and input/output fields containers must agree!");
#endif

    struct Functor {
      Kokkos::DynRankView<inoutFieldProperties...>  _inoutFields;
      Kokkos::DynRankView<inputFactorProperties...> _inputFactors;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<inoutFieldProperties...>  &inoutFields_,
              Kokkos::DynRankView<inputFactorProperties...> &inputFactors_)
        : _inoutFields(inoutFields_), _inputFactors(inputFactors_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) {
        const ordinal_type numFields = outputFields.dimension(1);
        const ordinal_type numPoints = outputFields.dimension(2);
        const ordinal_type dim1Tens  = outputFields.dimension(3);
        const ordinal_type dim2Tens  = outputFields.dimension(4);

        for(ordinal_type bf = 0; bf < numFields; ++bf)
          for(ordinal_type pt = 0; pt < numPoints; ++pt)
            for(ordinal_type iTens1 = 0; iTens1 < dim1Tens; ++iTens1)
              for(ordinal_type iTens2 = 0; iTens2 < dim2Tens; ++iTens2)
                inoutFields(cl, bf, pt, iTens1, iTens2)
                  = inoutFields(cl, bf, pt, iTens1, iTens2) * inputFactors(cl, bf);
      }
    };

    const size_type numCells = outputFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(inoutFactors, inputFields) );
  }

} // end namespace Intrepid2
