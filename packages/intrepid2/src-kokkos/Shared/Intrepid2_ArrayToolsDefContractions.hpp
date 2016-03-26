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

/** \file   Intrepid2_ArrayToolsDefContractions.hpp
    \brief  Definition file for contraction (integration) operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_CONTRACTIONS_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_CONTRACTIONS_HPP__

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...leftFieldProperties,
           class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  contractFieldField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                      const Kokkos::DynRankView<leftFieldProperties...>   leftFields,
                      const Kokkos::DynRankView<rightFieldProperties...>  rightFields,
                      const bool sumInto ) {
    struct Functor {
      /**/  Kokkos::DynRankView<outputFieldProperties...> _outputFields;
      const Kokkos::DynRankView<leftFieldProperties...>   _leftFields;
      const Kokkos::DynRankView<rightFieldProperties...>  _rightFields;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldProperties...> &outputFields_,
              Kokkos::DynRankView<leftFieldProperties...>   &leftFields_,
              Kokkos::DynRankView<rightFieldProperties...>  &rightFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _leftFields(leftFields_), _rightFields(rightFields_), _sumInto(sumInto_) {}
      
      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl) {
        const size_type numLeftInput  = _leftFields.dimension(1);
        const size_type numRightInput = _rightFields.dimension(1);
        const size_type numPoints     = _leftFields.dimension(2);
        const size_type dim1Tensor    = _leftFields.dimension(3);
        const size_type dim2Tensor    = _leftFields.dimension(4);

        if (_sumInto) {
          for (size_type lbf = 0; lbf < numLeftInput; ++lbf) 
            for (size_type rbf = 0; rbf < numRightInput; ++rbf) {
              value_type tmp(0);        
              for (size_type qp = 0; qp < numPoints; ++qp) 
                for (size_type iTens1 = 0; iTens1 < dim1Tensor; ++iTens1) 
                  for (size_type iTens2 = 0; iTens2 < dim2Tensor; ++iTens2) 
                    tmp += _leftFields(cl, lbf, qp, iTens1, iTens2)*_rightFields(cl, rbf, qp, iTens1, iTens2);
              _outputFields(cl, lbf, rbf) += tmp;
            }
        } else {
          for (size_type lbf = 0; lbf < numLeftInput; ++lbf) 
            for (size_type rbf = 0; rbf < numRightInput; ++rbf) {
              value_type tmp(0);        
              for (size_type qp = 0; qp < numPoints; ++qp) 
                for (size_type iTens1 = 0; iTens1 < dim1Tensor; ++iTens1) 
                  for (size_type iTens2 = 0; iTens2 < dim2Tensor; ++iTens2) 
                    tmp += _leftFields(cl, lbf, qp, iTens1, iTens2)*_rightInmput(cl, rbf, qp, iTens1, iTens2);
              _outputFields(cl, lbf, rbf) = tmp;
            }
        }
      }
    };
    const size_type numCells = leftFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(outputFields, leftFields, rightFields, sumInto) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  contractDataField( /**/  Kokkos::DynRankView<outputFieldProperties...>      outputFields,
                     const Kokkos::DynRankView<inputDataProperties...>   inputData,
                     const Kokkos::DynRankView<inputFieldsFieldProperties...> inputFields,
                     const bool sumInto ) {
    struct Functor {
      /**/  Kokkos::DynRankView<outputFieldProperties...>      _outputFields;
      const Kokkos::DynRankView<inputDataProperties...>   _inputData;
      const Kokkos::DynRankView<inputFieldsFieldProperties...> _inputFields;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldProperties...>      &outputFields_,
              Kokkos::DynRankView<inputDataProperties...>   &inputData_,
              Kokkos::DynRankView<inputFieldsFieldProperties...> &inputFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _inputData(inputData_), _inputFields(inputFields_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl) {
        const size_type numFields = _inputFields.dimension(1);
        const size_type numPoints = _inputFields.dimension(2);
        const size_type dim1Tens  = _inputFields.dimension(3);
        const size_type dim2Tens  = _inputFields.dimension(4);

        if (_sumInto) {
          for (size_type bf = 0; bf < numFields; ++bf) {
            value_type tmp(0);        
            for (size_type qp = 0; qp < numPoints; ++qp) 
              for (size_type iTens1 = 0; iTens1 < dim1Tensor; ++iTens1) 
                for (size_type iTens2 = 0; iTens2 < dim2Tensor; ++iTens2) 
                  tmp += _inputFields(cl, bf, qp, iTens1, iTens2)*_inputData(cl, qp, iTens1, iTens2);
            _outputFields(cl, bf) += tmp;
          }
        } else {
          for (size_type bf = 0; bf < numFields; ++bf) {
            value_type tmp(0);        
            for (size_type qp = 0; qp < numPoints; ++qp) 
              for (size_type iTens1 = 0; iTens1 < dim1Tensor; ++iTens1) 
                for (size_type iTens2 = 0; iTens2 < dim2Tensor; ++iTens2) 
                  tmp += _inputFields(cl, bf, qp, iTens1, iTens2)*_inputData(cl, qp, iTens1, iTens2);
            _outputFields(cl, bf) = tmp;
          }
        }
      }
    };
    const size_type numCells = inputData.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(outputFields, inputData, inputFields, sumInto) );
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::Internal::
  contractDataData( /**/  Kokkos::DynRankView<outputDataProperties...>          outputData,
                    const Kokkos::DynRankView<inputDataLeftFieldProperties...>  inputDataLeft,
                    const Kokkos::DynRankView<inputDataRightFieldProperties...> inputDataRight,
                    const bool sumInto ) {
    
    struct Functor {
      /**/  Kokkos::DynRankView<outputDataProperties...>          _outputData;
      const Kokkos::DynRankView<inputDataLeftFieldProperties...>  _inputDataLeft;
      const Kokkos::DynRankView<inputDataRightFieldProperties...> _inputDataRight;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputDataProperties...>          &outputData_,
              Kokkos::DynRankView<inputDataLeftFieldProperties...>  &inputDataLeft_,
              Kokkos::DynRankView<inputDataRightFieldProperties...> &inputDataRight_,
              const bool sumInto_) 
        : _outputData(outputData_), _inputDataLeft(inputDataLeft_), _inputDataRight(inputDataRight_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl) {
        size_type numPoints  = _inputDataLeft.dimension(1);
        size_type dim1Tensor = _inputDataLeft.dimension(2);
        size_type dim2Tensor = _inputDataLeft.dimension(3);

        value_type tmp(0);        
        for (size_type qp = 0; qp < numPoints; ++qp) 
          for (size_type iTens1 = 0; iTens1 < dim1Tensor; ++iTens1) 
            for (size_type iTens2 = 0; iTens2 < dim2Tensor; ++iTens2) 
              tmp += _inputDataLeft(cl, qp, iTens1, iTens2)*_inputDataRight(cl, qp, iTens1, iTens2);
        
        if (_sumInto) 
          _outputData(cl) += tmp;
        else 
          _outputData(cl) = tmp;
      }
    };
    const size_type numCells = inputDataLeft.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
    Kokkos::parallel_for( policy, Functor(outputData, inputDataLeft, inputDataRight, sumInto) );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...leftFieldProperties,
           class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldScalar( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( leftFields.rank() != 3,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( rightFields.rank() != 3, 
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 3,
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(0) != rightFields.dimension(0),
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(2) != rightFields.dimension(2),
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != rightFields.dimension(0),
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != leftFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(2) != rightFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...leftFieldProperties,
           class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldVector( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( leftFields.rank() != 4, 
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
    INTREPID2_TEST_FOR_ABORT( rightFields.rank() != 4,
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 3, 
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(0) != rightFields.dimension(0),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(2) != rightFields.dimension(2),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(3) != rightFields.dimension(3),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != rightFields.dimension(0),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != leftFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(2) != rightFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...leftFieldProperties,
           class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldTensor( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldProperties...>  rightFields,
                            const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( leftFields.rank()  != 5, 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
    INTREPID2_TEST_FOR_ABORT( rightFields.rank() != 5, 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 3,
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(0) != rightFields.dimension(0), 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(2) != rightFields.dimension(2), 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(3) != rightFields.dimension(3), 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( leftFields.dimension(4) != rightFields.dimension(4), 
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != rightFields.dimension(0),
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != leftFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(2) != rightFields.dimension(1),
                              ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataFieldScalar( /**/  Kokkos::DynRankView<outputFieldProperties...>  outputFields,
                           const Kokkos::DynRankView<inputDataProperties...>    inputData,
                           const Kokkos::DynRankView<intputFieldProperties...>  inputFields,
                           const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputFields.rank()  != 3, 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 2, 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 2, 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != inputData.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2) && 
                              inputData.dimension(1) != 1, 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputFields.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != inputFields.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  } 
  
  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataFieldVector( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<intputFieldProperties...> inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputFields.rank()  != 4, 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 3, 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 2, 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != inputData.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                              inputData.dimension(1) != 1, 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(3) != inputData.dimension(2), 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputFields.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != inputFields.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }
  
  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataFieldTensor( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<intputFieldProperties...> inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputFields.rank()  != 5, 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
    INTREPID2_TEST_FOR_ABORT( inputData.rank() != 4, 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
    INTREPID2_TEST_FOR_ABORT( outputFields.rank() != 2, 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(0) != inputData.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                              inputData.dimension(1) != 1, 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(3) != inputData.dimension(2), 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
    INTREPID2_TEST_FOR_ABORT( inputFields.dimension(4) != inputData.dimension(3), 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(0) != inputFields.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputFields.dimension(1) != inputFields.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataDataScalar( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank()  != 2, 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() != 2, 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 1, 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  } 
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataDataVector( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank()  != 3, 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() != 3, 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 1, 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) != inputDataRight.dimension(2), 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  contractDataDataTensor( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.rank()  != 4, 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
    INTREPID2_TEST_FOR_ABORT( inputDataRight.rank() != 4, 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
    INTREPID2_TEST_FOR_ABORT( outputData.rank() != 1, 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(2) != inputDataRight.dimension(2), 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( inputDataLeft.dimension(3) != inputDataRight.dimension(3), 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
    INTREPID2_TEST_FOR_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), 
                              ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }
  
}
#endif
