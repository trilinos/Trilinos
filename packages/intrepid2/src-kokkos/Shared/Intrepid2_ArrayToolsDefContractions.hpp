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
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::Internal::
  contractFieldField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                      const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                      const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                      const bool sumInto ) {
    typedef leftFieldValueType value_type;

    struct Functor {
      Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> _outputFields;
      Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   _leftFields;
      Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  _rightFields;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields_,
              Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields_,
              Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _leftFields(leftFields_), _rightFields(rightFields_), _sumInto(sumInto_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, lbf, rbf;
        Util::unrollIndex( cl, lbf, rbf, 
                _leftFields.dimension(0),
                _leftFields.dimension(1),
                iter );

        auto result = Kokkos::subdynrankview( _outputFields,  cl, lbf, rbf );

        const auto left  = Kokkos::subdynrankview( _leftFields,  cl, lbf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto right = Kokkos::subdynrankview( _rightFields, cl, rbf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        const size_type npts = left.dimension(0);
        const size_type iend = left.dimension(1);
        const size_type jend = left.dimension(2);

        value_type tmp(0);        
        for (size_type qp = 0; qp < npts; ++qp) 
          for (size_type i = 0; i < iend; ++i) 
            for (size_type j = 0; j < jend; ++j) 
              tmp += left(qp, i, j)*right(qp, i, j);

        result() = value_type(_sumInto)*result() + tmp;
      }
    };
    const size_type loopSize = leftFields.dimension(0)*leftFields.dimension(1)*rightFields.dimension(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(outputFields, leftFields, rightFields, sumInto) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::Internal::
  contractDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>      outputFields,
                     const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>        inputData,
                     const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...> inputFields,
                     const bool sumInto ) {
    typedef inputFieldValueType value_type;

    struct Functor {
      Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> _outputFields;
      Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   _inputData;
      Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  _inputFields;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields_,
              Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData_,
              Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _inputData(inputData_), _inputFields(inputFields_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, bf;
        Util::unrollIndex( cl, bf, 
                     _inputFields.dimension(0),
                     iter );
        
        auto result = Kokkos::subdynrankview( _outputFields, cl, bf );

        const auto field = Kokkos::subdynrankview( _inputFields, cl, bf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto data  = Kokkos::subdynrankview( _inputData,   cl,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        const size_type npts = field.dimension(0);
        const size_type iend = field.dimension(1);
        const size_type jend = field.dimension(2);

        value_type tmp(0);        

        if(data.dimension(1) != 1)
          for (size_type qp = 0; qp < npts; ++qp)
            for (size_type i = 0; i < iend; ++i)
              for (size_type j = 0; j < jend; ++j)
                tmp += field(qp, i, j) * data(qp, i, j);
        else
          for (size_type qp = 0; qp < npts; ++qp)
            for (size_type i = 0; i < iend; ++i)
              for (size_type j = 0; j < jend; ++j)
                tmp += field(qp, i, j) * data(qp, 0, j);


        result() = value_type(_sumInto)*result() + tmp;
      }
    };
    const size_type loopSize = inputFields.dimension(0)*inputFields.dimension(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(outputFields, inputData, inputFields, sumInto) );
  }



  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::Internal::
  contractDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          outputData,
                    const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                    const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                    const bool sumInto ) {
    typedef inputDataLeftValueType value_type;

    struct Functor {
      Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          _outputData;
      Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  _inputDataLeft;
      Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> _inputDataRight;
      const bool _sumInto; 

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          outputData_,
              Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft_,
              Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight_,
              const bool sumInto_) 
        : _outputData(outputData_), _inputDataLeft(inputDataLeft_), _inputDataRight(inputDataRight_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        const size_type cl = iter;
        
        auto result = Kokkos::subdynrankview( _outputData, cl );
        const auto left  = Kokkos::subdynrankview( _inputDataLeft,  cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto right = Kokkos::subdynrankview( _inputDataRight, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        size_type npts = left.dimension(0);
        size_type iend = left.dimension(1);
        size_type jend = left.dimension(2);

        value_type tmp(0);        
        for (size_type qp = 0; qp < npts; ++qp) 
          for (size_type i = 0; i < iend; ++i) 
            for (size_type j = 0; j < jend; ++j) 
              tmp += left(qp, i, j)*right(qp, i, j);
        result() = value_type(_sumInto)*result() + tmp;
      }
    };
    const size_type loopSize = inputDataLeft.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(outputData, inputDataLeft, inputDataRight, sumInto) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldScalar( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    { 
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.rank() != 3, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( rightFields.rank() != 3, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 3, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(0) != rightFields.dimension(0), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(2) != rightFields.dimension(2), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != rightFields.dimension(0), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != leftFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(2) != rightFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldVector( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.rank() != 4, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( rightFields.rank() != 4, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 3, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(0) != rightFields.dimension(0), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(2) != rightFields.dimension(2), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(3) != rightFields.dimension(3), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != rightFields.dimension(0), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != leftFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(2) != rightFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldTensor( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.rank()  != 5, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( rightFields.rank() != 5, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 3, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(0) != rightFields.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(2) != rightFields.dimension(2), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(3) != rightFields.dimension(3), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftFields.dimension(4) != rightFields.dimension(4), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != rightFields.dimension(0), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != leftFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(2) != rightFields.dimension(1), dbgInfo,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif

    ArrayTools<ExecSpaceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldScalar( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>  outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>    inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>   inputFields,
                           const bool sumInto ) {

#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.rank()  != 3, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.rank() != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(0) != inputData.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");

      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(1) != inputFields.dimension(2) && 
                                      inputData.dimension(1) != 1, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != inputFields.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");      
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != inputFields.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
      
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  } 

  
  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldVector( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.rank()  != 4, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.rank() != 3, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(0) != inputData.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(1) != inputFields.dimension(2) &&
                                      inputData.dimension(1) != 1, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(3) != inputData.dimension(2), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != inputFields.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != inputFields.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif

    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }
  


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldTensor( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.rank()  != 5, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.rank() != 4, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.rank() != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(0) != inputData.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputData.dimension(1) != inputFields.dimension(2) && inputData.dimension(1) != 1, dbgInfo,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(3) != inputData.dimension(2), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputFields.dimension(4) != inputData.dimension(3), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(0) != inputFields.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputFields.dimension(1) != inputFields.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif

    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }


  
  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataScalar( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.rank()  != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.rank() != 2, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.rank() != 1, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif

    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  } 

  
  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataVector( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.rank()  != 3, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.rank() != 3, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.rank() != 1, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(2) != inputDataRight.dimension(2), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }

  
  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataTensor( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.rank()  != 4, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataRight.rank() != 4, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.rank() != 1, dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(1) != inputDataRight.dimension(1), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(2) != inputDataRight.dimension(2), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inputDataLeft.dimension(3) != inputDataRight.dimension(3), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputData.dimension(0) != inputDataRight.dimension(0), dbgInfo, 
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }
  
}
#endif
