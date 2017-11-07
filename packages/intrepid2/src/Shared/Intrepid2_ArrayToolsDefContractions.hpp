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


    namespace FunctorArrayTools {
    /**
      \brief Functor to contractFieldField see Intrepid2::ArrayTools for more
    */ 
    template < typename outFieldViewType , typename leftFieldViewType , typename rightFieldViewType >
    struct F_contractFieldField{
      outFieldViewType     _outputFields;
      leftFieldViewType    _leftFields;
      rightFieldViewType   _rightFields;
      const bool _sumInto; 
      typedef typename outFieldViewType::value_type value_type;

      KOKKOS_INLINE_FUNCTION
      F_contractFieldField(outFieldViewType outputFields_,
              leftFieldViewType leftFields_,
              rightFieldViewType rightFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _leftFields(leftFields_), _rightFields(rightFields_), _sumInto(sumInto_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, lbf, rbf;
        unrollIndex( cl, lbf, rbf, 
                           _outputFields.dimension(0),
                           _outputFields.dimension(1),
                           _outputFields.dimension(2),
                           iter );

        auto result = Kokkos::subview( _outputFields,  cl, lbf, rbf );

        const auto left  = Kokkos::subview( _leftFields,  cl, lbf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto right = Kokkos::subview( _rightFields, cl, rbf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        const size_type npts = left.dimension(0);
        const ordinal_type iend = left.dimension(1);
        const ordinal_type jend = left.dimension(2);

        value_type tmp(0);        
        for (size_type qp = 0; qp < npts; ++qp) 
          for (ordinal_type i = 0; i < iend; ++i) 
            for (ordinal_type j = 0; j < jend; ++j) 
              tmp += left(qp, i, j)*right(qp, i, j);
        if (_sumInto)
          result() = result() + tmp;
        else
          result() = tmp;
      }
    };
    } //end namespace

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<SpT>::Internal::
  contractFieldField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                      const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                      const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                      const bool sumInto ) {

    typedef Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outFieldViewType;
    typedef Kokkos::DynRankView<leftFieldValueType,leftFieldProperties...> leftFieldViewType;
    typedef Kokkos::DynRankView<rightFieldValueType,rightFieldProperties...> rightFieldViewType;
    typedef FunctorArrayTools::F_contractFieldField<outFieldViewType, leftFieldViewType, rightFieldViewType> FunctorType;
    typedef typename ExecSpace< typename outFieldViewType::execution_space, SpT >::ExecSpaceType ExecSpaceType;


    const size_type loopSize = leftFields.dimension(0)*leftFields.dimension(1)*rightFields.dimension(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(outputFields, leftFields, rightFields, sumInto) );
  }


    namespace FunctorArrayTools {
    /**
      \brief Functor to contractDataField see Intrepid2::ArrayTools for more
    */ 
    template < typename outputFieldsViewType , typename inputDataViewType , typename inputFieldsViewType >
    struct F_contractDataField {
      outputFieldsViewType  _outputFields;
      inputDataViewType     _inputData;
      inputFieldsViewType   _inputFields;
      const bool _sumInto; 
      typedef typename outputFieldsViewType::value_type value_type;

      KOKKOS_INLINE_FUNCTION
      F_contractDataField(outputFieldsViewType outputFields_,
              inputDataViewType inputData_,
              inputFieldsViewType inputFields_,
              const bool sumInto_) 
        : _outputFields(outputFields_), _inputData(inputData_), _inputFields(inputFields_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~F_contractDataField() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, bf;
        unrollIndex( cl, bf, 
                           _inputFields.dimension(0),
                           _inputFields.dimension(1),
                           iter );
        
        auto result = Kokkos::subview( _outputFields, cl, bf );

        const auto field = Kokkos::subview( _inputFields, cl, bf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto data  = Kokkos::subview( _inputData,   cl,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        const size_type npts = field.dimension(0);
        const ordinal_type iend = field.dimension(1);
        const ordinal_type jend = field.dimension(2);

        value_type tmp(0);        

        if(_inputData.dimension(1) != 1)
          for (size_type qp = 0; qp < npts; ++qp)
            for (ordinal_type i = 0; i < iend; ++i)
              for (ordinal_type j = 0; j < jend; ++j)
                tmp += field(qp, i, j) * data(qp, i, j);
        else
          for (size_type qp = 0; qp < npts; ++qp)
            for (ordinal_type i = 0; i < iend; ++i)
              for (ordinal_type j = 0; j < jend; ++j)
                tmp += field(qp, i, j) * data(0, i, j);

        if (_sumInto)
          result() = result() + tmp;
        else
          result() = tmp;
      }
    };
    } //namespace

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<SpT>::Internal::
  contractDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>  outputFields,
                     const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>    inputData,
                     const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>   inputFields,
                     const bool sumInto ) {

    typedef Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>                 outputFieldsViewType;
    typedef Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>                   inputDataViewType;
    typedef Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>                  inputFieldsViewType;
    typedef FunctorArrayTools::F_contractDataField<outputFieldsViewType, inputDataViewType, inputFieldsViewType>  FunctorType;
    typedef typename ExecSpace< typename inputFieldsViewType::execution_space , SpT >::ExecSpaceType               ExecSpaceType;

    const size_type loopSize = inputFields.dimension(0)*inputFields.dimension(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(outputFields, inputData, inputFields, sumInto) );
  }


    namespace FunctorArrayTools {
    /**
      \brief Functor to contractDataData see Intrepid2::ArrayTools for more
    */ 
    template < typename outputDataViewType , typename inputDataLeftViewType , typename inputDataRightViewType >
    struct F_contractDataData {
      outputDataViewType  _outputData;
      inputDataLeftViewType   _inputDataLeft;
      inputDataRightViewType  _inputDataRight;
      const bool _sumInto; 
      typedef typename outputDataViewType::value_type value_type;

      KOKKOS_INLINE_FUNCTION
      F_contractDataData(outputDataViewType outputData_,
              inputDataLeftViewType inputDataLeft_,
              inputDataRightViewType inputDataRight_,
              const bool sumInto_) 
        : _outputData(outputData_), _inputDataLeft(inputDataLeft_), _inputDataRight(inputDataRight_), _sumInto(sumInto_) {}

      KOKKOS_INLINE_FUNCTION
      ~F_contractDataData() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        const size_type cl = iter;
        
        auto result = Kokkos::subview( _outputData, cl );
        const auto left  = Kokkos::subview( _inputDataLeft,  cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto right = Kokkos::subview( _inputDataRight, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        size_type npts = left.dimension(0);
        ordinal_type iend = left.dimension(1);
        ordinal_type jend = left.dimension(2);

        value_type tmp(0);        
        for (size_type qp = 0; qp < npts; ++qp) 
          for (ordinal_type i = 0; i < iend; ++i) 
            for (ordinal_type j = 0; j < jend; ++j) 
              tmp += left(qp, i, j)*right(qp, i, j);

        if (_sumInto)
          result() = result() + tmp;
        else
          result() = tmp;
      }
    };
    } //namespace

  template<typename SpT>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<SpT>::Internal::
  contractDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          outputData,
                    const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                    const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                    const bool sumInto ) {
    typedef Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>       outputDataViewType;
    typedef Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>    inputDataLeftViewType;
    typedef Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...>   inputDataRightViewType;
    typedef FunctorArrayTools::F_contractDataData<outputDataViewType, inputDataLeftViewType, inputDataRightViewType> FunctorType;
    typedef typename ExecSpace< typename inputDataLeftViewType::execution_space , SpT>::ExecSpaceType  ExecSpaceType;

    const size_type loopSize = inputDataLeft.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight, sumInto) );
  }



  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldScalar(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID2_DEBUG
    { 
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of the left input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( rightFields.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of right input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(2) != rightFields.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != leftFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(2) != rightFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");

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
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldVector(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.rank() != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of the left input argument must equal 4!");
      INTREPID2_TEST_FOR_EXCEPTION( rightFields.rank() != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of right input argument must equal 4!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(2) != rightFields.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(3) != rightFields.dimension(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != leftFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(2) != rightFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractFieldFieldTensor(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                            const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                            const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                            const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.rank()  != 5, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of the left input argument must equal 5!");
      INTREPID2_TEST_FOR_EXCEPTION( rightFields.rank() != 5, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of right input argument must equal 5!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Rank of output argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(2) != rightFields.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(3) != rightFields.dimension(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.dimension(4) != rightFields.dimension(4), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != rightFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != leftFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(2) != rightFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldScalar(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>  outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>    inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>   inputFields,
                           const bool sumInto ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank()  != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the fields input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of the data input argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");

      INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(2) &&
                                      inputData.dimension(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != inputFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");      
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != inputFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldVector(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank()  != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the fields input argument must equal 4!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of the data input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(2) &&
                                      inputData.dimension(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(3) != inputData.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != inputFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != inputFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataFieldTensor(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank()  != 5, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the fields input argument must equal 5!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of the data input argument must equal 4!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Rank of output argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(0) != inputData.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.dimension(1) != inputFields.dimension(2) && inputData.dimension(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(3) != inputData.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.dimension(4) != inputData.dimension(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(0) != inputFields.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.dimension(1) != inputFields.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataScalar(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank()  != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of the left input argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() != 2, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of right input argument must equal 2!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataVector( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank()  != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of the left input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() != 3, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of right input argument must equal 3!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) != inputDataRight.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
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
  void
  ArrayTools<ExecSpaceType>::
  contractDataDataTensor(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank()  != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of the left input argument must equal 4");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() != 4, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of right input argument must equal 4!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Rank of output argument must equal 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(1) != inputDataRight.dimension(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(2) != inputDataRight.dimension(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.dimension(3) != inputDataRight.dimension(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.dimension(0) != inputDataRight.dimension(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    }
#endif
    
    ArrayTools<ExecSpaceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }
  
}
#endif
