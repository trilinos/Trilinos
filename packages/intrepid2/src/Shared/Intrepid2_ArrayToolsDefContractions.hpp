// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
                           _outputFields.extent(0),
                           _outputFields.extent(1),
                           _outputFields.extent(2),
                           iter );

        const size_type npts = _leftFields.extent(2);
        const ordinal_type iend = _leftFields.extent(3);
        const ordinal_type jend = _leftFields.extent(4);

        _outputFields( cl, lbf, rbf ) *= (_sumInto ? 1 : 0); 
        for (size_type qp = 0; qp < npts; ++qp) 
          for (ordinal_type i = 0; i < iend; ++i) 
            for (ordinal_type j = 0; j < jend; ++j) 
              _outputFields( cl, lbf, rbf ) += _leftFields(cl, lbf, qp, i, j)*_rightFields(cl, rbf, qp, i, j);
      }
    };
    } //end namespace

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<DeviceType>::Internal::
  contractFieldField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                      const Kokkos::DynRankView<leftFieldValueType,  leftFieldProperties...>   leftFields,
                      const Kokkos::DynRankView<rightFieldValueType, rightFieldProperties...>  rightFields,
                      const bool sumInto ) {

    typedef Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outFieldViewType;
    typedef Kokkos::DynRankView<leftFieldValueType,leftFieldProperties...> leftFieldViewType;
    typedef Kokkos::DynRankView<rightFieldValueType,rightFieldProperties...> rightFieldViewType;
    typedef FunctorArrayTools::F_contractFieldField<outFieldViewType, leftFieldViewType, rightFieldViewType> FunctorType;

    const size_type loopSize = leftFields.extent(0)*leftFields.extent(1)*rightFields.extent(1);
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

      KOKKOS_DEFAULTED_FUNCTION
      ~F_contractDataField() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cl, bf;
        unrollIndex( cl, bf, 
                           _inputFields.extent(0),
                           _inputFields.extent(1),
                           iter );
        
        auto result = Kokkos::subview( _outputFields, cl, bf );

        const auto field = Kokkos::subview( _inputFields, cl, bf, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto data  = Kokkos::subview( _inputData,   cl,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        const size_type npts = field.extent(0);
        const ordinal_type iend = field.extent(1);
        const ordinal_type jend = field.extent(2);

        result() *= (_sumInto ? 1 : 0); 

        if(_inputData.extent(1) != 1)
          for (size_type qp = 0; qp < npts; ++qp)
            for (ordinal_type i = 0; i < iend; ++i)
              for (ordinal_type j = 0; j < jend; ++j)
                result() += field(qp, i, j) * data(qp, i, j);
        else
          for (size_type qp = 0; qp < npts; ++qp)
            for (ordinal_type i = 0; i < iend; ++i)
              for (ordinal_type j = 0; j < jend; ++j)
                result() += field(qp, i, j) * data(0, i, j);
      }
    };
    } //namespace

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::Internal::
  contractDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>  outputFields,
                     const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>    inputData,
                     const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>   inputFields,
                     const bool sumInto ) {

    typedef Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...>                 outputFieldsViewType;
    typedef Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>                   inputDataViewType;
    typedef Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>                  inputFieldsViewType;
    typedef FunctorArrayTools::F_contractDataField<outputFieldsViewType, inputDataViewType, inputFieldsViewType>  FunctorType;

    const size_type loopSize = inputFields.extent(0)*inputFields.extent(1);
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

      KOKKOS_DEFAULTED_FUNCTION
      ~F_contractDataData() = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        const size_type cl = iter;
        
        auto result = Kokkos::subview( _outputData, cl );
        const auto left  = Kokkos::subview( _inputDataLeft,  cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );
        const auto right = Kokkos::subview( _inputDataRight, cl, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() );

        size_type npts = left.extent(0);
        ordinal_type iend = left.extent(1);
        ordinal_type jend = left.extent(2);

        result() *= (_sumInto ? 1 : 0); 
        for (size_type qp = 0; qp < npts; ++qp) 
          for (ordinal_type i = 0; i < iend; ++i) 
            for (ordinal_type j = 0; j < jend; ++j) 
              result() += left(qp, i, j)*right(qp, i, j);
      }
    };
    } //namespace

  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::Internal::
  contractDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>          outputData,
                    const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                    const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                    const bool sumInto ) {
    typedef Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>       outputDataViewType;
    typedef Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>    inputDataLeftViewType;
    typedef Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...>   inputDataRightViewType;
    typedef FunctorArrayTools::F_contractDataData<outputDataViewType, inputDataLeftViewType, inputDataRightViewType> FunctorType;

    const size_type loopSize = inputDataLeft.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight, sumInto) );
  }



  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(2) != rightFields.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != leftFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(2) != rightFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldScalar): Second dimension of output container and first dimension of right input container must agree!");

    }
#endif

    ArrayTools<DeviceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 


  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(2) != rightFields.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(3) != rightFields.extent(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Third dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != leftFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(2) != rightFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldVector): Second dimension of output container and first dimension of right input container must agree!");
    }
#endif

    ArrayTools<DeviceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  } 


  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename leftFieldValueType,   class ...leftFieldProperties,
           typename rightFieldValueType,  class ...rightFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(2) != rightFields.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(3) != rightFields.extent(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Third dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( leftFields.extent(4) != rightFields.extent(4), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Fourth dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != rightFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != leftFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): First dimension of output container and first dimension of left input container must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(2) != rightFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractFieldFieldTensor): Second dimension of output container and first dimension of right input container must agree!");
    }
#endif

    ArrayTools<DeviceType>::Internal::contractFieldField( outputFields,
                                                             leftFields,
                                                             rightFields,
                                                             sumInto );
  }


  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");

      INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2) &&
                                      inputData.extent(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Second dimension of fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");      
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != inputFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldScalar): First dimensions (number of fields) of the fields input and output containers must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  } 

  
  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2) &&
                                      inputData.extent(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(3) != inputData.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Third dimension of the fields input container and second dimension of data input container (vector index) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != inputFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldVector): First dimensions of output container and fields input container (number of fields) must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }
  


  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2) && inputData.extent(1) != 1, std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree or first data dimension must be 1!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(3) != inputData.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Third dimension of the fields input container and second dimension of data input container (first tensor dimension) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(4) != inputData.extent(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Fourth dimension of the fields input container and third dimension of data input container (second tensor dimension) must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(0) != inputFields.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): Zeroth dimensions (numbers of integration domains) of the fields input and output containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputFields.extent(1) != inputFields.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataFieldTensor): First dimensions (number of fields) of output container and fields input container must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataField( outputFields,
                                                            inputData,
                                                            inputFields,
                                                            sumInto );
  }


  
  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataScalar): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  } 

  
  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) != inputDataRight.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Second dimensions (numbers of vector components) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataVector): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }

  
  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
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
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (number of integration domains) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): First dimensions (numbers of integration points) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(2) != inputDataRight.extent(2), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Second dimensions (first tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(3) != inputDataRight.extent(3), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Third dimensions (second tensor dimensions) of the left and right input containers must agree!");
      INTREPID2_TEST_FOR_EXCEPTION( outputData.extent(0) != inputDataRight.extent(0), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::contractDataDataTensor): Zeroth dimensions (numbers of integration domains) of the input and output containers must agree!");
    }
#endif
    
    ArrayTools<DeviceType>::Internal::contractDataData( outputData,
                                                           inputDataLeft,
                                                           inputDataRight,
                                                           sumInto );
  }
  
}
#endif
