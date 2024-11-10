// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayToolsDefScalar.hpp
    \brief  Definition file for scalar multiply operations of the array tools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_SCALAR_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_SCALAR_HPP__

namespace Intrepid2 {

  namespace FunctorArrayTools {
    /**
       \brief Functor for scalarMultiply see Intrepid2::ArrayTools for more
    */
    template<typename OutputViewType,
             typename inputLeftViewType,
             typename inputRightViewType,
             bool equalRank,
             bool reciprocal>
    struct F_scalarMultiply {
            OutputViewType _output;
      const inputLeftViewType _inputLeft;
      const inputRightViewType _inputRight;

      KOKKOS_INLINE_FUNCTION
      F_scalarMultiply(OutputViewType output_,
                       inputLeftViewType inputLeft_,
                       inputRightViewType inputRight_)
        : _output(output_),
          _inputLeft(inputLeft_),
          _inputRight(inputRight_) {}
      
      // DataData
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type pt) const {
        const auto val = _inputLeft(cl , pt%_inputLeft.extent(1));
        
        //const ordinal_type cp[2] = { cl, pt };
        //ViewAdapter<2,OutputViewType> out(cp, _output);
        auto out = Kokkos::subview(_output, cl, pt, Kokkos::ALL(), Kokkos::ALL());
        if (equalRank) {
          //const ViewAdapter<2,inputRightViewType> right(cp, _inputRight);
          const auto right = Kokkos::subview(_inputRight, cl, pt, Kokkos::ALL(), Kokkos::ALL());
          if (reciprocal) Kernels::inv_scalar_mult_mat(out, val, right);
          else            Kernels::    scalar_mult_mat(out, val, right);
        } else {
          //const ordinal_type p[1] = { pt };
          //const ViewAdapter<1,inputRightViewType> right(p, _inputRight);
          const auto right = Kokkos::subview(_inputRight, pt, Kokkos::ALL(), Kokkos::ALL());
          if (reciprocal) Kernels::inv_scalar_mult_mat(out, val, right);
          else            Kernels::    scalar_mult_mat(out, val, right);
        } 
      }
      
      // DataField
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type bf,
                      const ordinal_type pt) const {
        const auto val = _inputLeft(cl , pt%_inputLeft.extent(1));

        //const ordinal_type cfp[3] = { cl, bf, pt };
        //ViewAdapter<3,OutputViewType> out(cfp, _output);
        auto out = Kokkos::subview(_output, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL());
        if (equalRank) {
          //const ViewAdapter<3,inputRightViewType> right(cfp, _inputRight);          
          auto right = Kokkos::subview(_inputRight, cl, bf, pt, Kokkos::ALL(), Kokkos::ALL());
          if (reciprocal) Kernels::inv_scalar_mult_mat(out, val, right);
          else            Kernels::    scalar_mult_mat(out, val, right);          
        } else {
          //const ordinal_type fp[2] = { bf, pt };
          //const ViewAdapter<2,inputRightViewType> right(fp, _inputRight);          
          auto right = Kokkos::subview(_inputRight, bf, pt, Kokkos::ALL(), Kokkos::ALL());
          if (reciprocal) Kernels::inv_scalar_mult_mat(out, val, right);
          else            Kernels::    scalar_mult_mat(out, val, right);
        } 
      }
    };
  } 

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  ArrayTools<DeviceType>::
  scalarMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool reciprocal ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputData.rank() != 2, std::invalid_argument,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input data container must have rank 2.");

      if (outputFields.rank() <= inputFields.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 3, 4, or 5.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != inputFields.rank(), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input and output fields containers must have the same rank.");
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(0) != inputData.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions (number of integration domains) of the fields and data input containers must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(2) &&
                                  inputData.extent(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Second dimension of the fields input container and first dimension of data input container (number of integration points) must agree, or first data dimension must be 1!");
        for (size_type i=0;i<inputFields.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(i) != outputFields.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i) of outputFields");
        }
      }
      else {
        INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 2 || inputFields.rank() > 4, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Input fields container must have rank 2, 3, or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( outputFields.rank() != (inputFields.rank()+1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): The rank of the input fields container must be one less than the rank of the output fields container.");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(1) != inputFields.extent(1) &&
                                  inputData.extent(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): First dimensions of fields input container and data input container (number of integration points) must agree or first data dimension must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputData.extent(0) != outputFields.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataField): Zeroth dimensions of fields output container and data input containers (number of integration domains) must agree!");
        for (size_type i=0;i<inputFields.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputFields.extent(i) != outputFields.extent(i+1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataField): inputFields dimension (i) does not match to the dimension (i+1) of outputFields");
        }
      }
    }
#endif

    typedef Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFieldViewType;
    typedef Kokkos::DynRankView<inputDataValueType,inputDataProperties...> inputDataViewType;
    typedef Kokkos::DynRankView<inputFieldValueType,inputFieldProperties...> inputFieldViewType;
    
    using range_policy_type = Kokkos::MDRangePolicy
        < ExecSpaceType, Kokkos::Rank<3>, Kokkos::IndexType<ordinal_type> >;

    const range_policy_type policy( { 0, 0, 0 },
                                    { /*C*/ outputFields.extent(0), /*F*/ outputFields.extent(1), /*P*/ outputFields.extent(2) } );
    
    const bool equalRank = ( outputFields.rank() == inputFields.rank() );
    if (equalRank) 
      if (reciprocal) {
        typedef FunctorArrayTools::F_scalarMultiply<outputFieldViewType,inputDataViewType,inputFieldViewType,true,true> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputFields, inputData, inputFields) );
      } else {
        typedef FunctorArrayTools::F_scalarMultiply<outputFieldViewType,inputDataViewType,inputFieldViewType,true,false> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputFields, inputData, inputFields) );
      }
    else 
      if (reciprocal) {
        typedef FunctorArrayTools::F_scalarMultiply<outputFieldViewType,inputDataViewType,inputFieldViewType,false,true> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputFields, inputData, inputFields) );
      } else {
        typedef FunctorArrayTools::F_scalarMultiply<outputFieldViewType,inputDataViewType,inputFieldViewType,false,false> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputFields, inputData, inputFields) );
      }
  }


  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  ArrayTools<DeviceType>::
  scalarMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool reciprocal ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.rank() != 2, std::invalid_argument,
                                ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");

      if (outputData.rank() <= inputDataRight.rank()) {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 2 || inputDataRight.rank() > 4, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2, 3, or 4.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != inputDataRight.rank(), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(0) != inputDataLeft.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(1) &&
                                  inputDataLeft.extent(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree, or first dimension of the left data input container must be 1!");
        for (size_type i=0;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(i) != outputData.extent(i), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i) of outputData");
        }
      } else {
        INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.rank() < 1 || inputDataRight.rank() > 3, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1, 2, or 3.");
        INTREPID2_TEST_FOR_EXCEPTION( outputData.rank() != (inputDataRight.rank()+1), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(1) != inputDataRight.extent(0)  &&
                                  inputDataLeft.extent(1) != 1, std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
        INTREPID2_TEST_FOR_EXCEPTION( inputDataLeft.extent(0) != outputData.extent(0), std::invalid_argument,
                                  ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
        for (size_type i=0;i<inputDataRight.rank();++i) {
          INTREPID2_TEST_FOR_EXCEPTION( inputDataRight.extent(i) != outputData.extent(i+1), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::scalarMultiplyDataData): inputDataRight dimension (i) does not match to the dimension (i+1) of outputData");
        }
      }
    }
#endif

    typedef Kokkos::DynRankView<outputDataValueType,outputDataProperties...> outputDataViewType;
    typedef Kokkos::DynRankView<inputDataLeftValueType,inputDataLeftProperties...> inputDataLeftViewType;
    typedef Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRightViewType;

    using range_policy_type = Kokkos::MDRangePolicy
      < ExecSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;

    const range_policy_type policy( { 0, 0 },
                                    { /*C*/ outputData.extent(0), /*P*/ outputData.extent(1) } );

    const bool equalRank = ( outputData.rank() == inputDataRight.rank() );    
    if (equalRank) 
      if (reciprocal) {
        typedef FunctorArrayTools::F_scalarMultiply<outputDataViewType,inputDataLeftViewType,inputDataRightViewType,true,true> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight) );
      } else {
        typedef FunctorArrayTools::F_scalarMultiply<outputDataViewType,inputDataLeftViewType,inputDataRightViewType,true,false> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight) );
      }
    else 
      if (reciprocal) {
        typedef FunctorArrayTools::F_scalarMultiply<outputDataViewType,inputDataLeftViewType,inputDataRightViewType,false,true> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight) );
      } else {
        typedef FunctorArrayTools::F_scalarMultiply<outputDataViewType,inputDataLeftViewType,inputDataRightViewType,false,false> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputData, inputDataLeft, inputDataRight) );
      }
  }
  
} // end namespace Intrepid2

#endif
