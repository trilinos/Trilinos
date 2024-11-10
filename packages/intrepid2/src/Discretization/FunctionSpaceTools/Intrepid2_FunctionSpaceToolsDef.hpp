// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_FunctionSpaceToolsDef.hpp
    \brief  Definition file for the Intrepid2::FunctionSpaceTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_FUNCTIONSPACETOOLS_DEF_HPP__
#define __INTREPID2_FUNCTIONSPACETOOLS_DEF_HPP__

#include "Intrepid2_FunctorIterator.hpp"
#include "Intrepid2_TensorArgumentIterator.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace Intrepid2 {
  // ------------------------------------------------------------------------------------
  template<typename DeviceType>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HGRADtransformVALUE(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
                       const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
    if(output.rank() == input.rank()) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      for (size_type i=0;i< input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( (input.extent(i) != output.extent(i)), std::invalid_argument,
                                        ">>> ERROR (FunctionSpaceTools::HGRADtransformVALUE): Dimensions of input and output fields containers do not match.");
      }
    }
#endif
      RealSpaceTools<DeviceType>::clone(output, input);
    }
    else
      ArrayTools<DeviceType>::cloneFields(output, input);
  }

  template<typename DeviceType>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHGradDataFromPhysToRef(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
                       const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
    if(output.rank() == input.rank()) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      for (size_type i=0;i< input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( (input.extent(i) != output.extent(i)), std::invalid_argument,
                                        ">>> ERROR (FunctionSpaceTools::mapHGradDataFromPhysToRef): Dimensions of input and output fields containers do not match.");
      }
    }
#endif
      RealSpaceTools<DeviceType>::clone(output, input);
    }
    else
      ArrayTools<DeviceType>::cloneFields(output, input);
  }

  template<typename DeviceType>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHGradDataFromPhysSideToRefSide(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
                       const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
    if(output.rank() == input.rank()) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      for (size_type i=0;i< input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( (input.extent(i) != output.extent(i)), std::invalid_argument,
                                        ">>> ERROR (FunctionSpaceTools::mapHGradDataSideFromPhysToRefSide): Dimensions of input and output fields containers do not match.");
      }
    }
#endif
      RealSpaceTools<DeviceType>::clone(output, input);
    }
    else
      ArrayTools<DeviceType>::cloneFields(output, input);
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {
   /**
       \brief Functor for calculation HGRADtransformGRAD, see Intrepid2::FunctionSpaceTools for more
   */
  }
  
  template<typename DeviceType>
  template<typename OutputValViewType,
           typename JacobianInverseViewType,
           typename InputValViewType>
  void
  FunctionSpaceTools<DeviceType>::
  HGRADtransformGRAD(       OutputValViewType       outputVals,
                      const JacobianInverseViewType jacobianInverse,
                      const InputValViewType        inputVals ) {
    return HCURLtransformVALUE(outputVals, jacobianInverse, inputVals);
  }
  
  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HCURLtransformVALUE(       Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                       const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                       const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals ) {
    ArrayTools<DeviceType>::matvecProductDataField(outputVals, jacobianInverse, inputVals, 'T');
  }

  template<typename DeviceType>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianValueType,        class ...jacobianProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHCurlDataFromPhysToRef(       Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                              const Kokkos::DynRankView<jacobianValueType,       jacobianProperties...>        jacobian,
                              const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals ) {
    ArrayTools<DeviceType>::matvecProductDataData(outputVals, jacobian, inputVals, 'T');
  }


  namespace FunctorFunctionSpaceTools {

  template<typename outViewType,
  typename inputViewType,
  typename metricViewType
  >
  struct F_negativeWeighted2dInputCrossK {
    outViewType output_;
    const inputViewType input_;
    const metricViewType metricTensorDet_;

    KOKKOS_INLINE_FUNCTION
    F_negativeWeighted2dInputCrossK( outViewType output,
        const inputViewType input,
        const metricViewType metricTensorDet)
    : output_(output), input_(input), metricTensorDet_(metricTensorDet){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type ic) const {
      for (size_t pt=0; pt < input_.extent(1); pt++) {
        auto measure = std::sqrt(metricTensorDet_(ic,pt));
        output_(ic,pt,0) = - measure*input_(ic,pt,1);
        output_(ic,pt,1) = measure*input_(ic,pt,0);
      }
    }
  };
  }

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename tangentsValueType,    class ...tangentsProperties,
           typename metricTensorInvValueType,    class ...metricTensorInvProperties,
           typename metricTensorDetValueType, class ...metricTensorDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHCurlDataCrossNormalFromPhysSideToRefSide(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<tangentsValueType,   tangentsProperties...>    tangents,
                      const Kokkos::DynRankView<metricTensorInvValueType,metricTensorInvProperties...> metricTensorInv,
                      const Kokkos::DynRankView<metricTensorDetValueType,metricTensorDetProperties...> metricTensorDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    auto work = Kokkos::DynRankView<outputValValueType,  outputValProperties...>("work", outputVals.extent(0), outputVals.extent(1),outputVals.extent(2));
    ArrayTools<DeviceType>::matvecProductDataData(outputVals, tangents, inputVals, 'T');
    typename DeviceType::execution_space().fence();
    ArrayTools<DeviceType>::matvecProductDataData(work, metricTensorInv, outputVals);
    typename DeviceType::execution_space().fence();
    using FunctorType = FunctorFunctionSpaceTools::F_negativeWeighted2dInputCrossK<decltype(outputVals),decltype(work),decltype(metricTensorDet)>;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,outputVals.extent(0)), FunctorType(outputVals, work, metricTensorDet) );
  }


  namespace FunctorFunctionSpaceTools {

  template<typename outViewType,
  typename inputViewType,
  typename metricViewType
  >
  struct F_weighedInput {
    outViewType output_;
    const inputViewType input_;
    const metricViewType metricTensorDet_;
    const double scaling_;

    KOKKOS_INLINE_FUNCTION
    F_weighedInput( outViewType output,
        const inputViewType input,
        const metricViewType metricTensorDet,
        const double scaling)
    : output_(output), input_(input), metricTensorDet_(metricTensorDet), scaling_(scaling){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type ic) const {
      for (size_t pt=0; pt < input_.extent(1); pt++) {
        auto measure = std::sqrt(metricTensorDet_(ic,pt));
        output_.access(ic,pt) = scaling_ * measure * input_.access(ic,pt);
      }
    }
  };
  }

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHCurlDataCrossNormalFromPhysSideToRefSide(
            Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> metricTensorDet,
      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
      using FunctorType = FunctorFunctionSpaceTools::F_weighedInput<decltype(outputVals),decltype(inputVals),decltype(metricTensorDet)>;
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,outputVals.extent(0)), FunctorType(outputVals, inputVals, metricTensorDet, -1.0) );
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HCURLtransformCURL(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    if(jacobian.data()==NULL || jacobian.extent(2)==2) //2D case
      return HVOLtransformVALUE(outputVals, jacobianDet, inputVals);
    else
      return HDIVtransformVALUE(outputVals, jacobian, jacobianDet, inputVals);
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HCURLtransformCURL(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.rank() == 4, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HCURLtransformCURL): Output rank must have rank 3.\n If these are 3D fields, then use the appropriate overload of this function.");
    }
#endif
    return HVOLtransformVALUE(outputVals, jacobianDet, inputVals);
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HGRADtransformCURL(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(3)!=2, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformCURL):\n output field is 3D by the function is meant for 2D fields");
    }
#endif
    return HDIVtransformVALUE(outputVals, jacobian, jacobianDet, inputVals);
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HDIVtransformVALUE(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<DeviceType>::matvecProductDataField(outputVals, jacobian, inputVals, 'N');
    ArrayTools<DeviceType>::scalarMultiplyDataField(outputVals, jacobianDet, outputVals, true);
  }

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianInverseValueType,    class ...jacobianInverseProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHDivDataFromPhysToRef(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianInverseValueType,   jacobianInverseProperties...>    jacobianInv,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<DeviceType>::matvecProductDataData(outputVals, jacobianInv, inputVals, 'N');
    typename DeviceType::execution_space().fence();
    ArrayTools<DeviceType>::scalarMultiplyDataData(outputVals, jacobianDet, outputVals, false);
  }



  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHDivDataDotNormalFromPhysSideToRefSide(
            Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> metricTensorDet,
      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    using FunctorType = FunctorFunctionSpaceTools::F_weighedInput<decltype(outputVals),decltype(inputVals),decltype(metricTensorDet)>;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,outputVals.extent(0)), FunctorType(outputVals, inputVals, metricTensorDet, 1.0) );
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HDIVtransformDIV(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                    const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                    const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    return HVOLtransformVALUE(outputVals, jacobianDet, inputVals);
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  HVOLtransformVALUE(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<DeviceType>::scalarMultiplyDataField(outputVals, jacobianDet, inputVals, true);
  }
  
  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  mapHVolDataFromPhysToRef(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                            const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                            const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<DeviceType>::scalarMultiplyDataData(outputVals, jacobianDet, inputVals, false);
  }



  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename leftValueValueType,   class ...leftValueProperties,
           typename rightValueValueType,  class ...rightValueProperties>
  void
  FunctionSpaceTools<DeviceType>::
  integrate(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<leftValueValueType,  leftValueProperties...>   leftValues,
             const Kokkos::DynRankView<rightValueValueType, rightValueProperties...>  rightValues,
             const bool sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( leftValues.rank() < 2 ||
                                    leftValues.rank() > 4, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::integrate): Left data must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_EXCEPTION( outputValues.rank() < 1 ||
                                    outputValues.rank() > 3, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::integrate): Output values must have rank 1, 2 or 3.");
    }    
#endif
    
    const ordinal_type outRank  = outputValues.rank();
    const ordinal_type leftRank = leftValues.rank();
    const ordinal_type mode = outRank*10 + leftRank;

    switch (mode) {
      // DataData
    case 12:
      ArrayTools<DeviceType>::contractDataDataScalar( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;
    case 13:
      ArrayTools<DeviceType>::contractDataDataVector( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;
    case 14:
      ArrayTools<DeviceType>::contractDataDataTensor( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;

      // DataField
    case 22:
      ArrayTools<DeviceType>::contractDataFieldScalar( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;
    case 23:
      ArrayTools<DeviceType>::contractDataFieldVector( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;
    case 24:
      ArrayTools<DeviceType>::contractDataFieldTensor( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;

      // FieldField
    case 33:
      ArrayTools<DeviceType>::contractFieldFieldScalar( outputValues,
                                                 leftValues,
                                                 rightValues,
                                                 sumInto );
      break;
    case 34:
      ArrayTools<DeviceType>::contractFieldFieldVector( outputValues,
                                                 leftValues,
                                                 rightValues,
                                                 sumInto );
      break;
    case 35:
      ArrayTools<DeviceType>::contractFieldFieldTensor( outputValues,
                                                 leftValues,
                                                 rightValues,
                                                 sumInto );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 1 || outRank > 3, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::integrate): outRank must be 1,2, or 3.");
      INTREPID2_TEST_FOR_EXCEPTION( leftRank < 2 || leftRank > 4, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::integrate): leftRank must be 1,2, 3 or 4.");
    }
    }
  }

  // ------------------------------------------------------------------------------------  
  namespace FunctorFunctionSpaceTools {
   /**
       \brief Functor for calculation of cell measure, see Intrepid2::FunctionSpaceTools for more
   */
    template<typename outputValViewType,
             typename inputDetViewType,
             typename inputWeightViewType>
    struct F_computeCellMeasure {
            outputValViewType   _outputVals;
      const inputDetViewType    _inputDet;
      const inputWeightViewType    _inputWeight;
      
      KOKKOS_INLINE_FUNCTION
      F_computeCellMeasure( outputValViewType outputVals_,
                            inputDetViewType inputDet_,
                            inputWeightViewType inputWeight_)
        : _outputVals(outputVals_), 
          _inputDet(inputDet_),
          _inputWeight(inputWeight_) {}

      typedef ordinal_type value_type;
      
//      KOKKOS_INLINE_FUNCTION
//      void init( value_type &dst ) const {
//        dst = false;
//      }
      
//      KOKKOS_INLINE_FUNCTION
//      void join( volatile value_type &dst,
//                 const volatile value_type &src ) const {
//       dst |= src;
//      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cl, value_type &dst) const {
        // negative jacobian check
        const bool hasNegativeDet = (_inputDet(cl, 0) < 0.0);
        dst |= hasNegativeDet;

        // make it absolute
        const auto sign = (hasNegativeDet ? -1.0 : 1.0);
        const ordinal_type pt_end = _outputVals.extent(1);
        for (ordinal_type pt=0;pt<pt_end;++pt) 
          _outputVals(cl, pt) = sign*_inputDet(cl, pt)*_inputWeight(pt);
      }
    };
  }

  template<typename DeviceType>
  template<typename OutputValViewType,
           typename InputDetViewType,
           typename InputWeightViewType>
  bool
  FunctionSpaceTools<DeviceType>::
  computeCellMeasure(       OutputValViewType   outputVals,
                      const InputDetViewType    inputDet,
                      const InputWeightViewType inputWeights ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( rank(inputDet)     != 2 || 
                                    rank(inputWeights) != 1 || 
                                    rank(outputVals)   != 2, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Ranks are not compatible.");
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(0) != inputDet.extent(0), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Cell dimension does not match.");
      INTREPID2_TEST_FOR_EXCEPTION( inputDet.extent(1)      !=  outputVals.extent(1) ||
                                    inputWeights.extent(0) !=  outputVals.extent(1), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Point dimension does not match.");
    }
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(outputVals)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inputDet)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inputWeights)::memory_space>::accessible;
    static_assert(are_accessible, "FunctionSpaceTools<DeviceType>::computeCellMeasure(..): input/output views' memory spaces are not compatible with DeviceType");

    using FunctorType = FunctorFunctionSpaceTools::F_computeCellMeasure
                     <decltype(outputVals),decltype(inputDet),decltype(inputWeights)>;
    
    const ordinal_type C = inputDet.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    
    typename FunctorType::value_type hasNegativeDet = false;
    Kokkos::parallel_reduce( policy, FunctorType(outputVals, inputDet, inputWeights), hasNegativeDet );
    
    return hasNegativeDet;
  } 

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightProperties,
           typename scratchValueType,     class ...scratchProperties>
  void 
  FunctionSpaceTools<DeviceType>::
  computeFaceMeasure(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<inputJacValueType,  inputJacProperties...>     inputJac,
                      const Kokkos::DynRankView<inputWeightValueType,inputWeightProperties...> inputWeights,
                      const ordinal_type                   whichFace,
                      const shards::CellTopology  parentCell,
                      const Kokkos::DynRankView<scratchValueType,    scratchProperties...>     scratch ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inputJac.rank() != 4, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Input Jacobian container must have rank 4.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Scratch view imust have rank 1.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.size() < inputJac.size(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Scratch storage must be greater than or equal to inputJac's one.");
#endif

    // face normals (reshape scratch)
    // Kokkos::DynRankView<scratchValueType,scratchProperties...> faceNormals(scratch.data(), 
    //                                                                        inputJac.extent(0), 
    //                                                                        inputJac.extent(1), 
    //                                                                        inputJac.extent(2));
    auto vcprop = Kokkos::common_view_alloc_prop(scratch);
    //typedef Kokkos::DynRankView<scratchValueType, typename decltype(scratch)::memory_space> viewType;
    typedef Kokkos::DynRankView<scratchValueType, DeviceType> viewType;
    viewType faceNormals(Kokkos::view_wrap(scratch.data(), vcprop),
                         inputJac.extent(0),
                         inputJac.extent(1),
                         inputJac.extent(2));

    // compute normals
    CellTools<DeviceType>::getPhysicalFaceNormals(faceNormals, inputJac, whichFace, parentCell);

    // compute lenghts of normals
    RealSpaceTools<DeviceType>::vectorNorm(outputVals, faceNormals, NORM_TWO);

    // multiply with weights
    ArrayTools<DeviceType>::scalarMultiplyDataData(outputVals, outputVals, inputWeights);
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightProperties,
           typename scratchValueType,     class ...scratchProperties>
  void 
  FunctionSpaceTools<DeviceType>::
  computeEdgeMeasure(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>  outputVals,
                      const Kokkos::DynRankView<inputJacValueType,  inputJacProperties...>   inputJac,
                      const Kokkos::DynRankView<inputWeightValueType,inputWeightProperties...> inputWeights,
                      const ordinal_type                   whichEdge,
                      const shards::CellTopology  parentCell,
                      const Kokkos::DynRankView<scratchValueType,    scratchProperties...>    scratch ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( (inputJac.rank() != 4), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Input Jacobian container must have rank 4.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Scratch view must have a rank 1.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.size() < inputJac.size(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Scratch storage must be greater than or equal to inputJac'one.");
#endif

    // edge tangents (reshape scratch)
    // Kokkos::DynRankView<scratchValueType,scratchProperties...> edgeTangents(scratch.data(), 
    //                                                                         inputJac.extent(0), 
    //                                                                         inputJac.extent(1), 
    //                                                                         inputJac.extent(2));
    auto vcprop = Kokkos::common_view_alloc_prop(scratch);
    //typedef Kokkos::DynRankView<scratchValueType, typename decltype(scratch)::memory_space> viewType;
    typedef Kokkos::DynRankView<scratchValueType, DeviceType> viewType;
    viewType edgeTangents(Kokkos::view_wrap(scratch.data(), vcprop),
                         inputJac.extent(0),
                         inputJac.extent(1),
                         inputJac.extent(2));

    // compute normals
    CellTools<DeviceType>::getPhysicalEdgeTangents(edgeTangents, inputJac, whichEdge, parentCell);

    // compute lenghts of tangents
    RealSpaceTools<DeviceType>::vectorNorm(outputVals, edgeTangents, NORM_TWO);

    // multiply with weights
    ArrayTools<DeviceType>::scalarMultiplyDataData(outputVals, outputVals, inputWeights);
  }

  // ------------------------------------------------------------------------------------  

  template<typename DeviceType>
  template<typename outputValValueType,    class ...outputValProperties,
           typename inputMeasureValueType, class ...inputMeasureProperties,
           typename inputValValueType,     class ...inputValProperties>
  void
  FunctionSpaceTools<DeviceType>::
  multiplyMeasure(       Kokkos::DynRankView<outputValValueType,   outputValProperties...>    outputVals,
                   const Kokkos::DynRankView<inputMeasureValueType,inputMeasureProperties...> inputMeasure,
                   const Kokkos::DynRankView<inputValValueType,    inputValProperties...>     inputVals ) {
    scalarMultiplyDataField( outputVals, 
                             inputMeasure, 
                             inputVals );
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<DeviceType>::
  scalarMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool reciprocal ) {
    ArrayTools<DeviceType>::scalarMultiplyDataField( outputFields,
                                              inputData, 
                                              inputFields, 
                                              reciprocal );    
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputDataValuetype,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<DeviceType>::
  scalarMultiplyDataData(       Kokkos::DynRankView<outputDataValuetype,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool reciprocal ) {
    ArrayTools<DeviceType>::scalarMultiplyDataData( outputData,
                                             inputDataLeft, 
                                             inputDataRight, 
                                             reciprocal );
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<DeviceType>::
  dotMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                        const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                        const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    ArrayTools<DeviceType>::dotMultiplyDataField( outputFields,
                                           inputData, 
                                           inputFields );
  } 
  
  // ------------------------------------------------------------------------------------
  
  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<DeviceType>::
  dotMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                       const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {
    ArrayTools<DeviceType>::dotMultiplyDataData( outputData,
                                          inputDataLeft, 
                                          inputDataRight );
  }

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<DeviceType>::
  vectorMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 3:
    case 4:
      ArrayTools<DeviceType>::crossProductDataField( outputFields,
                                              inputData, 
                                              inputFields );
      break;
    case 5:
      ArrayTools<DeviceType>::outerProductDataField( outputFields,
                                              inputData, 
                                              inputFields );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 3 && outRank > 5, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataField): Output container must have rank 3, 4 or 5.");
    }
    }
  } 

  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<DeviceType>::
  vectorMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {
    const auto outRank = outputData.rank();
    switch (outRank) {
    case 2:
    case 3:
      ArrayTools<DeviceType>::crossProductDataData( outputData,
                                             inputDataLeft, 
                                             inputDataRight );
      break;
    case 4:
      ArrayTools<DeviceType>::outerProductDataData( outputData,
                                             inputDataLeft, 
                                             inputDataRight );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 2 && outRank > 4, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataData): Output container must have rank 2, 3 or 4.");
    }
    }    
  } 
  
  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<DeviceType>::
  tensorMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const char transpose ) {

    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 4:
      ArrayTools<DeviceType>::matvecProductDataField( outputFields,
                                               inputData, 
                                               inputFields, 
                                               transpose );
      break;
    case 5:
      ArrayTools<DeviceType>::matmatProductDataField( outputFields,
                                               inputData, 
                                               inputFields, 
                                               transpose );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 4 && outRank > 5, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
    }
    }
  } 
  
  // ------------------------------------------------------------------------------------

  template<typename DeviceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<DeviceType>::
  tensorMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const char transpose ) {
    const auto outRank = outputData.rank();
    switch (outRank) {
    case 3:
      ArrayTools<DeviceType>::matvecProductDataData( outputData,
                                              inputDataLeft, 
                                              inputDataRight, 
                                              transpose );
      break;
    case 4:
      ArrayTools<DeviceType>::matmatProductDataData( outputData,
                                              inputDataLeft, 
                                              inputDataRight, 
                                              transpose );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 4 && outRank > 5, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
    }
    }
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {

   /**
       \brief Functor for applyLeftFieldSigns, see Intrepid2::FunctionSpaceTools for more
   */
    template<typename inoutOperatorViewType,
             typename fieldSignViewType>
    struct F_applyLeftFieldSigns {
            inoutOperatorViewType _inoutOperator;
      const fieldSignViewType     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      F_applyLeftFieldSigns( inoutOperatorViewType inoutOperator_,
                             fieldSignViewType     fieldSigns_ )
        : _inoutOperator(inoutOperator_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type nlbf = _inoutOperator.extent(1);
        const ordinal_type nrbf = _inoutOperator.extent(2);

        for (ordinal_type lbf=0;lbf<nlbf;++lbf)
          for (ordinal_type rbf=0;rbf<nrbf;++rbf)
            _inoutOperator(cl, lbf, rbf) *= _fieldSigns(cl, lbf);
      }
    };
  }

  template<typename DeviceType>
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<DeviceType>::
  applyLeftFieldSigns(       Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator,
                       const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input operator container must have rank 3.");
    INTREPID2_TEST_FOR_EXCEPTION( fieldSigns.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.extent(0) != fieldSigns.extent(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.extent(1) != fieldSigns.extent(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): First dimensions (number of left fields) of the operator and field signs containers must agree!");
#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inoutOperator)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(fieldSigns)::memory_space>::accessible;
    static_assert(are_accessible, "FunctionSpaceTools<DeviceType>::applyLeftFieldSigns(..): input/output views' memory spaces are not compatible with DeviceType");

    using FunctorType = FunctorFunctionSpaceTools::F_applyLeftFieldSigns
      /**/           <decltype(inoutOperator),decltype(fieldSigns)>;

    const ordinal_type C = inoutOperator.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    Kokkos::parallel_for( policy, FunctorType(inoutOperator, fieldSigns) );
  } 

  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {
   /**
       \brief Functor for applyRightFieldSigns, see Intrepid2::FunctionSpaceTools for more
   */
    template<typename inoutOperatorViewType,
             typename fieldSignViewType>
    struct F_applyRightFieldSigns {
      /**/  inoutOperatorViewType _inoutOperator;
      const fieldSignViewType     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      F_applyRightFieldSigns( inoutOperatorViewType inoutOperator_,
                              fieldSignViewType     fieldSigns_ )
        : _inoutOperator(inoutOperator_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type nlbf = _inoutOperator.extent(1);
        const ordinal_type nrbf = _inoutOperator.extent(2);

        for (ordinal_type lbf=0;lbf<nlbf;++lbf)
          for (ordinal_type rbf=0;rbf<nrbf;++rbf)
            _inoutOperator(cl, lbf, rbf) *= _fieldSigns(cl, rbf);
      }
    };
  }

  template<typename DeviceType>
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<DeviceType>::
  applyRightFieldSigns(       Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator,
                        const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input operator container must have rank 3.");
    INTREPID2_TEST_FOR_EXCEPTION( fieldSigns.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.extent(0) != fieldSigns.extent(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.extent(2) != fieldSigns.extent(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Second dimension of the operator container and first dimension of the field signs container (number of right fields) must agree!");
#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inoutOperator)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(fieldSigns)::memory_space>::accessible;
    static_assert(are_accessible, "FunctionSpaceTools<DeviceType>::applyRightFieldSigns(..): input/output views' memory spaces are not compatible with DeviceType");

    using FunctorType = FunctorFunctionSpaceTools::F_applyRightFieldSigns
                     <decltype(inoutOperator),decltype(fieldSigns)>;

    const ordinal_type C = inoutOperator.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    Kokkos::parallel_for( policy, FunctorType(inoutOperator, fieldSigns) );
  } 
  
  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {
   /**
       \brief Functor for applyFieldSigns, see Intrepid2::FunctionSpaceTools for more
   */
    template<typename inoutFunctionViewType,
             typename fieldSignViewType>
    struct F_applyFieldSigns {
            inoutFunctionViewType _inoutFunction;
      const fieldSignViewType     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      F_applyFieldSigns( inoutFunctionViewType inoutFunction_,
                         fieldSignViewType     fieldSigns_)
        : _inoutFunction(inoutFunction_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type nbfs = _inoutFunction.extent(1);
        const ordinal_type npts = _inoutFunction.extent(2);
        const ordinal_type iend = _inoutFunction.extent(3);
        const ordinal_type jend = _inoutFunction.extent(4);

        for (ordinal_type bf=0;bf<nbfs;++bf) 
          for (ordinal_type pt=0;pt<npts;++pt)
            for (ordinal_type i=0;i<iend;++i) 
              for (ordinal_type j=0;j<jend;++j) 
                _inoutFunction(cl, bf, pt, i, j) *= _fieldSigns(cl, bf);
      }
    };
  }

  template<typename DeviceType>
  template<typename inoutFunctionValueType, class ...inoutFunctionProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<DeviceType>::
  applyFieldSigns(       Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...> inoutFunction,
                   const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.rank() < 2 || inoutFunction.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input function container must have rank 2, 3, 4, or 5.");
    INTREPID2_TEST_FOR_EXCEPTION( fieldSigns.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.extent(0) != fieldSigns.extent(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the function and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.extent(1) != fieldSigns.extent(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): First dimensions (number of fields) of the function and field signs containers must agree!");
    
#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inoutFunction)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(fieldSigns)::memory_space>::accessible;
    static_assert(are_accessible, "FunctionSpaceTools<DeviceType>::applyFieldSigns(..): input/output views' memory spaces are not compatible with DeviceType");

    using FunctorType = FunctorFunctionSpaceTools::F_applyFieldSigns
                     <decltype(inoutFunction),decltype(fieldSigns)>;

    const ordinal_type C = inoutFunction.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    Kokkos::parallel_for( policy, FunctorType(inoutFunction, fieldSigns) );
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {

   /**
       \brief Functor to evaluate functions, see Intrepid2::FunctionSpaceTools for more
   */
    template<typename outputPointViewType,
             typename inputCoeffViewType,
             typename inputFieldViewType>
    struct F_evaluate {
            outputPointViewType _outputPointVals;
      const inputCoeffViewType  _inputCoeffs;
      const inputFieldViewType  _inputFields;
      
      KOKKOS_INLINE_FUNCTION
      F_evaluate( outputPointViewType outputPointVals_,
                  inputCoeffViewType  inputCoeffs_,
                  inputFieldViewType  inputFields_ )
        : _outputPointVals(outputPointVals_), _inputCoeffs(inputCoeffs_), _inputFields(inputFields_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type nbfs = _inputFields.extent(1);
        const ordinal_type npts = _inputFields.extent(2);

        const ordinal_type iend = _inputFields.extent(3);
        const ordinal_type jend = _inputFields.extent(4);
        
        for (ordinal_type bf=0;bf<nbfs;++bf) 
          for (ordinal_type pt=0;pt<npts;++pt)
            for (ordinal_type i=0;i<iend;++i) 
              for (ordinal_type j=0;j<jend;++j) 
                _outputPointVals(cl, pt, i, j) += _inputCoeffs(cl, bf) * _inputFields(cl, bf, pt, i, j);
      }
    };
  }

  template<typename DeviceType>
  template<typename outputPointValueType, class ...outputPointProperties,
           typename inputCoeffValueType,  class ...inputCoeffProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<DeviceType>::
  evaluate(       Kokkos::DynRankView<outputPointValueType,outputPointProperties...> outputPointVals,
            const Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>  inputCoeffs,
            const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inputFields.rank() < 3 || inputFields.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Input fields container must have rank 3, 4, or 5.");
    INTREPID2_TEST_FOR_EXCEPTION( inputCoeffs.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Input coefficient container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( outputPointVals.rank() != (inputFields.rank()-1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Output values container must have rank one less than the rank of the input fields container.");
    INTREPID2_TEST_FOR_EXCEPTION( inputCoeffs.extent(0) != inputFields.extent(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the coefficient and fields input containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inputCoeffs.extent(1) != inputFields.extent(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): First dimensions (number of fields) of the coefficient and fields input containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( outputPointVals.extent(0) != inputFields.extent(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the input fields container and the output values container must agree!");
    for (size_type i=1;i<outputPointVals.rank();++i)
      INTREPID2_TEST_FOR_EXCEPTION( outputPointVals.extent(i) != inputFields.extent(i+1), std::invalid_argument, 
                                    ">>> ERROR (FunctionSpaceTools::evaluate): outputPointVals dimension(i) does not match to inputFields dimension(i+1).");
#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(outputPointVals)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inputCoeffs)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(inputFields)::memory_space>::accessible;
    static_assert(are_accessible, "FunctionSpaceTools<DeviceType>::evaluate(..): input/output views' memory spaces are not compatible with DeviceType");

    using FunctorType = FunctorFunctionSpaceTools::F_evaluate
                     <decltype(outputPointVals),decltype(inputCoeffs),decltype(inputFields)>;
    
    const ordinal_type C = inputFields.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    Kokkos::parallel_for( policy, FunctorType(outputPointVals, inputCoeffs, inputFields) );
  }
  
  // ------------------------------------------------------------------------------------
  
} // end namespace Intrepid2

#endif
