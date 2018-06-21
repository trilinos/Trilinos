// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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

/** \file   Intrepid2_FunctionSpaceToolsDef.hpp
    \brief  Definition file for the Intrepid2::FunctionSpaceTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_FUNCTIONSPACETOOLS_DEF_HPP__
#define __INTREPID2_FUNCTIONSPACETOOLS_DEF_HPP__

namespace Intrepid2 {

  // ------------------------------------------------------------------------------------
  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  FunctionSpaceTools<SpT>::
  HGRADtransformVALUE(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
                       const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
    ArrayTools<SpT>::cloneFields(output, input);
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorFunctionSpaceTools {
   /**
       \brief Functor for calculation HGRADtransformGRAD, see Intrepid2::FunctionSpaceTools for more
   */
    template <typename outputViewType,
              typename jacInverseViewType,
              typename inputViewType,
              ordinal_type spaceDim>
    struct F_HGRADtransformGRAD {
            outputViewType     _output;
      const jacInverseViewType  _jacInverse;
      const inputViewType _input;

      // output CPDD, left CPDD or PDD, right FPD
      KOKKOS_INLINE_FUNCTION
      F_HGRADtransformGRAD(outputViewType     output_,
                           jacInverseViewType  jacInverse_,
                           inputViewType input_)
        : _output(output_), 
          _jacInverse(jacInverse_), 
          _input(input_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type bf,
                      const ordinal_type pt) const {
        /* */ auto y = Kokkos::subview(_output,     cl, bf, pt, Kokkos::ALL());
        const auto A = Kokkos::subview(_jacInverse, cl,     pt, Kokkos::ALL(), Kokkos::ALL());
        const auto x = Kokkos::subview(_input,      bf,     pt, Kokkos::ALL());
        
        if (spaceDim == 2) {
          Kernels::Serial::matvec_trans_product_d2( y, A, x );
        } else {
          Kernels::Serial::matvec_trans_product_d3( y, A, x );
        }
      }
    };
  }
  
  template<typename SpT>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HGRADtransformGRAD(       Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                      const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                      const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputVals.rank()       != 3 ||
                                    jacobianInverse.rank() != 4 ||
                                    outputVals.rank()      != 4, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): Ranks are not compatible.");
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(0) != jacobianInverse.extent(0), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): Cell dimension does not match.");
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(1) != inputVals.extent(0), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): Field dimension does not match.");
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(2)      != inputVals.extent(1) || 
                                    jacobianInverse.extent(1) != inputVals.extent(1), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): Point dimension does not match.");
      const auto spaceDim = outputVals.extent(3);
      INTREPID2_TEST_FOR_EXCEPTION( jacobianInverse.extent(2) != spaceDim || 
                                    jacobianInverse.extent(3) != spaceDim || 
                                    inputVals.extent(2)       != spaceDim , std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): matvec dimensions are not compatible.");
    }
#endif
    ArrayTools<SpT>::matvecProductDataField(outputVals, jacobianInverse, inputVals, 'T');

    // this modification is for 2d and 3d (not 1d)
    // this is an attempt to measure the overhead of subview of dynrankview. 

    // typedef       Kokkos::DynRankView<outputValValueType,outputValProperties...> outputViewType;
    // typedef const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacInverseViewType;
    // typedef const Kokkos::DynRankView<inputValValueType,inputValProperties...>  inputViewType;

    // const ordinal_type 
    //   C = outputVals.extent(0),
    //   F = outputVals.extent(1),
    //   P = outputVals.extent(2);

    // using range_policy_type = Kokkos::Experimental::MDRangePolicy
    //   < SpT, Kokkos::Experimental::Rank<3>, Kokkos::IndexType<ordinal_type> >;
    // range_policy_type policy( { 0, 0, 0 },
    //                           { C, F, P } );

    // const ordinal_type spaceDim = inputVals.extent(2);
    // switch (spaceDim) {
    // case 2: {
    //   typedef FunctorFunctionSpaceTools::F_HGRADtransformGRAD<outputViewType, jacInverseViewType, inputViewType, 2> FunctorType;
    //   Kokkos::parallel_for( policy, FunctorType(outputVals, jacobianInverse, inputVals) );
    //   break;
    // }
    // case 3: {
    //   typedef FunctorFunctionSpaceTools::F_HGRADtransformGRAD<outputViewType, jacInverseViewType, inputViewType, 3> FunctorType;
    //   Kokkos::parallel_for( policy, FunctorType(outputVals, jacobianInverse, inputVals) );
    //   break;
    // }
    // default: {
    //   INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
    //                                 ">>> ERROR (FunctionSpaceTools::HGRADtransformGRAD): spaceDim is not 2 or 3.");
    //   break;
    // }
    // }
  }
  
  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HCURLtransformVALUE(       Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                       const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                       const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals ) {
    ArrayTools<SpT>::matvecProductDataField(outputVals, jacobianInverse, inputVals, 'T');
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HCURLtransformCURL(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<SpT>::matvecProductDataField(outputVals, jacobian, inputVals, 'N');
    ArrayTools<SpT>::scalarMultiplyDataField(outputVals, jacobianDet, outputVals, true);
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HDIVtransformVALUE(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<SpT>::matvecProductDataField(outputVals, jacobian, inputVals, 'N');
    ArrayTools<SpT>::scalarMultiplyDataField(outputVals, jacobianDet, outputVals, true);
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HDIVtransformDIV(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                    const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                    const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<SpT>::scalarMultiplyDataField(outputVals, jacobianDet, inputVals, true);
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<SpT>::
  HVOLtransformVALUE(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<SpT>::scalarMultiplyDataField(outputVals, jacobianDet, inputVals, true);
  }
  
  // ------------------------------------------------------------------------------------

  template<typename SpT>  
  template<typename outputValueValueType, class ...outputValueProperties,
           typename leftValueValueType,   class ...leftValueProperties,
           typename rightValueValueType,  class ...rightValueProperties>
  void
  FunctionSpaceTools<SpT>::
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
      ArrayTools<SpT>::contractDataDataScalar( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;
    case 13:
      ArrayTools<SpT>::contractDataDataVector( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;
    case 14:
      ArrayTools<SpT>::contractDataDataTensor( outputValues,
                                               leftValues,
                                               rightValues,
                                               sumInto );
      break;

      // DataField
    case 22:
      ArrayTools<SpT>::contractDataFieldScalar( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;
    case 23:
      ArrayTools<SpT>::contractDataFieldVector( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;
    case 24:
      ArrayTools<SpT>::contractDataFieldTensor( outputValues,
                                                leftValues,
                                                rightValues,
                                                sumInto );
      break;

      // FieldField
    case 33:
      ArrayTools<SpT>::contractFieldFieldScalar( outputValues,
                                                 leftValues,
                                                 rightValues,
                                                 sumInto );
      break;
    case 34:
      ArrayTools<SpT>::contractFieldFieldVector( outputValues,
                                                 leftValues,
                                                 rightValues,
                                                 sumInto );
      break;
    case 35:
      ArrayTools<SpT>::contractFieldFieldTensor( outputValues,
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

  template<typename SpT>  
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputDetValueType,    class ...inputDetProperties,
           typename inputWeightValueType, class ...inputWeightProperties>
  bool
  FunctionSpaceTools<SpT>::
  computeCellMeasure(       Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<inputDetValueType,   inputDetProperties...>    inputDet,
                      const Kokkos::DynRankView<inputWeightValueType,inputWeightProperties...> inputWeights ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inputDet.rank()     != 2 || 
                                    inputWeights.rank() != 1 || 
                                    outputVals.rank()   != 2, std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Ranks are not compatible.");
      INTREPID2_TEST_FOR_EXCEPTION( outputVals.extent(0) != inputDet.extent(0), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Cell dimension does not match.");
      INTREPID2_TEST_FOR_EXCEPTION( inputDet.extent(1)      !=  outputVals.extent(1) ||
                                    inputWeights.extent(0) !=  outputVals.extent(1), std::invalid_argument,
                                    ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Point dimension does not match.");
    }
#endif
    typedef          Kokkos::DynRankView<outputValValueType,  outputValProperties...>         outputValViewType;
    typedef          Kokkos::DynRankView<inputDetValueType,   inputDetProperties...>          inputDetViewType;
    typedef          Kokkos::DynRankView<inputWeightValueType,inputWeightProperties...>       inputWeightViewType;
    typedef          FunctorFunctionSpaceTools::F_computeCellMeasure
                     <outputValViewType,inputDetViewType,inputWeightViewType> FunctorType;
    
    const ordinal_type C = inputDet.extent(0);
    Kokkos::RangePolicy<SpT,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    
    typename FunctorType::value_type hasNegativeDet = false;
    Kokkos::parallel_reduce( policy, FunctorType(outputVals, inputDet, inputWeights), hasNegativeDet );
    
    return hasNegativeDet;
  } 

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightProperties,
           typename scratchValueType,     class ...scratchProperties>
  void 
  FunctionSpaceTools<SpT>::
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
    INTREPID2_TEST_FOR_EXCEPTION( scratch.span() < inputJac.span(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Scratch storage must be greater than or equal to inputJac's one.");
#endif

    // face normals (reshape scratch)
    // Kokkos::DynRankView<scratchValueType,scratchProperties...> faceNormals(scratch.data(), 
    //                                                                        inputJac.extent(0), 
    //                                                                        inputJac.extent(1), 
    //                                                                        inputJac.extent(2));
    auto vcprop = Kokkos::common_view_alloc_prop(scratch);
    //typedef Kokkos::DynRankView<scratchValueType, typename decltype(scratch)::memory_space> viewType;
    typedef Kokkos::DynRankView<scratchValueType, SpT> viewType;
    viewType faceNormals(Kokkos::view_wrap(scratch.data(), vcprop),
                         inputJac.extent(0),
                         inputJac.extent(1),
                         inputJac.extent(2));

    // compute normals
    CellTools<SpT>::getPhysicalFaceNormals(faceNormals, inputJac, whichFace, parentCell);

    // compute lenghts of normals
    RealSpaceTools<SpT>::vectorNorm(outputVals, faceNormals, NORM_TWO);

    // multiply with weights
    ArrayTools<SpT>::scalarMultiplyDataData(outputVals, outputVals, inputWeights);
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightProperties,
           typename scratchValueType,     class ...scratchProperties>
  void 
  FunctionSpaceTools<SpT>::
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
    INTREPID2_TEST_FOR_EXCEPTION( scratch.span() < inputJac.span(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Scratch storage must be greater than or equal to inputJac'one.");
#endif

    // edge tangents (reshape scratch)
    // Kokkos::DynRankView<scratchValueType,scratchProperties...> edgeTangents(scratch.data(), 
    //                                                                         inputJac.extent(0), 
    //                                                                         inputJac.extent(1), 
    //                                                                         inputJac.extent(2));
    auto vcprop = Kokkos::common_view_alloc_prop(scratch);
    //typedef Kokkos::DynRankView<scratchValueType, typename decltype(scratch)::memory_space> viewType;
    typedef Kokkos::DynRankView<scratchValueType, SpT> viewType;
    viewType edgeTangents(Kokkos::view_wrap(scratch.data(), vcprop),
                         inputJac.extent(0),
                         inputJac.extent(1),
                         inputJac.extent(2));

    // compute normals
    CellTools<SpT>::getPhysicalEdgeTangents(edgeTangents, inputJac, whichEdge, parentCell);

    // compute lenghts of tangents
    RealSpaceTools<SpT>::vectorNorm(outputVals, edgeTangents, NORM_TWO);

    // multiply with weights
    ArrayTools<SpT>::scalarMultiplyDataData(outputVals, outputVals, inputWeights);
  }

  // ------------------------------------------------------------------------------------  

  template<typename SpT>
  template<typename outputValValueType,    class ...outputValProperties,
           typename inputMeasureValueType, class ...inputMeasureProperties,
           typename inputValValueType,     class ...inputValProperteis>
  void
  FunctionSpaceTools<SpT>::
  multiplyMeasure(       Kokkos::DynRankView<outputValValueType,   outputValProperties...>    outputVals,
                   const Kokkos::DynRankView<inputMeasureValueType,inputMeasureProperties...> inputMeasure,
                   const Kokkos::DynRankView<inputValValueType,    inputValProperteis...>     inputVals ) {
    scalarMultiplyDataField( outputVals, 
                             inputMeasure, 
                             inputVals );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<SpT>::
  scalarMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool reciprocal ) {
    ArrayTools<SpT>::scalarMultiplyDataField( outputFields, 
                                              inputData, 
                                              inputFields, 
                                              reciprocal );    
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputDataValuetype,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<SpT>::
  scalarMultiplyDataData(       Kokkos::DynRankView<outputDataValuetype,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const bool reciprocal ) {
    ArrayTools<SpT>::scalarMultiplyDataData( outputData, 
                                             inputDataLeft, 
                                             inputDataRight, 
                                             reciprocal );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<SpT>::
  dotMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                        const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                        const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    ArrayTools<SpT>::dotMultiplyDataField( outputFields, 
                                           inputData, 
                                           inputFields );
  } 
  
  // ------------------------------------------------------------------------------------
  
  template<typename SpT>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<SpT>::  
  dotMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                       const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {
    ArrayTools<SpT>::dotMultiplyDataData( outputData, 
                                          inputDataLeft, 
                                          inputDataRight );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<SpT>::  
  vectorMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 3:
    case 4:
      ArrayTools<SpT>::crossProductDataField( outputFields, 
                                              inputData, 
                                              inputFields );
      break;
    case 5:
      ArrayTools<SpT>::outerProductDataField( outputFields, 
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

  template<typename SpT>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<SpT>::  
  vectorMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight ) {
    const auto outRank = outputData.rank();
    switch (outRank) {
    case 2:
    case 3:
      ArrayTools<SpT>::crossProductDataData( outputData, 
                                             inputDataLeft, 
                                             inputDataRight );
      break;
    case 4:
      ArrayTools<SpT>::outerProductDataData( outputData, 
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

  template<typename SpT>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<SpT>::  
  tensorMultiplyDataField(       Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties...>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const char transpose ) {

    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 4:
      ArrayTools<SpT>::matvecProductDataField( outputFields, 
                                               inputData, 
                                               inputFields, 
                                               transpose );
      break;
    case 5:
      ArrayTools<SpT>::matmatProductDataField( outputFields, 
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

  template<typename SpT>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<SpT>::  
  tensorMultiplyDataData(       Kokkos::DynRankView<outputDataValueType,    outputDataProperties...>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties...>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties...> inputDataRight,
                          const char transpose ) {
    const auto outRank = outputData.rank();
    switch (outRank) {
    case 3:
      ArrayTools<SpT>::matvecProductDataData( outputData, 
                                              inputDataLeft, 
                                              inputDataRight, 
                                              transpose );
      break;
    case 4:
      ArrayTools<SpT>::matmatProductDataData( outputData, 
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

  template<typename SpT>  
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<SpT>::  
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

    typedef          Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...>        inoutOperatorViewType;
    typedef          Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>            fieldSignViewType;
    typedef          FunctorFunctionSpaceTools::F_applyLeftFieldSigns
      /**/           <inoutOperatorViewType,fieldSignViewType>                                     FunctorType;
    typedef typename ExecSpace<typename inoutOperatorViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

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

  template<typename SpT>  
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<SpT>::  
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

    typedef          Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...>        inoutOperatorViewType;
    typedef          Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>            fieldSignViewType;
    typedef          FunctorFunctionSpaceTools::F_applyRightFieldSigns
                     <inoutOperatorViewType,fieldSignViewType>                                     FunctorType;
    typedef typename ExecSpace<typename inoutOperatorViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

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

  template<typename SpT>  
  template<typename inoutFunctionValueType, class ...inoutFunctionProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<SpT>::  
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

    typedef          Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...>        inoutFunctionViewType;
    typedef          Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>            fieldSignViewType;
    typedef          FunctorFunctionSpaceTools::F_applyFieldSigns
                     <inoutFunctionViewType,fieldSignViewType>                                     FunctorType;
    typedef typename ExecSpace<typename inoutFunctionViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

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

  template<typename SpT>    
  template<typename outputPointValueType, class ...outputPointProperties,
           typename inputCoeffValueType,  class ...inputCoeffProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<SpT>::  
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

    typedef          Kokkos::DynRankView<outputPointValueType,outputPointProperties...>         outputPointValViewType;
    typedef          Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>          inputCoeffViewType;
    typedef          Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>          inputFieldViewType; 
    typedef          FunctorFunctionSpaceTools::F_evaluate
                     <outputPointValViewType,inputCoeffViewType,inputFieldViewType>             FunctorType;
    typedef typename ExecSpace<typename inputCoeffViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;
    
    const ordinal_type C = inputFields.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, C);
    Kokkos::parallel_for( policy, FunctorType(outputPointVals, inputCoeffs, inputFields) );
  }
  
  // ------------------------------------------------------------------------------------
  
} // end namespace Intrepid2

#endif
