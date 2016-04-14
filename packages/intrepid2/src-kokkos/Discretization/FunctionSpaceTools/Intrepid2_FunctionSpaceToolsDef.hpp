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

/** \file   Intrepid_FunctionSpaceToolsDef.hpp
    \brief  Definition file for the Intrepid2::FunctionSpaceTools class.
    \author Created by P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/


namespace Intrepid2 {

  
  template<typename ExecSpaceType>
  template<typename outputValValueType, class ...outputValProperties,
           typename inputValValueType,  class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HGRADtransformVALUE( /**/  Kokkos::DynRankView<outputValValueType,outputValProperties> outputVals,
                       const Kokkos::DynRankView<inputValValueType, inputValProperties>  inputVals ) {
    ArrayTools<ExecSpaceType>::cloneFields(outputVals, inputVals);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HGRADtransformGRAD( /**/  Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                      const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                      const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outputVals, jacobianInverse, inputVals, transpose);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,       class ...outputValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inputValValueType,        class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HCURLtransformVALUE( /**/  Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                       const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                       const Kokkos::DynRankView<inputValValueType,       inputValProperties...>        inputVals,
                       const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outputVals, jacobianInverse, inputVals, transpose);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HCURLtransformCURL( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outputVals, jacobian, inputVals, transpose);
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outputVals, jacobianDet, outputVals, true);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HDIVtransformVALUE( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outputVals, jacobian, inputVals, transpose);
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outputVals, jacobianDet, outputVals, true);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HDIVtransformDIV( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                    const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                    const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outputVals, jacobianDet, inputVals, true);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inputValValueType,    class ...inputValProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  HVOLtransformVALUE( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>   outputVals,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inputValValueType,   inputValProperties...>    inputVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outputVals, jacobianDet, inputVals, true);
  }
  

  template<typename ExecSpaceType>  
  template<typename outputValueValueType, class ...outputValueProperties,
           typename leftValueValueType,   class ...leftValueProperties,
           typename rightValueValueType,  class ...rightValueProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  integrate( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties> outputValues,
             const Kokkos::DynRankView<leftValueValueType,  leftValueProperties>   leftValues,
             const Kokkos::DynRankView<rightValueValueType, rightValueProperties>  rightValues,
             const bool         sumInto ) {
    
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
      // *** DataData
    case 12:
      ArrayTools<ExecSpaceType>::contractDataDataScalar( outputValues,
                                                         leftValues,
                                                         rightValues,
                                                         sumInto );
      break;
    case 13:
      ArrayTools<ExecSpaceType>::contractDataDataVector( outputValues,
                                                         leftValues,
                                                         rightValues,
                                                         sumInto );
      break;
    case 14:
      ArrayTools<ExecSpaceType>::contractDataDataTensor( outputValues,
                                                         leftValues,
                                                         rightValues,
                                                         sumInto );
      break;

      // *** DataField
    case 22:
      ArrayTools<ExecSpaceType>::contractDataFieldScalar( outputValues,
                                                          leftValues,
                                                          rightValues,
                                                          sumInto );
      break;
    case 23:
      ArrayTools<ExecSpaceType>::contractDataFieldVector( outputValues,
                                                          leftValues,
                                                          rightValues,
                                                          sumInto );
      break;
    case 24:
      ArrayTools<ExecSpaceType>::contractDataFieldTensor( outputValues,
                                                          leftValues,
                                                          rightValues,
                                                          sumInto );
      break;

      // *** FieldField
    case 32:
      ArrayTools<ExecSpaceType>::contractFieldFieldScalar( outputValues,
                                                           leftValues,
                                                           rightValues,
                                                           sumInto );
      break;
    case 33:
      ArrayTools<ExecSpaceType>::contractFieldFieldVector( outputValues,
                                                           leftValues,
                                                           rightValues,
                                                           sumInto );
      break;
    case 34:
      ArrayTools<ExecSpaceType>::contractFieldFieldTensor( outputValues,
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
  

  template<typename ExecSpaceType>  
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputDetValueType,    class ...inputDetPropertes,
           typename inputWeightValueType, class ...inputWeightPropertes>
  void
  FunctionSpaceTools<ExecSpaceType>::
  computeCellMeasure( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>  outputVals,
                      const Kokkos::DynRankView<inputDetValueType,   inputDetPropertes...>    inputDet,
                      const Kokkos::DynRankView<inputWeightValueType,inputWeightPropertes...> inputWeights ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( (inDet.rank() != 1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Input determinants container must have rank 1.");
    
#endif
    struct Functor {
      typedef bool value_type;

      Kokkos::DynRankView<outputValValueType,  outputValProperties...>  _outputVals;
      Kokkos::DynRankView<inputDetValueType,   inputDetPropertes...>    _inputDet;
      Kokkos::DynRankView<inputWeightValueType,inputWeightPropertes...> _inputWeights;
      
      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputValValueType,  outputValProperties...>  outputVals_,
              Kokkos::DynRankView<inputDetValueType,   inputDetPropertes...>    inputDet_,
              Kokkos::DynRankView<inputWeightValueType,inputWeightPropertes...> inputWeights_)
        : _outputVals(outputVals_), _inputDet(inputDet_), _inputWeights(inputWeights_) {}

      KOKKOS_INLINE_FUNCTION
      void init(value_type &dst ) const {
        dst = false;
      }

      KOKKOS_INLINE_FUNCTION
      void join( volatile value_type &dst,
                 const volatile value_type &src ) const {
        dst |= src;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cell, value_type &dst) const {
        dst |= (_inDet(cell) < 0.0);
      }
    };
    const size_type loopSize = inDet.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

    // instead flip the sign of negative input, it reports an error
    bool hasNegativeDet = false;
    Kokkos::parallel_reduce( policy, Functor(outputVals, inputDet, inputWeights), hasNegativeDet );
    
    INTREPID2_TEST_FOR_EXCEPTION( hasNegativeDet, std::runtime_error,
                                  ">>> ERROR (FunctionSpaceTools::computeCellMeasure): inputDet has negative values");
    
    ArrayTools<ExecSpaceType>::scalarMultiplyDataData(outputVals, inDet, inWeights);
  } 


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightPropertes,
           typename scratchValueType,     class ...scratchProperties>
  void 
  FunctionSpaceTools<ExecSpaceType>::
  computeFaceMeasure(/**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>  outputVals,
                     const Kokkos::DynRankView<inputJacValueTypem,  inputJacProperties...>   inputJac,
                     const Kokkos::DynRankView<inputWeightValueType,inputWeightPropertes...> inputWeights,
                     const int                   whichFace,
                     const shards::CellTopology &parentCell,
                     const Kokkos::DynRankView<scratchValueType,    scratchProperties...>   scratch) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inJac.rank() != 4, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Input Jacobian container must have rank 4.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Scratch space always has a unit rank.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.span() >= inJac.span(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Scratch space must be greater than inJac.");
#endif

    // face normals (reshape scratch)
    Kokkos::DynRankView<scratchValueType,scratchProperties...> faceNormals(scratch.data(), 
                                                                           inJac.dimension(0), 
                                                                           inJac.dimension(1), 
                                                                           inJac.dimension(2));
    
    // compute normals
    CellTools<ExecSpaceType>::getPhysicalFaceNormals(faceNormals, inJac, whichFace, parentCell);

    // compute lenghts of normals
    RealSpaceTools<ExecSpaceType>::vectorNorm(outputVals, faceNormals, NORM_TWO);

    // multiply with weights
    ArrayTools<ExecSpaceType>::scalarMultiplyDataData(outputVals, outputVals, inWeights);
  }


  template<typename ExecSpaceType>
  template<typename outputValValueType,   class ...outputValProperties,
           typename inputJacValueType,    class ...inputJacProperties,
           typename inputWeightValueType, class ...inputWeightPropertes,
           typename scratchValueType,     class ...scratchPropertes>
  void 
  FunctionSpaceTools<ExecSpaceType>::
  computeEdgeMeasure( /**/  Kokkos::DynRankView<outputValValueType,  outputValProperties...>  outputVals,
                      const Kokkos::DynRankView<inputJacValueTypem,  inputJacProperties...>   inputJac,
                      const Kokkos::DynRankView<inputWeightValueType,inputWeightPropertes...> inputWeights,
                      const int                   whichEdge,
                      const shards::CellTopology &parentCell,
                      const Kokkos::DynRankView<scratchValueType,    scratchProperties...>    scratch ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Input Jacobian container must have rank 4.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Scratch space always has a unit rank.");
    INTREPID2_TEST_FOR_EXCEPTION( scratch.span() >= inJac.span(), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Scratch space must be greater than inJac.");
#endif

    // edge tangents (reshape scratch)
    Kokkos::DynRankView<scratchValueType,scratchProperties...> edgeTangents(scratch.data(), 
                                                                            inJac.dimension(0), 
                                                                            inJac.dimension(1), 
                                                                            inJac.dimension(2));
    
    // compute normals
    CellTools<ExecSpaceType>::getPhysicalEdgeTangents(edgeTangents, inJac, whichEdge, parentCell);

    // compute lenghts of tangents
    RealSpaceTools<ExecSpaceType>::vectorNorm(outputVals, edgeTangents, NORM_TWO);

    // multiply with weights
    ArrayTools<ExecSpaceType>::scalarMultiplyDataData<Scalar>(outputVals, outputVals, inWeights);
  }

  
  template<typename ExecSpaceType>
  template<typename outputValValueType,    class ...outputValProperties,
           typename inputMeasureValueType, class ...inputMeasureProperties,
           typename inputValValueType,     class ...inputValProperteis>
  void
  FunctionSpaceTools<ExecSpaceType>::
  multiplyMeasure( /**/  Kokkos::DynRankView<outputValValueType,      outputValProperties...>       outputVals,
                   const Kokkos::DynRankView<inputMeasureValueType,inputMeasureProperties...> inputMeasure,
                   const Kokkos::DynRankView<inputValValueType,    inputValProperteis...>     inputVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField( outputVals, 
                                                        inputMeasure, 
                                                        inputVals );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataPropertes,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  scalarMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties...> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataPropertes...>    inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields,
                           const bool reciprocal ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField( outputFields, 
                                                        inputData, 
                                                        inputFields, 
                                                        reciprocal );
  }


  template<typename ExecSpaceType>
  template<typename outputDataValuetype,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  scalarMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValuetype,    outputDataProperties>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties> inputDataRight,
                          const bool reciprocal ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataData( outputData, 
                                                       inputDataLeft, 
                                                       inputDataRight, 
                                                       reciprocal );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::
  dotMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties> outputFields,
                        const Kokkos::DynRankView<inputDataValueType,  inputDataProperties>   inputData,
                        const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties>  inputFields ) {
    ArrayTools<ExecSpaceType>::dotMultiplyDataField( outputFields, 
                                                     inputData, 
                                                     inputFields );
  } 
  
  
  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  dotMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties>     outputData,
                       const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties>  inputDataLeft,
                       const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties> inputDataRight ) {
    ArrayTools<ExecSpaceType>::dotMultiplyDataData( outputData, 
                                                    inputDataLeft, 
                                                    inputDataRight );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  vectorMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties>  inputFields ) {
    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 3:
    case 4:
      ArrayTools<ExecSpaceType>::crossProductDataField( outputFields, 
                                                        inputData, 
                                                        inputFields );
      break;
    case 5:
      ArrayTools<ExecSpaceType>::outerProductDataField( outputFields, 
                                                        inputData, 
                                                        inputFields );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 3 && outRank > 5, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataField): Output container must have rank 3, 4 or 5.");
    }
    }
  } 


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  vectorMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties> inputDataRight ) {
    const auto outRank = outputData.rank();
    switch (outRank) {
    case 2:
    case 3:
      ArrayTools<ExecSpaceType>::crossProductDataData( outputData, 
                                                       inputDataLeft, 
                                                       inputDataRight );
      break;
    case 4:
      ArrayTools<ExecSpaceType>::outerProductDataData( outputData, 
                                                       inputDataLeft, 
                                                       inputDataRight );
      break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( outRank < 2 && outRank > 4, std::runtime_error,
                                    ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataData): Output container must have rank 2, 3 or 4.");
    }
    }    
  } 
  

  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  tensorMultiplyDataField( /**/  Kokkos::DynRankView<outputFieldValueType,outputFieldProperties> outputFields,
                           const Kokkos::DynRankView<inputDataValueType,  inputDataProperties>   inputData,
                           const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties>  inputFields,
                           const char transpose = 'N') {

    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 4:
      ArrayTools<ExecSpaceType>::matvecProductDataField( outputFields, 
                                                         inputData, 
                                                         inputFields, 
                                                         transpose );
      break;
    case 5:
      ArrayTools<ExecSpaceType>::matmatProductDataField( outputFields, 
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
  


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  tensorMultiplyDataData( /**/  Kokkos::DynRankView<outputDataValueType,    outputDataProperties>     outputData,
                          const Kokkos::DynRankView<inputDataLeftValueType, inputDataLeftProperties>  inputDataLeft,
                          const Kokkos::DynRankView<inputDataRightValueType,inputDataRightProperties> inputDataRight,
                          const char transpose = 'N' ) {
    const auto outRank = outputFields.rank();
    switch (outRank) {
    case 3:
      ArrayTools<ExecSpaceType>::matvecProductDataData( outputData, 
                                                        inputDataLeft, 
                                                        inputDataRight, 
                                                        transpose );
      break;
    case 4:
      ArrayTools<ExecSpaceType>::matmatProductDataData( outputData, 
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


  template<typename ExecSpaceType>  
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  applyLeftFieldSigns( /**/  Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator,
                       const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input operator container must have rank 3.");
    INTREPID2_TEST_FOR_EXCEPTION( fieldSigns.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.dimension(0) != fieldSigns.dimension(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.dimension(1) != fieldSigns.dimension(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): First dimensions (number of left fields) of the operator and field signs containers must agree!");
#endif

    struct Functor {
      Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> _inoutOperator;
      Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator_,
              Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns_)
        : _inoutOperator(inoutOperator_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cell) const {
        const size_type nlbf = _inoutOperator.dimension(1);
        const size_type nrbf = _inoutOperator.dimension(2);

        for (size_type lbf=0;lbf<nlbf;++lbf)
          for (size_type rbf=0;rbf<nrbf;++rbf)
            _inoutOperator(cell, lbf, rbf) *= _fieldSigns(cell, lbf);
      }
    };
    const size_type loopSize = inoutOperator.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(inoutOperator, fieldSigns) );
  } 
  

  template<typename ExecSpaceType>  
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::  
  applyRightFieldSigns( /**/  Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator,
                        const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( getrank(inoutOperator) != 3, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input operator container must have rank 3.");
    INTREPID2_TEST_FOR_EXCEPTION( getrank(fieldSigns) != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.dimension(0) != fieldSigns.dimension(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutOperator.dimension(2) != fieldSigns.dimension(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Second dimension of the operator container and first dimension of the field signs container (number of right fields) must agree!");
#endif

    struct Functor {
      Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> _inoutOperator;
      Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator_,
              Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns_)
        : _inoutOperator(inoutOperator_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cell) const {
        const size_type nlbf = _inoutOperator.dimension(1);
        const size_type nrbf = _inoutOperator.dimension(2);

        for (size_type lbf=0;lbf<nlbf;++lbf)
          for (size_type rbf=0;rbf<nrbf;++rbf)
            _inoutOperator(cell, lbf, rbf) *= _fieldSigns(cell, rbf);
      }
    };
    const size_type loopSize = inoutOperator.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(inoutOperator, fieldSigns) );
  } 
  

  template<typename ExecSpaceType>  
  template<typename inoutFunctionValueType, class ...inoutFunctionProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  applyFieldSigns( /**/  Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...> inoutFunction,
                   const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.rank() < 2 || inoutFunction.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input function container must have rank 2, 3, 4, or 5.");
    INTREPID2_TEST_FOR_EXCEPTION( fieldSigns.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.dimension(0) != fieldSigns.dimension(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the function and field signs containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inoutFunction.dimension(1) != fieldSigns.dimension(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::applyFieldSigns): First dimensions (number of fields) of the function and field signs containers must agree!");
    
#endif
    
    struct Functor {
      Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...> _inoutFunction;
      Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     _fieldSigns;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...> inoutFunction_,
              Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns_)
        : _inoutFunction(inoutFunction_), _fieldSigns(fieldSigns_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cell) const {
        const size_type nbfs = _inoutFunction.dimension(1);
        const size_type npts = _inoutFunction.dimension(2);
        const size_type iend = _inoutFunction.dimension(3);
        const size_type jend = _inoutFunction.dimension(4);

        for (size_type bf=0;bf<nbfs;++bf) 
          for (size_type pt=0;pt<npts;++pt)
            for (size_type i=0;i<iend;++i) 
              for (size_type j=0;j<jend;++j) 
                inoutFunction(cell, bf, pt, i, j) *= fieldSigns(cell, bf);
      }
    };
    const size_type loopSize = inoutFunction.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(inoutFunction, fieldSigns) );
  }


  template<typename ExecSpaceType>    
  template<typename outputPointValueType, class ...outputPointProperties,
           typename inputCoeffValueType,  class ...inputCoeffProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  void
  FunctionSpaceTools<ExecSpaceType>::  
  evaluate( /**/  Kokkos::DynRankView<outputPointValueType,outputPointProperties...> outputPointVals,
            const Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>  inputCoeffs,
            const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inFields.rank() < 3 || inFields.rank() > 5, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Input fields container must have rank 3, 4, or 5.");
    INTREPID2_TEST_FOR_EXCEPTION( inCoeffs.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Input coefficient container must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( outPointVals.rank() != (inFields.rank()-1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Output values container must have rank one less than the rank of the input fields container.");
    INTREPID2_TEST_FOR_EXCEPTION( inCoeffs.dimension(0) != inFields.dimension(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the coefficient and fields input containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( inCoeffs.dimension(1) != inFields.dimension(1), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): First dimensions (number of fields) of the coefficient and fields input containers must agree!");
    INTREPID2_TEST_FOR_EXCEPTION( outPointVals.dimension(0) != inFields.dimension(0), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the input fields container and the output values container must agree!");
    for (auto i=1;i<outPointVals.rank();++i) 
      INTREPID2_TEST_FOR_EXCEPTION( outPointVals.dimension(i) != inFields.dimension(i+1), std::invalid_argument, 
                                    ">>> ERROR (FunctionSpaceTools::evaluate): outPointVals dimension(i) does not match to inFields dimension(i+1).");
                                    errmsg );
#endif
    
    struct Functor {
      Kokkos::DynRankView<outputPointValueType,outputPointProperties...> _outputPointVals;
      Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>  _inputCoeffs;
      Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  _inputFields;
      
      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<outputPointValueType,outputPointProperties...> outputPointVals_,
              Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>  inputCoeffs_,
              Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields_)
        : _outputPointVals(outputPointVals_), _inputCoeffs(inputCoeffs_), _inputFields(inputFields_) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type cell) const {
        const size_type nbfs = _inoutFunction.dimension(1);
        const size_type npts = _inoutFunction.dimension(2);
        const size_type iend = _inoutFunction.dimension(3);
        const size_type jend = _inoutFunction.dimension(4);
        
        for (size_type bf=0;bf<nbfs;++bf) 
          for (size_type pt=0;pt<npts;++pt)
            for (size_type i=0;i<iend;++i) 
              for (size_type j=0;j<jend;++j) 
                outPointVals(cell, pt, i, j) += inCoeffs(cell, bf) * inFields(cell, bf, pt, i, j);
      }
    };
    const size_type loopSize = inputFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(outputPointVals, inputCoeffs, inputFields) );
  }

} // end namespace Intrepid2
