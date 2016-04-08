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
  template<typename outValValueType, class ...outValProperties,
           typename inValValueType,  class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HGRADtransformVALUE( /**/  Kokkos::DynRankView<outValValueType,outValProperties> outVals,
                       const Kokkos::DynRankView<inValValueType, inValProperties>  inVals ) {
    ArrayTools<ExecSpaceType>::cloneFields(outVals, inVals);
  }



  template<typename ExecSpaceType>
  template<typename outValValueType,          class ...outValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inValValueType,           class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HGRADtransformGRAD( /**/  Kokkos::DynRankView<outValValueType,         outValProperties...>          outVals,
                      const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                      const Kokkos::DynRankView<inValValueType,          inValProperties...>           inVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outVals, jacobianInverse, inVals, transpose);
  }

  template<typename ExecSpaceType>
  template<typename outValValueType,          class ...outValProperties,
           typename jacobianInverseValueType, class ...jacobianInverseProperties,
           typename inValValueType,           class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HCURLtransformVALUE( /**/  Kokkos::DynRankView<outValValueType,         outValProperties...>          outVals,
                       const Kokkos::DynRankView<jacobianInverseValueType,jacobianInverseProperties...> jacobianInverse,
                       const Kokkos::DynRankView<inValValueType,          inValProperties...>           inVals,
                       const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outVals, jacobianInverse, inVals, transpose);
  }



  template<typename ExecSpaceType>
  template<typename outValValueType,      class ...outValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inValValueType,       class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HCURLtransformCURL( /**/  Kokkos::DynRankView<outValValueType,     outValProperties...>      outVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inValValueType,      inValProperties...>       inVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outVals, jacobian, inVals, transpose);
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outVals, jacobianDet, outVals, true);
  }


  template<typename ExecSpaceType>
  template<typename outValValueType,      class ...outValProperties,
           typename jacobianValueType,    class ...jacobianProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inValValueType,       class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HDIVtransformVALUE( /**/  Kokkos::DynRankView<outValValueType,     outValProperties...>      outVals,
                      const Kokkos::DynRankView<jacobianValueType,   jacobianProperties...>    jacobian,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inValValueType,      inValProperties...>       inVals,
                      const char transpose = 'T' ) {
    ArrayTools<ExecSpaceType>::matvecProductDataField(outVals, jacobian, inVals, transpose);
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outVals, jacobianDet, outVals, true);
  }



  template<typename ExecSpaceType>
  template<typename outValValueType,      class ...outValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inValValueType,       class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HDIVtransformDIV( /**/  Kokkos::DynRankView<outValValueType,     outValProperties...>      outVals,
                    const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                    const Kokkos::DynRankView<inValValueType,      inValProperties...>       inVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outVals, jacobianDet, inVals, true);
  }



  template<typename ExecSpaceType>
  template<typename outValValueType,      class ...outValProperties,
           typename jacobianDetValueType, class ...jacobianDetProperties,
           typename inValValueType,       class ...inValProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  HVOLtransformVALUE( /**/  Kokkos::DynRankView<outValValueType,     outValProperties...>      outVals,
                      const Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...> jacobianDet,
                      const Kokkos::DynRankView<inValValueType,      inValProperties...>       inVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField(outVals, jacobianDet, inVals, true);
  }
  

  template<typename ExecSpaceType>  
  template<typename outputValueValueType, class ...outputValueProperties,
           typename leftValueValueType,   class ...leftValueProperties,
           typename rightValueValueType,  class ...rightValueProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  integrate( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties> outputValues,
             const Kokkos::DynRankView<leftValueValueType,  leftValueProperties>   leftValues,
             const Kokkos::DynRankView<rightValueValueType, rightValueProperties>  rightValues,
             const ECompEngine  compEngine,
             const bool         sumInto ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( leftValues.rank() < 2 ||
                                      leftValues.rank() > 4, dbgInfo,
                                      ">>> ERROR (FunctionSpaceTools::integrate): Left data must have rank 2, 3 or 4.");
      INTREPID2_TEST_FOR_DEBUG_ABORT( outputValues.rank() < 1 ||
                                      outputValues.rank() > 3, dbgInfo,
                                      ">>> ERROR (FunctionSpaceTools::integrate): Output values must have rank 1, 2 or 3.");
      INTREPID2_TEST_FOR_DEBUG_ABORT( !isValidCompEngine(compEngine),
                                      ">>> ERROR (FunctionSpaceTools::integrate): CompEngine is not valid.");

#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return;
#endif
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
      INTREPID2_TEST_FOR_ABORT( outRank < 1 || outRank > 3,
                                ">>> ERROR (FunctionSpaceTools::integrate): outRank must be 1,2, or 3.");
      INTREPID2_TEST_FOR_ABORT( leftRank < 2 || leftRank > 4,
                                ">>> ERROR (FunctionSpaceTools::integrate): leftRank must be 1,2, 3 or 4.");
    }
    }
  }

  //   template <class Scalar,class ArrayOut,class ArrayDet>
  //   struct computeCellMeasure_Abs {
  //     ArrayOut outVals;
  //     ArrayDet inDet;
  //     typedef typename conditional_eSpace<ArrayOut>::execution_space execution_space;
  //     computeCellMeasure_Abs (ArrayOut outVals_, ArrayDet inDet_) :
  //       outVals (outVals_),inDet (inDet_)
  //     {}


  //     KOKKOS_INLINE_FUNCTION
  //     void operator () (const index_type cell) const {
  //       if (inDet(cell,0) < 0.0) {
  //         for (index_type point=0; point<static_cast<index_type>(outVals.dimension(1)); point++) {
  //           outVals(cell, point) *= -1.0;
  //         }
  //       }
 
  //     }
  //   };

  //   template<class Scalar, class ArrayOut, class ArrayDet, class ArrayWeights>
  //   inline void FunctionSpaceTools::computeCellMeasure(ArrayOut             & outVals,
  //                                                      const ArrayDet       & inDet,
  //                                                      const ArrayWeights   & inWeights) {
  // #ifdef HAVE_INTREPID2_DEBUG

  //     TEUCHOS_TEST_FOR_EXCEPTION( (inDet.rank() != 2), std::invalid_argument,
  //                                 ">>> ERROR (FunctionSpaceTools::computeCellMeasure): Input determinants container must have rank 2.");

  // #endif

  //     ArrayTools::scalarMultiplyDataData<Scalar>(outVals, inDet, inWeights);
  //     // must use absolute value of inDet, so flip sign where needed

  //     Kokkos::parallel_for (outVals.dimension(0), computeCellMeasure_Abs<Scalar, ArrayWrapper<Scalar,ArrayOut, Rank<ArrayOut >::value, false>,ArrayDet> (outValsWrap,inDet));
 

  //   } 


  //   template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  //   void FunctionSpaceTools::computeFaceMeasure(ArrayOut                   & outVals,
  //                                               const ArrayJac             & inJac,
  //                                               const ArrayWeights         & inWeights,
  //                                               const int                    whichFace,
  //                                               const shards::CellTopology & parentCell) {

  // #ifdef HAVE_INTREPID2_DEBUG

  //     TEUCHOS_TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
  //                                 ">>> ERROR (FunctionSpaceTools::computeFaceMeasure): Input Jacobian container must have rank 4.");

  // #endif

  //     // temporary storage for face normals
  //     FieldContainer<Scalar> faceNormals(inJac.dimension(0), inJac.dimension(1), inJac.dimension(2));

  //     // compute normals
  //     CellTools<Scalar>::getPhysicalFaceNormals(faceNormals, inJac, whichFace, parentCell);

  //     // compute lenghts of normals
  //     RealSpaceTools<Scalar>::vectorNorm(outVals, faceNormals, NORM_TWO);

  //     // multiply with weights
  //     ArrayTools::scalarMultiplyDataData<Scalar>(outVals, outVals, inWeights);
  //   }


  //   template<class Scalar, class ArrayOut, class ArrayJac, class ArrayWeights>
  //   void FunctionSpaceTools::computeEdgeMeasure(ArrayOut                   & outVals,
  //                                               const ArrayJac             & inJac,
  //                                               const ArrayWeights         & inWeights,
  //                                               const int                    whichEdge,
  //                                               const shards::CellTopology & parentCell) {

  // #ifdef HAVE_INTREPID2_DEBUG

  //     TEUCHOS_TEST_FOR_EXCEPTION( (inJac.rank() != 4), std::invalid_argument,
  //                                 ">>> ERROR (FunctionSpaceTools::computeEdgeMeasure): Input Jacobian container must have rank 4.");

  // #endif

  //     // temporary storage for edge tangents
  //     FieldContainer<Scalar> edgeTangents(inJac.dimension(0), inJac.dimension(1), inJac.dimension(2));

  //     // compute normals
  //     CellTools<Scalar>::getPhysicalEdgeTangents(edgeTangents, inJac, whichEdge, parentCell);

  //     // compute lenghts of tangents
  //     RealSpaceTools<Scalar>::vectorNorm(outVals, edgeTangents, NORM_TWO);

  //     // multiply with weights
  //     ArrayTools::scalarMultiplyDataData<Scalar>(outVals, outVals, inWeights);

  //   }

  
  template<typename ExecSpaceType>
  template<typename outValValueType,       class ...outValProperties,
           typename inputMeasureValueType, class ...inputMeasureProperties,
           typename inputValValueType,     class ...inputValProperteis>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::
  multiplyMeasure( /**/  Kokkos::DynRankView<outValValueType,      outValProperties...>       outVals,
                   const Kokkos::DynRankView<inputMeasureValueType,inputMeasureProperties...> inputMeasure,
                   const Kokkos::DynRankView<inputValValueType,    inputValProperteis...>     inputVals ) {
    ArrayTools<ExecSpaceType>::scalarMultiplyDataField( outVals, 
                                                        inputMeasure, 
                                                        inputVals );
  }


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataPropertes,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
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
  KOKKOS_INLINE_FUNCTION
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
  KOKKOS_INLINE_FUNCTION
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
  KOKKOS_INLINE_FUNCTION
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
  KOKKOS_INLINE_FUNCTION
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
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 3) && (outRank != 4) && (outRank != 5)), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataField): Output container must have rank 3, 4 or 5.");
    }
    
  } 


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
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
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 2) && (outRank != 3) && (outRank != 4)), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::vectorMultiplyDataData): Output container must have rank 2, 3 or 4.");
    }    
  } 


  template<typename ExecSpaceType>
  template<typename outputFieldValueType, class ...outputFieldProperties,
           typename inputDataValueType,   class ...inputDataProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
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
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 4) && (outRank != 5)), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
    }
  } 
  


  template<typename ExecSpaceType>
  template<typename outputDataValueType,     class ...outputDataProperties,
           typename inputDataLeftValueType,  class ...inputDataLeftProperties,
           typename inputDataRightValueType, class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
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
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ((outRank != 4) && (outRank != 5)), std::invalid_argument,
                                  ">>> ERROR (FunctionSpaceTools::tensorMultiplyDataField): Output container must have rank 4 or 5.");
    }
  }


  template<typename ExecSpaceType>  
  template<typename inoutOperatorValueType, class ...inoutOperatorProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::  
  applyLeftFieldSigns( /**/  Kokkos::DynRankView<inoutOperatorValueType,inoutOperatorProperties...> inoutOperator,
                       const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.rank() != 3), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input operator container must have rank 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Input field signs container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyLeftFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
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
    const size_type iend = inoutOperator.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
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
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(inoutOperator) != 3), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input operator container must have rank 3.");
    TEUCHOS_TEST_FOR_EXCEPTION( (getrank(fieldSigns) != 2), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Input field signs container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyRightFieldSigns): Zeroth dimensions (number of cells) of the operator and field signs containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutOperator.dimension(2) != fieldSigns.dimension(1) ), std::invalid_argument,
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
    const size_type iend = inoutOperator.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(inoutOperator, fieldSigns) );
  } 
  



  template<typename ExecSpaceType>  
  template<typename inoutFunctionValueType, class ...inoutFunctionProperties,
           typename fieldSignValueType,     class ...fieldSignProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::  
  applyFieldSigns( /**/  Kokkos::DynRankView<inoutFunctionValueType,inoutFunctionProperties...> inoutFunction,
                   const Kokkos::DynRankView<fieldSignValueType,    fieldSignProperties...>     fieldSigns ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    
    TEUCHOS_TEST_FOR_EXCEPTION( ((inoutFunction.rank() < 2) || (inoutFunction.rank() > 5)), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input function container must have rank 2, 3, 4, or 5.");
    TEUCHOS_TEST_FOR_EXCEPTION( (fieldSigns.rank() != 2), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Input field signs container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutFunction.dimension(0) != fieldSigns.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::applyFieldSigns): Zeroth dimensions (number of integration domains) of the function and field signs containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inoutFunction.dimension(1) != fieldSigns.dimension(1) ), std::invalid_argument,
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
    const size_type iend = inoutFunction.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(inoutFunction, fieldSigns) );
  }




  template<typename ExecSpaceType>    
  template<typename outputPointValueType, class ...outputPointProperties,
           typename inputCoeffValueType,  class ...inputCoeffProperties,
           typename inputFieldValueType,  class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  void
  FunctionSpaceTools<ExecSpaceType>::  
  evaluate( /**/  Kokkos::DynRankView<outputPointValueType,outputPointProperties...> outputPointVals,
            const Kokkos::DynRankView<inputCoeffValueType, inputCoeffProperties...>  inputCoeffs,
            const Kokkos::DynRankView<inputFieldValueType, inputFieldProperties...>  inputFields ) {
    
#ifdef HAVE_INTREPID2_DEBUG
    
    TEUCHOS_TEST_FOR_EXCEPTION( ((inFields.rank() < 3) || (inFields.rank() > 5)), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): Input fields container must have rank 3, 4, or 5.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.rank() != 2), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): Input coefficient container must have rank 2.");
    TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.rank() != inFields.rank()-1), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): Output values container must have rank one less than the rank of the input fields container.");
    TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the coefficient and fields input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (inCoeffs.dimension(1) != inFields.dimension(1) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): First dimensions (number of fields) of the coefficient and fields input containers must agree!");
    TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.dimension(0) != inFields.dimension(0) ), std::invalid_argument,
                                ">>> ERROR (FunctionSpaceTools::evaluate): Zeroth dimensions (number of cells) of the input fields container and the output values container must agree!");
    for (int i=1; i<outPointVals.rank(); i++) {
      std::string errmsg  = ">>> ERROR (FunctionSpaceTools::evaluate): Dimensions ";
      errmsg += (char)(48+i);
      errmsg += " and ";
      errmsg += (char)(48+i+1);
      errmsg += " of the output values and input fields containers must agree!";
      TEUCHOS_TEST_FOR_EXCEPTION( (outPointVals.dimension(i) != inFields.dimension(i+1)), std::invalid_argument, errmsg );
    }
    
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
    const size_type iend = inputFields.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(outputPointVals, inputCoeffs, inputFields) );
  }

} // end namespace Intrepid2
