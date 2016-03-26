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

/** \file   Intrepid2_ArrayToolsDefTensor.hpp
    \brief  Definition file for tensor multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_TENSOR_HPP__

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  crossProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {

#ifdef HAVE_INTREPID_DEBUG
    const char errmsg[] = ">>> ERROR (ArrayTools::crossProductDataField):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P,D);    
     *      (2) inputFields(C,F,P,D) or (F,P,D);   
     *      (3) outputFields(C,F,P,D) in 3D, or (C,F,P) in 2D
     */
    // (1) inputData is (C, P, D) and 2 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputData, 3, 3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 2, 2,3), 
				std::invalid_argument, errmsg);
    // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputFields, 3,4), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputFields, getrank(inputFields)-1, 2,3), 
				std::invalid_argument, errmsg);
    // (3) outputFields is (C,F,P,D) in 3D and (C,F,P) in 2D => rank = inputData.dimension(2) + 1
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputFields, inputData.dimension(2)+1, inputData.dimension(2)+1), 
				std::invalid_argument, errmsg); 
    /*
     *   Dimension cross-checks:
     *      (1) inputData    vs. inputFields
     *      (2) outputFields vs. inputData
     *      (3) outputFields vs. inputFields 
     *
     *   Cross-check (1):
     */
    if( getrank(inputFields) == 4) {
      // inputData(C,P,D) vs. inputFields(C,F,P,D): dimensions C, P, D must match 
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputData, 0,1,2,  inputFields, 0,2,3),
				  std::invalid_argument, errmsg);
    }
    else{
      // inputData(C,P,D) vs. inputFields(F,P,D): dimensions P, D must match 
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputData, 1,2,  inputFields, 1,2),
				  std::invalid_argument, errmsg);      
    }
    /* 
     *  Cross-check (2): 
     */
    if(inputData.dimension(2) == 2) {
      //  in 2D: outputFields(C,F,P) vs. inputData(C,P,D): dimensions C,P must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, 0,2,  inputData, 0,1),
				  std::invalid_argument, errmsg);
    }
    else{
      // in 3D: outputFields(C,F,P,D) vs. inputData(C,P,D): dimensions C,P,D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, 0,2,3,  inputData, 0,1,2),
				  std::invalid_argument, errmsg);
    }
    /* 
     *  Cross-check (3): 
     */
    if(inputData.dimension(2) == 2) {
      // In 2D:
      if(getrank(inputFields) == 4){
        //  and rank-4 inputFields: outputFields(C,F,P) vs. inputFields(C,F,P,D): dimensions C,F,P must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, 0,1,2,  inputFields, 0,1,2),
				    std::invalid_argument, errmsg);
      }
      else{
        //  and rank-3 inputFields: outputFields(C,F,P) vs. inputFields(F,P,D): dimensions F,P must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, 1,2,  inputFields, 0,1),
				    std::invalid_argument, errmsg);
      }
    }
    else{
      // In 3D:
      if(getrank(inputFields) == 4){
        //  and rank-4 inputFields: outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C,F,P,D must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields,  inputFields),
				    std::invalid_argument, errmsg);
      }
      else{
        // and rank-3 inputFields: outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F,P,D must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, 1,2,3,  inputFields, 0,1,2),
				    std::invalid_argument, errmsg);
      }
    }
#endif  

    // body
  }
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  crossProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::crossProductDataData):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P,D);    
     *      (2) inputDataRight(C,P,D) or (P,D);   
     *      (3) outputData(C,P,D) in 3D, or (C,P) in 2D
     */
    // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataLeft, 3,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 2, 2,3), 
				std::invalid_argument, errmsg);
    // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataRight, 2,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataRight, getrank(inputDataRight)-1,  2,3), 
				std::invalid_argument, errmsg);
    // (3) outputData is (C,P,D) in 3D and (C,P) in 2D => rank = inputDataLeft.dimension(2)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputData, inputDataLeft.dimension(2), inputDataLeft.dimension(2)), 
				std::invalid_argument, errmsg); 
    /*
     *   Dimension cross-checks:
     *      (1) inputDataLeft vs. inputDataRight
     *      (2) outputData    vs. inputDataLeft
     *      (3) outputData    vs. inputDataRight 
     *
     *   Cross-check (1):
     */
    if( getrank(inputDataRight) == 3) {
      // inputDataLeft(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match 
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputDataLeft, inputDataRight),
				  std::invalid_argument, errmsg);
    }
    // inputDataLeft(C, P,D) vs. inputDataRight(P,D): dimensions P, D must match
    else{
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputDataLeft, 1,2, inputDataRight, 0,1),
				  std::invalid_argument, errmsg);      
    }
    /* 
     *  Cross-check (2): 
     */
    if(inputDataLeft.dimension(2) == 2){
      // in 2D: outputData(C,P) vs. inputDataLeft(C,P,D): dimensions C, P must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputDataLeft, 0,1,  outputData, 0,1),
				  std::invalid_argument, errmsg);
    }
    else{
      // in 3D: outputData(C,P,D) vs. inputDataLeft(C,P,D): all dimensions C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputDataLeft, outputData),
				  std::invalid_argument, errmsg);
    }
    /* 
     *  Cross-check (3): 
     */
    if(inputDataLeft.dimension(2) == 2) {
      // In 2D:
      if(getrank(inputDataRight) == 3){
        //  and rank-3 inputDataRight: outputData(C,P) vs. inputDataRight(C,P,D): dimensions C,P must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData, 0,1,  inputDataRight, 0,1),
				    std::invalid_argument, errmsg);
      }
      else{
        //  and rank-2 inputDataRight: outputData(C,P) vs. inputDataRight(P,D): dimension P must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData, 1,  inputDataRight, 0),
				    std::invalid_argument, errmsg);
      }
    }
    else{
      // In 3D:
      if(getrank(inputDataRight) == 3){
        //  and rank-3 inputDataRight: outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C,P,D must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData,  inputDataRight),
				    std::invalid_argument, errmsg);
      }
      else{
        //  and rank-2 inputDataRight: outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData, 1,2,  inputDataRight, 0,1),
				    std::invalid_argument, errmsg);
      }
    }
#endif  
    // body
    
  }
  
  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  outerProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                         const Kokkos::DynRankView<inputDataProperties...>   inputData,
                         const Kokkos::DynRankView<intputFieldProperties...> inputFields ) {
  
#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::outerProductDataField):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P,D);    
     *      (2) inputFields(C,F,P,D) or (F,P,D);   
     *      (3) outputFields(C,F,P,D,D)
     */
    // (1) inputData is (C, P, D) and 2 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputData,  3,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 2,  2,3), 
				std::invalid_argument, errmsg);
    // (2) inputFields is (C, F, P, D) or (F, P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputFields, 3,4), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputFields, getrank(inputFields)-1,  2,3), 
				std::invalid_argument, errmsg);
    // (3) outputFields is (C,F,P,D,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputFields,  5,5), std::invalid_argument, errmsg);      
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputFields, 3,  2,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputFields, 4,  2,3), 
				std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputData    vs. inputFields
     *      (2) outputFields vs. inputData
     *      (3) outputFields vs. inputFields 
     *
     *   Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P,D): dimensions C, P, D must match
     */
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
						       outputFields, 0,2,3,4,
						       inputData,    0,1,2,2),
				std::invalid_argument, errmsg);    
    /*
     *   Cross-checks (1,3):
     */
    if( getrank(inputFields) == 4) {
      // Cross-check (1): inputData(C,P,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 inputData,    0,1,2, 
							 inputFields,  0,2,3),
				  std::invalid_argument, errmsg);  
      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D): dimensions C, F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputFields, 0,1,2,3,4,
							 inputFields,  0,1,2,3,3),
				  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputData(C,P,D) vs. inputFields(F,P,D): dimensions  P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 inputData,    1,2, 
							 inputFields,  1,2),
				  std::invalid_argument, errmsg);      
      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D): dimensions F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputFields, 1,2,3,4, 
							 inputFields,  0,1,2,2),
				  std::invalid_argument, errmsg);
    }
#endif  
    // body
  }
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  outerProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                        const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                        const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::outerProductDataData):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P,D);    
     *      (2) inputDataRight(C,P,D) or (P,D);   
     *      (3) outputData(C,P,D,D)
     */
    // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataLeft,  3,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 2,  2,3), 
				std::invalid_argument, errmsg);
    // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataRight,  2,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataRight, getrank(inputDataRight)-1,  2,3), 
				std::invalid_argument, errmsg);
    // (3) outputData is (C,P,D,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputData, 4, 4), std::invalid_argument, errmsg);      
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputData, 2,  2,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputData, 3,  2,3), 
				std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputDataLeft vs. inputDataRight
     *      (2) outputData    vs. inputDataLeft
     *      (3) outputData    vs. inputDataRight 
     *
     *   Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P,D): dimensions C, P, D must match
     */
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
						       outputData,    0,1,2,3,
						       inputDataLeft, 0,1,2,2),
				std::invalid_argument, errmsg);    
    /*
     *   Cross-checks (1,3):
     */
    if( getrank(inputDataRight) == 3) {
      // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(C,P,D):  all dimensions  C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, inputDataLeft, inputDataRight),
				  std::invalid_argument, errmsg);  
      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D): dimensions C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputData,     0,1,2,3,
							 inputDataRight, 0,1,2,2),
				  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(P,D): dimensions  P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 inputDataLeft,  1,2, 
							 inputDataRight, 0,1),
				  std::invalid_argument, errmsg);      
      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D): dimensions P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputData,     1,2,3, 
							 inputDataRight, 0,1,1),
				  std::invalid_argument, errmsg);
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputFieldProperties,
           class ...inputDataProperties,
           class ...inputFieldProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>
  matvecProductDataField( /**/  Kokkos::DynRankView<outputFieldProperties...> outputFields,
                          const Kokkos::DynRankView<inputDataProperties...>   inputData,
                          const Kokkos::DynRankView<intputFieldProperties...> inputFields,
                          const char transpose ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::matvecProductDataField):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
     *      (2) inputFields(C,F,P,D) or (F,P,D);   
     *      (3) outputFields(C,F,P,D)
     */
    // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputData,  2,4), 
				std::invalid_argument, errmsg);
    if(getrank(inputData) > 2) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 2,  1,3), 
				  std::invalid_argument, errmsg);
    }
    if(getrank(inputData) == 4) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 3,  1,3), 
				  std::invalid_argument, errmsg);
    }
    // (2) inputFields is (C, F, P, D) or (F, P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputFields, 3,4), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputFields, getrank(inputFields)-1,  1,3), 
				std::invalid_argument, errmsg);
    // (3) outputFields is (C,F,P,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputFields,  4,4), std::invalid_argument, errmsg);      
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputFields, 3,  1,3), 
				std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputData    vs. inputFields
     *      (2) outputFields vs. inputData
     *      (3) outputFields vs. inputFields
     *     
     *   Cross-check (2): outputFields(C,F,P,D) vs. inputData(C,P), (C,P,D) or (C,P,D,D): 
     *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
     *   data is specified (P>1). Do not check P dimensions with constant data, i.e., when P=1 in
     *   inputData(C,1,...)
     */
    if(inputData.dimension(1) > 1){ // check P dimension if P>1 in inputData
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 2,
							 inputData,    1),
				  std::invalid_argument, errmsg);    
    }
    if(getrank(inputData) == 2) { // inputData(C,P) -> C match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,
							 inputData,    0),
				  std::invalid_argument, errmsg);    
    }
    if(getrank(inputData) == 3){ // inputData(C,P,D) -> C, D match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,3,
							 inputData,    0,2),
				  std::invalid_argument, errmsg);    
      
    }
    if(getrank(inputData) == 4){ // inputData(C,P,D,D) -> C, D, D match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,3,3,
							 inputData,    0,2,3),
				  std::invalid_argument, errmsg);
    }
    /*
     *   Cross-checks (1,3):
     */
    if(getrank(inputFields) == 4) {      
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match      
      if(inputData.dimension(1) > 1){ // check P dimension if P>1 in inputData
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							   inputFields,  2,
							   inputData,    1),
				    std::invalid_argument, errmsg);    
      }      
      if(getrank(inputData) == 2){ // inputData(C,P) -> C match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0, 
							   inputFields,  0),
				    std::invalid_argument, errmsg);  
      }
      if(getrank(inputData) == 3){  // inputData(C,P,D) -> C, D match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0,2, 
							   inputFields,  0,3),
				    std::invalid_argument, errmsg);  
      }
      if(getrank(inputData) == 4){   // inputData(C,P,D,D) -> C, D, D match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0,2,3, 
							   inputFields,  0,3,3),
				    std::invalid_argument, errmsg);  
      }
      // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C, F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, inputFields),
				  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D): dimensions  P, D must match
      if(inputData.dimension(1) > 1){ // check P if P>1 in inputData 
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    1, 
							   inputFields,  1),
				    std::invalid_argument, errmsg);    
      }
      if(getrank(inputData) == 3){
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    2, 
							   inputFields,  2),
				    std::invalid_argument, errmsg);    
      }
      if(getrank(inputData) == 4){
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    2,3, 
							   inputFields,  2,2),
				    std::invalid_argument, errmsg);            
      }
      // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputFields, 1,2,3, 
							 inputFields,  0,1,2),
				  std::invalid_argument, errmsg);
    }
#endif
    // body
  } 

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  matvecProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>    outputData,
                         const Kokkos::DynRankView<inputDataLeftProperties...> inputDataLeft,
                         const Kokkos::DynRankView<intputFieldProperties...>   inputDataRight,
                         const char transpose ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::matvecProductDataData):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data   
     *      (2) inputDataRight(C,P,D) or (P,D);   
     *      (3) outputData(C,P,D)
     */
    // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataLeft,  2,4), 
                                std::invalid_argument, errmsg);
    if(getrank(inputDataLeft) > 2) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 2,  1,3), 
                                  std::invalid_argument, errmsg);
    }
    if(getrank(inputDataLeft) == 4) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 3,  1,3), 
                                  std::invalid_argument, errmsg);
    }
    // (2) inputDataRight is (C, P, D) or (P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataRight,  2,3), 
                                std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataRight, getrank(inputDataRight)-1,  1,3), 
                                std::invalid_argument, errmsg);
    // (3) outputData is (C,P,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputData, 3,3), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputData, 2,  1,3), 
                                std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputDataLeft vs. inputDataRight
     *      (2) outputData    vs. inputDataLeft
     *      (3) outputData    vs. inputDataRight 
     *
     *   Cross-check (2): outputData(C,P,D) vs. inputDataLeft(C,P), (C,P,D) or (C,P,D,D):
     *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
     *   data is specified (P>1). Do not check P dimensions with constant left data, i.e., when P=1 in
     *   inputDataLeft(C,1,...)
     */
    if(inputDataLeft.dimension(1) > 1){ // check P dimension if P>1 in inputDataLeft
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,     1,
                                                         inputDataLeft,  1),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 2){  // inputDataLeft(C,P): check C
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,
                                                         inputDataLeft, 0),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check C and D
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,2,
                                                         inputDataLeft, 0,2),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check C and D
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,2,2,
                                                         inputDataLeft, 0,2,3),
                                  std::invalid_argument, errmsg);    
    }
    /*
     *   Cross-checks (1,3):
     */
    if( getrank(inputDataRight) == 3) {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D):  dimensions  C, P, D must match
      if(inputDataLeft.dimension(1) > 1){ // check P dimension if P>1 in inputDataLeft
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                           inputDataLeft,  1,
                                                           inputDataRight, 1),
                                    std::invalid_argument, errmsg);    
      }      
      if(getrank(inputDataLeft) == 2){  // inputDataLeft(C,P): check C
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0, 
                                                           inputDataRight, 0),
                                    std::invalid_argument, errmsg);  
      }      
      if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check C and D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0,2, 
                                                           inputDataRight, 0,2),
                                    std::invalid_argument, errmsg);  
      }      
      if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check C and D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0,2,3, 
                                                           inputDataRight, 0,2,2),
                                    std::invalid_argument, errmsg);  
      }
      
      // Cross-check (3): outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData, inputDataRight),
                                  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D): dimensions  P, D must match
      if(inputDataLeft.dimension(1) > 1){ // check P if P>1 in inputData 
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                           inputDataLeft,  1,
                                                           inputDataRight, 0),
                                    std::invalid_argument, errmsg);    
      }
      if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  2, 
                                                           inputDataRight, 1),
                                    std::invalid_argument, errmsg); 
      }
      if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check D      
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  2,3, 
                                                           inputDataRight, 1,1),
                                    std::invalid_argument, errmsg); 
      }
      // Cross-check (3): outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                         outputData,     1,2, 
                                                         inputDataRight, 0,1),
                                  std::invalid_argument, errmsg);
    }
#endif
    // body
  } 
  
  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>::
  matmatProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                         const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                         const char transpose  ) {

#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::matmatProductDataField):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputData(C,P), (C,P,D) or (C,P,D,D);   P=1 is admissible to allow multiply by const. data
     *      (2) inputFields(C,F,P,D,D) or (F,P,D,D);   
     *      (3) outputFields(C,F,P,D,D)
     */
    // (1) inputData is (C,P), (C, P, D) or (C, P, D, D) and 1 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputData,  2,4), 
				std::invalid_argument, errmsg);
    if(getrank(inputData) > 2) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 2,  1,3), 
				  std::invalid_argument, errmsg);
    }
    if(getrank(inputData) == 4) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputData, 3,  1,3), 
				  std::invalid_argument, errmsg);
    }
    // (2) inputFields is (C,F,P,D,D) or (F,P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputFields, 4,5), std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputFields, getrank(inputFields)-1,  1,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputFields, getrank(inputFields)-2,  1,3), 
				std::invalid_argument, errmsg);
    // (3) outputFields is (C,F,P,D,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputFields,  5,5), std::invalid_argument, errmsg);      
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputFields, 3,  1,3), 
				std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputFields, 4,  1,3), 
				std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputData    vs. inputFields
     *      (2) outputFields vs. inputData
     *      (3) outputFields vs. inputFields 
     *
     *   Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P), (C,P,D) or (C,P,D,D): 
     *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
     *   data is specified (P>1). Do not check P dimensions with constant data, i.e., when P=1 in
     *   inputData(C,1,...)
     */
    if(inputData.dimension(1) > 1){ // check P dimension if P>1 in inputData
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 2,
							 inputData,    1),
				  std::invalid_argument, errmsg);    
    }
    if(getrank(inputData) == 2) { // inputData(C,P) -> C match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,
							 inputData,    0),
				  std::invalid_argument, errmsg);    
    }
    if(getrank(inputData) == 3){ // inputData(C,P,D) -> C, D match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,3,4,
							 inputData,    0,2,2),
				  std::invalid_argument, errmsg);    
      
    }
    if(getrank(inputData) == 4){ // inputData(C,P,D,D) -> C, D, D match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							 outputFields, 0,3,4,
							 inputData,    0,2,3),
				  std::invalid_argument, errmsg);
    }
    /*
     *   Cross-checks (1,3):
     */
    if( getrank(inputFields) == 5) {
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D,D):  dimensions  C, P, D must match
      if(inputData.dimension(1) > 1){ // check P dimension if P>1 in inputData
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
							   inputData,    1,
							   inputFields,  2),
				    std::invalid_argument, errmsg);    
      }      
      if(getrank(inputData) == 2){ // inputData(C,P) -> C match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0, 
							   inputFields,  0),
				    std::invalid_argument, errmsg);  
      }
      if(getrank(inputData) == 3){  // inputData(C,P,D) -> C, D match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0,2,2, 
							   inputFields,  0,3,4),
				    std::invalid_argument, errmsg);  
      }
      if(getrank(inputData) == 4){   // inputData(C,P,D,D) -> C, D, D match
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    0,2,3, 
							   inputFields,  0,3,4),
				    std::invalid_argument, errmsg);  
      }
      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D,D): all dimensions C, F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputFields, inputFields),
				  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D,D): dimensions  P, D must match
      if(inputData.dimension(1) > 1){ // check P if P>1 in inputData 
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    1, 
							   inputFields,  1),
				    std::invalid_argument, errmsg);    
      }
      if(getrank(inputData) == 3){
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    2,2, 
							   inputFields,  2,3),
				    std::invalid_argument, errmsg);    
      }
      if(getrank(inputData) == 4){
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							   inputData,    2,3, 
							   inputFields,  2,3),
				    std::invalid_argument, errmsg);            
      }
      // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D,D): dimensions F, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
							 outputFields, 1,2,3,4, 
							 inputFields,  0,1,2,3),
				  std::invalid_argument, errmsg);
    }
#endif
    // body
  }

  template<typename ExecSpaceType>
  template<class ...outputDataProperties,
           class ...inputDataLeftProperties,
           class ...inputDataRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  ArrayTools<ExecSpaceType>
  matmatProductDataData( /**/  Kokkos::DynRankView<outputDataProperties...>     outputData,
                         const Kokkos::DynRankView<inputDataLeftProperties...>  inputDataLeft,
                         const Kokkos::DynRankView<inputDataRightProperties...> inputDataRight,
                         const char transpose  ) {
    
#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::matmatProductDataData):";
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) inputDataLeft(C,P), (C,P,D) or (C,P,D,D); P=1 is admissible to allow multiply by const. left data   
     *      (2) inputDataRight(C,P,D,D) or (P,D,D);   
     *      (3) outputData(C,P,D,D)
     */
    // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required  
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataLeft,  2,4), 
                                std::invalid_argument, errmsg);
    if(getrank(inputDataLeft) > 2) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 2,  1,3), 
                                  std::invalid_argument, errmsg);
    }
    if(getrank(inputDataLeft) == 4) {
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataLeft, 3,  1,3), 
                                  std::invalid_argument, errmsg);
    }
    // (2) inputDataRight is (C,P,D,D) or (P,D,D) and 1 <= (D=dimension(rank-1),(rank-2)) <= 3 is required. 
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, inputDataRight,  3,4), 
                                std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataRight, getrank(inputDataRight)-1,  1,3), 
                                std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, inputDataRight, getrank(inputDataRight)-2,  1,3), 
                                std::invalid_argument, errmsg);
    // (3) outputData is (C,P,D,D)
    TEUCHOS_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, outputData, 4, 4), std::invalid_argument, errmsg);      
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputData, 2,  1,3), 
                                std::invalid_argument, errmsg);
    TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, outputData, 3,  1,3), 
                                std::invalid_argument, errmsg);
    /*
     *   Dimension cross-checks:
     *      (1) inputDataLeft vs. inputDataRight
     *      (2) outputData    vs. inputDataLeft
     *      (3) outputData    vs. inputDataRight 
     *
     *   Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P), (C,P,D) or (C,P,D,D):
     *   dimensions C, and D must match in all cases, dimension P must match only when non-constant
     *   data is specified (P>1). Do not check P dimensions with constant left data, i.e., when P=1 in
     *   inputDataLeft(C,1,...)
     */
    if(inputDataLeft.dimension(1) > 1){ // check P dimension if P>1 in inputDataLeft
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,     1,
                                                         inputDataLeft,  1),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 2){  // inputDataLeft(C,P): check C
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,
                                                         inputDataLeft, 0),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check C and D
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,2,3,
                                                         inputDataLeft, 0,2,2),
                                  std::invalid_argument, errmsg);    
    }
    if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check C and D
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                         outputData,    0,2,3,
                                                         inputDataLeft, 0,2,3),
                                  std::invalid_argument, errmsg);    
    }
    /*
     *   Cross-checks (1,3):
     */
    if( getrank(inputDataRight) == 4) {
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D,D):  dimensions  C, P, D must match
      if(inputDataLeft.dimension(1) > 1){ // check P dimension if P>1 in inputDataLeft
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                           inputDataLeft,  1,
                                                           inputDataRight, 1),
                                    std::invalid_argument, errmsg);    
      }      
      if(getrank(inputDataLeft) == 2){  // inputDataLeft(C,P): check C
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0, 
                                                           inputDataRight, 0),
                                    std::invalid_argument, errmsg);  
      }      
      if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check C and D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0,2,2, 
                                                           inputDataRight, 0,2,3),
                                    std::invalid_argument, errmsg);  
      }      
      if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check C and D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  0,2,3, 
                                                           inputDataRight, 0,2,3),
                                    std::invalid_argument, errmsg);  
      }
      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D,D): all dimensions C, P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, outputData, inputDataRight),
                                  std::invalid_argument, errmsg);
    }
    else{
      // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D,D): dimensions  P, D must match
      if(inputDataLeft.dimension(1) > 1){ // check P if P>1 in inputData 
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg,
                                                           inputDataLeft,  1,
                                                           inputDataRight, 0),
                                    std::invalid_argument, errmsg);    
      }
      if(getrank(inputDataLeft) == 3){   // inputDataLeft(C,P,D): check D
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  2,2, 
                                                           inputDataRight, 1,2),
                                    std::invalid_argument, errmsg); 
      }
      if(getrank(inputDataLeft) == 4){   // inputDataLeft(C,P,D,D): check D      
        TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                           inputDataLeft,  2,3, 
                                                           inputDataRight, 1,2),
                                    std::invalid_argument, errmsg); 
      }
      // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D,D): dimensions P, D must match
      TEUCHOS_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, 
                                                         outputData,     1,2,3, 
                                                         inputDataRight, 0,1,2),
                                  std::invalid_argument, errmsg);
    }
#endif
    //body
  }
  
} // end namespace Intrepid
#endif
