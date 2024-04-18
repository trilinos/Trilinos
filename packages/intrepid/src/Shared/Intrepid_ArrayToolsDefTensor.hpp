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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_ArrayToolsDefTensor.hpp
    \brief  Definition file for tensor multiply operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {
  
  template<class Scalar, 
           class ArrayOutFields, 
           class ArrayInData, 
           class ArrayInFields>
  void ArrayTools::crossProductDataField(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields){
#ifdef HAVE_INTREPID_DEBUG
    std::string errmsg = ">>> ERROR (ArrayTools::crossProductDataField):";
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
ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);    
ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);    
ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value, true>inputFieldsWrap(inputFields);    

    // 3D cross product
    if(inputData.dimension(2) == 3) {
      
      // inputFields is (C,F,P,D)
      if(getrank(inputFields) == 4){
        
        for(size_t cell = 0; cell < static_cast<size_t>(outputFields.dimension(0)); cell++){
          for(size_t field = 0; field < static_cast<size_t>(outputFields.dimension(1)); field++){
            for(size_t point = 0; point < static_cast<size_t>(outputFields.dimension(2)); point++){
              // 
              outputFieldsWrap(cell, field, point, 0) = \
                inputDataWrap(cell, point, 1)*inputFieldsWrap(cell, field, point, 2) - 
                inputDataWrap(cell, point, 2)*inputFieldsWrap(cell, field, point, 1); 
              // 
              outputFieldsWrap(cell, field, point, 1) = \
                inputDataWrap(cell, point, 2)*inputFieldsWrap(cell, field, point, 0) - 
                inputDataWrap(cell, point, 0)*inputFieldsWrap(cell, field, point, 2); 
              // 
              outputFieldsWrap(cell, field, point, 2) = \
                inputDataWrap(cell, point, 0)*inputFieldsWrap(cell, field, point, 1) - 
                inputDataWrap(cell, point, 1)*inputFieldsWrap(cell, field, point, 0); 
            }// point
          }// field
        } // cell
      }// rank = 4
      // inputFields is (F,P,D)
      else if(getrank(inputFields) == 3){
        
        for(size_t cell = 0; cell < static_cast<size_t>(outputFields.dimension(0)); cell++){
          for(size_t field = 0; field < static_cast<size_t>(outputFields.dimension(1)); field++){
            for(size_t point = 0; point < static_cast<size_t>(outputFields.dimension(2)); point++){
              // 
              outputFieldsWrap(cell, field, point, 0) = \
		inputDataWrap(cell, point, 1)*inputFieldsWrap(field, point, 2) - 
		inputDataWrap(cell, point, 2)*inputFieldsWrap(field, point, 1); 
              // 
              outputFieldsWrap(cell, field, point, 1) = \
                inputDataWrap(cell, point, 2)*inputFieldsWrap(field, point, 0) - 
                inputDataWrap(cell, point, 0)*inputFieldsWrap(field, point, 2); 
              // 
              outputFieldsWrap(cell, field, point, 2) = \
                inputDataWrap(cell, point, 0)*inputFieldsWrap(field, point, 1) - 
                inputDataWrap(cell, point, 1)*inputFieldsWrap(field, point, 0); 
            }// point
          }// field
        } // cell
      }// rank = 3
      else{
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				    ">>> ERROR (ArrayTools::crossProductDataField): inputFields rank 3 or 4 required.")
	  }    
    }
    // 2D cross product
    else if(inputData.dimension(2) == 2){
      
      // inputFields is (C,F,P,D)
      if(getrank(inputFields) == 4){
        
        for(size_t cell = 0; cell < static_cast<size_t>(outputFields.dimension(0)); cell++){
          for(size_t field = 0; field < static_cast<size_t>(outputFields.dimension(1)); field++){
            for(size_t point = 0; point < static_cast<size_t>(outputFields.dimension(2)); point++){
              outputFieldsWrap(cell, field, point) = \
                inputDataWrap(cell, point, 0)*inputFieldsWrap(cell, field, point, 1) - 
                inputDataWrap(cell, point, 1)*inputFieldsWrap(cell, field, point, 0); 
            }// point
          }// field
        } // cell
      }// rank = 4
      // inputFields is (F,P,D)
      else if(getrank(inputFields) == 3) {
        
        for(size_t cell = 0; cell < static_cast<size_t>(outputFields.dimension(0)); cell++){
          for(size_t field = 0; field < static_cast<size_t>(outputFields.dimension(1)); field++){
            for(size_t point = 0; point < static_cast<size_t>(outputFields.dimension(2)); point++){
              outputFieldsWrap(cell, field, point) = \
                inputDataWrap(cell, point, 0)*inputFieldsWrap(field, point, 1) - 
                inputDataWrap(cell, point, 1)*inputFieldsWrap(field, point, 0); 
            }// point
          }// field
        } // cell
      }// rank = 3
    }
    // Error: wrong dimension
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::crossProductDataField): spatial dimension 2 or 3 required.")
	}
  }
  
  
  
  template<class Scalar, 
           class ArrayOutData, 
           class ArrayInDataLeft, 
           class ArrayInDataRight>
  void ArrayTools::crossProductDataData(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight){


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


ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);    
ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);    
ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value, true>inputDataRightWrap(inputDataRight); 
    // 3D cross product
    if(inputDataLeft.dimension(2) == 3) {
      
      // inputDataRight is (C,P,D)
      if(getrank(inputDataRight) == 3){
        
        for(size_t cell = 0; cell < (size_t)inputDataLeft.dimension(0); cell++){
          for(size_t point = 0; point < (size_t)inputDataLeft.dimension(1); point++){
            // 
            outputDataWrap(cell, point, 0) = \
	      inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(cell, point, 2) - 
	      inputDataLeftWrap(cell, point, 2)*inputDataRightWrap(cell, point, 1); 
            // 
            outputDataWrap(cell, point, 1) = \
              inputDataLeftWrap(cell, point, 2)*inputDataRightWrap(cell, point, 0) - 
              inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(cell, point, 2); 
            // 
            outputDataWrap(cell, point, 2) = \
              inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(cell, point, 1) - 
              inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(cell, point, 0); 
          }// point
        } // cell
      }// rank = 3
       // inputDataRight is (P,D)
      else if(getrank(inputDataRight) == 2){
        
        for(size_t cell = 0; cell < (size_t)inputDataLeft.dimension(0); cell++){
          for(size_t point = 0; point < (size_t)inputDataLeft.dimension(1); point++){
            // 
            outputDataWrap(cell, point, 0) = \
	      inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(point, 2) - 
	      inputDataLeftWrap(cell, point, 2)*inputDataRightWrap(point, 1); 
            // 
            outputDataWrap(cell, point, 1) = \
              inputDataLeftWrap(cell, point, 2)*inputDataRightWrap(point, 0) - 
              inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(point, 2); 
            // 
            outputDataWrap(cell, point, 2) = \
              inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(point, 1) - 
              inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(point, 0); 
          }// point
        } // cell
      }// rank = 2
      else{
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				    ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight rank 2 or 3 required.")
	  }    
    }
    // 2D cross product
    else if(inputDataLeft.dimension(2) == 2){
      
      // inputDataRight is (C,P,D)
      if(getrank(inputDataRight) == 3){
        
        for(size_t cell = 0; cell < (size_t)inputDataLeft.dimension(0); cell++){
	  for(size_t point = 0; point < (size_t)inputDataLeft.dimension(1); point++){
	    outputDataWrap(cell, point) = \
	      inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(cell, point, 1) - 
	      inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(cell, point, 0); 
	  }// point
        } // cell
      }// rank = 3
       // inputDataRight is (P,D)
      else if(getrank(inputDataRight) == 2) {
        
        for(size_t cell = 0; cell < (size_t)inputDataLeft.dimension(0); cell++){
	  for(size_t point = 0; point < (size_t)inputDataLeft.dimension(1); point++){
	    outputDataWrap(cell, point) = \
	      inputDataLeftWrap(cell, point, 0)*inputDataRightWrap(point, 1) - 
	      inputDataLeftWrap(cell, point, 1)*inputDataRightWrap(point, 0); 
	  }// point
        } // cell
      }// rank = 2
    }
    // Error: wrong dimension
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::crossProductDataData): spatial dimension 2 or 3 required.")
	}
  }
  
  
  
  template<class Scalar, 
           class ArrayOutFields, 
           class ArrayInData, 
           class ArrayInFields>
  void ArrayTools::outerProductDataField(ArrayOutFields &       outputFields,
                                         const ArrayInData &    inputData,
                                         const ArrayInFields &  inputFields){

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
ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);    
ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);    
ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value, true>inputFieldsWrap(inputFields);  
    
    // inputFields is (C,F,P,D)
 if(getrank(inputFields) == 4){
      
      for(size_t cell = 0; cell < (size_t)outputFields.dimension(0); cell++){
        for(size_t field = 0; field < (size_t)outputFields.dimension(1); field++){
          for(size_t point = 0; point < (size_t)outputFields.dimension(2); point++){
            for(size_t row = 0; row < (size_t)outputFields.dimension(3); row++){
              for(size_t col = 0; col < (size_t)outputFields.dimension(4); col++){
                outputFieldsWrap(cell, field, point, row, col) = \
                  inputDataWrap(cell, point, row)*inputFieldsWrap(cell, field, point, col);
              }// col
            }// row
          }// point
        }// field
      } // cell
    }// rank = 4
     // inputFields is (F,P,D)
 else if(getrank(inputFields) == 3){
      
      for(size_t cell = 0; cell < (size_t)outputFields.dimension(0); cell++){
        for(size_t field = 0; field < (size_t)outputFields.dimension(1); field++){
          for(size_t point = 0; point < (size_t)outputFields.dimension(2); point++){
            for(size_t row = 0; row < (size_t)outputFields.dimension(3); row++){
              for(size_t col = 0; col < (size_t)outputFields.dimension(4); col++){
                outputFieldsWrap(cell, field, point, row, col) = \
                  inputDataWrap(cell, point, row)*inputFieldsWrap(field, point, col);
              }// col
            }// row
          }// point
        }// field
      } // cell
    }// rank = 3
    else{
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::outerProductDataField): inputFields rank 3 or 4 required.")
	}    
  }
  
  
  
  template<class Scalar, 
           class ArrayOutData, 
           class ArrayInDataLeft, 
           class ArrayInDataRight>
  void ArrayTools::outerProductDataData(ArrayOutData &            outputData,
                                        const ArrayInDataLeft &   inputDataLeft,
                                        const ArrayInDataRight &  inputDataRight){

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

ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);    
ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);    
ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value, true>inputDataRightWrap(inputDataRight); 

    // inputDataRight is (C,P,D)
 if(getrank(inputDataRight) == 3){
      
      for(size_t cell = 0; cell <(size_t) inputDataLeft.dimension(0); cell++){
        for(size_t point = 0; point <(size_t) inputDataLeft.dimension(1); point++){
          for(size_t row = 0; row <(size_t) inputDataLeft.dimension(2); row++){
            for(size_t col = 0; col <(size_t) inputDataLeft.dimension(2); col++){
              
              outputDataWrap(cell, point, row, col) = \
                inputDataLeftWrap(cell, point, row)*inputDataRightWrap(cell, point, col); 
            }// col
          }// row
        }// point
      } // cell
    }// rank = 3
     // inputDataRight is (P,D)
 else if(getrank(inputDataRight) == 2){
      
      for(size_t cell = 0; cell <(size_t) inputDataLeft.dimension(0); cell++){
        for(size_t point = 0; point <(size_t) inputDataLeft.dimension(1); point++){
          for(size_t row = 0; row <(size_t) inputDataLeft.dimension(2); row++){
            for(size_t col = 0; col <(size_t) inputDataLeft.dimension(2); col++){
              // 
              outputDataWrap(cell, point, row, col) = \
		inputDataLeftWrap(cell, point, row)*inputDataRightWrap(point, col); 
            } // col
          } // row
        } // point
      } // cell
    }// rank = 2
    else{
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::crossProductDataData): inputDataRight rank 2 or 3 required.")
	}    
  }
  template<class Scalar, 
           class ArrayOutFields, 
           class ArrayInData, 
           class ArrayInFields>
  void ArrayTools::matvecProductDataField(ArrayOutFields &       outputFields,
                                          const ArrayInData &    inputData,
                                          const ArrayInFields &  inputFields,
                                          const char              transpose){

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

   ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);
   ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);
   ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value,true>inputFieldsWrap(inputFields);	
    size_t dataRank   = getrank(inputData);
    size_t numDataPts = inputData.dimension(1);
    size_t inRank     = getrank(inputFields);    
    size_t numCells   = outputFields.dimension(0);
    size_t numFields  = outputFields.dimension(1);
    size_t numPoints  = outputFields.dimension(2);
    size_t matDim     = outputFields.dimension(3);
    /*********************************************************************************************
     *                              inputFields is (C,F,P,D)                                     *
     *********************************************************************************************/
    if(inRank == 4){
      if(numDataPts != 1){  // non-constant data
        
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for( size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, point)*inputFieldsWrap(cell, field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, point, row)*inputFieldsWrap(cell, field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, point, row, col)*inputFieldsWrap(cell, field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, point, col, row)*inputFieldsWrap(cell, field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      }
      else{  // constant data case
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, 0)*inputFieldsWrap(cell, field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, 0, row)*inputFieldsWrap(cell, field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, 0, row, col)*inputFieldsWrap(cell, field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, 0, col, row)*inputFieldsWrap(cell, field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      } // end constant data case
    } // inputFields rank 4
    /*********************************************************************************************
     *                              inputFields is (F,P,D)                                       *
     *********************************************************************************************/
    else if(inRank == 3) {
      if(numDataPts != 1){  // non-constant data
        
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, point)*inputFieldsWrap(field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
		    inputDataWrap(cell, point, row)*inputFieldsWrap(field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, point, row, col)*inputFieldsWrap(field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
			inputDataWrap(cell, point, col, row)*inputFieldsWrap(field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      }
      else{  // constant data case
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
                    inputDataWrap(cell, 0)*inputFieldsWrap(field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  outputFieldsWrap(cell, field, point, row) = \
                    inputDataWrap(cell, 0, row)*inputFieldsWrap(field, point, row);
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
                        inputDataWrap(cell, 0, row, col)*inputFieldsWrap(field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    outputFieldsWrap(cell, field, point, row) = 0.0;
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row) += \
                        inputDataWrap(cell, 0, col, row)*inputFieldsWrap(field, point, col);
		    }// col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      } // end constant data case
    } // inputFields rank 3
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::matvecProductDataField): inputFields rank 3 or 4 required.")      
	}// rank error
  } 
  
  
 
   template<class Scalar, 
           class ArrayOutData, 
           class ArrayInDataLeft, 
           class ArrayInDataRight>
  void ArrayTools::matvecProductDataData(ArrayOutData &            outputData,
                                         const ArrayInDataLeft &   inputDataLeft,
                                         const ArrayInDataRight &  inputDataRight,
                                         const char                transpose){

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

 ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);
 ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);
 ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value, true>inputDataRightWrap(inputDataRight);
    size_t dataLeftRank   = getrank(inputDataLeft);
    size_t numDataLeftPts = inputDataLeft.dimension(1);
    size_t dataRightRank  = getrank(inputDataRight);    
    size_t numCells       = outputData.dimension(0);
    size_t numPoints      = outputData.dimension(1);
    size_t matDim         = outputData.dimension(2);
    
    /*********************************************************************************************
     *                              inputDataRight is (C,P,D)                                   *
     *********************************************************************************************/
    if(dataRightRank == 3){
      if(numDataLeftPts != 1){  // non-constant left data
        
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, point)*inputDataRightWrap(cell, point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, point, row)*inputDataRightWrap(cell, point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, point, row, col)*inputDataRightWrap(cell, point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, point, col, row)*inputDataRightWrap(cell, point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      }
      else{  // constant data case
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, 0)*inputDataRightWrap(cell, point, row);
	      } // Row-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, 0, row)*inputDataRightWrap(cell, point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, 0, row, col)*inputDataRightWrap(cell, point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, 0, col, row)*inputDataRightWrap(cell, point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      } // end constant data case
    } // inputDataRight rank 4
    /*********************************************************************************************
     *                              inputDataRight is (P,D)                                     *
     *********************************************************************************************/
    else if(dataRightRank == 2) {
      if(numDataLeftPts != 1){  // non-constant data
        
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, point)*inputDataRightWrap(point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, point, row)*inputDataRightWrap(point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, point, row, col)*inputDataRightWrap(point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, point, col, row)*inputDataRightWrap(point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      }
      else{  // constant data case
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, 0)*inputDataRightWrap(point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		outputDataWrap(cell, point, row) = \
                  inputDataLeftWrap(cell, 0, row)*inputDataRightWrap(point, row);
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, 0, row, col)*inputDataRightWrap(point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  outputDataWrap(cell, point, row) = 0.0;
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row) += \
                      inputDataLeftWrap(cell, 0, col, row)*inputDataRightWrap(point, col);
		  }// col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      } // end constant inputDataLeft case
    } // inputDataRight rank 2
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight rank 2 or 3 required.")      
	}// rank error
  } 
  
 
  
  template<class Scalar, 
           class ArrayOutFields, 
           class ArrayInData, 
           class ArrayInFields>
  void ArrayTools::matmatProductDataField(ArrayOutFields &       outputFields,
                                          const ArrayInData &    inputData,
                                          const ArrayInFields &  inputFields,
                                          const char             transpose){
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
    ArrayWrapper<Scalar,ArrayOutFields, Rank<ArrayOutFields >::value, false>outputFieldsWrap(outputFields);
    ArrayWrapper<Scalar,ArrayInData, Rank<ArrayInData >::value, true>inputDataWrap(inputData);
    ArrayWrapper<Scalar,ArrayInFields, Rank<ArrayInFields >::value, true>inputFieldsWrap(inputFields);


    size_t dataRank   = getrank(inputData);
    size_t numDataPts = inputData.dimension(1);
    size_t inRank     = getrank(inputFields);    
    size_t numCells   = outputFields.dimension(0);
    size_t numFields  = outputFields.dimension(1);
    size_t numPoints  = outputFields.dimension(2);
    size_t matDim     = outputFields.dimension(3);
    
    /*********************************************************************************************
     *                              inputFields is (C,F,P,D,D)                                     *
     *********************************************************************************************/
    if(inRank == 5){
      if(numDataPts != 1){  // non-constant data
        
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, point)*inputFieldsWrap(cell, field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, point, row)*inputFieldsWrap(cell, field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, point, row, i)*inputFieldsWrap(cell, field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, point, i, row)*inputFieldsWrap(cell, field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      }
      else{  // constant data case
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, 0)*inputFieldsWrap(cell, field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, 0, row)*inputFieldsWrap(cell, field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, 0, row, i)*inputFieldsWrap(cell, field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, 0, i, row)*inputFieldsWrap(cell, field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      } // end constant data case
    } // inputFields rank 5
    /**********************************************************************************************
     *                              inputFields is (F,P,D,D)                                     *
     *********************************************************************************************/
    else if(inRank == 4) {
      if(numDataPts != 1){  // non-constant data
        
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, point)*inputFieldsWrap(field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, point, row)*inputFieldsWrap(field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, point, row, i)*inputFieldsWrap(field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, point, i, row)*inputFieldsWrap(field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      }
      else{  // constant data case
        switch(dataRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, 0)*inputFieldsWrap(field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t field = 0; field < numFields; field++) {
	      for(size_t point = 0; point < numPoints; point++) {
		for(size_t row = 0; row < matDim; row++) {
		  for(size_t col = 0; col < matDim; col++) {
		    outputFieldsWrap(cell, field, point, row, col) = \
                      inputDataWrap(cell, 0, row)*inputFieldsWrap(field, point, row, col);
		  }// Col-loop
		} // Row-loop
	      } // P-loop
	    } // F-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, 0, row, i)*inputFieldsWrap(field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t field = 0; field < numFields; field++){
		for(size_t point = 0; point < numPoints; point++){
		  for(size_t row = 0; row < matDim; row++){
		    for(size_t col = 0; col < matDim; col++){
		      outputFieldsWrap(cell, field, point, row, col) = 0.0;
		      for(size_t i = 0; i < matDim; i++){
			outputFieldsWrap(cell, field, point, row, col) += \
                          inputDataWrap(cell, 0, i, row)*inputFieldsWrap(field, point, i, col);
		      }// i
		    } // col
		  } //row
		}// point
	      }// field
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataField): inputData rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      } // end constant data case
    } // inputFields rank 4
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::matmatProductDataField): inputFields rank 4 or 5 required.")      
	}// rank error
  }


   template<class Scalar, 
           class ArrayOutData, 
           class ArrayInDataLeft, 
           class ArrayInDataRight>
  void ArrayTools::matmatProductDataData(ArrayOutData &            outputData,
                                         const ArrayInDataLeft &   inputDataLeft,
                                         const ArrayInDataRight &  inputDataRight,
                                         const char                transpose){

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


 ArrayWrapper<Scalar,ArrayOutData, Rank<ArrayOutData >::value, false>outputDataWrap(outputData);
 ArrayWrapper<Scalar,ArrayInDataLeft, Rank<ArrayInDataLeft >::value, true>inputDataLeftWrap(inputDataLeft);
 ArrayWrapper<Scalar,ArrayInDataRight, Rank<ArrayInDataRight >::value, true>inputDataRightWrap(inputDataRight);

    size_t dataLeftRank   = getrank(inputDataLeft);
    size_t numDataLeftPts = inputDataLeft.dimension(1);
    size_t dataRightRank  = getrank(inputDataRight);    
    size_t numCells       = outputData.dimension(0);
    size_t numPoints      = outputData.dimension(1);
    size_t matDim         = outputData.dimension(2);
    
    /*********************************************************************************************
     *                              inputDataRight is (C,P,D,D)                                 *
     *********************************************************************************************/
    if(dataRightRank == 4){
      if(numDataLeftPts != 1){  // non-constant data
        
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, point)*inputDataRightWrap(cell, point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, point, row)*inputDataRightWrap(cell, point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, point, row, i)*inputDataRightWrap(cell, point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, point, i, row)*inputDataRightWrap(cell, point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputData rank
      }
      else{  // constant data case
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, 0)*inputDataRightWrap(cell, point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, 0, row)*inputDataRightWrap(cell, point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, 0, row, i)*inputDataRightWrap(cell, point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, 0, i, row)*inputDataRightWrap(cell, point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      } // end constant data case
    } // inputDataRight rank 4
    /**********************************************************************************************
     *                              inputDataRight is (P,D,D)                                    *
     *********************************************************************************************/
    else if(dataRightRank == 3) {
      if(numDataLeftPts != 1){  // non-constant data
        
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, point)*inputDataRightWrap(point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, point, row)*inputDataRightWrap(point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, point, row, i)*inputDataRightWrap(point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, point, i, row)*inputDataRightWrap(point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      }
      else{  // constant data case
        switch(dataLeftRank){
	case 2:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, 0)*inputDataRightWrap(point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 3:
	  for(size_t cell = 0; cell < numCells; cell++) {
	    for(size_t point = 0; point < numPoints; point++) {
	      for(size_t row = 0; row < matDim; row++) {
		for(size_t col = 0; col < matDim; col++) {
		  outputDataWrap(cell, point, row, col) = \
		    inputDataLeftWrap(cell, 0, row)*inputDataRightWrap(point, row, col);
		}// Col-loop
	      } // Row-loop
	    } // P-loop
	  }// C-loop
	  break;
            
	case 4:
	  if ((transpose == 'n') || (transpose == 'N')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, 0, row, i)*inputDataRightWrap(point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } // no transpose
	  else if ((transpose == 't') || (transpose == 'T')) {
	    for(size_t cell = 0; cell < numCells; cell++){
	      for(size_t point = 0; point < numPoints; point++){
		for(size_t row = 0; row < matDim; row++){
		  for(size_t col = 0; col < matDim; col++){
		    outputDataWrap(cell, point, row, col) = 0.0;
		    for(size_t i = 0; i < matDim; i++){
		      outputDataWrap(cell, point, row, col) += \
			inputDataLeftWrap(cell, 0, i, row)*inputDataRightWrap(point, i, col);
		    }// i
		  } // col
		} //row
	      }// point
	    }// cell
	  } //transpose
	  else {
	    TEUCHOS_TEST_FOR_EXCEPTION( !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ), std::invalid_argument,
					">>> ERROR (ArrayTools::matmatProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
	  }
	  break;
            
	default:
	  TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ), std::invalid_argument,
				      ">>> ERROR (ArrayTools::matmatProductDataData): inputDataLeft rank 2, 3 or 4 required.")      
	    } // switch inputDataLeft rank
      } // end constant data case
    } // inputDataRight rank 3
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
				  ">>> ERROR (ArrayTools::matmatProductDataData): inputDataRight rank 3 or 4 required.")      
	}// rank error
  }
  
  
  
  
} // end namespace Intrepid






















#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

