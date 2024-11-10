  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecRight<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,-1>{
    matmatProductDataDataTempSpecRight(ArrayOutData& outputData,
				       ArrayInDataLeft& inputDataLeft,
				       ArrayInDataRight& inputDataRight,
				       const char                transpose){
   
      ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Rank<ArrayInDataLeft>::value,Rank<ArrayInDataRight>::value>(outputData, inputDataLeft, inputDataRight,transpose);		
								
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecRight<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,3>{
    matmatProductDataDataTempSpecRight(ArrayOutData& outputData,
				       ArrayInDataLeft& inputDataLeft,
				       ArrayInDataRight& inputDataRight,
				       const char                transpose){								
      ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Rank<ArrayInDataLeft>::value,3>(outputData, inputDataLeft, inputDataRight,transpose);		
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecRight<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,4>{
    matmatProductDataDataTempSpecRight(ArrayOutData& outputData,
				       ArrayInDataLeft& inputDataLeft,
				       ArrayInDataRight& inputDataRight,
				       const char                transpose){
									
      ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Rank<ArrayInDataLeft>::value,4>(outputData, inputDataLeft, inputDataRight,transpose);										
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,-1,-1>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataLeftRank   = inputDataLeft.rank();
    
      int dataRightRank  = inputDataRight.rank();   							
      switch(dataLeftRank){
      case 2:
	if(dataRightRank==3){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 2,3>(outputData, inputDataLeft, inputDataRight,transpose);										
	
	}else if(dataRightRank==4){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 2,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	}
	break;
      case 3:
	if(dataRightRank==3){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 3,3>(outputData, inputDataLeft, inputDataRight,transpose);										
	
	}else if(dataRightRank==4){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 3,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	}
	break;		
      case 4:
	if(dataRightRank==3){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 4,3>(outputData, inputDataLeft, inputDataRight,transpose);										

	}else if(dataRightRank==4){
	  ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 4,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	}
	break;		
		
      }								
    }
								
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,2,-1>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
   
   
      int dataRightRank  = inputDataRight.rank();    
      switch(dataRightRank){
      case 3:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 2,3>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;
		
      case 4:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 2,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;					
									
									
      }							
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,2,3>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);
      if(numDataLeftPts != 1){			
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, point)*inputDataRight(point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop 						
      }else{
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, 0)*inputDataRight(point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop
		
      }									
									
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,2,4>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);

      if(numDataLeftPts != 1){
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, point)*inputDataRight(cell, point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop							
      }else{
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, 0)*inputDataRight(cell, point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop
      }							
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,3,-1>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataRightRank  = inputDataRight.rank();    
      switch(dataRightRank){
      case 3:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 3,3>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;
		
      case 4:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 3,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;					
									
									
      }						
									
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,3,3>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
									
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);
      if(numDataLeftPts != 1){			
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, point, row)*inputDataRight(point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop		      							
      }else{
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, 0, row)*inputDataRight(point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop		           
				
      }						
																		
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,3,4>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
									
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);
      if(numDataLeftPts != 1){			
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, point, row)*inputDataRight(cell, point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop						
      }else{
	for(int cell = 0; cell < numCells; cell++) {
	  for(int point = 0; point < numPoints; point++) {
	    for( int row = 0; row < matDim; row++) {
	      for( int col = 0; col < matDim; col++) {
		outputData(cell, point, row, col) = \
		  inputDataLeft(cell, 0, row)*inputDataRight(cell, point, row, col);
	      }// Col-loop
	    } // Row-loop
	  } // P-loop
	}// C-loop
      }						
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,4,-1>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){

      int dataRightRank  = inputDataRight.rank();    
      switch(dataRightRank){
      case 3:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 4,3>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;
		
      case 4:
	ArrayTools::matmatProductDataDataTempSpecLeft<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, 4,4>(outputData, inputDataLeft, inputDataRight,transpose);										

	break;					
									
									
      }							
									
									
									
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,4,3>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);									
      if(numDataLeftPts != 1){
	if ((transpose == 'n') || (transpose == 'N')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, point, row, i)*inputDataRight(point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} // no transpose
	else if ((transpose == 't') || (transpose == 'T')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, point, i, row)*inputDataRight(point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} //transpose								
      }else{
	if ((transpose == 'n') || (transpose == 'N')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, 0, row, i)*inputDataRight(point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} // no transpose
	else if ((transpose == 't') || (transpose == 'T')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, 0, i, row)*inputDataRight(point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell	
			
	}						
      }							
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  struct ArrayTools::matmatProductDataDataTempSpecLeft<Scalar, ArrayOutData, ArrayInDataLeft, ArrayInDataRight,4,4>{
    matmatProductDataDataTempSpecLeft(ArrayOutData& outputData,
				      ArrayInDataLeft& inputDataLeft,
				      ArrayInDataRight& inputDataRight,
				      const char                transpose){
      int dataLeftRank   = inputDataLeft.rank();
      int numDataLeftPts = inputDataLeft.dimension(1);
      int dataRightRank  = inputDataRight.rank();    
      int numCells       = outputData.dimension(0);
      int numPoints      = outputData.dimension(1);
      int matDim         = outputData.dimension(2);
      if(numDataLeftPts != 1){				
	if ((transpose == 'n') || (transpose == 'N')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, point, row, i)*inputDataRight(cell, point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} // no transpose
	else if ((transpose == 't') || (transpose == 'T')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, point, i, row)*inputDataRight(cell, point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} //transpose 							
      }else{
	if ((transpose == 'n') || (transpose == 'N')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, 0, row, i)*inputDataRight(cell, point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} // no transpose
	else if ((transpose == 't') || (transpose == 'T')) {
	  for(int cell = 0; cell < numCells; cell++){
	    for(int point = 0; point < numPoints; point++){
	      for(int row = 0; row < matDim; row++){
		for(int col = 0; col < matDim; col++){
		  outputData(cell, point, row, col) = 0.0;
		  for(int i = 0; i < matDim; i++){
		    outputData(cell, point, row, col) += \
		      inputDataLeft(cell, 0, i, row)*inputDataRight(cell, point, i, col);
		  }// i
		} // col
	      } //row
	    }// point
	  }// cell
	} //transpose	
									
      }						
    }
  };
  template<class Scalar, class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight>
  void ArrayTools::matmatProductDataDataTemp(ArrayOutData &            outputData,
					     const ArrayInDataLeft &   inputDataLeft,
					     const ArrayInDataRight &  inputDataRight,
					     const char                transpose){
    ArrayTools::matmatProductDataDataTempSpecRight<Scalar,ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Rank<ArrayInDataRight>::value>(outputData, inputDataLeft, inputDataRight,transpose);		


 								  
  }
  

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

