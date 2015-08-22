#ifndef INTREPID_ARRAYTOOLSDEFSCALAR_KOKKOS_HPP
#define INTREPID_ARRAYTOOLSDEFSCALAR_KOKKOS_HPP

#ifdef HAVE_INTREPID_KOKKOSCORE

#include "Intrepid_ArrayToolsDefScalar.hpp"

namespace Intrepid{


//Functors for parallel_for rank 2_2
//Below you can find the switch staments and paralle fors that call
//The functors.
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const2_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const2_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(i,j)/inputDataLeft(i,j);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const2_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const2_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(i,j)*inputDataLeft(i,j);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const2_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const2_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(i,j)/inputDataLeft(i,0);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const2_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const2_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(i,j)*inputDataLeft(i,0);
     }
  }
};
//End Functors for parallel_for rank 2_2

//Partially specialized implementation of scalarMultiplyDataData with rank 2_2 
	  
	template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,2,2>{
		scalarMultiplyDataData2Kokkos(ArrayOutData & outputData,
                                ArrayInDataLeft & inputDataLeft,
                                ArrayInDataRight & inputDataRight,
                                const bool           reciprocal){
	
#ifdef HAVE_INTREPID_DEBUG
  			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 2)), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != inputDataLeft.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first dimension of the left data input container must be 1!");
    		for (int i=0; i<inputDataRight.rank(); i++) {
      			std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimension ";
      			errmsg += (char)(48+i);
      			errmsg += " of the right input and output data containers must agree!";
      			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    		}
  
 
#endif
  
 // get sizes

	int numCells       = outputData.dimension(0);
	int numDataPoints  = inputDataLeft.dimension(1);
   
     if (numDataPoints != 1) {
          if (reciprocal) {
          	  Kokkos::parallel_for(numCells,PFor__Recip_Not_Const2_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));

          }
         else {
         Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const2_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
              }
		}
		else{
		   if (reciprocal) {
		Kokkos::parallel_for(numCells,PFor__Recip__Const2_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
            Kokkos::parallel_for(numCells,PFor_Not_Recip__Const2_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
              }
            }   	
		
		}
	};
	
//This is a test that uses nested paralleism
//Functors for parallel_for rank 3_3
#if defined( KOKKOS_HAVE_CXX11 )
#define VECTOR_LENGTH 32
typedef Kokkos::TeamVectorPolicy<8> team_policy;
typedef team_policy::member_type team_member ;


template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const3_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const3_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member & thread) const {
	 int i = thread.league_rank(); 
	 thread.team_par_for(outputData.dimension_1(),[&] (const int& j){
		 thread.vector_par_for(outputData.dimension_2(),[&] (const int& k){
			 outputData(i,j,k) = inputDataRight(i,j,k)/inputDataLeft(i,j);
			 });		  	
		 });
   /* for(int j = 0; j < outputData.dimension_1(); j++){
    	for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(i,j,k)/inputDataLeft(i,j);
      	}
     }*/
  }
};
#endif
/*
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const3_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const3_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
	  
    for(int j = 0; j < outputData.dimension_1(); j++){
    	for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(i,j,k)/inputDataLeft(i,j);
      	}
     }
  }
};
*/

template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const3_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const3_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(i,j,k)*inputDataLeft(i,j);
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const3_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const3_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(i,j,k)/inputDataLeft(i,0);
      	}
     }
  }
};

template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const3_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const3_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(i,j,k)*inputDataLeft(i,0);
      	}
     }
  }
};	
//End Functors for parallel_for rank 3_3

//Partially specialized implementation of scalarMultiplyDataData with rank 3_3 	
		
template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
	struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,3,3>{	
			scalarMultiplyDataData2Kokkos(ArrayOutData & outputData,
                                ArrayInDataLeft & inputDataLeft,
                                ArrayInDataRight & inputDataRight,
                                const bool           reciprocal){
		#ifdef HAVE_INTREPID_DEBUG
  			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 3)), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 3.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank of 3.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != inputDataLeft.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first dimension of the left data input container must be 1!");
    		for (int i=0; i<inputDataRight.rank(); i++) {
      			std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimension ";
      			errmsg += (char)(48+i);
      			errmsg += " of the right input and output data containers must agree!";
      			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    		}
  
 
#endif

 // get sizes

	int numCells       = outputData.dimension(0);
	int numDataPoints  = inputDataLeft.dimension(1);
  
      if (numDataPoints != 1) {
 
         if (reciprocal) {
#if defined( KOKKOS_HAVE_CXX11 )
			  const team_policy policy( numCells , 2 );
              Kokkos::parallel_for(policy,PFor__Recip_Not_Const3_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, inputDataRight));
#endif
          }
          else {
    
              Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const3_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		else{
	 
          if (reciprocal) {
      
         	  Kokkos::parallel_for(numCells,PFor__Recip__Const3_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
       
              Kokkos::parallel_for(numCells,PFor_Not_Recip__Const3_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		}
	};	
	
	
//Functors for parallel_for rank 4_4
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const4_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const4_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
    	for(int k = 0; k < outputData.dimension_2(); k++){
    		for(int l = 0; l < outputData.dimension_3(); l++){
      		outputData(i,j,k,l) = inputDataRight(i,j,k,l)/inputDataLeft(i,j);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const4_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const4_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      			outputData(i,j,k,l) = inputDataRight(i,j,k,l)*inputDataLeft(i,j);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const4_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const4_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      				outputData(i,j,k,l) = inputDataRight(i,j,k,l)/inputDataLeft(i,0);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const4_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const4_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      					outputData(i,j,k,l) = inputDataRight(i,j,k,l)/inputDataLeft(i,0);
      		}
      	}
     }
  }
};	
//End Functors for parallel_for rank 4_4

	
//Partially specialized implementation of scalarMultiplyDataData with rank 4_4 		

template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
	struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,4,4>{
	scalarMultiplyDataData2Kokkos(ArrayOutData&  outputData,
                            ArrayInDataLeft&  inputDataLeft,
                            ArrayInDataRight& inputDataRight,
                            const bool           reciprocal){
		#ifdef HAVE_INTREPID_DEBUG
  			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 4)), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 4.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input and output data containers must have the same rank of 4.");
    		TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(0) != inputDataLeft.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions (number of integration domains) of the left and right data input containers must agree!");
    		TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(1) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): First dimensions of the left and right data input containers (number of integration points) must agree or first dimension of the left data input container must be 1!");
    		for (int i=0; i<inputDataRight.rank(); i++) {
      			std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimension ";
      			errmsg += (char)(48+i);
      			errmsg += " of the right input and output data containers must agree!";
      			TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i)), std::invalid_argument, errmsg );
    		}
  
 
#endif 
 // get sizes

	int numCells       = outputData.dimension(0);
	int numDataPoints  = inputDataLeft.dimension(1);

      if (numDataPoints != 1) {
          if (reciprocal) {
              Kokkos::parallel_for(numCells,PFor__Recip_Not_Const4_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const4_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		else{
		
           if (reciprocal) {
      
              Kokkos::parallel_for(numCells,PFor__Recip__Const4_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip__Const4_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		}
		
};
	

//Functors for parallel_for rank 1_2	
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const1_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const1_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(j)/inputDataLeft(i,j);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const1_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const1_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(j)*inputDataLeft(i,j);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const1_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const1_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(j)/inputDataLeft(i,0);
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const1_2 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const1_2(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
      outputData(i,j) = inputDataRight(j)*inputDataLeft(i,0);
     }
  }
};
//End Functors for parallel_for rank 1_2

//Partially specialized implementation of scalarMultiplyDataData with rank 1_2 		
    
     	template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
	struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,1,2>{
		scalarMultiplyDataData2Kokkos(ArrayOutData&  outputData,
                            ArrayInDataLeft&  inputDataLeft,
                            ArrayInDataRight& inputDataRight,
                            const bool           reciprocal){
			#ifdef HAVE_INTREPID_DEBUG
  				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 1.");
    			TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( inputDataLeft.dimension(0) != outputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
    			for (int i=0; i<inputDataRight.rank(); i++) {
      				std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimensions ";
      				errmsg += (char)(48+i);
      				errmsg += " and ";
      				errmsg += (char)(48+i+1);
      				errmsg += " of the right input and output data containers must agree!";
      				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i+1)), std::invalid_argument, errmsg );
    					}
  			#endif
  			
 // get sizes  			
	int numCells       = outputData.dimension(0);
	int numDataPoints  = inputDataLeft.dimension(1);
 
         if (numDataPoints != 1) {
         if (reciprocal) {
              Kokkos::parallel_for(numCells,PFor__Recip_Not_Const1_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const1_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		else{
          if (reciprocal) {
         	  Kokkos::parallel_for(numCells,PFor__Recip__Const1_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip__Const1_2<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
      }
  };
      	
//Functors for parallel_for rank 2_3
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const2_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const2_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
	 
    for(int j = 0; j < outputData.dimension_1(); j++){
    	for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(j,k)/inputDataLeft(i,j);
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const2_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const2_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(j,k)*inputDataLeft(i,j);
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const2_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const2_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(j,k)/inputDataLeft(i,0);
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const2_3 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const2_3(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		outputData(i,j,k) = inputDataRight(j,k)*inputDataLeft(i,0);
      	}
     }
  }
};	
//End Functors for parallel_for rank 2_3

//Partially specialized implementation of scalarMultiplyDataData with rank 2_3

 template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,2,3>{
      		scalarMultiplyDataData2Kokkos(ArrayOutData&  outputData,
                            ArrayInDataLeft&  inputDataLeft,
                            ArrayInDataRight& inputDataRight,
                            const bool           reciprocal){
								
#ifdef HAVE_INTREPID_DEBUG
  				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 2) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 2.");
    			TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( inputDataLeft.dimension(0) != outputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
    			for (int i=0; i<inputDataRight.rank(); i++) {
      				std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimensions ";
      				errmsg += (char)(48+i);
      				errmsg += " and ";
      				errmsg += (char)(48+i+1);
      				errmsg += " of the right input and output data containers must agree!";
      				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i+1)), std::invalid_argument, errmsg );
    			}
#endif
 // get sizes
  int numCells       = outputData.dimension(0);
  int numDataPoints  = inputDataLeft.dimension(1);
      if (numDataPoints != 1) {
         if (reciprocal) {
              Kokkos::parallel_for(numCells,PFor__Recip_Not_Const2_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const2_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		else{
          if (reciprocal) {
         	  Kokkos::parallel_for(numCells,PFor__Recip__Const2_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip__Const2_3<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
 
      }
  };
   
//Functors for parallel_for rank 3_4
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip_Not_Const3_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip_Not_Const3_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
    	for(int k = 0; k < outputData.dimension_2(); k++){
    		for(int l = 0; l < outputData.dimension_3(); l++){
      		outputData(i,j,k,l) = inputDataRight(j,k,l)/inputDataLeft(i,j);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip_Not_Const3_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip_Not_Const3_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      			outputData(i,j,k,l) = inputDataRight(j,k,l)*inputDataLeft(i,j);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor__Recip__Const3_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor__Recip__Const3_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      				outputData(i,j,k,l) = inputDataRight(j,k,l)/inputDataLeft(i,0);
      		}
      	}
     }
  }
};
template<class ViewType, class ViewType1,class ViewType2>
struct PFor_Not_Recip__Const3_4 {
  ViewType outputData;
  typename ViewType1::const_type inputDataLeft;
  typename ViewType2::const_type inputDataRight;
  PFor_Not_Recip__Const3_4(ViewType outputData_, ViewType1 inputDataLeft_, 
              ViewType2 inputDataRight_):outputData(outputData_),inputDataLeft(inputDataLeft_),inputDataRight(inputDataRight_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
        for(int k = 0; k < outputData.dimension_2(); k++){
      		for(int l = 0; l < outputData.dimension_3(); l++){
      					outputData(i,j,k,l) = inputDataRight(j,k,l)/inputDataLeft(i,0);
      		}
      	}
     }
  }
};	
//End Functors for parallel_for rank 3_4     

//Partially specialized implementation of scalarMultiplyDataData with rank 3_4
 
   template<class ArrayOutData, class ArrayInDataLeft, class ArrayInDataRight, class Layout, class MemorySpace>
struct ArrayTools::scalarMultiplyDataData2Kokkos<ArrayOutData, ArrayInDataLeft, ArrayInDataRight, Layout, MemorySpace,3,4>{
   		scalarMultiplyDataData2Kokkos(ArrayOutData&  outputData,
                            ArrayInDataLeft&  inputDataLeft,
                            ArrayInDataRight& inputDataRight,
                            const bool           reciprocal){
#ifdef HAVE_INTREPID_DEBUG
  				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataLeft.rank() != 2), std::invalid_argument,
                      ">>> ERROR (ArrayTools::scalarMultiplyDataData): Left input data container must have rank 2.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.rank() != 3) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Right input data container must have rank 3.");
    			TEUCHOS_TEST_FOR_EXCEPTION( (outputData.rank() != inputDataRight.rank()+1), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): The rank of the right input data container must be one less than the rank of the output data container.");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( (inputDataRight.dimension(0) != inputDataLeft.dimension(1)) && (inputDataLeft.dimension(1) != 1) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimension of the right input data container and first dimension of the left data input container (number of integration points) must agree or first dimension of the left data input container must be 1!");
    			TEUCHOS_TEST_FOR_EXCEPTION( ( inputDataLeft.dimension(0) != outputData.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (ArrayTools::scalarMultiplyDataData): Zeroth dimensions of data output and left data input containers (number of integration domains) must agree!");
    			for (int i=0; i<inputDataRight.rank(); i++) {
      				std::string errmsg  = ">>> ERROR (ArrayTools::scalarMultiplyDataData): Dimensions ";
      				errmsg += (char)(48+i);
      				errmsg += " and ";
      				errmsg += (char)(48+i+1);
      				errmsg += " of the right input and output data containers must agree!";
      				TEUCHOS_TEST_FOR_EXCEPTION( (inputDataRight.dimension(i) != outputData.dimension(i+1)), std::invalid_argument, errmsg );
    					}
#endif


  int numCells       = outputData.dimension(0);
  int numDataPoints  = inputDataLeft.dimension(1);


  
  if (numDataPoints != 1) {
         if (reciprocal) {
              Kokkos::parallel_for(numCells,PFor__Recip_Not_Const3_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip_Not_Const3_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight >(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
		else{
          if (reciprocal) {
         	  Kokkos::parallel_for(numCells,PFor__Recip__Const3_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
          else {
              Kokkos::parallel_for(numCells,PFor_Not_Recip__Const3_4<ArrayOutData,ArrayInDataLeft,ArrayInDataRight>(outputData, inputDataLeft, 
              inputDataRight));
          }
		}
      } 
  };  
  
	}//end namespace 
#endif	
#endif
