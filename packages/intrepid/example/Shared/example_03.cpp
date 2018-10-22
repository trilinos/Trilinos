#include "Intrepid_ArrayTools.hpp"
#include <iostream>
#ifdef INTREPID_OLD_KOKKOS_CODE
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_ScalarTraits.hpp"




int main(){
#ifdef INTREPID_OLD_KOKKOS_CODE

   Kokkos::initialize();
  //initialize viewsto random values
	const int v=100,x=300,y=40,z=40;
	Kokkos::View<double****> inputview1("X",v,x,y,z);
	Kokkos::View<double****> inputview2("Y",v,x,y,z);
	Kokkos::View<double****> outputview2("Z",v,x,y,z);
	Intrepid::FieldContainer<double> inview2FieldContainer(v, x, y, z);
	Intrepid::FieldContainer<double> inview1FieldContainer(v, x, y, z);
    Intrepid::FieldContainer<double> outview2FieldContainer(v, x, y, z);
  //These are the wrapper structures that are used to avoid compiletime rank issues
    ArrayWrapper<double,Kokkos::View<double****>, Rank<Kokkos::View<double****> >::value,false>inputview1wrap(inputview1);
	ArrayWrapper<double,Kokkos::View<double****>, Rank<Kokkos::View<double****> >::value,false>outputview2wrap(outputview2);
	ArrayWrapper<double,Intrepid::FieldContainer<double>,Rank<Intrepid::FieldContainer<double> >::value,false>inputfieldcontainer2wrap(inview2FieldContainer);
	ArrayWrapper<double,Kokkos::View<double****>, Rank<Kokkos::View<double****> >::value,false>inputview2wrap(inputview2);
  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
        
  //test getrank partial template specialization for kokkos views and field containers
//	std::cout <<"Rankstuff: "<<getrank(inputview1)<<","<<getrank(inview2FieldContainer)<<std::endl;
	
//	int numDataLeftPts = inputview1.dimension(1);
    int numCells       = outputview2.dimension(0);
    int numPoints      = outputview2.dimension(1);
    int matDim         = outputview2.dimension(2);
//    std::cout <<numCells<<","<<numPoints<<","<<matDim<<std::endl; 
  //setup field container using values from kokkos view and set output to zeros
    for(int cell = 0; cell < numCells; cell++) {
        for(int point = 0; point < numPoints; point++) {
            for( int row = 0; row < matDim; row++) {
                for( int col = 0; col < matDim; col++) {
					inview2FieldContainer(cell, point, row, col)=inputview2(cell, point, row, col);
					inview1FieldContainer(cell, point, row, col)=inputview1(cell, point, row, col);					
					outputview2wrap(cell, point, row, col)=0.0;
					outview2FieldContainer(cell, point, row, col)=0.0;
                }// Col-loop
            } // Row-loop
        } // P-loop
    }// C-loop	*/
 	Kokkos::fence();
	Kokkos::Impl::Timer structviewstimer;
   //example with two kokkos views and structs	
    for(int cell = 0; cell < numCells; cell++) {
       for(int point = 0; point < numPoints; point++) {
          for(int row = 0; row < matDim; row++) {
             for(int col = 0; col < matDim; col++) {                   
                   outputview2wrap(cell,point,row,col)= inputview1wrap(cell,point,row,col)*inputview2wrap(cell, point, row, col);
                }// Col-loop
            } // Row-loop
        } // P-loop
    }// C-loop	
    Kokkos::fence();
    //double timestructviews = structviewstimer.seconds();
//    std::cout <<std::setprecision(9)<<"Time for structviews"<<timestructviews<<std::endl;
    

    Kokkos::fence();
    Kokkos::Impl::Timer rawviewstimer;
    
   //example with two kokkos views without structs
	for(int cell = 0; cell < numCells; cell++) {
       for(int point = 0; point < numPoints; point++) {
          for( int row = 0; row < matDim; row++) {
             for( int col = 0; col < matDim; col++) {
                    outputview2(cell, point, row, col) =
                    inputview1(cell,point, row, col)*inputview2(cell, point, row, col);
                }// Col-loop
            } // Row-loop
        } // P-loop
    }// C-loop
	 	
 	Kokkos::fence();	
  //  double timerawviews = rawviewstimer.seconds();
 //   std::cout <<"Time for rawviews"<<timerawviews<<std::endl;


 /*           
   //example with kokkos view and intrepid field container    
    for(int cell = 0; cell < numCells; cell++) {
       for(int point = 0; point < numPoints; point++) {
          for(int row = 0; row < matDim; row++) {
             for( int col = 0; col < matDim; col++) {
                    outputview2wrap(cell, point, row, col) =                    inputview1wrap(cell,point, row, col)*inputfieldcontainer2wrap(cell, point, row, col);
                }// Col-loop
            } // Row-loop
        } // P-loop
    }// C-loop	
  //Test passing first element of field container

*/

         Kokkos::fence();
 	Kokkos::Impl::Timer fieldcontainertimer;
    for(int cell = 0; cell < numCells; cell++) {
       for(int point = 0; point < numPoints; point++) {
          for(int row = 0; row < matDim; row++) {
             for( int col = 0; col < matDim; col++) {
                    outview2FieldContainer(cell, point, row, col) = 
                    inview1FieldContainer(cell,point, row, col)*inview2FieldContainer(cell, point, row, col);
                }// Col-loop
            } // Row-loop
        } // P-loop
    }// C-loop	
  

    Kokkos::fence();
//double timerfieldcontainermultiply=fieldcontainertimer.seconds();

//std::cout <<"FieldContainerTimer: "<<timerfieldcontainermultiply <<std::endl;

   Kokkos::finalize();
  
#endif 
  
	
	return 0;
}
