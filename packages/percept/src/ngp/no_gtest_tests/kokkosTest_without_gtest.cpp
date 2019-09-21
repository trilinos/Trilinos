#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>


#include <Kokkos_Core.hpp>

#include <iostream>

//#include "simple_mesh_adapt_unstructured_cuda.cpp"


#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda   ExecSpaceNoG ;
  typedef Kokkos::CudaSpace   MemSpaceNoG ;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpaceNoG ;
  typedef Kokkos::OpenMP     MemSpaceNoG ;
#else
  typedef Kokkos::Serial   ExecSpaceNoG ;
  typedef Kokkos::HostSpace   MemSpaceNoG ;
#endif


struct viewWrapper
{
	Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG> theView;
	
	viewWrapper(int _length, Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG>  passedView);
	
	void initializeView();

	int length;
};

viewWrapper::viewWrapper(int _length, Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG>  passedView)
{	
	length = _length;
	theView = passedView;
}

void viewWrapper::initializeView()
{
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceNoG>(0,length), KOKKOS_LAMBDA (int i){
		theView(i) = 5.7;
		printf("%g\n",theView(i));
	});


	/*
	double sum=0;
	Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceNoG>(0,length), KOKKOS_LAMBDA (int i,double & interimSum){
		interimSum+= theView(i); //complains that an illegal memory access has occured
	},sum);

	std::cout << sum << std::endl;
	*/
}


//Kokkos::View<double**[2], Kokkos::LayoutRight, MemSpace>::HostMirror host_inputmesh, host_outputmesh; 
//host_inputmesh  = Kokkos::create_mirror_view(inputmesh);  

void initializeViewFunctional(Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG> & passedView)
{
	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceNoG>(0,20), KOKKOS_LAMBDA (int i){
		passedView(i) = 5.7;
		printf("%g\n",passedView(i));
	});

	
}


int main(int argc, char **argv)
{	
	Kokkos::initialize();

	#ifdef KOKKOS_ENABLE_CUDA
		std::cout << "Running Test with Kokkos::Cuda execution space and Kokkos::CudaSpace memory space" <<std::endl;
	#elif defined(KOKKOS_ENABLE_OPENMP)
		std::cout << "Running Test with Kokkos::OpenMP execution space and Kokkos::OpenMP memory space" <<std::endl;
	#else
		std::cout << "Running Test with Kokkos::Serial execution space and Kokkos::HostSpace memory space" <<std::endl;
	#endif

	Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG> passedView("passedView", 20);

	viewWrapper wrapper(20, passedView);

	wrapper.initializeView(); //if I comment out this function use and just initialize the view from within this main function it works.


	//Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceNoG>(0,20), KOKKOS_LAMBDA (int i){
	//	passedView(i) = 5.7;
	//	printf("%g\n",passedView(i));
	//});

	//initializeViewFunctional(passedView); //works this way, too.
	//initializeViewFunctional(viewWrapper.theView); //doesn't like this


	/*
	Kokkos::View<double*, Kokkos::LayoutRight, MemSpaceNoG>::HostMirror hostView;

	hostView = Kokkos::create_mirror_view(passedView);

	Kokkos::deep_copy(hostView, passedView);

	for(int i=0;i<20;i++)
		std::cout<< hostView(i) << std::endl;

	
	*/
	Kokkos::finalize();
	return 0;
	
}


//ISSUES:
/*
	WITHIN SIERRA:

	I can build and run a simple unstrucutured mesh example that will run WITH LAMBDAS. The catch is that I have to use functions outside of the class to run the PF's otherwise it will complain about illegal memory access at run time

	I can build and run a simple class or struct that essentially is a view with an member initializer function that uses a lambda. However, if I try and run the initializer function it complains about illegal memory access.
		Again, if I approach this purely functionally, I can initialize the view.



	OUTSIDE SIERRA:

	I can get a simple case in which I declare a view to build and run. However, as soon as I try and add a PF with a lambda, it errors during the build
		We altered some build options, now this is working
*/



