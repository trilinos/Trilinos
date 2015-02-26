

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <Kokkos_Random.hpp>

template<class ViewType>
struct ComputeFunctor{
  ViewType outputData;
  const  ViewType inputData1;
  const ViewType inputData2;
  ComputeFunctor (ViewType outputData_, const ViewType inputData1_, 
               const ViewType inputData2_):outputData(outputData_),inputData1(inputData1_),inputData2(inputData2_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < outputData.dimension_1(); j++){
          for(int k = 0; k < outputData.dimension_2(); k++){
                      outputData(i,j,k) = inputData1(i,j,k)*inputData2(i,j,k)+inputData2(i,j,k);
          }
     }
  }
};


//typedef Kokkos::TeamVectorPolicy<32> team_policy;
//typedef team_policy::member_type team_member ;
#ifdef KOKKOS_HAVE_CXX11
/*
template<class ViewType>
struct NestedComputeFunctor {
  ViewType outputData;
  const ViewType inputData1;
  const ViewType inputData2;
  NestedComputeFunctor(ViewType outputData_, const ViewType inputData1_, 
              const ViewType inputData2_):outputData(outputData_),inputData1(inputData1_),inputData2(inputData2_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member & thread) const {
	 int i = thread.league_rank(); 
	 thread.team_par_for(outputData.dimension_1(),[&] (const int& j){
		 thread.vector_par_for(outputData.dimension_2(),[&] (const int& k){
			 outputData(i,j,k) = inputData1(i,j,k)*inputData2(i,j,k)+inputData2(i,j,k);
			 });		  	
		 });

  }
};
*/
#endif

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
int main(){
	const int loop2=8,loop3=3;
	int myints[] ={10,25,50,100,250,500,1000,2500,5000};
        std::vector<int>loop1(myints, myints + sizeof(myints) / sizeof(int) );
        //std::vector<int>loop1={100, 250, 500, 1000, 2500, 5000};
	int hyperthreads=0;
	std::vector<double>regularfortimevector;
	std::vector<double>parallelfortimevector;
	std::vector<double>nestedparallelfortimevector;
#if defined( KOKKOS_HAVE_CUDA )
  //Kokkos::Cuda::host_mirror_device_type::initialize();
  Kokkos::initialize();
  hyperthreads =256;
  std::cout << "CUDA device has been initialized" <<std::endl;
#else

 #if defined( KOKKOS_HAVE_OPENMP )
 int num_threads = 4;

 if (Kokkos::hwloc::available()) {
      std::cout <<"hwloc"<<std::endl;
      num_threads = Kokkos::hwloc::get_available_numa_count()
                    * Kokkos::hwloc::get_available_cores_per_numa()
                    * Kokkos::hwloc::get_available_threads_per_core()
                    ;
   }
  Kokkos::OpenMP::initialize( num_threads);
  std::cout << "OpenMP device has been initialized" <<std::endl;
  std::cout <<"   number of threads= " << num_threads<<std::endl;
  std::cout << "available threads: " << omp_get_max_threads() << std::endl;
 #endif
	 if ( Kokkos::hwloc::available() ) {
    std::cout << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
        << "] x CORE["    << Kokkos::hwloc::get_available_cores_per_numa()
        << "] x HT["      << Kokkos::hwloc::get_available_threads_per_core()
        << "] )"
        << std::endl ;
       hyperthreads = Kokkos::hwloc::get_available_threads_per_core();
  }
#endif
std::cout << "hyperthreads= " <<hyperthreads<< std::endl;


for(int itt=0;itt<loop1.size();itt++){
	//create kokkos views    
	Kokkos::View<double***> inputview1("X",loop1[itt],loop2,loop3);
	Kokkos::View<double***> inputview2("Y",loop1[itt],loop2,loop3);
	Kokkos::View<double***> outputview2("Z",loop1[itt],loop2,loop3);
        Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
        Kokkos::fill_random(inputview1,rand_pool64,100);	
        Kokkos::fill_random(inputview2,rand_pool64,100);


	//generate random data
	/*
	for(int i=0;i<loop1[itt];i++){
		for(int j=0;j<loop2;j++){
			for(int k=0;k<loop3;k++){
				inputview1(i,j,k)=fRand(0.0,1.0);
				inputview2(i,j,k)=fRand(0.0,1.0);
				outputview2(i,j,k)=fRand(0.0,1.0);
			}
		}
	}*/	
std::cout <<"loopvalues test: "<<1000000/loop1[itt]<<"\n";
	//regular for loop time	
	Kokkos::Impl::Timer RegularForTimer;
/*	for(int w=0;w<1000000/loop1[itt];w++){
   	for(int i=0;i<loop1[itt];i++){
		for(int j=0;j<loop2;j++){
			for(int k=0;k<loop3;k++){
				outputview2(i,j,k)=inputview1(i,j,k)*inputview2(i,j,k)+inputview1(i,j,k);
			}
		}
	}	
}
*/
   	Kokkos::fence();  
  double regularfortime = RegularForTimer.seconds();
  regularfortimevector.push_back(regularfortime);
  std::cout <<"Regular For Time: "<<regularfortime<<"\n\n";
  
  //parallel for loop 
  Kokkos::Impl::Timer ParallelForTimer;
  for(int w=0;w<1000000/loop1[itt];w++){
     Kokkos::parallel_for(loop1[itt],ComputeFunctor<Kokkos::View<double***> >(outputview2, inputview1, inputview2));
  }
  Kokkos::fence();  
  double parallelfortime = ParallelForTimer.seconds();
  parallelfortimevector.push_back(parallelfortime);
   std::cout <<"Parallel For Time: "<<parallelfortime<<"\n\n";
   
#ifdef KOKKOS_HAVE_CXX11
  Kokkos::Impl::Timer NestedParallelForTimer;
  //nested parallel for loop
  //const team_policy policy( loop1[itt] , hyperthreads);
  for(int w=0;w<1000000/loop1[itt];w++){
		
    //Kokkos::parallel_for(policy,NestedComputeFunctor<Kokkos::View<double***> >(outputview2, inputview1,  inputview2));
 }
    double nestedparallelfortime = NestedParallelForTimer.seconds();
    nestedparallelfortimevector.push_back(nestedparallelfortime);
   std::cout <<"Nested Parallel For Time: "<<nestedparallelfortime<<"\n\n";
#endif
}

 std::cout << std::endl << "finalize" << std::endl;

#if defined( KOKKOS_HAVE_CUDA )
  Kokkos::finalize();
  //Kokkos::Cuda::host_mirror_device_type::finalize();
#else

#if defined( KOKKOS_HAVE_OPENMP )
  Kokkos::OpenMP::finalize();
#endif
#endif
std::ofstream ofs ("intelphiperformance.csv", std::ofstream::out);
for(int i=0;i<loop1.size();i++){
	ofs<<loop1[i]<<","<<regularfortimevector[i]<<","<<parallelfortimevector[i]<<","<<nestedparallelfortimevector[i]<<"\n";
}
ofs.close();
	return 0;
	}
