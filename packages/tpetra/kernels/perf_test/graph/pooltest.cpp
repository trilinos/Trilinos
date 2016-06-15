#include <iostream>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"





typedef long idx;

template <typename mpool_type, typename ExecutionSpace>
struct MemoryPoolTest{
  mpool_type my_memory_pool ;

  KOKKOS_INLINE_FUNCTION
  size_t get_thread_id(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type & teamMember) const{
#if defined( KOKKOS_HAVE_SERIAL )
  if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
    return 0;
  }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
  if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
    return Kokkos::Threads::hardware_thread_id();
  }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
  if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
    return Kokkos::OpenMP::hardware_thread_id();
  }
#endif

    return teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
  }

  MemoryPoolTest (mpool_type my_memory_pool_):
    my_memory_pool(my_memory_pool_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type & teamMember) const {
    volatile idx * myData = NULL;
    size_t tid = this->get_thread_id(teamMember);

    int trial = 0;
    while (myData == NULL){
      ++trial;
      Kokkos::single(Kokkos::PerThread(teamMember),[&] (volatile idx * &memptr) {
        memptr = (volatile idx * )this->my_memory_pool.allocate_chunk(tid);

      }, myData);

    }


    for (int i = 0; i < 1000; ++i);



    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      printf("me:%ld lr:%d ts:%d tr:%d, Had Memory location:%ld with chunk_index:%ld in this many tries:%d\n",
          tid, (int) teamMember.league_rank(), (int) teamMember.team_size(),   (int)  teamMember.team_rank(),
          myData,my_memory_pool.get_chunk_index((idx *) myData), trial );
      this->my_memory_pool.release_chunk((idx *) myData);
    });

    myData = NULL;
  }
};

template <typename MyExecSpace1>
void run_test(int TEAMSIZE, int VECTORSIZE){




  typedef Kokkos::TeamPolicy<MyExecSpace1> team_policy_t1 ;
  typedef typename KokkosKernels::Experimental::Util::UniformMemoryPool<MyExecSpace1, idx> simple_pool1;

  int chunk_size = 32;
  int space_1_concurrency  = MyExecSpace1::concurrency() / VECTORSIZE;

  simple_pool1 space1_one2one_mp (space_1_concurrency, chunk_size, -1, KokkosKernels::Experimental::Util::OneThread2OneChunk);
  simple_pool1 space1_many2one_mp (4, chunk_size, -1, KokkosKernels::Experimental::Util::ManyThread2OneChunk);

  MemoryPoolTest<simple_pool1, MyExecSpace1> t1 (space1_one2one_mp);
  MemoryPoolTest<simple_pool1, MyExecSpace1> t2 (space1_many2one_mp);


  printf("\n\nRunning One2One\n");
  printf("Memory Pool\n");
  space1_one2one_mp.print_memory_pool();


  Kokkos::parallel_for( team_policy_t1(16 , TEAMSIZE, VECTORSIZE), t1);
  MyExecSpace1::fence();

  printf("\n\nRunning Many2One\n");
  printf("Memory Pool\n");
  space1_many2one_mp.print_memory_pool();

  Kokkos::parallel_for( team_policy_t1(16, TEAMSIZE, VECTORSIZE), t2);
  MyExecSpace1::fence();


}

int  main (int  argc, char ** argv){
  //Kokkos::initialize(argc, argv);

  int threads = 16;
  if (argc > 1){
    threads = atoi(argv[1]);
  }

#if defined( KOKKOS_HAVE_SERIAL )
  printf("\n\n\nSERIAL RUN\n");
  Kokkos::Serial::initialize();
  Kokkos::Serial::print_configuration(std::cout);
  run_test<Kokkos::Serial>(1,1);
  Kokkos::Serial::finalize();
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
  printf("\n\n\nTHREADS RUN\n");
  Kokkos::Threads::initialize(threads);
  Kokkos::Threads::print_configuration(std::cout);
  run_test<Kokkos::Threads>(1,1);
  Kokkos::Threads::finalize();
#endif

#if defined( KOKKOS_HAVE_OPENMP )
  printf("\n\n\nOpenMP RUN\n");
  Kokkos::OpenMP::initialize(threads);
  Kokkos::OpenMP::print_configuration(std::cout);
  run_test<Kokkos::OpenMP>(1,1);
  Kokkos::OpenMP::finalize();
#endif

#if defined( KOKKOS_HAVE_CUDA )
  printf("\n\n\nCuda RUN\n");

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
  Kokkos::Cuda::print_configuration(std::cout);
  run_test<Kokkos::Cuda>(8,32);
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
#endif






  //openmp_one2one_mp.print_memory_pool();
  Kokkos::finalize();
  return 0;
}
