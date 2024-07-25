/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

// This is the example that drives Adelus to solve for a dense linear system
// with one RHS vector:
//1. Determine the matrix and RHS vector distribution
//2. Randomly generate the assigned matrix on each MPI process and the reference
//   solution vector
//3. Compute the assigned RHS vector on MPI processes
//4. Launch Adelus
//5. Gather results
//6. Compare the returned solution vector with the reference vector

#include <iostream>
#include <fstream>
#include <iterator>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_gemv.hpp>
#include "Adelus.hpp"

int main( int argc, char* argv[] )
{
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  int rank, size;

  unsigned int N       = 20;// number of rows (columns)
  unsigned int nrepeat = 5; // number of repeats of the test
  int nprocs_row       = 1; // number of processors per row
  int nrhs             = 1; // number of RHS
  int nptile           = 1; // number of processors per node
  int diagdom          = 0; // whether diagonal is dominant

  struct timeval begin, end;
  double secs;

  MPI_Init (&argc, &argv);                           // starts MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);             // get current process id
  MPI_Comm_size (MPI_COMM_WORLD, &size);             // get number of processes
  MPI_Get_processor_name(processor_name, &name_len); // get name of the processor

  // Print off a hello world message
  printf("Hello world from node %s, rank %d out of %d processors\n", processor_name, rank, size);

  MPI_Barrier (MPI_COMM_WORLD);

  // Read command line arguments
  for ( int i = 0; i < argc; i++ ) {
    if ( strcmp( argv[ i ], "-N" ) == 0 ) {
      N = atol( argv[ ++i ] );
      printf( "  User N is %ld\n", N );
    }
    else if ( strcmp( argv[ i ], "-npr" ) == 0 ) {
      nprocs_row = atoi( argv[ ++i ] );
      if(rank==0) printf( "  User number of processors per row is %d\n", nprocs_row );
    }
    else if ( strcmp( argv[ i ], "-nrhs" ) == 0 ) {
      nrhs = atoi( argv[ ++i ] );
      if(rank==0) printf( "  User number of RHS vectors is %d\n", nrhs );
    }
    else if ( strcmp( argv[ i ], "-nrepeat" ) == 0 ) {
      nrepeat = atoi( argv[ ++i ] );
      if(rank==0) printf( "  User number of repetitions is %d\n", nrepeat );
    }
    else if ( strcmp( argv[ i ], "-nptile" ) == 0 ) {
      nptile = atoi( argv[ ++i ] );
      if(rank==0) printf( "  User number of processors per node is %d\n", nptile );
    }
    else if ( strcmp( argv[ i ], "-diag" ) == 0 ) {
      diagdom = atoi( argv[ ++i ] );
      if(rank==0) printf( "  User whether diagonal is dominant %d\n", diagdom );
    }
    else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "-help" ) == 0 ) ) {
      if(rank==0) printf( "  Options:\n" );
      if(rank==0) printf( "  -N <int>      :   determines number of rows (and columns) (default: 20)\n" );
      if(rank==0) printf( "  -npr <int>    :   number of processors per row (default: 1)\n" );
      if(rank==0) printf( "  -nrhs <int>   :   number of RHS vectors (default: 1)\n" );
      if(rank==0) printf( "  -nrepeat <int>:   number of repetitions (default: 5)\n" );
      if(rank==0) printf( "  -nptile <int> :   number of processors per node (default: 1)\n" );
      if(rank==0) printf( "  -diag <int>   :   whether matrix diagonal is dominant (0: rand, 1: dominant, default: 0)\n" );
      if(rank==0) printf( "  -help (-h)    :   print this message\n\n" );
      exit( 1 );
    }
  }

  printf("  Rank %d, N = %d, number of processors per row = %d, number of RHS = %d\n", rank, N, nprocs_row, nrhs);

  MPI_Barrier (MPI_COMM_WORLD);

  // 1. Determine distribution
  int matrix_size = N;
  int my_rows;
  int my_cols;
  int my_first_row;
  int my_first_col;
  int my_rhs;
  int my_row;
  int my_col;

  int my_rows_max;
  int my_cols_max;

  Adelus::GetDistribution( MPI_COMM_WORLD,
                           nprocs_row, matrix_size, nrhs,
                           my_rows, my_cols, my_first_row, my_first_col,
                           my_rhs, my_row, my_col );

  MPI_Allreduce( &my_rows, &my_rows_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce( &my_cols, &my_cols_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  fprintf(stderr,"GetDistribution -- rank %d, nprocs %d, my_rows %d, my_cols %d, my_first_row %d, my_first_col %d, my_rhs %d, my_row %d, my_col %d, my_rows_max %d, my_cols_max %d\n", rank,size,my_rows,my_cols,my_first_row,my_first_col,my_rhs,my_row,my_col,my_rows_max,my_cols_max);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  int gpu_count;
#ifdef KOKKOS_ENABLE_CUDA
  cudaGetDeviceCount ( &gpu_count );
#else
  hipGetDeviceCount ( &gpu_count );
#endif
  Kokkos::InitializationSettings args;
  args.set_num_threads(1);
  printf("Rank %d, Before Kokkos initialization, GPU %d/%d\n", rank, args.get_device_id(), gpu_count);

  MPI_Barrier (MPI_COMM_WORLD);

  Kokkos::initialize( args );
#else
  MPI_Barrier (MPI_COMM_WORLD);

  Kokkos::initialize( argc, argv );
#endif
  {
    char* env_var;
    env_var = getenv ("OMP_NUM_THREADS");
    if (env_var!=NULL) printf("Rank %d, the current OMP_NUM_THREADS is: %s\n",rank, env_var);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    int device_id;
#ifdef KOKKOS_ENABLE_CUDA
    cudaGetDevice ( &device_id );
#else
    hipGetDevice ( &device_id );
#endif
  printf("Rank %d, After Kokkos initialization, GPU %d/%d\n", rank, device_id, gpu_count);
#endif

    typedef Kokkos::LayoutLeft Layout; //Note: matrix is stored column-wise (column-major, LayoutLeft)

    // Allocate matrices and vectors on device
#ifdef DREAL
    typedef Kokkos::View<double**, Layout>  ViewMatrixType;
    typedef Kokkos::View<double*,  Layout>  ViewVectorType;
    typedef Kokkos::View<double**, Layout, Kokkos::HostSpace>  ViewMatrixType_Host;
    typedef Kokkos::View<double*,  Layout, Kokkos::HostSpace>  ViewVectorType_Host;
#elif defined(SREAL)
    typedef Kokkos::View<float**, Layout>  ViewMatrixType;
    typedef Kokkos::View<float*,  Layout>  ViewVectorType;
    typedef Kokkos::View<float**, Layout, Kokkos::HostSpace>  ViewMatrixType_Host;
    typedef Kokkos::View<float*,  Layout, Kokkos::HostSpace>  ViewVectorType_Host;
#elif defined(SCPLX)
    typedef Kokkos::View<Kokkos::complex<float>**, Layout>  ViewMatrixType;
    typedef Kokkos::View<Kokkos::complex<float>*,  Layout>  ViewVectorType;
    typedef Kokkos::View<Kokkos::complex<float>**, Layout, Kokkos::HostSpace>  ViewMatrixType_Host;
    typedef Kokkos::View<Kokkos::complex<float>*,  Layout, Kokkos::HostSpace>  ViewVectorType_Host;
#else
    typedef Kokkos::View<Kokkos::complex<double>**, Layout>  ViewMatrixType;
    typedef Kokkos::View<Kokkos::complex<double>*,  Layout>  ViewVectorType;
    typedef Kokkos::View<Kokkos::complex<double>**, Layout, Kokkos::HostSpace>  ViewMatrixType_Host;
    typedef Kokkos::View<Kokkos::complex<double>*,  Layout, Kokkos::HostSpace>  ViewVectorType_Host;
#endif
    typedef typename ViewMatrixType::device_type::execution_space execution_space;
    typedef typename ViewMatrixType::device_type::memory_space memory_space;
    typedef typename ViewMatrixType::value_type ScalarA;

    printf("Rank %d, ViewMatrixType execution_space %s, memory_space %s, value_type %s\n", rank, typeid(execution_space).name(), typeid(memory_space).name(), typeid(ScalarA).name());

    ViewMatrixType my_A( "my_A", my_rows, my_cols + my_rhs + 6 );//store the local matrix and RHS, padding 6 to prevent possible seg fault during factor and backsolve
    ViewVectorType my_B( "my_B", my_rows );                      //store the partial RHS vector on device
    ViewVectorType X0  ( "X0", N );                              //store the reference solution vector on device

    // Create host views
    ViewMatrixType::HostMirror h_my_A = Kokkos::create_mirror( my_A );     //backup data for multiple runs
#if defined(HOSTPTR_API) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    ViewMatrixType::HostMirror h_my_A_hptr = Kokkos::create_mirror( my_A );//in place of my_A, for testing host pointer API with CUDA or HIP enabled
#endif
    ViewVectorType_Host h_X ( "h_X",  N );                                 //store the final solution vector on host
    ViewVectorType_Host h_X0( "h_X0", N );                                 //store the referencen solution vector on host for error checking

    //2. Randomly generate the assigned matrix on each MPI process and the reference solution vector
    uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + rank;
    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

    // Reference solution (assumming one solution vector)
    gettimeofday( &begin, NULL );
    if(rank == 0) {//generate on rank 0
      Kokkos::fill_random(X0,rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());
      Kokkos::fence();//Note: need a fence here to make sure "fill_random" completes before broadcast
    }
    MPI_Bcast(reinterpret_cast<char *>(X0.data()), N*sizeof(ScalarA), MPI_CHAR, 0, MPI_COMM_WORLD);//broadcast X0 to other nodes
    Kokkos::deep_copy( h_X0, X0 );//copy from device to host
    gettimeofday( &end, NULL );
    printf( "  Rank %d generates reference solution vector: time( %lf s )\n", rank, 1.0 *    ( end.tv_sec  - begin.tv_sec  ) +
                                                                                    1.0e-6 * ( end.tv_usec - begin.tv_usec ));

    // Impedance matrix A
    gettimeofday( &begin, NULL );
    Kokkos::fill_random(my_A, rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());

    if (diagdom) {
      //Adjust diagonal
      Kokkos::parallel_for( Kokkos::RangePolicy<execution_space>(0,N), KOKKOS_LAMBDA ( const int& i ) {
#ifdef DREAL
        ScalarA tmp = 310.0 + (double)(rank*10);
#elif defined(SREAL)
        ScalarA tmp = 310.0 + (float)(rank*10);
#elif defined(SCPLX)
        ScalarA tmp(310.0 + (float)(rank*10),250.0 + (float)(rank*10));
#else
        ScalarA tmp(310.0 + (double)(rank*10),250.0 + (double)(rank*10));
#endif
        if ( ( (i-(my_first_row-1)) >= 0 ) && ( (i-(my_first_col-1)) >= 0 ) &&
             ( (i-(my_first_row-1)) < my_rows ) && ( (i-(my_first_col-1)) < my_cols ) )
          my_A(i-(my_first_row-1),i-(my_first_col-1)) = my_A(i-(my_first_row-1),i-(my_first_col-1)) + tmp;
      });
    }

    Kokkos::deep_copy( h_my_A, my_A );//copy from device (local) to host (local)
    gettimeofday( &end, NULL );
    printf( "  Rank %d generates impedance matrix: time( %lf s )\n", rank, 1.0 *    ( end.tv_sec  - begin.tv_sec  ) +
                                                                           1.0e-6 * ( end.tv_usec - begin.tv_usec ));

    //3. Compute the assigned RHS vector on MPI processes
    // RHS B = A*X0 (assumming 1 RHS vector)
    gettimeofday( &begin, NULL );
    ScalarA alpha = 1.0;
    ScalarA beta  = 0.0;

    auto sub_my_A = subview(my_A, Kokkos::ALL(), Kokkos::make_pair(0, my_cols));
    auto sub_X0   = subview(X0,   Kokkos::make_pair(my_first_col-1, my_first_col-1 + my_cols));

    KokkosBlas::gemv("N",alpha,sub_my_A,sub_X0,beta,my_B); Kokkos::fence();

    if (nprocs_row==1)
      Kokkos::deep_copy( subview(my_A,Kokkos::ALL(),my_cols), my_B );
    else {
      MPI_Comm rowcomm;//define a row based communicator
      MPI_Comm_split(MPI_COMM_WORLD, my_row, my_col, &rowcomm);
      MPI_Reduce(my_B.data(), my_A.data()+(unsigned long long int)my_rows*(unsigned long long int)my_cols, my_rows, ADELUS_MPI_DATA_TYPE, MPI_SUM, 0, rowcomm);
    }
    Kokkos::deep_copy( h_my_A, my_A );//copy from device (local) to host (local)

    gettimeofday( &end, NULL );
    printf( "  Rank %d generates RHS vector: time( %lf s )\n", rank, 1.0 *    ( end.tv_sec  - begin.tv_sec  ) +
                                                                     1.0e-6 * ( end.tv_usec - begin.tv_usec ));

    MPI_Barrier (MPI_COMM_WORLD);

    //4. Launch Adelus
    if(rank == 0) {
      printf("Using ADELUS\n");
#ifdef KOKKOS_ENABLE_CUDA
      printf("CUDA is enabled\n");
#endif
#ifdef KOKKOS_ENABLE_HIP
      printf("HIP is enabled\n");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
      printf("OPENMP is enabled\n");
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
      printf("Found MKL TPL\n");
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
      printf("Found BLAS TPL\n");
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
      printf("Found cuBLAS TPL\n");
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
      printf("Found rocBLAS TPL\n");
#endif
    }

    // Create handle
    Adelus::AdelusHandle<typename ViewMatrixType::value_type, execution_space, memory_space>
      ahandle(0, MPI_COMM_WORLD, matrix_size, nprocs_row, nrhs );

    double time = 0.0;

    MPI_Barrier (MPI_COMM_WORLD);

    for ( int repeat = 0; repeat < nrepeat; repeat++ ) {
#if defined(HOSTPTR_API) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
      Kokkos::deep_copy( h_my_A_hptr, h_my_A ); //restore the orig. matrix and RHS
#else
      Kokkos::deep_copy( my_A, h_my_A ); //restore the orig. matrix and RHS
#endif

      Kokkos::fence();

      MPI_Barrier (MPI_COMM_WORLD);

      gettimeofday( &begin, NULL );

#ifdef KKVIEW_API
      Adelus::FactorSolve (ahandle, my_A, &secs);
#endif
#if defined(DEVPTR_API) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
      Adelus::FactorSolve_devPtr (ahandle, reinterpret_cast<ADELUS_DATA_TYPE *>(my_A.data()), my_rows, my_cols, my_rhs, &matrix_size, &nprocs_row, &nrhs, &secs);
#endif
#if defined(HOSTPTR_API) && !(defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))//KOKKOS_ENABLE_OPENMP
      Adelus::FactorSolve_hostPtr (ahandle, reinterpret_cast<ADELUS_DATA_TYPE *>(my_A.data()), my_rows, my_cols, my_rhs, &matrix_size, &nprocs_row, &nrhs, &secs);
#endif
#if defined(HOSTPTR_API) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
      Adelus::FactorSolve_hostPtr (ahandle, reinterpret_cast<ADELUS_DATA_TYPE *>(h_my_A_hptr.data()), my_rows, my_cols, my_rhs, &matrix_size, &nprocs_row, &nrhs, &secs);
#endif

      Kokkos::fence();

      gettimeofday( &end, NULL );

      MPI_Barrier (MPI_COMM_WORLD);

      // Calculate time
      double time_iter = 1.0 *    ( end.tv_sec  - begin.tv_sec  ) +
                         1.0e-6 * ( end.tv_usec - begin.tv_usec );

      printf( "  Real runs on rank %d: repeat ( %d ) time( %lf s )\n", rank, repeat, time_iter);

      time += time_iter;
    }

#if defined(HOSTPTR_API) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    Kokkos::deep_copy( my_A, h_my_A_hptr );
#endif

    //5. Gather results (solution vector) (assuming nrhs=1)
    if (size==1) {
      Kokkos::deep_copy( h_X, subview(my_A,Kokkos::ALL(),my_cols) );//copy directly from device (local) to host (global)
    }
    else {// size>1
      if (nprocs_row==size) {
        if (rank == 0) {
          Kokkos::deep_copy( h_X, subview(my_A,Kokkos::ALL(),my_cols) );
        }
      }
      else {//nprocs_row<size
        memset( h_X.data(), 0, N*sizeof(ScalarA) );
        if (my_rhs > 0) {
            auto sub_my_A = subview(my_A, Kokkos::ALL(), my_cols);
            auto sub_h_X  = subview(h_X,  Kokkos::make_pair(my_first_row-1, my_first_row-1+my_rows));
            Kokkos::deep_copy( sub_h_X, sub_my_A );
        }
        MPI_Allreduce( MPI_IN_PLACE, h_X.data(), N, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);
      }
    }

    //6. Compare the returned solution vector with the reference vector on host
    if(rank == 0) {
      for (unsigned int i=0; i<N; i++) {
        if ( abs(h_X(i) - h_X0(i)) > 1e-8 ) {
#if defined(SCPLX) || defined(ZCPLX)
          printf( "    Error: result( %.12lf + j%.12lf ) != solution( %.12lf + j%.12lf ) at (%ld)\n", h_X(i).real(), h_X(i).imag(), h_X0(i).real(), h_X0(i).imag(), i );
#else
          printf( "    Error: result( %.12lf ) != solution( %.12lf ) at (%ld)\n", h_X(i), h_X0(i), i );
#endif
        }
      }
    }

    double Gbytes    = 1.0e-9 * double( sizeof(ScalarA) * ( N*N + N + N ) );
    double my_Gbytes = 1.0e-9 * double( sizeof(ScalarA) * ( my_rows*(my_cols + my_rhs + 6) + my_rows + N ) );

    // Print results (problem size, time).
    printf( "  Rank %d: N ( %d ), diag dominant ( %d ), nrepeat ( %d ), total problem( %g MB ), my problem( %g MB ), totaltime( %g s ), avg_time( %g s (totaltime/nrepeat))\n", rank, N, diagdom, nrepeat, Gbytes * 1000, my_Gbytes * 1000, time, time/nrepeat);
  }
  Kokkos::finalize();

  MPI_Finalize();

  return 0;
}
