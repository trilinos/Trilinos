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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <Adelus.hpp>

int main(int argc, char *argv[])
{
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  int rank, size;

  int  myrows;
  int  mycols;
  int  myfirstrow;
  int  myfirstcol;
  int  myrhs;
  int  my_row;
  int  my_col;
  int  matrix_size;
  int  nprocs_per_row;
  int  nptile = 1; // number of processors per node
  int  numrhs = 1;

  double mflops;

  MPI_Comm rowcomm, colcomm;

  static int buf[4];

  int i, m, k;

  int mlen;   // Message length for input data

  unsigned int seed= 10;

  double secs;

  double eps;

  double othird;

  double four_thirds = 4./3.;

  double tempc;

  int result=1;

  // Enroll into MPI

  MPI_Init(&argc,&argv);                             /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);             /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);             /* get number of processes */
  MPI_Get_processor_name(processor_name, &name_len); /* get name of the processor */

  // Initialize Input buffer

  for(i=0;i<4;i++) buf[i]=-1;

  std::cout << "proc " << rank << " (" << processor_name << ") is alive of " << size << " Processors" << std::endl;

  if( rank == 0 ) {
    // Check for commandline input

    if (argc > 1) {
      // argv[1] should be size of matrix
      buf[0] = atoi(argv[1]);
      if (argc > 2) {
        // argv[2] should be #procs per row
        buf[1] = atoi(argv[2]);
        // argv[3] should be #procs per node
        buf[2] = atoi(argv[3]);
        // argv[4] should be #rhs
        buf[3] = atoi(argv[4]);
      }
      else {
        // default is 1, but sqrt(p) would be better
        buf[1] = 1; buf[2] = 1; buf[3] = 1;
      }
    }
    else {
      // Input Data about matrix and distribution

      if (buf[0] < 0) {
        std::cout << "Enter size of matrix " << std::endl;
        std::cin >> buf[0];
      }
      if (buf[1] < 0) {
        std::cout << "Enter number of processors to which each row is assigned "  << std::endl;
        std::cin >> buf[1];
      }
      if (buf[2] < 0) {
        std::cout << "Enter number of processors per node "  << std::endl;
        std::cin >> buf[2];
      }
      if (buf[3] < 0) {
        std::cout << "Enter number of rhs vectors "  << std::endl;
        std::cin >> buf[3];
      }
    }
  }

  // Send the initilization data to each processor
  mlen = 4*sizeof(int);

  MPI_Bcast(reinterpret_cast<char *>(buf), mlen, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Set the values where needed

  matrix_size = buf[0];

  nprocs_per_row = buf[1];

  nptile = buf[2];

  numrhs = buf[3];

  if( rank == 0 ) {
    std::cout << " Matrix Size " << matrix_size << std::endl;
    std::cout << " Processors in a row  "  << nprocs_per_row << std::endl;
    std::cout << " Processors in a node  " << nptile << std::endl;
    std::cout << " Number of RHS vectors " << numrhs << std::endl;
  }

  if( rank == 0) {
    std::cout << " ---- Building Adelus solver ----" << std::endl;
  }

  // Get Info to build the matrix on a processor

  Adelus::GetDistribution( MPI_COMM_WORLD,
                           nprocs_per_row, matrix_size, numrhs,
                           myrows, mycols, myfirstrow, myfirstcol,
                           myrhs, my_row, my_col );

  // Define new communicators: rowcomm and colcomm

  MPI_Comm_split(MPI_COMM_WORLD,my_row,my_col,&rowcomm);
  MPI_Comm_split(MPI_COMM_WORLD,my_col,my_row,&colcomm);

  std::cout << " ------ PARALLEL Distribution Info for : ---------" <<std::endl;

  std::cout << "   Processor  " << rank << std::endl
       << "    my rows  " << myrows << std::endl
       << "    my cols  " << mycols << std::endl
       << "    my rhs  " << myrhs << std::endl
       << "    my first col  " << myfirstcol  << std::endl
       << "    my first row  " << myfirstrow << std::endl
       << "    my_row  " << my_row << std::endl
       << "    num procs row   " << nprocs_per_row << std::endl
       << "    my_col  " << my_col << std::endl;

  // Adelus example using the Kokkos Views
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  int gpu_count;
#ifdef KOKKOS_ENABLE_CUDA
  cudaGetDeviceCount ( &gpu_count );
#else
  hipGetDeviceCount ( &gpu_count );
#endif
  if (nptile > gpu_count) {
    if( rank == 0 ) {
      std::cout << "Request more GPUs than the number of GPUs available "
                << "to MPI processes (requested: " << nptile
                << " vs. available: " << gpu_count
                << "). Exit without test." << std::endl;
    }
    MPI_Finalize() ;
    return 0;
  }

  Kokkos::InitializationSettings args;
  args.set_num_threads(1);
  std::cout << "   Processor  " << rank << " (" << processor_name << "), GPU: "
            << args.get_device_id() << "/" << gpu_count << std::endl;
  Kokkos::initialize( args );
#else
  Kokkos::initialize( argc, argv );
#endif
  {
  // Local size -- myrows  * (mycols + myrhs)

  using Layout = Kokkos::LayoutLeft;
#if defined(KOKKOS_ENABLE_CUDA)
  using TestSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
  using TestSpace = Kokkos::HIPSpace;
#else
  using TestSpace = Kokkos::HostSpace;
#endif
#ifdef DREAL
  using ViewMatrixType         = Kokkos::View<double**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<double*,  Layout, Kokkos::HostSpace>;
  using ViewMatrixType_Host    = Kokkos::View<double**, Layout, Kokkos::HostSpace>;
  using ViewNrmVectorType_Host = Kokkos::View<double*,  Layout, Kokkos::HostSpace>;
#elif defined(SREAL)
  using ViewMatrixType         = Kokkos::View<float**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<float*,  Layout, Kokkos::HostSpace>;
  using ViewMatrixType_Host    = Kokkos::View<float**, Layout, Kokkos::HostSpace>;
  using ViewNrmVectorType_Host = Kokkos::View<float*,  Layout, Kokkos::HostSpace>;
#elif defined(SCPLX)
  using ViewMatrixType         = Kokkos::View<Kokkos::complex<float>**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<Kokkos::complex<float>*,  Layout, Kokkos::HostSpace>;
  using ViewMatrixType_Host    = Kokkos::View<Kokkos::complex<float>**, Layout, Kokkos::HostSpace>;
  using ViewNrmVectorType_Host = Kokkos::View<float*,  Layout, Kokkos::HostSpace>;
#else
  using ViewMatrixType         = Kokkos::View<Kokkos::complex<double>**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<Kokkos::complex<double>*,  Layout, Kokkos::HostSpace>;
  using ViewMatrixType_Host    = Kokkos::View<Kokkos::complex<double>**, Layout, Kokkos::HostSpace>;
  using ViewNrmVectorType_Host = Kokkos::View<double*,  Layout, Kokkos::HostSpace>;
#endif

  using ViewIntType_Host= Kokkos::View<int*, Layout, Kokkos::HostSpace>;

  using execution_space = typename ViewMatrixType::device_type::execution_space;
  using memory_space    = typename ViewMatrixType::device_type::memory_space;
  using ScalarA         = typename ViewMatrixType::value_type;

  printf("Rank %d, ViewMatrixType execution_space %s, memory_space %s, value_type %s\n",rank, typeid(execution_space).name(), typeid(memory_space).name(), typeid(ScalarA).name());

  ViewMatrixType A( "A", myrows, mycols );
  ViewMatrixType B( "B", myrows, myrhs + 6 );

  ViewMatrixType::HostMirror h_A = Kokkos::create_mirror( A );
  ViewMatrixType::HostMirror h_B = Kokkos::create_mirror( B );

  // Some temp arrays

  ViewVectorType_Host temp  ( "temp", myrows );

  ViewVectorType_Host temp2 ( "temp2", myrows );

  ViewMatrixType_Host rhs   ( "rhs", matrix_size, numrhs );

  ViewMatrixType_Host temp3 ( "temp3", matrix_size, numrhs );

  ViewMatrixType_Host temp4 ( "temp4", matrix_size, numrhs );

  ViewMatrixType_Host tempp ( "tempp", matrix_size, numrhs );

  ViewMatrixType_Host temp22( "temp22", matrix_size, numrhs );

  ViewNrmVectorType_Host rhs_nrm( "rhs_nrm", numrhs );

  ViewNrmVectorType_Host m_nrm  ( "m_nrm", numrhs );

  ViewIntType_Host h_permute( "h_permute", matrix_size);// Permutation array for factor and solve done independently

  // Set Random values

  if( rank == 0 )
    std::cout << " ****   Setting Random Matrix    ****" << std::endl;

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed+rank);
  Kokkos::fill_random(A, rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());

  Kokkos::deep_copy( h_A, A );

  // Now Create the RHS

  if( rank == 0 )
    std::cout << " ****   Creating RHS   ****" << std::endl;

  // Sum the portion of the row that I have

  for (k= 0; k < myrows; k++) {
    temp(k) = 0;
    for (m=0; m < mycols; m++) {
     temp(k) = temp(k) + h_A(k,m);
    }
  }

  // Sum from all processes and distribute the result back to all processes in rowcomm

  MPI_Allreduce(temp.data(), temp2.data(), myrows, ADELUS_MPI_DATA_TYPE, MPI_SUM, rowcomm);

  // Find the location of my RHS in the global RHS

  int *nrhs_procs_rowcomm;
  int my_rhs_offset = 0;

  nrhs_procs_rowcomm  = (int*)malloc( nprocs_per_row * sizeof(int));
  MPI_Allgather(&myrhs, 1, MPI_INT, nrhs_procs_rowcomm, 1, MPI_INT, rowcomm);//gather numbers of rhs of other processes

  for (i=0; i<my_col; i++) {
    my_rhs_offset += nrhs_procs_rowcomm[i];
  }

  if( rank == 0 )
    std::cout << " ****   Backing up RHS   ****" << std::endl;

  // Now put the RHS in the appropriate position

  if( myrhs > 0 ) {
    for (k = 0; k < myrhs; k++) {
#if defined(DREAL) || defined(ZCPLX)
      ScalarA scal_factor = static_cast<double>(my_rhs_offset+k+1);
#else
      ScalarA scal_factor = static_cast<float>(my_rhs_offset+k+1);
#endif
      auto cur_rhs_vec_1d = subview(h_B,Kokkos::ALL(),k);
      Kokkos::deep_copy( cur_rhs_vec_1d, temp2 );
      KokkosBlas::scal(cur_rhs_vec_1d,scal_factor,cur_rhs_vec_1d);
    }
    for (k = 0; k < numrhs; k++) {
#if defined(DREAL) || defined(ZCPLX)
      ScalarA scal_factor = static_cast<double>(k+1);
#else
      ScalarA scal_factor = static_cast<float>(k+1);
#endif
      auto cur_rhs_vec_1d = subview(rhs,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows),k);
      Kokkos::deep_copy( cur_rhs_vec_1d, temp2 );
      KokkosBlas::scal(cur_rhs_vec_1d,scal_factor,cur_rhs_vec_1d);
    }
  }

  // Globally Sum the RHS needed for testing later

  MPI_Allreduce(rhs.data(), temp4.data(), matrix_size*numrhs, ADELUS_MPI_DATA_TYPE, MPI_SUM, colcomm);

  // Pack back into RHS

  Kokkos::deep_copy( rhs, temp4 );

  KokkosBlas::nrm2(rhs_nrm, rhs);

  Kokkos::deep_copy( B, h_B );


  // Create handle
  Adelus::AdelusHandle<typename ViewMatrixType::value_type, execution_space, memory_space>
    ahandle(0, MPI_COMM_WORLD, matrix_size, nprocs_per_row, numrhs );

  // Now Factor the matrix

  if( rank == 0 )
    std::cout << " ****   Beginning Matrix Factor   ****" << std::endl;

  Adelus::Factor (ahandle, A, h_permute, &secs);

  if( rank == 0) {
    std::cout << " ----  Factor time  ----   " << secs << "  in secs. " << std::endl;

    mflops = 2./3.*pow(matrix_size,3.)/secs/1000000.;

    std::cout << " *****   MFLOPS   *****  " << mflops << std::endl;
  }

  // Call Solve (1st time)

  if( rank == 0 )
    std::cout << " ****   Beginning Matrix Solve (1st)   ****" << std::endl;

  Adelus::Solve (ahandle, A, B, h_permute, &secs);

  if( rank == 0)
    std::cout << " ----  Solution time (1st)  ----   " << secs << "  in secs. " << std::endl;

  // Restore the orig. RHS for testing Adelus::Solve() on a pre-computed LU factorization
  Kokkos::deep_copy( B, h_B );

  // Call Solve (2nd time)
  if( rank == 0 )
    std::cout << " ****   Beginning Matrix Solve (2nd)   ****" << std::endl;

  Adelus::Solve (ahandle, A, B, h_permute, &secs);

  if( rank == 0)
    std::cout << " ----  Solution time (2nd)  ----   " << secs << "  in secs. " << std::endl;

  // Now Check the Solution

  Kokkos::deep_copy( h_B, B );


  // Pack the Answer into the apropriate position

  if ( myrhs > 0 ) {
    Kokkos::deep_copy( subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows),
                                     Kokkos::make_pair(my_rhs_offset, my_rhs_offset + myrhs)),
                       subview(h_B,Kokkos::ALL(),Kokkos::make_pair(0, myrhs)) );
  }

  // All processors get the answer

  MPI_Allreduce(tempp.data(), temp22.data(), matrix_size*numrhs, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);

  // Perform the Matrix vector product

  ScalarA alpha = 1.0;
  ScalarA beta  = 0.0;

  KokkosBlas::gemm("N", "N", alpha,
                   subview(h_A,Kokkos::ALL(),Kokkos::make_pair(0, mycols)),
                   subview(temp22,Kokkos::make_pair(myfirstcol - 1, myfirstcol - 1 + mycols),Kokkos::ALL()),
                   beta, subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows),Kokkos::ALL()));

  MPI_Allreduce(tempp.data(), temp3.data(), matrix_size*numrhs, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);

  if( rank == 0) {
    std::cout <<  "======================================" << std::endl;
    std::cout << " ---- Error Calculation ----" << std::endl;

    ScalarA alpha_ = -1.0;

    KokkosBlas::axpy(alpha_, rhs, temp3);//temp3=temp3-rhs

    KokkosBlas::nrm2(m_nrm, temp3);
  }

  // Machine epsilon Calculation

  othird = four_thirds - 1.;

  tempc = othird + othird + othird;

  eps = fabs(tempc-1.0);

  if ( rank == 0 ) {
	std::cout << "   Machine eps  " << eps  << std::endl;

    std::cout << "   Threshold = " << eps*1e4  << std::endl;

    for (k = 0; k < numrhs; k++) {
      std::cout << "   Solution " << k << ":   ||Ax - b||_2 = " << m_nrm(k) << std::endl;

      std::cout << "   Solution " << k << ":   ||b||_2 = " << rhs_nrm(k) << std::endl;

      std::cout << "   Solution " << k << ":   ||Ax - b||_2 / ||b||_2  = " << m_nrm(k)/rhs_nrm(k)  << std::endl;

      if ( m_nrm(k)/rhs_nrm(k)  > (eps*1e4)) {
        std::cout << " ****   Solution " << k << " Fails   ****" <<  std::endl;
        result = 1;
        break;
      }
      else {
        std::cout << " ****   Solution " << k << " Passes   ****" << std::endl;
        result = 0;
      }
    }
    std::cout <<  "======================================" << std::endl;
  }

  MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);

  free(nrhs_procs_rowcomm);

  }
  Kokkos::finalize();

  MPI_Finalize() ;

  return (result);
}
