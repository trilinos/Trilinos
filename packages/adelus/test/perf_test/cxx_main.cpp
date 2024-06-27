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

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Adelus.hpp"

int main(int argc, char *argv[])
{
  using Teuchos::MpiComm;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::StackedTimer;

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

  double mflops;

  MPI_Comm rowcomm;

  static int buf[3];
  int numrhs;

  int i, m, k;

  int mlen;   // Message length for input data

  unsigned int seed= 10;

  double secs;

  double eps;

  double othird;

  double four_thirds = 4./3.;

  double tempc;

  double rhs_nrm, m_nrm;

  int result=1;

  // Enroll into MPI

  MPI_Init(&argc,&argv);                             /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);             /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);             /* get number of processes */
  MPI_Get_processor_name(processor_name, &name_len); /* get name of the processor */

  // Make a Comm.  This one happens to wrap MPI_COMM_WORLD.
  RCP<const Comm<int> > comm = rcp (new MpiComm<int> (MPI_COMM_WORLD));

  // Initialize Input buffer

  for(i=0;i<3;i++) buf[i]=-1;

  std::cout << "proc " << rank << " (" << processor_name << ") is alive of " << size << " Processors" << std::endl;

  if( rank == 0 ) {
    // Check for command-line input
    if (argc > 1) {
      // argv[1] should be size of matrix
      buf[0] = atoi(argv[1]);
      if (argc > 2) {
        // argv[2] should be #procs per row
        buf[1] = atoi(argv[2]);
      }
      else {
        // default is 1, but sqrt(p) would be better
        buf[1] = 1;
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
    }
  }

  /* Send the initilization data to each processor    */
  mlen = 3*sizeof(int);

  MPI_Bcast(reinterpret_cast<char *>(buf), mlen, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Set the values where needed

  matrix_size = buf[0];

  nprocs_per_row = buf[1];

  if( rank == 0 ) {
    std::cout << " Matrix Size " << matrix_size << std::endl;
    std::cout << " Processors in a row  "  << nprocs_per_row << std::endl;
  }

  // Example for 1 RHS

  numrhs = 1;

  if( rank == 0) {
    std::cout << " ---- Building Adelus solver ----" << std::endl;
  }

  // Get Info to build the matrix on a processor

  Adelus::GetDistribution( MPI_COMM_WORLD,
                           nprocs_per_row, matrix_size, numrhs,
                           myrows, mycols, myfirstrow, myfirstcol,
                           myrhs, my_row, my_col );

  //   Define a new communicator

  MPI_Comm_split(MPI_COMM_WORLD,my_row,my_col,&rowcomm);

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
  Kokkos::InitializationSettings args;
  args.set_num_threads(1);
  std::cout << "   Processor  " << rank << " (" << processor_name << "), GPU: "
            << args.get_device_id() << "/" << gpu_count << std::endl;
  Kokkos::initialize( args );
#else
  Kokkos::initialize( argc, argv );
#endif
  {
  //  Local size -- myrows  * (mycols + myrhs)

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
#elif defined(SREAL)
  using ViewMatrixType         = Kokkos::View<float**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<float*,  Layout, Kokkos::HostSpace>;
#elif defined(SCPLX)
  using ViewMatrixType         = Kokkos::View<Kokkos::complex<float>**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<Kokkos::complex<float>*,  Layout, Kokkos::HostSpace>;
#else
  using ViewMatrixType         = Kokkos::View<Kokkos::complex<double>**, Layout, TestSpace>;
  using ViewVectorType_Host    = Kokkos::View<Kokkos::complex<double>*,  Layout, Kokkos::HostSpace>;
#endif
  using execution_space = typename ViewMatrixType::device_type::execution_space;
  using memory_space    = typename ViewMatrixType::device_type::memory_space;
  using ScalarA         = typename ViewMatrixType::value_type;

  printf("Rank %d, ViewMatrixType execution_space %s, memory_space %s, value_type %s\n",rank, typeid(execution_space).name(), typeid(memory_space).name(), typeid(ScalarA).name());

  ViewMatrixType A( "A", myrows, mycols + myrhs + 6 );

  ViewMatrixType::HostMirror h_A = Kokkos::create_mirror( A );

  // Some temp arrays

  ViewVectorType_Host temp  ( "temp", myrows );

  ViewVectorType_Host temp2 ( "temp2", myrows );

  ViewVectorType_Host rhs   ( "rhs", matrix_size );

  ViewVectorType_Host temp3 ( "temp3", matrix_size );

  ViewVectorType_Host temp4 ( "temp4", matrix_size );

  ViewVectorType_Host tempp ( "tempp", matrix_size );

  ViewVectorType_Host temp22( "temp22", matrix_size );

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

  // Sum to Processor 0

  MPI_Allreduce(temp.data(), temp2.data(), myrows, ADELUS_MPI_DATA_TYPE, MPI_SUM, rowcomm);

  if( rank == 0 )
    std::cout << " ****   Packing RHS in Matrix   ****" << std::endl;

  // Now put the RHS in the appropriate position

  if( myrhs > 0 ) {
    Kokkos::deep_copy( subview(h_A,Kokkos::ALL(),mycols), temp2 );
    Kokkos::deep_copy( subview(rhs,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)), temp2 );
  }

  // Globally Sum the RHS needed for testing later

  MPI_Allreduce(rhs.data(), temp4.data(), matrix_size, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);

  // Pack back into RHS

  Kokkos::deep_copy( rhs, temp4 );

  rhs_nrm = KokkosBlas::nrm2(rhs);

  Kokkos::deep_copy( subview(A,Kokkos::ALL(),mycols), subview(h_A,Kokkos::ALL(),mycols) );

  // Create handle
  Adelus::AdelusHandle<typename ViewMatrixType::value_type, execution_space, memory_space>
    ahandle(0, MPI_COMM_WORLD, matrix_size, nprocs_per_row, numrhs );

  // Now Solve the Problem
  RCP<StackedTimer> timer = rcp(new StackedTimer("Adelus: total"));
  TimeMonitor::setStackedTimer(timer);

  if( rank == 0 )
    std::cout << " ****   Beginning Matrix Solve   ****" << std::endl;

  Adelus::FactorSolve (ahandle, A, &secs);

  if( rank == 0) {
    std::cout << " ----  Solution time  ----   " << secs << "  in secs. " << std::endl;

    mflops = 2./3.*pow(matrix_size,3.)/secs/1000000.;

    std::cout << " *****   MFLOPS   *****  " << mflops << std::endl;
  }

  timer->stopBaseTimer();
  StackedTimer::OutputOptions options;
  options.print_warnings = false;
  options.output_minmax  = true;
  timer->report(std::cout, comm, options);

  std::string testBaseName = "Adelus Factor and Solve ";
  auto xmlOut = timer->reportWatchrXML(testBaseName + std::to_string(size) + " ranks " + std::to_string(nprocs_per_row) + " procs_per_row", comm);
  if(rank == 0)
  {
    if(xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
  }

  // Now Check the Solution

  Kokkos::deep_copy( subview(h_A,Kokkos::ALL(),mycols), subview(A,Kokkos::ALL(),mycols) );

  // Pack the Answer into the apropriate position

  if ( myrhs > 0) {
    Kokkos::deep_copy( subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)), subview(h_A,Kokkos::ALL(),mycols) );
  }

  // All processors get the answer

  MPI_Allreduce(tempp.data(), temp22.data(), matrix_size, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);

  // perform the Matrix vector product

  ScalarA alpha = 1.0;
  ScalarA beta  = 0.0;

  KokkosBlas::gemv("N", alpha, subview(h_A,Kokkos::ALL(),Kokkos::make_pair(0, mycols)),
                               subview(temp22,Kokkos::make_pair(myfirstcol - 1, myfirstcol - 1 + mycols)),
                         beta, subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)));

  MPI_Allreduce(tempp.data(), temp3.data(), matrix_size, ADELUS_MPI_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD);

  if( rank == 0) {
    std::cout <<  "======================================" << std::endl;
    std::cout << " ---- Error Calculation ----" << std::endl;

    ScalarA alpha_ = -1.0;

    KokkosBlas::axpy(alpha_,rhs,temp3);//temp3=temp3-rhs

    m_nrm = KokkosBlas::nrm2(temp3);
  }

  // Machine epsilon Calculation

  othird = four_thirds - 1.;

  tempc = othird + othird + othird;

  eps = fabs(tempc-1.0);

  if ( rank == 0 ) {
	std::cout << "   Machine eps  " << eps  << std::endl;
  }

  if ( rank == 0 ) {

    std::cout << "   ||Ax - b||_2 = " << m_nrm << std::endl;

    std::cout << "   ||b||_2 = " << rhs_nrm << std::endl;

    std::cout << "   ||Ax - b||_2 / ||b||_2  = " << m_nrm/rhs_nrm  << std::endl;

    std::cout << "   Threshold = " << eps*1e4  << std::endl;

    if ( m_nrm/rhs_nrm  > (eps*1e4)) {
      std::cout << " ****    Solution Fails   ****" <<  std::endl;
      result = 1;
    }
    else {
      std::cout << " ****   Solution Passes   ****" << std::endl;
      result = 0;
    }
    std::cout <<  "======================================" << std::endl;
  }

  MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);

  }
  Kokkos::finalize();

  MPI_Finalize() ;

  return result;
}
