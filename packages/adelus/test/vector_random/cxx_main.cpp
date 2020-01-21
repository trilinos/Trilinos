/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <Pliris.hpp>
//#include "flops.h"

#ifdef DREAL
#define MPI_VALUE_TYPE MPI_DOUBLE 
#endif
#ifdef SREAL
#define MPI_VALUE_TYPE MPI_FLOAT
#endif
#ifdef ZCPLX
#define MPI_VALUE_TYPE MPI_DOUBLE_COMPLEX
#endif
#ifdef SCPLX
#define MPI_VALUE_TYPE MPI_COMPLEX 
#endif

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

  double mflops;

  MPI_Comm rowcomm;

  static int buf[4];
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
      }
      else
        // default is 1, but sqrt(p) would be better
        buf[1] = 1;
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
    std::cout << " Matrix Size "<< buf[0] <<"\n";
    std::cout << " Processors in a row  " <<  buf[1] << "\n";
  }

  /* Send the initilization data to each processor    */
  mlen = 4*sizeof(int);

  MPI_Bcast(reinterpret_cast<char *>(buf), mlen, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Set the values where needed

  matrix_size = buf[0];

  nprocs_per_row  = buf[1];

  // Example for 1 RHS

  numrhs = 1;

  if( rank == 0) {
    std::cout << " ---- Building Pliris solver ----" << std::endl;
  }

  // Get Info to build the matrix on a processor

  Pliris::GetDistribution( &nprocs_per_row,
                           &matrix_size,
                           &numrhs,
                           &myrows,
                           &mycols,
                           &myfirstrow,
                           &myfirstcol,
                           &myrhs,
                           &my_row,
                           &my_col );

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

  // Pliris example using the Kokkos Views
  Kokkos::initialize( argc, argv );
  {
  //  Local size -- myrows  * (mycols + myrhs)
  
  typedef Kokkos::LayoutLeft Layout;
  typedef Kokkos::View<Kokkos::complex<double>**, Layout>  ViewMatrixType;
  typedef Kokkos::View<Kokkos::complex<double>*,  Layout, Kokkos::HostSpace>  ViewVectorType_Host;
  typedef typename ViewMatrixType::device_type::execution_space execution_space;
  typedef typename ViewMatrixType::device_type::memory_space memory_space;
  printf("Rank %d, ViewMatrixType execution_space %s, memory_space %s\n",rank,typeid(execution_space).name(),typeid(memory_space).name());

  typedef typename ViewMatrixType::value_type ScalarA;

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
  //Kokkos::fill_random(A, rand_pool,ScalarA(10));
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

  MPI_Allreduce(temp.data(), temp2.data(), myrows, MPI_VALUE_TYPE, MPI_SUM, rowcomm);

  if( rank == 0 )
    std::cout << " ****   Packing RHS in Matrix   ****" << std::endl;

  // Now put the RHS in the appropriate position

  if( myrhs > 0 ) {
    Kokkos::deep_copy( subview(h_A,Kokkos::ALL(),mycols), temp2 );
    Kokkos::deep_copy( subview(rhs,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)), temp2 );
  }

  // Globally Sum the RHS needed for testing later

  MPI_Allreduce(rhs.data(), temp4.data(), matrix_size, MPI_VALUE_TYPE, MPI_SUM, MPI_COMM_WORLD);

  // Pack back into RHS

  Kokkos::deep_copy( rhs, temp4 );

  rhs_nrm = KokkosBlas::nrm2(rhs);

  Kokkos::deep_copy( subview(A,Kokkos::ALL(),mycols), subview(h_A,Kokkos::ALL(),mycols) );

  // Now Solve the Problem

  if( rank == 0 )
    std::cout << " ****   Beginning Matrix Solve   ****" << std::endl;

  Pliris::FactorSolve (A, myrows, mycols, &matrix_size, &nprocs_per_row, &numrhs, &secs);

  if( rank == 0) {
    std::cout << " ----  Solution time  ----   " << secs << "  in secs. " << std::endl;

    mflops = 2./3.*pow(matrix_size,3.)/secs/1000000.;

    std::cout << " *****   MFLOPS   *****  " << mflops << std::endl;
  }

  // Now Check the Solution

  Kokkos::deep_copy( subview(h_A,Kokkos::ALL(),mycols), subview(A,Kokkos::ALL(),mycols) );

  // Pack the Answer into the apropriate position

  if ( myrhs > 0) {
    Kokkos::deep_copy( subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)), subview(h_A,Kokkos::ALL(),mycols) );
  }

  // All processors get the answer

  MPI_Allreduce(tempp.data(), temp22.data(), matrix_size, MPI_VALUE_TYPE, MPI_SUM, MPI_COMM_WORLD);
  
  // perform the Matrix vector product

  ScalarA alpha = 1.0;
  ScalarA beta  = 0.0;

  KokkosBlas::gemv("N", alpha, subview(h_A,Kokkos::ALL(),Kokkos::make_pair(0, mycols)),
                               subview(temp22,Kokkos::make_pair(myfirstcol - 1, myfirstcol - 1 + mycols)),
                         beta, subview(tempp,Kokkos::make_pair(myfirstrow - 1, myfirstrow - 1 + myrows)));
	
  MPI_Allreduce(tempp.data(), temp3.data(), matrix_size, MPI_VALUE_TYPE, MPI_SUM, MPI_COMM_WORLD);

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

    if ( m_nrm/rhs_nrm  > (eps*1e4))
      std::cout << " ****    Solution Fails   ****" <<  std::endl;
    else
      std::cout << " ****   Solution Passes   ****" << std::endl;

    std::cout <<  "======================================" << std::endl;
  }

  }
  Kokkos::finalize();

  MPI_Finalize() ;

  return (0);
}
