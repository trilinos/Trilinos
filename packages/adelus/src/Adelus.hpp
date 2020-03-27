/*
//@HEADER
// ************************************************************************
//
//                        Adelus v. 1.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// 3. Neither the name of NTESS nor the names of the contributors may be
// used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL NTESS OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Vinh Dang (vqdang@sandia.gov)
//                    Joseph Kotulski (jdkotul@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#pragma once

#include <Kokkos_View.hpp>
#include <Adelus_defines.h>
#include <Adelus_xlu_solve.hpp>
#include <Adelus_distribute.hpp>
#include <mpi.h>

// Adelus: provides the functionality to interface to a dense LU solver

namespace Adelus {

  /// Adelus GetDistirbution
  /// Gives the distribution information that is required by the dense solver
  
  /// \param nprocs_row_ (In)        - number of processors for a row
  /// \param number_of_unknowns (In) - order of the dense matrix
  /// \param nrhs_ (In)              - number of right hand sides
  /// \param my_rows_ (Out)          - number of rows of the matrix on this processor
  /// \param my_cols_ (Out)          - number of columns of the matrix on this processor
  /// \param my_first_row_ (Out)     - first (global) row number on this processor (array starts at index 1)
  /// \param my_first_col_ (Out)     - first (global) column number on this processor (array starts at index 1)
  /// \param my_rhs_ (Out)           - number of right hand sides on this processor
  /// \param my_row (Out)            - row number in processor mesh, 0 to the  number of processors for a column -1
  /// \param my_col (Out)            - column number in processor mesh, 0 to the  number of processors for a row -1
    
  inline
  int GetDistribution( int* nprocs_row_,
                       int* number_of_unknowns,
                       int* nrhs_,
                       int* my_rows_,
                       int* my_cols_,
                       int* my_first_row_,
                       int* my_first_col_,
                       int* my_rhs_,
                       int* my_row,
                       int* my_col ) {
    // This function echoes the multiprocessor distribution of the matrix

    distmat_(nprocs_row_,
             number_of_unknowns,
             nrhs_,
             my_rows_,
             my_cols_,
             my_first_row_,
             my_first_col_,
             my_rhs_,
             my_row,
             my_col);

    return(0);

  }

  /// Adelus FactorSolve
  /// Factors and solves the dense matrix

  /// \param AA (InOut)       -- Kokkos View that has the matrix and rhs packed (Note: matrix and rhs are overwritten)
  /// \param my_rows_ (In)    -- number of rows of the matrix on this processor
  /// \param my_cols_ (In)    -- number of columns of the matrix on this processor
  /// \param matrix_size (In) -- order of the dense matrix
  /// \param num_procsr (In)  -- number of processors for a row
  /// \param num_rhs (In)     -- number of right hand sides
  /// \param secs (Out)       -- factor and solve time in seconds
    
  template<class ZDView>
  inline
  void FactorSolve( ZDView AA,
                    int my_rows_,
                    int my_cols_,
                    int* matrix_size,
                    int* num_procsr,
                    int* num_rhs,
                    double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
#ifdef PRINT_STATUS
    printf("FactorSolve (Kokkos View interface) in rank %d -- my_rows %u , my_cols %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

  }

#ifdef ZCPLX
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  inline
  void FactorSolve_devPtr( ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<double> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  inline
  void FactorSolve_hostPtr( ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<double> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace > AA_Internal_dev;

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
#endif

#ifdef DREAL
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  inline
  void FactorSolve_devPtr( ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<double *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  inline
  void FactorSolve_hostPtr( ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<double *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace > AA_Internal_dev;

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
#endif

#ifdef SCPLX
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  inline
  void FactorSolve_devPtr( ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<float> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  inline
  void FactorSolve_hostPtr( ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<float> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace > AA_Internal_dev;

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
#endif

#ifdef SREAL
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  inline
  void FactorSolve_devPtr( ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<float *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  inline
  void FactorSolve_hostPtr( ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<float *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace > AA_Internal_dev;

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
#endif

}

