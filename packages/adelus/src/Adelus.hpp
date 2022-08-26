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

#ifndef __ADELUS_HPP__
#define __ADELUS_HPP__

#pragma once

#include <Kokkos_Core.hpp>
#include <Adelus_defines.h>
#include <Adelus_distribute.hpp>
#include <Adelus_xlu_solve.hpp>
#include <Adelus_x_factor.hpp>
#include <Adelus_x_solve.hpp>

#include <mpi.h>

// Adelus: provides the functionality to interface to a dense LU solver

namespace Adelus {

  /// Adelus GetDistirbution
  /// Gives the distribution information that is required by the dense solver
  
  /// \param comm (In)               - communicator that Adelus runs on
  /// \param nprocs_row (In)         - number of processors for a row
  /// \param number_of_unknowns (In) - order of the dense matrix
  /// \param nrhs (In)               - number of right hand sides
  /// \param my_rows (Out)           - number of rows of the matrix on this processor
  /// \param my_cols (Out)           - number of columns of the matrix on this processor
  /// \param my_first_row (Out)      - first (global) row number on this processor (array starts at index 1)
  /// \param my_first_col (Out)      - first (global) column number on this processor (array starts at index 1)
  /// \param my_rhs (Out)            - number of right hand sides on this processor
  /// \param my_row (Out)            - row number in processor mesh, 0 to the  number of processors for a column -1
  /// \param my_col (Out)            - column number in processor mesh, 0 to the  number of processors for a row -1
    
  inline
  int GetDistribution( MPI_Comm comm,
                       const int nprocs_row,
                       const int number_of_unknowns,
                       const int nrhs,
                       int& my_rows,
                       int& my_cols,
                       int& my_first_row,
                       int& my_first_col,
                       int& my_rhs,
                       int& my_row,
                       int& my_col ) {
    // This function echoes the multiprocessor distribution of the matrix

    distmat_(comm,
             nprocs_row,
             number_of_unknowns,
             nrhs,
             my_rows,
             my_cols,
             my_first_row,
             my_first_col,
             my_rhs,
             my_row,
             my_col);

    return(0);

  }

  /// Adelus FactorSolve
  /// Factors and solves the dense matrix

  /// \param ahandle (In)     -- handle that contains metadata needed by the Adelus solver
  /// \param AA (InOut)       -- Kokkos View that has the matrix and rhs packed in this processor
  ///                            (Note: matrix and rhs are overwritten)
  /// \param secs (Out)       -- factor and solve time in seconds
    
  template<class HandleType, class ZRHSViewType>
  inline
  void FactorSolve( HandleType& ahandle,
                    ZRHSViewType& AA,
                    double* secs ) {
	
#ifdef PRINT_STATUS
    printf("FactorSolve (Kokkos View interface) in rank %d\n", ahandle.get_myrank());
#endif

    lusolve_(ahandle, AA, secs);

  }

  /// Adelus Factor
  /// Factors the dense matrix for later solve

  /// \param ahandle (In)     -- handle that contains metadata needed by the Adelus solver
  /// \param AA (InOut)       -- Kokkos View that has the matrix in this processor (Note: matrix is overwritten)
  /// \param permute (In)     -- Kokkos View that has the global pivot vector
  /// \param secs (Out)       -- factor and solve time in seconds

  template<class HandleType, class ZViewType, class PViewType>
  inline
  void Factor( HandleType& ahandle,
               ZViewType& AA,
               PViewType& permute,
               double* secs ) {

#ifdef PRINT_STATUS
    printf("Factor (Kokkos View interface) in rank %d\n", ahandle.get_myrank());
#endif

    lu_(ahandle, AA, permute, secs);

  }

  /// Adelus Solve
  /// Solves the previously factored dense matrix for provided RHS

  /// \param ahandle (In)     -- handle that contains metadata needed by the Adelus solver
  /// \param AA (In)          -- Kokkos View that has the LU-factorized matrix
  /// \param BB (InOut)       -- Kokkos View that has the rhs and solution (Note: rhs are overwritten)
  /// \param permute (In)     -- Kokkos View that has the global pivot vector
  /// \param secs (Out)       -- factor and solve time in seconds

  template<class HandleType, class ZViewType,
           class RHSViewType, class PViewType>
  inline
  void Solve( HandleType& ahandle,
              ZViewType& AA,
              RHSViewType& BB,
              PViewType& permute,
              double* secs ) {

#ifdef PRINT_STATUS
    printf("Solve (Kokkos View interface) in rank %d\n", ahandle.get_myrank());
#endif

    solve_(ahandle, AA, BB, permute, secs);

  }

#ifdef ZCPLX
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  template<class HandleType>
  inline
  void FactorSolve_devPtr( HandleType& ahandle,
                           ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace,
#else
                         Kokkos::Experimental::HIPSpace,
#endif
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<double> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  template<class HandleType>
  inline
  void FactorSolve_hostPtr( HandleType& ahandle,
                            ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<double> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) 
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace> AA_Internal_dev;
#else
                         Kokkos::Experimental::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(ahandle, AA_i_dev, secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }
#endif

#ifdef DREAL
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  template<class HandleType>
  inline
  void FactorSolve_devPtr( HandleType& ahandle,
                           ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace,
#else
                         Kokkos::Experimental::HIPSpace,
#endif
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<double *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  template<class HandleType>
  inline
  void FactorSolve_hostPtr( HandleType& ahandle,
                            ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<double *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) 
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace> AA_Internal_dev;
#else
                         Kokkos::Experimental::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(ahandle, AA_i_dev, secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }
#endif

#ifdef SCPLX
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  template<class HandleType>
  inline
  void FactorSolve_devPtr( HandleType& ahandle,
                           ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace,
#else
                         Kokkos::Experimental::HIPSpace,
#endif
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<float> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  template<class HandleType>
  inline
  void FactorSolve_hostPtr( HandleType& ahandle,
                            ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<float> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace> AA_Internal_dev;
#else
                         Kokkos::Experimental::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(ahandle, AA_i_dev, secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }
#endif

#ifdef SREAL
  /// Adelus FactorSolve_devPtr
  /// Matrix and rhs are packed and passed as device pointer

  template<class HandleType>
  inline
  void FactorSolve_devPtr( HandleType& ahandle,
                           ADELUS_DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace,
#else
                         Kokkos::Experimental::HIPSpace,
#endif
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<float *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#ifdef PRINT_STATUS
    printf("FactorSolve_devPtr (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }

  /// Adelus FactorSolve_hostPtr
  /// Matrix and rhs are packed and passed as host pointer

  template<class HandleType>
  inline
  void FactorSolve_hostPtr( HandleType& ahandle,
                            ADELUS_DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs ) {

    { // Note: To avoid segmentation fault when FactorSolve is called multiple times with the unmanaged View, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<float *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
#ifdef KOKKOS_ENABLE_CUDA
                         Kokkos::CudaSpace> AA_Internal_dev;
#else
                         Kokkos::Experimental::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(ahandle, AA_i_dev, secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", ahandle.get_myrank(), my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    lusolve_(ahandle, AA_i, secs);
#endif
    }
  }
#endif

}//namespace Adelus

#endif
