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

#ifndef __ADELUS_HPP__
#define __ADELUS_HPP__

#pragma once

#include <Kokkos_Core.hpp>
#include <Adelus_global_comm.hpp>
#include <Adelus_defines.h>
#include <Adelus_distribute.hpp>
#include <Adelus_xlu_solve.hpp>
#include <Adelus_x_factor.hpp>
#include <Adelus_x_solve.hpp>

#include <mpi.h>

// Adelus: provides the functionality to interface to a dense LU solver

namespace Adelus {

  /// Adelus GetDistribution
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

  /// Adelus GetDistribution (old interface)
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

    distmat_( get_global_comm(),
              *nprocs_row_,
              *number_of_unknowns,
              *nrhs_,
              *my_rows_,
              *my_cols_,
              *my_first_row_,
              *my_first_col_,
              *my_rhs_,
              *my_row,
              *my_col );

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

  /// Adelus FactorSolve (old interface)
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

    MPI_Comm_rank(get_global_comm(), &rank);

#ifdef PRINT_STATUS
    printf("FactorSolve (Kokkos View interface) in rank %d -- my_rows %u , my_cols %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, *matrix_size, *num_procsr, *num_rhs);
#endif

    using value_type = typename ZDView::value_type;
    using execution_space = typename ZDView::device_type::execution_space;
    using memory_space    = typename ZDView::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
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
                         Kokkos::HIPSpace,
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
                         Kokkos::HIPSpace> AA_Internal_dev;
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

  /// Adelus FactorSolve_hostPtr (old interface)
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

    MPI_Comm_rank(get_global_comm(), &rank) ;

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
                         Kokkos::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    using value_type = typename AA_Internal_dev::value_type;
    using execution_space = typename AA_Internal_dev::device_type::execution_space;
    using memory_space    = typename AA_Internal_dev::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i_dev, secs );

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    using value_type = typename AA_Internal::value_type;
    using execution_space = typename AA_Internal::device_type::execution_space;
    using memory_space    = typename AA_Internal::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i, secs );
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
                         Kokkos::HIPSpace,
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
                         Kokkos::HIPSpace> AA_Internal_dev;
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

  /// Adelus FactorSolve_hostPtr (old interface)
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

    MPI_Comm_rank(get_global_comm(), &rank) ;

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
                         Kokkos::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    using value_type = typename AA_Internal_dev::value_type;
    using execution_space = typename AA_Internal_dev::device_type::execution_space;
    using memory_space    = typename AA_Internal_dev::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i_dev, secs );

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    using value_type = typename AA_Internal::value_type;
    using execution_space = typename AA_Internal::device_type::execution_space;
    using memory_space    = typename AA_Internal::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i, secs );
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
                         Kokkos::HIPSpace,
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
                         Kokkos::HIPSpace> AA_Internal_dev;
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

  /// Adelus FactorSolve_hostPtr (old interface)
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

    MPI_Comm_rank(get_global_comm(), &rank) ;

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
                         Kokkos::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    using value_type = typename AA_Internal_dev::value_type;
    using execution_space = typename AA_Internal_dev::device_type::execution_space;
    using memory_space    = typename AA_Internal_dev::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i_dev, secs );

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float complex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    using value_type = typename AA_Internal::value_type;
    using execution_space = typename AA_Internal::device_type::execution_space;
    using memory_space    = typename AA_Internal::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i, secs );
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
                         Kokkos::HIPSpace,
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
                         Kokkos::HIPSpace> AA_Internal_dev;
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

  /// Adelus FactorSolve_hostPtr (old interface)
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

    MPI_Comm_rank(get_global_comm(), &rank) ;

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
                         Kokkos::HIPSpace> AA_Internal_dev;
#endif

    AA_Internal_dev AA_i_dev( "AA_i_dev", my_rows_, my_cols_ + my_rhs_ + 6 );

#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with CUDA solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    Kokkos::deep_copy( AA_i_dev, AA_i );

    using value_type = typename AA_Internal_dev::value_type;
    using execution_space = typename AA_Internal_dev::device_type::execution_space;
    using memory_space    = typename AA_Internal_dev::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i_dev, secs );

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
#ifdef PRINT_STATUS
    printf("FactorSolve_hostPtr with host solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);
#endif

    using value_type = typename AA_Internal::value_type;
    using execution_space = typename AA_Internal::device_type::execution_space;
    using memory_space    = typename AA_Internal::device_type::memory_space;

    Adelus::AdelusHandle<value_type, execution_space, memory_space>
      ahandle(0, get_global_comm(), *matrix_size, *num_procsr, *num_rhs );
    lusolve_( ahandle, AA_i, secs );
#endif
    }
  }
#endif

}//namespace Adelus

#endif
