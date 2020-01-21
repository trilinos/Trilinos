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

#pragma once

#include <Kokkos_View.hpp>
#include <defines.h>
#include <xlu_solve.hpp>
#include <distribute.hpp>
#include <mpi.h>

//! Pliris: An Obect-Oriented Interface to a Dense LU Solver.
/*! The Pliris class : Provides the functionality to interface to a dense LU
*/


namespace Pliris {



  //! Pliris RHS Set
  /*! Associates an already defined Epetra_MultiVector (or Epetra_Vector) as the right-hand-side of
      the linear system.
   */

  ////int SetRHS(linearAlgebra::VectorBlock *X);

  //! Pliris Matrix Set
  /*! Associates an already defined  VectorBlock as the matrix (column ordered) of
      the linear system.
   */
  ////int SetMatrix( linearAlgebra::Matrix  *A);


  //! Pliris GetDistirbution
  /*! Gives the distribution information that is required by the dense solver
   */

    /*!
     \param nprocs_row(In)   - number of processors for a row
     \param number_of_unknowns(In)  - order of the dense matrix
     \param nrhs(In)  - number of right hand sides
     \param my_rows(Out)  - number of rows of the matrix on this processor
     \param my_cols (Out)  - number of columns of the matrix on this processor
     \param  my_first_row(Out)  - first (global) row number on this processor (array starts at index 1)
     \param  my_first_col (Out)  - first (global) column number on this processor (array starts at index 1)
     \param  my_rhs(Out)  - number of Right hand sides on this processor
     \param my_row(Out)  - row number in processor mesh, 0 to the  number of processors for a column -1
     \param my_col(Out)  - column  number in processor mesh, 0 to the  number of processors for a row -1
    */

  int GetDistribution( int* nprocs_row_,
                       int* number_of_unknowns,
                       int* nrhs_,
                       int* my_rows_,
                       int* my_cols_,
                       int* my_first_row_,
                       int* my_first_col_,
                       int* my_rhs_,
                       int* my_row,
                       int* my_col ){
    // This function echoes the multiprocessor distribution of the matrix

    distmat_(  nprocs_row_,
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

  //! Pliris FactorSolve
  /*! Factors and solves the dense matrix
   */

    /*!
     \param A(InOut) --  Epetra Vector that has the matrix and rhs packed( Note: matrix is overwritten)
     \param my_rows(In) --  number of rows of the matrix on this processor
     \param my_cols(In) --  number of columns of the matrix on this processor
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param num_rhs(In) --  number of right hand sides
     \param secs(Out) -- factor and solve time in seconds
    */


  //! Pliris FactorSolve
  /*! Factors and solves the dense matrix
   */

    /*!
     \param AA(In) --  Epetra Serial Dense Vector that has the matrix and rhs packed
     \param my_rows(In) --  number of rows of the matrix on this processor
     \param my_cols(In) --  number of columns of the matrix on this processor
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param num_rhs(In) --  number of right hand sides
     \param secs(Out) -- factor and solve time in seconds
    */

  template<class ZDView>
  void FactorSolve( ZDView AA,
                    int my_rows_,
                    int my_cols_,
                    int* matrix_size,
                    int* num_procsr,
                    int* num_rhs,
                    double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ; 

    printf("FactorSolve (Kokkos View interface) in rank %d -- my_rows %u , my_cols %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

  }

#ifdef ZCPLX
  void FactorSolve_devPtr( DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<double>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<double> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

    printf("FactorSolve_devPtr (dcomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
  void FactorSolve_hostPtr( DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
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

    printf("FactorSolve_hostPtr with CUDA solve (dcomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
    printf("FactorSolve_hostPtr with host solve (dcomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

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
  void FactorSolve_devPtr( DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<double *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

    printf("FactorSolve_devPtr (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
  void FactorSolve_hostPtr( DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
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

    printf("FactorSolve_hostPtr with CUDA solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
    printf("FactorSolve_hostPtr with host solve (double pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

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
  void FactorSolve_devPtr( DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<Kokkos::complex<float>**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<Kokkos::complex<float> *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

    printf("FactorSolve_devPtr (scomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
  void FactorSolve_hostPtr( DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
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

    printf("FactorSolve_hostPtr with CUDA solve (scomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
    printf("FactorSolve_hostPtr with host solve (scomplex pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

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
  void FactorSolve_devPtr( DATA_TYPE* AA,
                           int my_rows_,
                           int my_cols_,
                           int my_rhs_,
                           int* matrix_size,
                           int* num_procsr,
                           int* num_rhs,
                           double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
#ifdef KOKKOS_ENABLE_CUDA 
    typedef Kokkos::View<float**,
                         Kokkos::LayoutLeft,
                         Kokkos::CudaSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > AA_Internal;

    AA_Internal AA_i(reinterpret_cast<float *>(AA), my_rows_, my_cols_ + my_rhs_ + 6);

    printf("FactorSolve_devPtr (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
  void FactorSolve_hostPtr( DATA_TYPE* AA,
                            int my_rows_,
                            int my_cols_,
                            int my_rhs_,
                            int* matrix_size,
                            int* num_procsr,
                            int* num_rhs,
                            double* secs){
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    { // Note: To avoid segmentation fault when FactorSolve is called mulyiple times with the unmanaged view, it's safest to make sure unmanaged View falls out of scope before freeing its memory.
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

    printf("FactorSolve_hostPtr with CUDA solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    Kokkos::deep_copy( AA_i_dev, AA_i );

    lusolve_(AA_i_dev,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);

    Kokkos::deep_copy( AA_i, AA_i_dev );
#else//OpenMP
    printf("FactorSolve_hostPtr with host solve (float pointer interface) in rank %d -- my_rows %u , my_cols %u, my_rhs %u , matrix_size %u, num_procs_per_row %u, num_rhs %u\n", rank, my_rows_, my_cols_, my_rhs_, *matrix_size, *num_procsr, *num_rhs);

    lusolve_(AA_i,
             matrix_size,
             num_procsr,
             num_rhs,
             secs);
#endif
    }
  }
#endif

  //! Pliris Factor
  /*! Factors the dense matrix
   */

    /*!
     \param A(In) --  Epetra Vector that has the matrix packed
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param permute(In) -- permutation matrix
     \param secs(Out) -- factor and solve time in seconds
    */


   ////nt Factor( linearAlgebra::Matrix* A,
   ////           int* matrix_size,
   ////           int* num_procsr,
   ////           int* permute,
   ////           double* secs);


  //! Pliris Solve
  /*! Solves the previously factored dense matrix
   */

    /*!
     \param permute(In) -- permutation matrix
     \param num_rhs(In) -- factor and solve time in seconds
     \ Note that the matrix has been previously factored by Factor
     \ The RHS has been set by SetRHS(Epetra_MultiVector * B)
     \ On output the result is in the RHS
    */

    ////int Solve(int* permute,
    ////          int* num_rhs);

}

