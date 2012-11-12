/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#include "Tpetra_CrsMatrixSolveOp.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "Tpetra_CrsMatrixSolveOp_def.hpp"
// need this to instantiate CrsMatrix::solve()
#include "Tpetra_CrsMatrix_def.hpp"

#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOSCLASSIC_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
#  include <Kokkos_OpenMPNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE)
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

#define INST_SER_GO_SCALAR(GO,SCALAR) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(SCALAR,SCALAR,int,GO,Kokkos::SerialNode)
#ifdef HAVE_KOKKOSCLASSIC_TBB
#define INST_TBB_GO_SCALAR(GO,SCALAR) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(SCALAR,SCALAR,int,GO,Kokkos::TBBNode)
#else
#define INST_TBB_GO_SCALAR(GO,SCALAR)
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#define INST_OMP_GO_SCALAR(GO,SCALAR) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(SCALAR,SCALAR,int,GO,Kokkos::OpenMPNode)
#else
#define INST_OMP_GO_SCALAR(GO,SCALAR)
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#define INST_TPI_GO_SCALAR(GO,SCALAR) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(SCALAR,SCALAR,int,GO,Kokkos::TPINode)
#else
#define INST_TPI_GO_SCALAR(GO,SCALAR)
#endif

#define INST_SER_GO_SCALAR_SCALAR(GO,S1,S2) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(S1,S2,int,GO,Kokkos::SerialNode)
#ifdef HAVE_KOKKOSCLASSIC_TBB
#define INST_TBB_GO_SCALAR_SCALAR(GO,S1,S2) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(S1,S2,int,GO,Kokkos::TBBNode)
#else
#define INST_TBB_GO_SCALAR_SCALAR(GO,S1,S2)
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#define INST_OMP_GO_SCALAR_SCALAR(GO,S1,S2) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(S1,S2,int,GO,Kokkos::OpenMPNode)
#else
#define INST_OMP_GO_SCALAR_SCALAR(GO,S1,S2)
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#define INST_TPI_GO_SCALAR_SCALAR(GO,S1,S2) TPETRA_CRSMATRIX_SOLVEOP_INSTANT(S1,S2,int,GO,Kokkos::TPINode)
#else
#define INST_TPI_GO_SCALAR_SCALAR(GO,S1,S2)
#endif

#define INST_CPU_GO_SCALAR(GO,SCALAR) \
  INST_SER_GO_SCALAR(GO,SCALAR) \
  INST_TBB_GO_SCALAR(GO,SCALAR) \
  INST_OMP_GO_SCALAR(GO,SCALAR) \
  INST_TPI_GO_SCALAR(GO,SCALAR)

#ifdef HAVE_TPETRA_INST_INT_LONG
#define INST_CPU_SCALAR(SCALAR) \
        INST_CPU_GO_SCALAR(int,SCALAR) \
        INST_CPU_GO_SCALAR(long,SCALAR)
#else
#define INST_CPU_SCALAR(SCALAR) \
        INST_CPU_GO_SCALAR(int,SCALAR)
#endif

namespace Tpetra {

  //
  // instantiate all single-scalar implementations; these are needed internally by CrsMatrix
  //

#if defined(HAVE_TPETRA_INST_DOUBLE)
  INST_CPU_SCALAR(double)
#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,double,int,int,Kokkos::ThrustGPUNode)
#ifdef HAVE_TPETRA_INST_INT_LONG
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,double,int,long,Kokkos::ThrustGPUNode)
#endif
#endif
#endif // double

#if defined(HAVE_TPETRA_INST_FLOAT)
  INST_CPU_SCALAR(float)
#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,int,Kokkos::ThrustGPUNode)
#ifdef HAVE_TPETRA_INST_INT_LONG
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,long,Kokkos::ThrustGPUNode)
#endif
#endif 
#endif // float

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  INST_CPU_SCALAR(std::complex<double>)
// not yet supported by cusparse
// #if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE)
//     TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,std::complex<double>,int,int,Kokkos::ThrustGPUNode)
// #endif
#endif // complex double

#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  INST_CPU_SCALAR(std::complex<float>)
// not yet supported by cusparse
// #if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT)
//     TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,std::complex<float>,int,int,Kokkos::ThrustGPUNode)
// #endif
#endif // complex float

  // cross scalar applications
  // double-float and float-double
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
  INST_TPI_GO_SCALAR_SCALAR(int,float,double)
  INST_TPI_GO_SCALAR_SCALAR(int,double,float)
#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,float,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,double,int,int,Kokkos::ThrustGPUNode)
#ifdef HAVE_TPETRA_INST_INT_LONG
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,float,int,long,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,double,int,long,Kokkos::ThrustGPUNode)
#endif
#endif
#endif

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
