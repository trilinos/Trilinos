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

#include "Tpetra_CrsMatrix.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

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

#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"

namespace Tpetra {

//
// Explicit instantiations for LocalOrdinal = int, GlobalOrdinal = int
//
#if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,double,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,double,int,int,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::TBBNode)
#endif // HAVE_KOKKOSCLASSIC_TBB

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::TPINode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,double,int,int,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::TPINode)
#endif // HAVE_KOKKOSCLASSIC_THREADPOOL

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,double,int,int,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::OpenMPNode)
#endif // HAVE_KOKKOSCLASSIC_OPENMP

#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,double,int,int,Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, int, Kokkos::ThrustGPUNode)
#endif // HAVE_KOKKOSCLASSIC_CUSPARSE && HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
#endif // DOUBLE

#if defined(HAVE_TPETRA_INST_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,float,int,int,Kokkos::SerialNode)
#  if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,double,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,float,int,int,Kokkos::SerialNode)
#  endif // if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::TBBNode)
#  if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,double,int,int,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,float,int,int,Kokkos::TBBNode)
#  endif // if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::TPINode)
#  if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,double,int,int,Kokkos::TPINode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,float,int,int,Kokkos::TPINode)
#  endif // if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::OpenMPNode)
#  if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,double,int,int,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,float,int,int,Kokkos::OpenMPNode)
#  endif // if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::OpenMPNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::ThrustGPUNode)
#  if defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(float,double,int,int,Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,float,int,int,Kokkos::ThrustGPUNode)
#  endif // if defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, int, Kokkos::ThrustGPUNode)
#endif
#endif // FLOAT

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,int,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,int,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,int,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::OpenMPNode)
#endif
//
// No std::complex<double> on GPU support for now.
//
//#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
//  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,int,Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, int, Kokkos::ThrustGPUNode)
//#endif
#endif // COMPLEX DOUBLE

#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,int,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,int,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,int,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::OpenMPNode)
#endif
//
// No std::complex<float> on GPU support for now.
//
//#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
//  TPETRA_CRSMATRIX_INSTANT(double,int,int,Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, int, Kokkos::ThrustGPUNode)
//#endif
#endif // COMPLEX FLOAT

//
// Explicit instantiations for LocalOrdinal = int, GlobalOrdinal = long
//
#ifdef HAVE_TPETRA_INST_INT_LONG
#if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::TBBNode)
#endif // HAVE_KOKKOSCLASSIC_TBB

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::TPINode)
#endif // HAVE_KOKKOSCLASSIC_THREADPOOL

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::OpenMPNode)
#endif // HAVE_KOKKOSCLASSIC_OPENMP

#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(double, int, long, Kokkos::ThrustGPUNode)
#endif // HAVE_KOKKOSCLASSIC_CUSPARSE && HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
#endif // DOUBLE

#if defined(HAVE_TPETRA_INST_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(float,int,long,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(float,int,long,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(float,int,long,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(float,int,long,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::OpenMPNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(float,int,long,Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::ThrustGPUNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(float, int, long, Kokkos::ThrustGPUNode)
#endif
#endif // FLOAT

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,long,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,long,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,long,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,long,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::OpenMPNode)
#endif
//
// No std::complex<double> on GPU support for now.
//
//#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
//  TPETRA_CRSMATRIX_INSTANT(std::complex<double>,int,long,Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<double>, int, long, Kokkos::ThrustGPUNode)
//#endif
#endif // COMPLEX DOUBLE

#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,long,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::SerialNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,long,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::TBBNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,long,Kokkos::TPINode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::TPINode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_CRSMATRIX_INSTANT(std::complex<float>,int,long,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::OpenMPNode)
#endif
//
// No std::complex<float> on GPU support for now.
//
//#if defined(HAVE_KOKKOSCLASSIC_CUSPARSE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
//  TPETRA_CRSMATRIX_INSTANT(double,int,long,Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::ThrustGPUNode)
//  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(std::complex<float>, int, long, Kokkos::ThrustGPUNode)
//#endif
#endif // COMPLEX FLOAT
#endif


} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
