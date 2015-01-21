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

#include "Tpetra_RowMatrixTransposer.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_RowMatrixTransposer_def.hpp"

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // mfh 21 Jul 2013: Unfortunately, the Tpetra explicit template
  // instantiation (ETI) system doesn't currently let us control the
  // Scalar types for a specific Node; it attempts to instantiate for
  // all supported Scalar types.  This is a problem for ThrustGPUNode
  // with CUSPARSEOps, since the latter only supports Scalar = void,
  // float, or double.  I deal with this as follows:
  //
  // 1. I only instantiate all supported Scalar, LocalOrdinal, and
  //    GlobalOrdinal combinations for non-GPU Node types.
  //
  // 2. With Node=ThrustGPUNode, I manually instantiate
  //    RowMatrixTransposer only for the supported Scalar types.  In
  //    this case, that's Scalar = float and double.  (Scalar = void
  //    is only meaningful for CUSPARSEOps; it's a trick to deal with
  //    the fact that Tpetra::CrsGraph has to refer to
  //    CUSPARSEOps<Scalar, ...>, but doesn't have a sensible Scalar
  //    type other than void.)

  // Let Tpetra's ETI system handle non-GPU types.
  TPETRA_INSTANTIATE_SLGN_NOGPU(TPETRA_ROWMATRIXTRANSPOSER_INSTANT)

  // mfh 21 Jul 2013: Since we're rolling explicit instantiation by
  // hand here, we don't have to use the typedefs.  (The typedefs work
  // around the fact that macros do not deal well with commas or
  // spaces in their arguments.  We're not using macros to instantiate
  // the classes here, so we don't need the typedefs.)

  // FIXME (mfh 21 Jul 2013) The instantiations below do not include
  // GlobalOrdinal = unsigned long (CMake BOOL option
  // Tpetra_INST_UNSIGNED_LONG) and other combinations of LocalOrdinal
  // and GlobalOrdinal types enabled by Tpetra's ETI system.  However,
  // it does avoid link errors for common cases.

  // FIXME (mfh 22 Jul 2013) There must be a more concise way to ask
  // if ThrustGPUNode with CUSPARSEOps (the default sparse kernels for
  // ThrustGPUNode) works.  However, this works for now.
#if defined(HAVE_KOKKOSCLASSIC_CUDA) && defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_KOKKOSCLASSIC_CUSPARSE)

  //
  // LocalOrdinal = int, GlobalOrdinal = int
  //

#  if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
  template class RowMatrixTransposer<float, int, int, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  template class RowMatrixTransposer<double, int, int, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT)
  template class RowMatrixTransposer<std::complex<float>, int, int, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE)
  template class RowMatrixTransposer<std::complex<double>, int, int, KokkosClassic::ThrustGPUNode>;
#  endif

  //
  // LocalOrdinal = int, GlobalOrdinal = unsigned int
  //

#  if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
  template class RowMatrixTransposer<float, int, unsignedint, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  template class RowMatrixTransposer<double, int, unsignedint, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT)
  template class RowMatrixTransposer<std::complex<float>, int, unsignedint, KokkosClassic::ThrustGPUNode>;
#  endif
#  if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE)
  template class RowMatrixTransposer<std::complex<double>, int, unsignedint, KokkosClassic::ThrustGPUNode>;
#  endif

  //
  // LocalOrdinal = int, GlobalOrdinal = long
  //

#  ifdef HAVE_TPETRA_INST_INT_LONG
#    if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT)
  template class RowMatrixTransposer<float, int, long, KokkosClassic::ThrustGPUNode>;
#    endif
#    if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  template class RowMatrixTransposer<double, int, long, KokkosClassic::ThrustGPUNode>;
#    endif
#    if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT)
  template class RowMatrixTransposer<std::complex<float>, int, long, KokkosClassic::ThrustGPUNode>;
#    endif
#    if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE)
  template class RowMatrixTransposer<std::complex<double>, int, long, KokkosClassic::ThrustGPUNode>;
#    endif
#  endif // HAVE_TPETRA_INST_INT_LONG

#endif // defined(HAVE_KOKKOSCLASSIC_CUDA) && defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_KOKKOSCLASSIC_CUSPARSE)

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
