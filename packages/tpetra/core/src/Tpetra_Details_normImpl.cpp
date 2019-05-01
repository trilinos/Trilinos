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

#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)

#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Tpetra_Details_normImpl_decl.hpp"
#include "Tpetra_Details_normImpl_def.hpp"
#include "TpetraCore_ETIHelperMacros.h"


#if defined(HAVE_TPETRA_INST_INT_INT)
// don't need to do anything; Scalar=int is already added
# define TPETRA_DETAILS_NORMIMPL_INSTANT_INT( NODE )
#else
# define TPETRA_DETAILS_NORMIMPL_INSTANT_INT( NODE ) \
  TPETRA_DETAILS_NORMIMPL_INSTANT( int, NODE )
#endif

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SN_NO_ORDINAL_SCALAR( TPETRA_DETAILS_NORMIMPL_INSTANT )
  TPETRA_INSTANTIATE_GN( TPETRA_DETAILS_NORMIMPL_INSTANT )
#ifndef HAVE_TPETRA_INST_INT_INT
  TPETRA_INSTANTIATE_N( TPETRA_DETAILS_NORMIMPL_INSTANT_INT )
#endif // NOT HAVE_TPETRA_INST_INT_INT

#ifdef HAVE_TPETRA_INST_CUDA

  using cuda_host_mirror_device_type =
    Kokkos::Device<Kokkos::DefaultHostExecutionSpace,
                   Kokkos::CudaUVMSpace>;

#define TPETRA_DETAILS_NORMIMPL_INSTANT_CUDAHOSTMIRROR( S ) \
  TPETRA_DETAILS_NORMIMPL_INSTANT( S, cuda_host_mirror_device_type )

  TPETRA_INSTANTIATE_S_NO_ORDINAL_SCALAR( TPETRA_DETAILS_NORMIMPL_INSTANT_CUDAHOSTMIRROR )
  TPETRA_INSTANTIATE_G( TPETRA_DETAILS_NORMIMPL_INSTANT_CUDAHOSTMIRROR )
#if ! defined(HAVE_TPETRA_INST_INT_INT)
  TPETRA_DETAILS_NORMIMPL_INSTANT_CUDAHOSTMIRROR( int )
#endif

#endif // HAVE_TPETRA_INST_CUDA

} // namespace Tpetra

#endif // Whether we should build this specialization
