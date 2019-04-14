// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "TpetraCore_config.h"
#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TPETRA_INST_CUDA)

// mfh 13 Apr 2019: I can't just use Stokhos' system for Tpetra ETI,
// because that only works for Tpetra ETI macros that take four
// arguments (SCALAR, LO, GO, NODE).  The macro for ETI of
// Tpetra::Details::normImpl only takes two arguments: (SCALAR, NODE).

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "Tpetra_Details_normImpl_decl.hpp"
#include "Tpetra_Details_normImpl_def.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#define TPETRA_DETAILS_NORMIMPL_INSTANT_N( NODE ) \
  using scalar_type_##NODE = Sacado::MP::Vector<Stokhos::StaticFixedStorage<int, double, 16, NODE::execution_space> >; \
  TPETRA_DETAILS_NORMIMPL_INSTANT( scalar_type_##NODE, NODE )

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_DETAILS_NORMIMPL_INSTANT_N( Kokkos_Compat_KokkosCudaWrapperNode )

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION && HAVE_TPETRA_INST_CUDA
