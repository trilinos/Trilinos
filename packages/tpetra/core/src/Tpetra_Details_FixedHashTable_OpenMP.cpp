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

#include "Tpetra_Details_FixedHashTable_decl.hpp"

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(KOKKOS_ENABLE_OPENMP)

#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_Details_FixedHashTable_def.hpp"

namespace Tpetra {
namespace Details {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> openmp_device_type;

#define TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_OPENMP( LO, GO ) \
  TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT( LO, GO, openmp_device_type )

  TPETRA_INSTANTIATE_LG( TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_OPENMP )

  // mfh 26 Sep 2015: Make sure that the {KeyType = LO, ValueType =
  // LO} and {KeyType = int, ValueType = LO} specializations get
  // built, since Directory needs them.

  // KeyType = int doesn't get built if GO = int is disabled.
#ifndef HAVE_TPETRA_INST_INT_INT
#  define TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_OPENMP_INT( LO ) \
  TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT( LO, int, openmp_device_type )

  TPETRA_INSTANTIATE_L( TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_OPENMP_INT )

  // FIXME (mfh 26 Sep 2015) Once it becomes possible to disable LO =
  // int, add an instantiation here for {KeyType = LO, ValueType =
  // LO}.  However, this case has likely already been covered above,
  // since disabling LO = int leaves LO = GO = long or long long as
  // the most reasonable options.  (It would be silly to use LO = long
  // and GO = long long if sizeof(long) = sizeof(long long).)

#endif // HAVE_TPETRA_INST_INT_INT

} // namespace Details
} // namespace Tpetra

#endif // defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(KOKKOS_ENABLE_OPENMP)
