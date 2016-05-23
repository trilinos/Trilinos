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

// Include here just the bare minimum needed for the outer #if.
// Only include more once we know we need the test.
#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_OPENMP)

#include "Tpetra_Test_CrsMatrix_WithGraph.hpp"

namespace Tpetra {
namespace Test {

//
// INSTANTIATIONS
//

TPETRA_ETI_MANGLING_TYPEDEFS()

// Declare a colon- and comma-free typedef, to avoid macro issues.
typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node_type;

#define UNIT_TEST_GROUP_OPENMP( SCALAR, LO, GO ) \
  UNIT_TEST_GROUP( SCALAR, LO, GO, openmp_node_type )

TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_OPENMP )

} // namespace Test
} // namespace Tpetra

#endif // defined(HAVE_TPETRA_OPENMP)

