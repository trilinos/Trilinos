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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Experimental_BlockView.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using std::endl;
  typedef Teuchos::Array<int>::size_type size_type;

  //
  // UNIT TESTS
  //

  // The "little" blocks and vectors do not depend on Tpetra's
  // GlobalOrdinal type.  At some point, it would make sense to
  // include the Kokkos device type as well, but for now we omit it,
  // since LittleBlock and LittleVector as yet live in host memory.

  // Make sure that LittleBlock and LittleVector have the correct
  // typedefs and other stuff.
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LittleBlockVector, Typedefs, ST )
  {
    typedef Tpetra::Experimental::LittleBlock<ST> blk_type;
    typedef Tpetra::Experimental::LittleVector<ST> vec_type;

    static_assert (std::is_same<typename blk_type::array_layout, Kokkos::LayoutLeft>::value ||
                   std::is_same<typename blk_type::array_layout, Kokkos::LayoutRight>::value ||
                   std::is_same<typename blk_type::array_layout, Kokkos::LayoutStride>::value,
                   "LittleBlock needs to have the array_layout typedef with one of the following "
                   "values: LayoutLeft, LayoutRight, or LayoutStride (all in the Kokkos namespace).");
    static_assert (blk_type::rank == 2, "LittleBlock must have rank 2.");
    static_assert (vec_type::rank == 1, "LittleVector must have rank 1.");
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LittleBlockVector, Typedefs, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_GROUP )

} // namespace (anonymous)


