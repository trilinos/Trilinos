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
  // GlobalOrdinal type.  This is why we only include three template
  // parameters: Scalar (ST) and LocalOrdinal (LO).  At some point, it
  // would make sense to include Node as well, but for now we omit it,
  // since LittleBlock and LittleVector as yet live in host memory.

  // Test small dense block LU factorization and solve, with an easy
  // problem (the identity matrix).
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, SolveIdentity, ST, LO )
  {
    typedef Tpetra::Experimental::LittleBlock<ST, LO> block_type;
    typedef Tpetra::Experimental::LittleVector<ST, LO> vec_type;
    const ST zero = static_cast<ST> (0.0);
    const ST one = static_cast<ST> (1.0);
    const LO minBlockSize = 1; // 1x1 "blocks" should also work
    const LO maxBlockSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<ST> blockPool (maxBlockSize * maxBlockSize);
    // Memory pool for the LittleVector instances (x and b).
    Teuchos::Array<ST> vecPool (maxBlockSize * 2);
    // Memory pool for the pivot vector.
    Teuchos::Array<int> ipivPool (maxBlockSize);

    for (LO blockSize = minBlockSize; blockSize <= maxBlockSize; ++blockSize) {
      block_type A (blockPool (0, blockSize*blockSize).getRawPtr (),
                    blockSize, 1, blockSize);
      Teuchos::ArrayView<ST> x_view = vecPool (0, blockSize);
      vec_type x (x_view.getRawPtr (), blockSize, 1);
      Teuchos::ArrayView<ST> b_view = vecPool (blockSize, blockSize);
      vec_type b (b_view.getRawPtr (), blockSize, 1);
      Teuchos::ArrayView<int> ipiv = ipivPool (0, blockSize);

      A.fill (zero);
      for (LO i = 0; i < blockSize; ++i) {
        A(i,i) = one;
        b(i) = static_cast<ST> (i + 1);
        x(i) = b(i); // copy of right-hand side on input
        ipiv[i] = 0;
      }

      int info = 0;
      A.factorize (ipiv.getRawPtr (), info);

      TEST_EQUALITY_CONST( info, 0 );
      if (info == 0) {
        A.solve (x, ipiv.getRawPtr ());
      }

      // Re-fill b, in case A.solve brokenly clobbered it.
      for (LO i = 0; i < blockSize; ++i) {
        b(i) = static_cast<ST> (i + 1);
      }

      TEST_COMPARE_ARRAYS( x_view, b_view );
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, SolveIdentity, SCALAR, LOCAL_ORDINAL )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SL_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)


