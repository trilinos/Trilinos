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
#include "Tpetra_BlockView.hpp"
#include "Teuchos_BLAS.hpp"
#include <vector>

namespace {
  using std::endl;
  typedef int LO;

  template<class IST, class LayoutType>
  void
  testGEMV (Teuchos::FancyOStream& out,
            bool& success,
            const LO maxBlkSize)
  {
    using Tpetra::GEMV;
    typedef Kokkos::View<IST**, LayoutType, Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> blk_type;
    typedef Kokkos::View<IST*, LayoutType, Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> vec_type;
    typedef Kokkos::Details::ArithTraits<IST> KAT;

    Teuchos::OSTab tab2 (out);
    const IST zero = KAT::zero ();
    const IST one = KAT::one ();
    //const IST two = one + one;
    // Temporary space for the blocks.  For now, we only exercise
    // HostSpace blocks.
    std::vector<IST> blkMem (maxBlkSize * maxBlkSize);
    std::vector<IST> vecMem (3 * maxBlkSize);

    for (LO blkSize = 1; blkSize < maxBlkSize; ++blkSize) {
      // For now, we only test square blocks.  Change the definitions
      // of numRows and numCols below, and change the temporary
      // allocations above, if you wish to generalize this test to
      // non-square blocks.
      const LO numRows = blkSize;
      const LO numCols = blkSize;

      blk_type A (&blkMem[0], numRows, numCols);
      vec_type x (&vecMem[0], numCols);
      vec_type y (&vecMem[blkSize], numRows);
      vec_type y_expected (&vecMem[2*blkSize], numRows);
      for (LO i = 0; i < numRows; ++i) {
        y_expected(i) = zero; // to be revised below
      }

      IST curVal = one;
      for (LO i = 0; i < numRows; ++i) {
        for (LO j = 0; j < numCols; ++j) {
          A(i,j) = curVal;
          y_expected(i) += curVal;
          curVal = curVal + one;
        }
        y(i) = zero;
      }
      for (LO j = 0; j < numCols; ++j) {
        x(j) = one;
      }

      GEMV (one, A, x, y);
      for (LO i = 0; i < numRows; ++i) {
        TEST_EQUALITY( y(i), y_expected(i) );
      }

      GEMV (-one, A, x, y);
      for (LO i = 0; i < numRows; ++i) {
        TEST_EQUALITY( y(i), zero );
      }
    }
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( LittleBlockOps, GEMV, ST )
  {
    typedef Kokkos::Details::ArithTraits<ST> KAT;
    typedef typename KAT::val_type IST;

    out << "Test GEMV" << endl;
    Teuchos::OSTab tab1 (out);

    // Test blocks of dimensions 1 x 1, 2 x 2, ...,
    // maxBlkSize x maxBlkSize.
    const LO maxBlkSize = 31;

    // Test LayoutRight (row-major) Views.
    {
      out << "Test LayoutRight (row-major) Views" << endl;
      Teuchos::OSTab tab2 (out);
      testGEMV<IST, Kokkos::LayoutRight> (out, success, maxBlkSize);
    }

    // Test LayoutLeft (column-major) Views.
    {
      out << "Test LayoutLeft (column-major) Views" << endl;
      Teuchos::OSTab tab2 (out);
      testGEMV<IST, Kokkos::LayoutLeft> (out, success, maxBlkSize);
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( LittleBlockOps, GEMV, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_S( UNIT_TEST_GROUP )

} // namespace (anonymous)


