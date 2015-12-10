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
      std::cerr << "Factor A for blockSize = " << blockSize << std::endl;
      A.factorize (ipiv.getRawPtr (), info);

      TEST_EQUALITY_CONST( info, 0 );
      if (info == 0) {
        std::cerr << "Solve: blockSize = " << blockSize << std::endl;
        A.solve (x, ipiv.getRawPtr ());
      }
      std::cerr << "Done with factor and solve" << std::endl;

      // Re-fill b, in case A.solve brokenly clobbered it.
      for (LO i = 0; i < blockSize; ++i) {
        b(i) = static_cast<ST> (i + 1);
      }

      TEST_COMPARE_ARRAYS( x_view, b_view );
    }
  }


  // Test small dense block QR factorization and solve, with an easy
  // problem (the identity matrix).
  //
  // FIXME (mfh 17 Sep 2015) Right now, this only tests whether
  // calling GEQRF compiles.
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, GEQRF, ST, LO )
  {
    using std::endl;
    typedef Tpetra::Experimental::LittleBlock<ST, LO> block_type;
    const ST zero = static_cast<ST> (0.0);
    const ST one = static_cast<ST> (1.0);
    const LO minBlockSize = 1; // 1x1 "blocks" should also work
    const LO maxBlockSize = 32;

    typename Tpetra::Details::GetLapackType<ST>::lapack_type lapack;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<ST> blockPool (maxBlockSize * maxBlockSize);
    // Memory pool for LAPACK's temporary workspace.  It might need to
    // be resized in the loop below, because LAPACK may want more
    // workspace than just the minimum (the number of columns, which
    // is what the BLAS 2 QR factorization GEQR2 requires).  This must
    // have length at least one, for the workspace query.
    Teuchos::Array<ST> workPool (std::max (1, maxBlockSize));
    // Memory pool for the TAU output array.
    Teuchos::Array<ST> tauPool (maxBlockSize);

    for (LO blockSize = minBlockSize; blockSize <= maxBlockSize; ++blockSize) {
      block_type A (blockPool (0, blockSize*blockSize).getRawPtr (),
                    blockSize, 1, blockSize);

      // Fill A with the identity matrix.
      A.fill (zero);
      for (LO i = 0; i < blockSize; ++i) {
        A(i,i) = one;
      }

      Teuchos::ArrayView<ST> tauView = tauPool (0, blockSize);

      // Workspace query.
      Teuchos::ArrayView<ST> workView = workPool (0, std::max (1, blockSize));
      int lda = blockSize;
      int lwork = -1;
      int info = 0;
      out << "Workspace query" << endl;
      lapack.GEQRF (blockSize, blockSize, A.getRawPtr (), lda,
                    tauView.getRawPtr (),
                    workView.getRawPtr (), lwork, &info);

      TEST_EQUALITY_CONST( info, 0 );
      if (info != 0) {
        continue; // workspace query failed; skip the rest
      }
      lwork = static_cast<int> (Teuchos::ScalarTraits<ST>::real (workView[0]));
      TEST_ASSERT( lwork >= 0 );
      if (lwork < 0) {
        continue; // workspace query failed; skip the rest
      }

      if (workPool.size () < static_cast< decltype (workPool.size ()) > (lwork)) {
        workPool.resize (lwork);
      }
      workView = workPool (0, lwork);
      out << "Workspace size: " << lwork << endl;

      out << "Factor A for blockSize = " << blockSize << endl;
      lapack.GEQRF (blockSize, blockSize, A.getRawPtr (), lda,
                    tauView.getRawPtr (),
                    workView.getRawPtr (), lwork, &info);
      TEST_EQUALITY_CONST( info, 0 );

      out << "Done with factoring A" << endl;
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExpBlockView, SWAP, ST )
  {
    if (! Teuchos::ScalarTraits<ST>::isOrdinal) { // skip integer types
      Teuchos::BLAS<int, ST> blas;
      const int n = 6;
      Teuchos::Array<ST> x (n), y (n), x_cpy (n), y_cpy (n);
      int incx, incy;

      incx = 1;
      incy = 1;
      for (int i = 0; i < n; ++i) {
        x[i] = static_cast<ST> (i + 1);
        x_cpy[i] = static_cast<ST> (i + 1);
        y[i] = 2 * static_cast<ST> (i + 1);
        y_cpy[i] = 2 * static_cast<ST> (i + 1);
      }
      blas.SWAP (n, x.getRawPtr (), incx, y.getRawPtr (), incy);
      TEST_COMPARE_ARRAYS( x, y_cpy );
      TEST_COMPARE_ARRAYS( y, x_cpy );

      // FIXME (mfh 16 Sep 2015) Fix the negative and strided INCX and
      // INCY cases.
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExpBlockView, LAPY2, ST )
  {
    if (! Teuchos::ScalarTraits<ST>::isOrdinal) { // skip integer types
      typename Tpetra::Details::GetLapackType<ST>::lapack_type lapack;
      // Rough tolerance for rounding errors.  LAPY2 uses a different
      // formula, so I expect it to commit different rounding error.
      const auto tol = 10.0 * Teuchos::ScalarTraits<ST>::eps ();
      ST x, y;

      x = 2.0;
      y = 3.0;
      auto correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      auto lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = -2.0;
      y = 3.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = 2.0;
      y = -3.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = 0.0;
      y = 3.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = 5.0;
      y = 0.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = 0.0;
      y = -3.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );

      x = -5.0;
      y = 0.0;
      correctResult = Teuchos::ScalarTraits<ST>::squareroot (x*x + y*y);
      lapy2Result = lapack.LAPY2 (x, y);
      TEST_FLOATING_EQUALITY( correctResult, lapy2Result, tol );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExpBlockView, LARFGP, ST )
  {
    if (! Teuchos::ScalarTraits<ST>::isOrdinal) { // skip integer types
      typename Tpetra::Details::GetLapackType<ST>::lapack_type lapack;

      const ST zero = Teuchos::ScalarTraits<ST>::zero ();
      const ST one = Teuchos::ScalarTraits<ST>::one ();
      // Rough tolerance for rounding errors.
      const auto tol = 10.0 * Teuchos::ScalarTraits<ST>::eps ();
      const int n = 2;
      const int incx = 1;
      ST alpha, x[2], tau;

      alpha = 0.0;
      x[0] = 0.0;
      x[1] = 0.0;
      tau = 0.0;

      lapack.LARFG (n, &alpha, x, incx, &tau);
      TEST_FLOATING_EQUALITY( tau, zero, tol );

      alpha = 0.0;
      x[0] = 1.0;
      x[1] = 0.0;
      tau = 1.0;

      lapack.LARFG (n, &alpha, x, incx, &tau);
      TEST_FLOATING_EQUALITY( tau, one, tol );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, SCAL, ST, LO )
  {
    typedef typename Kokkos::Details::ArithTraits<ST>::val_type IST; // "impl_scalar_type"

    typedef Tpetra::Experimental::LittleBlock<IST, LO> blk_type;
    typedef Tpetra::Experimental::LittleVector<IST, LO> vec_type;
    const IST zero = static_cast<IST> (0.0);
    const IST one = static_cast<IST> (1.0);
    const IST two = one + one;
    const IST three = two + one;
    const IST five = three + two;
    const IST six = three + three;
    const IST ten = five + five;
    const LO minBlkSize = 1; // 1x1 "blks" should also work
    const LO maxBlkSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blkPool (2 * maxBlkSize * maxBlkSize);
    // Memory pool for the LittleVector instances.
    Teuchos::Array<IST> vecPool (2 * maxBlkSize);

    for (LO blkSize = minBlkSize; blkSize <= maxBlkSize; ++blkSize) {
      blk_type A1 (blkPool (0, blkSize*blkSize).getRawPtr (),
                   blkSize, 1, blkSize);
      blk_type A2 (blkPool (blkSize*blkSize, blkSize*blkSize).getRawPtr (),
                   blkSize, 1, blkSize);
      vec_type x1 (vecPool (0, blkSize).getRawPtr (),
                   blkSize, 1);
      vec_type x2 (vecPool (blkSize, blkSize).getRawPtr (),
                   blkSize, 1);

      // A1 == A2 and x1 == x2.  We will use SCAL on A1 and x1, and
      // use conventional loops on A2 and x2, then compare the
      // results.  The numbers are small enough that the test need not
      // worry about rounding error.  We use a different value for
      // each entry, in order to catch possible layout bugs (e.g.,
      // mixing up the strides).
      IST curVecVal = one;
      IST curBlkVal = one;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          A1(i,j) = curBlkVal;
          A2(i,j) = curBlkVal;
          curBlkVal += one;
        }
        x1(j) = curVecVal;
        x2(j) = curVecVal;
        curVecVal += one;
      }

      Tpetra::Experimental::SCAL (two, A1);
      Tpetra::Experimental::SCAL (three, x1);

      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          A2(i,j) *= two;
        }
        x2(j) *= three;
      }

      bool blksEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          if (A1(i,j) != A2(i,j)) {
            blksEq = false;
            break;
          }
        }
      }
      TEST_ASSERT( blksEq );

      bool vecsEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        if (x1(j) != x2(j)) {
          vecsEq = false;
          break;
        }
      }
      TEST_ASSERT( vecsEq );
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, COPY, ST, LO )
  {
    typedef typename Kokkos::Details::ArithTraits<ST>::val_type IST; // "impl_scalar_type"

    typedef Tpetra::Experimental::LittleBlock<IST, LO> blk_type;
    typedef Tpetra::Experimental::LittleVector<IST, LO> vec_type;
    const IST zero = static_cast<IST> (0.0);
    const IST one = static_cast<IST> (1.0);
    const IST two = one + one;
    const IST three = two + one;
    const IST five = three + two;
    const IST six = three + three;
    const IST ten = five + five;
    const LO minBlkSize = 1; // 1x1 "blks" should also work
    const LO maxBlkSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blkPool (3 * maxBlkSize * maxBlkSize);
    // Memory pool for the LittleVector instances.
    Teuchos::Array<IST> vecPool (3 * maxBlkSize);

    for (LO blkSize = minBlkSize; blkSize <= maxBlkSize; ++blkSize) {
      blk_type A1 (blkPool (0, blkSize*blkSize).getRawPtr (),
                   blkSize, 1, blkSize);
      blk_type A2 (blkPool (blkSize*blkSize, blkSize*blkSize).getRawPtr (),
                   blkSize, 1, blkSize);
      blk_type A3 (blkPool (2*blkSize*blkSize, blkSize*blkSize).getRawPtr (),
                   blkSize, 1, blkSize);
      vec_type x1 (vecPool (0, blkSize).getRawPtr (),
                   blkSize, 1);
      vec_type x2 (vecPool (blkSize, blkSize).getRawPtr (),
                   blkSize, 1);
      vec_type x3 (vecPool (2*blkSize, blkSize).getRawPtr (),
                   blkSize, 1);

      // A1 == A2 and x1 == x2.  We will use COPY to copy A1 into A3
      // and x1 into A3, then compare the result against A2 resp. x2.
      IST curVecVal = one;
      IST curBlkVal = one;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          A1(i,j) = curBlkVal;
          A2(i,j) = curBlkVal;
          A3(i,j) = zero;
          curBlkVal += one;
        }
        x1(j) = curVecVal;
        x2(j) = curVecVal;
        x3(j) = zero;
        curVecVal += one;
      }

      Tpetra::Experimental::COPY (A1, A3);
      Tpetra::Experimental::COPY (x1, x3);

      bool blksEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          if (A2(i,j) != A3(i,j)) {
            blksEq = false;
            break;
          }
        }
      }
      TEST_ASSERT( blksEq );

      bool vecsEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        if (x2(j) != x3(j)) {
          vecsEq = false;
          break;
        }
      }
      TEST_ASSERT( vecsEq );
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, SolveIdentity, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, SCAL, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, COPY, SCALAR, LOCAL_ORDINAL )

#define UNIT_TEST_GROUP2( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExpBlockView, SWAP, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExpBlockView, LAPY2, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExpBlockView, LARFGP, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  //
  // FIXME (mfh 17 Sep 2015) Fails for __float128!
  //
  TPETRA_INSTANTIATE_SL_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

  // FIXME (mfh 17 Sep 2015) Define ETI / test macros for real Scalar
  // types only.  Note that in LAPACK, _LAPY2 only exists for _ = S, D
  // -- thus, real-valued Scalar types.

  //TPETRA_INSTANTIATE_S_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP2 )

#ifdef TPETRA_INST_FLOAT
  UNIT_TEST_GROUP2( float )
#endif // TPETRA_INST_FLOAT

#ifdef TPETRA_INST_DOUBLE
  UNIT_TEST_GROUP2( double )
#endif // TPETRA_INST_DOUBLE

#ifdef TPETRA_INST_FLOAT128
  UNIT_TEST_GROUP2( __float128 )
#endif // TPETRA_INST_FLOAT128

} // namespace (anonymous)


