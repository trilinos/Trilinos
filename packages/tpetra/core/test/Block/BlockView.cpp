// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_BlockView.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef HAVE_TPETRA_INST_FLOAT128
#  include "Teuchos_Details_Lapack128.hpp"
#endif // HAVE_TPETRA_INST_FLOAT128
#ifdef HAVE_TPETRA_INST_LONG_DOUBLE
#  include "Teuchos_Details_LapackLongDouble.hpp"
#endif // HAVE_TPETRA_INST_LONG_DOUBLE

#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using std::endl;
  typedef Teuchos::Array<int>::size_type size_type;

  /// \brief Return the Teuchos::LAPACK specialization corresponding
  ///   to the given Scalar type.
  ///
  /// The reason this exists is the same reason why the
  /// impl_scalar_type typedef in Tpetra::MultiVector may differ from
  /// its Scalar template parameter.  For example, Scalar =
  /// std::complex<T> corresponds to impl_scalar_type =
  /// Kokkos::complex<T>.  The latter has no Teuchos::LAPACK
  /// specialization, so we have to map it back to std::complex<T>.
  template<class Scalar>
  struct GetLapackType {
    typedef Scalar lapack_scalar_type;
    typedef Teuchos::LAPACK<int, Scalar> lapack_type;
  };

  template<class T>
  struct GetLapackType<Kokkos::complex<T> > {
    typedef std::complex<T> lapack_scalar_type;
    typedef Teuchos::LAPACK<int, std::complex<T> > lapack_type;
  };

#ifdef HAVE_TPETRA_INST_FLOAT128
  template<>
  struct GetLapackType<__float128> {
    typedef __float128 lapack_scalar_type;
    // Use the Lapack128 class we declared above to implement the
    // linear algebra operations needed for small dense blocks and
    // vectors.
    typedef Teuchos::Details::Lapack128 lapack_type;
  };
#endif // HAVE_TPETRA_INST_FLOAT128

  //
  // UNIT TESTS
  //

  // The "little" blocks and vectors do not depend on Tpetra's
  // GlobalOrdinal type.  This is why we only include three template
  // parameters: Scalar (ST) and LocalOrdinal (LO).  For now, we just
  // test with the default execution and memory spaces, but at some
  // point, it would make sense to test over all the enabled spaces.

  // This example tests the factorization routine with a matrix
  // that does not require partial pivoting
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, Factor, ST, LO )
  {
    using Teuchos::Array;
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;
    typedef Teuchos::LAPACK<LO, ST> lapack_type;
    typedef Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace> block_type;
    typedef Kokkos::View<LO*, Kokkos::HostSpace> int_vec_type;
    typedef Kokkos::View<IST*, Kokkos::HostSpace> scalar_vec_type;

    const auto tol = 10.0 * Kokkos::ArithTraits<IST>::eps ();

    TEST_ASSERT( Kokkos::is_initialized () );
    if (! Kokkos::is_initialized ()) {
      return; // don't bother to continue
    }

    // Create a matrix
    block_type A("A",3,3);
    A(0,0) = 10;
    A(0,1) = 2;
    A(0,2) = 3;
    A(1,0) = 4;
    A(1,1) = 20;
    A(1,2) = 5;
    A(2,0) = 6;
    A(2,1) = 7;
    A(2,2) = 30;

    // Create the pivot vector
    int_vec_type ipiv("ipiv",3);

    // Compute the true factorization
    block_type true_A("trueA",3,3);
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        true_A(i,j) = A(i,j);
    int_vec_type true_piv("trueipiv",3);
    LO info;

    lapack_type lapackOBJ;
    lapackOBJ.GETRF (3, 3, reinterpret_cast<ST*> (true_A.data ()), 3, true_piv.data(), &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute our factorization
    Tpetra::GETF2<block_type,int_vec_type>(A,ipiv,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the two solutions
    Teuchos::ArrayView<ST> ptr1 (reinterpret_cast<ST*> (A.data ()), 9);
    Teuchos::ArrayView<ST> ptr2 (reinterpret_cast<ST*> (true_A.data ()), 9);
    TEST_COMPARE_FLOATING_ARRAYS( ptr1, ptr2, tol );
    TEST_COMPARE_ARRAYS( ipiv, true_piv );

    // Create a RHS
    scalar_vec_type rhs("rhs",3);
    rhs(0) = 100.0;
    rhs(1) = 200.0;
    rhs(2) = 300.0;

    scalar_vec_type true_rhs("truerhs",3);
    for(int i=0; i<3; i++)
      true_rhs(i) = rhs(i);

    // Compute the true solution
    lapackOBJ.GETRS ('n', 3, 1, reinterpret_cast<ST*> (true_A.data ()),
                     3, true_piv.data (),
                     reinterpret_cast<ST*> (true_rhs.data ()), 3,
                     &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute our solution
    Tpetra::GETRS<block_type,int_vec_type,scalar_vec_type>("n",A,ipiv,rhs,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the solutions
    Teuchos::ArrayView<ST> ptr3 (reinterpret_cast<ST*> (rhs.data ()), 3);
    Teuchos::ArrayView<ST> ptr4 (reinterpret_cast<ST*> (true_rhs.data ()), 3);
    TEST_COMPARE_FLOATING_ARRAYS( ptr3, ptr4, tol );

    // Compute the inverse
    scalar_vec_type work("work",3);
    Tpetra::GETRI<block_type,int_vec_type,scalar_vec_type>(A,ipiv,work,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute the true inverse
    lapackOBJ.GETRI (3, reinterpret_cast<ST*> (true_A.data ()), 3,
                     true_piv.data (),
                     reinterpret_cast<ST*> (work.data ()), 3, &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the inverses
    TEST_COMPARE_FLOATING_ARRAYS( ptr1, ptr2, tol );
  }


  // This example tests the factorization routine with a matrix
  // that requires partial pivoting
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, FactorPivot, ST, LO )
  {
    using Teuchos::Array;
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;
    typedef Teuchos::LAPACK<LO, ST> lapack_type;
    typedef Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace> block_type;
    typedef Kokkos::View<LO*, Kokkos::HostSpace> int_vec_type;
    typedef Kokkos::View<IST*, Kokkos::HostSpace> scalar_vec_type;

    const auto tol = 10.0 * Kokkos::ArithTraits<IST>::eps ();

    TEST_ASSERT( Kokkos::is_initialized () );
    if (! Kokkos::is_initialized ()) {
      return; // don't bother to continue
    }

    // Create a matrix
    block_type A("A",3,3);
    A(2,0) = 10.0;
    A(2,1) = 2.0;
    A(2,2) = 3.0;
    A(0,0) = 4.0;
    A(0,1) = 20.0;
    A(0,2) = 5.0;
    A(1,0) = 6.0;
    A(1,1) = 7.0;
    A(1,2) = 30.0;

    // Create the pivot vector
    int_vec_type ipiv("ipiv",3);

    // Compute the true factorization
    block_type true_A("trueA",3,3);
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        true_A(i,j) = A(i,j);
    int_vec_type true_piv("trueipiv",3);
    LO info;

    lapack_type lapackOBJ;
    lapackOBJ.GETRF (3, 3, reinterpret_cast<ST*> (true_A.data ()),
                     3, true_piv.data (), &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute our factorization
    Tpetra::GETF2<block_type,int_vec_type>(A,ipiv,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the two solutions
    Teuchos::ArrayView<ST> ptr1 (reinterpret_cast<ST*> (A.data ()), 9);
    Teuchos::ArrayView<ST> ptr2 (reinterpret_cast<ST*> (true_A.data ()), 9);
    TEST_COMPARE_FLOATING_ARRAYS( ptr1, ptr2, tol );
    TEST_COMPARE_ARRAYS( ipiv, true_piv );

    // Create a RHS
    scalar_vec_type rhs("rhs",3);
    rhs(0) = 100.0;
    rhs(1) = 200.0;
    rhs(2) = 300.0;

    scalar_vec_type true_rhs("truerhs",3);
    for(int i=0; i<3; i++)
      true_rhs(i) = rhs(i);

    // Compute the true solution
    lapackOBJ.GETRS ('n', 3, 1, reinterpret_cast<ST*> (true_A.data ()),
                     3, true_piv.data (),
                     reinterpret_cast<ST*> (true_rhs.data ()), 3,
                     &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute our solution
    Tpetra::GETRS<block_type,int_vec_type,scalar_vec_type>("n",A,ipiv,rhs,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the solutions
    Teuchos::ArrayView<ST> ptr3 (reinterpret_cast<ST*> (rhs.data ()), 3);
    Teuchos::ArrayView<ST> ptr4 (reinterpret_cast<ST*> (true_rhs.data ()), 3);
    TEST_COMPARE_FLOATING_ARRAYS( ptr3, ptr4, tol );

    // Compute the inverse
    scalar_vec_type work("work",3);
    Tpetra::GETRI<block_type,int_vec_type,scalar_vec_type>(A,ipiv,work,info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compute the true inverse
    lapackOBJ.GETRI (3, reinterpret_cast<ST*> (true_A.data ()), 3,
                     true_piv.data (),
                     reinterpret_cast<ST*> (work.data ()), 3, &info);
    TEST_EQUALITY_CONST( info, 0 );

    // Compare the inverses
    TEST_COMPARE_FLOATING_ARRAYS( ptr1, ptr2, tol );
  }


  // Test small dense block LU factorization and solve, with an easy
  // problem (the identity matrix).
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, SolveIdentity, ST, LO )
  {
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;
    typedef Kokkos::View<IST**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> block_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vec_type;
    typedef Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> piv_type;
    const IST zero = static_cast<IST> (0.0);
    const IST one = static_cast<IST> (1.0);
    const LO minBlockSize = 1; // 1x1 "blocks" should also work
    const LO maxBlockSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blockPool (maxBlockSize * maxBlockSize);
    // Memory pool for the LittleVector instances (x and b).
    Teuchos::Array<IST> vecPool (maxBlockSize * 2);
    // Memory pool for the pivot vector.
    Teuchos::Array<int> ipivPool (maxBlockSize);

    for (LO blockSize = minBlockSize; blockSize <= maxBlockSize; ++blockSize) {
      block_type A (blockPool (0, blockSize*blockSize).getRawPtr (), blockSize, blockSize);
      Teuchos::ArrayView<IST> x_view = vecPool (0, blockSize);
      vec_type x (x_view.getRawPtr (), blockSize);
      Teuchos::ArrayView<IST> b_view = vecPool (blockSize, blockSize);
      vec_type b (b_view.getRawPtr (), blockSize);
      piv_type ipiv (ipivPool.getRawPtr (), blockSize);

      Tpetra::FILL (A, zero); // assign zero to each entry
      for (LO i = 0; i < blockSize; ++i) {
        A(i,i) = one;
        b(i) = static_cast<IST> (i + 1);
        x(i) = b(i); // copy of right-hand side on input
        ipiv(i) = 0;
      }

      int info = 0;
      std::cerr << "Factor A for blockSize = " << blockSize << std::endl;
      Tpetra::GETF2 (A, ipiv, info);
      TEST_EQUALITY_CONST( info, 0 );
      if (info == 0) {
        std::cerr << "Solve: blockSize = " << blockSize << std::endl;
        Tpetra::GETRS ("N", A, ipiv, x, info);
      }
      std::cerr << "Done with factor and solve" << std::endl;

      // Re-fill b, in case A.solve brokenly clobbered it.
      for (LO i = 0; i < blockSize; ++i) {
        b(i) = static_cast<IST> (i + 1);
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
    typedef typename Kokkos::ArithTraits<ST>::val_type IST;
    typedef Kokkos::View<IST**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> block_type;
    const IST zero = static_cast<IST> (0.0);
    const IST one = static_cast<IST> (1.0);
    const LO minBlockSize = 1; // 1x1 "blocks" should also work
    const LO maxBlockSize = 32;

    typename GetLapackType<IST>::lapack_type lapack; // ??? does IST (Kokkos::complex) work here?

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blockPool (maxBlockSize * maxBlockSize);
    // Memory pool for LAPACK's temporary workspace.  It might need to
    // be resized in the loop below, because LAPACK may want more
    // workspace than just the minimum (the number of columns, which
    // is what the BLAS 2 QR factorization GEQR2 requires).  This must
    // have length at least one, for the workspace query.
    Teuchos::Array<IST> workPool (std::max (1, maxBlockSize));
    // Memory pool for the TAU output array.
    Teuchos::Array<IST> tauPool (maxBlockSize);

    for (LO blockSize = minBlockSize; blockSize <= maxBlockSize; ++blockSize) {
      block_type A (blockPool (0, blockSize*blockSize).getRawPtr (), blockSize, blockSize);

      // Fill A with the identity matrix.
      Tpetra::FILL (A, zero); // assign zero to each entry
      for (LO i = 0; i < blockSize; ++i) {
        A(i,i) = one;
      }

      Teuchos::ArrayView<IST> tauView = tauPool (0, blockSize);

      // Workspace query.
      Teuchos::ArrayView<IST> workView = workPool (0, std::max (1, blockSize));
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
      lwork = static_cast<int> (Kokkos::ArithTraits<IST>::real (workView[0]));
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
      typename GetLapackType<ST>::lapack_type lapack;
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
      typename GetLapackType<ST>::lapack_type lapack;

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
    typedef typename Kokkos::ArithTraits<ST>::val_type IST; // "impl_scalar_type"
    typedef Kokkos::View<IST**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> blk_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vec_type;

    //const IST zero = static_cast<IST> (0.0); // unused
    const IST one = static_cast<IST> (1.0);
    const IST two = one + one;
    const IST three = two + one;
    //const IST five = three + two; // unused
    //const IST six = three + three; // unused
    //const IST ten = five + five; // unused
    const LO minBlkSize = 1; // 1x1 "blks" should also work
    const LO maxBlkSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blkPool (2 * maxBlkSize * maxBlkSize);
    // Memory pool for the LittleVector instances.
    Teuchos::Array<IST> vecPool (2 * maxBlkSize);

    for (LO blkSize = minBlkSize; blkSize <= maxBlkSize; ++blkSize) {
      blk_type A1 (blkPool (0, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A2 (blkPool (blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      vec_type x1 (vecPool (0, blkSize).getRawPtr (), blkSize);
      vec_type x2 (vecPool (blkSize, blkSize).getRawPtr (), blkSize);

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

      Tpetra::SCAL (two, A1);
      Tpetra::SCAL (three, x1);

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
    typedef typename Kokkos::ArithTraits<ST>::val_type IST; // "impl_scalar_type"
    typedef Kokkos::View<IST**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> blk_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vec_type;

    const IST zero = static_cast<IST> (0.0);
    const IST one = static_cast<IST> (1.0);
    //const IST two = one + one; // unused
    //const IST three = two + one; // unused
    //const IST five = three + two; // unused
    //const IST six = three + three; // unused
    //const IST ten = five + five; // unused
    const LO minBlkSize = 1; // 1x1 "blks" should also work
    const LO maxBlkSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blkPool (3 * maxBlkSize * maxBlkSize);
    // Memory pool for the LittleVector instances.
    Teuchos::Array<IST> vecPool (3 * maxBlkSize);

    for (LO blkSize = minBlkSize; blkSize <= maxBlkSize; ++blkSize) {
      blk_type A1 (blkPool (0, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A2 (blkPool (blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A3 (blkPool (2*blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      vec_type x1 (vecPool (0, blkSize).getRawPtr (), blkSize);
      vec_type x2 (vecPool (blkSize, blkSize).getRawPtr (), blkSize);
      vec_type x3 (vecPool (2*blkSize, blkSize).getRawPtr (), blkSize);

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

      Tpetra::COPY (A1, A3);
      Tpetra::COPY (x1, x3);

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


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ExpBlockView, AXPY, ST, LO )
  {
    typedef typename Kokkos::ArithTraits<ST>::val_type IST; // "impl_scalar_type"
    typedef Kokkos::View<IST**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> blk_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vec_type;

    //const IST zero = static_cast<IST> (0.0); // unused
    const IST one = static_cast<IST> (1.0);
    const IST two = one + one;
    const IST three = two + one;
    const IST five = three + two;
    const LO minBlkSize = 1; // 1x1 "blks" should also work
    const LO maxBlkSize = 32;

    // Memory pool for the LittleBlock instances.
    Teuchos::Array<IST> blkPool (4 * maxBlkSize * maxBlkSize);
    // Memory pool for the LittleVector instances.
    Teuchos::Array<IST> vecPool (4 * maxBlkSize);

    for (LO blkSize = minBlkSize; blkSize <= maxBlkSize; ++blkSize) {
      blk_type A1 (blkPool (0, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A2 (blkPool (blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A3 (blkPool (2*blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      blk_type A4 (blkPool (3*blkSize*blkSize, blkSize*blkSize).getRawPtr (), blkSize, blkSize);
      vec_type x1 (vecPool (0, blkSize).getRawPtr (), blkSize);
      vec_type x2 (vecPool (blkSize, blkSize).getRawPtr (), blkSize);
      vec_type x3 (vecPool (2*blkSize, blkSize).getRawPtr (), blkSize);
      vec_type x4 (vecPool (3*blkSize, blkSize).getRawPtr (), blkSize);

      // Compare AXPY(alpha, A1, A3) and AXPY(alpha, x1, x3) with the
      // manual equivalent of AXPY(alpha, A2, A4) resp. AXPY(alpha,
      // x2, x4).
      IST curVecVal = one;
      IST curBlkVal = one;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          A1(i,j) = curBlkVal;
          A2(i,j) = curBlkVal;
          A3(i,j) = three; // just something different
          A4(i,j) = three;
          curBlkVal += one;
        }
        x1(j) = curVecVal;
        x2(j) = curVecVal;
        x3(j) = three; // just something different
        x4(j) = three; // just something different
        curVecVal += one;
      }

      const IST alpha = five;

      Tpetra::AXPY (alpha, A1, A3); // A3 := A3 + alpha*A1
      Tpetra::AXPY (alpha, x1, x3); // x3 := x3 + alpha*x1

      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          A4(i,j) = A4(i,j) + alpha * A2(i,j);
        }
        x4(j) = x4(j) + alpha * x2(j);
      }

      bool blksEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        for (LO i = 0; i < blkSize; ++i) {
          if (A3(i,j) != A4(i,j)) {
            blksEq = false;
            break;
          }
        }
      }
      TEST_ASSERT( blksEq );

      bool vecsEq = true;
      for (LO j = 0; j < blkSize; ++j) {
        if (x3(j) != x4(j)) {
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, Factor, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, FactorPivot, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, SolveIdentity, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, SCAL, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, COPY, SCALAR, LOCAL_ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ExpBlockView, AXPY, SCALAR, LOCAL_ORDINAL )

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

#ifdef TPETRA_INST_LONG_DOUBLE
  UNIT_TEST_GROUP2( long double )
#endif // TPETRA_INST_LONG_DOUBLE

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);
  Kokkos::initialize (argc, argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  Kokkos::finalize ();
  return errCode;
}
