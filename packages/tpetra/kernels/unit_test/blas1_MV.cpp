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

#include <Kokkos_Blas1_MV.hpp>
#include <Teuchos_Comm.hpp>
#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#else
#  include <Teuchos_DefaultSerialComm.hpp>
#endif // HAVE_MPI

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;

namespace { // (anonymous)

enum EWhichNorm {
  TPETRAKERNELS_TEST_NORM_INF = 0,
  TPETRAKERNELS_TEST_NORM_ONE = 1,
  TPETRAKERNELS_TEST_NORM_TWO = 2,
  TPETRAKERNELS_TEST_NORM_INVALID = 3
};

template<class RV, class XMV>
void
kokkosNorm (const RV& norms, const XMV& X, const EWhichNorm whichNorm)
{
  if (whichNorm == TPETRAKERNELS_TEST_NORM_INF) {
    KokkosBlas::nrmInf (norms, X);
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_ONE) {
    KokkosBlas::nrm1 (norms, X);
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_TWO) {
    KokkosBlas::nrm2_squared (norms, X);
  }
}

template<class RV, class XMV>
void
kokkosNorm (const RV& norms, const XMV& X, const size_t j,
            const EWhichNorm whichNorm)
{
  if (whichNorm == TPETRAKERNELS_TEST_NORM_INF) {
    KokkosBlas::nrmInf (norms, X, j);
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_ONE) {
    KokkosBlas::nrm1 (norms, X, j);
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_TWO) {
    KokkosBlas::nrm2_squared (norms, X, j);
  }
}

} // namespace (anonymous)

template<class Scalar, class Layout, class Device>
bool
testFill (std::ostream& out, const bool prvSuccess)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  bool curSuccess = true;

  out << "Testing KokkosBlas::fill" << endl;

  const size_type numRows = 10;
  const size_type numCols = 3;

  mv_type X ("X", numRows, numCols);
  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);

  out << "  Test fill with zero" << endl;
  KokkosBlas::fill (X, ATS::zero ());
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != ATS::zero ()) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j) << " != " << ATS::zero () << endl;
      }
    }
  }

  out << "  Test fill with one" << endl;
  KokkosBlas::fill (X, ATS::one ());
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != ATS::one ()) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j) << " != " << ATS::one () << endl;
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess && prvSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testSum (std::ostream& out, const bool prvSuccess)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  // sum() uses the Device's preferred Layout for the output array.
  typedef Kokkos::View<Scalar*, Device> sums_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  bool curSuccess = true;

  out << "Testing KokkosBlas::sum" << endl;

  const size_type numRows = 4;
  const size_type numCols = 3;

  mv_type X ("X", numRows, numCols);
  sums_type r ("r", numCols);

  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);
  typename sums_type::HostMirror r_h = Kokkos::create_mirror_view (r);

  out << "  Test that the sum of zeros is zero" << endl;
  KokkosBlas::fill (X, ATS::zero ());
  KokkosBlas::sum (r, X);
  Kokkos::deep_copy (r_h, r);
  for (size_type j = 0; j < numCols; ++j) {
    if (r_h(j) != ATS::zero ()) {
      curSuccess = false;
      out << "    FAILED: r_h(" << j << ") = " << r_h(j)
          << " != " << ATS::zero () << endl;
    }
  }

  // Tetractys test.
  out << "  Test that the sum of [1, 2, 3, 4] is 10" << endl;
  const Scalar ONE = ATS::one ();
  const Scalar TWO = ONE + ONE;
  const Scalar THREE = TWO + ONE;
  const Scalar FOUR = THREE + ONE;
  const Scalar TEN = ONE + TWO + THREE + FOUR;
  for (size_type j = 0; j < numCols; ++j) {
    X_h(0,j) = ONE;
    X_h(1,j) = TWO;
    X_h(2,j) = THREE;
    X_h(3,j) = FOUR;
  }
  Kokkos::deep_copy (X, X_h);
  KokkosBlas::sum (r, X);
  KokkosBlas::fill (X, ATS::zero ());
  KokkosBlas::sum (r, X);
  Kokkos::deep_copy (r_h, r);
  for (size_type j = 0; j < numCols; ++j) {
    if (r_h(j) != TEN) {
      curSuccess = false;
      out << "    FAILED: r_h(" << j << ") = " << r_h(j)
          << " != " << TEN << endl;
    }
  }

#ifdef KOKKOS_HAVE_CXX11
  if (numCols > 1) {
    out << "  Repeat previous test, one column at a time" << endl;
    // Make sure that we get the same result one column at a time, for a
    // single vector (1-D Views), as we get when processing all the
    // columns of the multivector at once (2-D Views).
    for (size_type j = 0; j < numCols; ++j) {
      auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
      auto r_j = Kokkos::subview (r, j);
      KokkosBlas::sum (r_j, X_j);
    }
    Kokkos::deep_copy (r_h, r);
    for (size_type j = 0; j < numCols; ++j) {
      if (r_h(j) != TEN) {
        curSuccess = false;
        out << "    FAILED: r_h(" << j << ") = " << r_h(j)
            << " != " << TEN << endl;
      }
    }
  }
#endif // KOKKOS_HAVE_CXX11

  // Make sure that sum() and nrm1() are different, by changing the
  // tetractys test slightly.
  out << "  Test that the sum of [-1, 2, -3, 4] is 2" << endl;
  for (size_type j = 0; j < numCols; ++j) {
    X_h(0,j) = -ONE;
    X_h(1,j) = TWO;
    X_h(2,j) = -THREE;
    X_h(3,j) = FOUR;
  }
  Kokkos::deep_copy (X, X_h);
  KokkosBlas::sum (r, X);
  Kokkos::deep_copy (r_h, r);
  for (size_type j = 0; j < numCols; ++j) {
    if (r_h(j) != TWO) {
      curSuccess = false;
      out << "    FAILED: r_h(" << j << ") = " << r_h(j)
          << " != " << TWO << endl;
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess && prvSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testAnyNorm (std::ostream& out, const EWhichNorm whichNorm, const bool prvSuccess)
{
  using Kokkos::subview;
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> MV;
  typedef typename MV::size_type size_type;
  // typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  typedef Kokkos::Details::InnerProductSpaceTraits<Scalar> IPT;
  typedef typename IPT::mag_type mag_type;
  typedef Kokkos::View<mag_type*, Device> norms_type;
  typedef std::pair<size_t, size_t> pair_type;
  typedef Kokkos::Details::ArithTraits<mag_type> ATM;
  bool curSuccess = true;

  const size_type numRows = 10;
  const size_type numCols = 3;
  norms_type X_norms ("X_norms", numCols);
  typename norms_type::HostMirror X_norms_h = Kokkos::create_mirror_view (X_norms);

  out << "Testing properties common to all norms, for norm ";
  if (whichNorm == TPETRAKERNELS_TEST_NORM_INF) {
    out << "inf";
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_ONE) {
    out << "1";
  } else if (whichNorm == TPETRAKERNELS_TEST_NORM_TWO) {
    out << "2";
  }
  out << endl;

  // Norm of an MV with zero rows is zero.
  out << "  Test that the norm(s) of a(n M)V with zero rows is zero" << endl;
  {
    MV X_empty ("X_empty", 0, 1);
    kokkosNorm<norms_type, MV> (subview (X_norms, pair_type (0, 1)), X_empty, whichNorm);

    Kokkos::deep_copy (X_norms_h, X_norms);
    if (X_norms_h(0) != ATM::zero ()) {
      curSuccess = false;
      out << "    FAILED" << endl;
    }
  }

  // Norm of an MV with zero rows is zero, no matter how many columns it has.
  out << "  Test that the norm(s) of an 0 x " << numCols << " MV is/are zero" << endl;
  {
    MV X_empty ("X_empty", 0, numCols);
    kokkosNorm<norms_type, MV> (subview (X_norms, pair_type (0, numCols)), X_empty, whichNorm);

    Kokkos::deep_copy (X_norms_h, X_norms);
    if (X_norms_h(0) != ATM::zero ()) {
      out << "    FAILED" << endl;
    }
  }

  // Degenerate MVs (with zero rows, columns, or both) must be valid input.
  // They shouldn't make KokkosBlas throw or segfault.
  out << "  Test that degenerate MVs don't cause a throw or segfault" << endl;
  {
    MV X0 ("X0", 0, 0);
    MV X1 ("X1", numRows, 0);
    // This causes a throw below, because it creates a 1 x 1 View for some reason.
    //
    //norms_type X_norms_degenerate = subview (X_norms, pair_type (0, 0));
    norms_type X_norms_degenerate ("X_norms_degenerate", 0);

    try {
      kokkosNorm<norms_type, MV> (X_norms_degenerate, X0, whichNorm);
      kokkosNorm<norms_type, MV> (X_norms_degenerate, X1, whichNorm);
    } catch (std::exception& e) {
      curSuccess = false;
      out << "    FAILED: Norm threw exception: " << e.what () << endl;
    }
  }

  // Norm of a nondegenerate MV full of zeros is always zero.
  out << "  Test that the norms of a nondegenerate MV are always zero" << endl;
  {
    MV X1 ("X1", numRows, 1);
    MV X2 ("X2", numRows, numCols);

    kokkosNorm<norms_type, MV> (subview (X_norms, pair_type (0, 1)), X1, whichNorm);
    Kokkos::deep_copy (X_norms_h, X_norms);
    if (X_norms_h(0) != ATM::zero ()) {
      curSuccess = false;
      out << "    FAILED: Norm of a nondegenerate single-column MV full of zeros was not zero" << endl;
    }

    kokkosNorm<norms_type, MV> (subview (X_norms, pair_type (0, numCols)), X2, whichNorm);
    Kokkos::deep_copy (X_norms_h, X_norms);
    for (size_type j = 0; j < numCols; ++j) {
      if (X_norms_h(j) != ATM::zero ()) {
        curSuccess = false;
        out << "    FAILED: Norm of a nondegenerate multicolumn MV full of zeros was not zero" << endl;
        break;
      }
    }
  }
  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess && prvSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testNormInf (std::ostream& out, const bool prvSuccess)
{
  using Kokkos::subview;
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> MV;
  typedef typename MV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  typedef Kokkos::Details::InnerProductSpaceTraits<Scalar> IPT;
  typedef typename IPT::mag_type mag_type;
  typedef Kokkos::View<mag_type*, Device> norms_type;
  typedef Kokkos::Details::ArithTraits<mag_type> ATM;
  bool curSuccess = true;

  const size_type numRows = 10;
  const size_type numCols = 3;
  norms_type X_norms ("X_norms", numCols);
  typename norms_type::HostMirror X_norms_h =
    Kokkos::create_mirror_view (X_norms);

  out << "Testing inf-norm" << endl;

  MV X ("X", numRows, numCols);
  typename MV::HostMirror X_host = Kokkos::create_mirror_view (X);
  norms_type norms ("norms", numCols);
  typename norms_type::HostMirror norms_host = Kokkos::create_mirror_view (norms);

  // Test that inf-norm picks out the max value in each column.
  out << "  Test that inf-norm picks out the max value" << endl;
  {
    KokkosBlas::fill (X_host, ATS::zero ());
    for (size_type j = 0; j < numCols; ++j) {
      const size_type i = (j+1) % numRows;
      X_host(i,j) = (j % 2 == size_type (0)) ? ATS::one () : -ATS::one ();
    }
    Kokkos::deep_copy (X, X_host);
    KokkosBlas::nrmInf (norms, X);
    Kokkos::deep_copy (norms_host, norms);

    for (size_type j = 0; j < numCols; ++j) {
      if (norms_host(j) != ATM::one ()) {
        curSuccess = false;
        out << "    FAILED: norms_host(j) = " << norms_host(j)
            << " != " << ATM::one () << endl;
      }
    }
  }

  // Test that inf-norm [-2,1] is 2.
  out << "  Test that inf-norm [-2,1] is 2" << endl;
  {
    MV X1 ("X1", 2, 1);
    norms_type norms1 ("norms1", 1);
    typename MV::HostMirror X1_host = Kokkos::create_mirror_view (X1);
    X1_host(0,0) = -(ATS::one () + ATS::one ());
    X1_host(1,0) = +ATS::one ();
    Kokkos::deep_copy (X1, X1_host);

    KokkosBlas::nrmInf (norms1, X1);
    typename norms_type::HostMirror norms1_host =
      Kokkos::create_mirror_view (norms1);
    Kokkos::deep_copy (norms1_host, norms1);

    if (norms1_host(0) != ATM::one () + ATM::one ()) {
      curSuccess = false;
      out << "    FAILED: norms_host(0) = " << norms_host(0)
          << " != " << (ATM::one () + ATM::one ()) << endl;
    }
  }

  // Test that inf-norm of a column of ones is just one.
  // This distinguishes inf-norm from the other norms.
  out << "  Test that inf-norm of ones is just one" << endl;
  {
    KokkosBlas::fill (X, ATS::one ());
    KokkosBlas::nrmInf (norms, X);
    Kokkos::deep_copy (norms_host, norms);

    for (size_type j = 0; j < numCols; ++j) {
      if (norms_host(j) != ATM::one ()) {
        curSuccess = false;
        out << "    FAILED: norms_host(j) = " << norms_host(j)
            << " != " << ATM::one () << endl;
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess && prvSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testNorm1 (std::ostream& out, const bool prvSuccess)
{
  using Kokkos::subview;
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> MV;
  typedef typename MV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  typedef Kokkos::Details::InnerProductSpaceTraits<Scalar> IPT;
  typedef typename IPT::mag_type mag_type;
  typedef Kokkos::View<mag_type*, Device> norms_type;
  typedef Kokkos::Details::ArithTraits<mag_type> ATM;
  bool curSuccess = true;

  const size_type numRows = 10;
  const size_type numCols = 3;
  norms_type X_norms ("X_norms", numCols);
  typename norms_type::HostMirror X_norms_h =
    Kokkos::create_mirror_view (X_norms);

  out << "Testing 1-norm" << endl;

  MV X ("X", numRows, numCols);
  typename MV::HostMirror X_host = Kokkos::create_mirror_view (X);
  norms_type norms ("norms", numCols);
  typename norms_type::HostMirror norms_host = Kokkos::create_mirror_view (norms);

  // Test that the 1-norm of a column of -ones, is the number of rows
  // in the column.
  out << "  Test that the 1-norm of -ones, is the number of rows" << endl;
  {
    KokkosBlas::fill (X, -ATS::one ());
    KokkosBlas::nrm1 (norms, X);
    Kokkos::deep_copy (norms_host, norms);
    const mag_type expectedNorm = static_cast<mag_type> (numRows);
    for (size_type j = 0; j < numCols; ++j) {
      if (norms_host(j) != expectedNorm) {
        curSuccess = false;
        out << "    FAILED: norms_host(j) = " << norms_host(j)
            << " != " << expectedNorm << endl;
      }
    }
  }

  // Test that the 1-norm of [1, 2, 3, 4] is 10.
  // This distinguishes the 1-norm from the other norms.
  out << "  Test that the 1-norm of [1, 2, 3, 4] is 10" << endl;
  {
    const mag_type ONE = ATM::one ();
    const mag_type TWO = ONE + ONE;
    const mag_type THREE = TWO + ONE;
    const mag_type FOUR = THREE + ONE;
    const mag_type TEN = ONE + TWO + THREE + FOUR;

    MV X1 ("X1", 4, 1);
    typename MV::HostMirror X1_host = Kokkos::create_mirror_view (X1);
    X1_host(0,0) = ATS::one ();
    X1_host(1,0) = X1_host(0,0) + ATS::one ();
    X1_host(2,0) = X1_host(1,0) + ATS::one ();
    X1_host(3,0) = X1_host(2,0) + ATS::one ();
    Kokkos::deep_copy (X1, X1_host);

    norms_type norms1 ("norms1", 1);
    typename norms_type::HostMirror norms1_host =
      Kokkos::create_mirror_view (norms1);

    KokkosBlas::nrm1 (norms1, X1);
    Kokkos::deep_copy (norms1_host, norms1);

    if (norms1_host(0) != TEN) {
      curSuccess = false;
      out << "    FAILED: norms1_host(0) = " << norms1_host(0)
          << " != " << TEN << endl;
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess && prvSuccess;
}

template<class Scalar, class Layout, class Device>
bool
testMV (std::ostream& out, const bool prvSuccess)
{
  bool curSuccess = true;

  curSuccess = testFill<Scalar, Layout, Device> (out, prvSuccess);
  for (int wn = static_cast<int> (TPETRAKERNELS_TEST_NORM_INF);
       wn < static_cast<int> (TPETRAKERNELS_TEST_NORM_INVALID);
       ++wn) {
    const EWhichNorm whichNorm = static_cast<EWhichNorm> (wn);
    curSuccess = testAnyNorm<Scalar, Layout, Device> (out, whichNorm, curSuccess);
  }
  curSuccess = testNorm1<Scalar, Layout, Device> (out, curSuccess);
  curSuccess = testNormInf<Scalar, Layout, Device> (out, curSuccess);
  curSuccess = testSum<Scalar, Layout, Device> (out, curSuccess);
  return curSuccess && prvSuccess;
}


template<class Scalar, class Device>
bool
testOverLayouts (std::ostream& out, const bool prvSuccess)
{
  using std::endl;
  out << endl << "Testing Scalar = " << typeid (Scalar).name () << endl;
  bool curSuccess = true;

  out << endl << "Testing LayoutLeft" << endl;
  curSuccess = testMV<Scalar, Kokkos::LayoutLeft, Device> (out, curSuccess);
  out << endl << "Testing LayoutRight" << endl;
  curSuccess = testMV<Scalar, Kokkos::LayoutRight, Device> (out, curSuccess);

  return curSuccess && prvSuccess;
}

template<class Device>
bool
testOverScalarsAndLayouts (std::ostream& out, const bool prvSuccess)
{
  const bool testComplex = true;
  bool curSuccess = true;

  curSuccess = testOverLayouts<double, Device> (out, curSuccess);
  curSuccess = testOverLayouts<float, Device> (out, curSuccess);
  curSuccess = testOverLayouts<int, Device> (out, curSuccess);

  if (testComplex) {
    curSuccess = testOverLayouts<Kokkos::complex<float>, Device> (out, curSuccess);
    curSuccess = testOverLayouts<Kokkos::complex<double>, Device> (out, curSuccess);
  }
  return curSuccess && prvSuccess;
}

bool
testOverScalarsAndLayoutsAndDevices (std::ostream& out,
                                     const bool prvSuccess)
{
  using std::endl;
  bool curSuccess = true;

#ifdef KOKKOS_HAVE_SERIAL
  {
    typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> device_type;
    if (Kokkos::Serial::is_initialized ()) {
      out << endl << "Testing Serial" << endl;
      curSuccess = testOverScalarsAndLayouts<device_type> (out, curSuccess);
    } else {
      out << endl << "Serial NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
  {
    typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> device_type;
    if (Kokkos::OpenMP::is_initialized ()) {
      out << endl << "Testing OpenMP" << endl;
      curSuccess = testOverScalarsAndLayouts<device_type> (out, curSuccess);
    } else {
      out << endl << "OpenMP NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
  {
    typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> device_type;
    if (Kokkos::Threads::is_initialized ()) {
      out << endl << "Testing Threads" << endl;
      curSuccess = testOverScalarsAndLayouts<device_type> (out, curSuccess);
    } else {
      out << endl << "Threads NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
  {
    typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> device_type;
    if (Kokkos::Cuda::is_initialized ()) {
      out << endl << "Testing Cuda,CudaSpace" << endl;
      curSuccess = testOverScalarsAndLayouts<device_type> (out, curSuccess);
    } else {
      out << endl << "Cuda NOT initialized; skipping test" << endl;
    }
  }
  {
    typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> device_type;
    if (Kokkos::Cuda::is_initialized ()) {
      out << endl << "Testing Cuda,CudaUVMSpace" << endl;
      curSuccess = testOverScalarsAndLayouts<device_type> (out, curSuccess);
    } else {
      out << endl << "Cuda NOT initialized; skipping test" << endl;
    }
  }
#endif // KOKKOS_HAVE_CUDA

  return curSuccess && prvSuccess;
}


int
main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  Kokkos::initialize (argc, argv);

#ifdef HAVE_MPI
  RCP<const Comm<int> > comm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
  RCP<const Comm<int> > comm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI
  const int myRank = comm->getRank ();

  bool success = true;

  success = testOverScalarsAndLayoutsAndDevices (std::cerr, success);
  if (success) {
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
  } else {
    if (myRank == 0) {
      cout << "End Result: TEST FAILED" << endl;
    }
  }
  Kokkos::finalize ();
  return EXIT_SUCCESS;
}
