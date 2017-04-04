/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <Kokkos_Blas1_MV.hpp>

#ifndef KOKKOSKERNELS_BLAS1_MV_UNITTESTS_HPP
#define KOKKOSKERNELS_BLAS1_MV_UNITTESTS_HPP

#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace { // (anonymous)

enum EWhichNorm {
  KOKKOSKERNELS_TEST_NORM_INF = 0,
  KOKKOSKERNELS_TEST_NORM_ONE = 1,
  KOKKOSKERNELS_TEST_NORM_TWO = 2,
  KOKKOSKERNELS_TEST_NORM_INVALID = 3
};

template<class RV, class XMV>
void
kokkosNorm (const RV& norms, const XMV& X, const EWhichNorm whichNorm)
{
  if (whichNorm == KOKKOSKERNELS_TEST_NORM_INF) {
    KokkosBlas::nrmInf (norms, X);
  } else if (whichNorm == KOKKOSKERNELS_TEST_NORM_ONE) {
    KokkosBlas::nrm1 (norms, X);
  } else if (whichNorm == KOKKOSKERNELS_TEST_NORM_TWO) {
    KokkosBlas::nrm2_squared (norms, X);
  }
}

} // namespace (anonymous)


namespace KokkosBlas {
namespace Impl {

template<class Scalar, class Layout, class Device>
bool
testFill (std::ostream& out, const int theNumCols, const bool oneCol)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  bool curSuccess = true;

  out << "Testing KokkosBlas::fill" << endl;

  const size_type numRows = 10;
  const size_type numCols = static_cast<size_type> (theNumCols);

  mv_type X ("X", numRows, numCols);
  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);

  out << "  Test fill with zero" << endl;
  KokkosBlas::fill (X, ATS::zero ());
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != ATS::zero ()) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = "
            << X_h(i,j) << " != " << ATS::zero () << endl;
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
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << ATS::one () << endl;
      }
    }
  }

  if (oneCol) {
    out << "  Repeat previous test, one column at a time" << endl;

    out << "    Test fill with zero" << endl;
    for (size_type j = 0; j < numCols; ++j) {
      auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
      KokkosBlas::fill (X_j, ATS::zero ());
      Kokkos::deep_copy (X_h, X);
      for (size_type i = 0; i < numRows; ++i) {
        if (X_h(i,j) != ATS::zero ()) {
          curSuccess = false;
          out << "    FAILED: X_h(" << i << "," << j << ") = "
              << X_h(i,j) << " != " << ATS::zero () << endl;
        }
      }
    }

    out << "    Test fill with one" << endl;
    for (size_type j = 0; j < numCols; ++j) {
      auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
      KokkosBlas::fill (X_j, ATS::one ());
      Kokkos::deep_copy (X_h, X);
      for (size_type i = 0; i < numRows; ++i) {
        if (X_h(i,j) != ATS::one ()) {
          curSuccess = false;
          out << "    FAILED: X_h(" << i << "," << j << ") = "
              << X_h(i,j) << " != " << ATS::one () << endl;
        }
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testDot (std::ostream& out, const int theNumCols)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  typedef Kokkos::Details::InnerProductSpaceTraits<Scalar> IPT;
  typedef typename IPT::dot_type dot_type;
  typedef Kokkos::View<dot_type, Layout, Device> zero_d_result_type;
  typedef Kokkos::View<dot_type*, Layout, Device> one_d_result_type;

  bool curSuccess = true;

  out << "Testing KokkosBlas::dot" << endl;

  const size_type numRows = 10;
  const size_type numCols = static_cast<size_type> (theNumCols);

  mv_type X ("X", numRows, numCols);
  mv_type Y ("Y", numRows, numCols);
  one_d_result_type R ("R", numCols);
  one_d_result_type R_expected ("R_expected", numCols);
  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);
  typename mv_type::HostMirror Y_h = Kokkos::create_mirror_view (Y);
  typename one_d_result_type::HostMirror R_h = Kokkos::create_mirror_view (R);
  typename one_d_result_type::HostMirror R_expected_h ("R_expected_h", numCols);

  const Scalar ONE = ATS::one ();
  const Scalar TWO = ONE + ONE;

  KokkosBlas::fill (R, Kokkos::Details::ArithTraits<dot_type>::zero ());
  Scalar curVal = ONE;
  for (size_type j = 0; j < numCols; ++j) {
    auto X_h_j = Kokkos::subview (X_h, Kokkos::ALL (), j);
    auto Y_h_j = Kokkos::subview (Y_h, Kokkos::ALL (), j);
    KokkosBlas::fill (X_h_j, curVal);
    KokkosBlas::fill (Y_h_j, TWO * curVal);

    dot_type curExpectedVal = Kokkos::Details::ArithTraits<dot_type>::zero ();
    for (size_type i = 0; i < numRows; ++i) {
      curExpectedVal += IPT::dot (curVal, TWO * curVal);
    }
    R_expected_h[j] = curExpectedVal;
    curVal += ONE;
  }
  Kokkos::deep_copy (X, X_h);
  Kokkos::deep_copy (Y, Y_h);

  out << "  Test 2-D dot 2-D -> 1-D case" << endl;
  KokkosBlas::dot (R, X, Y);
  Kokkos::deep_copy (R_h, R);

  for (size_type j = 0; j < numCols; ++j) {
    // We chose both the values and the number of terms in each sum
    // such that rounding error should not occur.  Thus, we may
    // compare exactly.
    if (R_expected_h[j] != R_h[j]) {
      curSuccess = false;
    }
  }

  out << "  Test 1-D dot 1-D -> 0-D case" << endl;
  zero_d_result_type r ("r");
  typename zero_d_result_type::HostMirror r_h = Kokkos::create_mirror_view (r);
  for (size_type j = 0; j < numCols; ++j) {
    auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
    auto Y_j = Kokkos::subview (Y, Kokkos::ALL (), j);
    KokkosBlas::dot (r, X_j, Y_j);
    Kokkos::deep_copy (r_h, r);
    if (r_h () != R_expected_h[j]) {
      curSuccess = false;
    }
  }

  out << "  Test 2-D dot 1-D -> 1-D case" << endl;
  one_d_result_type R2 ("R2", numCols);
  typename one_d_result_type::HostMirror R2_h = Kokkos::create_mirror_view (R2);
  curVal = ONE;
  for (size_type j = 0; j < numCols; ++j) {
    // auto X_h_j = Kokkos::subview (X_h, Kokkos::ALL (), j);
    // auto Y_h_j = Kokkos::subview (Y_h, Kokkos::ALL (), j);
    // KokkosBlas::fill (X_h_j, curVal);
    // KokkosBlas::fill (Y_h_j, TWO * curVal);

    dot_type curExpectedVal = Kokkos::Details::ArithTraits<dot_type>::zero ();
    for (size_type i = 0; i < numRows; ++i) {
      curExpectedVal += IPT::dot (curVal, TWO);
    }
    R_expected_h[j] = curExpectedVal;
    curVal += ONE;
  }
  auto Y_0 = Kokkos::subview (Y, Kokkos::ALL (), 0);
  KokkosBlas::dot (R2, X, Y_0);
  Kokkos::deep_copy (R2_h, R2);
  for (size_type j = 0; j < numCols; ++j) {
    if (R2_h[j] != R_expected_h[j]) {
      curSuccess = false;
    }
  }

  out << "  Test 1-D dot 2-D -> 1-D case" << endl;
  Kokkos::deep_copy (R2, Kokkos::Details::ArithTraits<dot_type>::zero ());
  const Scalar THREE = TWO + ONE; // need this, else same as previous test
  curVal = THREE;
  for (size_type j = 0; j < numCols; ++j) {
    // auto X_h_j = Kokkos::subview (X_h, Kokkos::ALL (), j);
    // auto Y_h_j = Kokkos::subview (Y_h, Kokkos::ALL (), j);
    // KokkosBlas::fill (X_h_j, curVal);
    // KokkosBlas::fill (Y_h_j, TWO * curVal);

    dot_type curExpectedVal = Kokkos::Details::ArithTraits<dot_type>::zero ();
    for (size_type i = 0; i < numRows; ++i) {
      curExpectedVal += IPT::dot (THREE, TWO * curVal);
    }
    R_expected_h[j] = curExpectedVal;
    curVal += ONE;
  }
  auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
  KokkosBlas::dot (R2, X_0, Y);
  Kokkos::deep_copy (R2_h, R2);
  for (size_type j = 0; j < numCols; ++j) {
    if (R2_h[j] != R_expected_h[j]) {
      curSuccess = false;
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}



template<class Scalar, class Layout, class Device>
bool
testAxpby (std::ostream& out, const int theNumCols)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  bool curSuccess = true;

  out << "Testing KokkosBlas::axpby" << endl;

  const size_type numRows = 10;
  const size_type numCols = static_cast<size_type> (theNumCols);

  mv_type X ("X", numRows, numCols);
  //typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);
  mv_type Y ("Y", numRows, numCols);
  typename mv_type::HostMirror Y_h = Kokkos::create_mirror_view (Y);

  const Scalar ZERO = ATS::zero ();
  const Scalar ONE = ATS::one ();
  const Scalar TWO = ONE + ONE;
  const Scalar THREE = TWO + ONE;
  const Scalar FOUR = THREE + ONE;
  const Scalar FIVE = FOUR + ONE;
  const Scalar EIGHT = FOUR + FOUR;
  const Scalar TEN = FIVE + FIVE;
  const Scalar TWELVE = TEN + TWO;
  const Scalar TWENTY_TWO = TWELVE + TEN;

  Scalar alpha;
  Scalar beta;

  KokkosBlas::fill (X, ONE);
  KokkosBlas::fill (Y, TWO);
  alpha = ZERO;
  beta = ZERO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != ZERO) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << ZERO << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = ONE;
  beta = ZERO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != FOUR) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << FOUR << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = -ONE;
  beta = ZERO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != -FOUR) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << -FOUR << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = TWO;
  beta = ZERO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != EIGHT) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << EIGHT << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = ZERO;
  beta = ONE;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != FIVE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << FIVE << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = ZERO;
  beta = -ONE;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != -FIVE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << -FIVE << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = ZERO;
  beta = TWO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != TEN) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << TEN << endl;
      }
    }
  }

  KokkosBlas::fill (X, TWO);
  KokkosBlas::fill (Y, THREE);
  alpha = ONE;
  beta = ONE;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != FIVE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << FIVE << endl;
      }
    }
  }

  KokkosBlas::fill (X, TWO);
  KokkosBlas::fill (Y, THREE);
  alpha = -ONE;
  beta = ONE;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != ONE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << ONE << endl;
      }
    }
  }

  KokkosBlas::fill (X, TWO);
  KokkosBlas::fill (Y, THREE);
  alpha = ONE;
  beta = -ONE;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != -ONE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << -ONE << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = THREE;
  beta = ZERO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != TWELVE) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << TWELVE << endl;
      }
    }
  }

  KokkosBlas::fill (X, FOUR);
  KokkosBlas::fill (Y, FIVE);
  alpha = THREE;
  beta = TWO;
  out << "  Test axpby(" << alpha << ", X, " << beta << ", Y)" << endl;

  KokkosBlas::axpby (alpha, X, beta, Y);
  Kokkos::deep_copy (Y_h, Y);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (Y_h(i,j) != TWENTY_TWO) {
        curSuccess = false;
        out << "    FAILED: Y_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << TWENTY_TWO << endl;
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testScal1Arg (std::ostream& out, const int theNumCols)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  bool curSuccess = true;

  out << "Testing KokkosBlas::axpby" << endl;

  const size_type numRows = 4;
  const size_type numCols = static_cast<size_type> (theNumCols);

  mv_type X ("X", numRows, numCols);
  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view (X);

  const Scalar ZERO = ATS::zero ();
  const Scalar ONE = ATS::one ();
  const Scalar TWO = ONE + ONE;
  Scalar alpha;

  // Set up the test problem.
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      X_h(i,j) = static_cast<Scalar> (i + 1);
    }
  }
  Kokkos::deep_copy (X, X_h);

  // Test with alpha == 0.
  alpha = ZERO;
  out << "  Test scal(X, " << alpha << ", X)" << endl;
  KokkosBlas::scal (X, alpha, X);
  // Compare against the right answer.
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != ZERO) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << ZERO << endl;
      }
    }
  }

  // Set up the test problem.
  for (size_type j = 0; j < numCols; ++j) {
    if (numRows > static_cast<size_type> (0)) {
      X_h(0,j) = ONE;
    }
    for (size_type i = 1; i < numRows; ++i) {
      X_h(i,j) = X_h(i-1,j) + ONE;
    }
  }
  Kokkos::deep_copy (X, X_h);

  // Test with alpha == 1.
  alpha = ONE;
  out << "  Test scal(X, " << alpha << ", X)" << endl;
  KokkosBlas::scal (X, alpha, X);
  // Compare against the right answer.
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != static_cast<Scalar> (i + 1)) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << static_cast<Scalar> (i+1) << endl;
      }
    }
  }

  // Set up the test problem.
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      X_h(i,j) = static_cast<Scalar> (i + 1);
    }
  }
  Kokkos::deep_copy (X, X_h);

  // Test with alpha == -1.
  alpha = -ONE;
  out << "  Test scal(X, " << alpha << ", X)" << endl;
  KokkosBlas::scal (X, alpha, X);
  // Compare against the right answer.
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != -static_cast<Scalar> (i + 1)) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << static_cast<Scalar> (i+1) << endl;
      }
    }
  }

  // Set up the test problem.
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      X_h(i,j) = static_cast<Scalar> (i + 1);
    }
  }
  Kokkos::deep_copy (X, X_h);

  // Test with alpha == 2.
  alpha = TWO;
  out << "  Test scal(X, " << alpha << ", X)" << endl;
  KokkosBlas::scal (X, alpha, X);
  // Compare against the right answer.
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      if (X_h(i,j) != TWO * static_cast<Scalar> (i + 1)) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << TWO * static_cast<Scalar> (i+1) << endl;
      }
    }
  }

  // Set up the test problem.
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      X_h(i,j) = static_cast<Scalar> (i + 1);
    }
  }
  Kokkos::deep_copy (X, X_h);

  // Test with alpha == -2.
  alpha = -TWO;
  out << "  Test scal(X, " << alpha << ", X)" << endl;
  KokkosBlas::scal (X, alpha, X);
  // Compare against the right answer.
  Kokkos::deep_copy (X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      const Scalar should = -TWO * static_cast<Scalar> (i + 1);
      if (X_h(i,j) != should) {
        curSuccess = false;
        out << "    FAILED: X_h(" << i << "," << j << ") = " << X_h(i,j)
            << " != " << should << endl;
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}

template<class Scalar, class Layout, class Device>
bool
testSum (std::ostream& out, const int theNumCols, const bool oneCol)
{
  using std::endl;
  typedef Kokkos::View<Scalar**, Layout, Device> mv_type;
  // sum() uses the Device's preferred Layout for the output array.
  typedef Kokkos::View<Scalar*, Device> sums_type;
  typedef typename mv_type::size_type size_type;
  typedef Kokkos::Details::ArithTraits<Scalar> ATS;
  const Scalar ZERO = ATS::zero ();
  const Scalar ONE = ATS::one ();
  const Scalar TWO = ONE + ONE;
  const Scalar THREE = TWO + ONE;
  const Scalar FOUR = THREE + ONE;
  const Scalar TEN = ONE + TWO + THREE + FOUR;
  bool curSuccess = true;

  out << "Testing KokkosBlas::sum" << endl;

  const size_type numRows = 4;
  const size_type numCols = static_cast<size_type> (theNumCols);

  mv_type X ("X", numRows, numCols);
  sums_type r ("r", numCols);

  if (X.dimension_0 () != numRows || X.dimension_1 () != numCols) {
    out << "  FAILED: X should be " << numRows << " x " << numCols
        << ", but is " << X.dimension_0 () << " x " << X.dimension_1 ()
        << " instead!" << endl;
    return false;
  }
  if (r.dimension_0 () != numCols) {
    out << "  FAILED: r should be " << numRows << " x 1, but is "
        << r.dimension_0 () << " x 1 instead!" << endl;
    return false;
  }

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

  if (oneCol) {
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
      if (r_h(j) != ZERO) {
        curSuccess = false;
        out << "    FAILED: r_h(" << j << ") = " << r_h(j)
            << " != " << ZERO << endl;
      }
    }
  }

  {
    out << "  Test the sum of ones" << endl;
    Scalar should = ZERO;
    for (size_type i = 0; i < numRows; ++i) {
      should += ONE;
    }
    KokkosBlas::fill (X, ONE);
    KokkosBlas::sum (r, X);
    Kokkos::deep_copy (r_h, r);
    for (size_type j = 0; j < numCols; ++j) {
      if (r_h(j) != should) {
        curSuccess = false;
        out << "    FAILED: r_h(" << j << ") = " << r_h(j)
            << " != " << should << endl;
      }
    }
  }

  if (oneCol) {
    out << "  Repeat previous test, one column at a time" << endl;
    Scalar should = ZERO;
    for (size_type i = 0; i < numRows; ++i) {
      should += ONE;
    }
    for (size_type j = 0; j < numCols; ++j) {
      auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
      auto r_j = Kokkos::subview (r, j);
      KokkosBlas::fill (X_j, ONE);
      KokkosBlas::sum (r_j, X_j);
    }
    Kokkos::deep_copy (r_h, r);
    for (size_type j = 0; j < numCols; ++j) {
      if (r_h(j) != should) {
        curSuccess = false;
        out << "    FAILED: r_h(" << j << ") = " << r_h(j)
            << " != " << should << endl;
      }
    }
  }

  // Tetractys test.
  out << "  Test that the sum of [1, 2, 3, 4] is 10" << endl;
  for (size_type j = 0; j < numCols; ++j) {
    X_h(0,j) = ONE;
    X_h(1,j) = TWO;
    X_h(2,j) = THREE;
    X_h(3,j) = FOUR;
  }
  Kokkos::deep_copy (X, X_h);
  KokkosBlas::sum (r, X);
  Kokkos::deep_copy (r_h, r);
  for (size_type j = 0; j < numCols; ++j) {
    if (r_h(j) != TEN) {
      curSuccess = false;
      out << "    FAILED: r_h(" << j << ") = " << r_h(j)
          << " != " << TEN << endl;
    }
  }

  if (oneCol) {
    out << "  Repeat previous test, one column at a time" << endl;
    KokkosBlas::fill (r, ATS::zero ()); // reset r
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

  if (oneCol) {
    out << "  Repeat previous test, one column at a time" << endl;
    KokkosBlas::fill (r, ATS::zero ()); // reset r
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
      if (r_h(j) != TWO) {
        curSuccess = false;
        out << "    FAILED: r_h(" << j << ") = " << r_h(j)
            << " != " << TWO << endl;
      }
    }
  }

  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testAnyNorm (std::ostream& out, const EWhichNorm whichNorm,
             const int theNumCols)
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
  const size_type numCols = static_cast<size_type> (theNumCols);
  norms_type X_norms ("X_norms", numCols);
  typename norms_type::HostMirror X_norms_h = Kokkos::create_mirror_view (X_norms);

  out << "Testing properties common to all norms, for norm ";
  if (whichNorm == KOKKOSKERNELS_TEST_NORM_INF) {
    out << "inf";
  } else if (whichNorm == KOKKOSKERNELS_TEST_NORM_ONE) {
    out << "1";
  } else if (whichNorm == KOKKOSKERNELS_TEST_NORM_TWO) {
    out << "2";
  }
  out << endl;

  if (numCols != 1) {
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
      out << "    FAILED: Norm of a nondegenerate single-column MV "
        "full of zeros was not zero" << endl;
    }

    kokkosNorm<norms_type, MV> (subview (X_norms, pair_type (0, numCols)), X2, whichNorm);
    Kokkos::deep_copy (X_norms_h, X_norms);
    for (size_type j = 0; j < numCols; ++j) {
      if (X_norms_h(j) != ATM::zero ()) {
        curSuccess = false;
        out << "    FAILED: Norm of a nondegenerate multicolumn MV "
          "full of zeros was not zero" << endl;
        break;
      }
    }
  }
  if (curSuccess) {
    out << "  SUCCESS" << endl;
  } else {
    out << "  FAILURE" << endl;
  }
  return curSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testNormInf (std::ostream& out, const int theNumCols)
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
  const size_type numCols = static_cast<size_type> (theNumCols);
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
  return curSuccess;
}


template<class Scalar, class Layout, class Device>
bool
testNorm1 (std::ostream& out, const int theNumCols)
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
  const size_type numCols = static_cast<size_type> (theNumCols);
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
  return curSuccess;
}

template<class Scalar, class Layout, class Device>
bool
testMV (std::ostream& out, const int numCols, const bool oneCol)
{
  bool success = true;
  bool curSuccess = true;

  curSuccess = testFill<Scalar, Layout, Device> (out, numCols, oneCol);
  success = success && curSuccess;

  curSuccess = testDot<Scalar, Layout, Device> (out, numCols);
  success = success && curSuccess;

  // The compiler frowns upon looping with enums, so we cast to int.
  for (int wn = static_cast<int> (KOKKOSKERNELS_TEST_NORM_INF);
       wn < static_cast<int> (KOKKOSKERNELS_TEST_NORM_INVALID);
       ++wn) {
    const EWhichNorm whichNorm = static_cast<EWhichNorm> (wn);
    curSuccess = testAnyNorm<Scalar, Layout, Device> (out, whichNorm, numCols);
    success = success && curSuccess;
  }
  curSuccess = testNorm1<Scalar, Layout, Device> (out, numCols);
  success = success && curSuccess;
  curSuccess = testNormInf<Scalar, Layout, Device> (out, numCols);
  success = success && curSuccess;
  curSuccess = testSum<Scalar, Layout, Device> (out, numCols, oneCol);
  success = success && curSuccess;
  curSuccess = testScal1Arg<Scalar, Layout, Device> (out, numCols);
  success = success && curSuccess;

  return success;
}

template<class Scalar, class Device>
bool
testOverLayouts (std::ostream& out, const int numCols, const bool oneCol)
{
  using std::endl;
  out << endl << "Testing Scalar = " << typeid (Scalar).name () << endl;
  bool curSuccess = true;
  bool success = true;

  out << endl << "Testing LayoutLeft" << endl;
  curSuccess = testMV<Scalar, Kokkos::LayoutLeft, Device> (out, numCols, oneCol);
  success = success && curSuccess;
  out << endl << "Testing LayoutRight" << endl;
  curSuccess = testMV<Scalar, Kokkos::LayoutRight, Device> (out, numCols, oneCol);
  success = success && curSuccess;

  return success;
}

template<class Device>
bool
testOverScalarsAndLayouts (std::ostream& out,
                           const char deviceName[],
                           const int numCols,
                           const bool oneCol,
                           const bool testComplex)
{
  using std::endl;
  typedef Device device_type;
  bool success = true;

  if (device_type::execution_space::is_initialized ()) {
    out << endl << "Testing Device = " << deviceName << endl;
    bool curSuccess = true;

    curSuccess = testOverLayouts<double, device_type> (out, numCols, oneCol);
    success = success && curSuccess;
    curSuccess = testOverLayouts<float, device_type> (out, numCols, oneCol);
    success = success && curSuccess;
    curSuccess = testOverLayouts<int, device_type> (out, numCols, oneCol);
    success = success && curSuccess;

    if (testComplex) {
      curSuccess = testOverLayouts<Kokkos::complex<float>, device_type> (out, numCols, oneCol);
      success = success && curSuccess;
      curSuccess = testOverLayouts<Kokkos::complex<double>, device_type> (out, numCols, oneCol);
      success = success && curSuccess;
    }
  }
  else {
    out << "Device = " << deviceName << " NOT initialized; skipping test" << endl;
  }

  return success;
}

//
// Declarations of full specializations of testOverScalarsAndLayouts.
//
// We _define_ each full specialization in its own .cpp file.  This
// reduces compilation time.  All full specializations are the same as
// the generic version of testOverScalarsAndLayouts above, though one
// could make them differ if necessary (for example, some Scalar types
// might only work on certain Devices).
//

#ifdef KOKKOS_HAVE_SERIAL
template<>
bool
testOverScalarsAndLayouts<Kokkos::Device<Kokkos::Serial,
                                         Kokkos::HostSpace> > (std::ostream& out,
                                                               const char deviceName[],
                                                               const int numCols,
                                                               const bool oneCol,
                                                               const bool testComplex);
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
bool
testOverScalarsAndLayouts<Kokkos::Device<Kokkos::OpenMP,
                                         Kokkos::HostSpace> > (std::ostream& out,
                                                               const char deviceName[],
                                                               const int numCols,
                                                               const bool oneCol,
                                                               const bool testComplex);
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
bool
testOverScalarsAndLayouts<Kokkos::Device<Kokkos::Threads,
                                         Kokkos::HostSpace> > (std::ostream& out,
                                                               const char deviceName[],
                                                               const int numCols,
                                                               const bool oneCol,
                                                               const bool testComplex);
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
bool
testOverScalarsAndLayouts<Kokkos::Device<Kokkos::Cuda,
                                         Kokkos::CudaSpace> > (std::ostream& out,
                                                               const char deviceName[],
                                                               const int numCols,
                                                               const bool oneCol,
                                                               const bool testComplex);
template<>
bool
testOverScalarsAndLayouts<Kokkos::Device<Kokkos::Cuda,
                                         Kokkos::CudaUVMSpace> > (std::ostream& out,
                                                                  const char deviceName[],
                                                                  const int numCols,
                                                                  const bool oneCol,
                                                                  const bool testComplex);
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOSKERNELS_BLAS1_MV_UNITTESTS_HPP
