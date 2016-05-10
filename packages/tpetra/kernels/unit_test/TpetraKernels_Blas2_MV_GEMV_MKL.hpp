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

#ifndef TPETRAKERNELS_BLAS2_MV_GEMV_MKL_HPP
#define TPETRAKERNELS_BLAS2_MV_GEMV_MKL_HPP

#include "Kokkos_Blas1_MV.hpp"
#include "Kokkos_Blas2_MV.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include <sstream>
#include <stdexcept>

#ifdef HAVE_TPETRAKERNELS_MKL
#  include "mkl.h"
#endif // HAVE_TPETRAKERNELS_MKL

namespace { // (anonymous)

using std::endl;

// Serial version of gemv (dense matrix-vector multiply).
template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType = typename YViewType::non_const_value_type,
         class BetaCoeffType = typename YViewType::non_const_value_type,
         class IndexType = typename AViewType::size_type>
struct SerialGEMV {
  static void
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    using Kokkos::Details::ArithTraits;
    typedef typename std::decay<decltype (A(0,0)) >::type A_value_type;
    typedef typename std::decay<decltype (y(0)) >::type y_value_type;

    const IndexType A_numRows = A.dimension_0 ();
    const IndexType A_numCols = A.dimension_1 ();

    const bool conjugate = trans[0] == 'C' || trans[0] == 'c' ||
      trans[0] == 'H' || trans[0] == 'h';
    const bool transpose = trans[0] != 'N' && trans[0] != 'n';

    if (transpose) {
      for (IndexType j = 0; j < A_numCols; ++j) {
        y_value_type y_j = ArithTraits<y_value_type>::zero ();

        if (alpha != ArithTraits<AlphaCoeffType>::zero ()) {
          for (IndexType i = 0; i < A_numRows; ++i) {
            if (conjugate) {
              y_j += alpha * ArithTraits<A_value_type>::conj (A(i, j)) * x(i);
            }
            else {
              y_j += alpha * A(i, j) * x(i);
            }
          }
        }

        const y_value_type y_j_orig =
          beta == ArithTraits<BetaCoeffType>::zero () ?
          ArithTraits<y_value_type>::zero () :
          beta * y[j];
        y[j] = y_j + y_j_orig;
      }
    }
    else { // not transpose
      for (IndexType i = 0; i < A_numRows; ++i) {
        y_value_type y_i = ArithTraits<y_value_type>::zero ();

        if (alpha != ArithTraits<AlphaCoeffType>::zero ()) {
          for (IndexType j = 0; j < A_numCols; ++j) {
            y_i += alpha * A(i, j) * x(j);
          }
        }

        if (beta == ArithTraits<BetaCoeffType>::zero ()) {
          y(i) = y_i;
        }
        else {
          y(i) = beta * y(i) + y_i;
        }
      }
    }
  }
};

// Serial version of gemv (dense matrix-vector multiply).
template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType>
void
gemv (const char trans[],
      const AlphaCoeffType& alpha,
      const AViewType& A,
      const XViewType& x,
      const BetaCoeffType& beta,
      const YViewType& y)
{
  const bool transpose = trans[0] != 'N' && trans[0] != 'n';
  if ((transpose && (A.dimension_1 () != y.dimension_0 () ||
                     A.dimension_0 () != x.dimension_0 ())) ||
      (! transpose && (A.dimension_0 () != y.dimension_0 () ||
                       A.dimension_1 () != x.dimension_0 ()))) {
    std::ostringstream os;
    os << "gemv: Dimensions don't match. A is " << A.dimension_0 ()
       << " x " << A.dimension_1 () << ", x is " << x.dimension_0 ()
       << ", and y is " << y.dimension_0 () << ".";
    throw std::invalid_argument (os.str ());
  }
  SerialGEMV<AViewType, XViewType, YViewType, AlphaCoeffType,
    BetaCoeffType>::gemv (trans, alpha, A, x, beta, y);
}

//
// Run a single test of KokkosBlas::Impl::tryMklGemv.
//
template<class Scalar, class Layout, class Device>
bool
testGemvOne (std::ostream& curOut,
             std::string indent,
             const typename Kokkos::View<Scalar*, Layout, Device>::size_type numRows,
             const typename Kokkos::View<Scalar*, Layout, Device>::size_type numCols,
             const char trans,
             const Scalar& alpha,
             const Scalar& beta,
             const bool printOnlyOnFailure,
             const bool debug)
{
  using Kokkos::Details::ArithTraits;
  using Kokkos::Details::InnerProductSpaceTraits;
  typedef Kokkos::View<Scalar**, Layout, Device> matrix_type;
  typedef Kokkos::View<Scalar*, Layout, Device> vector_type;
  typedef typename vector_type::size_type size_type;
  typedef ArithTraits<Scalar> ATS;
  typedef typename InnerProductSpaceTraits<Scalar>::mag_type mag_type;
  const Scalar one = ATS::one ();
  const Scalar two = one + one;
  const Scalar three = two + one;
  const Scalar four = two + two;
  const bool transpose = (trans != 'N' && trans != 'n');

  bool curSuccess = true;
  indent = indent + "  ";

  std::ostringstream os;
  // Gather non-debug output into the above ostringstream first.  This
  // makes the test output easier to read.  Always print debug output
  // immediately, since that lets us catch unexpected program
  // termination.
  os << "numRows: " << numRows
     << ", numCols: " << numCols
     << ", trans: '" << trans << "'"
     << ", alpha: " << alpha
     << ", beta: " << beta;

  Scalar curVal;

  matrix_type A ("A", numRows, numCols);
  auto A_h = Kokkos::create_mirror_view (A);
  curVal = one / two;
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      A_h(i,j) = curVal;
      curVal += (one / four);
    }
  }
  Kokkos::deep_copy (A, A_h);

  const size_type x_numEnt = transpose ? numRows : numCols;
  vector_type x ("x", x_numEnt);
  auto x_h = Kokkos::create_mirror_view (x);
  curVal = one / two;
  for (size_type j = 0; j < x_numEnt; ++j) {
    x_h(j) = curVal;
    curVal += (one / four);
  }
  Kokkos::deep_copy (x, x_h);

  const size_type y_numEnt = transpose ? numCols : numRows;
  vector_type y ("y", y_numEnt);
  auto y_h = Kokkos::create_mirror_view (y);

  curVal = three;
  for (size_type i = 0; i < y_numEnt; ++i) {
    y_h(i) = curVal;
    curVal += one;
  }
  Kokkos::deep_copy (y, y_h);

  if (debug) {
    // Debug output goes straight to the output stream, even if we're
    // supposed to be holding output until after the test runs.  This
    // is because the point of debug output is to catch unexpected
    // stopping of execution, like segfaults or thrown exceptions.
    curOut << indent << "Call KokkosBlas::gemv" << endl;
  }
  KokkosBlas::gemv (&trans, alpha, A, x, beta, y);

  if (debug) {
    curOut << indent << "Call comparison function" << endl;
  }
  vector_type y_copy ("y_copy", y_numEnt);
  auto y_copy_h = Kokkos::create_mirror_view (y_copy);

  curVal = three;
  for (size_type i = 0; i < y_numEnt; ++i) {
    y_copy_h(i) = curVal;
    curVal += one;
  }
  Kokkos::deep_copy (y_copy, y_copy_h);

  const bool callSucceeded =
    KokkosBlas::Impl::tryMklGemv (&trans, alpha, A, x, beta, y_copy);
  if (debug && ! callSucceeded) {
    curOut << indent << "FAILURE: KokkosBlas::Impl::tryMklGemv returned false" << endl;
  }
  curSuccess = curSuccess && callSucceeded;

  if (debug) {
    curOut << indent << "Compare KokkosBlas::Impl::tryMklGemv result to "
      "hand-rolled code result" << endl;
  }
  Kokkos::deep_copy (y_h, y);
  Kokkos::deep_copy (y_copy_h, y_copy);

  // Accumulate the max-norm difference.
  mag_type infNorm = ArithTraits<mag_type>::zero ();
  for (size_type j = 0; j < y_numEnt; ++j) {
    const mag_type diff =
      InnerProductSpaceTraits<Scalar>::norm (y_h(j) - y_copy_h(j));
    infNorm = std::max (infNorm, diff);
  }

  const mag_type tol = (transpose ? numRows : numCols) * ArithTraits<Scalar>::eps ();
  if (infNorm > tol) {
    curSuccess = false;
  }

  // It looks a little weird to compare tryMklGemv against the MKL.
  // However, this is a reasonable test for correct MKL CBLAS use
  // inside tryMklGemv.
  //
  // Compare against the MKL's CBLAS, but only for the array layouts
  // that the MKL supports (row major and column major).  We only test
  // nonzero numbers of rows and columns here, to avoid LDA confusion
  // for those cases.

  bool ranMklTest = false;
  mag_type infNormMkl = ArithTraits<mag_type>::zero ();
#ifdef HAVE_TPETRAKERNELS_MKL
  if (numRows != 0 && numCols != 0 &&
      (typeid (Layout) == typeid (Kokkos::LayoutLeft) ||
       typeid (Layout) == typeid (Kokkos::LayoutRight))) {
    const CBLAS_LAYOUT layoutMkl =
      typeid (Layout) == typeid (Kokkos::LayoutLeft) ?
      CblasColMajor : CblasRowMajor;
    const bool conjugate = (trans == 'C' || trans == 'c' || trans == 'H' || trans == 'h');

    CBLAS_TRANSPOSE transCblas = CblasNoTrans;
    if (transpose) {
      if (conjugate) {
        transCblas = CblasConjTrans;
      }
      else {
        transCblas = CblasTrans;
      }
    }

    const MKL_INT numRowsMkl = static_cast<MKL_INT> (numRows);
    const MKL_INT numColsMkl = static_cast<MKL_INT> (numCols);

    MKL_INT strides[8];
    A.stride (strides);
    const MKL_INT ldaMkl = (layoutMkl == CblasColMajor) ? strides[1] : strides[0];
    const MKL_INT incxMkl = 1;
    const MKL_INT incyMkl = 1;

    // Refill y_copy, so that we can use it in the test.
    curVal = three;
    for (size_type i = 0; i < y_numEnt; ++i) {
      y_copy_h(i) = curVal;
      curVal += one;
    }
    Kokkos::deep_copy (y_copy, y_copy_h);

    // Skip unsupported Scalar types.  We have to do reinterpret_cast
    // because this code must compile for every Scalar type, even if
    // it never gets called.  This also explains why we take the real
    // parts of alpha and beta.
    if (typeid (Scalar) == typeid (double)) {
      cblas_dgemv (layoutMkl, transCblas, numRowsMkl, numColsMkl,
                   ATS::real (alpha),
                   reinterpret_cast<const double*> (A.ptr_on_device ()), ldaMkl,
                   reinterpret_cast<const double*> (x.ptr_on_device ()),
                   incxMkl, ATS::real (beta),
                   reinterpret_cast<double*> (y_copy.ptr_on_device ()),
                   incyMkl);
      ranMklTest = true;
    }
    else if (typeid (Scalar) == typeid (float)) {
      cblas_sgemv (layoutMkl, transCblas, numRowsMkl, numColsMkl,
                   ATS::real (alpha),
                   reinterpret_cast<const float*> (A.ptr_on_device ()), ldaMkl,
                   reinterpret_cast<const float*> (x.ptr_on_device ()),
                   incxMkl, ATS::real (beta),
                   reinterpret_cast<float*> (y_copy.ptr_on_device ()),
                   incyMkl);
      ranMklTest = true;
    }
    else if (typeid (Scalar) == typeid (Kokkos::complex<float>)) {
      cblas_cgemv (layoutMkl, transCblas, numRowsMkl, numColsMkl,
                   reinterpret_cast<const void*> (&alpha),
                   reinterpret_cast<const void*> (A.ptr_on_device ()), ldaMkl,
                   reinterpret_cast<const void*> (x.ptr_on_device ()),
                   incxMkl, reinterpret_cast<const void*> (&beta),
                   reinterpret_cast<void*> (y_copy.ptr_on_device ()),
                   incyMkl);
      ranMklTest = true;
    }
    else if (typeid (Scalar) == typeid (Kokkos::complex<double>)) {
      cblas_zgemv (layoutMkl, transCblas, numRowsMkl, numColsMkl,
                   reinterpret_cast<const void*> (&alpha),
                   reinterpret_cast<const void*> (A.ptr_on_device ()), ldaMkl,
                   reinterpret_cast<const void*> (x.ptr_on_device ()),
                   incxMkl, reinterpret_cast<const void*> (&beta),
                   reinterpret_cast<void*> (y_copy.ptr_on_device ()),
                   incyMkl);
      ranMklTest = true;
    }

    if (ranMklTest) {
      if (debug) {
        curOut << indent << "Compare KokkosBlas::Impl::tryMklGemv result to "
          "MKL CBLAS result" << endl;
      }
      Kokkos::deep_copy (y_copy_h, y_copy);

      // Accumulate the max-norm difference.
      infNormMkl = ArithTraits<mag_type>::zero ();
      for (size_type j = 0; j < y_numEnt; ++j) {
        const mag_type diff =
          InnerProductSpaceTraits<Scalar>::norm (y_h(j) - y_copy_h(j));
        infNormMkl = std::max (infNorm, diff);
      }
    }
  }
#endif // HAVE_TPETRAKERNELS_MKL

  if (callSucceeded) {
    os << ", infNorm (hand): " << infNorm << (infNorm <= tol ? " <= " : " > ") << tol;
    if (ranMklTest) {
      os << ", infNorm (MKL): " << infNormMkl << (infNormMkl <= tol ? " <= " : " > ") << tol;
    }
  }
  else {
    os << ", tryMklGemv returned false (failed)";
  }

  if (curSuccess) {
    if (! printOnlyOnFailure) {
      curOut << indent << "SUCCESS: {" << os.str () << "}" << endl;
    }
  }
  else if (! curSuccess) {
    curOut << indent << "FAILURE: {" << os.str () << "}" << endl;
  }

  return curSuccess;
}

template<class Scalar, class Layout, class Device>
bool
testGemvAll (std::ostream& out,
             std::string indent,
             const bool printOnlyOnFailure,
             const bool debug)
{
  using Kokkos::Details::ArithTraits;
  typedef typename Kokkos::View<Scalar*, Layout, Device>::size_type size_type;
  bool success = true;
  bool curSuccess = true;

  const Scalar zero = ArithTraits<Scalar>::zero ();
  const Scalar one = ArithTraits<Scalar>::one ();
  const Scalar three = one + one + one;

  const int numCoeffValues = 5;
  Scalar coeffValues[5];
  coeffValues[0] = zero;
  coeffValues[1] = one;
  coeffValues[2] = -one;
  coeffValues[3] = three;
  coeffValues[4] = -three;

  const int numNumRowsValues = 7;
  size_type numRowsValues[7];
  numRowsValues[0] = 0;
  numRowsValues[1] = 1;
  numRowsValues[2] = 3;
  numRowsValues[3] = 20;
  numRowsValues[4] = 100;
  numRowsValues[5] = 200;
  numRowsValues[6] = 500;

  const int numNumColsValues = 6;
  size_type numColsValues[6];
  numColsValues[0] = 0;
  numColsValues[1] = 1;
  numColsValues[2] = 2;
  numColsValues[3] = 5;
  numColsValues[4] = 11;
  numColsValues[5] = 20;

  const int numTransValues = 3;
  char transValues[3];
  transValues[0] = 'N';
  transValues[1] = 'T';
  transValues[2] = 'C';

  for (int alphaCoeff = 0; alphaCoeff < numCoeffValues; ++alphaCoeff) {
    const Scalar alpha = coeffValues[alphaCoeff];

    for (int betaCoeff = 0; betaCoeff < numCoeffValues; ++betaCoeff) {
      const Scalar beta = coeffValues[betaCoeff];

      for (int numRowsInd = 0; numRowsInd < numNumRowsValues; ++numRowsInd) {
        const size_type numRows = numRowsValues[numRowsInd];

        for (int numColsInd = 0; numColsInd < numNumColsValues; ++numColsInd) {
          const size_type numCols = numColsValues[numColsInd];

          for (int transInd = 0; transInd < numTransValues; ++transInd) {
            const char trans = transValues[transInd];

            curSuccess =
              testGemvOne<Scalar, Layout, Device> (out, indent, numRows,
                                                   numCols, trans, alpha, beta,
                                                   printOnlyOnFailure,
                                                   debug);
            success = success && curSuccess;
          }
        }
      }
    }
  }
  return success;
}

//
// Test all BLAS 2 operations, for a single combination of (Scalar,
// Layout, Device) template parameters.
//
template<class Scalar, class Layout, class Device>
bool
testMV (std::ostream& out,
        std::string indent,
        const bool printOnlyOnFailure,
        const bool debug)
{
  bool success = true;
  bool curSuccess = true;
  indent = indent + "  ";

  out << indent << "Test KokkosBlas::Impl::tryMklGemv" << endl;
  curSuccess =
    testGemvAll<Scalar, Layout, Device> (out, indent, printOnlyOnFailure, debug);
  success = success && curSuccess;

  return success;
}

//
// Test all BLAS 2 operations, for all Layouts in {LayoutLeft,
// LayoutRight}, and a single combination of (Scalar, Device) template
// parameters.
//
template<class Scalar, class Device>
bool
testLayouts (std::ostream& out,
             std::string indent,
             const bool testLayoutRight,
             const bool printOnlyOnFailure,
             const bool debug)
{
  using Kokkos::LayoutLeft;
  using Kokkos::LayoutRight;
  bool curSuccess = true;
  bool success = true;

  indent = indent + "  ";

  out << indent << "Testing Layout = LayoutLeft" << endl;
  curSuccess =
    testMV<Scalar, LayoutLeft, Device> (out, indent, printOnlyOnFailure, debug);
  success = success && curSuccess;

  if (testLayoutRight) {
    out << indent << "Testing Layout = LayoutRight" << endl;
    curSuccess =
      testMV<Scalar, LayoutRight, Device> (out, indent, printOnlyOnFailure, debug);
    success = success && curSuccess;
  }

  return success;
}

//
// Test all BLAS 2 operations, for all Scalar types in {double, float,
// int, complex<double>, complex<float>}, for all Layout types in
// {LayoutLeft, LayoutRight}, and a single Device template parameter.
//
template<class Device>
bool
testScalarsLayouts (std::ostream& out,
                    std::string indent,
                    const bool testFloat,
                    const bool testComplex,
                    const bool testLayoutRight,
                    const bool printOnlyOnFailure,
                    const bool debug)
{
  using std::endl;
  typedef typename Device::execution_space execution_space;
  bool success = true;
  indent = indent + "  ";

  if (execution_space::is_initialized ()) {
    bool curSuccess = true;

    out << indent << "Testing Scalar = double" << endl;
    curSuccess =
      testLayouts<double, Device> (out,
                                   indent,
                                   testLayoutRight,
                                   printOnlyOnFailure,
                                   debug);
    success = success && curSuccess;

    if (testFloat) {
      out << indent << "Testing Scalar = float" << endl;
      curSuccess =
        testLayouts<float, Device> (out,
                                    indent,
                                    testLayoutRight,
                                    printOnlyOnFailure,
                                    debug);
      success = success && curSuccess;
    }

    if (testComplex) {
      if (testFloat) { // this includes both real and complex
        out << indent << "Testing Scalar = Kokkos::complex<float>" << endl;
        curSuccess =
          testLayouts<Kokkos::complex<float>, Device> (out,
                                                       indent,
                                                       testLayoutRight,
                                                       printOnlyOnFailure,
                                                       debug);
        success = success && curSuccess;
      }

      out << indent << "Testing Scalar = Kokkos::complex<double>" << endl;
      curSuccess =
        testLayouts<Kokkos::complex<double>, Device> (out,
                                                      indent,
                                                      testLayoutRight,
                                                      printOnlyOnFailure,
                                                      debug);
      success = success && curSuccess;
    }
  }
  else {
    out << "Device NOT initialized; skipping test" << endl;
  }

  return success;
}

//
// Test all BLAS 2 operations, for all Scalar types in {double, float,
// int, complex<double>, complex<float>}, for all Layout types in
// {LayoutLeft, LayoutRight}, and for only the Serial and OpenMP
// execution spaces.  Don't test Threads, because those threads will
// interfere with the OpenMP threads used by the MKL.
//
bool
testScalarsLayoutsDevices (std::ostream& out,
                           std::string indent,
                           const bool testFloat,
                           const bool testComplex,
                           const bool testLayoutRight,
                           const bool printOnlyOnFailure,
                           const bool debug)
{
  bool success = true;
  bool curSuccess = true;

  indent = indent + "  ";
  out << indent << "Test all BLAS 2 operations" << endl;

#ifdef KOKKOS_HAVE_SERIAL
  {
    out << indent << "Test Serial, HostSpace" << endl;
    typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> device_type;
    curSuccess = testScalarsLayouts<device_type> (out, indent,
                                                  testFloat,
                                                  testComplex,
                                                  testLayoutRight,
                                                  printOnlyOnFailure,
                                                  debug);
    success = success && curSuccess;
  }
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
  {
    out << indent << "Test OpenMP, HostSpace" << endl;
    typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> device_type;
    curSuccess = testScalarsLayouts<device_type> (out, indent,
                                                  testFloat,
                                                  testComplex,
                                                  testLayoutRight,
                                                  printOnlyOnFailure,
                                                  debug);
    success = success && curSuccess;
  }
#endif // KOKKOS_HAVE_OPENMP

  return success;
}

} // namespace (anonymous)

#endif // TPETRAKERNELS_BLAS2_MV_GEMV_MKL_HPP
