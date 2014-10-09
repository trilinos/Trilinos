//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

/// \file Tsqr_CombineTest.hpp
/// \brief Test accuracy of TSQR::Combine
///
/// TSQR::Combine implements two structured QR factorizations: [R1;
/// R2] (where R1 and R2 are both square and upper triangular) and [R;
/// A] (where R is square and upper triangular, and A is a general
/// dense matrix).  This file implements accuracy tests for both
/// versions.
#ifndef __TSQR_Test_CombineTest_hpp
#define __TSQR_Test_CombineTest_hpp

#include <Tsqr_ConfigDefs.hpp>


namespace TSQR {
  namespace Test {

    /// \brief Test accuracy of TSQR::Combine.
    ///
    /// Test the accuracy of TSQR::Combine, and print the results to
    /// stdout.  Rather than template on Ordinal and Scalar, as with
    /// the full TSQR tests, we pick Ordinal=int and try four
    /// different Scalar types (the same four that LAPACK supports:
    /// double, float, complex<double>, complex<float> -- i.e.,
    /// S,D,C,Z).
    /// 
    /// \param numRows [in] When testing the [R; A] routines, the
    ///   number of rows in the cache block A.
    ///
    /// \param numCols [in] Number of columns in the test matrices.
    ///
    /// \param testReal [in] Whether or not to test TSQR::Combine for
    ///   real arithmetic.  For now, this means Scalar={float,double}.
    ///
    /// \param testComplex [in] Whether or not to test TSQR::Combine
    ///   for complex arithmetic.  For now, this means
    ///   Scalar={std::complex<float>, std::complex<double>}.
    ///
    /// \param printFieldNames [in] Whether to print field names (to
    ///   make machine parsing of results easier).
    ///
    /// \param simulateSequentialTsqr [in] Whether to use 
    ///   TSQR::Combine to simulate SequentialTsqr.
    ///
    /// \param debug [in] Whether to print (possibly verbose)
    ///   debugging output to stderr.
    ///
    /// \note This routine will only test Combine if at least one of
    ///   testReal and testComplex is true.
    void
    verifyCombine (const int numRows,
		   const int numCols, 
		   const bool testReal,
		   const bool testComplex,
		   const bool printFieldNames,
		   const bool simulateSequentialTsqr,
		   const bool debug);

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_CombineTest_hpp
