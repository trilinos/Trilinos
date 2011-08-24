//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
