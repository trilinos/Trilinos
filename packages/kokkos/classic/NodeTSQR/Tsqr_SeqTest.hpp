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

#ifndef __TSQR_Test_SeqTest_hpp
#define __TSQR_Test_SeqTest_hpp

#include <Tsqr_ConfigDefs.hpp>

#include <cstring> // size_t definition
#include <string>
#include <iostream>


namespace TSQR {
  namespace Test {

    /// \brief Test accuracy of SequentialTsqr.
    ///
    /// Test the accuracy of our sequential TSQR implementation
    /// (SequentialTsqr), on an nrows by ncols matrix, using the given
    /// cache size hint (in bytes).  Print the results to the given
    /// output stream out.
    void
    verifySeqTsqr (std::ostream& out,
		   const int nrows, 
		   const int ncols, 
		   const size_t cache_size_hint,
		   const bool test_complex_arithmetic,
		   const bool save_matrices,
		   const bool contiguous_cache_blocks,
		   const std::string& additionalFieldNames,
		   const std::string& additionalData,
		   const bool printFieldNames,
		   const bool human_readable = false,
		   const bool b_debug = false);

    /// \brief Test accuracy of LAPACK's QR factorization.
    ///
    /// Test the accuracy of LAPACK's QR factorization (_GEQRF +
    /// _ORGQR) on an nrows by ncols matrix, and print the results to
    /// the given output stream out.
    void
    verifyLapack (std::ostream& out,
		  const int nrows, 
		  const int ncols, 
		  const bool test_complex_arithmetic,
		  const std::string& additionalFieldNames,
		  const std::string& additionalData,
		  const bool printFieldNames,
		  const bool human_readable,
		  const bool b_debug = false);

    /// \brief Test performance of SequentialTsqr.
    ///
    /// Test the run time over ntrials trials of sequential TSQR, on
    /// an nrows by ncols matrix (using the given cache block size (in
    /// bytes)), and print the results to the given output stream out.
    ///
    /// \param human_readable [in] If true, print the benchmark
    ///   results to stdout in human-readable format.  Otherwise,
    ///   print them as two rows of comma-delimited ASCII, in an
    ///   abbreviated format suitable for automatic processing.
    void
    benchmarkSeqTsqr (std::ostream& out,
		      const int numRows,
		      const int numCols,
		      const int numTrials,
		      const size_t cacheSizeHint,
		      const bool contiguousCacheBlocks,
		      const bool testComplex,
		      const std::string& additionalFieldNames,
		      const std::string& additionalData,
		      const bool printFieldNames,
		      const bool humanReadable);

    /// \brief Test performance of LAPACK's QR factorization.
    ///
    /// Test the run time over numTrials trials of LAPACK QR (_GEQRF +
    /// _ORGQR), on a numRows by numCols matrix, and print the results
    /// to the given output stream out.
    ///
    /// \param humanReadable [in] If true, print the benchmark results
    ///   to out in human-readable format.  Otherwise, print them as
    ///   two rows of comma-delimited ASCII, in an abbreviated format
    ///   suitable for automatic processing.
    void
    benchmarkLapack (std::ostream& out,
		     const int numRows,
		     const int numCols,
		     const int numTrials,
		     const bool testComplex,
		     const std::string& additionalFieldNames,
		     const std::string& additionalData,
		     const bool printFieldNames,
		     const bool humanReadable);

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_SeqTest_hpp
