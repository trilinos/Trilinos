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

#include <MatrixMarket_Raw_Checker.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <algorithm>

using std::endl;

namespace {
  // Sample Matrix Market sparse matrix file.  We include this so we
  // can test without needing to read in a file.  Notice that all the
  // decimal floating-point values in this example can be represented
  // exactly in binary floating point.  This example has correct
  // syntax, so you won't need to use tolerant mode to parse it.
  const char sampleMatrixMarketFile[] =
    "%%MatrixMarket matrix coordinate real general\n"
    "5 5 10\n"
    "5 5 55.0\n"
    "4 4 44.0\n"
    "3 3 33.0\n"
    "2 2 22.0\n"
    "1 1 11.0\n"
    "4 5 45.0\n"
    "3 4 34.0\n"
    "2 3 23.0\n"
    "1 2 12.0\n"
    "1 5 15.0\n";

  // Given the three arrays of a CSR data structure, along with the
  // numbers of rows and columns, print the result to the given output
  // stream as a MatrixMarket file.
  template<class OrdinalType, class ScalarType>
  void
  csrToMatrixMarket (std::ostream& out,
                     Teuchos::ArrayView<OrdinalType> ptr,
                     Teuchos::ArrayView<OrdinalType> ind,
                     Teuchos::ArrayView<ScalarType> val,
                     const OrdinalType numRows,
                     const OrdinalType numCols)
  {
    using Teuchos::ArrayView;
    using std::endl;
    typedef ArrayView<OrdinalType>::size_type size_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "%%MatrixMarket matrix coordinate ";
    if (STS::isComplex) {
      out << "complex ";
    }
    else {
      out << "real ";
    }
    out << "general" << endl;
    out << numRows << " " << numCols << " " << ptr[numRows] << endl;
    for (OrdinalType rowIndex = 0; rowIndex < numRows; ++rowIndex) {
      for (OrdinalType k = ptr[rowIndex]; k < ptr[rowIndex+1]; ++k) {
        // Matrix Market files use 1-based row and column indices.
        out << (rowIndex+1) << " " << (ind[k]+1) << " ";
        if (STS::isComplex) {
          out << STS::real (val[k]) << " " << STS::imag (val[k]);
        }
        else {
          out << val[k];
        }
        out << endl;
      }
    }
  }


  // Return an RCP to a Kokkos Node instance.
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode() {
    throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
  }

  template<>
  Teuchos::RCP<Kokkos::SerialNode>
  getNode() {
    static Teuchos::RCP<Kokkos::SerialNode> node;
    if (node.is_null ()) {
      Teuchos::ParameterList defaultParams;
      node = Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
    }
    return node;
  }

#if defined(HAVE_KOKKOS_TBB)
  template<>
  Teuchos::RCP<Kokkos::TBBNode>
  getNode() {
    static Teuchos::RCP<Kokkos::TBBNode> node;
    if (node.is_null ()) {
      // "Num Threads" specifies the number of threads.  Defaults to an
      // automatically chosen value.
      Teuchos::ParameterList defaultParams;
      node = Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
    }
    return node;
  }
#endif // defined(HAVE_KOKKOS_TBB)

#if defined(HAVE_KOKKOS_TPL)
  template<>
  Teuchos::RCP<Kokkos::TPLNode>
  getNode() {
    static Teuchos::RCP<Kokkos::TBBNode> node;
    if (node.is_null ()) {
      // "Num Threads" specifies the number of threads.  Defaults to an
      // automatically chosen value.
      Teuchos::ParameterList defaultParams;
      node = Teuchos::rcp (new Kokkos::TPLNode (defaultParams));
    }
    return node;
  }
#endif // defined(HAVE_KOKKOS_TPL)
} // namespace (anonymous)

// Benchmark driver
int
main (int argc, char *argv[])
{
  using Teuchos::MatrixMarket::Raw::Checker;
  using Teuchos::MatrixMarket::Raw::Reader;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;
  typedef double scalar_type;
  typedef int ordinal_type;

  // Name of the Matrix Market sparse matrix file to read.
  std::string filename;
  // If true, just check the sparse matrix file.  Otherwise,
  // do a full conversion to CSR (compressed sparse row) format.
  bool checkOnly = true;
  // Whether to echo the sparse matrix to stdout after reading it
  // successfully.
  bool echo = false;
  // Whether to parse the Matrix Market file tolerantly.
  bool tolerant = false;
  // Verbosity of output
  bool verbose = false;
  // Whether to print debugging-level output
  bool debug = false;

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("filename", &filename,
                  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("checkOnly", "fullTest", &checkOnly,
                  "If true, just check the syntax of the input file.  "
                  "Otherwise, do a full test.");
  cmdp.setOption ("echo", "noecho", &echo,
                  "Whether to echo the sparse matrix contents to stdout "
                  "after reading it successfully.");
  cmdp.setOption ("tolerant", "strict", &tolerant,
                  "Whether to tolerate syntax errors in the Matrix Market file.");
  cmdp.setOption ("verbose", "quiet", &verbose,
                  "Print status output to stdout.");
  cmdp.setOption ("debug", "nodebug", &debug,
                  "Print possibly copious debugging output to stderr.");
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      std::cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
       parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
       std::invalid_argument, "Failed to parse command-line arguments.");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, the test passes trivially.
  bool success = true;
    if (checkOnly) {
      typedef Checker<scalar_type, ordinal_type> checker_type;
      checker_type checker (echo, tolerant, debug);

      RCP<const Comm<int> > comm = rcp (new SerialComm);
      if (filename != "") {
        success = success && checker.readFile (*comm, filename);
      }
      else {
        std::istringstream inStr (sampleMatrixMarketFile);
        success = success && checker.read (*comm, inStr);
      }
    }
    else {
      typedef Reader<scalar_type, ordinal_type> reader_type;
      reader_type reader (tolerant, debug);
      ArrayRCP<ordinal_type> ptr, ind;
      ArrayRCP<scalar_type> val;
      ordinal_type numRows, numCols;
      if (filename != "") {
        success = success && reader.readFile (ptr, ind, val, numRows, numCols, filename);
      }
      else {
        std::istringstream inStr (sampleMatrixMarketFile);
        success = success && reader.read (ptr, ind, val, numRows, numCols, inStr);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error, "Matrix Market "
        "reader failed to read the given file or input stream.");

      // Here's the fun part.  Output the CSR data to an output
      // stream.  Then read in the output stream.  The resulting
      // matrix should be exactly the same (unless the original file
      // had elements at the same location that were added together
      // with rounding error).
      std::ostringstream outStr;
      csrToMatrixMarket (outStr, ptr, ind, val, numRows, numCols);

      std::istringstream inStr (outStr.str ());
      ArrayRCP<ordinal_type> newptr, newind;
      ArrayRCP<scalar_type> newval;
      ordinal_type newNumRows, newNumCols;
      success = success && reader.read (newptr, newind, newval, newNumRows, newNumCols, inStr);
    }
  }

  if (success) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}



