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

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Tpetra_MultiVectorFiller.hpp"

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::cerr;
  using std::cout;
  using std::endl;

  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef Tpetra::MultiVector<>::local_ordinal_type LO;
#if defined (HAVE_TPETRA_INT_LONG_LONG)
  typedef long long GO;
#else
  typedef Tpetra::MultiVector<>::global_ordinal_type GO;
#endif
  typedef Tpetra::MultiVector<>::node_type NT;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();

  //
  // Read in command line arguments.
  //
  int unknownsPerNode = 20; // number of unknowns per process
  int unknownsPerElt = 3; // number of unknowns per (overlapping) element
  int numCols = 1;
  bool verbose = false;

  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("unknownsPerNode", &unknownsPerNode,
                  "Number of unknowns per process");
  cmdp.setOption ("unknownsPerElt", &unknownsPerElt,
                  "Number of unknowns per (overlapping) element.");
  cmdp.setOption ("numCols", &numCols,
                  "Number of columns in the multivector.  Must be positive.");
  cmdp.setOption ("verbose", "quiet", &verbose,
                  "Whether to print verbose output.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  const Teuchos::EVerbosityLevel verbLevel =
    verbose ? Teuchos::VERB_EXTREME : Teuchos::VERB_DEFAULT;
  RCP<FancyOStream> out = verbose ? getFancyOStream (rcpFromRef (std::cout)) :
    getFancyOStream (rcpFromRef (blackHole));

  RCP<ParameterList> nodeParams = parameterList ("Kokkos Node");

  // Run the test.
  bool succeeded = true;
  try {
    using Tpetra::Test::testMultiVectorFiller;
    testMultiVectorFiller<ST, LO, GO, NT> (comm, static_cast<size_t> (unknownsPerNode),
                                           static_cast<GO> (unknownsPerElt),
                                           static_cast<size_t> (numCols), out, verbLevel);
    succeeded = true;
  } catch (std::exception& e) {
    *out << "MultiVectorFiller test threw an exception:  " << e.what() << endl;
    succeeded = false;
  }

  const int localSuccess = succeeded ? 1 : 0;
  int globalSuccess = localSuccess;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  reduceAll<int, int> (*comm, REDUCE_MIN, localSuccess, outArg (globalSuccess));

  if (globalSuccess == 1) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


