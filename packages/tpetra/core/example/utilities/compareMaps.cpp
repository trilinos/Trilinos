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

#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>
#include <sstream>
#include <iostream>

// This program reads two Maps from a file.  The Maps must have the
// same process counts as the input communicator.
int
main (int argc, char* argv[])
{
  using Teuchos::broadcast;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::cout;
  using std::endl;
  typedef Tpetra::Map<> map_type;
  // Reader must be templated on a CrsMatrix specialization
  typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<> > reader_type;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int rootRank = 0;

    std::string mapFile1, mapFile2;
    bool tolerant = false;
    bool debug = false;
    {
      const bool throwExceptions = false;
      const bool recognizeAllOptions = true; // let Kokkos options through
      Teuchos::CommandLineProcessor cmdLineProc (throwExceptions, recognizeAllOptions);

      cmdLineProc.setOption ("mapFile1", &mapFile1, "Filename of first Map file, to read on Process 0");
      cmdLineProc.setOption ("mapFile2", &mapFile2, "Filename of second Map file, to read on Process 0");
      cmdLineProc.setOption ("tolerant", "strict", &tolerant, "Read Maps in tolerant mode");
      cmdLineProc.setOption ("debug", "release", &debug, "Enable debug output");

      auto result = cmdLineProc.parse (argc, argv);
      if (result != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
	if (result == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
	  return EXIT_SUCCESS; // printed help; now we're done
	}
	else if (result != Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION) {
	  return EXIT_FAILURE; // some error that's not just an unrecognized option
	}
      }
    }

    // On Process 0, test whether the files exist.
    int lclFileExists = 0;
    int gblFileExists = 0;
    if (myRank == rootRank) {
      std::ifstream f1 (mapFile1);
      lclFileExists = f1.good () ? 1 : 0;
      gblFileExists = lclFileExists;
      f1.close ();
    }
    broadcast<int, int> (*comm, rootRank, outArg (gblFileExists));
    if (gblFileExists != 1) {
      if (myRank == rootRank) {
	cerr << "Cannot read the first Map file \"" << mapFile1 << "\"!" << endl;
      }
      return EXIT_FAILURE;
    }

    if (myRank == rootRank) {
      std::ifstream f2 (mapFile2);
      lclFileExists = f2.good ();
      gblFileExists = lclFileExists;
      f2.close ();
    }
    broadcast<int, int> (*comm, rootRank, outArg (gblFileExists));
    if (gblFileExists != 1) {
      if (myRank == rootRank) {
	cerr << "Cannot read the second Map file \"" << mapFile2 << "\"!" << endl;
      }
      return EXIT_FAILURE;
    }

    if (debug && myRank == rootRank) {
      cerr << "Both Map files are readable" << endl;
    }

    RCP<const map_type> map1;
    RCP<const map_type> map2;

    int lclSuccess = 0;
    int gblSuccess = 0;
    std::ostringstream errStrm;
    try {
      map1 = reader_type::readMapFile (mapFile1, comm, tolerant, debug);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      errStrm << "Process " << myRank << ": Caught exception while reading "
	"first Map file: " << e.what () << endl;
    }
    catch (...) {
      errStrm << "Process " << myRank << ": Caught exception not a subclass of "
	"std::exception, while reading first Map file" << endl;
    }
    if (map1.is_null ()) {
      errStrm << "Process " << myRank << ": Attempt to read first Map file "
	"returned null" << endl;
      lclSuccess = 0;
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      Tpetra::Details::gathervPrint (std::cerr, errStrm.str (), *comm);
      return EXIT_FAILURE;
    }

    if (debug && myRank == rootRank) {
      cerr << "Read first Map file" << endl;
    }

    try {
      map2 = reader_type::readMapFile (mapFile2, comm, tolerant, debug);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      errStrm << "Process " << myRank << ": Caught exception while reading "
	"second Map file: " << e.what () << endl;
    }
    catch (...) {
      errStrm << "Process " << myRank << ": Caught exception not a subclass of "
	"std::exception, while reading second Map file" << endl;
    }
    if (map2.is_null ()) {
      errStrm << "Process " << myRank << ": Attempt to read second Map file "
	"returned null" << endl;
      lclSuccess = 0;
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      Tpetra::Details::gathervPrint (std::cerr, errStrm.str (), *comm);
      return EXIT_FAILURE;
    }

    if (debug && myRank == rootRank) {
      cerr << "Read second Map file" << endl;
    }

    // At this point, we have two valid Maps.  Let's print stuff about them.
    const bool compat = map1->isCompatible (*map2);
    const bool same = map1->isSameAs (*map2);
    // This is a local (per MPI process) property.  Thus, we have to
    // all-reduce to find out whether it's true on all processes
    // (which is actually what we want to know).
    const bool locallyFitted = map1->isLocallyFitted (*map2);
    const int locallyFittedInt = locallyFitted ? 1 : 0;
    int globallyFittedInt = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, locallyFittedInt,
			 outArg (globallyFittedInt));
    const bool globallyFitted = (globallyFittedInt == 1);

    if (myRank == rootRank) {
      cout << "Maps compatible? " << (compat ? "YES" : "NO") << endl
	   << "Maps same?       " << (same ? "YES" : "NO") << endl
	   << "Maps fitted?     " << (globallyFitted ? "YES" : "NO") << endl;
    }
  }
  return EXIT_SUCCESS;
}
