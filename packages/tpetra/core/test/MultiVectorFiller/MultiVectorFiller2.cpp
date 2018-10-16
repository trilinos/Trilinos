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

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_MultiVectorFiller.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <sstream>

namespace { // (anonymous)
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::CommRequest;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::wait;
  using Teuchos::waitAll;
  using std::endl;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::MultiVector<> MV;
  typedef Tpetra::global_size_t GST;
  typedef MV::scalar_type ST;
  typedef MV::local_ordinal_type LO;
  typedef MV::global_ordinal_type GO;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef MV::mag_type MT;
  typedef MV::execution_space execution_space;


  void
  syncErrs (bool& success,
            std::ostream& out,
            std::ostringstream& errStrm,
            int& lclExPass,
            int& gblSuccess,
            const Teuchos::Comm<int>& comm)
  {
    if (! success) {
      lclExPass = 0;
    }
    reduceAll<int, int> (comm, REDUCE_MIN, lclExPass, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      const int myRank = comm.getRank ();
      const int numProcs = comm.getSize ();

      ArrayRCP<size_t> sizeBuf (1);
      const std::string errStrmStr = errStrm.str ();
      const size_t errStrmStrSize = errStrmStr.size ();
      sizeBuf[0] = (lclExPass == 1) ? size_t (0) : errStrmStrSize;
      ArrayRCP<char> msgBuf;
      Array<RCP<CommRequest<int> > > requests (2);
      const int sizeTag = 42;
      const int msgTag = 43;

      if (myRank == 0 && lclExPass != 1) {
        out << errStrm.str ();
      }
      for (int p = 1; p < numProcs; ++p) {
        if (myRank == 0) {
          requests[0] = ireceive<int, size_t> (sizeBuf, p, sizeTag, comm);
          (void) wait<int> (comm, outArg (requests[0]));
          if (sizeBuf[0] != 0) {
            msgBuf.resize (sizeBuf[0] + 1); // size doesn't count the '\0'
            requests[1] = ireceive<int, char> (msgBuf, p, msgTag, comm);
            (void) wait<int> (comm, outArg (requests[1]));
            out << msgBuf.getRawPtr ();
          }
        }
        else if (myRank == p) {
          requests[0] = isend<int, size_t> (sizeBuf, 0, sizeTag, comm);
          (void) wait<int> (comm, outArg (requests[0]));
          if (lclExPass != 1) {
            msgBuf.resize (errStrmStrSize + 1);
            std::copy (errStrmStr.begin (), errStrmStr.end (), msgBuf.begin ());
            msgBuf[errStrmStrSize] = '\0';
            requests[1] = isend<int, char> (msgBuf, 0, sizeTag, comm);
            (void) wait<int> (comm, outArg (requests[1]));
          }
        }
      }

      if (myRank == 0) {
        out << "FAILED test on some process; exiting early" << endl;
      }
    }
  }

#define TPETRA_MVFILLER_SYNC_ERRS() do { \
  syncErrs (success, out, errStrm, lclExPass, gblSuccess, *comm); \
  if (gblSuccess != 1) { \
    return; \
  } \
} while (0)

#define TPETRA_MVFILLER_SYNC_OUT(msg) do { \
  if (debug) { \
    if (comm->getRank () == 0) { \
      std::cerr << msg << std::endl; \
    } \
  } else { \
    out << ">>> " << msg << std::endl; \
  } \
} while (0)


  // Test that MultiVectorFiller::clear actually clears out the object.
  TEUCHOS_UNIT_TEST( MultiVectorFiller, FillClearFill)
  {
    auto comm = Tpetra::getDefaultComm ();
    int gblSuccess = 1;
    //int lclSuccess = 1;
    int lclExPass = 1; // whether the calling process did NOT threw an exception
    const bool debug = true;

    Teuchos::OSTab tab0 (out);
    TPETRA_MVFILLER_SYNC_OUT("MultiVectorFiller double fill test");
    Teuchos::OSTab tab1 (out);

    TPETRA_MVFILLER_SYNC_OUT("Create nonoverlapping Map");

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const LO myNumRows = 5;
    auto nonoverlappingMap =
      rcp (new map_type (INV, static_cast<size_t> (myNumRows),
                         indexBase, comm));

    TPETRA_MVFILLER_SYNC_OUT("Create overlapping Map");

    // Make an overlapping Map, using the overlap pattern of a 1-D
    // Poisson equation (where the rows are distributed using the
    // above, nonoverlapping Map).  Put "remotes" at the end.
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (numProcs < 1 || myRank < 0 || myRank >= numProcs, std::logic_error,
       "Your Teuchos::Comm instance is totally messed up.  comm->getRank() = "
       << myRank << " and comm->getSize() = " << numProcs << ".");

    LO myNumInds = myNumRows;
    if (numProcs > 1) { // my process has at least one neighbor
      if (myRank == 0 || myRank == numProcs - 1) {
        myNumInds = myNumRows + 1; // one neighbor
      } else { // middle process; two neighbors
        myNumInds = myNumRows + 2;
      }
    }
    Teuchos::Array<GO> myInds (myNumInds);
    for (LO i_lcl = 0; i_lcl < myNumRows; ++i_lcl) {
      myInds[i_lcl] = nonoverlappingMap->getGlobalElement (i_lcl);
    }
    if (numProcs > 1) {
      if (myRank == 0) { // leftmost process
        myInds[myNumRows] = myInds[myNumRows-1] + static_cast<GO> (1);
      } else if (myRank == numProcs - 1) { // rightmost process
        myInds[myNumRows] = myInds[0] - static_cast<GO> (1);
      }
      else { // middle process; two neighbors (left and right)
        myInds[myNumRows] = myInds[0] - static_cast<GO> (1);
        myInds[myNumRows+1] = myInds[myNumRows-1] + static_cast<GO> (1);
      }
    }
    // Make the overlapping Map.
    auto overlappingMap = rcp (new map_type (INV, myInds (), indexBase, comm));

    // Number of columns in all the MultiVectors here.
    const size_t numCols = 3;

    TPETRA_MVFILLER_SYNC_OUT("First fill pass");

    // Fill all the entries of the overlapping Filler with ones.
    // While doing so, fill the entries of the overlapping MultiVector
    // with ones.  The two should represent the same data with the
    // same distribution.  We'll check this below.
    Tpetra::MultiVectorFiller<MV> filler (nonoverlappingMap, numCols);
    MV X_overlap (overlappingMap, numCols);
    X_overlap.sync<MV::dual_view_type::t_host::memory_space> (); // sync to host
    X_overlap.modify<MV::dual_view_type::t_host::memory_space> (); // will modify on host

    std::ostringstream errStrm;
    try {
      if (overlappingMap->getNodeNumElements () != 0) {
        auto X_overlap_lcl = X_overlap.getLocalView<MV::dual_view_type::t_host::memory_space> ();
        Teuchos::Array<GO> rows (1);
        Teuchos::Array<ST> vals (1);
        for (GO i_gbl = overlappingMap->getMinGlobalIndex ();
             i_gbl <= overlappingMap->getMaxGlobalIndex (); ++i_gbl) {
          const LO i_lcl = overlappingMap->getLocalElement (i_gbl);
          rows[0] = i_gbl;
          vals[0] = STS::one ();
          for (size_t col = 0; col < numCols; ++col) {
            filler.sumIntoGlobalValues (rows (), col, vals ());
            X_overlap_lcl(i_lcl, col) = vals[0];
          }
        }
        X_overlap.sync<MV::dual_view_type::t_dev::memory_space> (); // sync to device
      }
    } catch (std::exception& e) {
      errStrm << "Proc " << myRank << " threw an exception while filling: "
              << e.what () << endl;
      lclExPass = 0;
    }
    TPETRA_MVFILLER_SYNC_ERRS();

    TPETRA_MVFILLER_SYNC_OUT("Done with first fill pass");

    // Make a MultiVector using the nonoverlapping Map, and assemble
    // into it, using the Filler.
    MV X_result (nonoverlappingMap, numCols);

    try {
      //filler.globalAssemble (X_result, false, debug);
      filler.globalAssemble (X_result);
    } catch (std::exception& e) {
      errStrm << "Proc " << myRank << " threw an exception during "
        "filler.globalAssemble: " << e.what () << endl;
      lclExPass = 0;
    }
    TPETRA_MVFILLER_SYNC_ERRS();

    TPETRA_MVFILLER_SYNC_OUT("Repeat first fill pass with Export");

    // Do the same thing, but use Export instead.
    Tpetra::Export<> exp (overlappingMap, nonoverlappingMap);
    MV X_result2 (nonoverlappingMap, numCols);
    X_result2.doExport (X_overlap, exp, Tpetra::ADD);

    // Compare the two MultiVectors using the 1-norm of their
    // difference.  The 1-norm helps to show how _many_ entries are
    // off.
    X_result2.update (-STS::one (), X_result, STS::one ());
    Kokkos::View<MT*, execution_space> norms ("norms", numCols);
    X_result2.norm1 (norms);
    for (size_t col = 0; col < numCols; ++col) {
      TEST_EQUALITY_CONST( norms[col], STS::zero () );
    }

    //
    // Test that MultiVectorFiller::clear() actually clears out the object.
    //
    TPETRA_MVFILLER_SYNC_OUT("Call MultiVectorFiller::clear");

    try {
      filler.clear ();
    } catch (std::exception& e) {
      errStrm << "Proc " << myRank << " threw an exception during "
        "filler.clear(): " << e.what () << endl;
      lclExPass = 0;
    }
    TPETRA_MVFILLER_SYNC_ERRS();

    TPETRA_MVFILLER_SYNC_OUT("Second fill pass");

    try {
      X_overlap.putScalar (STS::zero ());
      X_overlap.sync<MV::dual_view_type::t_host::memory_space> (); // sync to host
      X_overlap.modify<MV::dual_view_type::t_host::memory_space> (); // will modify on host
      if (overlappingMap->getNodeNumElements () != 0) {
        auto X_overlap_lcl = X_overlap.getLocalView<MV::dual_view_type::t_host::memory_space> ();
        Teuchos::Array<GO> rows (1);
        Teuchos::Array<ST> vals (1);
        for (GO i_gbl = overlappingMap->getMinGlobalIndex ();
             i_gbl <= overlappingMap->getMaxGlobalIndex (); ++i_gbl) {
          const LO i_lcl = overlappingMap->getLocalElement (i_gbl);
          rows[0] = i_gbl;
          vals[0] = STS::one ();
          for (size_t col = 0; col < numCols; ++col) {
            filler.sumIntoGlobalValues (rows (), col, vals ());
            X_overlap_lcl(i_lcl, col) = vals[0];
          }
        }
        X_overlap.sync<MV::dual_view_type::t_dev::memory_space> (); // sync to device
      }
    } catch (std::exception& e) {
      errStrm << "Proc " << myRank << " threw an exception while filling "
        "(2nd pass): " << e.what () << endl;
      lclExPass = 0;
    }
    TPETRA_MVFILLER_SYNC_ERRS();

    TPETRA_MVFILLER_SYNC_OUT("Done with second fill pass");

    // Make a MultiVector using the nonoverlapping Map, and assemble
    // into it, using the Filler.
    X_result.putScalar (STS::zero ());
    try {
      //filler.globalAssemble (X_result, false, debug);
      filler.globalAssemble (X_result);
    } catch (std::exception& e) {
      errStrm << "Proc " << myRank << " threw an exception during "
        "filler.globalAssemble (2nd pass): " << e.what () << endl;
      lclExPass = 0;
    }
    TPETRA_MVFILLER_SYNC_ERRS();

    TPETRA_MVFILLER_SYNC_OUT("Done with second globalAssemble");

    // Do the same thing, but use the Export instead.
    X_result2.putScalar (STS::zero ());
    X_result2.doExport (X_overlap, exp, Tpetra::ADD);

    // Compare the two MultiVectors, using the inf-norm of the difference.
    X_result2.update (-STS::one (), X_result, STS::one ());
    X_result2.normInf (norms);
    for (size_t col = 0; col < numCols; ++col) {
      TEST_EQUALITY_CONST( norms[col], STS::zero () );
    }
  }

} // namespace (anonymous)
