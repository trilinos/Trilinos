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

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_ETIHelperMacros.h>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include "Kokkos_DefaultNode.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

namespace {

  template<class MapType>
  class Test {
  public:
    typedef double ST;
    typedef Teuchos::ScalarTraits<double>::magnitudeType MT;

    typedef typename MapType::local_ordinal_type LO;
    typedef typename MapType::global_ordinal_type GO;
    typedef typename MapType::node_type NT;

    typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;
    typedef typename Teuchos::ArrayView<const GO>::size_type size_type;

    static Teuchos::RCP<const MapType>
    createTestMap (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   const Teuchos::RCP<NT>& node)
    {
      using Tpetra::createNonContigMapWithNode;
      using Teuchos::ArrayView;

      // Elements of the Map owned by Proc 0.  No other processes in
      // the communicator own anything.  Don't put this in an inner
      // scope, since it has to persist up to Map creation.
      const GO procZeroEltList[] = {0, 1, 2,
                                    12, 13, 14,
                                    6, 7, 8,
                                    18, 19, 20,
                                    3, 4, 5,
                                    15, 16, 17,
                                    9, 10, 11,
                                    21, 22, 23};
      const size_type numProcZeroElts = 24;

      ArrayView<const GO> eltList;
      if (comm->getRank () == 0) {
        eltList = ArrayView<const GO> (procZeroEltList, numProcZeroElts);
      }
      else  {
        eltList = ArrayView<const GO> (NULL, 0);
      }
      return createNonContigMapWithNode<LO, GO, NT> (eltList, comm, node);
    }

    static Teuchos::RCP<const MV>
    createTestVector (const Teuchos::RCP<const MapType>& map)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // We write in row major form, rather than the column major form
      // that MultiVector's constructor requires, so that it's easy for
      // humans to read.  This makes us convert to column major below.
      const ST rawDataRowMajor[] = {
        1.0000,         0,         0,    0.5000,         0,         0,
             0,    1.0000,         0,   -0.5000,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,    0.5000,         0,         0,
             0,    1.0000,         0,    0.5000,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0,
        1.0000,         0,         0,         0,         0,         0,
             0,    1.0000,         0,         0,         0,         0,
             0,         0,    1.0000,         0,         0,         0};

      const size_type numRows = 24;
      const size_type numCols = 6;
      const size_type numRawData = as<size_type> (numRows * numCols);

      // Convert from row major to column major.
      Array<ST> rawDataColMajor (numRawData);
      for (size_type j = 0; j < numCols; ++j) {
        for (size_type i = 0; i < numRows; ++i) {
          rawDataColMajor[i + j*numRows] = rawDataRowMajor[j + i*numCols];
        }
      }

      ArrayView<const ST> data;
      size_t LDA;
      if (map->getComm ()->getRank () == 0) {
        data = rawDataColMajor ();
        LDA = numRows;
      }
      else {
        data = ArrayView<const ST> (NULL, 0);
        LDA = 0;
      }
      return rcp (new MV (map, data, LDA, numCols));
    }

    static std::string
    writeMultiVectorToString (const MV& X) {
      using Teuchos::rcpFromRef;
      typedef Tpetra::CrsMatrix<ST, LO, GO, NT> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;

      std::ostringstream os;
      writer_type::writeDense (os, rcpFromRef (X));
      return os.str ();
    }

    static Teuchos::RCP<const MV>
    readMultiVectorFromString (const std::string& s,
                               Teuchos::RCP<const MapType> map)
    {
      typedef Tpetra::CrsMatrix<ST, LO, GO, NT> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

      std::istringstream in (s);
      // map needs to be a nonconst RCP to const MapType, because
      // readDense would set map to a computed Map on output if if
      // were null on input.
      return reader_type::readDense (in, map->getComm (), map->getNode (), map);
    }

    // Ensure that X and Y have the same dimensions, and that the
    // corresponding columns of X and Y have 2-norms that differ by no
    // more than a prespecified tolerance (that accounts for rounding
    // errors in computing the 2-norm).
    static void
    assertMultiVectorsEqual (const Teuchos::RCP<const MV>& X,
                             const Teuchos::RCP<const MV>& Y)
    {
      using Teuchos::Array;
      using Teuchos::as;
      using Teuchos::RCP;
      using std::cerr;
      using std::cout;
      using std::endl;
      typedef Teuchos::ScalarTraits<ST> STS;
      typedef Teuchos::ScalarTraits<MT> STM;

      TEUCHOS_TEST_FOR_EXCEPTION(X->getGlobalLength() != Y->getGlobalLength(),
                                 std::logic_error, "Y has a different number of rows than X.");
      TEUCHOS_TEST_FOR_EXCEPTION(X->getNumVectors() != Y->getNumVectors(),
                                 std::logic_error, "Y has a different number of columns than X.");

      Tpetra::global_size_t numRows = X->getGlobalLength();
      const size_t numVecs = X->getNumVectors();
      Array<MT> X_norm2 (numVecs);
      Array<MT> Y_norm2 (numVecs);
      X->norm2 (X_norm2 ());
      Y->norm2 (Y_norm2 ());

      // For the relative tolerance, I'm using the typical heuristic: a
      // (fudge factor) times (machine epsilon) times sqrt(the number of
      // floating-point numbers involved in computing the norm for one
      // column).  Our output routine is careful to use enough digits,
      // so the input matrix shouldn't be that much different.
      const MT tol = as<MT> (10) *
        STS::magnitude (STS::eps ()) *
        STM::squareroot (as<MT> (numRows));
      Array<size_t> badColumns;
      for (size_t j = 0; j < numVecs; ++j) {
        // If the norm of the current column of X is zero, use the
        // absolute difference; otherwise use the relative difference.
        if ((X_norm2[j] == STM::zero() && STS::magnitude (Y_norm2[j]) > tol) ||
            STS::magnitude (X_norm2[j] - Y_norm2[j]) > tol) {
          badColumns.push_back (j);
        }
      }

      if (badColumns.size() > 0) {
        const size_t numBad = badColumns.size();
        std::ostringstream os;
        std::copy (badColumns.begin(), badColumns.end(),
                   std::ostream_iterator<size_t> (os, " "));
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Column"
          << (numBad != 1 ? "s" : "") << " [" << os.str() << "] of X and Y have "
          "norms that differ relatively by more than " << tol << ".");
      }
    }

    static void
    testPermutedMultiVectorOutput (Teuchos::FancyOStream& out,
                                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                   const Teuchos::RCP<NT>& node)
    {
      using Teuchos::getFancyOStream;
      using Teuchos::FancyOStream;
      using Teuchos::RCP;
      using Teuchos::rcpFromRef;
      using std::endl;

      //RCP<FancyOStream> out = getFancyOStream (rcpFromRef (std::cerr));

      RCP<const MapType> map = createTestMap (comm, node);
      RCP<const MV> X_orig = createTestVector (map);

      // Test that writing out the multivector and reading it back in
      // doesn't change the multivector.
      const std::string outStr = writeMultiVectorToString (*X_orig);
      RCP<const MV> X_outIn = readMultiVectorFromString (outStr, map);

      // Before we finalize the test, print out the multivectors in
      // different ways, so we can debug any test failures.
      out << "Here is the original multivector:" << endl;
      X_orig->describe (out, Teuchos::VERB_EXTREME);
      out << "Here is the printout from the original multivector:" << endl
           << outStr << endl;
      out << "Here is the multivector read in from the printout:" << endl;
      X_outIn->describe (out, Teuchos::VERB_EXTREME);

      assertMultiVectorsEqual (X_orig, X_outIn);
    }
  };

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Tpetra_MatrixMarket, MultiVector_Output_Perm, LO, GO )
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef Kokkos::DefaultNode::DefaultNodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  // Get the default communicator.
  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();

  if (myRank == 0) {
    out << "Test with " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;
  }

  // This test doesn't make much sense if there are multiple MPI processes.
  // We let it pass trivially in that case.
  if (numProcs != 1) {
    out << "Number of processes in world is not one; test passes trivially." << endl;
    return;
  }

  // Get a Kokkos Node instance.  It would be nice if we could pass in
  // parameters here, but threads don't matter for this test; it's a
  // test for distributed-memory capabilities.
  if (myRank == 0) {
    out << "Creating Kokkos Node of type " << TypeNameTraits<NT>::name () << endl;
  }
  RCP<NT> node = Kokkos::DefaultNode::getDefaultNode();
  // Run the actual test.
  //RCP<const map_type> map = Test<map_type>::createTestMap (comm, node);
  Test<map_type>::testPermutedMultiVectorOutput (out, comm, node);
}

//////////////////////////////////////////////////////////////////////
// INSTANTIATE THE TEMPLATED UNIT TESTS
//////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Tpetra_MatrixMarket, MultiVector_Output_Perm, LO, GO )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP )

} // namespace (anonymous)
