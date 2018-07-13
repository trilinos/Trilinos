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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {

  using std::endl;
  using std::string;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceDomainMapAndImporter, LO, GO, Scalar, Node )
  {
    // Based on the FullTriDiag tests...

    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef ScalarTraits<MT> STM;

    const size_t ONE  = OrdinalTraits<size_t>::one();
    const size_t ZERO = OrdinalTraits<GO>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,ONE,comm);

    // RCP<FancyOStream> fos = Teuchos::fancyOStream(rcp(&std::cout,false));

    /* Create the following matrix:
    0  [2 1       ]   [2 1]
    1  [1 4 1     ]   [1 2] + [2 1]
    2  [  1 4 1   ]           [1 2] +
    3  [    1     ] =
       [       4 1]
   n-1 [       1 2]
    */

    MAT A(map,4);
    MAT B(map,4);
    A.setObjectLabel("The Matrix");
    A.setObjectLabel("The Other Matrix");
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*STS::one(),STS::one(),static_cast<Scalar>(2)*STS::one()));
      Array<GO> cols(tuple<GO>(myImageID,myImageID + 1));
      A.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      A.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
      B.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      B.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
    }
    A.fillComplete();
    B.fillComplete();

    // Build a one-process Map.
    if (comm->getSize() > 1) {
      MT norm = STM::zero ();

      // we know the map is contiguous...
      const size_t NumMyElements = (comm->getRank () == 0) ?
        A.getDomainMap ()->getGlobalNumElements () : 0;
      RCP<const Map<LO,GO,Node> > NewMap =
        rcp (new Map<LO,GO,Node> (INVALID, NumMyElements, ZERO, comm));
      RCP<const Tpetra::Import<LO,GO,Node> > NewImport =
        rcp (new Import<LO,GO,Node> (NewMap, A.getColMap ()));

      B.replaceDomainMapAndImporter (NewMap, NewImport);

      // Fill a random vector on the original map
      Vector<Scalar,LO,GO,Node> AVecX(A.getDomainMap());
      AVecX.randomize();

      // Import this vector to the new domainmap
      Vector<Scalar,LO,GO,Node> BVecX(B.getDomainMap());
      Tpetra::Import<LO,GO,Node> TempImport(A.getDomainMap(),NewMap); // (source,target)
      BVecX.doImport(AVecX,TempImport,Tpetra::ADD);

      // Now do some multiplies
      Vector<Scalar,LO,GO,Node> AVecY(A.getRangeMap());
      Vector<Scalar,LO,GO,Node> BVecY(B.getRangeMap());
      A.apply(AVecX,AVecY);
      B.apply(BVecX,BVecY);

      BVecY.update (-STS::one (), AVecY, STS::one ());
      norm = BVecY.norm2();

      out << "Residual 2-norm: " << norm << endl
          << "Residual 1-norm: " << BVecY.norm1 () << endl
          << "Residual Inf-norm: " << BVecY.normInf () << endl;

      // Macros don't like spaces, so we put the test outside the
      // macro.  Use <= rather than <, so the test passes even if
      // Scalar is an integer type.
      //
      // FIXME (mfh 10 Mar 2013) We should pick the tolerance relative
      // to the Scalar type and problem size.
      const bool normSmallEnough = norm <= as<MT> (1e-10);
      TEST_EQUALITY ( normSmallEnough, true );
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceDomainMapAndImporter, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}


