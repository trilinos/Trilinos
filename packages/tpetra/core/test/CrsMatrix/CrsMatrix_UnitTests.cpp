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
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_getNumDiags.hpp"
#include "Tpetra_Details_residual.hpp"

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Tpetra::createContigMapWithNode;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::NO_TRANS;
  //using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using std::endl;
  typedef Tpetra::global_size_t GST;

#define STD_TESTS(matrix) \
  { \
    using Teuchos::outArg; \
    RCP<const Comm<int> > STCOMM = matrix.getComm(); \
    ArrayView<const GO> STMYGIDS = matrix.getRowMap()->getNodeElementList(); \
    ArrayView<const LO> loview; \
    ArrayView<const Scalar> sview; \
    size_t STMAX = 0; \
    for (size_t STR=0; STR < matrix.getNodeNumRows(); ++STR) { \
      const size_t numEntries = matrix.getNumEntriesInLocalRow(STR); \
      TEST_EQUALITY( numEntries, matrix.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
      matrix.getLocalRowView(STR,loview,sview); \
      TEST_EQUALITY( static_cast<size_t>(loview.size()), numEntries ); \
      TEST_EQUALITY( static_cast<size_t>( sview.size()), numEntries ); \
      STMAX = std::max( STMAX, numEntries ); \
    } \
    TEST_EQUALITY( matrix.getNodeMaxNumRowEntries(), STMAX ); \
    GST STGMAX; \
    Teuchos::reduceAll<int,GST>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
    TEST_EQUALITY( matrix.getGlobalMaxNumRowEntries(), STGMAX ); \
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, BadCalls, LO, GO, Scalar, Node )
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef RCP<const Tpetra::Map<LO,GO,Node> > RCPMap;
    typedef Teuchos::ScalarTraits<Mag> MT;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    // create the zero matrix
    RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > zero;
    {
      RCPMap map  = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
      MV mv(map,1);
      zero = rcp( new MAT(map,0,TPETRA_DEFAULT_PROFILE_TYPE) );
      TEST_THROW(zero->apply(mv,mv), std::runtime_error);
#   if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS)
      // throw exception because we required increased allocation
      TEST_THROW(zero->insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(0),tuple<Scalar>(ST::one())), std::runtime_error);
#   endif
      TEST_ASSERT( zero->getProfileType() == TPETRA_DEFAULT_PROFILE_TYPE );
      zero->fillComplete();
    }
    STD_TESTS((*zero));
    TEST_EQUALITY_CONST( zero->getRangeMap() == zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( zero->getFrobeniusNorm(), MT::zero() );
    const RCPMap drmap = zero->getDomainMap();
    {
      MV mv1(drmap,1), mv2(drmap,2), mv3(drmap,3);
      TEST_THROW(zero->apply(mv2,mv1), std::runtime_error); // MVs have different number of vectors
      TEST_THROW(zero->apply(mv2,mv3), std::runtime_error); // MVs have different number of vectors
    }
    // test that our assumptions on the maps are correct:
    // that is, that badmap is not equal to the range, domain, row or colum map of the matrix
    const RCPMap badmap = createContigMapWithNode<LO,GO,Node>(INVALID,1,comm);
    TEST_EQUALITY_CONST( badmap != zero->getRowMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getColMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getRangeMap(),  true );
    TEST_EQUALITY_CONST( *badmap != *zero->getRowMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getColMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getRangeMap(),  true );
    // now test the multivector against the matrix operators
    // Bugzilla bug #5247
    {
      MV mvbad(badmap,1);
#ifdef HAVE_TPETRA_DEBUG
      const Scalar ONE = ST::one(), ZERO = ST::zero();
      MV mvcol(zero->getColMap(),1);
      MV mvrow(zero->getRowMap(),1);
      TEST_THROW(zero->localApply(mvcol,mvbad,  NO_TRANS,ONE,ZERO), std::runtime_error); // bad output map
      TEST_THROW(zero->localApply(mvbad,mvrow,  NO_TRANS,ONE,ZERO), std::runtime_error); // bad input map
      TEST_THROW(zero->localApply(mvbad,mvcol,CONJ_TRANS,ONE,ZERO), std::runtime_error); // bad output map
      TEST_THROW(zero->localApply(mvrow,mvbad,CONJ_TRANS,ONE,ZERO), std::runtime_error); // bad input map
#endif
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TheEyeOfTruth, LO, GO, Scalar, Node )
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Tpetra::Map<LO,GO,Node> > map =
      createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    // create the identity matrix
    GO base = numLocal*myImageID;
    RCP<Tpetra::RowMatrix<Scalar,LO,GO,Node> > eye;
    {
      RCP<MAT> eye_crs = rcp(new MAT(map,numLocal,TPETRA_DEFAULT_PROFILE_TYPE));
      for (size_t i=0; i<numLocal; ++i) {
        eye_crs->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<Scalar>(ST::one()));
      }
      TEST_ASSERT( eye_crs->getProfileType() == TPETRA_DEFAULT_PROFILE_TYPE );
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (*eye), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (*eye), static_cast<LO> (numLocal) );
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);
    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, SimpleEigTest, LO, GO, Scalar, Node )
  {
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    const size_t ONE = Teuchos::OrdinalTraits<size_t>::one();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    RCP<const Tpetra::Map<LO,GO,Node> > map =
      createContigMapWithNode<LO,GO,Node>(INVALID,ONE,comm);
    // create a multivector ones(n,1)
    MV ones(map,ONE,false), threes(map,ONE,false);
    ones.putScalar(ST::one());
    /* create the following matrix:
       [2 1           ]
       [1 1 1         ]
       [  1 1 1       ]
       [   . . .      ]
       [     . . .    ]
       [       . . .  ]
       [         1 1 1]
       [           1 2]
     this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
    */
    size_t myNNZ;
    MAT A(map,3);
    if (myImageID == 0) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*ST::one(), ST::one()));
      Array<GO> cols(tuple<GO>(myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(ST::one(), static_cast<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID-1,myImageID));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else {
      myNNZ = 3;
      Array<Scalar> vals(3,ST::one());
      Array<GO> cols(tuple<GO>(myImageID-1, myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.getGlobalNumEntries()   , static_cast<size_t>(3*numImages-2));
    TEST_EQUALITY(A.getNodeNumEntries()       , myNNZ);
    TEST_EQUALITY(A.getGlobalNumRows()       , static_cast<size_t>(numImages));
    TEST_EQUALITY_CONST(A.getNodeNumRows()     , ONE);
    TEST_EQUALITY(A.getNodeNumCols()           , myNNZ);
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (A), static_cast<GO> (numImages));
    TEST_EQUALITY_CONST( Tpetra::Details::getLocalNumDiags (A), static_cast<LO> (ONE) );
    TEST_EQUALITY(A.getGlobalMaxNumRowEntries() , (numImages > 2 ? 3 : 2));
    TEST_EQUALITY(A.getNodeMaxNumRowEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getDomainMap()), true);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getRangeMap()) , true);
    // test the action
    threes.randomize();
    A.apply(ones,threes);
    // now, threes should be 3*ones
    threes.update(static_cast<Scalar>(-3)*ST::one(),ones,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    threes.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ZeroMatrix, LO, GO, Scalar, Node )
  {
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Tpetra::Map<LO,GO,Node> > map =
      createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    // create the zero matrix
    MAT zero(map,0);
    zero.fillComplete();
    //
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    mvres.putScalar(1);
    zero.apply(mvrand,mvres);
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ImbalancedRowMatrix, LO, GO, Scalar, Node )
  {
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocalRows = 10;
    // compute the number of entries in the long rows required to trigger the imbalanced row apply
    // this isn't quite the minimum but it works: numLocalRows * 5 + Behavior::rowImbalanceThreshold()
    const size_t numLocalColumns = 1 + (5 + 1.5 * Tpetra::Details::Behavior::rowImbalanceThreshold()) / (1.0 - 1.0 / numLocalRows);
    const size_t numVecs = 2;
    const int numRanks = comm->getSize();
    RCP<const Tpetra::Map<LO,GO,Node> > rowMap =
      createContigMapWithNode<LO,GO,Node>(numLocalRows * numRanks, numLocalRows, comm);
    RCP<const Tpetra::Map<LO,GO,Node> > colMap =
      createContigMapWithNode<LO,GO,Node>(numLocalColumns * numRanks, numLocalColumns, comm);
    //Create a matrix that is dense in the last row on each proc, but has 5 entries in every other row
    Teuchos::Array<size_t> entriesPerRow(numLocalRows);
    for(size_t i = 0; i < numLocalRows - 1; i++)
      entriesPerRow[i] = 5;
    entriesPerRow[numLocalRows - 1] = numLocalColumns;
    MAT imba(rowMap, colMap, entriesPerRow());
    //Insert 1 as value in all entries.
    Teuchos::Array<Scalar> vals(numLocalColumns, ST::one());
    Teuchos::Array<LO> cols(numLocalColumns);
    for(size_t i = 0; i < numLocalColumns; i++)
      cols[i] = i;
    auto valsView = vals();
    auto colsView = cols();
    for(size_t i = 0; i < numLocalRows - 1; i++)
      imba.insertLocalValues(i, colsView.view(0, 5), valsView.view(0, 5));
    imba.insertLocalValues(numLocalRows - 1, colsView, valsView);
    imba.fillComplete(colMap, rowMap);
    //auto fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //imba.describe(*fancy,Teuchos::VERB_EXTREME);
    auto domainMap = imba.getDomainMap();
    //Input vector: elements are all 1s
    MV w(domainMap, numVecs, false);
    w.putScalar(ST::one());
    //Output vector: zero initially
    //note: rowMap is non-overlapping so this is OK
    MV v(rowMap, numVecs, true);
    //Do the apply. If cuSPARSE is enabled, the merge path algorithm will be used. Otherwise, this tests the fallback.
    imba.apply(w, v);
    Teuchos::Array<Scalar> vvals(numLocalRows * numVecs);
    v.get1dCopy(vvals(), numLocalRows);
    //Each output value should be equal to the number of entries in the row.
    //These are smallish integers so they should be represented exactly in 32- or 64-bit floating point
    for(size_t vec = 0; vec < numVecs; vec++)
    {
      size_t vecOffset = vec * numLocalRows;
      for(size_t i = 0; i < numLocalRows - 1; i++)
      {
        TEST_EQUALITY(static_cast<Mag>(5.0) * ST::one(), vvals[vecOffset + i]);
      }
      TEST_EQUALITY(static_cast<Mag>(numLocalColumns) * ST::one(), vvals[vecOffset + numLocalRows - 1]);
    }
    if(numVecs != 1)
    {
      //Now run again, but on single vectors only (rank-1)
      auto wcol = w.getVector(0);
      auto vcol = v.getVectorNonConst(0);
      vcol->putScalar(ST::zero());
      imba.apply(*wcol, *vcol);
      vcol->get1dCopy(vvals(), numLocalRows);
      for(size_t i = 0; i < numLocalRows - 1; i++)
      {
        TEST_EQUALITY(static_cast<Mag>(5.0) * ST::one(), vvals[i]);
      }
      TEST_EQUALITY(static_cast<Mag>(numLocalColumns) * ST::one(), vvals[numLocalRows - 1]);
      //Finally, test residual.
      V res(rowMap);
      //Here, have A*wcol = vcol. This means the residual of A, wcol, and vcol should be 0.
      Tpetra::Details::residual(imba, *wcol, *vcol, res);
      Teuchos::Array<Scalar> resVals(numLocalRows);
      res.get1dCopy(resVals(), numLocalRows);
      for(size_t i = 0; i < numLocalRows; i++)
      {
        TEST_EQUALITY(ST::zero(), resVals[i]);
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TheEyeOfTruth,  LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ZeroMatrix,     LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ImbalancedRowMatrix, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, BadCalls,       LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, SimpleEigTest,  LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
