// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {

  using std::endl;
  using std::string;

  using Teuchos::as;
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
  using Tpetra::Export;
  using Tpetra::global_size_t;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createVector;
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
  template <typename Scalar, typename LO, typename GO, typename Node>
  void applyAndCheckResult(
    Tpetra::CrsMatrix<Scalar, LO, GO, Node> &A,
    Tpetra::CrsMatrix<Scalar, LO, GO, Node> &B,
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node> > &newMap,
    Teuchos::FancyOStream &out,
    bool &success
  )
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    MT norm = ScalarTraits<MT>::zero ();

    /* Fill a random vector on the original map */
    Tpetra::Vector<Scalar,LO,GO,Node> vecX(A.getDomainMap());
    vecX.randomize();

    /* Now do some multiplies */
    Tpetra::Vector<Scalar,LO,GO,Node> AVecY(A.getRangeMap());
    Tpetra::Vector<Scalar,LO,GO,Node> BVecY(B.getRangeMap());
    A.apply(vecX,AVecY);
    B.apply(vecX,BVecY);

    /* Export BVecY to the original range map for comparisons */
    Tpetra::Vector<Scalar,LO,GO,Node> BVecYOrig(A.getRangeMap());
    Tpetra::Export<LO,GO,Node> TempExport(newMap, A.getRangeMap());
    BVecYOrig.doExport(BVecY,TempExport,Tpetra::INSERT);

    BVecYOrig.update (-STS::one (), AVecY, STS::one ());
    norm = BVecYOrig.norm2();

    std::cout << "Residual 2-norm: " << norm << std::endl
              << "Residual 1-norm: " << BVecYOrig.norm1 () << std::endl
              << "Residual Inf-norm: " << BVecYOrig.normInf () << std::endl;

    const bool normSmallEnough = (norm <= as<MT> (1e-10));
    TEST_EQUALITY ( normSmallEnough, true );
  }

  //////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceRangeMap, LO, GO, Scalar, Node )
  {
    // Based on the FullTriDiag tests...

    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> STS;

    const size_t ONE  = OrdinalTraits<size_t>::one();
    const size_t ZERO = OrdinalTraits<GO>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map =
              createContigMapWithNode<LO,GO,Node>(INVALID,ONE,comm);

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
    B.setObjectLabel("The Other Matrix");
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*STS::one(),
                                       STS::one(),
                                       static_cast<Scalar>(2)*STS::one()));
      Array<GO> cols(tuple<GO>(myImageID,myImageID + 1));
      A.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      A.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
      B.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      B.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
    }
    A.fillComplete();
    B.fillComplete();


    // Build a one-process Map.
    // we know the map is contiguous...
    const size_t NumMyElements = (comm->getRank () == 0) ?
      A.getRangeMap ()->getGlobalNumElements () : 0;
    RCP<const Map<LO,GO,Node> > NewMap =
      rcp (new Map<LO,GO,Node> (INVALID, NumMyElements, ZERO, comm));

    B.replaceRangeMap (NewMap);

    applyAndCheckResult(A, B, NewMap, out, success);
  }

  //////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceRangeMapAndExporter, LO, GO, Scalar, Node )
  {
    // Based on the FullTriDiag tests...

    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> STS;

    const size_t ONE  = OrdinalTraits<size_t>::one();
    const size_t ZERO = OrdinalTraits<GO>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map =
              createContigMapWithNode<LO,GO,Node>(INVALID,ONE,comm);

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
    B.setObjectLabel("The Other Matrix");
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*STS::one(),
                                       STS::one(),
                                       static_cast<Scalar>(2)*STS::one()));
      Array<GO> cols(tuple<GO>(myImageID,myImageID + 1));
      A.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      A.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
      B.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      B.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
    }
    A.fillComplete();
    B.fillComplete();

    // we know the map is contiguous...
    const size_t NumMyElements = (comm->getRank () == 0) ?
      A.getRangeMap ()->getGlobalNumElements () : 0;
    RCP<const Map<LO,GO,Node> > NewMap =
      rcp (new Map<LO,GO,Node> (INVALID, NumMyElements, ZERO, comm));
    RCP<const Tpetra::Export<LO,GO,Node> > NewExport =
      rcp (new Export<LO,GO,Node> (B.getRowMap (), NewMap));

    B.replaceRangeMapAndExporter (NewMap, NewExport);

    applyAndCheckResult(A, B, NewMap, out, success);
  }

  //////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, RangeMapEqualsRowMap, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> STS;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map =
              createContigMapWithNode<LO,GO,Node>(INVALID,2,comm);

    /* Create the following matrix:
    0  [2 1       ]   [2 1]
    1  [1 2       ]   [1 2]
    2  [    2 1   ]         + [2 1]
    3  [    1 2   ] =         [1 2] +
       [       2 1]
  2n-1 [       1 2]
    */

    MAT A(map,2);
    MAT B(map,2);
    A.setObjectLabel("The Matrix");
    B.setObjectLabel("The Other Matrix");
    Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*STS::one(),
                                     STS::one(),
                                     static_cast<Scalar>(2)*STS::one()));
    Array<GO> cols(tuple<GO>(2*myImageID,2*myImageID + 1));
    A.insertGlobalValues(2*myImageID  ,cols(),vals(0,2)); // insert [2 1]
    A.insertGlobalValues(2*myImageID+1,cols(),vals(1,2)); // insert [1 2]
    B.insertGlobalValues(2*myImageID  ,cols(),vals(0,2)); // insert [2 1]
    B.insertGlobalValues(2*myImageID+1,cols(),vals(1,2)); // insert [1 2]
    A.fillComplete();
    B.fillComplete();

    // Use the row map to exercise the check for Range Map == Row Map
    // in replaceRangeMap; By construction, Column Map is 1-to-1
    RCP<const Map<LO,GO,Node> > NewMap = B.getRowMap();

    B.replaceRangeMap (NewMap);

    applyAndCheckResult(A, B, NewMap, out, success);
  }
//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceRangeMap, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceRangeMapAndExporter, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, RangeMapEqualsRowMap, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
