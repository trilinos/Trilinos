// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MixedScalarMultiplyOp.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_getNumDiags.hpp"
#include "Tpetra_Details_residual.hpp"

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


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MixedScalarMultiplyOp, Apply, LO, GO, Scalar, Node )
  {
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Tpetra::Operator<Scalar,LO,GO,Node> OP;
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
    MV ones(map,ONE,false), threes1(map,ONE,false), threes2(map,ONE,false);
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
    MAT A(map,3);
    if (myImageID == 0) {
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*ST::one(), ST::one()));
      Array<GO> cols(tuple<GO>(myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      Array<Scalar> vals(tuple<Scalar>(ST::one(), static_cast<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID-1,myImageID));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else {
      Array<Scalar> vals(3,ST::one());
      Array<GO> cols(tuple<GO>(myImageID-1, myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    A.fillComplete();

    // Now, wrap A into a MixedScalarMultiplyOp. We generally build
    // only with a single scalar type, in which case we test wrapping
    // a scalar operator as a scalar operator.
    // But if we do have double and float, we can convert A to a float
    // matrix and then wrap it back to double.
    RCP<OP> mixOpA;
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
    if (typeid(Scalar) == typeid(double)) {
      typedef Tpetra::MixedScalarMultiplyOp<double,float,LO,GO,Node> MixedOp;
      mixOpA = Teuchos::rcp_dynamic_cast<OP>(rcp(new MixedOp(A.template convert<float>())),true);
    } else
#endif
    {
      typedef Tpetra::MixedScalarMultiplyOp<Scalar,Scalar,LO,GO,Node> MixedOp;
      mixOpA = rcp(new MixedOp(rcpFromRef(A)));
    }

    // test the properties of mixOpA
    TEST_EQUALITY_CONST(mixOpA->getDomainMap()->isSameAs(*A.getDomainMap()), true);
    TEST_EQUALITY_CONST(mixOpA->getRangeMap()->isSameAs(*A.getRangeMap()), true);

    // test the action
    threes1.randomize();
    threes2.randomize();
    A.apply(ones,threes1);
    mixOpA->apply(ones,threes2);

    // Both vectors should now be all threes.
    // Check threes2 against threes1.
    threes2.update(-ST::one(),threes1,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    threes2.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }

    // now, threes1 should be 3*ones
    threes1.update(static_cast<Scalar>(-3)*ST::one(),ones,ST::one());
    threes1.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

//
// INSTANTIATIONS
//

// FIXME_SYCL
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MixedScalarMultiplyOp, Apply,  LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

}
