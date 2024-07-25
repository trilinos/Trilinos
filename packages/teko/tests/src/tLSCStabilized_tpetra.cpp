// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tLSCStabilized_tpetra.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teko_Utilities.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
//
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  0  0 ]
//      [ -1  1  0  0 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;
using namespace Teko::NS;

void tLSCStabilized_tpetra::initializeTest() {
  tolerance_ = 1.0e-13;

  comm = GetComm_tpetra();
}

int tLSCStabilized_tpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                   int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tLSCStabilized_tpetra";

  status = test_diagonal(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal\" ... PASSED", "   \"diagonal\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
  /*
     status = test_diagonalNotSym(verbosity,failstrm);
     Teko_TEST_MSG_tpetra(stdstrm,1,"   \"diagonalNotSym\" ... PASSED","   \"diagonalNotSym\" ...
     FAILED"); allTests &= status; failcount += status ? 0 : 1; totalrun++;

     status = test_strategy(verbosity,failstrm);
     Teko_TEST_MSG_tpetra(stdstrm,1,"   \"strategy\" ... PASSED","   \"strategy\" ... FAILED");
     allTests &= status;
     failcount += status ? 0 : 1;
     totalrun++;
  */

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tLSCStabilized_tpetra...PASSED",
                         "tLSCStabilized_tpetra...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tLSCStabilized_tpetra...FAILED");
  }

  return failcount;
}

bool tLSCStabilized_tpetra::test_diagonal(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  ST vec[2];
  ST diff = 0.0;

  // build 4x4 matrix with block 2x2 diagonal subblocks
  //
  //            [ 1 0 7 0 ]
  // [ F G ] =  [ 0 2 0 8 ]
  // [ D C ]    [ 5 0 3 0 ]
  //            [ 0 6 0 4 ]
  //

  vec[0]     = 1.0;
  vec[1]     = 2.0;
  LinearOp F = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  LinearOp G = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  LinearOp D = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 3.0;
  vec[1]     = 4.0;
  LinearOp C = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]        = 0.030303030303030;
  vec[1]        = 0.02205882352941;
  LinearOp iBBt = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]       = 0.026041666666667;
  vec[1]       = 0.041666666666667;
  LinearOp aiD = Teko::Test::DiagMatrix_tpetra(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);

  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(iF, iBBt, aiD, Teuchos::null));
  RCP<Thyra::PreconditionerBase<ST> > prec = Thyra::prec<ST>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<ST> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));
  // construct a couple of vectors
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  Tpetra::Vector<ST, LO, GO, NT> ef(map), eg(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = BlockVector(ef, eg, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);

  ef.replaceGlobalValue(0, 0.407268709825528);
  ef.replaceGlobalValue(1, 1.560553633217993);
  eg.replaceGlobalValue(0, -0.058181244260790);
  eg.replaceGlobalValue(1, -0.265138408304498);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [-2 4 7 9]
  ea.replaceGlobalValue(0, -2.0);
  ea.replaceGlobalValue(1, 4.0);
  eb.replaceGlobalValue(0, 7.0);
  eb.replaceGlobalValue(1, 9.0);

  ef.replaceGlobalValue(0, 0.850880968778696);
  ef.replaceGlobalValue(1, 5.181660899653979);
  eg.replaceGlobalValue(0, -0.407268709825528);
  eg.replaceGlobalValue(1, -0.795415224913495);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [1 0 0 -5]
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, -5.0);

  ef.replaceGlobalValue(0, 1.000000000000000);
  ef.replaceGlobalValue(1, -1.767589388696655);
  eg.replaceGlobalValue(0, 0.000000000000000);
  eg.replaceGlobalValue(1, 0.441897347174164);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [4 -4 6 12]
  ea.replaceGlobalValue(0, 4.0);
  ea.replaceGlobalValue(1, -4.0);
  eb.replaceGlobalValue(0, 6.0);
  eb.replaceGlobalValue(1, 12.0);

  ef.replaceGlobalValue(0, 6.443612258953168);
  ef.replaceGlobalValue(1, 2.242214532871971);
  eg.replaceGlobalValue(0, -0.349087465564738);
  eg.replaceGlobalValue(1, -1.060553633217993);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  return allPassed;
}

bool tLSCStabilized_tpetra::test_diagonalNotSym(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  ST vec[2];
  ST diff = 0.0;

  // build 4x4 matrix with block 2x2 diagonal subblocks
  //
  //            [ 1 0 7 0 ]
  // [ F G ] =  [ 0 2 0 8 ]
  // [ D C ]    [ 5 0 3 0 ]
  //            [ 0 6 0 4 ]
  //

  vec[0]     = 1.0;
  vec[1]     = 2.0;
  LinearOp F = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  LinearOp G = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  LinearOp D = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 3.0;
  vec[1]     = 4.0;
  LinearOp C = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]        = 0.030303030303030;
  vec[1]        = 0.02205882352941;
  LinearOp iBBt = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]       = 0.026041666666667;
  vec[1]       = 0.041666666666667;
  LinearOp aiD = Teko::Test::DiagMatrix_tpetra(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);

  RCP<InverseLibrary> invLib  = InverseLibrary::buildFromStratimikos();
  RCP<InverseFactory> invFact = invLib->getInverseFactory("Ifpack2");

  RCP<InvLSCStrategy> lscStrat = rcp(new InvLSCStrategy(invFact));
  // lscStrat->setSymmetric(false);
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(lscStrat));
  RCP<Thyra::PreconditionerBase<ST> > prec = Thyra::prec<ST>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<ST> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));
  // construct a couple of vectors
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  Tpetra::Vector<ST, LO, GO, NT> ef(map), eg(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = BlockVector(ef, eg, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);

  ef.replaceGlobalValue(0, 0.407268709825528);
  ef.replaceGlobalValue(1, 1.560553633217993);
  eg.replaceGlobalValue(0, -0.058181244260790);
  eg.replaceGlobalValue(1, -0.265138408304498);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [-2 4 7 9]
  ea.replaceGlobalValue(0, -2.0);
  ea.replaceGlobalValue(1, 4.0);
  eb.replaceGlobalValue(0, 7.0);
  eb.replaceGlobalValue(1, 9.0);

  ef.replaceGlobalValue(0, 0.850880968778696);
  ef.replaceGlobalValue(1, 5.181660899653979);
  eg.replaceGlobalValue(0, -0.407268709825528);
  eg.replaceGlobalValue(1, -0.795415224913495);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [1 0 0 -5]
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, -5.0);

  ef.replaceGlobalValue(0, 1.000000000000000);
  ef.replaceGlobalValue(1, -1.767589388696655);
  eg.replaceGlobalValue(0, 0.000000000000000);
  eg.replaceGlobalValue(1, 0.441897347174164);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  // test vector [4 -4 6 12]
  ea.replaceGlobalValue(0, 4.0);
  ea.replaceGlobalValue(1, -4.0);
  eb.replaceGlobalValue(0, 6.0);
  eb.replaceGlobalValue(1, 12.0);

  ef.replaceGlobalValue(0, 6.443612258953168);
  ef.replaceGlobalValue(1, 2.242214532871971);
  eg.replaceGlobalValue(0, -0.349087465564738);
  eg.replaceGlobalValue(1, -1.060553633217993);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  TEST_ASSERT((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_,
              "   tLSCStabilized_tpetra::test_diagonal "
                  << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff
                  << " <= " << tolerance_ << ")\n"
                  << "      " << Print("x", x) << "      " << Print("y", y) << "      "
                  << Print("z", z));

  return allPassed;
}

bool tLSCStabilized_tpetra::test_strategy(int verbosity, std::ostream& os) {
  std::vector<LO> indicies(2);
  std::vector<ST> row0(2);
  int sz = 5;
  ST vec[5];

  bool status    = false;
  bool allPassed = true;

  vec[0]     = 1.0;
  vec[1]     = 2.0;
  vec[2]     = 3.0;
  vec[3]     = 4.0;
  vec[4]     = 5.0;
  LinearOp F = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  vec[2]     = 9.0;
  vec[3]     = 10.0;
  vec[4]     = 11.0;
  LinearOp G = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  vec[2]     = 7.0;
  vec[3]     = 8.0;
  vec[4]     = 9.0;
  LinearOp D = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]     = 3.0;
  vec[1]     = 4.0;
  vec[2]     = 5.0;
  vec[3]     = 6.0;
  vec[4]     = 7.0;
  LinearOp C = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]      = 1.0;
  vec[1]      = 1.0 / 2.0;
  vec[2]      = 1.0 / 3.0;
  vec[3]      = 1.0 / 4.0;
  vec[4]      = 1.0 / 5.0;
  LinearOp iF = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]           = 0.091304347826087;
  vec[1]           = 0.090517241379310;
  vec[2]           = 0.087646076794658;
  vec[3]           = 0.084000000000000;
  vec[4]           = 0.080152671755725;
  LinearOp iBQBtmC = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]       = 0.020202020202020;
  vec[1]       = 0.032323232323232;
  vec[2]       = 0.040404040404040;
  vec[3]       = 0.046176046176046;
  vec[4]       = 0.050505050505051;
  LinearOp aiD = Teko::Test::DiagMatrix_tpetra(sz, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);

  comm                                    = GetComm_tpetra();
  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(sz, 0, comm));

  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  RCP<InverseLibrary> invLib  = InverseLibrary::buildFromStratimikos();
  RCP<InverseFactory> invFact = invLib->getInverseFactory("Ifpack2");

  // build Mass matrix
  vec[0]        = 3.0;
  vec[1]        = 4.0;
  vec[2]        = 5.0;
  vec[3]        = 6.0;
  vec[4]        = 7.0;
  LinearOp mass = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]           = 1.0 / 3.0;
  vec[1]           = 1.0 / 4.0;
  vec[2]           = 1.0 / 5.0;
  vec[3]           = 1.0 / 6.0;
  vec[4]           = 1.0 / 7.0;
  LinearOp invMass = Teko::Test::DiagMatrix_tpetra(sz, vec);

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1.2e-2);
  tester.show_all_tests(true);
  std::stringstream ss;
  Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss), "      |||");

  Teko::BlockedLinearOp blkA = Teko::toBlockedLinearOp(A);

  // build preconditioner
  vec[0]       = 1.0;
  vec[1]       = 0.5;
  vec[2]       = 1.0 / 3.0;
  vec[3]       = 0.25;
  vec[4]       = 0.2;
  LinearOp p00 = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]       = 0.368351759561589;
  vec[1]       = 0.325933832979017;
  vec[2]       = 0.295436133965709;
  vec[3]       = 0.272240115440115;
  vec[4]       = 0.253891252128534;
  LinearOp p01 = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]       = 0;
  vec[1]       = 0;
  vec[2]       = 0;
  vec[3]       = 0;
  vec[4]       = 0;
  LinearOp p10 = Teko::Test::DiagMatrix_tpetra(sz, vec);

  vec[0]       = -0.052621679937370;
  vec[1]       = -0.081483458244754;
  vec[2]       = -0.098478711321903;
  vec[3]       = -0.108896046176046;
  vec[4]       = -0.115405114603879;
  LinearOp p11 = Teko::Test::DiagMatrix_tpetra(sz, vec);
  LinearOp P   = Thyra::block2x2(p00, p01, p10, p11);

  // Kluge to get around problem with Anasazi
  // Teko::computeSpectralRad(Thyra::multiply(invMass,F),5e-2,false,3)/3.0;
  // Teko::computeSpectralRad(Thyra::multiply(invMass,F),5e-2,false,3)/3.0;

  // build inverse strategy
  {
    bool result;
    Teko::NS::LSCPrecondState state;
    Teko::NS::InvLSCStrategy iStrat(invFact, mass, false);
    iStrat.setEigSolveParam(3);
    Teko::NS::LSCPreconditionerFactory factory(Teuchos::rcpFromRef(iStrat));
    LinearOp prec = factory.buildPreconditionerOperator(blkA, state);

    // test inverse mass
    ss.str("");
    result = tester.compare(*invMass, *iStrat.getInvMass(blkA, state), Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLSCStabilized_tpetra::test_strategy " << toString(status)
                            << " : Comparing mass operators");
    if (not result || verbosity >= 10) os << ss.str();

    // test inverse F
    ss.str("");
    result = tester.compare(*iF, *iStrat.getInvF(blkA, state), Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLSCStabilized_tpetra::test_strategy " << toString(status)
                            << " : Comparing F operators");
    if (not result || verbosity >= 10) os << ss.str();

    // test inverse B*Q*Bt-gamma*C
    ss.str("");
    result = tester.compare(*iBQBtmC, *iStrat.getInvBQBt(blkA, state), Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLSCStabilized_tpetra::test_strategy " << toString(status)
                            << " : Comparing B*Q*Bt-C operators");
    if (not result || verbosity >= 10) os << ss.str();

    // test alpha*inv(D)
    ss.str("");

#ifndef TEKO_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D

    // result = tester.compare( *aiD, *iStrat.getInvAlphaD(blkA,state), &fos );
    result =
        tester.compare(*aiD, *iStrat.getOuterStabilization(blkA, state), Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLSCStabilized_tpetra::test_strategy " << toString(status)
                            << " : Comparing alpha*inv(D) operators");
    if (not result || verbosity >= 10) os << ss.str();

#endif

    // test full op
    ss.str("");
    result = tester.compare(*P, *prec, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLSCStabilized_tpetra::test_strategy " << toString(status)
                            << " : Comparing full op");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
