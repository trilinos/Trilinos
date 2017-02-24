#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "BelosThyraAdapter.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"

// Typedefs

typedef Tpetra::Vector<>::scalar_type Scalar;
typedef Tpetra::Vector<>::local_ordinal_type LO;
typedef Tpetra::Vector<>::global_ordinal_type GO;
typedef Tpetra::Vector<>::node_type Node;

typedef Teuchos::ScalarTraits<Scalar> ST;
typedef Tpetra::Map<LO,GO,Node> Map;
typedef Tpetra::Vector<Scalar,LO,GO,Node> TV;
typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> CRSM;
typedef Thyra::VectorSpaceBase<Scalar> TVSB;
typedef Thyra::VectorBase<Scalar> TVB;
typedef Thyra::MultiVectorBase<Scalar> TMVB;
typedef Thyra::LinearOpBase<Scalar> TLOB;
typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node> TOVE;
typedef typename TMV::mag_type mag_type;

#if defined HAVE_TPETRACORE_CUDA
#define NUM_LOCAL 100
#else
#define NUM_LOCAL 100
#endif
const std::size_t numLocalElements = NUM_LOCAL;

// Tpetra matrix creation

Teuchos::RCP<CRSM> createTpetraOp(const Teuchos::RCP<const Map>& map)
{
  Teuchos::RCP<CRSM> mat = Teuchos::rcp(new CRSM(map, 2));

  for (std::size_t i = 0; i < map->getNodeNumElements(); ++i) {
    GO row = map->getGlobalElement(static_cast<LO>(i));
    Scalar val = static_cast<Scalar>(row+1);
    if (row == map->getGlobalNumElements()-1) {
      mat->insertGlobalValues(row,
                              Teuchos::tuple(row),
                              Teuchos::tuple(val));
    } else {
      mat->insertGlobalValues(row,
                              Teuchos::tuple(row, row+1),
                              Teuchos::tuple(val, -1.0*val));
    }
  }

  mat->fillComplete();
  return mat;
}

//Routines for checking solution

void checkMultiVectors(const Teuchos::RCP<TMV>& a,
                       const Teuchos::RCP<TMV>& b,
                       const Teuchos::ArrayView<mag_type>& expectedNorms,
                       Teuchos::FancyOStream& out,
                       bool& success)
{
  TEUCHOS_ASSERT(a->getNumVectors() == b->getNumVectors());
  TEUCHOS_ASSERT(a->getNumVectors() == expectedNorms.size());
  Scalar one = ST::one();
  b->update(-1.0*one, *a, one);
  std::vector<ST::magnitudeType> norms(expectedNorms.size());
  b->norm2(Teuchos::arrayViewFromVector(norms));

  ST::magnitudeType tol = 1.0e-14;
  ST::magnitudeType zero = ST::magnitude(ST::zero());
  for (auto norm = norms.begin(); norm != norms.end(); ++norm) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*norm, zero, tol, out, success);
  }
  a->norm2(Teuchos::arrayViewFromVector(norms));
  auto normIter = norms.begin();
  auto expNormIter = expectedNorms.begin();
  for (; normIter != norms.end(); ++normIter, ++expNormIter) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*normIter, *expNormIter, tol, out, success);
  }
}

// Unit tests

TEUCHOS_UNIT_TEST(Belos_LOWS, Factory)
{
  // Create Belos solve strategy with no prec
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  // Create the LOWSF
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory
    = builder.createLinearSolveStrategy("");
  Teuchos::RCP<Thyra::BelosLinearOpWithSolveFactory<Scalar> > factoryBelos
    = Teuchos::rcp_dynamic_cast<Thyra::BelosLinearOpWithSolveFactory<Scalar> >(factory);
  TEUCHOS_TEST_EQUALITY(nonnull(factoryBelos), true, out, success);

  // Create the LOWS
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > lows = factory->createOp();
  Teuchos::RCP<Thyra::BelosLinearOpWithSolveFactory<Scalar> > lowsBelos
    = Teuchos::rcp_dynamic_cast<Thyra::BelosLinearOpWithSolveFactory<Scalar> >(factory);
  TEUCHOS_TEST_EQUALITY(nonnull(lowsBelos), true, out, success);

  // Create a new solve strategy with Ifpack2 prec
  builder.unsetParameterList();
  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::Ifpack2PreconditionerFactory<CRSM> Impl;
  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
  p->set("Preconditioner Type", "Ifpack2");
  Teuchos::ParameterList& ifpackList = p->sublist("Preconditioner Types").sublist("Ifpack2");
  ifpackList.set("Prec Type", "ILUT");
  builder.setParameterList(p);

  // Create the LOWSF
  factory = builder.createLinearSolveStrategy("");
  factoryBelos = Teuchos::rcp_dynamic_cast<Thyra::BelosLinearOpWithSolveFactory<Scalar> >(factory);
  TEUCHOS_TEST_EQUALITY(nonnull(factoryBelos), true, out, success);

  // Create the LOWS
  lows = factory->createOp();
  lowsBelos = Teuchos::rcp_dynamic_cast<Thyra::BelosLinearOpWithSolveFactory<Scalar> >(factory);
  TEUCHOS_TEST_EQUALITY(nonnull(lowsBelos), true, out, success);
  Teuchos::RCP<Impl> precFactory = Teuchos::rcp_dynamic_cast<Impl>(factory->getPreconditionerFactory());
  TEUCHOS_TEST_EQUALITY(nonnull(precFactory), true, out, success);
}

TEUCHOS_UNIT_TEST(Belos_LOWS, FwdOp_Vec)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create Belos solve strategy with no prec
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  // Create the LOWS
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory
    = builder.createLinearSolveStrategy("");
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > lows = factory->createOp();

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);

  // Initialize the LOWS
  Thyra::initializeOp<Scalar>(*factory, op, lows.ptr());
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->domain()), true, out, success);
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->range()), true, out, success);

  // Create vector objects
  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);

  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  // Apply the forward op
  Thyra::apply<Scalar>(*lows, Thyra::NOTRANS, *x, y.ptr(), ST::one(), ST::one());

  // Check for correct result
  if (comm->getRank() == comm->getSize()-1) {
    LO row = static_cast<LO>(numLocalElements-1);
    Scalar val = static_cast<Scalar>(numGlobalElements);
    x_tpetra->sumIntoLocalValue(row, val);
  }
  mag_type val = static_cast<mag_type>(static_cast<Scalar>(numGlobalElements)
    *ST::squareroot(static_cast<Scalar>(1.0) + 3.0/numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_LOWS, FwdOp_MultiVec)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create Belos solve strategy with no prec
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  // Create the LOWS
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory
    = builder.createLinearSolveStrategy("");
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > lows = factory->createOp();

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);

  // Initialize the LOWS
  Thyra::initializeOp<Scalar>(*factory, op, lows.ptr());
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->domain()), true, out, success);
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->range()), true, out, success);

  // Create vector objects
  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);

  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  auto scales = Teuchos::tuple(one, 2.0*one, 3.0*one, 4.0*one, 5.0*one);
  x_tpetra->scale(scales);
  y_tpetra->scale(scales);

  // Apply the forward op
  Thyra::apply<Scalar>(*lows, Thyra::NOTRANS, *x, y.ptr(), ST::one(), ST::one());

  // Check for correct result
  if (comm->getRank() == comm->getSize()-1) {
    LO row = static_cast<LO>(numLocalElements-1);
    for (std::size_t col = 0; col < numCols; ++col) {
      Scalar val = static_cast<Scalar>((col+1)*numGlobalElements);
      x_tpetra->sumIntoLocalValue(row, col, val);
    }
  }
  mag_type val = static_cast<mag_type>(static_cast<Scalar>(numGlobalElements)
    *ST::squareroot(static_cast<Scalar>(1.0) + 3.0/numGlobalElements));
  auto ans = Teuchos::tuple(val, 2.0*val, 3.0*val, 4.0*val, 5.0*val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_LOWS, InvOp_Vec_NoPrec)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create Belos solve strategy with no prec
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", numGlobalElements+1);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", numGlobalElements+1);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  // Create the LOWS
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory
    = builder.createLinearSolveStrategy("");
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > lows = factory->createOp();

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);

  // Initialize the LOWS
  Thyra::initializeOp<Scalar>(*factory, op, lows.ptr());
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->domain()), true, out, success);
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->range()), true, out, success);

  // Create vector objects
  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);

  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  // Setup random initial guess and RHS with known answer
  x_tpetra->randomize(static_cast<Scalar>(-1.0), static_cast<Scalar>(1.0));
  y_tpetra->putScalar(ST::zero());
  if (comm()->getRank() == comm->getSize()-1) {
    LO row = static_cast<LO>(numLocalElements-1);
    Scalar val = static_cast<Scalar>(numGlobalElements);
    y_tpetra->sumIntoLocalValue(row, val);
  }

  // Create solve criteria
  Thyra::SolveCriteria<Scalar> criteria;
  criteria.requestedTol = 1.0e-8;
  criteria.solveMeasureType = Thyra::SolveMeasureType(Thyra::SOLVE_MEASURE_NORM_RESIDUAL,
    Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL);

  // Apply the inverse op
  Thyra::SolveStatus<Scalar> status;
  status = Thyra::solve<Scalar>(*lows, Thyra::NOTRANS, *y, x.ptr(), Teuchos::constPtr(criteria));
  //x_tpetra->describe(*Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);
  //y_tpetra->describe(*Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);
  TEUCHOS_ASSERT(status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);

  // Check for correct result
  x_tpetra->putScalar(ST::one());
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(numGlobalElements)));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_LOWS, InvOp_MultiVec_NoPrec)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create Belos solve strategy with no prec
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", numGlobalElements+1);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", numGlobalElements+1);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  // Create the LOWS
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory
    = builder.createLinearSolveStrategy("");
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > lows = factory->createOp();

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);

  // Initialize the LOWS
  Thyra::initializeOp<Scalar>(*factory, op, lows.ptr());
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->domain()), true, out, success);
  TEUCHOS_TEST_EQUALITY(space->isCompatible(*lows->range()), true, out, success);

  // Create vector objects
  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);

  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  // Setup random initial guess and RHS with known answer
  x_tpetra->randomize(static_cast<Scalar>(-1.0), static_cast<Scalar>(1.0));
  y_tpetra->putScalar(ST::zero());
  if (comm()->getRank() == comm->getSize()-1) {
    LO row = static_cast<LO>(numLocalElements-1);
    for (std::size_t col = 0; col < numCols; ++col) {
      Scalar val = static_cast<Scalar>((col+1)*numGlobalElements);
      y_tpetra->sumIntoLocalValue(row, col, val);
    }
  }

  // Create solve criteria
  Thyra::SolveCriteria<Scalar> criteria;
  criteria.requestedTol = 1.0e-8;
  criteria.solveMeasureType = Thyra::SolveMeasureType(Thyra::SOLVE_MEASURE_NORM_RESIDUAL,
    Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL);

  // Apply the inverse op
  Thyra::SolveStatus<Scalar> status;
  status = Thyra::solve<Scalar>(*lows, Thyra::NOTRANS, *y, x.ptr(), Teuchos::constPtr(criteria));
  TEUCHOS_ASSERT(status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);

  // Check for correct result
  /*x_tpetra->putScalar(ST::one());
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(numGlobalElements)));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);*/
}

