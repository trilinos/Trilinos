// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include <cstdlib>

using Teuchos::RCP;
using Teuchos::rcp;
typedef Tpetra::Map<int, int> map_type;
typedef map_type::local_ordinal_type LO;
typedef map_type::global_ordinal_type GO;
typedef double ST;
typedef Tpetra::CrsMatrix<ST, LO, GO> crs_matrix_type;
typedef Tpetra::Operator<ST, LO, GO> op_type;
typedef Tpetra::Vector<ST, LO, GO> vec_type;

Teuchos::RCP<const crs_matrix_type>
create_matrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  const Tpetra::global_size_t numGlobalElements = 100;

  const GO indexBase = 0;
  auto map = rcp (new map_type (numGlobalElements, indexBase, comm));

  const LO numMyElements = map->getLocalNumElements ();
  auto myGlobalElements = map->getLocalElementList ();
  auto A = rcp (new crs_matrix_type (map, 1));

  for (LO lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const GO gblInd = map->getGlobalElement (lclRow);
    const ST val = static_cast<ST> (1.0) / static_cast<ST> (gblInd + 1);
    A->insertGlobalValues (gblInd,
                           Teuchos::tuple (gblInd),
                           Teuchos::tuple (val));
  }
  A->fillComplete ();
  return A;
}

int
main (int argc, char *argv[])
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;
  int lclSuccess = 1;
  int gblSuccess = 0;

  Teuchos::GlobalMPISession session (&argc, &argv, NULL);
  auto out = Teuchos::VerboseObjectBase::getDefaultOStream ();
  const auto comm = Teuchos::DefaultComm<int>::getComm ();
  const int myRank = comm->getRank ();

  *out << "Creating matrix" << endl;
  RCP<const crs_matrix_type> A;
  try {
    A = create_matrix (comm);
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": create_matrix(comm) threw an "
      "exception: " << e.what () << endl;
    cerr << os.str ();
  }
  catch (...) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": create_matrix(comm) threw an "
      "exception not a subclass of std::exception." << endl;
    cerr << os.str ();
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "create_matrix(comm) threw an exception on some process."
          << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  // Make sure that A is nonnull on all processes.
  if (A.is_null ()) {
    lclSuccess = 0;
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "The result of create_matrix(comm) is null on at least one "
        "process." << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  *out << "Creating vectors" << endl;
  vec_type b (A->getRangeMap ());
  b.putScalar (1.0);
  vec_type x (A->getDomainMap ());
  x.putScalar (0.0);

  *out << "Creating Stratimikos linear solver builder" << endl;
  Stratimikos::DefaultLinearSolverBuilder builder;
  auto p = Teuchos::rcp (new Teuchos::ParameterList ());
  try {
    builder.setParameterList (p);
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": builder.setParameterList(p) threw an "
      "exception: " << e.what () << endl;
    cerr << os.str ();
  }
  catch (...) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": builder.setParameterList(p) threw an "
      "exception not a subclass of std::exception." << endl;
    cerr << os.str ();
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "builder.setParameterList(p) threw an exception on some process."
          << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  *out << "Calling builder.createLinearSolveStrategy" << endl;
  auto lowsFactory = builder.createLinearSolveStrategy ("");
  lowsFactory->setVerbLevel (Teuchos::VERB_LOW);

  *out << "Calling Thyra::createConstLinearOp" << endl;
  const op_type& opA = *A;
  auto thyraA = Thyra::createConstLinearOp (Teuchos::rcpFromRef (opA));

  // using Teuchos::rcp_implicit_cast;
  // auto thyraA = Thyra::createConstLinearOp (rcp_implicit_cast<const op_type> (A));

  // Make sure that thyraA is nonnull on all processes.
  if (thyraA.is_null ()) {
    lclSuccess = 0;
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "The result of Thyra::createConstLinearOp is null on at least "
        "one process." << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  *out << "Creating Thyra Ifpack2 factory" << endl;
  RCP<Thyra::PreconditionerFactoryBase<double> > factory =
    rcp (new Thyra::Ifpack2PreconditionerFactory<crs_matrix_type> ());

  *out << "Creating Ifpack2 preconditioner using factory" << endl;
  typedef Thyra::PreconditionerBase<ST> prec_type;
  RCP<prec_type> prec;
  try {
    prec = factory->createPrec ();
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": factory->createPrec() threw an "
      "exception: " << e.what () << endl;
    cerr << os.str ();
  }
  catch (...) {
    lclSuccess = 0;
    std::ostringstream os;
    os << "Proc " << myRank << ": factory->createPrec() threw an "
      "exception not a subclass of std::exception." << endl;
    cerr << os.str ();
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "factory->createPrec() threw an exception on some process."
          << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  // Make sure that prec is nonnull on all processes.
  if (prec.is_null ()) {
    lclSuccess = 0;
  }
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    if (myRank == 0) {
      *out << "The result of factory->createPrec() is null on at least one "
        "process." << endl;
    }
    *out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  Thyra::initializePrec (*factory, thyraA, prec.ptr ()); // segfault!
  
  // If this test starts failing, please check ifpack2/adapters/thyra/Thyra_Ifpack2PreconditionerFactory_def.hpp
  // and corresponding spot in stratimikos/adapters/belos/tpetra/Thyra_BelosTpetraPreconditionerFactory_def.hpp
  // See PR #11222 for most recent adjustments to these files. -JAL

  *out << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}
