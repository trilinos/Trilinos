// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_UnitTestAdditiveSchwarzOverlap.cpp
// \brief Unit test exercising different overlap with AdditiveSchwarz.

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_Core.hpp>


namespace {

// This example exercises the fix for Bug 5963.  It illustrates the
// "Add" combine mode with overlap 0.  This example only works with
// exactly 2 MPI processes.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(AdditiveSchwarz, AddCombineMode, ScalarType, LocalOrdinalType, GlobalOrdinalType)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::toString;
  using std::endl;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef Tpetra::Map<>::node_type node_type;

  typedef Tpetra::global_size_t GST;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vec_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  typedef Ifpack2::AdditiveSchwarz<row_matrix_type> global_solver_type;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)

  TEUCHOS_TEST_FOR_EXCEPTION(
    numProcs != 2, std::logic_error, "This test must run with exactly two MPI "
    "processes.  Please fix by changing CMakeLists.txt in this directory.");

  out << "Creating Maps" << endl;

  //
  // Each of the two processes gets two rows.
  //
  const size_t localNumRows = 2;
  const GST globalNumRows = comm->getSize () * localNumRows;
  const global_ordinal_type indexBase = 0;

  RCP<const map_type> rowMap (new map_type (globalNumRows, localNumRows, indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating and filling matrix A" << endl;

  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, 3));

  Array<scalar_type> val (3);
  Array<global_ordinal_type> ind (3);
  val[0] = as<scalar_type> (-1);
  val[1] = as<scalar_type> (2);
  val[2] = as<scalar_type> (-1);

  if (myRank == 0) {
    // Row GID 0: Insert [2 -1] into column GIDs [0 1]
    ind[1] = indexBase;
    ind[2] = indexBase + 1;
    A->insertGlobalValues (indexBase, ind.view (1, 2), val.view (1, 2));

    // Row GID 1: Insert [-1 2 -1] into column GIDs [0 1 2]
    ind[0] = indexBase;
    ind[1] = indexBase + 1;
    ind[2] = indexBase + 2;
    A->insertGlobalValues (indexBase+1, ind.view (0, 3), val.view (0, 3));
  }
  else if (myRank == 1) {
    // Row GID 2: Insert [-1 2 -1] into column GIDs [1 2 3]
    ind[0] = indexBase + 1;
    ind[1] = indexBase + 2;
    ind[2] = indexBase + 3;
    A->insertGlobalValues (indexBase+2, ind.view (0, 3), val.view (0, 3));

    // Row GID 3: Insert [-1 2] into column GIDs [2 3]
    ind[0] = indexBase + 2;
    ind[1] = indexBase + 3;
    A->insertGlobalValues (indexBase+3, ind.view (0, 2), val.view (0, 2));
  }

  out << "Calling fillComplete on A" << endl;
  A->fillComplete (domMap, ranMap);

  out << "Creating b for the test problem" << endl;
  vec_type b (ranMap);
  {
    ArrayRCP<scalar_type> b_view = b.get1dViewNonConst ();
    if (myRank == 0) {
      b_view[0] = as<scalar_type> (1);
      b_view[1] = as<scalar_type> (4);
    }
    else if (myRank == 1) {
      b_view[0] = as<scalar_type> (9);
      b_view[1] = as<scalar_type> (16);
    }
  }

  out << "Creating solver" << endl;
  global_solver_type solver (A);

  out << "Setting parameters: ADD combine mode, overlap 0" << endl;
  {
    Teuchos::ParameterList pl;
    pl.set ("schwarz: combine mode", "ADD");
    pl.set ("schwarz: overlap level", 0);
    pl.set ("inner preconditioner name", "DENSE");
    solver.setParameters (pl);
  }
  out << "Calling initialize" << endl;
  solver.initialize ();

  out << "Calling compute" << endl;
  solver.compute ();

  out << "Calling apply" << endl;
  vec_type x_add (ranMap);
  solver.apply (b, x_add);

  // Solution should be [2 3 34/3 41/3] to within a very small tolerance.
  vec_type x_add_expected (ranMap);
  {
    ArrayRCP<scalar_type> view = x_add_expected.get1dViewNonConst ();
    if (myRank == 0) {
      view[0] = as<scalar_type> (2);
      view[1] = as<scalar_type> (3);
    }
    else if (myRank == 1) {
      view[0] = as<scalar_type> (34) / as<scalar_type> (3);
      view[1] = as<scalar_type> (41) / as<scalar_type> (3);
    }
  }
  // Check computed solution against expected solution.
  {
    const magnitude_type denom = x_add_expected.norm2 ();
    x_add.update (STS::one (), x_add_expected, -STS::one ());
    const magnitude_type numer = x_add.norm2 ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      (numer/denom > 10 * STS::eps ()), std::runtime_error,
      "Relative error is " << numer/denom
      << " > 10*eps = " << 10*STS::eps() << ".");
  }
}


// This example exercises the fix for Bug 5963.  It illustrates the
// "Zero" combine mode with overlap 1.  This example only works with
// exactly 2 MPI processes.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(AdditiveSchwarz, ZeroCombineMode, ScalarType, LocalOrdinalType, GlobalOrdinalType)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::toString;
  using std::endl;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef Tpetra::Map<>::node_type node_type;

  typedef Tpetra::global_size_t GST;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vec_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  typedef Ifpack2::AdditiveSchwarz<row_matrix_type> global_solver_type;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)

  TEUCHOS_TEST_FOR_EXCEPTION(
    numProcs != 2, std::logic_error, "This test must run with exactly two MPI "
    "processes.  Please fix by changing CMakeLists.txt in this directory.");

  out << "Creating Maps" << endl;

  //
  // Each of the two processes gets two rows.
  //
  const size_t localNumRows = 2;
  const GST globalNumRows = comm->getSize () * localNumRows;
  const global_ordinal_type indexBase = 0;

  RCP<const map_type> rowMap (new map_type (globalNumRows, localNumRows, indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating and filling matrix A" << endl;

  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, 3));

  Array<scalar_type> val (3);
  Array<global_ordinal_type> ind (3);
  val[0] = as<scalar_type> (-1);
  val[1] = as<scalar_type> (2);
  val[2] = as<scalar_type> (-1);

  if (myRank == 0) {
    // Row GID 0: Insert [2 -1] into column GIDs [0 1]
    ind[1] = indexBase;
    ind[2] = indexBase + 1;
    A->insertGlobalValues (indexBase, ind.view (1, 2), val.view (1, 2));

    // Row GID 1: Insert [-1 2 -1] into column GIDs [0 1 2]
    ind[0] = indexBase;
    ind[1] = indexBase + 1;
    ind[2] = indexBase + 2;
    A->insertGlobalValues (indexBase+1, ind.view (0, 3), val.view (0, 3));
  }
  else if (myRank == 1) {
    // Row GID 2: Insert [-1 2 -1] into column GIDs [1 2 3]
    ind[0] = indexBase + 1;
    ind[1] = indexBase + 2;
    ind[2] = indexBase + 3;
    A->insertGlobalValues (indexBase+2, ind.view (0, 3), val.view (0, 3));

    // Row GID 3: Insert [-1 2] into column GIDs [2 3]
    ind[0] = indexBase + 2;
    ind[1] = indexBase + 3;
    A->insertGlobalValues (indexBase+3, ind.view (0, 2), val.view (0, 2));
  }

  out << "Calling fillComplete on A" << endl;
  A->fillComplete (domMap, ranMap);

  out << "Creating b for the test problem" << endl;
  vec_type b (ranMap);
  {
    ArrayRCP<scalar_type> b_view = b.get1dViewNonConst ();
    if (myRank == 0) {
      b_view[0] = as<scalar_type> (1);
      b_view[1] = as<scalar_type> (4);
    }
    else if (myRank == 1) {
      b_view[0] = as<scalar_type> (9);
      b_view[1] = as<scalar_type> (16);
    }
  }

  out << "Creating solver" << endl;
  global_solver_type solver (A);

  out << "Setting parameters: ADD combine mode, overlap 0" << endl;
  {
    Teuchos::ParameterList pl;
    pl.set ("schwarz: combine mode", "ADD");
    pl.set ("schwarz: overlap level", 0);
    pl.set ("inner preconditioner name", "DENSE");
    solver.setParameters (pl);
  }
  out << "Calling initialize" << endl;
  solver.initialize ();

  out << "Calling compute" << endl;
  solver.compute ();

  out << "Calling apply" << endl;
  vec_type x_add (ranMap);
  solver.apply (b, x_add);

  out << "Checking solution" << endl;

  // Solution should be [2 3 34/3 41/3] to within a very small tolerance.
  vec_type x_add_expected (ranMap);
  {
    ArrayRCP<scalar_type> view = x_add_expected.get1dViewNonConst ();
    if (myRank == 0) {
      view[0] = as<scalar_type> (2);
      view[1] = as<scalar_type> (3);
    }
    else if (myRank == 1) {
      view[0] = as<scalar_type> (34) / as<scalar_type> (3);
      view[1] = as<scalar_type> (41) / as<scalar_type> (3);
    }
  }
  // Check computed solution against expected solution.
  {
    const magnitude_type denom = x_add_expected.norm2 ();
    x_add.update (STS::one (), x_add_expected, -STS::one ());
    const magnitude_type numer = x_add.norm2 ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      (numer/denom > 10 * STS::eps ()), std::runtime_error,
      "Relative error is " << numer/denom
      << " > 10*eps = " << 10*STS::eps() << ".");
  }

  out << "Setting parameters: ZERO combine mode, overlap 1" << endl;
  {
    Teuchos::ParameterList pl;
    pl.set ("schwarz: combine mode", "ZERO");
    pl.set ("schwarz: overlap level", 1);
    pl.set ("inner preconditioner name", "DENSE");
    solver.setParameters (pl);
  }
  out << "Calling initialize" << endl;
  solver.initialize ();

  out << "Calling compute" << endl;
  solver.compute ();

  out << "Calling apply" << endl;
  vec_type x_zero (ranMap);
  solver.apply (b, x_zero);

  out << "Checking solution" << endl;

  // Solution should be [5 9 19 35/2] to within a very small tolerance.
  vec_type x_zero_expected (ranMap);
  {
    ArrayRCP<scalar_type> view = x_zero_expected.get1dViewNonConst ();
    if (myRank == 0) {
      view[0] = as<scalar_type> (5);
      view[1] = as<scalar_type> (9);
    }
    else if (myRank == 1) {
      view[0] = as<scalar_type> (19);
      view[1] = as<scalar_type> (35) / as<scalar_type> (2);
    }
  }
  // Check computed solution against expected solution.
  {
    const magnitude_type denom = x_zero_expected.norm2 ();
    x_zero.update (STS::one (), x_zero_expected, -STS::one ());
    const magnitude_type numer = x_zero.norm2 ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      (numer/denom > 10 * STS::eps ()), std::runtime_error,
      "Relative error is " << numer/denom
      << " > 10*eps = " << 10*STS::eps() << ".");
  }
}


// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AdditiveSchwarz, ZeroCombineMode, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

