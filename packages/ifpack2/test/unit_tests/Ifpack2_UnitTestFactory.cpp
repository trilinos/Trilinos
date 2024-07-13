// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestFactory.cpp

\brief Ifpack2 Unit test for the Factory template, and some basic tests
for preconditioners it produces.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Factory.hpp>

#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>

#include <Ifpack2_BlockRelaxation.hpp>
#include <type_traits>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void check_precond_basics(Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& precond, Teuchos::FancyOStream& out, bool& success)
{
  TEST_EQUALITY( precond->getNumInitialize(), 0);
  TEST_EQUALITY( precond->getNumCompute(), 0);
  TEST_EQUALITY( precond->getNumApply(), 0);

  precond->initialize();
  precond->compute();

  TEST_EQUALITY( precond->getNumInitialize(), 1);
  TEST_EQUALITY( precond->getNumCompute(), 1);
}


template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void check_precond_apply(Teuchos::RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& precond, Teuchos::FancyOStream& out, bool& success)
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  MV x(precond->getDomainMap(),1);
  MV y(precond->getRangeMap(),1);

  x.putScalar(1);
  y.putScalar(0);

  TEST_NOTHROW(precond->apply(x, y));
}


//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Factory, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using std::endl;
  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> prec_type;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  Teuchos::OSTab tab0 (out);
  out << "Ifpack2::Factory Test0" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Factory factory;
  RCP<prec_type> prec_ilut = factory.create<row_matrix_type> ("ILUT", crsmatrix);
  TEST_EQUALITY(prec_ilut != Teuchos::null, true);

  check_precond_basics(prec_ilut, out, success);

  RCP<prec_type> prec_riluk = factory.create<row_matrix_type> ("RILUK", crsmatrix);
  TEST_EQUALITY(prec_riluk != Teuchos::null, true);

  check_precond_basics(prec_riluk, out, success);

  // Avoid unsupported complex types for now. Can be remove on next KokkosKernels integration
  if constexpr (std::is_arithmetic_v<Scalar>) {
    RCP<prec_type> prec_mdf = factory.create<row_matrix_type> ("MDF", crsmatrix);
    TEST_EQUALITY(prec_mdf != Teuchos::null, true);

    check_precond_basics(prec_mdf, out, success);
  }

  RCP<prec_type> prec_relax = factory.create<row_matrix_type> ("RELAXATION", crsmatrix);
  TEST_EQUALITY(prec_relax != Teuchos::null, true);

  check_precond_basics(prec_relax, out, success);

  RCP<prec_type> prec_diag = factory.create<row_matrix_type> ("DIAGONAL", crsmatrix);
  TEST_EQUALITY(prec_diag != Teuchos::null, true);

  check_precond_basics(prec_diag, out, success);

  RCP<prec_type> prec_cheby = factory.create<row_matrix_type> ("CHEBYSHEV", crsmatrix);
  TEST_EQUALITY(prec_cheby != Teuchos::null, true);

  check_precond_basics(prec_cheby, out, success);
}


//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Factory, BlockCrs, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using std::endl;
  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> prec_type;
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  // typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  // typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;


  Teuchos::OSTab tab0 (out);
  out << "Ifpack2::Factory: Test BlockCrs" << endl;
  Teuchos::OSTab tab1 (out);

  const int num_rows_per_proc = 5;
  const int blockSize = 3;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    tif_utest::create_tridiag_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  RCP<block_crs_matrix_type> bcrsmatrix =
    Teuchos::rcp_const_cast<block_crs_matrix_type, const block_crs_matrix_type> (tif_utest::create_block_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize));

  RCP<const row_matrix_type> rowmatrix = bcrsmatrix;


  Ifpack2::Factory factory;
  RCP<prec_type> prec_relax = factory.create<row_matrix_type> ("RELAXATION", rowmatrix);
  TEST_EQUALITY(prec_relax != Teuchos::null, true);
  check_precond_basics(prec_relax, out, success);
  check_precond_apply(prec_relax, out, success);

  // Basic block relaxation tests
  prec_relax = factory.create<row_matrix_type> ("BLOCKRELAXATION", rowmatrix);
  TEST_EQUALITY(prec_relax != Teuchos::null, true);
  check_precond_basics(prec_relax, out, success);
  check_precond_apply(prec_relax, out, success);

  // Block-Tridiagonal
  {
    Teuchos::ParameterList params;
    params.set("relaxation: container", "BlockTriDi");
    params.set("relaxation: type", "MT Split Jacobi");
    params.set("partitioner: type", "linear");
    params.set("partitioner: local parts", num_rows_per_proc);

    prec_relax = factory.create<row_matrix_type> ("BLOCKRELAXATION", rowmatrix);
    TEST_EQUALITY(prec_relax != Teuchos::null, true);
    prec_relax->setParameters(params);
    check_precond_basics(prec_relax, out, success);
    check_precond_apply(prec_relax, out, success);
  }


}


#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Factory, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Factory, BlockCrs, Scalar, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

