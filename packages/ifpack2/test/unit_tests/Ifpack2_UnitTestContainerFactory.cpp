// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_UnitTestContainerFactory.cpp
/// \brief Unit test for Ifpack2::Details::createContainer.
///
/// Ifpack2::Container<MatrixType> is an implementation detail of
/// Ifpack2::BlockRelaxation.  Each Container implements one or more
/// intraprocess subdomain solves, in either (block) Jacobi or (block)
/// Gauss-Seidel fashion.

#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_ContainerFactory.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Ifpack2_BandedContainer.hpp"
#include "Ifpack2_TriDiContainer.hpp"
#include "Ifpack2_Details_DenseSolver.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include <type_traits>

namespace { // (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(ContainerFactory, TestTypesAndInput, SC, LO, GO)
{
  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::endl;
  typedef Tpetra::Map<>::node_type NT;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  const SC ONE = Teuchos::ScalarTraits<SC>::one ();

  out << "Ifpack2::Details::createContainer (\"ContainerFactory\"): "
    "Test return type as a function of input" << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create test matrix A" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  // Create a nonzero diagonal matrix with a single row per process.
  // We won't actually do anything with it; we just need to give the
  // ContainerFactory a nontrivial matrix.

  const LO lclNumRows = 3;
  const GO gblNumRows = comm->getSize () * lclNumRows;
  const GO indexBase = 0;

  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<size_t> (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> colMap = rowMap;
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  // For a diagonal matrix, we can use the row Map as the column Map.
  const size_t maxNumEntPerRow = 1;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (rowMap, colMap, maxNumEntPerRow));
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    A->insertLocalValues (lclRow, tuple (lclRow), tuple (ONE));
  }
  A->fillComplete (domMap, ranMap);

  //Get importer (OK if null, won't actually be used in parallel apply())
  auto import = A->getGraph()->getImporter();

  out << "Set up input of createContainer" << endl;

  Teuchos::Array<Teuchos::Array<LO> > lclRows (1);
  lclRows[0].resize(1);
  lclRows[0][0] = 0;
  
  typedef Ifpack2::Container<row_matrix_type> container_type;
  {
    out << "Try creating a DenseContainer" << endl;
    auto c = Ifpack2::ContainerFactory<row_matrix_type>::build ("Dense", A, lclRows, import, false);
    static_assert (std::is_same<decltype (c), RCP<container_type> >::value,
                   "The return type of createContainer must be RCP<Container<row_matrix_type> >.");

    auto c_dense = Teuchos::rcp_dynamic_cast<Ifpack2::DenseContainer<row_matrix_type, SC> > (c);
    TEUCHOS_ASSERT( ! c_dense.is_null () );
  }
  {
    out << "Try creating a TriDiContainer" << endl;
    auto c = Ifpack2::ContainerFactory<row_matrix_type>::build ("TriDi", A, lclRows, import, false);
    auto c_tridi = Teuchos::rcp_dynamic_cast<Ifpack2::TriDiContainer<row_matrix_type, SC> > (c);
    TEUCHOS_ASSERT( ! c_tridi.is_null () );
  }
  {
    out << "Try creating a BandedContainer" << endl;
    auto c = Ifpack2::ContainerFactory<row_matrix_type>::build ("Banded", A, lclRows, import, false);
    auto c_band = Teuchos::rcp_dynamic_cast<Ifpack2::BandedContainer<row_matrix_type, SC> > (c);
    TEUCHOS_ASSERT( ! c_band.is_null () );
  }
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ContainerFactory, TestTypesAndInput, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

