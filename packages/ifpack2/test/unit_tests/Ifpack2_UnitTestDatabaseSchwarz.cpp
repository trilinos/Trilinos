// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestDatabaseSchwarz.cpp

\brief Ifpack2 Unit test for DatabaseSchwarz.
*/

#include "Ifpack2_config.h"
#include "Teuchos_UnitTestHarness.hpp"

// Xpetra / Galeri
#ifdef HAVE_IFPACK2_XPETRA
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"
#endif // HAVE_IFPACK2_XPETRA

#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_DatabaseSchwarz.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"

namespace { // (anonymous)

using Tpetra::global_size_t;
using Teuchos::RCP;
using std::endl;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2DatabaseSchwarz, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  const Scalar one = STS::one ();
  const Scalar two = STS::one () + STS::one ();

  out << "Ifpack2::DatabaseSchwarz: Test0" << endl;
  Teuchos::OSTab tab1 (out);

  global_size_t num_rows_per_proc = 5;

  out << "Creating row Map and CrsMatrix" << endl;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap);

  out << "Creating DatabaseSchwarz instance" << endl;

  Ifpack2::DatabaseSchwarz<row_matrix_type> prec (crsmatrix);
  Teuchos::ParameterList params, zlist;

  out << "Filling in ParameterList for DatabaseSchwarz" << endl;

  // GH: The test matrix has nnzPerRow = [2 3 3 3 2] (probably a 1D Laplacian),
  // so we set a small patch size of 3, which results in 2 patches that compress down to a database size of 1
  params.set ("database schwarz: patch size", static_cast<typename Ifpack2::DatabaseSchwarz<row_matrix_type>::local_ordinal_type>(3));
  params.set ("database schwarz: patch tolerance", 1e-6);
  params.set ("database schwarz: skip database", false);
  params.set ("database schwarz: print database summary", true);

  out << "Setting DatabaseSchwarz's parameters" << endl;

  TEST_NOTHROW(prec.setParameters(params));

  out << "Testing domain and range Maps of DatabaseSchwarz" << endl;

  // FIXME (mfh 26 Jul 2015) The domain and range Maps of the
  // preconditioner don't have to be the same object; they only need
  // to be the same in the sense of Tpetra::Map::isSameAs().
  //
  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();
  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();
  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  out << "Calling DatabaseSchwarz's initialize()" << endl;
  prec.initialize();

  out << "Calling DatabaseSchwarz's compute()" << endl;
  prec.compute();

  MV x (rowmap, 2), y (rowmap, 2), z (rowmap, 2);
  x.putScalar (one);

  out << "Calling DatabaseSchwarz's apply()" << endl;
  prec.apply(x, y);

  // The solution should now be full of 1/2s
  z.putScalar (one / two);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*STS::eps ());
}



#  define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
     TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2DatabaseSchwarz, Test0, Scalar, LocalOrdinal,GlobalOrdinal) 

typedef Tpetra::MultiVector<>::scalar_type default_scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type default_local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type default_global_ordinal_type;

UNIT_TEST_GROUP_SCALAR_ORDINAL(default_scalar_type, default_local_ordinal_type, default_global_ordinal_type)

} // namespace (anonymous)

