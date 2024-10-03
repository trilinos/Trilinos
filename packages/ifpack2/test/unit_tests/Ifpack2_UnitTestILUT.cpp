// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestILUT.cpp

\brief Ifpack2 Unit test for the ILUT template.
*/


#include <iostream>
#include <type_traits>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Version.hpp>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_ILUT.hpp>


namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2ILUT, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::ILUT<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: drop tolerance", 0.0);

  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();

  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  Kokkos::View<STM*, Kokkos::HostSpace> norms("Initial norms", x.getNumVectors());
  Kokkos::View<STM*, Kokkos::HostSpace> lastNorms("previous norms", x.getNumVectors());

  x.norm2(norms);
  std::cout << "||x_init||=" << norms[0] << std::endl;
  y.norm2(norms);
  std::cout << "||y||=" << norms[0] << std::endl;

  prec.apply(x, y);

  y.norm2(norms);
  std::cout << "||y_final||=" << norms[0] << std::endl;

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2ILUT, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix3<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::ILUT<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: ilut level-of-fill", 6.0);
  params.set("fact: drop tolerance", 0.0);

  TEST_NOTHROW(prec.setParameters(params));

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  crsmatrix->apply(x,y);
  prec.apply(y, x);

  Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();

  //x should be full of 1's now.

  Teuchos::ArrayRCP<Scalar> ones(num_rows_per_proc*2, 1);

  TEST_COMPARE_FLOATING_ARRAYS(xview, ones(), 2*Teuchos::ScalarTraits<Scalar>::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2ILUT, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 10;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_banded_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 5);

  Ifpack2::ILUT<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: ilut level-of-fill", 0.5);
  params.set("fact: drop tolerance", 0.0);

  // Fill < 1.0 not allowed
  TEST_THROW(prec.setParameters(params), std::runtime_error);

  size_t numFill0, numFill1, numFill2;
  params.set("fact: ilut level-of-fill", 1.0);
  TEST_NOTHROW(prec.setParameters(params));
  prec.initialize();
  prec.compute();
  numFill0 = prec.getLocalNumEntries();
  TEST_EQUALITY(numFill0, crsmatrix->getLocalNumEntries());

  params.set("fact: ilut level-of-fill", 1.5);
  TEST_NOTHROW(prec.setParameters(params));
  prec.initialize();
  prec.compute();
  numFill1 = prec.getLocalNumEntries();

  params.set("fact: ilut level-of-fill", 2.0);
  TEST_NOTHROW(prec.setParameters(params));
  prec.initialize();
  prec.compute();
  numFill2 = prec.getLocalNumEntries();

  TEST_ASSERT(numFill0 < numFill1);
  TEST_ASSERT(numFill1 < numFill2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2ILUT, Test3, Scalar, LocalOrdinal, GlobalOrdinal)
{
  const global_size_t num_rows_per_proc = 10;
  auto rowMap = tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node> (num_rows_per_proc);
  auto A = tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (rowMap, 5);

  using prec_type = Ifpack2::ILUT<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >;
  using mag_type = decltype (prec_type (A).getAbsoluteThreshold ());

  // Check that setParameters does not change any parameters until it
  // checks that all parameters are valid.
  {
    prec_type prec (A);
    {
      Teuchos::ParameterList params;
      params.set ("fact: ilut level-of-fill", 1.5); // valid, not default
      params.set ("fact: absolute threshold", 0.3); // valid, not default
      params.set ("fact: relative threshold", 0.25); // valid, not default
      params.set ("fact: relax value", 0.20); // valid, not default
      params.set ("fact: drop tolerance", 0.15); // valid, not default
      TEST_NOTHROW( prec.setParameters (params) );
    }
    {
      Teuchos::ParameterList params;
      params.set ("fact: ilut level-of-fill", 0.9); // invalid
      params.set ("fact: absolute threshold", 0.25); // valid, not default
      params.set ("fact: relative threshold", 0.2); // valid, not default
      params.set ("fact: relax value", 0.15); // valid, not default
      params.set ("fact: drop tolerance", 0.1); // valid, not default

      TEST_THROW( prec.setParameters (params), std::runtime_error );
      // Fill level is double, not mag_type, because it depends on
      // LocalOrdinal, not on Scalar.
      TEST_EQUALITY( prec.getLevelOfFill (), double (1.5) ); // this is properly double
      TEST_EQUALITY( prec.getAbsoluteThreshold (), mag_type (0.3) );
      TEST_EQUALITY( prec.getRelativeThreshold (), mag_type (0.25) );
      TEST_EQUALITY( prec.getRelaxValue (), mag_type (0.20) );
      TEST_EQUALITY( prec.getDropTolerance (), mag_type (0.15) );
    }
  }

  // Repeat above test, but set parameters using mag_type, not double
  // (except for "fact: ilut level-of-fill", which is always double).
  if (! std::is_same<mag_type, double>::value) {
    prec_type prec (A);
    {
      Teuchos::ParameterList params;
      params.set ("fact: ilut level-of-fill", 1.5); // valid, not default
      params.set ("fact: absolute threshold", mag_type (0.3)); // valid, not default
      params.set ("fact: relative threshold", mag_type (0.25)); // valid, not default
      params.set ("fact: relax value", mag_type (0.20)); // valid, not default
      params.set ("fact: drop tolerance", mag_type (0.15)); // valid, not default
      TEST_NOTHROW( prec.setParameters (params) );
    }
    {
      Teuchos::ParameterList params;
      params.set ("fact: ilut level-of-fill", 0.9); // invalid
      params.set ("fact: absolute threshold", mag_type (0.25)); // valid, not default
      params.set ("fact: relative threshold", mag_type (0.2)); // valid, not default
      params.set ("fact: relax value", mag_type (0.15)); // valid, not default
      params.set ("fact: drop tolerance", mag_type (0.1)); // valid, not default

      TEST_THROW( prec.setParameters (params), std::runtime_error );
      // Fill level is double, not mag_type, because it depends on
      // LocalOrdinal, not on Scalar.
      TEST_EQUALITY( prec.getLevelOfFill (), double (1.5) );
      TEST_EQUALITY( prec.getAbsoluteThreshold (), mag_type (0.3) );
      TEST_EQUALITY( prec.getRelativeThreshold (), mag_type (0.25) );
      TEST_EQUALITY( prec.getRelaxValue (), mag_type (0.20) );
      TEST_EQUALITY( prec.getDropTolerance (), mag_type (0.15) );
    }
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2ILUT, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2ILUT, Test1, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2ILUT, Test2, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2ILUT, Test3, Scalar, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)
