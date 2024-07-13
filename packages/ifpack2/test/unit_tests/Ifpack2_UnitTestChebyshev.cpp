// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestChebyshev.cpp

\brief Ifpack2 Unit test for the Chebyshev template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Chebyshev.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Chebyshev, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Chebyshev<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  //for the diagonal matrix, set the following parameters to make
  //the preconditioner do an exact solve with a 1-degree polynomial:
  Teuchos::ParameterList params;
  params.set("chebyshev: degree", 1);
  Scalar lambdamin = zero;
  Scalar lambdamax = one;
  Scalar eigratio = one / 0.9;
  params.set("chebyshev: min eigenvalue", lambdamin);
  params.set("chebyshev: max eigenvalue", lambdamax);
  params.set("chebyshev: ratio eigenvalue", eigratio);

  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > mtx_dom_map_ptr = crsmatrix->getDomainMap();
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > mtx_rng_map_ptr = crsmatrix->getRangeMap();

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > prec_dom_map_ptr = prec.getDomainMap();
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > prec_rng_map_ptr = prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  prec.applyMat(x, y);

  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

    //Since crsmatrix is a diagonal matrix with 2 on the diagonal,
    //y should be full of 2's now.

    Teuchos::ArrayRCP<Scalar> twos(num_rows_per_proc*2, 2);

    TEST_COMPARE_FLOATING_ARRAYS(yview, twos(), Teuchos::ScalarTraits<Scalar>::eps());
  }

  prec.apply(x, y);

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType trial_tol = 1.e-13;
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol = std::max(trial_tol, Teuchos::ScalarTraits<Scalar>::eps());

  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), tol);
  }

  //If I now increase the degree of the polynomial to 4 the solve won't be
  //exact, but it should still be within a tol of 1.e-4 for this trivial data.
  params.set("chebyshev: degree", 4);
  prec.setParameters(params);
  prec.apply(x, y);

  tol = 1.e-4;

  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), tol);
  }

  crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, -one);
  Scalar n = Teuchos::as<Scalar>(rowmap->getGlobalNumElements());
  Scalar expectedLambdaMax = one-std::cos(Teuchos::ScalarTraits<Scalar>::pi()*n/(n+1));

  prec.setMatrix(crsmatrix);

  params.set("debug", true);
  params.remove("chebyshev: max eigenvalue");
  params.remove("chebyshev: min eigenvalue");

  params.set("eigen-analysis: type", "power method");
  params.set("chebyshev: eigenvalue max iterations",30);
  prec.setParameters(params);
  prec.compute();

  TEST_FLOATING_EQUALITY(prec.getLambdaMaxForApply(),expectedLambdaMax,4.5e-2);

  params.set("eigen-analysis: type", "cg");
  params.set("chebyshev: eigenvalue max iterations",10);
  prec.setParameters(params);
  prec.compute();

  TEST_FLOATING_EQUALITY(prec.getLambdaMaxForApply(),expectedLambdaMax,4.5e-2);
}

#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Chebyshev, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

