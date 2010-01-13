// ***********************************************************************
// 
//      Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Tifpack_UnitTestChebyshev.cpp

\brief Tifpack Unit test for the Chebyshev template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tifpack_Version.hpp>
#include <iostream>

#ifdef HAVE_TIFPACK_QD
#include <qd/dd_real.h>
#endif

#include <Tifpack_UnitTestHelpers.hpp>
#include <Tifpack_Chebyshev.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(TifpackChebyshev, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Tifpack::Version();
  out << "Tifpack::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Tifpack::Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node> prec(crsmatrix);

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

  TEUCHOS_TEST_NOTHROW(prec.SetParameters(params), out, success);

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();

  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();

  TEUCHOS_TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr, out, success );
  TEUCHOS_TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr, out, success );

  prec.Initialize();
  prec.Compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  prec.apply(x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //Since crsmatrix is a diagonal matrix with 2 on the diagonal,
  //y should be full of 2's now.

  Teuchos::ArrayRCP<Scalar> twos(num_rows_per_proc*2, 2);

  TEST_COMPARE_FLOATING_ARRAYS(yview, twos(), Teuchos::ScalarTraits<Scalar>::eps());

  prec.applyInverse(x, y);

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType trial_tol = 1.e-13;
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol = std::max(trial_tol, Teuchos::ScalarTraits<Scalar>::eps());
  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), tol);

  //If I now increase the degree of the polynomial to 4 the solve won't be
  //exact, but it should still be within a tol of 1.e-4 for this trivial data.
  params.set("chebyshev: degree", 4);
  prec.SetParameters(params);
  prec.applyInverse(x, y);

  tol = 1.e-4;
  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), tol);
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TifpackChebyshev, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
typedef long int LongInt;
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, int, LongInt)

#ifdef HAVE_TIFPACK_QD
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

