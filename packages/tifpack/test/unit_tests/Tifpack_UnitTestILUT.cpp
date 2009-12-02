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


/*! \file Tifpack_UnitTestILUT.cpp

\brief Tifpack Unit test for the ILUT template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tifpack_Version.hpp>
#include <iostream>

#include <Tifpack_UnitTestHelpers.hpp>
#include <Tifpack_ILUT.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(TifpackILUT, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Tifpack::Version();
  out << "Tifpack::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Tifpack::ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node> prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: drop tolerance", 1.0);

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

  prec.applyInverse(x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TifpackILUT, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)

}//namespace <anonymous>

