// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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


/*! \file Ifpack2_UnitTestILUT.cpp

\brief Ifpack2 Unit test for the ILUT template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>
#include <Ifpack2_ILUT.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  int overlapLevel=0;
  Ifpack2::AdditiveSchwarz<CrsType,Ifpack2::ILUT<CrsType > > prec(crsmatrix,overlapLevel);
  Teuchos::ParameterList params, zlist;
  zlist.set("order_method","rcm"); 
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: drop tolerance", 0.0);
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set("schwarz: use reordering",true);
  params.set("schwarz: reordering list",zlist);
#else
  params.set("schwarz: use reordering",false);
#endif

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

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2), z(rowmap,2);
  x.putScalar(1);
  prec.apply(x, y);

  // The solution should now be full of 1/2s 
  z.putScalar(0.5);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}


//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) return;

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  int overlapLevel=3;
  Ifpack2::AdditiveSchwarz<CrsType,Ifpack2::ILUT<CrsType > > prec(crsmatrix,overlapLevel);
  Teuchos::ParameterList params, zlist;
  zlist.set("order_method","rcm"); 
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: drop tolerance", 0.0);
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  params.set("schwarz: use reordering",true);
  params.set("schwarz: reordering list",zlist);
#else
  params.set("schwarz: use reordering",false);
#endif

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

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
  x.putScalar(1);
  prec.apply(x, y);

  // The solution should now be full of 1/2s 
  z.putScalar(0.5);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
}



//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsType;
  
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  // Don't run in serial.
  if(rowmap->getComm()->getSize()==1) return;

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix2<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  int overlapLevel=0;
  Ifpack2::AdditiveSchwarz<CrsType,Ifpack2::ILUT<CrsType > > prec(crsmatrix,overlapLevel);
  Teuchos::ParameterList params, zlist;
  params.set("fact: ilut level-of-fill", 6.0);
  params.set("fact: drop tolerance", 0.0);
  params.set("schwarz: use reordering",false);
  TEST_NOTHROW(prec.setParameters(params));


  prec.initialize();
  prec.compute();


  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,1), y(rowmap,1), z(rowmap,1);
  x.putScalar(1);
  crsmatrix->apply(x,z);
  prec.apply(z, y);

  // The solution should now be full of 1s 
  z.putScalar(1.0);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

  TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());



}



#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test0, Scalar, LocalOrdinal,GlobalOrdinal)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test1, Scalar, LocalOrdinal,GlobalOrdinal)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2AdditiveSchwarz, Test2, Scalar, LocalOrdinal,GlobalOrdinal)  


UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

