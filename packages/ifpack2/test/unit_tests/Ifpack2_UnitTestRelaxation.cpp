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


/*! \file Ifpack2_UnitTestRelaxation.cpp

\brief Ifpack2 Unit test for the Relaxation template.
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
#include <Ifpack2_Relaxation.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Relaxation<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");

  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();

  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>* prec_rng_map_ptr = &*prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr);
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  prec.applyMat(x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //Since crsmatrix is a diagonal matrix with 2 on the diagonal,
  //y should be full of 2's now.

  Teuchos::ArrayRCP<Scalar> twos(num_rows_per_proc*2, 2);

  TEST_COMPARE_FLOATING_ARRAYS(yview, twos(), Teuchos::ScalarTraits<Scalar>::eps());

  prec.apply(x, y);

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
}

// Test apply() with x == y.
// When x == y, apply() need to create internally an auxiliary vector, Xcopy.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Relaxation<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2);
  x.putScalar(1);

  prec.apply(x, x);

  //y should be full of 0.5's now.
  {
    Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    TEST_COMPARE_FLOATING_ARRAYS(xview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

// Test apply() with x != y  but x and y pointing to the same memory location.
// Tpetra Multivector public constructors are always copying input data so it is harder to reach such case than with Ifpack/Epetra.
// Nevertheless, it is still possible to create two different Tpetra vectors pointing to the same memory location by using MultiVector::subView().
// I can't imagine anyone trying to do this but... in this case, apply() need also to create internally an auxiliary vector, Xcopy.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Relaxation<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> y(rowmap,2);
  y.putScalar(1);
  Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xrcp = y.subView(Teuchos::Range1D(0,1));
  const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & x = *xrcp;

  TEST_INEQUALITY(&x, &y);                                               // vector x and y are different
  TEST_EQUALITY(&(x.get2dView()[0][0]), &(y.get2dView()[0][0]));         // vector x and y are pointing to the same memory location (such test only works if num of local elements != 0)
  TEST_EQUALITY(x.getLocalMV().getValues(), y.getLocalMV().getValues()); // another way to test if x and y are pointing to the same memory location

  prec.apply(x, y);

  //y should be full of 0.5's now.
  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

  // Test apply() with a partially "null" x and y. In parallel, it is possible that some nodes do not have any local elements.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test3, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 0;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  if(!comm->getRank()) num_rows_per_proc=5;


  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Relaxation<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), num_rows_per_proc);

  TEST_NOTHROW(x.getLocalMV().getValues());
  TEST_NOTHROW(y.getLocalMV().getValues());

  TEST_NOTHROW(prec.apply(x, y));
}

  // Test apply() to make sure the L1 methods work
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test4, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Relaxation<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  params.set("relaxation: use l1",true);
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), 5);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), 5);

  TEST_NOTHROW(x.getLocalMV().getValues());
  TEST_NOTHROW(y.getLocalMV().getValues());

  TEST_NOTHROW(prec.apply(x, y));
}

// Check that (symmetric) Gauss-Seidel works if there are some MPI processes with zero rows of the matrix.
// Test contributed by Jonathan Hu on 04 Jan 2013.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, SymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << endl;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  if (comm->getSize () == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  // Create the Kokkos Node instance.  This ensures that the test will
  // work even if Node is not the default Kokkos Node type.
  RCP<Node> node;
  { // All Node constructors demand a ParameterList input.
    ParameterList junk;
    node = rcp (new Node (junk));
  }

  // The number of rows of the matrix and vectors owned by the calling
  // process.  One process (Proc 0) owns zero rows, and the other
  // process(es) own a nonzero amount of rows.  This is the point of
  // the test -- to ensure that (symmetric) Gauss-Seidel works if some
  // (but not all) processes have zero rows.
  LO nelPerProc;
  if (comm->getRank () == 0) {
    nelPerProc = 0;
  }
  else {
    nelPerProc = 100;
  }

  const global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
  RCP<const map_type> rowMap (new map_type (INVALID, nelPerProc, 0, comm, node));
  RCP<const crs_matrix_type> crsMatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowMap);

  // We don't really need prec to be a pointer here, but it's
  // convenient for putting the constructor call in a TEST_NOTHROW
  // macro.  Otherwise the declaration of prec would be in an inner
  // scope and we couldn't use it below.
  RCP<Ifpack2::Relaxation<crs_matrix_type> > prec;
  {
    TEST_NOTHROW( prec = rcp (new Ifpack2::Relaxation<crs_matrix_type> (crsMatrix)) );
  }

  ParameterList params;
  params.set ("relaxation: type", "Symmetric Gauss-Seidel");
  prec->setParameters (params);
  TEST_NOTHROW( prec->initialize () );
  TEST_NOTHROW( prec->compute () );

  mv_type X (rowMap, 2);
  X.putScalar (STS::zero ());
  mv_type B (rowMap, 2);
  B.putScalar (STS::one ());

  prec->apply (B, X);
}


#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test1, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test2, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test3, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test4, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, SymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>


