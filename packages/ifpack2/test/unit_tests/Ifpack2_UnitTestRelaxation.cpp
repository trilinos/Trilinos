/*
HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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

#include <Tpetra_Experimental_BlockMultiVector.hpp>

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



// Check that local (symmetric) Gauss-Seidel works if there are some MPI processes with zero rows of the matrix.
// Test contributed by Jonathan Hu on 04 Jan 2013.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, LocalSymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal)
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

  Teuchos::ArrayRCP<LO> rowIndices(crsMatrix->getNodeNumRows());
  for(size_t i=0; i < crsMatrix->getNodeNumRows(); i++)
    rowIndices[i] = Teuchos::as<LO>(i);
  params.set("relaxation: local smoothing indices",rowIndices);

  prec->setParameters (params);
  TEST_NOTHROW( prec->initialize () );
  TEST_NOTHROW( prec->compute () );

  mv_type X (rowMap, 2);
  X.putScalar (STS::zero ());
  mv_type B (rowMap, 2);
  B.putScalar (STS::one ());

  prec->apply (B, X);
}


// Test apply() on a NotCrsMatirx with a partially "null" x and y. In parallel, it is possible that some nodes do not have any local elements.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, NotCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 0;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  if(!comm->getRank()) num_rows_per_proc=5;


  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<crs_matrix_type> crsmatrix = Teuchos::rcp_const_cast<crs_matrix_type,const crs_matrix_type>(tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap));

  Teuchos::RCP<tif_utest::NotCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > notcrsmatrix = Teuchos::rcp(new tif_utest::NotCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsmatrix));

  Ifpack2::Relaxation<crs_matrix_type > prec(notcrsmatrix);

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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestDiagonalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  const int num_rows_per_proc = 5;
  const int blockSize = 3;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<block_crs_matrix_type> bcrsmatrix =
      Teuchos::rcp_const_cast<block_crs_matrix_type,const block_crs_matrix_type>(tif_utest::create_block_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph();

  Ifpack2::Relaxation<block_crs_matrix_type > prec(bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock(*crsgraph->getRowMap(),blockSize,1), yBlock(*crsgraph->getRowMap(),blockSize,1);
  MV x = xBlock.getMultiVectorView();
  MV y = yBlock.getMultiVectorView();
  x.putScalar(1);

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);

  TEST_NOTHROW(x.getLocalMV().getValues());
  TEST_NOTHROW(y.getLocalMV().getValues());

  TEST_NOTHROW(prec.apply(x, y));

  const Scalar exactSol = 0.2;

  for (int k = 0; k < num_rows_per_proc; ++k)
  {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j)
    {
      TEST_FLOATING_EQUALITY(yb[j],exactSol,1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestLowerTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<block_crs_matrix_type> bcrsmatrix =
      Teuchos::rcp_const_cast<block_crs_matrix_type,const block_crs_matrix_type>(tif_utest::create_triangular_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,true>(crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph();

  Ifpack2::Relaxation<block_crs_matrix_type > prec(bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Gauss-Seidel");
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock(*crsgraph->getRowMap(),blockSize,1), yBlock(*crsgraph->getRowMap(),blockSize,1);
  MV x = xBlock.getMultiVectorView();
  MV y = yBlock.getMultiVectorView();
  x.putScalar(1);

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);

  TEST_NOTHROW(x.getLocalMV().getValues());
  TEST_NOTHROW(y.getLocalMV().getValues());

  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.5;
  exactSol[1] = -0.25;
  exactSol[2] = 0.625;

  for (int k = 0; k < num_rows_per_proc; ++k)
  {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j)
    {
      TEST_FLOATING_EQUALITY(yb[j],exactSol[k],1e-14);
    }
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestUpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<block_crs_matrix_type> bcrsmatrix =
      Teuchos::rcp_const_cast<block_crs_matrix_type,const block_crs_matrix_type>(tif_utest::create_triangular_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,false>(crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph();

  Ifpack2::Relaxation<block_crs_matrix_type > prec(bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Symmetric Gauss-Seidel");
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock(*crsgraph->getRowMap(),blockSize,1), yBlock(*crsgraph->getRowMap(),blockSize,1);
  MV x = xBlock.getMultiVectorView();
  MV y = yBlock.getMultiVectorView();
  x.putScalar(1);

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);

  TEST_NOTHROW(x.getLocalMV().getValues());
  TEST_NOTHROW(y.getLocalMV().getValues());

  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.625;
  exactSol[1] = -0.25;
  exactSol[2] = 0.5;

  for (int k = 0; k < num_rows_per_proc; ++k)
  {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j)
    {
      TEST_FLOATING_EQUALITY(yb[j],exactSol[k],1e-14);
    }
  }

}


#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test1, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test2, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test3, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test4, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, SymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, LocalSymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, NotCrsMatrix, Scalar, LocalOrdinal,GlobalOrdinal) \

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

# define UNIT_TEST_GROUP_LGN( Scalar, LocalOrdinal, GlobalOrdinal ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestDiagonalBlockCrsMatrix, Scalar, LocalOrdinal,GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestLowerTriangularBlockCrsMatrix, Scalar, LocalOrdinal,GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestUpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal,GlobalOrdinal) \
    UNIT_TEST_GROUP_LGN(double, int, int)

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>


