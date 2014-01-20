/*
//@HEADER
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************


/*! \file Ifpack2_UnitTestRILUK.cpp

\brief Ifpack2 Unit test for the RILUK template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_RILUK.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUK, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  if (rowmap->getComm()->getSize() > 1) {
    out << std::endl;
    out << "This test may only be run in serial." << std::endl;
    return;
  }

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::RILUK<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  int fill_level = 1;
  params.set("fact: iluk level-of-fill", fill_level);
  params.set("fact: iluk level-of-overlap", 0);

  TEST_NOTHROW(prec.setParameters(params));

  TEST_EQUALITY( prec.getLevelOfFill(), fill_level);

  prec.initialize();
  //trivial tests to insist that the preconditioner's domain/range maps are
  //the same as those of the matrix:
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& mtx_dom_map = *crsmatrix->getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& mtx_rng_map = *crsmatrix->getRangeMap();

  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& prec_dom_map = *prec.getDomainMap();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& prec_rng_map = *prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map.isSameAs(mtx_dom_map), true );
  TEST_EQUALITY( prec_rng_map.isSameAs(mtx_rng_map), true );

  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  prec.apply(x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUK, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

  out << "Testing Ifpack2::RILUK" << endl;

  const global_size_t num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node> (num_rows_per_proc);

  if (rowmap->getComm ()->getSize () > 1) {
    out << endl << "This test may only be run in serial or with a single MPI process." << endl;
    return;
  }

  out << "Creating matrix" << endl;
  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix2<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  out << "Creating preconditioner" << endl;
  Ifpack2::RILUK<crs_matrix_type> prec (crsmatrix);

  out << "Setting preconditioner's parameters" << endl;

  Teuchos::ParameterList params;
  params.set ("fact: iluk level-of-fill", 1);
  params.set ("fact: iluk level-of-overlap", 0);

  TEST_NOTHROW(prec.setParameters(params));

  out << "Calling initialize() and compute()" << endl;
  prec.initialize();
  prec.compute();

  out << "Creating test problem" << endl;
  MV x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  //apply the matrix to x, putting A*x in y
  crsmatrix->apply(x,y);

  out << "Calling preconditioner's apply()" << endl;

  //apply the preconditioner to y, putting ~A^-1*y in x
  //(this should set x back to 1's)
  prec.apply(y, x);

  out << "Checking result" << endl;

  Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();

  //x should be full of 1's now.

  Teuchos::ArrayRCP<Scalar> ones(num_rows_per_proc*2, 1);

  TEST_COMPARE_FLOATING_ARRAYS(xview, ones(), 2*Teuchos::ScalarTraits<Scalar>::eps());

  // Now test alpha != 1 and beta == 0.
  Scalar alpha = Teuchos::as<Scalar> (2);
  Scalar beta = Teuchos::as<Scalar> (0);
  out << "Testing apply() for alpha = " << alpha << " and beta = " << beta << endl;
  MV x_copy (x);
  MV y_copy (y);
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be 2*y now.
  MV y_times_2 (y);
  y_times_2.scale (Teuchos::as<Scalar> (2));
  TEST_COMPARE_FLOATING_ARRAYS(y_times_2.get1dView (), y_copy.get1dView (), 10*Teuchos::ScalarTraits<Scalar>::eps ());

  // Now test alpha == 0 and beta == 0.
  alpha = Teuchos::as<Scalar> (0);
  beta = Teuchos::as<Scalar> (0);
  out << "Testing apply() for alpha = " << alpha << " and beta = " << beta << endl;
  x_copy = x;
  y_copy = y;
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be zero now.
  MV y_zero (y.getMap (), 2);
  y_zero.putScalar (Teuchos::as<Scalar> (0));
  TEST_COMPARE_FLOATING_ARRAYS(y_zero.get1dView (), y_copy.get1dView (), Teuchos::ScalarTraits<Scalar>::zero ());

  // Now test alpha == 0 and beta == -1.
  alpha = Teuchos::as<Scalar> (0);
  beta = Teuchos::as<Scalar> (-1);
  out << "Testing apply() for alpha = " << alpha << " and beta = " << beta << endl;
  x_copy = x;
  y_copy = y;
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be -y now.
  MV y_neg (y);
  y_neg.update (Teuchos::as<Scalar> (0), y, Teuchos::as<Scalar> (-1));
  TEST_COMPARE_FLOATING_ARRAYS(y_neg.get1dView (), y_copy.get1dView (), Teuchos::as<Scalar> (0));

  out << "Done with test" << endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUK, FillLevel, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Test that ILU(k) computes correct factors in serial for fill levels 0 to 5.
  // This test does nothing in parallel.
  // 1) Create a banded matrix A with bandwidth equal to lof+2.
  //    Nonzero off-diagonals are subdiagonals +/-1 and +/-(lof+2).
  //    The matrix has 4's on the main diagonal, -1's on off-diagonals.
  // 2) Compute ILU(lof) of A.
  // 3) The incomplete factors should be equal to those of an exact LU decomposition without pivoting.
  //    Note that Ifpack2 stores the inverse of the diagonal separately and scales U.
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> multivector_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                map_type;
  typedef Teuchos::ScalarTraits<Scalar>                               TST;

  using Teuchos::RCP;

  Tpetra::DefaultPlatform::DefaultPlatformType& platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();

  if (comm->getSize() > 1) {
    out << std::endl;
    out << "This test is only meaningful in serial." << std::endl;
    return;
  }

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 10;

  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  for (GlobalOrdinal lof=0; lof<6; ++lof) {

    RCP<const crs_matrix_type > crsmatrix = tif_utest::create_banded_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap,lof+2);
    //std::string aFile = "A_bw=" + Teuchos::toString(lof+2) + ".mm";
    //RCP<crs_matrix_type> crsmatrix = reader_type::readSparseFile (aFile, comm, platform.getNode());
    //crsmatrix->describe(out,Teuchos::VERB_EXTREME);
    Ifpack2::RILUK<crs_matrix_type> prec(crsmatrix);

    Teuchos::ParameterList params;
    params.set("fact: iluk level-of-fill", lof);
    params.set("fact: iluk level-of-overlap", 0);

    prec.setParameters(params);
    prec.initialize();
    prec.compute();
    //extract incomplete factors
    const crs_matrix_type &iL = prec.getL();
    const crs_matrix_type &iU = prec.getU();
    const multivector_type &iD = prec.getD();

    ////if (lof==0)
    //{
    //  Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile("check_A.mm",crsmatrix);
    //  std::string outfile = "check_iL_bw=" + Teuchos::toString(lof+2) + ".mm";
    //  Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(outfile,rcpFromRef(iL));
    //  outfile = "check_iU_bw=" + Teuchos::toString(lof+2) + ".mm";
    //  Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(outfile,rcpFromRef(iU));
    //  outfile = "check_iD_bw=" + Teuchos::toString(lof+2) + ".mm";
    //  Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(outfile,rcpFromRef(iD));
    //}

    //read L,U, and D factors from file
    std::string lFile = "Lfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string uFile = "Ufactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string dFile = "Dfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    out << "reading " << lFile << ", " << uFile << ", " << dFile << std::endl;
    RCP<crs_matrix_type> L = reader_type::readSparseFile (lFile, comm, platform.getNode());
    RCP<crs_matrix_type> U = reader_type::readSparseFile (uFile, comm, platform.getNode());
    RCP<const map_type> rm = U->getRowMap();
    RCP<multivector_type> D = reader_type::readVectorFile (dFile, comm, platform.getNode(), rm);

    //compare factors
    out << "bandwidth = " << lof+2 << ", lof = " << lof << std::endl;
    D->update(TST::one(),iD,-TST::one());
    RCP<crs_matrix_type> matdiff = Tpetra::MatrixMatrix::add(1.,false,*L,-1.,false,iL);
    magnitudeType mag = matdiff->getFrobeniusNorm();
    TEST_EQUALITY(mag < 1e-12, true);
    out << "||L - iL||_fro = " << mag << std::endl;
    matdiff = Tpetra::MatrixMatrix::add(1.,false,*U,-1.,false,iU);
    mag = matdiff->getFrobeniusNorm();
    TEST_EQUALITY(mag < 1e-12, true);
    out << "||U - iU||_fro = " << mag << std::endl;
    Teuchos::Array<magnitudeType> norms(1);
    D->norm2(norms);
    TEST_EQUALITY(norms[0] < 1e-12, true);
    out << "||inverse(D) - inverse(iD)||_2 = " << norms[0] << std::endl;

  } //for (GlobalOrdinal lof=0; lof<6; ++lof)

} //unit test FillLevel()

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUK, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUK, Test1, Scalar, LocalOrdinal,GlobalOrdinal)

TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUK, FillLevel, double, int, int)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

