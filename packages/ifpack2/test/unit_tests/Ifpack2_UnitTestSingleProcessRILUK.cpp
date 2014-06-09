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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************


/*! \file Ifpack2_UnitTestRILUK.cpp

\brief Ifpack2 single-process unit tests for the RILUK template.
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
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
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
  MV x_copy = Tpetra::createCopy (x);
  MV y_copy = Tpetra::createCopy (y);
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be 2*y now.
  MV y_times_2 = Tpetra::createCopy (y);
  y_times_2.scale (Teuchos::as<Scalar> (2));
  TEST_COMPARE_FLOATING_ARRAYS(y_times_2.get1dView (), y_copy.get1dView (), 10*Teuchos::ScalarTraits<Scalar>::eps ());

  // Now test alpha == 0 and beta == 0.
  alpha = Teuchos::as<Scalar> (0);
  beta = Teuchos::as<Scalar> (0);
  out << "Testing apply() for alpha = " << alpha << " and beta = " << beta << endl;
  Tpetra::deep_copy (x_copy, x);
  Tpetra::deep_copy (y_copy, y);
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be zero now.
  MV y_zero (y.getMap (), 2);
  y_zero.putScalar (Teuchos::as<Scalar> (0));
  TEST_COMPARE_FLOATING_ARRAYS(y_zero.get1dView (), y_copy.get1dView (), Teuchos::ScalarTraits<Scalar>::zero ());

  // Now test alpha == 0 and beta == -1.
  alpha = Teuchos::as<Scalar> (0);
  beta = Teuchos::as<Scalar> (-1);
  out << "Testing apply() for alpha = " << alpha << " and beta = " << beta << endl;
  x_copy = Tpetra::createCopy (x);
  y_copy = Tpetra::createCopy (y);
  crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
  // y_copy should be -y now.
  MV y_neg = Tpetra::createCopy (y);
  y_neg.update (Teuchos::as<Scalar> (0), y, Teuchos::as<Scalar> (-1));
  TEST_COMPARE_FLOATING_ARRAYS(y_neg.get1dView (), y_copy.get1dView (), Teuchos::as<Scalar> (0));

  out << "Done with test" << endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, FillLevel, Scalar, LocalOrdinal, GlobalOrdinal)
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
    out << "||L - iL||_fro = " << mag << std::endl;
    TEST_EQUALITY(mag < 1e-12, true);
    out << std::endl;
    matdiff = Tpetra::MatrixMatrix::add(1.,false,*U,-1.,false,iU);
    mag = matdiff->getFrobeniusNorm();
    out << "||U - iU||_fro = " << mag << std::endl;
    TEST_EQUALITY(mag < 1e-12, true);
    out << std::endl;
    Teuchos::Array<magnitudeType> norms(1);
    D->norm2(norms);
    out << "||inverse(D) - inverse(iD)||_2 = " << norms[0] << std::endl;
    TEST_EQUALITY(norms[0] < 1e-12, true);
    out << std::endl;

  } //for (GlobalOrdinal lof=0; lof<6; ++lof)

} //unit test FillLevel()

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, IgnoreRowMapGIDs, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Test that ILU(k) ignores ordering of GIDs in the matrix rowmap.  This test is a virtual duplicate
  // of the test "FillLevel", with the exception that the row map GIDs are permuted.
  // This test is associated with bug#6033.
  //
  // This test does nothing in parallel.
  //
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
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();

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

  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  //Create a permuted row map.  The first entry is the same as the original row map,
  //the remainder are in descending order.
  Teuchos::ArrayView<const GlobalOrdinal> GIDs = rowMap->getNodeElementList();
  Teuchos::Array<GlobalOrdinal> permutedGIDs(GIDs.size());
  Teuchos::Array<GlobalOrdinal> origToPerm(GIDs.size());
  permutedGIDs[0] = GIDs[0];
  origToPerm[0] = 0;
  for (GlobalOrdinal i=1; i<GIDs.size(); ++i) {
    permutedGIDs[i] = GIDs[GIDs.size()-i];
    origToPerm[GIDs[GIDs.size()-i]] = i;
  }
  const LocalOrdinal indexBase = 0;
  Teuchos::RCP<const map_type> permRowMap = Teuchos::rcp(new map_type(INVALID, permutedGIDs(), indexBase, comm));

  for (GlobalOrdinal lof=0; lof<6; ++lof) {

    RCP<const crs_matrix_type > crsmatrix = tif_utest::create_banded_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,lof+2);

    //Copy the banded matrix into a new matrix with permuted row map GIDs.
    //This matrix will have the sparsity pattern as the original matrix.
    RCP<crs_matrix_type> permutedMatrix = Teuchos::rcp(new crs_matrix_type(permRowMap, 5));
    Teuchos::Array<GlobalOrdinal> Inds(5);
    Teuchos::Array<GlobalOrdinal> pInds(5);
    Teuchos::Array<Scalar>        Vals(5);
    Teuchos::Array<Scalar>        pVals(5);
    size_t numEntries;
    for (global_size_t i=0; i<num_rows_per_proc; ++i) {
      crsmatrix->getGlobalRowCopy(i,Inds(),Vals(),numEntries);
      pInds.resize(numEntries);
      pVals.resize(numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        pInds[j] = origToPerm[Inds[j]];
        pVals[j] = Vals[j];
      }
      permutedMatrix->insertGlobalValues(origToPerm[i],pInds(),pVals());
    }
    permutedMatrix->fillComplete();

    Ifpack2::RILUK<crs_matrix_type> prec(Teuchos::as< RCP<const crs_matrix_type> >(permutedMatrix));

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

    //read L,U, and D factors from file
    std::string lFile = "Lfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string uFile = "Ufactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string dFile = "Dfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    out << "reading " << lFile << ", " << uFile << ", " << dFile << std::endl;
    RCP<crs_matrix_type> L = reader_type::readSparseFile (lFile, comm, platform.getNode());
    RCP<crs_matrix_type> U = reader_type::readSparseFile (uFile, comm, platform.getNode());
    RCP<const map_type> rm = U->getRowMap();

    //Compare factors.  We can't use the Frobenius norm, as it uses GIDs.
    //Instead, we use the trick of multiply by the same random vector and comparing the
    //two norm of the results.  One of the random vectors is based on the original rowmap,
    //the other on the permuted row map.  Both contain the same random entries.
    multivector_type randVec(rowMap,1);
    randVec.randomize();
    multivector_type permRandVec(permRowMap,1);
    Teuchos::ArrayRCP<const Scalar> data  = randVec.getData(0);
    Teuchos::ArrayRCP<Scalar> pdata = permRandVec.getDataNonConst(0);
    for (global_size_t i=0; i<num_rows_per_proc; ++i)
      pdata[i] = data[i];
    data = pdata = Teuchos::null;

    out << "bandwidth = " << lof+2 << ", lof = " << lof << std::endl;
    multivector_type permResult(permRowMap,1);
    iL.apply(permRandVec,permResult);
    Teuchos::Array<magnitudeType> n1(1);
    permResult.norm2(n1);

    multivector_type result(rowMap,1);
    L->apply(randVec,result);
    Teuchos::Array<magnitudeType> n2(1);
    result.norm2(n2);

    out << "||L*randvec||_2 - ||iL*randvec||_2 = " << n1[0]-n2[0] << std::endl;
    TEST_EQUALITY(n1[0]-n2[0] < 1e-7, true);
    out << std::endl;

    iU.apply(permRandVec,permResult);
    permResult.norm2(n1);

    U->apply(randVec,result);
    result.norm2(n2);

    out << "||U*randvec||_2 - ||iU*randvec||_2 = " << n1[0]-n2[0] << std::endl;
    TEST_EQUALITY(n1[0]-n2[0] < 1e-7, true);
    out << std::endl;

    RCP<multivector_type> D = reader_type::readVectorFile (dFile, comm, platform.getNode(), rm);
    D->update(TST::one(),iD,-TST::one());
    Teuchos::Array<magnitudeType> norms(1);
    D->norm2(norms);
    out << "||inverse(D) - inverse(iD)||_2 = " << norms[0] << std::endl;
    TEST_EQUALITY(norms[0] < 1e-7, true);
    out << std::endl;

  } //for (GlobalOrdinal lof=0; lof<6; ++lof)


} //unit test IgnoreRowMapGIDs()

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, TestGIDConsistency, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Test that ILU(k) throws an exception if the ordering of the GIDs
  // in the row Map is not the same as the ordering of the local GIDs
  // in the column Map.  The ILU(k) setup and algorithm assumes this
  // for the moment.

  // 25April2014 JJH: The local filter appears to fix the column Map in parallel so that it's
  //                  consistently ordered with the row Map.  In otherwords, I can't get this
  //                  test to fail in parallel.  So this check is only necessary in serial.

  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::global_size_t GST;

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  if (comm->getSize () > 1) {
    out << endl << "This test only runs in serial." << endl;
    return;
  }
  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl;

  // Create matrix: 5 rows per process, 3 entries per row
  const LO indexBase = 0;
  GST numRowsPerProc = 5;
  RCP<map_type> rowMap =
    rcp (new map_type (INVALID, numRowsPerProc, indexBase, comm));

  // Create a column Map with the same GIDs at the row Map, but in
  // permuted order.  The first entry is the same as the row Map, the
  // remainder are in descending order.
  Teuchos::ArrayView<const GO> rowGIDs = rowMap->getNodeElementList ();
  Teuchos::Array<GO> colElements (rowGIDs.size ());
  colElements[0] = rowGIDs[0];
  for (GO i = 1; i < rowGIDs.size (); ++i) {
    colElements[i] = rowGIDs[rowGIDs.size () - i];
  }

  RCP<const map_type> colMap =
    rcp (new map_type (INVALID, colElements (), indexBase, comm));
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (rowMap, colMap, 3));

  // Construct a nondiagonal matrix.  It's not tridiagonal because of
  // process boundaries.
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one ();
  const Scalar two = one + one;
  Teuchos::Array<GO> col (3);
  Teuchos::Array<Scalar> val (3);
  size_t numLocalElts = rowMap->getNodeNumElements ();
  for (LO l_row = 0; static_cast<size_t> (l_row) < numLocalElts; ++l_row) {
    const GO g_row = rowMap->getGlobalElement (l_row);
    size_t i = 0;
    col[i] = g_row;
    val[i++] = two;
    if (l_row>0) {
      col[i] = rowMap->getGlobalElement (l_row - 1);
      val[i++] = -one;
    }
    if (static_cast<size_t> (l_row) < numLocalElts - 1) {
      col[i] = rowMap->getGlobalElement (l_row + 1);
      val[i++] = -one;
    }
    A->insertGlobalValues (g_row, col (0, i), val (0, i));
  }
  A->fillComplete ();

  RCP<const crs_matrix_type> constA = A;
  Ifpack2::RILUK<crs_matrix_type> prec (constA);

  Teuchos::ParameterList params;
  GO lof = 1;
  params.set ("fact: iluk level-of-fill", lof);
  params.set ("fact: iluk level-of-overlap", 0);

  prec.setParameters (params);
  TEST_THROW( prec.initialize (), std::runtime_error);

} //unit test TestGIDConsistency()

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, Test1, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, IgnoreRowMapGIDs, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, TestGIDConsistency, Scalar, LocalOrdinal, GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
//FIXME JJH 5/5/2014 Tpetra::MatrixMatrix::add is not instantiated for float, so instantiate FillLevel test separately
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, FillLevel, double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

