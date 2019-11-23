/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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


/*! \file Ifpack2_UnitTestRILUK.cpp

\brief Ifpack2 single-process unit tests for the RILUK template.
*/


#include "Teuchos_ConfigDefs.hpp"
#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <iostream>

#include "Tpetra_Core.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_RILUK.hpp"

namespace { // (anonymous)

using Tpetra::global_size_t;
typedef tif_utest::Node Node;


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RILUKSingleProcess, FillLevel, Scalar, LO, GO)
{
  // Test that ILU(k) computes correct factors in serial for fill levels 0 to 5.
  // This test does nothing in parallel.
  // 1) Read matrix A known to exhibit need for resizing with StaticProfile
  // 2) Compute ILU(lof) of A.
  // 3) Compare to ILU(lof) created with DynamicProfile

  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node>                multivector_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>                  crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>                  row_matrix_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type>         reader_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef Tpetra::Map<LO,GO,Node>                               map_type;
  typedef Teuchos::ScalarTraits<Scalar>                         TST;

  out << "Ifpack2::RILUK: FillLevel" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  if (comm->getSize() > 1) {
    out << endl << "This test is only meaningful in serial." << endl;
    return;
  }

  for (GO lof=0; lof<6; ++lof) {

    std::string aFile = "6309_input.mm";
    RCP<const crs_matrix_type> crsmatrix = reader_type::readSparseFile (aFile, comm);
    Ifpack2::RILUK<row_matrix_type> prec(crsmatrix);

    Teuchos::ParameterList params;
    params.set("fact: iluk level-of-fill", lof);
    params.set("fact: iluk level-of-overlap", 0);

    const double overalloc = 1.+(1./(lof+1)); // non-default only to test resize
    params.set("fact: iluk overalloc", overalloc);

    prec.setParameters(params);
    prec.initialize();
    prec.compute();
    //extract incomplete factors
    const crs_matrix_type &iL = prec.getL();
    const crs_matrix_type &iU = prec.getU();
    const multivector_type &iD = prec.getD();

#ifdef INCLUDE_TEST
    //read L,U, and D factors from file
    std::string lFile = "Lfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string uFile = "Ufactor_bw" + Teuchos::toString(lof+2) + ".mm";
    std::string dFile = "Dfactor_bw" + Teuchos::toString(lof+2) + ".mm";
    out << "reading " << lFile << ", " << uFile << ", " << dFile << std::endl;
    RCP<crs_matrix_type> L = reader_type::readSparseFile (lFile, comm);
    RCP<crs_matrix_type> U = reader_type::readSparseFile (uFile, comm);
    RCP<const map_type> rm = U->getRowMap();
    RCP<multivector_type> D = reader_type::readVectorFile (dFile, comm, rm);

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
#endif
  } //for (GO lof=0; lof<6; ++lof)
} // unit test FillLevel()


//
// Instantiate and run unit tests
//

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RILUKSingleProcess, FillLevel, SC, LO, GO )

// FIXME (21 Oct 2015) There was a FIXME here a while back about
// matrix-matrix add not getting instantiated for Scalar != double.
// Need to fix that to make this test work for Scalar != double.

#ifdef HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
     UNIT_TEST_GROUP_SC_LO_GO( double, LO, GO )
#else // NOT HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_LO_GO( LO, GO )
#endif // HAVE_TPETRA_INST_DOUBLE

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

IFPACK2_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

} // namespace (anonymous)

