/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "TpetraExt_TripleMatrixMultiply.hpp"

#include <algorithm>  //for replace()

namespace Tpetra
{
namespace TripleProductTesting
{
  using Tpetra::MatrixMarket::Reader;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using std::string;

//file names will have 'X' replaced by mat
//mat can be R, A, or P
template<typename MatrixType, typename MapType>
RCP<MatrixType> loadMatrix(char mat, string matFile, string domainFile, string rangeFile, string rowFile, string colFile, RCP<const Comm<int>>& comm)
{
  std::replace(matFile.begin(), matFile.end(), 'X', mat);
  std::replace(domainFile.begin(), domainFile.end(), 'X', mat);
  std::replace(rangeFile.begin(), rangeFile.end(), 'X', mat);
  std::replace(rowFile.begin(), rowFile.end(), 'X', mat);
  std::replace(colFile.begin(), colFile.end(), 'X', mat);
  RCP<const MapType> domainmap = Reader<MatrixType>::readMapFile (domainFile, comm);
  RCP<const MapType> rangemap  = Reader<MatrixType>::readMapFile (rangeFile, comm);
  RCP<const MapType> rowmap    = Reader<MatrixType>::readMapFile (rowFile, comm);
  RCP<const MapType> colmap    = Reader<MatrixType>::readMapFile (colFile, comm);
  return Reader<MatrixType>::readSparseFile (matFile, rowmap, colmap, domainmap, rangemap);
}

//These tests use the same R/A/P matrices as in MatrixMatrix_UnitTests, but
//use regular CRS multiply to get the "correct" Ac for checking
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_TripleProduct, RAP, SC, LO, GO, NT)
{
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> MatrixType;
  typedef Tpetra::Map<LO, GO, NT> MapType;
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  string matFile = "matrices/rapX.mtx";
  string domainFile = "matrices/rapX_domainmap.mtx";
  string rangeFile = "matrices/rapX_rangemap.mtx";
  string rowFile = "matrices/rapX_rowmap.mtx";
  string colFile = "matrices/rapX_colmap.mtx";
  auto R = loadMatrix<MatrixType, MapType>('R', matFile, domainFile, rangeFile, rowFile, colFile, comm);
  auto A = loadMatrix<MatrixType, MapType>('A', matFile, domainFile, rangeFile, rowFile, colFile, comm);
  auto P = loadMatrix<MatrixType, MapType>('P', matFile, domainFile, rangeFile, rowFile, colFile, comm);
  //Compute RAP using the triple product kernel
  MatrixType Ac(R->getRowMap(), 1);
  RCP<ParameterList> params = rcp(new ParameterList);
  params->set<size_t>("openmp: ltg thread max", 16);
  Tpetra::TripleMatrixMultiply::MultiplyRAP(*R, false, *A, false, *P, false, Ac, true, "RAP test", params);
  //Compute RAP using two separate (presumed correct) multiplications
  MatrixType AcGold(R->getRowMap(), 1);
  MatrixType AP(P->getRowMap(), 1);
  Tpetra::MatrixMatrix::Multiply(*A, false, *P, false, AP);
  Tpetra::MatrixMatrix::Multiply(*R, false, AP, false, AcGold);
  //Find the difference between Ac and AcGold, which should be very small
  SC plus1(1);
  SC minus1(-1);
  auto resid = Tpetra::MatrixMatrix::add(plus1, false, Ac, minus1, false, AcGold);
  TEST_COMPARE(resid->getFrobeniusNorm(), <=, 1e-7);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_TripleProduct, PTAP, SC, LO, GO, NT)
{
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> MatrixType;
  typedef Tpetra::Map<LO, GO, NT> MapType;
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  string matFile = "matrices/rapX.mtx";
  string domainFile = "matrices/rapX_domainmap.mtx";
  string rangeFile = "matrices/rapX_rangemap.mtx";
  string rowFile = "matrices/rapX_rowmap.mtx";
  string colFile = "matrices/rapX_colmap.mtx";
  auto A = loadMatrix<MatrixType, MapType>('A', matFile, domainFile, rangeFile, rowFile, colFile, comm);
  auto P = loadMatrix<MatrixType, MapType>('P', matFile, domainFile, rangeFile, rowFile, colFile, comm);
  //Compute RAP using the triple product kernel
  RCP<MapType> AcMap = rcp(new MapType(P->getGlobalNumCols(), 0, comm));
  MatrixType Ac(AcMap, 1);
  RCP<ParameterList> params = rcp(new ParameterList);
  params->set<size_t>("openmp: ltg thread max", 16);
  //because &R == &P, MultiplyRAP will call the P^T*A*P kernel
  Tpetra::TripleMatrixMultiply::MultiplyRAP(*P, true, *A, false, *P, false, Ac, true, "RAP test", params);
  //Compute RAP using two separate (presumed correct) multiplications
  MatrixType AcGold(AcMap, 1);
  MatrixType AP(P->getRowMap(), 1);
  Tpetra::MatrixMatrix::Multiply(*A, false, *P, false, AP);
  Tpetra::MatrixMatrix::Multiply(*P, true, AP, false, AcGold);
  //Find the difference between Ac and AcGold, which should be very small
  SC plus1(1);
  SC minus1(-1);
  auto resid = Tpetra::MatrixMatrix::add(plus1, false, Ac, minus1, false, AcGold);
  TEST_COMPARE(resid->getFrobeniusNorm(), <=, 1e-7);
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )			\
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_TripleProduct, RAP, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_TripleProduct, PTAP, SC, LO, GO, NT) \

  TPETRA_ETI_MANGLING_TYPEDEFS()
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )

}}  //namespace Tpetra::TripleProductTesting

