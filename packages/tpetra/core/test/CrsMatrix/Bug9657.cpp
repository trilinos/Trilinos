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
// ************************************************************************
// @HEADER
*/


#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"


namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug9657, InsertValues,
                                  Scalar, LO, GO, Node)
{
// Test for issue #9657
// Build a matrix; fillComplete it.
// Then resumeFill and Insert values.  
// fillComplete
// The host row pointers should change.  

  using map_t = Tpetra::Map<>;
  using matrix_t = Tpetra::CrsMatrix<Scalar>;
  using vector_t = Tpetra::Vector<Scalar>;

  auto comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  int nrows = 50;
  int maxNzPerRow = 5;

  Teuchos::RCP<const map_t> map = rcp(new map_t(nrows, 0, comm));
  matrix_t Amat(map, maxNzPerRow);

  Teuchos::Array<GO> cols(maxNzPerRow);
  Teuchos::Array<Scalar> vals(maxNzPerRow, 1.);

  RCP<ParameterList> params = parameterList ();
  params->set ("Optimize Storage", false);

  // Initialize matrix 
  for (size_t i = 0; i < map->getNodeNumElements(); i++) {
    GO gid = map->getGlobalElement(i);
    int nz = 0;
    cols[nz++] = gid;
    if (gid+1<nrows) {cols[nz++] = gid+1;}
    if (gid-1>=0)    {cols[nz++] = gid-1;}
    Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
  }
  Amat.fillComplete();

  std::cout << me << " of " << np << " After first fillComplete \n"
              << "  nrows     " << Amat.getNodeNumRows() << "\n"
              << "  nnz       " << Amat.getNodeNumEntries() << "\n"
              << "  maxPerRow " << Amat.getNodeMaxNumRowEntries() << "\n"
              << "  norm      " << Amat.getFrobeniusNorm() << "\n"
              << "  optimized " << Amat.isStorageOptimized() << "\n"
              << std::endl;

  auto oldLocalGraph = Amat.getCrsGraph()->getLocalGraphHost();
  auto oldRowPtrs = oldLocalGraph.row_map;

  Amat.resumeFill();
std::cout << "KDDKDD ONE" << std::endl;
  for (size_t i = 0; i < map->getNodeNumElements(); i++) {
    GO gid = map->getGlobalElement(i);
    int nz = 0;
    if (gid+2<nrows) {cols[nz++] = gid+2;}
    if (gid-2>=0)    {cols[nz++] = gid-2;}
std::cout << "KDDKDD " << gid << " inserting " << nz << std::endl;
    if (nz > 0) 
      Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
  }
std::cout << "KDDKDD TWO" << std::endl;
  Amat.fillComplete();
std::cout << "KDDKDD THREE" << std::endl;

  auto newLocalGraph = Amat.getCrsGraph()->getLocalGraphHost();
  auto newRowPtrs = newLocalGraph.row_map;

  // We inserted more values; rowPtrs should change
  TEST_ASSERT(newRowPtrs.data() != oldRowPtrs.data()); 

  std::cout << me << " of " << np << " After second fillComplete \n"
              << "  nrows     " << Amat.getNodeNumRows() << "\n"
              << "  nnz       " << Amat.getNodeNumEntries() << "\n"
              << "  maxPerRow " << Amat.getNodeMaxNumRowEntries() << "\n"
              << "  norm      " << Amat.getFrobeniusNorm() << "\n"
              << "  optimized " << Amat.isStorageOptimized() << "\n"
              << std::endl;

}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug9657, InsertValues, SCALAR, LO, GO, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
}
