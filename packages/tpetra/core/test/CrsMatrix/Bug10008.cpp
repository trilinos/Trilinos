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
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>

#include <iostream>
#include <vector>

namespace { // anonymous

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug10008, SC, LO, GO, NT)
{
  // Demonstrate that importAndFillComplete fails for one-to-many 
  // transfers (entire matrix initially on one rank, redistributed to all ranks)
  auto comm = Tpetra::getDefaultComm();
  const int me = comm->getRank();
  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  typedef Tpetra::global_size_t GST;
  const GST dummy = Teuchos::OrdinalTraits<GST>::invalid();

  typedef Tpetra::Map<LO,GO,NT> map_t;
  typedef Tpetra::Import<LO,GO,NT> import_t;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NT> matrix_t;

  // Map with everything on rank zero
  const GO nRows = 25;
  const int nPerRow = 5;
  Teuchos::RCP<const map_t> mapOne =
           rcp(new map_t(dummy, (me == 0 ? nRows : 0), 0, comm));

  // Matrix using mapOne -- all entries on rank 0
  matrix_t matrixOne(mapOne, nPerRow);
  
  if (me == 0) {
    Teuchos::Array<GO> gids(nPerRow);
    Teuchos::Array<SC> values(nPerRow);
    for (int i = 0; i < nPerRow; i++) values[i] = 1.;
    for (GO i = 0; i < nRows; i++) {
      // five-point stencil
      int nnz = 0;

      gids[nnz++] = i;

      if (((i-1) >= 0) && (((i-1)%nPerRow) < (i%nPerRow)))
        gids[nnz++] = i-1; // west nbor

      if (((i+1) < nRows) && (((i+1) % nPerRow) > (i%nPerRow)))
        gids[nnz++] = i+1; // east nbor

      if ((i-nPerRow) >= 0)
        gids[nnz++] = i-nPerRow; // south nbor

      if ((i+nPerRow) < nRows) 
       gids[nnz++] = i+nPerRow; // north nbor

      matrixOne.insertGlobalValues(i, gids(0, nnz), values(0, nnz));
    } 
  }
  matrixOne.fillComplete();
  // matrixOne.describe(foo, Teuchos::VERB_EXTREME);

  // Map across all processors
  Teuchos::RCP<const map_t> mapMany = Teuchos::rcp (new map_t (nRows, 0, comm));

  // Importer from one to many 
  import_t importer(mapOne, mapMany);
  // importer.describe(foo, Teuchos::VERB_EXTREME);

  // Import to spread matrix across all processors
  Teuchos::RCP<matrix_t> matrixMany;

  // #10008 reports that this call to importAndFillComplete throws and hangs
  matrixOne.importAndFillComplete(matrixMany, importer, mapMany, mapMany);
  // matrixMany->describe(foo, Teuchos::VERB_EXTREME);
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug10008, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // anonymous
