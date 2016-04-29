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

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"

bool testReverse = true;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
    "test-reverse", "test-forward", &testReverse,
    "Test reverse or forward communication." );
}

//
// Test forward/reverse communication in reference to issue #227
//

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, ReverseCommunication, LO, GO, NT ) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO, GO, NT> Tpetra_Map;
  typedef Tpetra::Import<LO, GO, NT> Tpetra_Import;
  typedef Tpetra::Export<LO, GO, NT> Tpetra_Export;
  typedef Tpetra::Vector<double,LO,GO,NT> Tpetra_Vector;

  const Tpetra::global_size_t INVALID =
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<NT> node = Tpetra::TestingUtilities::getNode<NT>();

  // This test requires exactly 2 processors
  const int numProc = comm->getSize();
  TEUCHOS_ASSERT(numProc == 2);

  // create Maps
  const int proc = comm->getRank();
  Teuchos::Array<GO> nonoverlap_gids, overlap_gids;
  nonoverlap_gids.resize(4);
  if (proc == 0) {
    for (int i=0; i<4; ++i)
      nonoverlap_gids[i] = i;
    overlap_gids.resize(8);
    for (int i=0; i<8; ++i)
      overlap_gids[i] = i;
  }
  else {
    for (int i=0; i<4; ++i)
      nonoverlap_gids[i] = 4+i;
    overlap_gids = nonoverlap_gids;
  }
  RCP<const Tpetra_Map> nonoverlap_map =
    rcp(new Tpetra_Map(INVALID, nonoverlap_gids(), 0, comm, node));
  RCP<const Tpetra_Map> overlap_map =
    rcp(new Tpetra_Map(INVALID, overlap_gids(), 0, comm, node));

  // create import and export objects
  Tpetra_Import importer(nonoverlap_map, overlap_map);
  Tpetra_Export exporter(overlap_map, nonoverlap_map);

  // create vectors
  Tpetra_Vector nonoverlap_vector(nonoverlap_map);
  Tpetra_Vector overlap_vector(overlap_map);
  nonoverlap_vector.putScalar(1.0);
  overlap_vector.putScalar(0.0);

  if (testReverse)
    overlap_vector.doImport(nonoverlap_vector, exporter, Tpetra::REPLACE);
  else
    overlap_vector.doImport(nonoverlap_vector, importer, Tpetra::REPLACE);

  if (proc == 0) {
    Teuchos::ArrayRCP<const double> overlap_entries = overlap_vector.getData();
    for (int i=4; i<8; ++i)
      TEST_EQUALITY( overlap_entries[i], 1.0 );
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_3( LO, GO, NT )                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, ReverseCommunication, LO, GO, NT )
TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_LGN( UNIT_TEST_3 )
