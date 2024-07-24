// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"


namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug9116, LargeImport,
                                  Scalar, LO, GO, Node)
{
// Test for issue #9116
// Build a one-to-one map and a map with large overlap
// Fill one-to-one vector with known values
// Import from one-to-one vector to large overlap vector
// Check overlap vector for expected values
// The bug in #9116 should cause this test to 
// pass with TPETRA_ASSUME_GPU_AWARE_MPI=1 and 
// fail with TPETRA_ASSUME_GPU_AWARE_MPI=0.

  using map_t = Tpetra::Map<>;
  using vector_t = Tpetra::Vector<Scalar>;
  using lno_t = typename map_t::local_ordinal_type;
  using gno_t = typename map_t::global_ordinal_type;
  using import_t = Tpetra::Import<lno_t, gno_t>;

  auto comm = Tpetra::getDefaultComm();
  int np = comm->getSize();

  if (np == 1) {
    std::cout << "This test is useful only when number of MPI ranks > 1"
              << std::endl;
    std::cout << "TEST PASSED" << std::endl;
    TEST_ASSERT(true);
    return;
  }

  size_t len = 4000000;

  // One-to-one (oto) map across all processors
  Teuchos::RCP<const map_t> map_oto = rcp(new map_t(len, 0, comm));
  size_t myLen = map_oto->getLocalNumElements();

  // create one-to-one Vector with entries = GIDs
  vector_t x_oto(map_oto);
  {
    auto x_h = x_oto.getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i = 0; i < myLen; i++)
      x_h(i,0) = map_oto->getGlobalElement(i);
  }

  // Create overlap map with large overlap 
  Teuchos::Array<gno_t> elts(myLen*2);
  
  for (size_t i = 0; i < myLen; i++)
    elts[i] = map_oto->getGlobalElement(i);
  for (size_t i = myLen; i < myLen*2; i++)
    elts[i] = (elts[i-1] + 1) % len;

  // Create overlap vector
  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
  Teuchos::RCP<const map_t> map_olap = rcp(new map_t(dummy, elts(), 0, comm));
  vector_t x_olap(map_olap);

  // Import from one-to-one to overlap vector
  import_t importer(map_oto, map_olap);
  x_olap.doImport(x_oto, importer, Tpetra::INSERT);

  // Check results; each vector entry should equal its GID
  size_t ierr = 0;
  {
    auto x_h = x_olap.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (size_t i = 0; i < map_olap->getLocalNumElements(); i++) {
      if (x_h(i,0) != map_olap->getGlobalElement(i)) {
//      std::cout << comm->getRank() << " of " << np << ": x_olap[ " << i 
//                << "] " << x_h(i, 0) << " != "
//                << map_olap->getGlobalElement(i) << " expected" << std::endl;
        ierr++;
        
      }
    }
  }

  size_t gerr;
  Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (gerr) std::cout << "TEST FAILED with " << gerr << " errors" << std::endl;
  else std::cout << "TEST PASSED" << std::endl;

  TEST_ASSERT(gerr == 0);
}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug9116, LargeImport, SCALAR, LO, GO, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
}
