// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Directory.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_UnitTestHarness.hpp"

/// Test Tpetra::DistributedContiguousDirectory, particularly in cases
/// where some processors have no IDs in the provided Map.
/// Exercises issues seen in #7223
///
/// Each test creates a distributed, contigous map.  
/// It then creates a DistributedContiguousDirectory.
/// It looks up all of the Map's IDs in the directory and checks that 
/// the Directory reports them on the correct processor.
/// It also looks up some IDs that are not in the Map and checks that
/// the Directory returns processor -1 for them.

namespace {

/////////////////////////////////////////////////////////////////////////////
// Utility functions used in all tests to check results

template <typename GO>
static int checkProc(GO gid, int proc, int correct) {
  // Report error if found processor is incorrect
  if (proc != correct) {
    std::cout << "\nError: ID " << gid << " on proc " << proc
              << "; should be on proc " << correct << std::endl;
    return 1;
  }
  return 0;
}

template <typename LO, typename GO>
static int checkLid(GO gid, LO lid, LO correct) {
  //  Report error if found LID is incorrec2
  if (lid != correct) {
    std::cout << "\nError: ID " << gid << " has LID " << lid
              << "; should be " << correct << std::endl;
    return 1;
  }
  return 0;
}

template <typename GO>
static int checkNotFound(GO gid, int proc) {
  // Report error if a valid processor is returned for an ID 
  // that is not in the map
  if (proc != -1) {
    std::cout << "\nError: ID " << gid << " on proc " << proc
              << "; should be unfound -1 " << std::endl;
    return 1;
  }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, Uniform, LO, GO )
{
  // Test -- three IDs on every processor
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  if (me == 0)
    std::cout << "\nTest 0:  three IDs per processor" << std::endl;
  size_t nMyIds = 3;

  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, 0, comm));
  Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

  size_t nIds = map->getGlobalNumElements();
  Teuchos::Array<GO> gids(nIds+1);
  Teuchos::Array<int> procs(nIds+1);
  Teuchos::Array<LO> lids(nIds+1);

  // Find location of all entries in Map, plus one bad entry nIds
  for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i);

  dir->getEntries(*map, gids(), procs(), lids(), true);

  int ierr = 0;
  for (size_t i = 0; i < nIds; i++) {
    ierr += checkProc(gids[i], procs[i], i/nMyIds);
    ierr += checkLid(gids[i], lids[i], LO(i%nMyIds));
  }
  ierr += checkNotFound(gids[nIds], procs[nIds]);

  TEST_EQUALITY_CONST(ierr, 0);
}

/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, Replicate7223, LO, GO )
{
  // Test -- two IDs on every even-numbered processor
  // When run on seven processors, this test exhibits the same behavior
  // as seen in #7332
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  if (me == 0) 
    std::cout << "\nTest 1:  one ID on each even-numbered processor; "
              << "when run on seven ranks, reproduces #7332"
              << std::endl;

  size_t nMyIds = (!(me % 2) ? 1 : 0);

  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, 0, comm));
  Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

  size_t nIds = map->getGlobalNumElements();
  Teuchos::Array<GO> gids(nIds+1);
  Teuchos::Array<int> procs(nIds+1);
  Teuchos::Array<LO> lids(nIds+1);

    // Find location of all entries in Map, plus one bad entry nIds
  for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i);
  
  dir->getEntries(*map, gids(), procs(), lids(), true);
  
  int ierr = 0;
  for (size_t i = 0; i < nIds; i++) {
    ierr += checkProc(gids[i], procs[i], 2*i);
    ierr += checkLid(gids[i], lids[i], LO(0));
  }
  ierr += checkNotFound(gids[nIds], procs[nIds]);

  TEST_EQUALITY_CONST(ierr, 0);
}

/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, EvenProcs, LO, GO )
{
  // Test -- two IDs on every even-numbered processor
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  if (me == 0)
    std::cout << "\nTest 2:  two IDs on each even-numbered processor"
              << std::endl;

  size_t nMyIds = (!(me % 2) ? 2 : 0);

  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, 0, comm));
  Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

  size_t nIds = map->getGlobalNumElements();
  Teuchos::Array<GO> gids(nIds+1);
  Teuchos::Array<int> procs(nIds+1);
  Teuchos::Array<LO> lids(nIds+1);

  // Find location of all entries in Map, plus one bad entry nIds
  for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i);
  
  dir->getEntries(*map, gids(), procs(), lids(), true);
  
  int ierr = 0;
  for (size_t i = 0; i < nIds; i++) {
    ierr += checkProc(gids[i], procs[i], i-i%2);
    ierr += checkLid(gids[i], lids[i], LO(i%2));
  }
  ierr += checkNotFound(gids[nIds], procs[nIds]);

  TEST_EQUALITY_CONST(ierr, 0);
}

/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, EvenProcsOffset, LO, GO )
{
  // Test -- two IDs on every even-numbered processor, starting with GID 72
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  size_t offset = 72;  // must be an even number for error check at end
  if (me == 0) 
    std::cout << "\nTest 3:  two IDs on each even-numbered processor; "
              << "min GID is " << offset << std::endl;
  size_t nMyIds = (!(me % 2) ? 2 : 0);

  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, offset, comm));
  Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

  size_t nIds = map->getGlobalNumElements();
  Teuchos::Array<GO> gids(nIds+2);
  Teuchos::Array<int> procs(nIds+2);
  Teuchos::Array<LO> lids(nIds+2);

  // Find location of all entries in Map, 
  // plus two bad entries offset-1 and nIds+offset
  for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i+offset);
  gids[nIds+1] = offset-1;

  dir->getEntries(*map, gids(), procs(), lids(), true);

  int ierr = 0;
  for (size_t i = 0; i < nIds; i++) {
    ierr += checkProc(gids[i], procs[i], i - i%2);
    ierr += checkLid(gids[i], lids[i], LO(i%2));
  }
  for (size_t i = nIds; i <= nIds+1; i++) 
    ierr += checkNotFound(gids[i], procs[i]);

  TEST_EQUALITY_CONST(ierr, 0);
}
/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, OddProcsOffset, LO, GO )
{
  // Test -- two IDs on every odd-numbered processor, starting with GID 72
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  size_t offset = 72;  // must be an even number for error check at end
  if (me == 0) 
    std::cout << "\nTest 3:  two IDs on each odd-numbered processor; "
              << "min GID is " << offset << std::endl;
  size_t nMyIds = ((me % 2) ? 2 : 0);

  auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, offset, comm));
  Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

  size_t nIds = map->getGlobalNumElements();
  Teuchos::Array<GO> gids(nIds+2);
  Teuchos::Array<int> procs(nIds+2);
  Teuchos::Array<LO> lids(nIds+2);

  // Find location of all entries in Map, 
  // plus two bad entries offset-1 and nIds+offset
  for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i+offset);
  gids[nIds+1] = offset-1;

  dir->getEntries(*map, gids(), procs(), lids(), true);

  int ierr = 0;
  for (size_t i = 0; i < nIds; i++) {
    ierr += checkProc(gids[i], procs[i], i + 1 - i%2);
    ierr += checkLid(gids[i], lids[i], LO(i%2));
  }
  for (size_t i = nIds; i <= nIds+1; i++) 
    ierr += checkNotFound(gids[i], procs[i]);

  TEST_EQUALITY_CONST(ierr, 0);
}
  
/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, OneProcOnly, LO, GO )
{
  // Test -- only one processor has IDs; cycle among the processors
  using map_t = Tpetra::Map<LO, GO>;
  using node_t = typename map_t::node_type;
  using dir_t = 
        Tpetra::Details::DistributedContiguousDirectory<LO, GO, node_t>;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int np = comm->getSize();
  int me = comm->getRank();

  int ierr = 0;
  for (int iter = 0; iter < np; iter++) {

    if (me == 0) 
      std::cout << "\nTest 4." << iter << ": all IDs on processor" << iter
                << std::endl;

    size_t nMyIds = (me == iter ? np : 0);
     
    auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nMyIds, 0, comm));
    Teuchos::RCP<const dir_t> dir = rcp(new dir_t(*map));

    size_t nIds = map->getGlobalNumElements();
    Teuchos::Array<GO> gids(nIds+1);
    Teuchos::Array<int> procs(nIds+1);
    Teuchos::Array<LO> lids(nIds+1);

    // Find location of all entries in Map, plus one bad entry nIds
    for (size_t i = 0; i <= nIds; i++) gids[i] = Teuchos::as<GO>(i);

    dir->getEntries(*map, gids(), procs(), lids(), true);

    for (size_t i = 0; i < nIds; i++) {
      ierr += checkProc(gids[i], procs[i], iter);
      ierr += checkLid(gids[i], lids[i], LO(i));
    }
    ierr += checkNotFound(gids[nIds], procs[nIds]);
  }
  TEST_EQUALITY_CONST(ierr, 0);
}
  
/////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, EmptyMap, LO, GO )
{
  // Test -- Empty directory:  no IDs on any processors
  int ierr = 0;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  if (me == 0) {
    std::cout << "\nTest 5:  Empty map; zero IDs on all processors" << std::endl;
    std::cout << "Empty maps are not considered distributed, so "
              << "DistributedContiguousDirectory does not apply." << std::endl;
  }
  TEST_EQUALITY_CONST(ierr, 0);
}

/////////////////////////////////////////////////////////////////////////////
//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Uniform, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Replicate7223, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, EvenProcs, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, EvenProcsOffset, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, OddProcsOffset, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, OneProcOnly, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, EmptyMap, LO, GO )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

}
