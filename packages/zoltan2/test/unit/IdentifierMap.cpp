// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
//
// Test the IdentifierMap class.

#include <string>
#include <ostream>
#include <iostream>
#include <exception>
#include <utility>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Zoltan2_IdentifierMap.hpp>
#include <TestAdapters.hpp>         // for TEST_FAIL macros

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;

// We're testing with user global Ids that don't necessarily
// define "<<(ostream)", so we do this with traits.
// TODO: Let's add stringify to InputTraits, to be
//   used in debugging output.

template <typename T>
struct UserIdTraits{
  static std::string &stringify(T val) {return std::string("INVALID");}
};

template<>
struct UserIdTraits<std::pair<int, int> >{
  static std::string stringify(std::pair<int, int> p) {
    ostringstream oss;
    oss << "pair(" << p.first << ", " << p.second << ")";
    return oss.str();
  }
};

template<>
struct UserIdTraits<long>{
  static std::string stringify(long val) {
    ostringstream oss;
    oss << val;
    return oss.str();
  }
};

template<>
struct UserIdTraits<int>{
  static std::string stringify(int val) {
    ostringstream oss;
    oss << val;
    return oss.str();
  }
};

template <typename IDMAP>
  void testIdMap( RCP<const Comm<int> > &comm,
    IDMAP *map, bool gnosAreGids, bool gnosAreConsecutive,
    ArrayRCP<typename IDMAP::gid_t> &gids, 
    ArrayRCP<typename IDMAP::lid_t> &lids, 
    ArrayRCP<typename IDMAP::gid_t> &remoteGids,
    bool verbose)
{
  typedef typename IDMAP::lno_t LNO;
  typedef typename IDMAP::gno_t GNO;
  typedef typename IDMAP::gid_t GID;
  typedef typename IDMAP::lid_t LID;

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  int fail = 0;

  if (map->gnosAreGids() != gnosAreGids)
    fail = 1;

  TEST_FAIL_AND_THROW(*comm, fail==0, "gnosAreGids")

  if (map->gnosAreConsecutive() != gnosAreConsecutive)
    fail = 1;

  TEST_FAIL_AND_THROW(*comm, fail==0, "consecutiveGids")

  // Get Zoltan2's global numbers given user global Ids

  size_t nLocalIds = gids.size();
  Array<GNO> z2Ids(nLocalIds);

  try {
    map->gidTranslate(gids(), z2Ids(), Zoltan2::TRANSLATE_APP_TO_LIB);
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "gidTranslate")

  if (verbose){
    comm->barrier();
    if (rank == 0)
      std::cout << "Zoltan2 GNOs = User GIDs: " << gnosAreGids << std::endl;
    for (int p=0; p < nprocs; p++){
      if (p == rank){
        std::cout << "Rank " << p << " gnos: ";
        for (size_t i=0; i < nLocalIds; i++){
          std::cout << z2Ids[i] << " ";
        }
        std::cout << std::endl;
        std::cout.flush();
      }
      comm->barrier();
    }
    comm->barrier();
    if (rank == 0){
      std::cout << "MIN GNO " << map->getMinimumGlobalId();
      std::cout << ", MAX GNO " << map->getMaximumGlobalId() << std::endl;
      std::cout.flush();
    }
    comm->barrier();
  }

  // Get Zoltan2's global numbers given user local Ids

  Array<GNO> z2Ids2(nLocalIds);

  try {
    map->lidTranslate(lids(), z2Ids2(), Zoltan2::TRANSLATE_APP_TO_LIB);
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "lidTranslate")

  for (size_t i=0; i < nLocalIds; i++){
    if (z2Ids2[i] != z2Ids[i]){
       fail = 1;
       break;
    }
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "lidTranslate results")

  // Get User's global Ids give Zoltan2's global numbers

  Array<GID> userGids(nLocalIds);

  try {
    map->gidTranslate(userGids(), z2Ids(), Zoltan2::TRANSLATE_LIB_TO_APP);
  }
  catch (std::exception &e){
    fail = 1;
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "gidTranslate 2")

  for (size_t i=0; i < nLocalIds; i++){
    if (userGids[i] != gids[i]){
       fail = 1;
       break;
    }
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "gidTranslate 2 results")

  if (nprocs > 1){
    // Get Zoltan2 global number and owner of some remote User global Ids
    size_t nRemoteIds = remoteGids.size();
    Array<GNO> remoteGno(nRemoteIds);
    Array<int> remoteProc(nRemoteIds);
  
    try {
      map->gidGlobalTranslate(remoteGids(), remoteGno(), remoteProc());
    }
    catch (std::exception &e){
      fail = 1;
    }

    TEST_FAIL_AND_THROW(*comm, fail==0, "gidGLobalTranslate")
  
    if (verbose){
      comm->barrier();
      for (int p=0; p < nprocs; p++){
        if (rank == 0)
          std::cout << "Global info obtained from map:" << std::endl;
        if (p == rank){
          std::cout << "Rank " << p << std::endl;
          for (LNO i=0; i < nRemoteIds; i++){
            std::cout << "  GID: ";
            std::cout << UserIdTraits<GID>::stringify(remoteGids[i]);
            std::cout << ", GNO " << remoteGno[i];
            std::cout << ", Owner " << remoteProc[i] << std::endl;
          }
          std::cout << std::endl;
          std::cout.flush();
        }
        comm->barrier();
      }
      comm->barrier();
    }
  }
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  RCP<Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  long numLocalObjects = 10;
  long numRemoteObjects = 3;   // numRemoteObjects < numLocalObjects
  bool verbose = true;
  bool consecutiveGids=true;
  bool gnosAreGids=true;

  // Test these cases:
  // 1. GIDs are non-consecutive ordinals
  // 2. GIDs are non-consecutive ordinals, but we ask IdentifierMap to
  //    map them to consecutive IDs
  // 3. GIDs are consecutive ordinals
  // 4. GIDs are not Teuchos Ordinals

  ArrayRCP<long> gids(new long [numLocalObjects], 0, numLocalObjects, true);
  ArrayRCP<long> remoteGids(new long [numRemoteObjects], 0, 
    numRemoteObjects, true);
  ArrayRCP<std::pair<int,int> > remoteGidPairs(
    new std::pair<int,int> [numRemoteObjects], 0, numRemoteObjects, true);
  ArrayRCP<int> lids(new int[numLocalObjects], 0, numLocalObjects, true);

  using Zoltan2::IdentifierMap;

  //////////////////////////////////////////////////////////
  //  Ids are non-consecutive ordinals.

  long base = 10000 * rank;
  int fail = 0;

  for (int i=0; i < numLocalObjects; i++){
    gids[i] = base + i;   
    lids[i] = i;
  }

  typedef IdentifierMap<int, long, int, long> mapLongGids_t;

  mapLongGids_t *idMap = NULL;

  try{
    idMap = new mapLongGids_t(comm, env, gids, lids, false);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1; 
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor first case", 1);

  if (nprocs > 1){
    int remoteProc = (rank ? rank-1 : nprocs-1);
    base = remoteProc * 10000;
    for (int i=0; i < numRemoteObjects; i++)
      remoteGids[i] = base + i;
  }

  // We're not asking IdentifierMap to create consecutive
  // IDs, so Zoltan2 GNOs will be the User's GIDs, and
  // we will not have consecutive GNOs.

  testIdMap(comm, idMap, gnosAreGids, !consecutiveGids, 
    gids, lids, remoteGids, verbose);

  delete idMap;

  //////////////////////////////////////////////////////////
  //  Ids are non-consecutive ordinals.  
  //  IdentifierMap is asked to map them to consecutive.

  try{
    idMap = new mapLongGids_t(comm, env, gids, lids, true); 
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1; 
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor second case", 1);

  // Because we're asking IdentifierMap to make the Zoltan2 GNOs
  // consecutive, the GNOs will not be the same as the user GIDs.
  // And because we specifically asked for consecutive GNOs, we
  // will have consecutive global Ids.

  testIdMap(comm, idMap, !gnosAreGids, consecutiveGids, 
    gids, lids, remoteGids, verbose);

  delete idMap;

  //////////////////////////////////////////////////////////
  //  Ids are consecutive ordinals.  

  base = rank * numLocalObjects;
  for (int i=0; i < numLocalObjects; i++){
    gids[i] = base + i;   
  }

  try{
    idMap = new mapLongGids_t(comm, env, gids, lids, false); 
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1; 
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor third case", 1);

  if (nprocs > 1){
    int remoteProc = (rank ? rank-1 : nprocs-1);
    base = remoteProc * numLocalObjects;
    for (int i=0; i < numRemoteObjects; i++)
      remoteGids[i] = base + i;
  }

  // Because the User GIDs are ordinals, the Zoltan2 GNOs will be
  // the User GIDs. And since the User GIDs are already consecutive,
  // the Zoltan2 GNOs are consecutive.

  testIdMap(comm, idMap, gnosAreGids, consecutiveGids, 
    gids, lids, remoteGids, verbose);

  delete idMap;


#if 0
  // TODO - there is a bug in the IdentifierMap constructor
  //   when GIDs are std::pair<int,int>
  //////////////////////////////////////////////////////////
  //  Ids are not ordinals.  

  ArrayRCP<std::pair<int,int> > nonOrdinalGids(
     new std::pair<int,int> [numLocalObjects],
     0, numLocalObjects, true);

  for (int i=0; i < numLocalObjects; i++){
    nonOrdinalGids[i] = std::pair<int, int>(rank, i);
  }

  typedef IdentifierMap<int, std::pair<int,int>, int, long> mapPairGids_t;

  mapPairGids_t *idMap2 = NULL;

  try{
    idMap2 = new mapPairGids_t(comm, env, nonOrdinalGids, lids, false); 
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1; 
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "constructor fourth case", 1);

  if (nprocs > 1){
    int remoteProc = (rank ? rank-1 : nprocs-1);
    base = remoteProc * numLocalObjects;
    for (int i=0; i < numRemoteObjects; i++)
      remoteGidPairs[i] = std::pair<int,int>(remoteProc,i);
  }

  // Because the User's GIDs are not Teuchos Ordinals, they
  // will not be used as Zoltan2 GNOs.  When Zoltan2 creates
  // the global Ids for the problem, it creates consecutive
  // Ids that begin at 0.

  testIdMap(comm, idMap2, !gnosAreGids, consecutiveGids, 
    nonOrdinalGids, lids, remoteGidPairs, verbose);

  delete idMap2;
#endif

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

