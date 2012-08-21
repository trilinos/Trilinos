// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Test the IdentifierMap class.
//
//   Test local IDs are implied, not supplied by app.

#include <Zoltan2_IdentifierMap.hpp>
#include <Zoltan2_TestHelpers.hpp>

#if 0
#include <string>
#include <ostream>
#include <iostream>
#include <exception>
#include <utility>
#endif

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>


using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;

template <typename T>
struct UserIdTraits{
  static std::string stringify(T val) {return std::string("INVALID");}
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
    ArrayRCP<typename IDMAP::gid_t> &remoteGids,
    bool verbose)
{
  typedef typename IDMAP::lno_t LNO;
  typedef typename IDMAP::gno_t GNO;
  typedef typename IDMAP::gid_t GID;

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
    fail = 2;
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

  // Get Zoltan2's global numbers given user local indices 

  Array<GNO> z2Ids2(nLocalIds);
  Array<LNO> indices(nLocalIds);

  for (LNO i=nLocalIds-1,j=0; i >= 0; i--,j++){
    indices[j] = i;
  }
   

  try {
    map->lnoTranslate(indices(), z2Ids2(), Zoltan2::TRANSLATE_APP_TO_LIB);
  }
  catch (std::exception &e){
    fail = 3;
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "lidTranslate")

  for (LNO i=nLocalIds-1, j=0; i >= 0; i--, j++){
    if (z2Ids2[j] != z2Ids[i]){
       fail = 4;
       break;
    }
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "lnoTranslate results")

  // Get User's global Ids give Zoltan2's global numbers

  Array<GID> userGids(nLocalIds);

  try {
    map->gidTranslate(userGids(), z2Ids(), Zoltan2::TRANSLATE_LIB_TO_APP);
  }
  catch (std::exception &e){
    fail = 5;
  }

  TEST_FAIL_AND_THROW(*comm, fail==0, "gidTranslate 2")

  for (size_t i=0; i < nLocalIds; i++){
    if (userGids[i] != gids[i]){
       fail = 6;
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
      fail = 7;
    }

    TEST_FAIL_AND_THROW(*comm, fail==0, "gidGLobalTranslate")
  
    if (verbose){
      comm->barrier();
      for (int p=0; p < nprocs; p++){
        if (rank == 0)
          std::cout << "Global info obtained from map:" << std::endl;
        if (p == rank){
          std::cout << "Rank " << p << std::endl;
          for (size_t i=0; i < nRemoteIds; i++){
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
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  lno_t numLocalObjects = 10;
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

  ArrayRCP<gno_t> gids(new gno_t [numLocalObjects], 0, numLocalObjects, true);
  ArrayRCP<gno_t> remoteGids(new gno_t [numRemoteObjects], 0, 
    numRemoteObjects, true);

  using Zoltan2::IdentifierMap;

  typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> UserTypes;

  //////////////////////////////////////////////////////////
  //  Ids are non-consecutive ordinals.

  gno_t base1 = 10000 * rank;
  gno_t base2 = base1 + 5000;
  int fail = 0;
  gno_t base = base1;

  for (lno_t i=0; i < numLocalObjects; i++){
    gids[i] = base + i;   
    if (i == numLocalObjects/2) base = base2;
  }

  typedef IdentifierMap<UserTypes> idmap_t;

  idmap_t *idMap = NULL;

  try{
    idMap = new idmap_t(env, comm, gids, false);
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
    gids, remoteGids, verbose);

  delete idMap;

  //////////////////////////////////////////////////////////
  //  Ids are non-consecutive ordinals.  
  //  IdentifierMap is asked to map them to consecutive.

  try{
    idMap = new idmap_t(env, comm, gids, true); 
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
    gids, remoteGids, verbose);

  delete idMap;

  //////////////////////////////////////////////////////////
  //  Ids are consecutive ordinals.  

  base = rank * numLocalObjects;
  for (lno_t i=0; i < numLocalObjects; i++){
    gids[i] = base + i;   
  }

  try{
    idMap = new idmap_t(env, comm, gids, false); 
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
    gids, remoteGids, verbose);

  delete idMap;

#if 0

  ArrayRCP<std::pair<int,int> > remoteGidPairs(
    new std::pair<int,int> [numRemoteObjects], 0, numRemoteObjects, true);

  // TODO - there is a bug in the IdentifierMap constructor
  //   when GIDs are std::pair<int,int>
  //////////////////////////////////////////////////////////
  //  Ids are not ordinals.  
  Zoltan2::BasicUserTypes<float, std::pair<int,int>, int, long> UserPairGids;

  ArrayRCP<std::pair<int,int> > nonOrdinalGids(
     new std::pair<int,int> [numLocalObjects],
     0, numLocalObjects, true);

  for (int i=0; i < numLocalObjects; i++){
    nonOrdinalGids[i] = std::pair<int, int>(rank, i);
  }

  typedef IdentifierMap<UserPairGids> mapPairGids_t;

  mapPairGids_t *idMap2 = NULL;

  try{
    idMap2 = new mapPairGids_t(env, comm, nonOrdinalGids, false); 
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
    nonOrdinalGids, remoteGidPairs, verbose);

  delete idMap2;
#endif

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

