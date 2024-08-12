// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Program to testing Zoltan2::findUniqueGids capability
// Input:  Vector of keys:  each key is an array with N entries
//         Result vector to be filled by findUniqueGids
// Output: Filled result vector


#include <iostream>
#include <vector>
#include <array>
#include <unordered_set>
#include <string>
#include <typeinfo>

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Directory_Impl.hpp>

namespace Zoltan2
{

template <typename key_t, typename gno_t>
size_t findUniqueGids(
  const std::vector<key_t> &keys,
  std::vector<gno_t> &gids,
  Teuchos::RCP<const Teuchos::Comm<int> > &comm
)
{
  // Compute the new GIDs
  const bool bUseLocalIDs = false;  // Local IDs not needed
  int debug_level = 0;

  typedef int lno_t; // unused

  typedef Zoltan2_Directory_Simple<key_t,lno_t,gno_t> directory_t;

  directory_t directory(comm, bUseLocalIDs, debug_level);

  directory.update(keys.size(), &keys[0], NULL, &gids[0], NULL,
    directory_t::Update_Mode::Replace);

  directory.remap_user_data_as_unique_gids();

  // Retrieve the global numbers and put in the result gids vector
  directory.find(keys.size(), &keys[0], NULL, &gids[0], NULL, NULL, false);

  // using ssize_t can be long* on clang and I get warnings here or in the original
  // findUniqueGids - using long long specifically is ok. Can we do an MPI type
  // as size_t safely or is a conversion necessary? TODO but this gives a clean
  // build result.
  typedef long long mpi_t;
  mpi_t nDDEntries = static_cast<mpi_t>(directory.node_map_size());
  mpi_t nUnique = 0;

  // TODO use Teuchos
#ifdef HAVE_MPI
  MPI_Allreduce(&nDDEntries, &nUnique, 1, MPI_LONG_LONG, MPI_SUM,
    Teuchos::getRawMpiComm(*comm));
#else
  MPI_Allreduce(&nDDEntries, &nUnique, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  return size_t(nUnique);
}

template<typename T>
struct type_name
{
  static const char* name() {
    std::cout << "You are missing a DECL_TYPE_NAME" << std::endl;
    return(NULL);
  }
};

#define DECL_TYPE_NAME(x) \
  template<> struct type_name<x> { static const char* name() {return #x;} }

DECL_TYPE_NAME(int);
DECL_TYPE_NAME(long long);

///////////////////////////////////////////////////////////////////////////
// Tests for correctness
static const std::string fail = "FAIL ";
static const std::string pass = "     ";

// Test for correct number of unique Gids
void checkNUnique(std::string &name, size_t nUniqueGids, size_t nExpected)
{
  if (nUniqueGids != nExpected)
    std::cout << fail << name
              << "nUniqueGids " << nUniqueGids << " != " << nExpected
              << std::endl;
}

// Test for correct maximum Gid
template <typename gno_t>
void checkMaxGid(
  std::string &name,
  std::vector<gno_t> &gids,
  gno_t maxExpected,
  Teuchos::RCP<const Teuchos::Comm<int> > &comm
)
{
  gno_t maxGid = 0, gmaxGid = 0;
  size_t len = gids.size();
  for (size_t i = 0; i < len; i++)
    if (gids[i] > maxGid) maxGid = gids[i];

  Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MAX, 1,
                                 &maxGid, &gmaxGid);
  if (gmaxGid != maxExpected)
    std::cout << fail << name
              << "max Gid " << gmaxGid << " != " << maxExpected
              << std::endl;
}

// Test for correct minimum Gid
template <typename gno_t>
void checkMinGid(
  std::string &name,
  std::vector<gno_t> &gids,
  gno_t minExpected,
  Teuchos::RCP<const Teuchos::Comm<int> > &comm
)
{
  gno_t minGid = std::numeric_limits<gno_t>::max(), gminGid;
  size_t len = gids.size();
  for (size_t i = 0; i < len; i++)
    if (gids[i] < minGid) minGid = gids[i];

  Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN, 1,
                                    &minGid, &gminGid);
  if (gminGid != minExpected)
    std::cout << fail << name
              << "min Gid " << gminGid << " != " << minExpected
              << std::endl;
}

// Test for number of locally unique Gids
template <typename gno_t>
void checkNLocallyUnique(
  std::string &name,
  std::vector<gno_t> &gids,
  size_t nExpected)
{
  size_t gidsLen = gids.size();
  std::unordered_set<gno_t> gidsSet(gidsLen);

  size_t nDups = 0;
  for (size_t i = 0; i < gidsLen; i++) {
    if (gidsSet.find(gids[i]) != gidsSet.end()) {
      // Gid is already found locally
      nDups++;
    }
    else
      gidsSet.insert(gids[i]);
  }
  size_t nUnique = gidsLen - nDups;
  if (nUnique != nExpected)
    std::cout << fail << name
              << "num locally unique Gids " << nUnique << " != " << nExpected
              << std::endl;
}

///////////////////////////////////////////////////////////////////////////

template <typename gno_t>
void test1(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 1:
  // Key has only one entry
  // Each proc has me+1 keys
  // Keys are in range [1,np]
  int me = comm->getRank();
  int np = comm->getSize();

  std::string name = std::string(" test1: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef std::array<gno_t, 1> zkey_t;
  typedef std::vector<zkey_t> keyvec_t;
  typedef std::vector<gno_t> gidvec_t;

  const size_t nKeys = me+1;
  keyvec_t keys(nKeys);
  gidvec_t gids(nKeys);

  for (size_t i = 0; i < nKeys; i++) {
    zkey_t k;
    k[0] = i+1;
    keys[i] = k;
  }

  size_t nUniqueGids = findUniqueGids<zkey_t, gno_t>(keys,gids,comm);

  // Test for correctness
  if (me == 0)
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(np));

  checkMaxGid(name, gids, gno_t(np-1), comm);

  checkMinGid(name, gids, gno_t(0), comm);

  checkNLocallyUnique(name, gids, nKeys);
}

///////////////////////////////////////////////////////////////////////////

template <typename gno_t>
void test2(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 2:
  // Key has two entries
  // Each proc has six keys
  // Three Keys are {rank, x} for x in {1, 2, 3}
  // Three Keys are {(rank+x)%np, x} for x in {1, 2, 3}
  // Each rank has three unique and three non-unique keys
  int me = comm->getRank();
  int np = comm->getSize();

  std::string name = std::string(" test2: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef std::array<gno_t, 2> zkey_t;
  typedef std::vector<zkey_t> keyvec_t;
  typedef std::vector<gno_t> gidvec_t;

  const size_t nKeys = 6;
  const size_t nKeysHalf = 3;
  keyvec_t keys(nKeys);
  gidvec_t gids(nKeys);

  for (size_t i = 0; i < nKeysHalf; i++) {
    zkey_t k;
    k[0] = gno_t(me);
    k[1] = gno_t(i+1);
    keys[i] = k;
  }
  for (size_t i = 0; i < nKeysHalf; i++) {
    zkey_t k;
    k[0] = gno_t((me+i+1)%np);
    k[1] = gno_t(i+1);
    keys[i+nKeysHalf] = k;
  }

  size_t nUniqueGids = findUniqueGids<zkey_t,gno_t>(keys,gids,comm);

  // Test for correctness
  if (me == 0)
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(nKeysHalf*np));

  checkMaxGid(name, gids, gno_t(nKeysHalf*np-1), comm);

  checkMinGid(name, gids, gno_t(0), comm);
}

///////////////////////////////////////////////////////////////////////////
template <typename gno_t>
void test3(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 3:
  // Key has three entries
  // Each proc has 2*np keys
  // np Keys are {x, x, x} for x in {0, 1, ..., np-1}
  // np Keys are {rank, rank, x} for x in {0, 1, ..., np-1}
  // Each proc has one locally duplicated key
  // Each proc contributes np unique keys
  int me = comm->getRank();
  int np = comm->getSize();

  std::string name = std::string(" test3: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef std::array<gno_t, 3> zkey_t;
  typedef std::vector<zkey_t> keyvec_t;
  typedef std::vector<gno_t> gidvec_t;

  const size_t nKeys = 2*np;
  const size_t nKeysHalf = np;
  keyvec_t keys(nKeys);
  gidvec_t gids(nKeys);

  for (size_t i = 0; i < nKeysHalf; i++) {
    zkey_t k;
    k[0] = gno_t(me);
    k[1] = gno_t(me);
    k[2] = gno_t(i);
    keys[i+nKeysHalf] = k;
  }
  for (size_t i = 0; i < nKeysHalf; i++) {
    zkey_t k;
    k[0] = gno_t(i);
    k[1] = gno_t(i);
    k[2] = gno_t(i);
    keys[i] = k;
  }

  size_t nUniqueGids = findUniqueGids<zkey_t,gno_t>(keys,gids,comm);

  // for (size_t i = 0; i < nKeys; i++)
  //   std::cout << me << " Key " << i << ": "
  //             << keys[i][0] << " " << keys[i][1] << " " << keys[i][2]
  //             << " GID " << gids[i]
  //             << std::endl;

  // Test for correctness
  if (me == 0)
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(np*np));

  checkMaxGid(name, gids, gno_t(np*np-1), comm);

  checkMinGid(name, gids, gno_t(0), comm);

  checkNLocallyUnique(name, gids, size_t(nKeys-1));
}

///////////////////////////////////////////////////////////////////////////

template <typename gno_t>
void test4(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 4:
  // Key has four entries
  // Each proc has (rank+1)%2 keys; odd-numbered ranks are empty
  // Keys are all identical {0, 1, 2, 3}
  int me = comm->getRank();

  std::string name = std::string(" test4: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef std::array<gno_t, 4> zkey_t;
  typedef std::vector<zkey_t> keyvec_t;
  typedef std::vector<gno_t> gidvec_t;

  const size_t nKeys = (me+1)%2;
  keyvec_t keys(nKeys);
  gidvec_t gids(nKeys);

  for (size_t i = 0; i < nKeys; i++) {
    zkey_t k;
    k[0] = gno_t(0);
    k[1] = gno_t(1);
    k[2] = gno_t(2);
    k[3] = gno_t(3);
    keys[i] = k;
  }

  size_t nUniqueGids = findUniqueGids<zkey_t,gno_t>(keys,gids,comm);

  // Test for correctness
  if (me == 0)
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(1));

  checkMaxGid(name, gids, gno_t(0), comm);

  checkMinGid(name, gids, gno_t(0), comm);

  checkNLocallyUnique(name, gids, (nKeys ? size_t(1): size_t(0)));
}

} // namespace Zoltan2

int main(int argc, char *argv[])
{
  Tpetra::ScopeGuard tscope(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  Zoltan2::test1<int>(comm);
  Zoltan2::test2<int>(comm);
  Zoltan2::test3<int>(comm);
  Zoltan2::test4<int>(comm);

  Zoltan2::test1<long long>(comm);
  Zoltan2::test2<long long>(comm);
  Zoltan2::test3<long long>(comm);
  Zoltan2::test4<long long>(comm);

  return 0;
}
