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

#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   
#include <Zoltan2_findUniqueGids.hpp>

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
DECL_TYPE_NAME(long);
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
  const Teuchos::Comm<int> &comm 
)
{
  gno_t maxGid = 0, gmaxGid = 0;
  size_t len = gids.size();
  for (size_t i = 0; i < len; i++)
    if (gids[i] > maxGid) maxGid = gids[i];

  Teuchos::reduceAll<int, gno_t>(comm, Teuchos::REDUCE_MAX, 1,
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
  const Teuchos::Comm<int> &comm 
)
{
  gno_t minGid = std::numeric_limits<gno_t>::max(), gminGid;
  size_t len = gids.size();
  for (size_t i = 0; i < len; i++)
    if (gids[i] < minGid) minGid = gids[i];

  Teuchos::reduceAll<int, gno_t>(comm, Teuchos::REDUCE_MIN, 1,
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

  size_t nUniqueGids = Zoltan2::findUniqueGids<zkey_t, gno_t>(keys,gids,*comm);

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(np));

  checkMaxGid(name, gids, gno_t(np-1), *comm);

  checkMinGid(name, gids, gno_t(0), *comm);

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

  size_t nUniqueGids = Zoltan2::findUniqueGids<zkey_t,gno_t>(keys,gids,*comm);

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(nKeysHalf*np));

  checkMaxGid(name, gids, gno_t(nKeysHalf*np-1), *comm);

  checkMinGid(name, gids, gno_t(0), *comm);
  
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

  size_t nUniqueGids = Zoltan2::findUniqueGids<zkey_t,gno_t>(keys,gids,*comm);

  // for (size_t i = 0; i < nKeys; i++)
  //   std::cout << me << " Key " << i << ": " 
  //             << keys[i][0] << " " << keys[i][1] << " " << keys[i][2] 
  //             << " GID " << gids[i]
  //             << std::endl;

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(np*np));

  checkMaxGid(name, gids, gno_t(np*np-1), *comm);

  checkMinGid(name, gids, gno_t(0), *comm);
  
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

  size_t nUniqueGids = Zoltan2::findUniqueGids<zkey_t,gno_t>(keys,gids,*comm);

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(1));

  checkMaxGid(name, gids, gno_t(0), *comm);

  checkMinGid(name, gids, gno_t(0), *comm);
  
  checkNLocallyUnique(name, gids, (nKeys ? size_t(1): size_t(0)));
}

///////////////////////////////////////////////////////////////////////////

template <typename gno_t>
void test5(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 5:  Same as Test 3 but using the Tpetra Multivector interface.
  // Key has three entries
  // Each proc has 2*np keys
  // np Keys are {x, x, x} for x in {0, 1, ..., np-1}
  // np Keys are {rank, rank, x} for x in {0, 1, ..., np-1}
  // Each proc has one locally duplicated key
  // Each proc contributes np unique keys
  int me = comm->getRank();
  int np = comm->getSize();

  std::string name = std::string(" test5: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef int lno_t;

  const size_t nVecs = 3;
  const size_t nKeys = 2*np;
  const size_t nKeysHalf = np;

  Tpetra::global_size_t gNEntries =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm),
                                      true);

  Tpetra::MultiVector<gno_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<gno_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i+nKeysHalf, 0, gno_t(me));
    keys.replaceLocalValue(i+nKeysHalf, 1, gno_t(me));
    keys.replaceLocalValue(i+nKeysHalf, 2, gno_t(i));
  }
  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i, 0, gno_t(i));
    keys.replaceLocalValue(i, 1, gno_t(i));
    keys.replaceLocalValue(i, 2, gno_t(i));
  }

  size_t nUniqueGids = Zoltan2::findUniqueGids<lno_t,gno_t>(keys,gids);

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(np*np));

  Teuchos::ArrayRCP<const gno_t> gidsData = gids.getData();
  std::vector<gno_t> gidsVec(nKeys);
  for (size_t i = 0; i < nKeys; i++) gidsVec[i] = gidsData[i];

  checkMaxGid(name, gidsVec, gno_t(np*np-1), *comm);

  checkMinGid(name, gidsVec, gno_t(0), *comm);

  checkNLocallyUnique(name, gidsVec, size_t(nKeys-1));
} 

///////////////////////////////////////////////////////////////////////////

template <typename gno_t>
void test6(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Test 6:  Same as Test 4 but using the Tpetra Multivector interface.
  // Key has four entries
  // Each proc has (rank+1)%2 keys; odd-numbered ranks are empty
  // Keys are all identical {0, 1, 2, 3}
  int me = comm->getRank();

  std::string name = std::string(" test6: ")
                   + std::string(type_name<gno_t>::name());
  if (me == 0) std::cout << "--------\n  Starting " << name << std::endl;

  typedef int lno_t;

  const size_t nVecs = 4;
  const size_t nKeys = (me+1)%2;

  Tpetra::global_size_t gNEntries =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm),
                                      true);

  Tpetra::MultiVector<gno_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<gno_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeys; i++) {
    keys.replaceLocalValue(i, 0, gno_t(0));
    keys.replaceLocalValue(i, 1, gno_t(1));
    keys.replaceLocalValue(i, 2, gno_t(2));
    keys.replaceLocalValue(i, 3, gno_t(3));
  }

  size_t nUniqueGids = Zoltan2::findUniqueGids<lno_t,gno_t>(keys,gids);

  // Test for correctness
  if (me == 0) 
    std::cout << " " << name << " nUniqueGids " << nUniqueGids << std::endl;

  checkNUnique(name, nUniqueGids, size_t(1));

  Teuchos::ArrayRCP<const gno_t> gidsData = gids.getData();
  std::vector<gno_t> gidsVec(nKeys);
  for (size_t i = 0; i < nKeys; i++) gidsVec[i] = gidsData[i];

  checkMaxGid(name, gidsVec, gno_t(0), *comm);

  checkMinGid(name, gidsVec, gno_t(0), *comm);
  
  checkNLocallyUnique(name, gidsVec, (nKeys ? size_t(1): size_t(0)));
} 

///////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

#ifdef HAVE_TPETRA_INT_INT
  test1<int>(comm);
  test2<int>(comm);
  test3<int>(comm);
  test4<int>(comm);
  test5<int>(comm);
  test6<int>(comm);
#else
  if (comm->getRank() == 0) 
    std::cout << "Skipping int tests because Tpetra is not build with "
              << "GO == int" << std::endl;
#endif

#ifdef HAVE_TPETRA_INT_LONG
  test1<long>(comm);
  test2<long>(comm);
  test3<long>(comm);
  test4<long>(comm);
  test5<long>(comm);
  test6<long>(comm);
#else
  if (comm->getRank() == 0) 
    std::cout << "Skipping long tests because Tpetra is not build with "
              << "GO == long " << std::endl;
#endif

#ifdef HAVE_TPETRA_INT_LONG_LONG
  test1<long long>(comm);
  test2<long long>(comm);
  test3<long long>(comm);
  test4<long long>(comm);
  test5<long long>(comm);
  test6<long long>(comm);
#else
  if (comm->getRank() == 0) 
    std::cout << "Skipping long long tests because Tpetra is not build with "
              << "GO == long long" << std::endl;
#endif

  return 0;
}
