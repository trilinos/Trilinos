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
// Input:  Multivector of keys:  each key has #vectors entries
//         Result vector to be filled by findUniqueGids
// Output: Filled result vector


#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <unordered_set>


template <typename scalar_t, typename lno_t, typename gno_t>
size_t findUniqueGids(
  Tpetra::MultiVector<scalar_t, lno_t, gno_t> &keys,
  Tpetra::Vector<scalar_t, lno_t, gno_t> &gids
)
{
  std::cout << "Placeholder function for compilation" << std::endl;
  return 0;
}

///////////////////////////////////////////////////////////////////////////
// Test for correct number of unique Gids

static const std::string fail = "FAIL ";

void checkNUnique(std::string &name, size_t nUniqueGids, size_t nExpected)
{  
  if (nUniqueGids != nExpected)
    std::cout << fail << name 
              << "nUniqueGids " << nUniqueGids << " != " << nExpected
              << std::endl;
}

// Test for correct maximum Gid
template <typename scalar_t, typename lno_t, typename gno_t>
void checkMaxGid(
  std::string &name, 
  Tpetra::Vector<scalar_t, lno_t, gno_t> &gids,
  scalar_t maxExpected)
{
  scalar_t maxGid = gids.normInf();
  if (maxGid != maxExpected)
    std::cout << fail << name 
              << "max Gid " << maxGid << " != " << maxExpected
              << std::endl;
}
  
// Test for correct minimum Gid
template <typename scalar_t, typename lno_t, typename gno_t>
void checkMinGid(
  std::string &name, 
  Tpetra::Vector<scalar_t, lno_t, gno_t> &gids,
  scalar_t minExpected)
{
  gids.scale(scalar_t(-1));
  scalar_t minGid = gids.normInf();
  if (minGid != minExpected)
    std::cout << fail << name 
              << "min Gid " << minGid << " != " << minExpected
              << std::endl;
  gids.scale(scalar_t(-1));
}

// Test for number of locally unique Gids 
template <typename scalar_t, typename lno_t, typename gno_t>
void checkNLocallyUnique(
  std::string &name, 
  Tpetra::Vector<scalar_t, lno_t, gno_t> &gids,
  size_t nExpected)
{
  size_t gidsLen = gids.getLocalLength();
  std::unordered_set<scalar_t> gidsSet(gidsLen);

  Teuchos::ArrayRCP<const scalar_t> gidsData = gids.getData();
  size_t nDups = 0;
  for (size_t i = 0; i < gidsLen; i++) {
    if (gidsSet.find(gidsData[i]) != gidsSet.end()) {
      // Gid is already found locally
      nDups++;
    }
    else 
      gidsSet.insert(gidsData[i]);
  }
  size_t nUnique = gidsLen - nDups;
  if (nUnique != nExpected) 
    std::cout << fail << name 
              << "num locally unique Gids " << nUnique << " != " << nExpected
              << std::endl;
}


///////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  int me = comm->getRank();
  int np = comm->getSize();

  {
  // Test 1:
  // Key has only one int entry
  // Each proc has me+1 keys
  // Keys are in range [1,np]

  std::string name = "test1: ";

  typedef int scalar_t;
  typedef int lno_t;
  typedef int gno_t;

  const size_t nVecs = 1;
  const size_t nKeys = me+1;

  Tpetra::global_size_t gNEntries = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm));

  Tpetra::MultiVector<scalar_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<scalar_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeys; i++)
    keys.replaceLocalValue(i, 0, i+1);

  size_t nUniqueGids = findUniqueGids<scalar_t, lno_t, gno_t>(keys, gids);

  // Test for correctness
  checkNUnique(name, nUniqueGids, size_t(np));

  checkMaxGid(name, gids, scalar_t(np-1));

  checkMinGid(name, gids, scalar_t(0));

  checkNLocallyUnique(name, gids, nKeys);
  }

  {
  // Test 2:
  // Key has two int entries
  // Each proc has six keys
  // Three Keys are {rank, x} for x in {1, 2, 3}
  // Three Keys are {(rank+x)%np, x} for x in {1, 2, 3}
  // Each rank has three unique and three non-unique keys

  std::string name = "test2: ";

  typedef int scalar_t;
  typedef int lno_t;
  typedef int gno_t;

  const size_t nVecs = 2;
  const size_t nKeys = 6;
  const size_t nKeysHalf = 3;

  Tpetra::global_size_t gNEntries = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm));

  Tpetra::MultiVector<scalar_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<scalar_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i, 0, me);
    keys.replaceLocalValue(i, 1, i+1);
  }
  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i+nKeysHalf, 0, (me+i+1)%np);
    keys.replaceLocalValue(i+nKeysHalf, 1, i+1);
  }

  size_t nUniqueGids = findUniqueGids<scalar_t, lno_t, gno_t>(keys, gids);

  // Test for correctness
  checkNUnique(name, nUniqueGids, size_t(nKeysHalf*np));

  checkMaxGid(name, gids, scalar_t(nKeysHalf*np-1));

  checkMinGid(name, gids, scalar_t(0));
  
  }

  {
  // Test 3:
  // Key has three long long entries
  // Each proc has 2*np keys
  // np Keys are {x, x, x} for x in {0, 1, ..., np-1}
  // np Keys are {rank, rank, x} for x in {0, 1, ..., np-1}
  // Each proc has one locally duplicated key
  // Each proc contributes np unique keys

  std::string name = "test3: ";

  typedef long long scalar_t;
  typedef int lno_t;
  typedef long long gno_t;

  const size_t nVecs = 3;
  const size_t nKeys = 2*np;
  const size_t nKeysHalf = np;

  Tpetra::global_size_t gNEntries = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm));

  Tpetra::MultiVector<scalar_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<scalar_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i+nKeysHalf, 0, scalar_t(me));
    keys.replaceLocalValue(i+nKeysHalf, 1, scalar_t(me));
    keys.replaceLocalValue(i+nKeysHalf, 2, scalar_t(i));
  }
  for (size_t i = 0; i < nKeysHalf; i++) {
    keys.replaceLocalValue(i, 0, scalar_t(i));
    keys.replaceLocalValue(i, 1, scalar_t(i));
    keys.replaceLocalValue(i, 2, scalar_t(i));
  }

  size_t nUniqueGids = findUniqueGids<scalar_t, lno_t, gno_t>(keys, gids);

  // Test for correctness
  checkNUnique(name, nUniqueGids, size_t(np*np));

  checkMaxGid(name, gids, scalar_t(np*np-1));

  checkMinGid(name, gids, scalar_t(0));
  
  checkNLocallyUnique(name, gids, size_t(nKeysHalf-1));
  }

  {
  // Test 4:
  // Key has four long long entries
  // Each proc has (rank+1)%2 keys; odd-numbered ranks are empty
  // Keys are all identical {0, 1, 2, 3}

  std::string name = "test4: ";

  typedef long long scalar_t;
  typedef int lno_t;
  typedef long long gno_t;

  const size_t nVecs = 4;
  const size_t nKeys = (me+1)%2;

  Tpetra::global_size_t gNEntries = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  typedef Tpetra::Map<lno_t, gno_t> map_t;
  Teuchos::RCP<const map_t> map = rcp(new map_t(gNEntries, nKeys, 0, comm));

  Tpetra::MultiVector<scalar_t, lno_t, gno_t> keys(map, nVecs);
  Tpetra::Vector<scalar_t, lno_t, gno_t> gids(map);

  for (size_t i = 0; i < nKeys; i++) {
    keys.replaceLocalValue(i, 0, scalar_t(0));
    keys.replaceLocalValue(i, 1, scalar_t(1));
    keys.replaceLocalValue(i, 2, scalar_t(2));
    keys.replaceLocalValue(i, 3, scalar_t(3));
  }

  size_t nUniqueGids = findUniqueGids<scalar_t, lno_t, gno_t>(keys, gids);

  // Test for correctness
  checkNUnique(name, nUniqueGids, size_t(1));

  checkMaxGid(name, gids, scalar_t(0));

  checkMinGid(name, gids, scalar_t(0));
  
  checkNLocallyUnique(name, gids, (nKeys ? size_t(1): size_t(0)));
  }
  return 0;
}
