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

/*! \file Zoltan2_findUniqueGids.hpp
    \brief Convert keys stored in std::vector to unique Gids stored in 
           std::vector.  
*/

#ifndef _ZOLTAN2_FINDUNIQUEGIDS_HPP_
#define _ZOLTAN2_FINDUNIQUEGIDS_HPP_

#include <Zoltan2_Standards.hpp>
#include <vector>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

#include <Zoltan2_TPLTraits.hpp>

#include <zoltan_dd.h>
#include <zoltan_dd_const.h>

namespace Zoltan2
{

template <typename gno_t>
size_t findUniqueGidsCommon(
  size_t num_keys,
  int num_gid,
  ZOLTAN_ID_PTR ddkeys,
  char *ddnewgids,
  MPI_Comm mpicomm
)
{
  int num_lid = 0;  // Local IDs not needed
  int debug_level = 0;
  int num_user = sizeof(gno_t);

  Zoltan_DD_Struct *dd = NULL;
  Zoltan_DD_Create(&dd, mpicomm, num_gid, num_lid, num_user, num_keys, 
                   debug_level);

  ZOLTAN_ID_PTR ddnotneeded = NULL;  // Local IDs not needed
  Zoltan_DD_Update(dd, ddkeys, ddnotneeded, ddnewgids, NULL, int(num_keys));

  //////////
  // Insert unique GIDs for DD entries in User data here.

  // Get value of first gid on this rank
  typedef long long mpi_t;
  mpi_t nDDEntries = (mpi_t)(dd->nodecnt);
  mpi_t firstIdx;
  MPI_Scan(&nDDEntries, &firstIdx, 1, MPI_LONG_LONG, MPI_SUM, mpicomm);
  firstIdx -= nDDEntries;  // do not include this rank's entries in prefix sum

  // Loop over all directory entries, updating their userdata with updated gid
  DD_NodeIdx cnt = 0;

  for (DD_NodeIdx i = 0; i < dd->nodelistlen; i++) {
    DD_Node *ptr = &(dd->nodelist[i]);
    if (!(ptr->free)) {
      char *userchar = (char*)(ptr->gid + (dd->gid_length + dd->lid_length));
      gno_t *newgid = (gno_t*) userchar;
      *newgid = gno_t(firstIdx + cnt);
      cnt++;
    }
  }

  ///////////
  // Retrieve the global numbers and put in the result gids vector
  Zoltan_DD_Find(dd, ddkeys, ddnotneeded, ddnewgids, NULL, int(num_keys), NULL);

  Zoltan_DD_Destroy(&dd);

  mpi_t nUnique = 0;
  MPI_Allreduce(&nDDEntries, &nUnique, 1, MPI_LONG_LONG, MPI_SUM, mpicomm);

  return size_t(nUnique);
}

////////////////////////////////////////////////////////////////////////////
template <typename lno_t, typename gno_t>
size_t findUniqueGids(
  Tpetra::MultiVector<gno_t, lno_t, gno_t> &keys,
  Tpetra::Vector<gno_t, lno_t, gno_t> &gids
)
{
  // Input:  Tpetra MultiVector of keys; key length = numVectors()
  //         May contain duplicate keys within a processor.
  //         May contain duplicate keys across processors.
  // Input:  Empty Tpetra Vector with same map for holding the results
  // Output: Filled gids vector, containing unique global numbers for
  //         each unique key.  Global numbers are in range [0,#UniqueKeys).

  size_t num_keys = keys.getLocalLength();
  size_t num_entries = keys.getNumVectors();

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpicomm = Teuchos::getRawMpiComm(*(keys.getMap()->getComm()));
#else
  // Zoltan's siMPI will be used here
  {
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
      int narg = 0;
      char **argv = NULL;
      MPI_Init(&narg, &argv);
    }
  }
  MPI_Comm mpicomm = MPI_COMM_WORLD;  // Will get MPI_COMM_WORLD from siMPI
#endif

  int num_gid = TPL_Traits<ZOLTAN_ID_PTR,gno_t>::NUM_ID * num_entries;
  int num_user = sizeof(gno_t);

  // Buffer the keys for Zoltan_DD
  Teuchos::ArrayRCP<const gno_t> *tmpKeyVecs =
           new Teuchos::ArrayRCP<const gno_t>[num_entries];
  for (size_t v = 0; v < num_entries; v++) tmpKeyVecs[v] = keys.getData(v);

  ZOLTAN_ID_PTR ddkeys = new ZOLTAN_ID_TYPE[num_gid * num_keys];
  size_t idx = 0;
  for (size_t i = 0; i < num_keys; i++) {
    for (size_t v = 0; v < num_entries; v++) {
      ZOLTAN_ID_PTR ddkey = &(ddkeys[idx]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(ddkey, tmpKeyVecs[v][i]);
      idx += TPL_Traits<ZOLTAN_ID_PTR,gno_t>::NUM_ID;
    }
  }
  delete [] tmpKeyVecs;

  // Allocate memory for the result
  char *ddnewgids = new char[num_user * num_keys];
  
  // Compute the new GIDs
  size_t nUnique = findUniqueGidsCommon<gno_t>(num_keys, num_gid,
                                               ddkeys, ddnewgids, mpicomm);

  // Copy the result into the output vector
  gno_t *result = (gno_t *)ddnewgids;
  for (size_t i = 0; i < num_keys; i++)
    gids.replaceLocalValue(i, result[i]);

  // Clean up
  delete [] ddkeys;
  delete [] ddnewgids;

  return nUnique;
}

////////////////////////////////////////////////////////////////////////////
template <typename key_t, typename gno_t>
size_t findUniqueGids(
  std::vector<key_t> &keys,
  std::vector<gno_t> &gids,
  const Teuchos::Comm<int> &comm
)
{
  // Input:  Vector of keys; key length = key_t.size()
  //         Each key must have the same size.  std::array<gno_t, N> is
  //         an example of a good key_t.
  //         May contain duplicate keys within a processor.
  //         May contain duplicate keys across processors.
  // Input:  Empty vector for holding the results
  // Output: Filled gids vector, containing unique global numbers for
  //         each unique key.  Global numbers are in range [0,#UniqueKeys).
  //
  // Note:  This code uses the Zoltan Distributed Directory to assign the
  //        unique global numbers.  Right now, it hacks into the Zoltan_DD
  //        data structures.  If we like this approach, we can add some
  //        elegance to the Zoltan_DD, allowing operations internal to the
  //        directory.

  size_t num_keys = keys.size();
  key_t dummy;
  size_t num_entries = dummy.size();

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpicomm = Teuchos::getRawMpiComm(comm);
#else
  // Zoltan's siMPI will be used here
  {
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
      int narg = 0;
      char **argv = NULL;
      MPI_Init(&narg, &argv);
    }
  }
  MPI_Comm mpicomm = MPI_COMM_WORLD;  // Will get MPI_COMM_WORLD from siMPI
#endif

  int num_gid = TPL_Traits<ZOLTAN_ID_PTR,gno_t>::NUM_ID * num_entries;
  int num_user = sizeof(gno_t);

  // Buffer the keys for Zoltan_DD
  ZOLTAN_ID_PTR ddkeys = new ZOLTAN_ID_TYPE[num_gid * num_keys];
  size_t idx = 0;
  for (size_t i = 0; i < num_keys; i++) {
    for (size_t v = 0; v < num_entries; v++) {
      ZOLTAN_ID_PTR ddkey = &(ddkeys[idx]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(ddkey, keys[i][v]);
      idx += TPL_Traits<ZOLTAN_ID_PTR,gno_t>::NUM_ID;
    }
  }

  // Allocate memory for the result
  char *ddnewgids = new char[num_user * num_keys];

  // Compute the new GIDs
  size_t nUnique = findUniqueGidsCommon<gno_t>(num_keys, num_gid,
                                               ddkeys, ddnewgids, mpicomm);

  // Copy the result into the output vector
  gno_t *result = (gno_t *)ddnewgids;
  for (size_t i = 0; i < num_keys; i++)
    gids[i] = result[i];

  // Clean up
  delete [] ddkeys;
  delete [] ddnewgids;

  return nUnique;
}


}                   // namespace Zoltan2
#endif
