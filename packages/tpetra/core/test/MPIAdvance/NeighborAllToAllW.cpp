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
// ************************************************************************
// @HEADER
*/

#include <iostream>
#include <random>
#include <vector>

#include <mpi.h>

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_Core.hpp"

#include <mpi_advance.h>

#include "print_vec.hpp"
#include "fake_alltoallw.hpp"

/*! \brief all ranks send and receive nothing

    Passes if it does not crash or error out
*/
void test_nothing(MPI_Comm comm, bool nullBufs, bool sameBufs,
                  bool emptyTypes, Teuchos::FancyOStream &out, bool &success) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  std::vector<int> sendcounts(size, 0);
  std::vector<int> recvcounts(size, 0);

  std::vector<int> senddispls(size, 0);
  std::vector<int> recvdispls(size, 0);
  
  // neighbor call uses different type
  std::vector<MPI_Aint> nbrsenddispls(size, 0);
  std::vector<MPI_Aint> nbrrecvdispls(size, 0);
  
  MPI_Datatype type;
  if (emptyTypes) {
    MPI_Type_contiguous(0, MPI_BYTE, &type);
    MPI_Type_commit(&type);
  } else {
    type = MPI_BYTE;
  }
  std::vector<MPI_Datatype> sendtypes(size, type);
  std::vector<MPI_Datatype> recvtypes(size, type);

  char stackByte0, stackByte1;
  void *sbuf, *rbuf;
  if (nullBufs) {
    sbuf = nullptr;
    rbuf = nullptr;
  } else {
    sbuf = &stackByte0;
    rbuf = &stackByte1;
  }

  if (sameBufs) {
    sbuf = rbuf;
  }

  // create MPIX communicator
  MPIX_Comm *mpixComm = nullptr;
  MPIX_Dist_graph_create_adjacent(
      comm, 0, /*indegree*/
      nullptr, /*sources*/
      nullptr, /*sourceweights*/
      0,       /*outdegree*/
      nullptr /*destinations*/, nullptr /*destweights*/, MPI_INFO_NULL /*info*/,
      0 /*reorder*/, &mpixComm);

  // reference implementation should be okay
  Fake_Alltoallw(sbuf, sendcounts.data(), senddispls.data(), sendtypes.data(), rbuf,
                   recvcounts.data(), recvdispls.data(), recvtypes.data(), comm);

  // MPI advance implementation
  MPIX_Neighbor_alltoallw(sbuf, sendcounts.data(), nbrsenddispls.data(), sendtypes.data(),
                          rbuf, recvcounts.data(), nbrrecvdispls.data(), recvtypes.data(),
                          mpixComm);

  MPIX_Comm_free(mpixComm);

  // we just require that we got this far
  success = true;

  if (emptyTypes) {
    MPI_Type_free(&type);
  }
}

/*! \brief all ranks send and receive some
 */
void test_random(MPI_Comm comm, int seed, Teuchos::FancyOStream &out,
                 bool &success) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // seeded rng so all ranks generate the same communication plan
  // plan [i * size + j] is count i sends to j
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(1, 100);
  std::vector<int> plan;
  std::vector<int> typeindplan;

  // types to choose from
  constexpr int numTypeOps = 4;
  MPI_Datatype typeoptions[numTypeOps] = {MPI_BYTE, MPI_SHORT, MPI_INT, MPI_LONG};
  int typesizes[numTypeOps];
  MPI_Type_size(MPI_BYTE, &(typesizes[0]));
  MPI_Type_size(MPI_SHORT, &(typesizes[1]));
  MPI_Type_size(MPI_INT, &(typesizes[2]));
  MPI_Type_size(MPI_LONG, &(typesizes[3]));
  std::uniform_int_distribution<int> typedist(0, numTypeOps - 1);
  
  for (int i = 0; i < size * size; ++i) {
    plan.push_back(dist(rng));
    typeindplan.push_back(typedist(rng));
  }

  // read my part of the plan
  std::vector<int> sendcounts, recvcounts, senddispls, recvdispls; // alltoallw
  std::vector<MPI_Datatype> sendtypes, recvtypes;
  std::vector<int> nbrsendcounts, nbrrecvcounts; // neighbor alltoallw
  std::vector<MPI_Aint> nbrsenddispls, nbrrecvdispls; // different types in neighbor call
  std::vector<MPI_Datatype> nbrsendtypes, nbrrecvtypes;
  std::vector<int> sources, sourceweights, destinations, destweights; // communicator

  // displs are in bytes
  std::uniform_int_distribution<int> soffsetdist(0, 150);
  rng.seed(seed + rank); // different seed -> different displs per rank
  int initsdispl = soffsetdist(rng);
  int sdispl = initsdispl;
  int nbrsdispl = initsdispl;
  for (int dest = 0; dest < size; ++dest) {
    // pick type
    int ind = typeindplan[rank * size + dest];
    MPI_Datatype mytype = typeoptions[ind];
    int mytypesize = typesizes[ind];
    
    senddispls.push_back(sdispl);
    int count = plan[rank * size + dest];
    sendcounts.push_back(count);
    sendtypes.push_back(mytype);
    
    // pick random displacement
    int offset = soffsetdist(rng);
    if (offset < mytypesize * count) { // prevent overwrites
      offset += mytypesize * count;
    }
    sdispl += offset;

    if (count > 0) {
      destinations.push_back(dest);
      destweights.push_back(count);
      nbrsendcounts.push_back(count);
      nbrsenddispls.push_back((MPI_Aint)nbrsdispl);
      nbrsendtypes.push_back(mytype);
      
      nbrsdispl += offset;
    }
  }

  // displs are in bytes
  std::uniform_int_distribution<int> roffsetdist(0, 150);
  rng.seed(seed + size + rank); // different seed -> different displs (also per rank)
  int initrdispl = roffsetdist(rng);
  int rdispl = initrdispl;
  int nbrrdispl = initrdispl;
  for (int source = 0; source < size; ++source) {
    // pick type
    int ind = typeindplan[source * size + rank];
    MPI_Datatype mytype = typeoptions[ind];
    int mytypesize = typesizes[ind];
    
    recvdispls.push_back(rdispl);
    int count = plan[source * size + rank];
    recvcounts.push_back(count);
    recvtypes.push_back(mytype);
    
    // pick random displacement
    int offset = roffsetdist(rng);
    if (offset < mytypesize * count) { // prevent overwrites
      offset += mytypesize * count;
    }
    rdispl += offset;

    if (count > 0) {
      sources.push_back(source);
      sourceweights.push_back(count);
      nbrrecvcounts.push_back(count);
      nbrrecvdispls.push_back((MPI_Aint)nbrrdispl);
      nbrrecvtypes.push_back(mytype);
      
      nbrrdispl += offset;
    }
  }


  print_vec(rank, "sources", sources);
  print_vec(rank, "sourceweights", sourceweights);
  print_vec(rank, "destinations", destinations);
  print_vec(rank, "destweights", destweights);

  // create MPIX communicator
  MPIX_Comm *mpixComm = nullptr;
  MPIX_Dist_graph_create_adjacent(
      comm, sources.size(), /*indegree*/
      sources.data(),       /*sources*/
      sourceweights.data(), /*sourceweights*/
      destinations.size(),  /*outdegree*/
      destinations.data() /*destinations*/, destweights.data() /*destweights*/,
      MPI_INFO_NULL /*info*/, 0 /*reorder*/, &mpixComm);

  // allocate send/recv bufs
  // displs are in bytes, so the displs are correct no matter types below.
  // (guarranteed by loops above about types)
  std::vector<char> sbuf(sdispl), exp(rdispl), act(rdispl);

  // fill send buf
  std::iota(sbuf.begin(), sbuf.end(), 0); // 0, 1, 2, ...

  // Use reference and MPI_Advance implementation to fill buffers
  Fake_Alltoallw(sbuf.data(), sendcounts.data(), senddispls.data(), sendtypes.data(),
                   exp.data(), recvcounts.data(), recvdispls.data(), recvtypes.data(),
                   comm);

  print_vec(rank, "nbrsendcounts", nbrsendcounts);
  print_vec(rank, "nbrsenddispls", nbrsenddispls);
  // print_vec(rank, "nbrsendtypes", nbrsendtypes);
  print_vec(rank, "nbrrecvcounts", nbrrecvcounts);
  print_vec(rank, "nbrrecvdispls", nbrrecvdispls);
  // print_vec(rank, "nbrrecvtypes", nbrrecvtypes);

  MPIX_Neighbor_alltoallw(sbuf.data(), nbrsendcounts.data(), nbrsenddispls.data(), nbrsendtypes.data(),
                          act.data(), nbrrecvcounts.data(), nbrrecvdispls.data(), nbrrecvtypes.data(),
                          mpixComm);

  MPIX_Comm_free(mpixComm);

  // two recv buffers should be the same
  for (int i = 0; i < rdispl; ++i) {
    TEST_ASSERT(exp[i] == act[i]);
  }
}

static MPI_Comm tpetra_default_comm_as_mpi_comm() {
  Teuchos::RCP<const Teuchos::Comm<int>> teuchosComm = Tpetra::getDefaultComm();

  const Teuchos::MpiComm<int> *teuchosMpiComm =
      dynamic_cast<const Teuchos::MpiComm<int> *>(teuchosComm.get());

  if (teuchosMpiComm) {
    return *(teuchosMpiComm->getRawMpiComm());
  } else {
    throw std::runtime_error(
        "couldn't convert tpetra::getDefaultComm() into MPI_Comm");
  }
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_emptytypes) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_null) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_nullemptytypes) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_same) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_sameemptytypes) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_nullsame) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_nothing_nullsameemptytypes) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborAllToAllW_random) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_random(comm, 42, out, success);
}
