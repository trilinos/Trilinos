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

/*! \file kokkosBlock.cpp
 *  \brief An example of partitioning global ids with Block.
 */

#include <Kokkos_Core.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Zoltan2_BasicKokkosIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

using Teuchos::Comm;
using Teuchos::RCP;

/*! \example kokkosBlock.cpp
 *  An example of the use of the Block algorithm to partition data.
 *  \todo error handling
 *  \todo write some examples that don't use teuchos
 */

int main(int narg, char *arg[]) {
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // For convenience, we'll use the Tpetra defaults for local/global ID types
  // Users can substitute their preferred local/global ID types
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef double scalar_t;

  ///////////////////////////////////////////////////////////////////////
  // Generate some input data.

  int localCount = 40 * (rank + 1);
  int totalCount = 20 * nprocs * (nprocs + 1);
  int targetCount = totalCount / nprocs;
  Kokkos::View<globalId_t*> globalIds("globalIds", localCount);

  if (rank == 0) {
    for (int i = 0, num = 40; i < nprocs ; i++, num += 40) {
      std::cout << "Rank " << i << " generates " << num << " ids." << std::endl;
    }
  }

  globalId_t offset = 0;
  for (int i = 1; i <= rank; i++) {
    offset += 40 * i;
  }

  for (int i = 0; i < localCount; i++) {
    globalIds(i) = offset++;
  }

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter with no weights

  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;
  typedef Zoltan2::BasicKokkosIdentifierAdapter<myTypes> inputAdapter_t;
  typedef typename Zoltan2::BasicKokkosIdentifierAdapter<myTypes>::weight_layout_t Layout;

  const int nWeights = 1;
  Kokkos::View<scalar_t **, Layout> weights("weights", localCount, nWeights);
  for (int index = 0; index < localCount; index++) {
    weights(index, 0) = 1; // Error check relies on uniform weights
  }

  inputAdapter_t ia(globalIds, weights);

  /////////////////////////////////////////////////////////////////////////
  // Create parameters for a Block problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  params.set("algorithm", "block");
  params.set("imbalance_tolerance", 1.1);
  params.set("num_global_parts", nprocs);

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem = 
      new Zoltan2::PartitioningProblem<inputAdapter_t>(&ia, &params);

  ///////////////////////////////////////////////////////////////////////
  // Solve the problem - do the partitioning

  problem->solve();

  ///////////////////////////////////////////////////////////////////////
  // Check and print the solution.
  // Count number of IDs assigned to each part; compare to targetCount

  Kokkos::View<globalId_t *> ids;
  ia.getIDsKokkosView(ids);

  Kokkos::View<int*> partCounts("partCounts", nprocs, 0);

  Kokkos::View<int*> globalPartCounts("globalPartCounts", nprocs);

  for (size_t i = 0; i < ia.getLocalNumIDs(); i++) {
    int pp = problem->getSolution().getPartListView()[i];
    std::cout << rank << " LID " << i << " GID " << ids(i)
              << " PART " << pp << std::endl;
    partCounts(pp)++;
  }

  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, nprocs, 
      &partCounts(0), &globalPartCounts(0));

  if (rank == 0) {
    int ierr = 0;
    for (int i = 0; i < nprocs; i++) {
      if (globalPartCounts(i) != targetCount) {
        std::cout << "FAIL: part " << i << " has " << globalPartCounts(i)
                  << " != " << targetCount << "; " << ++ierr << " errors"
                  << std::endl;
      }
    }
    if (ierr == 0) {
      std::cout << "PASS" << std::endl;
    }
  }

  delete problem;
}

