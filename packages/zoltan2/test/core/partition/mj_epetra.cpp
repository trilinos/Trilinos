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

/*! \file mj_epetra
    \brief Protects the epetra acess to Zoltan2_XpetraMultiVectorAdapter.hpp.
           Tests via a multijagged partitioning.
*/

#include <vector>
#include <numeric>
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_MultiVector.h>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  typedef int gid_t;

  const int N = 100; // num coords
  const int num_parts = 5;
  const int dim = 3;

  int rank = comm->getRank();

  // params
  Teuchos::ParameterList params("test params");
  params.set("algorithm", "multijagged");
  params.set("num_global_parts", num_parts);

  // create gids
  std::vector<gid_t> global_ids(N);
  std::iota(global_ids.begin(), global_ids.end(), rank * N);

#ifdef HAVE_MPI
  Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm epetra_comm;
#endif

  Epetra_BlockMap map(-1, N, 1, 0, epetra_comm);

  typedef Epetra_MultiVector mv_t;
  Teuchos::RCP<mv_t> mv = Teuchos::rcp(new mv_t(map, dim)); // all 0's

  std::vector<int> stride;
  std::vector<const double *> weights;
  Teuchos::RCP<const mv_t> cmv = Teuchos::rcp_const_cast<const mv_t>(mv);

  typedef Zoltan2::XpetraMultiVectorAdapter<mv_t> inputAdapter_t;

  Teuchos::RCP<inputAdapter_t> ia =
    Teuchos::rcp(new inputAdapter_t(cmv, weights, stride));

  Teuchos::RCP<Zoltan2::PartitioningProblem<inputAdapter_t>> problem =
    Teuchos::rcp(new Zoltan2::PartitioningProblem<inputAdapter_t>(ia.get(), &params));

  problem->solve();

  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
  Teuchos::RCP<quality_t> metricObject = Teuchos::rcp(new quality_t(
    ia.get(), &params, comm, &problem->getSolution()));

  int err = 0;
  if (comm->getRank() == 0) {
    metricObject->printMetrics(std::cout);
    double imb = metricObject->getObjectCountImbalance();
    if (imb <= 1.01)
      std::cout << "balance satisfied " << imb << std::endl;
    else {
      std::cout << "balance failed " << imb << std::endl;
      err++;
    }
  }

  if(rank == 0 && err == 0) {
    std::cout << "PASS" << std::endl;
  }

  return 0;
}

