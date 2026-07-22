// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

// Test for issue #2010:  no IDs provided to partitioner

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int fail=0, gfail=0;

  int rank = comm->getRank();

  zlno_t numMyIdentifiers = 0;      // no IDs provided

  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> mydata_t;
  typedef Zoltan2::BasicIdentifierAdapter<mydata_t> adapter_t;

  adapter_t *adapter = new adapter_t(numMyIdentifiers, NULL);

  Teuchos::ParameterList params("test parameters");
  params.set("num_global_parts", 4);
  params.set("algorithm", "block");
  params.set("partitioning_approach", "partition");
  
  Zoltan2::PartitioningProblem<adapter_t> problem(adapter, &params);

  problem.solve();

  Zoltan2::PartitioningSolution<adapter_t> solution = problem.getSolution();

  if (solution.getActualGlobalNumberOfParts() != 0) 
    fail = true;

  gfail = globalFail(*comm, fail);

  if (gfail){
    if (rank==0)
      std::cout << "FAIL GlobalNumberOfParts = " 
                << solution.getActualGlobalNumberOfParts() << std:: endl;
    return 1;
  }

  if (rank==0)
    std::cout << "PASS" << std:: endl;

  return 0;
}
