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
