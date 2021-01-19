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
#ifndef _ZOLTAN2_ALGQUOTIENT_HPP_
#define _ZOLTAN2_ALGQUOTIENT_HPP_

#include <Zoltan2_CommGraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>

/////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgQuotient.hpp
//  
// This algorithm partitions the CommGraphModel (communication graph model)
//   and transfers the solution (which is only stored on active ranks of 
//   CommGraphModel) to all MPI ranks. Please see model/Zoltan2_CommGraphModel.hpp
//   for more details.
/////////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter>
class AlgQuotient : public Algorithm<Adapter>
{
public:

  typedef CommGraphModel<typename Adapter::base_adapter_t> graphModel_t;
  typedef typename Adapter::part_t part_t;

  /*! AlgQuotient constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param model a graph
   *  \param algName the (inner) algorithm to partition the quotient graph
   */
  AlgQuotient(const RCP<const Environment> &env__,
              const RCP<const Comm<int> > &problemComm__,
              const RCP<graphModel_t> &model__,
	      const std::string algName__) :
    env(env__), problemComm(problemComm__), 
    model(model__)
  {
    
    if(algName__ == "parmetis")
      this->innerAlgorithm = rcp(new AlgParMETIS<Adapter, graphModel_t>(env,
									problemComm,
									model));
    else
      throw std::logic_error(algName__ + " is not implemented in AlgQuotient yet\n");
    
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  void migrateBack(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<graphModel_t > model;

  RCP<Algorithm<Adapter>> innerAlgorithm;            // algorithm to partition the quotient graph
  RCP<PartitioningSolution<Adapter>> quotientSolution;     // the solution stored on active ranks
};


/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgQuotient<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  if(model->getMigrated()){

    // Create a different object for pre-migration solution 
    PartitioningSolution<Adapter> *soln = NULL;
    try{
      soln = new PartitioningSolution<Adapter>(env, problemComm, 1);
    }
    Z2_FORWARD_EXCEPTIONS;
    quotientSolution = rcp(soln);

    try {
      this->innerAlgorithm->partition(quotientSolution);
    }
    Z2_FORWARD_EXCEPTIONS;

    // Migrate the solution 
    migrateBack(solution); 
  }
  else{

    try {
      this->innerAlgorithm->partition(solution);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

}

// Pre-condition: 
//   The partitioning solution in quotientSolution is only stored on active MPI ranks,
//   which is a single rank if commGraphModel has less than 1024 vertices. 
// Post-condition:
//   The partitioning solution in output parameter solution is distributed to all MPI ranks, 
//   so that each rank has a single entry of the solution.
template <typename Adapter>
void AlgQuotient<Adapter>::migrateBack(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
    int me = problemComm->getRank();  
    int nActiveRanks = model->getNumActiveRanks();
    int dRank = model->getDestinationRank();

    Teuchos::ArrayRCP<part_t> parts(1); // Each rank has a single vertex for now.
    RCP<CommRequest<int>> *requests = new RCP<CommRequest<int>>[1];

    
    requests[0] = Teuchos::ireceive<int, part_t>(*problemComm, parts, dRank);

    if(me < nActiveRanks){

      const part_t *qtntSlnView = quotientSolution->getPartListView();

      int sRank = model->getStartRank();
      int eRank = model->getEndRank();

      ArrayView<size_t> vtxdist; 
      model->getVertexDist(vtxdist);
      for(int i = sRank; i < eRank; i++)
	Teuchos::send<int, part_t>(*problemComm, 1, &qtntSlnView[i-sRank], i);

    }
    
    Teuchos::waitAll<int>(*problemComm, Teuchos::arrayView(requests, 1)); 
    
    solution->setParts(parts);
  
}

} // namespace Zoltan2

#endif
