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
// This algorithm
//   - creates a CommGraphModel (communication graph model)
//   - partitions the CommGraphModel using ParMETIS, and
//   - transfers the solution (which is only stored on active ranks of
//   CommGraphModel) to all MPI ranks.
//
// Please see model/Zoltan2_CommGraphModel.hpp for more details.
/////////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter>
class AlgQuotient : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef CommGraphModel<typename Adapter::base_adapter_t> graphModel_t;


  /*! AlgQuotient constructors
   *  \param env          parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter      the user's input adapter
   *
   *  The algorithm works only on a XpetraCrsGraphAdpater for now.
   */
  AlgQuotient(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const IdentifierAdapter<user_t> > &adapter__,
          const modelFlag_t& graphFlags_) :
    env(env__), problemComm(problemComm__), adapter(adapter__), graphFlags(graphFlags_)
  {
    std::string errStr = "cannot build CommGraphModel from IdentifierAdapter, ";
    errStr            += "AlgQuotient requires Graph Adapter";
    throw std::runtime_error(errStr);
  }

  AlgQuotient(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const VectorAdapter<user_t> > &adapter__,
          const modelFlag_t& graphFlags_) :
    env(env__), problemComm(problemComm__), adapter(adapter__), graphFlags(graphFlags_)
  {
    std::string errStr = "cannot build CommGraphModel from VectorAdapter, ";
    errStr            += "AlgQuotient requires Graph Adapter";
    throw std::runtime_error(errStr);
  }

  AlgQuotient(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const MatrixAdapter<user_t,userCoord_t> > &adapter__,
          const modelFlag_t& graphFlags_) :
    env(env__), problemComm(problemComm__), adapter(adapter__), graphFlags(graphFlags_)
  {
    std::string errStr = "cannot build CommGraphModel from MatrixAdapter, ";
    errStr            += "AlgQuotient has not been implemented for Matrix Adapter yet.";
    throw std::runtime_error(errStr);
  }

  AlgQuotient(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const MeshAdapter<user_t> > &adapter__,
          const modelFlag_t& graphFlags_) :
    env(env__), problemComm(problemComm__), adapter(adapter__), graphFlags(graphFlags_)
  {
    std::string errStr = "cannot build CommGraphModel from MeshAdapter, ";
    errStr            += "AlgQuotient has not been implemented for Mesh Adapter yet.";
    throw std::runtime_error(errStr);
  }

  AlgQuotient(const RCP<const Environment> &env__,
          const RCP<const Comm<int> > &problemComm__,
          const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__,
          const modelFlag_t& graphFlags_) :
    env(env__), problemComm(problemComm__), adapter(adapter__), graphFlags(graphFlags_)
  {
    this->innerAlgorithm =
        rcp(new AlgParMETIS<Adapter, graphModel_t>(env, problemComm, adapter, graphFlags));
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & pl)
  {
    pl.set("quotient_threshold", 1, "threshold for the number of vertices on the active ranks",
	   Environment::getAnyIntValidator());
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  void migrateBack(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const base_adapter_t> adapter;
  modelFlag_t graphFlags;

  RCP<Algorithm<Adapter>> innerAlgorithm;               // algorithm to partition the quotient graph
  RCP<PartitioningSolution<Adapter>> quotientSolution;  // the solution stored on the active ranks

};


/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgQuotient<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  HELLO;

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

  env->memory("Zoltan2-Quotient: After creating solution");

}

// Pre-condition:
//   The partitioning solution in quotientSolution is only stored on active MPI ranks,
//   which is a single rank if commGraphModel has less than 'threshold_' vertices.
// Post-condition:
//   The partitioning solution in output parameter 'solution' is distributed to all MPI ranks,
//   so that each rank gets a partitioning vector with the size of the number of local vertices,
//   and each entry in a goven rank should have the same part value.
template <typename Adapter>
void AlgQuotient<Adapter>::migrateBack(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  const auto model = rcp(new CommGraphModel<base_adapter_t>(
      this->adapter, this->env, this->problemComm));
  int me = problemComm->getRank();
  int nActiveRanks = model->getNumActiveRanks();
  int dRank = model->getDestinationRank();

  // Receive the (single entry) partitioning solution for this MPI rank
  Teuchos::ArrayRCP<part_t> parts(1);
  RCP<CommRequest<int>> *requests = new RCP<CommRequest<int>>[1];
  requests[0] = Teuchos::ireceive<int, part_t>(*problemComm, parts, dRank);
  if (me < nActiveRanks) {

    const part_t *qtntSlnView = quotientSolution->getPartListView();

    int sRank = model->getStartRank();
    int eRank = model->getEndRank();

    ArrayView<size_t> vtxdist;
    model->getVertexDist(vtxdist);
    for (int i = sRank; i < eRank; i++)
      Teuchos::send<int, part_t>(*problemComm, 1, &qtntSlnView[i - sRank], i);
    }
    Teuchos::waitAll<int>(*problemComm, Teuchos::arrayView(requests, 1));

    // Extend the single-entry solution for all local vertices
    size_t numLocalVertices = adapter->getLocalNumIDs();
    Teuchos::ArrayRCP<part_t> extendedParts(numLocalVertices);
    for(size_t i = 0; i < numLocalVertices; i++)
      extendedParts[i] = parts[0];

    // Set the extended partitioning solution
    solution->setParts(extendedParts);

}


} // namespace Zoltan2

////////////////////////////////////////////////////////////////////////


#endif

