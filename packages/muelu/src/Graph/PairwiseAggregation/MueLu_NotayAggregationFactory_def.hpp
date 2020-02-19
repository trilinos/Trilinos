
// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_NotayAggregationFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Types.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: pairwise: size");
    SET_VALID_ENTRY("aggregation: compute aggregate qualities");
    SET_VALID_ENTRY("aggregation: Dirichlet threshold");
#undef SET_VALID_ENTRY


   
    /*
    validParamList->getEntry("aggregation: ordering").setValidator(
	  rcp(new validatorType(Teuchos::tuple<std::string>("natural", "graph", "random", "cuthill-mckee"), "aggregation: ordering")));
    */    

    // general variables needed in AggregationFactory
    validParamList->set< RCP<const FactoryBase> >("A",           null, "Generating factory of the matrix");
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
    validParamList->set< RCP<const FactoryBase> >("AggregateQualities", null, "Generating factory for variable \'AggregateQualities\'");

    
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
    Input(currentLevel, "A");
    Input(currentLevel, "DofsPerNode");

    const ParameterList& pL = GetParameterList();

    if (pL.get<bool>("aggregation: compute aggregate qualities")) {
        Input(currentLevel, "AggregateQualities");
    }

    
  }

  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    const ParameterList& pL = GetParameterList();

    // Parameters
    int pairsize=2;
    if (pL.isParameter("aggregation: pairwise: size"))
      pairsize = pL.get<int>("aggregation: pairwise: size");
    TEUCHOS_TEST_FOR_EXCEPTION(pairsize != 2 && pairsize != 4 && pairsize != 8,
			       Exceptions::RuntimeError,
			       "NotayAggregationFactory::Build(): pairsize needs to be a power of two");


    RCP<const GraphBase> graph = Get< RCP<GraphBase> >(currentLevel, "Graph");
    RCP<const Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    // Setup aggregates & aggStat objects
    RCP<Aggregates> aggregates = rcp(new Aggregates(*graph));
    aggregates->setObjectLabel("UC");

    const LO numRows = graph->GetNodeNumVertices();

    // construct aggStat information
    std::vector<unsigned> aggStat(numRows, READY);


    const int DofsPerNode =  Get<int>(currentLevel,"DofsPerNode");
    TEUCHOS_TEST_FOR_EXCEPTION(DofsPerNode != 1, Exceptions::RuntimeError,
			       "Pairwise only supports one dof per node");
    
    // This follows the paper:
    // Notay, "Aggregation-based algebraic multigrid for convection-diffusion equations", SISC 34(3), pp. A2288-2316.


    // FIXME: Do the ordering
    

    // Get the party stated
    LO numNonAggregatedNodes;
    Build_InitialAggregation(pL,A,*aggregates,aggStat,numNonAggregatedNodes);
    

      
    
    // DO stuff
    
    aggregates->AggregatesCrossProcessors(false);
    Set(currentLevel, "Aggregates", aggregates);
    GetOStream(Statistics0) << aggregates->description() << std::endl;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build_InitialAggregation(const Teuchos::ParameterList& params,
                                                                                                    const RCP<const Matrix>& A,
                                                                                                    Aggregates& aggregates,
                                                                                                    std::vector<unsigned>& aggStat,
                                                                                                    LO& numNonAggregatedNodes) const {
    Monitor m(*this, "Build_InitialAggregation");
    using STS = Teuchos::ScalarTraits<Scalar>;
    using MT = typename STS::magnitudeType;
    using RealValuedVector = Xpetra::Vector<MT,LocalOrdinal,GlobalOrdinal,Node>;
      
    MT MT_ZERO = Teuchos::ScalarTraits<MT>::zero();
    SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();    
    MT MT_ONE  = Teuchos::ScalarTraits<MT>::one();
    MT MT_TWO  = MT_ONE + MT_ONE;
    LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    LO LO_ZERO = Teuchos::OrdinalTraits<LO>::zero();

    const MT kappa = STS::magnitude(as<SC>(params.get<double>("aggregation: Dirichlet threshold")));
    const MT kappa_init  = kappa / (kappa - MT_TWO);
    
    LO numRows = aggStat.size();    
    numNonAggregatedNodes = numRows;
    const int myRank  = A->getMap()->getComm()->getRank();
    
    // FIXME: Assumes 1 dof per node


    // Extract diagonal, rowsums, etc
    RCP<Vector> ghostedDiag = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
    RCP<Vector> ghostedRowSum = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDeletedRowsum(*A);
    RCP<RealValuedVector> ghostedAbsRowSum = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedAbsDeletedRowsum(*A);
    const ArrayRCP<const SC> D     = ghostedDiag->getData(0);
    const ArrayRCP<const SC> S     = ghostedRowSum->getData(0);
    const ArrayRCP<const MT> AbsRs = ghostedAbsRowSum->getData(0);
      
    // Aggregates stuff
    ArrayRCP<LO> vertex2AggId_rcp = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner_rcp   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    ArrayView<LO> vertex2AggId    = vertex2AggId_rcp();
    ArrayView<LO> procWinner      = procWinner_rcp();
     
    // Algorithm 4.2

    // 0,1 : Initialize: Flag boundary conditions
    // Modification: We assume symmetry here aij = aji
    
    for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getNodeNumElements()); ++row) {
      MT aii = STS::magnitude(D[row]);
      MT rowsum = AbsRs[row];
      
      if(aii >= kappa_init * rowsum) {
	aggStat[row] = IGNORED;
	numNonAggregatedNodes--;
      }
    }

    
    // FIXME: Add ordering here
    
    // 2 : Iteration
    LO current_idx = 0;
    LO aggIndex = LO_ZERO;
    while (numNonAggregatedNodes > 0 && current_idx < numRows) {
      // If we're aggregated already, skip this guy
      if(aggStat[current_idx] != READY) 
        current_idx++;

      MT best_mu = MT_ZERO;
      LO best_idx = LO_INVALID;
      
      size_t nnz = A->getNumEntriesInLocalRow(current_idx);
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A->getLocalRowView(current_idx, indices, vals);

      MT aii = STS::real(D[current_idx]);
      MT si  = STS::real(S[current_idx]);
      for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
        // Skip aggregated neighbors
        if(aggStat[indices[colID]] != READY) 
          continue;

        // Skips off-rank neighbors
        if(indices[colID] > numRows)
          continue;

	MT aij = STS::real(vals[colID]);
	MT ajj = STS::real(D[colID]);
	MT sj  = STS::real(S[colID]);
	// FIXME: Add isaggregated check
	if(aii - si + ajj - sj >= MT_ZERO) {
	  MT mu_top    = MT_TWO / ( MT_ONE / aii + MT_ONE / ajj);
	  MT mu_bottom =  - aij + MT_ONE / ( MT_ONE / (aii - sj) + MT_ONE / (ajj - sj) );
	  MT mu = mu_top / mu_bottom;
	  if (best_idx == LO_INVALID ||  mu < best_mu) {
	    best_mu  = mu;
	    best_idx = indices[colID];
	  }
	}	
      }// end column for loop
      
      if(best_idx == LO_INVALID) {
        // We found no potential node-buddy, so let's just make this a singleton        
        // NOTE: The behavior of what to do if you have no un-aggregated neighbors is not specified in
        // the paper        
        aggStat[current_idx] = ONEPT;
        vertex2AggId[current_idx] = aggIndex;
        procWinner[current_idx]   = myRank;
	numNonAggregatedNodes--;
        aggIndex++;
      }
      else {
        // We have a buddy, so aggregate, either as a singleton or as a pair, depending on mu
        if(best_mu <= kappa) { 
          LO mybuddy_idx = indices[best_idx];
          aggStat[current_idx] = AGGREGATED;
          aggStat[mybuddy_idx] = AGGREGATED;
          vertex2AggId[current_idx] = aggIndex;
          vertex2AggId[mybuddy_idx] = aggIndex;
          procWinner[current_idx]   = myRank;
          procWinner[mybuddy_idx]   = myRank;
          numNonAggregatedNodes-=2;
          aggIndex++;
        }
        else {
          aggStat[current_idx] = ONEPT;
          vertex2AggId[current_idx] = aggIndex;
          procWinner[current_idx]   = myRank;
          numNonAggregatedNodes--;
          aggIndex++;
        }
      }
    }// end Algorithm 4.2

    // update aggregate object
    aggregates.SetNumAggregates(aggIndex);
    
  }




} //namespace MueLu

#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_ */
