
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

#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> NotayAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
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

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Graph");
    Input(currentLevel, "A");
    Input(currentLevel, "DofsPerNode");

    const ParameterList& pL = GetParameterList();

    if (pL.get<bool>("aggregation: compute aggregate qualities")) {
        Input(currentLevel, "AggregateQualities");
    }

    
  }

  
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
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
    Build_InitialAggregation(pL,A,aggregates,aggStat,numNonAggregatedNodes);
    

      
    
    // DO stuff
    
    aggregates->AggregatesCrossProcessors(false);
    Set(currentLevel, "Aggregates", aggregates);
    GetOStream(Statistics0) << aggregates->description() << std::endl;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void NotayAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build_InitialAggregation(const Teuchos::ParameterList& params,
											    const Matrix& A,
											    Aggregates& aggregates,
											    std::vector<unsigned>& aggStat,
											    LO& numNonAggregatedNodes) const {
    FactoryMonitor m(*this, "Build_InitialAggregation", currentLevel);
    using STS = Teuchos::ScalarTraits<Scalar>;
    using MT = typename STS::magnitudeType;
    MT MT_ZERO = Teuchos::ScalarTraits<MT>::zero;
    SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero;    
    MT MT_ONE  = Teuchos::ScalarTraits<MT>::one;
    MT MT_TWO  = MT_ONE + MT_ONE;
    LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid;
    
    const MT kappa = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
    const MT kappa_init  = kappa / (kappa - MT_TWO);
    
    LO numRows = aggStat.size();    
    numNonAggregatedNodes = numRows;
    
    
    // FIXME: Assumes 1 dof per node


    // Extract diagonal
    RCP<Vector> ghostedDiag = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
    RCP<Vector> ghostedRowSum = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDeletedRowsum(*A);
    RCP<Vector> ghostedAbsRowSum = MueLu::Utilities<MT,LO,GO,NO>::GetMatrixOverlappedAbsDeletedRowsum(*A);
    const ArrayRCP<const SC> D     = ghostedDiag->getData(0);
    const ArrayRCP<const SC> S     = ghostedRowSum->getData(0);
    const ArrayRCP<const SC> AbsRs = ghostedAbsRowSum->getData(0);
           
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
    LO current_idx;
    for (current_idx = 0; current_idx < numRows; current_idx++)
      if(aggStat[current_idx] == READY)
	break;
    
    
    // 2 : Iteration
    while (numNonAggregatedNodes > 0) {
      MT best_mu = MT_ZERO;
      bool have_buddy = false;
      LO best_idx = LO_INVALID;
      
      size_t nnz = A->getNumEntriesInLocalRow(current_idx);
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A->getLocalRowView(row, indices, vals);

      MT aii = STS::real(D[current_idx]);
      MT si  = S[current_idx];
      for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
	MT aij = STS::real(vals[colID]);
	MT ajj = STS::real(D[colID]);
	MT sj  = S[colID];
	// FIXME: Add isaggregated check
	if(aii - si + ajj - sj >= MT_ZERO) {
	  MT mu_top    = MT_TWO / ( MT_ONE / aii + MT_ONE / ajj);
	  MT mu_bottom =  - aij + MT_ONE / ( MT_ONE / (aii - sj) + MT_ONE / (ajj - sj) );
	  MT mu = mu_top / mu_bottom;
	  if (best_idx == LO_INVALID ||  mymu < best_mu) {
	    best_mu  = my;
	    best_idx = row;
	  }
	}
	
      }
    }

    // next
  }




} //namespace MueLu

#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DEF_HPP_ */
