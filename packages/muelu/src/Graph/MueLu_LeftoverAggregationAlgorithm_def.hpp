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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_LEFTOVERAGGREGATIONALGORITHM_DEF_HPP
#define MUELU_LEFTOVERAGGREGATIONALGORITHM_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_LeftoverAggregationAlgorithm_decl.hpp"

#include "MueLu_Aggregates_decl.hpp" // MUELU_UNASSIGNED macro
#include "MueLu_Utilities_decl.hpp"  // sumAll macro
#include "MueLu_Graph.hpp"
#include "MueLu_UCAggregationCommHelper.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LeftoverAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::LeftoverAggregationAlgorithm():
    phase3AggCreation_(.5),
    minNodesPerAggregate_(1)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LeftoverAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AggregateLeftovers(Graph const &graph, Aggregates &aggregates) const {
    Monitor m(*this, "Leftovers");

    my_size_t nVertices = graph.GetNodeNumVertices();
    int exp_nRows    = aggregates.GetMap()->getNodeNumElements(); // Tentative fix... was previously exp_nRows = nVertices + graph.GetNodeNumGhost();
    int myPid        = graph.GetComm()->getRank();
    my_size_t nAggregates  = aggregates.GetNumAggregates();

    int minNodesPerAggregate = GetMinNodesPerAggregate();

    const RCP<const Map> nonUniqueMap = aggregates.GetMap();
    const RCP<const Map> uniqueMap    = graph.GetDomainMap(); // Q: DomainMap or RowMap??

    MueLu::UCAggregationCommHelper<LO,GO,NO,LMO> myWidget(uniqueMap, nonUniqueMap);

    RCP<Xpetra::Vector<double,LO,GO,NO> > distWeights = Xpetra::VectorFactory<double,LO,GO,NO>::Build(nonUniqueMap);

    // Aggregated vertices not "definitively" assigned to processors are
    // arbitrated by ArbitrateAndCommunicate(). There is some
    // additional logic to prevent losing root nodes in arbitration.
    {
      ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
      ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
      ArrayRCP<double>    weights     = distWeights->getDataNonConst(0);

      for (size_t i=0;i<nonUniqueMap->getNodeNumElements();i++) {
        if (procWinner[i] == MUELU_UNASSIGNED) {
          if (vertex2AggId[i] != MUELU_UNAGGREGATED) {
            weights[i] = 1.;
            if (aggregates.IsRoot(i)) weights[i] = 2.;
          }
        }
      }

      // views on distributed vectors are freed here.
    }

    myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
    // All tentatively assigned vertices are now definitive

    // Tentatively assign any vertex (ghost or local) which neighbors a root
    // to the aggregate associated with the root.
    {
      ArrayRCP<LO>       vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
      ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
      ArrayRCP<double>   weights      = distWeights->getDataNonConst(0);

      for (my_size_t i = 0; i < nVertices; i++) {
        if ( aggregates.IsRoot(i) && (procWinner[i] == myPid) ) {

          // neighOfINode is the neighbor node list of node 'i'.
          ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

          for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            int colj = *it;
            if (vertex2AggId[colj] == MUELU_UNAGGREGATED) {
              weights[colj]= 1.;
              vertex2AggId[colj] = vertex2AggId[i];
            }
          }
        }
      }

      // views on distributed vectors are freed here.
    }

    myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
    // All tentatively assigned vertices are now definitive

    // Record the number of aggregated vertices
    GO total_phase_one_aggregated = 0;
    {
      ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);

      GO phase_one_aggregated = 0;
      for (my_size_t i = 0; i < nVertices; i++) {
        if (vertex2AggId[i] != MUELU_UNAGGREGATED)
          phase_one_aggregated++;
      }

      sumAll(graph.GetComm(), phase_one_aggregated, total_phase_one_aggregated);

      GO local_nVertices = nVertices, total_nVertices = 0;
      sumAll(graph.GetComm(), local_nVertices, total_nVertices);

      /* Among unaggregated points, see if we can make a reasonable size    */
      /* aggregate out of it. We do this by looking at neighbors and seeing */
      /* how many are unaggregated and on my processor. Loosely,            */
      /* base the number of new aggregates created on the percentage of     */
      /* unaggregated nodes.                                                */

      ArrayRCP<double>    weights      = distWeights->getDataNonConst(0);

      double factor = 1.;
      factor = ((double) total_phase_one_aggregated)/((double)(total_nVertices + 1));
      factor = pow(factor, GetPhase3AggCreation());

      for (my_size_t i = 0; i < nVertices; i++) {
        if (vertex2AggId[i] == MUELU_UNAGGREGATED)
          {

            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);
            int rowi_N = neighOfINode.size();

            int nonaggd_neighbors = 0;
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int colj = *it;
              if (vertex2AggId[colj] == MUELU_UNAGGREGATED && colj < nVertices)
                nonaggd_neighbors++;
            }
            if (  (nonaggd_neighbors > minNodesPerAggregate) &&
                  (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
              {
                vertex2AggId[i] = (nAggregates)++;
                for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                  int colj = *it;
                  if (vertex2AggId[colj]==MUELU_UNAGGREGATED) {
                    vertex2AggId[colj] = vertex2AggId[i];
                    if (colj < nVertices) weights[colj] = 2.;
                    else                  weights[colj] = 1.;
                  }
                }
                aggregates.SetIsRoot(i);
                weights[i] = 2.;
              }
          }
      } // for (i = 0; i < nVertices; i++)

      // views on distributed vectors are freed here.
    }

    myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
    //All tentatively assigned vertices are now definitive

    if (IsPrint(Statistics1)) {
      GO Nphase1_agg = nAggregates;
      GO total_aggs;

      sumAll(graph.GetComm(), Nphase1_agg, total_aggs);

      GetOStream(Statistics1, 0) << "Phase 1 - nodes aggregated = " << total_phase_one_aggregated << std::endl;
      GetOStream(Statistics1, 0) << "Phase 1 - total aggregates = " << total_aggs << std::endl;

      GO i = nAggregates - Nphase1_agg;
      { GO ii; sumAll(graph.GetComm(),i,ii); i = ii; }
      GetOStream(Statistics1, 0) << "Phase 3 - additional aggregates = " << i << std::endl;
    }

    // Determine vertices that are not shared by setting Temp to all ones
    // and doing NonUnique2NonUnique(..., ADD). This sums values of all
    // local copies associated with each Gid. Thus, sums > 1 are shared.

    //         std::cout << "exp_nrows=" << exp_nRows << " (nVertices= " << nVertices << ", numGhost=" << graph.GetNodeNumGhost() << ")" << std::endl;
    //         std::cout << "nonUniqueMap=" << nonUniqueMap->getNodeNumElements() << std::endl;

    RCP<Xpetra::Vector<double,LO,GO,NO> > temp_ = Xpetra::VectorFactory<double,LO,GO,NO> ::Build(nonUniqueMap,false); //no need to zero out vector in ctor
    temp_->putScalar(1.);

    RCP<Xpetra::Vector<double,LO,GO,NO> > tempOutput_ = Xpetra::VectorFactory<double,LO,GO,NO> ::Build(nonUniqueMap);

    myWidget.NonUnique2NonUnique(*temp_, *tempOutput_, Xpetra::ADD);

    std::vector<bool> gidNotShared(exp_nRows);
    {
      ArrayRCP<const double> tempOutput = tempOutput_->getData(0);
      for (int i = 0; i < exp_nRows; i++) {
        if (tempOutput[i] > 1.)
          gidNotShared[i] = false;
        else
          gidNotShared[i] = true;
      }
    }

    // Phase 4.
    double nAggregatesTarget;
    nAggregatesTarget = ((double)  uniqueMap->getGlobalNumElements())* (((double) uniqueMap->getGlobalNumElements())/ ((double) graph.GetGlobalNumEdges()));

    GO nAggregatesLocal=nAggregates, nAggregatesGlobal; sumAll(graph.GetComm(), nAggregatesLocal, nAggregatesGlobal);

    LO minNAggs; minAll(graph.GetComm(), nAggregates, minNAggs);
    LO maxNAggs; maxAll(graph.GetComm(), nAggregates, maxNAggs);

    //
    // Only do this phase if things look really bad. THIS
    // CODE IS PRETTY EXPERIMENTAL
    //
#define MUELU_PHASE4BUCKETS 6
    if ((nAggregatesGlobal < graph.GetComm()->getSize()) &&
        (2.5*nAggregatesGlobal < nAggregatesTarget) &&
        (minNAggs ==0) && (maxNAggs <= 1)) {

      // Modify seed of the random algorithm used by temp_->randomize()
      {
        typedef Teuchos::ScalarTraits<double> scalarTrait; // temp_ is of type double.
        scalarTrait::seedrandom(static_cast<unsigned int>(myPid*2 + (int) (11*scalarTrait::random())));
        int k = (int)ceil( (10.*myPid)/graph.GetComm()->getSize());
        for (int i = 0; i < k+7; i++) scalarTrait::random();
        temp_->setSeed(static_cast<unsigned int>(scalarTrait::random()));
      }

      temp_->randomize();

      ArrayRCP<double> temp = temp_->getDataNonConst(0);

      // build a list of candidate root nodes (vertices not adjacent
      // to aggregated vertices)

      my_size_t nCandidates = 0;
      global_size_t nCandidatesGlobal;

      ArrayRCP<LO> candidates = Teuchos::arcp<LO>(nVertices+1);

      double priorThreshold = 0.;
      for (int kkk = 0; kkk < MUELU_PHASE4BUCKETS; kkk++) {

        {
          ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
          ArrayView<const LO> vertex2AggIdView = vertex2AggId();
          RootCandidates(nVertices, vertex2AggIdView, graph, candidates, nCandidates, nCandidatesGlobal);
          // views on distributed vectors are freed here.
        }

        double nTargetNewGuys =  nAggregatesTarget - nAggregatesGlobal;
        double threshold      =  priorThreshold + (1. - priorThreshold)*nTargetNewGuys/(nCandidatesGlobal + .001);

        threshold = (threshold*(kkk+1.))/((double) MUELU_PHASE4BUCKETS);
        priorThreshold = threshold;

        {
          ArrayRCP<LO>     vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
          ArrayRCP<double> weights      = distWeights->getDataNonConst(0);

          for (int k = 0; k < nCandidates; k++ ) {
            int i = candidates[k];
            if ((vertex2AggId[i] == MUELU_UNAGGREGATED) && (fabs(temp[i])  < threshold)) {
              // Note: priorThreshold <= fabs(temp[i]) <= 1

              // neighOfINode is the neighbor node list of node 'iNode'.
              ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

              if (neighOfINode.size() > minNodesPerAggregate) { //TODO: check if this test is exactly was we want to do
                int count = 0;
                for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                  int Adjacent    = *it;
                  // This might not be true if someone close to i
                  // is chosen as a root via fabs(temp[]) < Threshold
                  if (vertex2AggId[Adjacent] == MUELU_UNAGGREGATED){
                    count++;
                    vertex2AggId[Adjacent] = nAggregates;
                    weights[Adjacent] = 1.;
                  }
                }
                if (count >= minNodesPerAggregate) {
                  vertex2AggId[i] = nAggregates++;
                  weights[i] = 2.;
                  aggregates.SetIsRoot(i);
                }
                else { // undo things
                  for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                    int Adjacent    = *it;
                    if (vertex2AggId[Adjacent] == nAggregates){
                      vertex2AggId[Adjacent] = MUELU_UNAGGREGATED;
                      weights[Adjacent] = 0.;
                    }
                  }
                }
              }
            }
          }
          // views on distributed vectors are freed here.
        }
        myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
        // All tentatively assigned vertices are now definitive
        nAggregatesLocal=nAggregates;
        sumAll(graph.GetComm(), nAggregatesLocal, nAggregatesGlobal);

        // check that there are no aggregates sizes below minNodesPerAggregate

        aggregates.SetNumAggregates(nAggregates);

        RemoveSmallAggs(aggregates, minNodesPerAggregate, distWeights, myWidget);

        nAggregates = aggregates.GetNumAggregates();
      }   // one possibility
    }

    // Initialize things for Phase 5. This includes building the transpose
    // of the matrix ONLY for transposed rows that correspond to unaggregted
    // ghost vertices. Further, the transpose is only a local transpose.
    // Nonzero edges which exist on other processors are not represented.


    int observedNAgg=-1; //number of aggregates that contain vertices on this process

    {
      ArrayRCP<LO>       vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
      ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
      for(LO k = 0; k < vertex2AggId.size(); ++k )
        if(vertex2AggId[k]>observedNAgg)
          observedNAgg=vertex2AggId[k];
      observedNAgg++;
    }

    ArrayRCP<int> Mark = Teuchos::arcp<int>(exp_nRows+1);
    ArrayRCP<int> agg_incremented = Teuchos::arcp<int>(observedNAgg);
    ArrayRCP<int> SumOfMarks = Teuchos::arcp<int>(observedNAgg);

    for (int i = 0; i < exp_nRows; i++)   Mark[i] = MUELU_DISTONE_VERTEX_WEIGHT;
    for (int i = 0; i < agg_incremented.size(); i++) agg_incremented[i] = 0;
    for (int i = 0; i < SumOfMarks.size(); i++) SumOfMarks[i] = 0;

    // Grab the transpose matrix graph for unaggregated ghost vertices.
    //     a) count the number of nonzeros per row in the transpose
    std::vector<int> RowPtr(exp_nRows+1-nVertices);
    //{
    ArrayRCP<const LO> vertex2AggIdCst = aggregates.GetVertex2AggId()->getData(0);

    for (int i = nVertices; i < exp_nRows;  i++) RowPtr[i-nVertices] = 0;
    for (int i = 0; i < nVertices;  i++) {

      // neighOfINode is the neighbor node list of node 'iNode'.
      ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

      for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int j = *it;
        if ( (j >= nVertices) && (vertex2AggIdCst[j] == MUELU_UNAGGREGATED)){
          RowPtr[j-nVertices]++;
        }
      }
    }

    //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

    int iSum = 0, iTemp;
    for (int i = nVertices; i < exp_nRows;  i++) {
      iTemp = RowPtr[i-nVertices];
      RowPtr[i-nVertices] = iSum;
      iSum += iTemp;
    }
    RowPtr[exp_nRows-nVertices] = iSum;
    std::vector<LO> cols(iSum+1);

    //     c) Traverse matrix and insert entries in proper location.
    for (int i = 0; i < nVertices;  i++) {

      // neighOfINode is the neighbor node list of node 'iNode'.
      ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

      for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int j = *it;
        if ( (j >= nVertices) && (vertex2AggIdCst[j] == MUELU_UNAGGREGATED)){
          cols[RowPtr[j-nVertices]++] = i;
        }
      }
    }

    //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.
    for (int i = exp_nRows; i > nVertices;  i--)
      RowPtr[i-nVertices] = RowPtr[i-1-nVertices];
    RowPtr[0] = 0;

    // views on distributed vectors are freed here.
    vertex2AggIdCst = Teuchos::null;
    //}

    int bestScoreCutoff;
    int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

    // Stick unaggregated vertices into existing aggregates as described above.

    {
      int ncalls=0;

      for (int kk = 0; kk < 10; kk += 2) {
        bestScoreCutoff = thresholds[kk];

        ArrayRCP<LO> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
        ArrayRCP<const LO> procWinner = aggregates.GetProcWinner()->getData(0);
        ArrayRCP<double> weights       = distWeights->getDataNonConst(0);

        for (int i = 0; i < exp_nRows; i++) {

          if (vertex2AggId[i] == MUELU_UNAGGREGATED) {

            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const LO> neighOfINode;

            // Grab neighboring vertices which is either in graph for local ids
            // or sits in transposed fragment just constructed above for ghosts.
            if (i < nVertices) {
              neighOfINode = graph.getNeighborVertices(i);
            }
            else {
              LO *rowi_col = NULL, rowi_N;
              rowi_col = &(cols[RowPtr[i-nVertices]]);
              rowi_N   = RowPtr[i+1-nVertices] - RowPtr[i-nVertices];

              neighOfINode = ArrayView<const LO>(rowi_col, rowi_N);
            }
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];

              //Adjacent is aggregated and either I own the aggregate
              // or I could own the aggregate after arbitration.
              if ((AdjacentAgg != MUELU_UNAGGREGATED) &&
                  ((procWinner[Adjacent] == myPid) ||
                   (procWinner[Adjacent] == MUELU_UNASSIGNED))){
                SumOfMarks[AdjacentAgg] += Mark[Adjacent];
              }
            }
            int best_score = MUELU_NOSCORE;
            int best_agg = -1;
            int BestMark = -1;
            bool cannotLoseAllFriends=false; // Used to address possible loss of vertices in arbitration of shared nodes discussed above. (Initialized to false only to avoid a compiler warning).

            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];
              //Adjacent is unaggregated, has some value and no
              //other processor has definitively claimed him
              if ((AdjacentAgg != MUELU_UNAGGREGATED) &&
                  (SumOfMarks[AdjacentAgg] != 0) &&
                  ((procWinner[Adjacent] == myPid) ||
                   (procWinner[Adjacent] == MUELU_UNASSIGNED ))) {

                // first figure out the penalty associated with
                // AdjacentAgg having already been incremented
                // during this phase, then compute score.

                double penalty = (double) (INCR_SCALING*agg_incremented[AdjacentAgg]);
                if (penalty > MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]))
                  penalty = MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]);
                int score = SumOfMarks[AdjacentAgg]- ((int) floor(penalty));

                if (score > best_score) {
                  best_agg             = AdjacentAgg;
                  best_score           = score;
                  BestMark             = Mark[Adjacent];
                  cannotLoseAllFriends = false;

                  // This address issue mentioned above by checking whether
                  // Adjacent could be lost in arbitration. weight==0 means that
                  // Adjacent was not set during this loop of Phase 5 (and so it
                  // has already undergone arbitration). GidNotShared == true
                  // obviously implies that Adjacent cannot be lost to arbitration
                  if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                    cannotLoseAllFriends = true;
                }
                // Another vertex within current best aggregate found.
                // We should have (best_score == score). We need to see
                // if we can improve BestMark and cannotLoseAllFriends.
                else if (best_agg == AdjacentAgg) {
                  if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                    cannotLoseAllFriends = true;
                  if (Mark[Adjacent] > BestMark) BestMark = Mark[Adjacent];
                }
              }
            }
            // Clean up
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];
              if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
            }
            // Tentatively assign vertex to best_agg.
            if ( (best_score >= bestScoreCutoff) && (cannotLoseAllFriends)) {

              TEUCHOS_TEST_FOR_EXCEPTION(best_agg == -1 || BestMark == -1, MueLu::Exceptions::RuntimeError, "MueLu::UCAggregationFactory internal error"); // should never happen

              vertex2AggId[i] = best_agg;
              weights[i] = best_score;
              agg_incremented[best_agg]++;
              Mark[i] = (int) ceil(   ((double) BestMark)/2.);
            }
          }

          // views on distributed vectors are freed here.
        }

        vertex2AggId = Teuchos::null;
        procWinner   = Teuchos::null;
        weights      = Teuchos::null;

        ++ncalls;
        myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
        // All tentatively assigned vertices are now definitive
      }

      //       if (graph.GetComm()->getRank()==0)
      //         std::cout << "#calls to Arb&Comm=" << ncalls << std::endl;
    }

    // Phase 6: Aggregate remain unaggregated vertices and try at all costs
    //          to avoid small aggregates.
    //          One case where we can find ourselves in this situation
    //          is if all vertices vk adjacent to v have already been
    //          put in other processor's aggregates and v does not have
    //          a direct connection to a local vertex in any of these
    //          aggregates.

    int Nleftover = 0, Nsingle = 0;
    {

      ArrayRCP<LO> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
      ArrayRCP<double> weights       = distWeights->getDataNonConst(0);
      ArrayRCP<const LO> procWinner = aggregates.GetProcWinner()->getData(0);

      int count = 0;
      for (my_size_t i = 0; i < nVertices; i++) {
        if (vertex2AggId[i] == MUELU_UNAGGREGATED) {
          Nleftover++;

          // neighOfINode is the neighbor node list of node 'iNode'.
          ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

          // We don't want too small of an aggregate. So lets see if there is an
          // unaggregated neighbor that we can also put with this vertex

          vertex2AggId[i] = nAggregates;
          weights[i] = 1.;
          if (count == 0) aggregates.SetIsRoot(i);
          count++;
          for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            int j = *it;
            if ((j != i)&&(vertex2AggId[j] == MUELU_UNAGGREGATED)&&
                (j < nVertices)) {
              vertex2AggId[j] = nAggregates;
              weights[j] = 1.;
              count++;
            }
          }
          if ( count >= minNodesPerAggregate) {
            nAggregates++;
            count = 0;
          }
        }
      }

      // We have something which is under minNodesPerAggregate when
      if (count != 0) {
#ifdef FIXME
        // Can stick small aggregate with 0th aggregate?
        if (nAggregates > 0) {
          for (my_size_t i = 0; i < nVertices; i++) {
            if ((vertex2AggId[i] == nAggregates) && (procWinner[i] == myPid)) {
              vertex2AggId[i] = 0;
              aggregates.SetIsRoot(i,false);
            }
          }
        }
        else {
          Nsingle++;
          nAggregates++;
        }
#else
        // Can stick small aggregate with 0th aggregate?
        if (nAggregates > 0) {
          for (my_size_t i = 0; i < nVertices; i++) {
            // TW: This is not a real fix. This may produce ugly bad aggregates!
            // I removed the procWinner[i] == myPid check. it makes no sense to me since
            // it leaves vertex2AggId[i] == nAggregates -> crash in ComputeAggregateSizes().
            // Maybe it's better to add the leftovers to the last generated agg on the current proc.
            // The best solution would be to add them to the "next"/nearest aggregate, that may be
            // on an other processor
            if (vertex2AggId[i] == nAggregates) {
              vertex2AggId[i] = nAggregates-1; //0;
              aggregates.SetIsRoot(i,false);
            }
          }
        }
        else {
          Nsingle++;
          nAggregates++;
        }
#endif
      }

      // views on distributed vectors are freed here.
    }

    myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, false);

    if (IsPrint(Statistics1)) {
      GO total_Nsingle=0;   sumAll(graph.GetComm(), (GO)Nsingle,     total_Nsingle);
      GO total_Nleftover=0; sumAll(graph.GetComm(), (GO)Nleftover,   total_Nleftover);
      // GO total_aggs;        sumAll(graph.GetComm(), (GO)nAggregates, total_aggs);
      // GetOStream(Statistics1, 0) << "Phase 6 - total aggregates = " << total_aggs << std::endl;
      GetOStream(Statistics1, 0) << "Phase 6 - leftovers = " << total_Nleftover << " and singletons = " << total_Nsingle << std::endl;
    }

    aggregates.SetNumAggregates(nAggregates);

  } //AggregateLeftovers

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LeftoverAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RootCandidates(my_size_t nVertices, ArrayView<const LO> & vertex2AggId, Graph const &graph,
                      ArrayRCP<LO> &candidates, my_size_t &nCandidates, global_size_t &nCandidatesGlobal) const
  {
    nCandidates = 0;

    for (my_size_t i = 0; i < nVertices; i++ ) {
      if (vertex2AggId[i] == MUELU_UNAGGREGATED) {
        bool noAggdNeighbors = true;

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

        for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int adjacent    = *it;
          if (vertex2AggId[adjacent] != MUELU_UNAGGREGATED)
            noAggdNeighbors = false;
        }
        if (noAggdNeighbors == true) candidates[nCandidates++] = i;
      }
    }

    sumAll(graph.GetComm(), (GO)nCandidates, nCandidatesGlobal);

  } //RootCandidates

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int LeftoverAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RemoveSmallAggs(Aggregates& aggregates, int min_size,
                      RCP<Xpetra::Vector<double,LO,GO,NO> > & distWeights, const MueLu::UCAggregationCommHelper<LO,GO,NO,LMO> & myWidget) const {
    int myPid = aggregates.GetMap()->getComm()->getRank();

    LO nAggregates = aggregates.GetNumAggregates();

    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    //ArrayRCP<int> AggInfo = Teuchos::arcp<int>(nAggregates+1);
    //aggregates.ComputeAggSizes(AggInfo);
    ArrayRCP<LO> AggInfo = aggregates.ComputeAggregateSizes();

    ArrayRCP<double> weights = distWeights->getDataNonConst(0);

    // Make a list of all aggregates indicating New AggId
    // Use AggInfo array for this.

    LO NewNAggs = 0;
    for (LO i = 0; i < nAggregates; i++) {
      if ( AggInfo[i] < min_size) {
        AggInfo[i] =  MUELU_UNAGGREGATED;
      }
      else AggInfo[i] = NewNAggs++;
    }

    for (LO k = 0; k < size; k++ ) {
      if (procWinner[k] == myPid) {
        if (vertex2AggId[k] !=  MUELU_UNAGGREGATED) {
          vertex2AggId[k] = AggInfo[vertex2AggId[k]];
          weights[k] = 1.;
        }
        if (vertex2AggId[k] ==  MUELU_UNAGGREGATED)
          aggregates.SetIsRoot(k,false);
      }
    }
    nAggregates = NewNAggs;

    myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
    // All tentatively assigned vertices are now definitive

    // procWinner is not set correctly for aggregates which have
    // been eliminated
    for (LO i = 0; i < size; i++) {
      if (vertex2AggId[i] == MUELU_UNAGGREGATED)
        procWinner[i] = MUELU_UNASSIGNED;
    }
    aggregates.SetNumAggregates(nAggregates);

    return 0; //TODO
  } //RemoveSmallAggs

} //namespace MueLu

#endif // MUELU_LEFTOVERAGGREGATIONALGORITHM_DEF_HPP
