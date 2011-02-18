#include "Cthulhu_VectorFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AggregationOptions.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AggAlgorithm2Comm.hpp"

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <iostream>


UCAggregationFactory::UCAggregationFactory() {}

virtual UCAggregationFactory::~UCAggregationFactory() {}

Teuchos::RCP<Aggregates> UCAggregationFactory::Build(Graph const &graph, AggregationOptions const &options) const
{
  Teuchos::RCP<Aggregates> aggregates = CoarsenUncoupled(options,graph);
  std::string name = "UC_CleanUp";
  AggregateLeftOvers(options, *aggregates, name, graph);
  return aggregates;
}

RCP<MueLu::Aggregates<LO,GO> >
UCAggregationFactory::CoarsenUncoupled(MueLu::AggregationOptions const & aggOptions, MueLu::Graph<LO,GO> const & graph) const
{
  /* Create Aggregation object */
  const std::string name = "Uncoupled";
  int nAggregates = 0;
  RCP<MueLu::Aggregates<LO,GO> > aggregates = Teuchos::rcp(new MueLu::Aggregates<LO,GO>(graph, name));

  /* ============================================================= */
  /* aggStat indicates whether this node has been aggreated, and   */
  /* vertex2AggId stores the aggregate number where this node has  */
  /* been aggregated into.                                         */
  /* ============================================================= */

  Teuchos::ArrayRCP<NodeState> aggStat;
  const int nRows = graph.GetNodeNumVertices();
  if (nRows > 0) aggStat = Teuchos::arcp<NodeState>(nRows);
  for ( int i = 0; i < nRows; ++i ) aggStat[i] = READY;

  /* unused */
  // Teuchos::ArrayRCP<int> aggCntArray = Teuchos::arcp<int>(nRows+1);
  // for ( int i = 0; i <= nRows; ++i ) aggCntArray[i] = 0;

  /* ============================================================= */
  /* Phase 1  :                                                    */
  /*    for all nodes, form a new aggregate with its neighbors     */
  /*    if the number of its neighbors having been aggregated does */
  /*    not exceed a given threshold                               */
  /*    (aggOptions.GetMaxNeighAlreadySelected() = 0 ===> Vanek's scheme) */
  /* ============================================================= */

  /* some general variable declarations */   
  const int ordering = aggOptions.GetOrdering();
  Teuchos::ArrayRCP<int> randomVector;
  MueLu_Node       *nodeHead=NULL, *nodeTail=NULL, *newNode=NULL;
  MueLu_SuperNode  *aggHead=NULL, *aggCurrent=NULL, *supernode=NULL;
  /**/

  if ( ordering == 1 )       /* random ordering */
  {
    randomVector = Teuchos::arcp<int>(nRows);
    for (int i = 0; i < nRows; ++i) randomVector[i] = i;
      MueLu_RandomReorder(randomVector, *(graph.GetDomainMap()));
  } 
  else if ( ordering == 2 )  /* graph ordering */
  {
    newNode = new MueLu_Node;      
    newNode->nodeId = 0;
    nodeHead = newNode;
    nodeTail = newNode;
    newNode->next = NULL;
  }

  /* main loop */
  {
    int iNode  = 0;
    int iNode2 = 0;
    
    Teuchos::ArrayRCP<int> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0); // output only: contents ignored
    
    while (iNode2 < nRows)
      {
        /*------------------------------------------------------ */
        /* pick the next node to aggregate                       */
        /*------------------------------------------------------ */

        if      ( ordering == 0 ) iNode = iNode2++;
        else if ( ordering == 1 ) iNode = randomVector[iNode2++];
        else if ( ordering == 2 ) 
          {
            if ( nodeHead == NULL ) 
              {
                for ( int jNode = 0; jNode < nRows; ++jNode ) 
                  {
                    if ( aggStat[jNode] == READY )
                      { 
                        newNode = new MueLu_Node;
                        newNode->nodeId = jNode;
                        nodeHead = newNode;
                        nodeTail = newNode;
                        newNode->next = NULL;
                        break;
                      }
                  }
              }
            if ( nodeHead == NULL ) break;
            newNode = nodeHead;
            iNode = newNode->nodeId;
            nodeHead = newNode->next;
            delete newNode;
          }

        /*------------------------------------------------------ */
        /* consider further only if the node is in READY mode    */
        /*------------------------------------------------------ */

        if ( aggStat[iNode] == READY ) 
          {
            // neighOfINode is the neighbor node list of node 'iNode'.
            Teuchos::ArrayView<const int> neighOfINode = graph.getNeighborVertices(iNode);
            int length = neighOfINode.size();
          
            supernode = new MueLu_SuperNode;
            try {
              supernode->list = Teuchos::arcp<int>(length+1);
            } catch (std::bad_alloc&) {
              printf("Error:couldn't allocate memory for supernode! %d\n", length);
              exit(1);
            }

            supernode->maxLength = length;
            supernode->length = 1;
            supernode->list[0] = iNode;
          
            int selectFlag = 1;
            {
              /*--------------------------------------------------- */
              /* count the no. of neighbors having been aggregated  */
              /*--------------------------------------------------- */
            
              int count = 0;
              for (Teuchos::ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                {
                  int index = *it;
                  if ( index < nRows ) 
                    {
                      if ( aggStat[index] == READY || 
                           aggStat[index] == NOTSEL ) 
                        supernode->list[supernode->length++] = index;
                      else count++;
                    
                    }
                }
            
              /*--------------------------------------------------- */
              /* if there are too many neighbors aggregated or the  */
              /* number of nodes in the new aggregate is too few,   */
              /* don't do this one                                  */
              /*--------------------------------------------------- */
            
              if ( count > aggOptions.GetMaxNeighAlreadySelected() ) selectFlag = 0;
            }

            // Note: the supernode length is actually 1 more than the 
            //       number of nodes in the candidate aggregate. The 
            //       root is counted twice. I'm not sure if this is 
            //       a bug or a feature ... so I'll leave it and change
            //       < to <= in the if just below.

            if (selectFlag != 1 || 
                supernode->length <= aggOptions.GetMinNodesPerAggregate()) 
              {
                aggStat[iNode] = NOTSEL;
                delete supernode;
                if ( ordering == 2 ) /* if graph ordering */
                  {
                    for (Teuchos::ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                      {
                        int index = *it;
                        if  ( index < nRows && aggStat[index] == READY )
                          { 
                            newNode = new MueLu_Node;
                            newNode->nodeId = index;
                            newNode->next = NULL;
                            if ( nodeHead == NULL )
                              {
                                nodeHead = newNode;
                                nodeTail = newNode;
                              } else {
                              nodeTail->next = newNode;
                              nodeTail = newNode;
                            }
                          } 
                      } 
                  } 
              } 
            else 
              {
                aggregates->SetIsRoot(iNode);
                for ( int j = 0; j < supernode->length; ++j ) 
                  {
                    int jNode = supernode->list[j];
                    aggStat[jNode] = SELECTED;
                    vertex2AggId[jNode] = nAggregates;
                    if ( ordering == 2 ) /* if graph ordering */
                      {

                        Teuchos::ArrayView<const int> neighOfJNode = graph.getNeighborVertices(jNode);

                        for (Teuchos::ArrayView<const int>::const_iterator it = neighOfJNode.begin(); it != neighOfJNode.end(); ++it)
                          {
                            int index = *it;
                            if ( index < nRows && aggStat[index] == READY )
                              { 
                                newNode = new MueLu_Node;
                                newNode->nodeId = index;
                                newNode->next = NULL;
                                if ( nodeHead == NULL )
                                  {
                                    nodeHead = newNode;
                                    nodeTail = newNode;
                                  } else {
                                  nodeTail->next = newNode;
                                  nodeTail = newNode;
                                }
                              }
                          } 
                      } 
                  }
                supernode->next = NULL;
                supernode->index = nAggregates;
                if ( nAggregates == 0 ) 
                  {
                    aggHead = supernode;
                    aggCurrent = supernode;
                  } 
                else 
                  {
                    aggCurrent->next = supernode;
                    aggCurrent = supernode;
                  } 
                nAggregates++;
                // unused aggCntArray[nAggregates] = supernode->length;
              }
          }
      } // end of 'for'

    // views on distributed vectors are freed here.

  } // end of 'main loop'

  if ( ordering == 2 ) 
    {
      while ( nodeHead != NULL )
        {
          newNode = nodeHead;
          nodeHead = newNode->next;
          delete newNode;
        }
    }

  /* Update aggregate object */  
  aggregates->SetNumAggregates(nAggregates);

  /* Verbose */
  // TODO: replace AllReduce by Reduce to proc 0
  if (aggOptions.GetPrintFlag() < 7) { //FIXME
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();
    int myPid = comm->getRank();
    
    {
      int localReady=0, globalReady;
      
      // Compute 'localReady'
      for ( int i = 0; i < nRows; ++i ) 
        if ( aggStat[i] == READY ) localReady++;
      
      // Compute 'globalReady'
      sumAll(comm, localReady, globalReady);
      
      if (myPid == 0 && globalReady > 0)
        printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",globalReady);
    }
    
    {
      int localSelected=0, globalSelected;
      int globalNRows;
      
      // Compute 'localSelected'
      for ( int i = 0; i < nRows; ++i ) 
          if ( aggStat[i] == SELECTED ) localSelected++;
      
      // Compute 'globalSelected' and 'globalNRows'
      sumAll(comm, localSelected, globalSelected);
      sumAll(comm, nRows, globalNRows);
      
      if (myPid == 0)
        printf("Aggregation(UC) : Phase 1 - nodes aggregated = %d (%d)\n",globalSelected, globalNRows);
    }
    
    {
      int nAggregatesGlobal; 
      sumAll(comm, nAggregates, nAggregatesGlobal);
      if (myPid == 0)
        printf("Aggregation(UC) : Phase 1 - total aggregates = %d \n",nAggregatesGlobal);
    }
    
  } // if myPid == 0 ...
  
  /* ------------------------------------------------------------- */
  /* clean up                                                      */
  /* ------------------------------------------------------------- */

  aggCurrent = aggHead;
  while ( aggCurrent != NULL ) 
    {
      supernode = aggCurrent;
      aggCurrent = aggCurrent->next;
      delete supernode;
    }

  return aggregates;
} //CoarsenUncoupled

int UCAggregationFactory::AggregateLeftOvers(AggregationOptions const &aggOptions, 
          Aggregates &aggregates, std::string const & label, Graph const &graph) const
{
  int nVertices    = graph.GetNodeNumVertices();
  int exp_nRows    = nVertices + graph.GetNodeNumGhost();
  int myPid        = graph.GetComm()->getRank();
  double printFlag = aggOptions.GetPrintFlag();
  int nAggregates  = aggregates.GetNumAggregates();

  int minNodesPerAggregate = aggOptions.GetMinNodesPerAggregate();

  const RCP<const Cthulhu::Map<LO,GO> > nonUniqueMap = aggregates.GetMap();
  const RCP<const Cthulhu::Map<LO,GO> > uniqueMap = graph.GetDomainMap();

  AggAlgorithm2Comm myWidget(uniqueMap, nonUniqueMap);

  RCP<Cthulhu::Vector<double> > distWeights = Cthulhu::VectorFactory<double>::Build(nonUniqueMap);

  // Aggregated vertices not "definitively" assigned to processors are
  // arbitrated by MueLu_ArbitrateAndCommunicate(). There is some
  // additional logic to prevent losing root nodes in arbitration.
  {
    ArrayRCP<const int> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
    ArrayRCP<const int> procWinner   = aggregates.GetProcWinner()->getData(0);
    ArrayRCP<double>    weights      = distWeights->getDataNonConst(0);
    
    distWeights->putScalar(0.);
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

  myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
  distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

  // Tentatively assign any vertex (ghost or local) which neighbors a root
  // to the aggregate associated with the root.
  {
    ArrayRCP<int>       vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<const int> procWinner   = aggregates.GetProcWinner()->getData(0);
    ArrayRCP<double>    weights      = distWeights->getDataNonConst(0);

    for (int i = 0; i < nVertices; i++) { 
      if ( aggregates.IsRoot(i) && (procWinner[i] == myPid) ) {
        
        // neighOfINode is the neighbor node list of node 'i'.
        ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
        
        for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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

  myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
  distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

  // Record the number of aggregated vertices
  int total_phase_one_aggregated = 0;
  {
    ArrayRCP<int> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
      
    int phase_one_aggregated = 0;
    for (int i = 0; i < nVertices; i++) {
      if (vertex2AggId[i] != MUELU_UNAGGREGATED)
        phase_one_aggregated++;
    }
    
    sumAll(graph.GetComm(), phase_one_aggregated, total_phase_one_aggregated);

    int total_nVertices = 0;
    sumAll(graph.GetComm(), nVertices, total_nVertices);
    
    /* Among unaggregated points, see if we can make a reasonable size    */
    /* aggregate out of it. We do this by looking at neighbors and seeing */
    /* how many are unaggregated and on my processor. Loosely,            */
    /* base the number of new aggregates created on the percentage of     */
    /* unaggregated nodes.                                                */

    ArrayRCP<double>    weights      = distWeights->getDataNonConst(0);

    double factor = 1.;
    factor = ((double) total_phase_one_aggregated)/((double)(total_nVertices + 1));
    factor = pow(factor, aggOptions.GetPhase3AggCreation());
    
    for (int i = 0; i < nVertices; i++) {
      if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
        {
          
          // neighOfINode is the neighbor node list of node 'iNode'.
          ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
          int rowi_N = neighOfINode.size();
          
          int nonaggd_neighbors = 0;
          for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            int colj = *it;
            if (vertex2AggId[colj] == MUELU_UNAGGREGATED && colj < nVertices)
              nonaggd_neighbors++;
          }
          if (  (nonaggd_neighbors > minNodesPerAggregate) &&
                (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
            {
              vertex2AggId[i] = (nAggregates)++;
              for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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

  myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
  distWeights->putScalar(0.);//All tentatively assigned vertices are now definitive

  if ( printFlag < 7) { //FIXME
    int Nphase1_agg = nAggregates;
    int total_aggs;

    sumAll(graph.GetComm(), Nphase1_agg, total_aggs);

    if (myPid == 0) {
      printf("Aggregation(%s) : Phase 1 - nodes aggregated = %d \n",label.c_str(), total_phase_one_aggregated);
      printf("Aggregation(%s) : Phase 1 - total aggregates = %d\n", label.c_str(), total_aggs);
    }
    int i = nAggregates - Nphase1_agg;

    { int ii; sumAll(graph.GetComm(),i,ii); i = ii; }
    if ( myPid == 0 ) {
      printf("Aggregation(%s) : Phase 3 - additional aggregates = %d\n",label.c_str(), i);
    }
  }

  // Determine vertices that are not shared by setting Temp to all ones
  // and doing MueLu_NonUnique2NonUnique(..., ADD). This sums values of all
  // local copies associated with each Gid. Thus, sums > 1 are shared.

  RCP<Cthulhu::Vector<double> > temp_ = Cthulhu::VectorFactory<double>::Build(nonUniqueMap);

  RCP<Cthulhu::Vector<double> > tempOutput_ = Cthulhu::VectorFactory<double>::Build(nonUniqueMap);
  ArrayRCP<double> tempOutput = tempOutput_->getDataNonConst(0);

  temp_->putScalar(1.);  
  tempOutput_->putScalar(0.); 

  myWidget.MueLu_NonUnique2NonUnique(*temp_, *tempOutput_, Cthulhu::ADD);
   
  std::vector<bool> gidNotShared(exp_nRows);
  for (int i = 0; i < exp_nRows; i++) {
    if (tempOutput[i] > 1.) gidNotShared[i] = false; 
    else  gidNotShared[i] = true; 
  }

  // Phase 4. 
double nAggregatesTarget; int    nAggregatesGlobal, minNAggs, maxNAggs; nAggregatesTarget = ((double)  uniqueMap->getGlobalNumElements())* (((double) uniqueMap->getGlobalNumElements())/ ((double) graph.GetGlobalNumEdges())); sumAll(graph.GetComm(), nAggregates, nAggregatesGlobal); minAll(graph.GetComm(), nAggregates, minNAggs); maxAll(graph.GetComm(), nAggregates, maxNAggs);

  //
  // Only do this phase if things look really bad. THIS
  // CODE IS PRETTY EXPERIMENTAL 
  //
#define MUELU_PHASE4BUCKETS 6
  if ((nAggregatesGlobal < graph.GetComm()->getSize()) &&
      (2.5*nAggregatesGlobal < nAggregatesTarget) &&
      (minNAggs ==0) && (maxNAggs <= 1)) {

    //TODO TODO
//     Epetra_Util util;
//     util.SetSeed( (unsigned int) myPid*2 + (int) (11*rand()));
//     k = (int)ceil( (10.*myPid)/graph.GetComm()->getSize());
//     for (i = 0; i < k+7; i++) util.SetSeed( (unsigned int) util.RandomInt() );
//    temp_->SetSeed( (unsigned int) util.RandomInt() );
    temp_->randomize(); 

    ArrayRCP<double> temp = temp_->getDataNonConst(0);

    // build a list of candidate root nodes (vertices not adjacent
    // to aggregated vertices)

    int nCandidates = 0, nCandidatesGlobal;

    ArrayRCP<int> candidates = Teuchos::arcp<int>(nVertices+1);

    double priorThreshold = 0.;
    for (int kkk = 0 ; kkk < MUELU_PHASE4BUCKETS; kkk++) {

      {
        ArrayRCP<const int> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
        ArrayView<const int> vertex2AggIdView = vertex2AggId();
        MueLu_RootCandidates(nVertices, vertex2AggIdView, graph,
                             candidates, nCandidates, nCandidatesGlobal);
        // views on distributed vectors are freed here.
      }

      double nTargetNewGuys =  nAggregatesTarget - nAggregatesGlobal;
      double threshold      =  priorThreshold + (1. - priorThreshold)*nTargetNewGuys/(nCandidatesGlobal + .001);
   
      threshold = (threshold*(kkk+1.))/((double) MUELU_PHASE4BUCKETS);
      priorThreshold = threshold;

      {
        ArrayRCP<int>    vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
        ArrayRCP<double> weights      = distWeights->getDataNonConst(0);

        for (int k = 0; k < nCandidates ; k++ ) { 
          int i = candidates[k];                  
          if ((vertex2AggId[i] == MUELU_UNAGGREGATED) && (fabs(temp[i])  < threshold)) {
            // Note: priorThreshold <= fabs(temp[i]) <= 1

            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

            if (neighOfINode.size() > minNodesPerAggregate) { //TODO: check if this test is exactly was we want to do
              int count = 0;
              for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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
                for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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
      myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
      distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive
      sumAll(graph.GetComm(), nAggregates, nAggregatesGlobal);

      // check that there are no aggregates sizes below minNodesPerAggregate
     
      aggregates.SetNumAggregates(nAggregates);

      MueLu_RemoveSmallAggs(aggregates, minNodesPerAggregate, distWeights, myWidget);
      
      nAggregates = aggregates.GetNumAggregates();
    }   // one possibility
  }

  // Initialize things for Phase 5. This includes building the transpose
  // of the matrix ONLY for transposed rows that correspond to unaggregted
  // ghost vertices. Further, the transpose is only a local transpose. 
  // Nonzero edges which exist on other processors are not represented.

  ArrayRCP<int> Mark = Teuchos::arcp<int>(exp_nRows+1);
  ArrayRCP<int> agg_incremented = Teuchos::arcp<int>(nAggregates+1);
  ArrayRCP<int> SumOfMarks = Teuchos::arcp<int>(nAggregates+1);

  for (int i = 0; i < exp_nRows; i++)   Mark[i] = MUELU_DISTONE_VERTEX_WEIGHT;
  for (int i = 0; i < nAggregates; i++) agg_incremented[i] = 0;
  for (int i = 0; i < nAggregates; i++) SumOfMarks[i] = 0;

  // Grab the transpose matrix graph for unaggregated ghost vertices.
  //     a) count the number of nonzeros per row in the transpose
  std::vector<int> RowPtr(exp_nRows+1-nVertices);
  //{
  ArrayRCP<const int> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);

  for (int i = nVertices; i < exp_nRows ;  i++) RowPtr[i-nVertices] = 0;
  for (int i = 0; i < nVertices;  i++) {

    // neighOfINode is the neighbor node list of node 'iNode'.
    ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

    for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
      int j = *it;
      if ( (j >= nVertices) && (vertex2AggId[j] == MUELU_UNAGGREGATED)){
        RowPtr[j-nVertices]++;
      }
    }
  }

  //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

  int iSum = 0, iTemp;
  for (int i = nVertices; i < exp_nRows ;  i++) {
    iTemp = RowPtr[i-nVertices];
    RowPtr[i-nVertices] = iSum;
    iSum += iTemp;
  }
  RowPtr[exp_nRows-nVertices] = iSum;
  std::vector<int> cols(iSum+1);
   
  //     c) Traverse matrix and insert entries in proper location.
  for (int i = 0; i < nVertices;  i++) {

    // neighOfINode is the neighbor node list of node 'iNode'.
    ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);
    
    for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
      int j = *it;
      if ( (j >= nVertices) && (vertex2AggId[j] == MUELU_UNAGGREGATED)){
        cols[RowPtr[j-nVertices]++] = i;
      }
    }
  }

  //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.
  for (int i = exp_nRows; i > nVertices;  i--)
    RowPtr[i-nVertices] = RowPtr[i-1-nVertices];
  RowPtr[0] = 0;
    
  // views on distributed vectors are freed here.
  vertex2AggId = Teuchos::null;
  //}

  int bestScoreCutoff;
  int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

  // Stick unaggregated vertices into existing aggregates as described above. 
  bool cannotLoseAllFriends; // Used to address possible loss of vertices in 
  // arbitration of shared nodes discussed above.
  for (int kk = 0; kk < 10; kk += 2) {
    bestScoreCutoff = thresholds[kk];
    for (int i = 0; i < exp_nRows; i++) {

      ArrayRCP<int> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
      ArrayRCP<const int> procWinner = aggregates.GetProcWinner()->getData(0);
      ArrayRCP<double> weights       = distWeights->getDataNonConst(0);

      if (vertex2AggId[i] == MUELU_UNAGGREGATED) {

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const int> neighOfINode;

        // Grab neighboring vertices which is either in graph for local ids
        // or sits in transposed fragment just constructed above for ghosts.
        if (i < nVertices) {
          neighOfINode = graph.getNeighborVertices(i);
        }
        else {
          int *rowi_col = NULL, rowi_N;
          rowi_col = &(cols[RowPtr[i-nVertices]]);
          rowi_N   = RowPtr[i+1-nVertices] - RowPtr[i-nVertices];
          
          neighOfINode = ArrayView<const int>(rowi_col, rowi_N);
        }
        for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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
        int best_agg; // TODO: init ?
        int BestMark;
        for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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
         for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int Adjacent    = *it;
          int AdjacentAgg = vertex2AggId[Adjacent];
          if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
        }
        // Tentatively assign vertex to best_agg. 
        if ( (best_score >= bestScoreCutoff) && (cannotLoseAllFriends)) { 
          vertex2AggId[i] = best_agg;
          weights[i] = best_score;
          agg_incremented[best_agg]++;
          Mark[i] = (int) ceil(   ((double) BestMark)/2.);
        }
      }

      // views on distributed vectors are freed here.
    }

    myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
    distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive
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

    ArrayRCP<int> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<double> weights       = distWeights->getDataNonConst(0);
    ArrayRCP<const int> procWinner = aggregates.GetProcWinner()->getData(0);

    int count = 0;
    for (int i = 0; i < nVertices; i++) { 
      if ((vertex2AggId[i] == MUELU_UNAGGREGATED) ) {
        Nleftover++;

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

        // We don't want too small of an aggregate. So lets see if there is an
        // unaggregated neighbor that we can also put with this vertex

        vertex2AggId[i] = nAggregates;
        weights[i] = 1.;
        if (count == 0) aggregates.SetIsRoot(i);
        count++;
        for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
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
      // Can stick small aggregate with 0th aggregate?
      if (nAggregates > 0) {
        for (int i = 0; i < nVertices; i++) {
          if ((vertex2AggId[i] == nAggregates) && (procWinner[i] == myPid)){
            vertex2AggId[i] = 0;
            aggregates.SetIsRoot(i,false);
          }
        }
      }
      else {
        Nsingle++;
        nAggregates++;
      }
    }

    // views on distributed vectors are freed here.
  }

  myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, false);

  if (printFlag < 7) { //FIXME
    { int total_Nsingle=0;   sumAll(graph.GetComm(), Nsingle, total_Nsingle);     Nsingle = total_Nsingle; }
    { int total_Nleftover=0; sumAll(graph.GetComm(), Nleftover, total_Nleftover); Nleftover = total_Nleftover; }
    int total_aggs;
    sumAll(graph.GetComm(), nAggregates, total_aggs);
    if ( myPid == 0 ) {
      printf("Aggregation(%s) : Phase 3 - total aggregates = %d\n",label.c_str(), total_aggs);
      printf("Aggregation(%s) : Phase 6 - leftovers = %d and singletons = %d\n",label.c_str() ,Nleftover, Nsingle);
    }
  }

  aggregates.SetNumAggregates(nAggregates);

  return 0;
} //AggregateLeftOvers

int UCAggregationFactory::MueLu_RootCandidates(int nVertices, ArrayView<const int> & vertex2AggId, const Graph graph,
                         ArrayRCP<int> &candidates, int &nCandidates, int &nCandidatesGlobal) const
{
  nCandidates = 0;
 
  for (int i = 0; i < nVertices; i++ ) {
    if (vertex2AggId[i] == MUELU_UNAGGREGATED) {
      bool noAggdNeighbors = true;

      // neighOfINode is the neighbor node list of node 'iNode'.
      ArrayView<const int> neighOfINode = graph.getNeighborVertices(i);

      for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        int adjacent    = *it;
        if (vertex2AggId[adjacent] != MUELU_UNAGGREGATED) 
          noAggdNeighbors = false;
      }
      if (noAggdNeighbors == true) candidates[nCandidates++] = i;
    }
  }

  sumAll(graph.GetComm(), nCandidates, nCandidatesGlobal);

  return 0;
} //MueLu_RootCandidates

int UCAggregationFactory::MueLu_ComputeAggSizes(Aggregates & aggregates, ArrayRCP<int> & aggSizes) const
{
  int myPid = aggregates.GetMap()->getComm()->getRank();

  int nAggregates = aggregates.GetNumAggregates();

  ArrayRCP<int> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  ArrayRCP<int> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  int size = procWinner.size();

  for (int i = 0; i < nAggregates; i++) aggSizes[i] = 0;
  for (int k = 0; k < size; k++ ) {
    if (procWinner[k] == myPid) aggSizes[vertex2AggId[k]]++;
  }

  return 0; //TODO
} //MueLu_ComputeAggSizes

int UCAggregationFactory::MueLu_RemoveSmallAggs(Aggregates& aggregates, int min_size,
                          RCP<Cthulhu::Vector<double> > & distWeights, const AggAlgorithm2Comm & myWidget) const
{
  int myPid = aggregates.GetMap()->getComm()->getRank();
  
  int nAggregates = aggregates.GetNumAggregates();

  ArrayRCP<int> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
  ArrayRCP<int> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  int size = procWinner.size();

  ArrayRCP<int> AggInfo = Teuchos::arcp<int>(nAggregates+1);

  ArrayRCP<double> weights = distWeights->getDataNonConst(0);

  MueLu_ComputeAggSizes(aggregates, AggInfo);

  // Make a list of all aggregates indicating New AggId
  // Use AggInfo array for this.

  int NewNAggs = 0; 
  for (int i = 0; i < nAggregates; i++) {
    if ( AggInfo[i] < min_size) { 
      AggInfo[i] =  MUELU_UNAGGREGATED;
    }
    else AggInfo[i] = NewNAggs++;
  }

  for (int k = 0; k < size; k++ ) {
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

  myWidget.MueLu_ArbitrateAndCommunicate(*distWeights, aggregates, true);
  distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

  // procWinner is not set correctly for aggregates which have 
  // been eliminated
  for (int i = 0; i < size; i++) {
    if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
      procWinner[i] = MUELU_UNASSIGNED;
  }
  aggregates.SetNumAggregates(nAggregates);

  return 0; //TODO
} //MueLu_RemoveSmallAggs

int UCAggregationFactory::MueLu_RandomReorder(Teuchos::ArrayRCP<int> list, const Map &map) const
{

  TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "RandomReorder: TODO");

//   Epetra_Vector     RandVec(map);
//   Epetra_IntVector iRandVec(map);

//   double *ptr;
//   int  *iptr;

//   RandVec.Random(); RandVec.ExtractView(&ptr); iRandVec.ExtractView(&iptr);
//   for (int i=0; i <  map.NumMyElements(); ++i) iptr[i] = (int) (10000.*ptr[i]);
//   Epetra_Util::Sort(true,RandVec.getMap().NumMyElements(), iptr, 0,NULL,1,&list);

  return 0; 
} //MueLu_RandomReorder
