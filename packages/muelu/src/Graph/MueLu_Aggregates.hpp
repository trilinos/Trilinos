#ifndef MUELU_AGGREGATES_HPP
#define MUELU_AGGREGATES_HPP

#include <Teuchos_Describable.hpp>

#include "Xpetra_VectorFactory.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Graph.hpp"

#define MUELU_UNAGGREGATED  -1   /* indicates that a node is unassigned to  */
                                 /* any aggregate.                          */

#define MUELU_UNASSIGNED    -1   /* indicates a vertex is not yet claimed   */
                                 /* by a processor during aggregation.      */
                                 /* Note, it is possible at                 */
                                 /* this stage that some processors may have*/
                                 /* claimed their copy of a vertex for one  */
                                 /* of their aggregates.  However, some     */
                                 /* arbitration still needs to occur.       */
                                 /* The corresponding procWinner[]'s remain */
                                 /* as MUELU_UNASSIGNED until               */
                                 /* ArbitrateAndCommunicate() is            */
                                 /* invoked to arbitrate.                   */

/***************************************************************************** 
   Structure holding aggregate information. Right now, nAggregates, IsRoot,
   Vertex2AggId, procWinner are populated.  This allows us to look at a node
   and determine the aggregate to which it has been assigned and the id of the 
   processor that owns this aggregate. It is not so easy to determine vertices
   within the kth aggregate or the size of the kth aggregate. Thus, it might be
   useful to have a secondary structure which would be a rectangular CrsGraph 
   where rows (or vertices) correspond to aggregates and colunmns (or edges) 
   correspond to nodes. While not strictly necessary, it might be convenient.
 *****************************************************************************/

namespace MueLu {

   template <class LocalOrdinal  = int, 
             class GlobalOrdinal = LocalOrdinal, 
             class Node          = Kokkos::DefaultNode::DefaultNodeType, 
             class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
   class Aggregates : public Teuchos::Describable {

#include "MueLu_UseShortNamesOrdinal.hpp"

   public:
     
     Aggregates(const Graph & graph, const std::string & objectLabel = "");
     virtual ~Aggregates() {}
     
     inline LO GetNumAggregates()                 { return nAggregates_;         } // rename GetNumLocal ?
     inline void SetNumAggregates(LO nAggregates) { nAggregates_ = nAggregates;  }
     inline RCP<LOVector> & GetVertex2AggId()     { return vertex2AggId_;        } // LocalOrdinal because it's an array of local id
     inline RCP<LOVector> & GetProcWinner()       { return procWinner_;          }
     inline bool IsRoot(LO i)                     { return isRoot_[i];           } // Local
     inline void SetIsRoot(LO i, bool value=true) { isRoot_[i] = value;          } // Local
     
     inline const RCP<const Xpetra::Map<LO,GO> > GetMap() { return GetVertex2AggId()->getMap(); }

     /*! @brief Compute sizes of all the aggregates.

        - FIXME Is this dangerous, i.e., could the user change this?
     */
     Teuchos::ArrayRCP<LO> ComputeAggregateSizes() const
     {
       if (aggregateSizes_ == Teuchos::null)
       {
         aggregateSizes_ = Teuchos::ArrayRCP<LO>(nAggregates_);
         int myPid = vertex2AggId_->getMap()->getComm()->getRank();
         Teuchos::ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
         Teuchos::ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);
         LO size = procWinner.size();

         for (LO i = 0; i < nAggregates_; ++i) aggregateSizes_[i] = 0;
         for (LO k = 0; k < size; ++k ) {
           if (procWinner[k] == myPid) aggregateSizes_[vertex2AggId[k]]++;
         }
       }

       return aggregateSizes_;
     } //ComputeAggSizes

     /*! @brief Compute lookup table that provides DOFs belonging to a given aggregate.

         @param aggToRowMap aggToRowMap[i][j] is the jth local DOF in local aggregate i
     */
     void ComputeAggregateToRowMap(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LO> > &aggToRowMap) const {
       int myPid = vertex2AggId_->getMap()->getComm()->getRank();
       ArrayRCP<LO> procWinner   = procWinner_->getDataNonConst(0);
       ArrayRCP<LO> vertex2AggId = vertex2AggId_->getDataNonConst(0);

       ArrayRCP<LO> aggSizes = ComputeAggregateSizes();
       LO t=0;
       for (typename ArrayRCP<ArrayRCP<LO> >::iterator a2r=aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r)
         *a2r = ArrayRCP<LO>(aggSizes[t++]);
       ArrayRCP< LO > numDofs(nAggregates_,0);  //Track how many DOFS have been recorded so far
                                           //for each each aggregate in aggToRowMap.
       LO size = procWinner.size();
       for (LO k = 0; k < size; ++k ) {
         LO myAgg = vertex2AggId[k];
         if (procWinner[k] == myPid) {
           aggToRowMap[ myAgg ][ numDofs[myAgg] ] = k;
           ++(numDofs[myAgg]);
         }
       }

     } //AggregateToRowMap

   private:
    LO   nAggregates_;              /* Number of aggregates on this processor  */
    
    RCP<LOVector> vertex2AggId_;    /* vertex2AggId[k] gives a local id        */
                                    /* corresponding to the aggregate to which */
                                    /* local id k has been assigned.  While k  */
    RCP<LOVector> procWinner_;      /* is the local id on my processor (MyPID),*/
                                    /* vertex2AggId[k] is the local id on the  */
                                    /* processor which actually owns the       */
                                    /* aggregate. This owning processor has id */
                                    /* given by procWinner[k].                 */

    Teuchos::ArrayRCP<bool> isRoot_;/* IsRoot[i] indicates whether vertex i  */
                                    /* is a root node.                       */

    //! Array of sizes of each local aggregate.
    mutable Teuchos::ArrayRCP<LO> aggregateSizes_;

  }; //class Aggregates

  // Constructors to create aggregates.
   template <class LocalOrdinal ,
             class GlobalOrdinal,
             class Node         ,
            class LocalMatOps>
   Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Aggregates(const MueLu::Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & graph, const std::string & objectLabel)
  {
    
    setObjectLabel(objectLabel);
    
    nAggregates_  = 0;
    
    vertex2AggId_ = LOVectorFactory::Build(graph.GetImportMap());
    vertex2AggId_->putScalar(MUELU_UNAGGREGATED);
    
    procWinner_ = LOVectorFactory::Build(graph.GetImportMap());
    procWinner_->putScalar(MUELU_UNASSIGNED);
    
    isRoot_ = Teuchos::ArrayRCP<bool>(graph.GetImportMap()->getNodeNumElements());
    for (size_t i=0; i < graph.GetImportMap()->getNodeNumElements(); i++)
      isRoot_[i] = false;

  }

} //namespace MueLu

#define MUELU_AGGREGATES_SHORT
#endif //ifndef MUELU_AGGREGATES_HPP
