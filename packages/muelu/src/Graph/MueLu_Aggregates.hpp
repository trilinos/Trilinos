#ifndef MUELU_AGGREGATES_HPP
#define MUELU_AGGREGATES_HPP

#include <Teuchos_Describable.hpp>

#include "MueLu_Graph.hpp"

#define MUELOO_AGGR_READY    -11  /* indicates that a node is available to be*/
                                  /* selected as a root node of an aggregate */
#define MUELOO_AGGR_NOTSEL   -12  /* indicates that a node has been rejected */
                                  /* as a root node. This could perhaps be   */
                                  /* because if this node had been selected a*/
                                  /* small aggregate would have resulted.    */
#define MUELOO_AGGR_SELECTED -13  /* indicates that a node has been assigned */
                                  /* to an aggregate.                        */
#define MUELOO_UNAGGREGATED  -1   /* indicates that a node is unassigned to  */
                                  /* any aggregate.                          */
#define MUELOO_UNASSIGNED    -1   /* indicates a vertex is not yet claimed   */
                                  /* by a processor during aggregation.      */
                                  /* Note, it is possible at                 */
                                  /* this stage that some processors may have*/
                                  /* claimed their copy of a vertex for one  */
                                  /* of their aggregates.  However, some     */
                                  /* arbitration still needs to occur.       */
                                  /* The corresponding procWinner[]'s remain */
                                  /* as MUELOO_UNASSIGNED until              */
                                  /* MueLu_ArbitrateAndCommunicate() is     */
                                  /* invoked to arbitrate.                   */
#define MUELOO_NOSCORE       -100 /* indicates that a quality score has not  */
                                  /* yet been assigned when determining to   */
                                  /* which existing aggregate a vertex       */
                                  /* should be assigned.                     */
#define MUELOO_DISTONE_VERTEX_WEIGHT 100  /* Weights associated with all     */
                                  /* vertices that have a direct connection  */
                                  /* to the aggregate root.                  */
#define INCR_SCALING 3            /* Determines how much of a penalty should */
                                  /* be deduced from a score during Phase 5  */
                                  /* for each Phase 5 vertex already added   */
                                  /* to this aggregate. Specifically the     */
                                  /* penalty associated with aggregate y is  */
                                  /*   max (INCR_SCALING*NNewVtx,            */
                                  /*        UnpenalizedScore*(1-             */
                                  /*              MUELOO_PENALTYFACTOR))*/
                                  /* where NNewVtx is the number of phase 5  */
                                  /* vertices already assigned to y.         */
#define MUELOO_PENALTYFACTOR .30 /* determines maximum allowable        */
                                  /* percentage of a score that can be       */
                                  /* deducted from this score for having     */
                                  /* already enlargened an aggregate to      */
                                  /* which we are contemplated adding another*/
                                  /* vertex.  Should be between 0 and 1.     */

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

  class Aggregates : public Teuchos::Describable {
    
  public:
    
    Aggregates(Graph *Graph, const std::string & objectLabel);
    ~Aggregates();

    inline int GetNumAggregates()                      { return nAggregates_;            }
    inline void SetNumAggregates(int nAggregates)      { nAggregates_ = nAggregates;     }
    inline Epetra_IntVector* GetVertex2AggId()         { return vertex2AggId_;           }
    inline Epetra_IntVector* GetProcWinner() { return procWinner_;             }
    inline bool IsRoot(int i)                          { return isRoot_[i];              }
    inline void SetIsRoot(int vertex, bool value=true) { isRoot_[vertex] = value; }

  private:
    int   nAggregates_;             /* Number of aggregates on this processor  */
    
    Epetra_IntVector *vertex2AggId_;/* vertex2AggId[k] gives a local id        */
                                    /* corresponding to the aggregate to which */
                                    /* local id k has been assigned.  While k  */
    Epetra_IntVector *procWinner_;  /* is the local id on my processor (MyPID),*/
                                    /* vertex2AggId[k] is the local id on the  */
                                    /* processor which actually owns the       */
                                    /* aggregate. This owning processor has id */
                                    /* given by procWinner[k].                 */

    bool *isRoot_;                  /* IsRoot[i] indicates whether vertex i  */
                                    /* is a root node.                       */

    // Epetra_CrsGraph *agg2Vertex_;/* Currently not used                    */

  };

  // Constructors to create aggregates.
  Aggregates::Aggregates(Graph *graph, const std::string & objectLabel)
  {
    
    setObjectLabel(objectLabel);
    
    nAggregates_  = 0;
    
    vertex2AggId_ = new Epetra_IntVector(graph->GetImportMap());
    vertex2AggId_->PutValue(MUELOO_UNAGGREGATED);
    
    procWinner_ = new Epetra_IntVector(graph->GetImportMap());
    procWinner_->PutValue(MUELOO_UNASSIGNED);
    
    isRoot_ = new bool[graph->GetImportMap().NumMyElements()];
    for (int i=0; i < graph->GetImportMap().NumMyElements(); i++)
      isRoot_[i] = false;

    //agg2Vertex_   = NULL; /* Currently not used */
  }

  // Destructor for aggregates.
  Aggregates::~Aggregates()
  {
    if (isRoot_       != NULL) delete [] (isRoot_);
    if (procWinner_   != NULL) delete procWinner_;
    if (vertex2AggId_ != NULL) delete vertex2AggId_;
    // if (agg2Vertex_   != NULL) delete agg2Vertex_; /* Currently not used */
  }

}

#endif
