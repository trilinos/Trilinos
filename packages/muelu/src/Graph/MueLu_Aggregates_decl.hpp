#ifndef MUELU_AGGREGATES_DECL_HPP
#define MUELU_AGGREGATES_DECL_HPP

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Aggregates_fwd.hpp"

#include "MueLu_Graph_fwd.hpp"

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
****************************************************************************/

namespace MueLu {

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class Aggregates : public BaseClass {
#undef MUELU_AGGREGATES_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
     
    Aggregates(const Graph & graph);
    Aggregates(Aggregates & a); // special copy constructor (that skips amalgamation information)
    virtual ~Aggregates() { }
     
    LO GetNumAggregates() const           { return nAggregates_;        } // rename GetNumLocal ?
    void SetNumAggregates(LO nAggregates) { nAggregates_ = nAggregates; }
    RCP<LOVector> & GetVertex2AggIdNonConst()     { return vertex2AggId_;       } // LocalOrdinal because it's an array of local id
    RCP<LOVector> & GetProcWinnerNonConst()       { return procWinner_;         }
    const RCP<LOVector> & GetVertex2AggId() const { return vertex2AggId_;       } // LocalOrdinal because it's an array of local id
    const RCP<LOVector> & GetProcWinner() const   { return procWinner_;         }

    bool IsRoot(LO i) const               { return isRoot_[i];          } // Local
    void SetIsRoot(LO i, bool value=true) { isRoot_[i] = value;         } // Local
    
    const RCP<const Map> GetMap() const { return vertex2AggId_->getMap(); }

    /*! @brief Compute sizes of aggregates
     * 
     * returns the number of nodes in each aggregate in an array.
     * 
     */
    Teuchos::ArrayRCP<LO> ComputeAggregateSizes() const;

    //! @name Overridden from Teuchos::Describable 
    //@{
     
    //! Return a simple one-line description of this object.
    std::string description() const;
     
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const;

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

    //! Get global number of aggregates
    // This method is private because it is used only for printing and because with the current implementation, communication occurs each time this method is called.
    GO GetNumGlobalAggregates() const;
  };

} //namespace MueLu

#define MUELU_AGGREGATES_SHORT
#endif // MUELU_AGGREGATES_DECL_HPP
