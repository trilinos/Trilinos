#ifndef MUELU_AGGREGATES_DECL_HPP
#define MUELU_AGGREGATES_DECL_HPP

#include <Teuchos_Describable.hpp>

#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Utilities.hpp"

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

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class Aggregates : public BaseClass {

#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
     
    Aggregates(const Graph & graph);
    virtual ~Aggregates() ;
     
    inline LO GetNumAggregates() const           ; // rename GetNumLocal ?
    inline void SetNumAggregates(LO nAggregates) ;
    inline RCP<LOVector> & GetVertex2AggId()     ; // LocalOrdinal because it's an array of local id
    inline RCP<LOVector> & GetProcWinner()       ;
    inline bool IsRoot(LO i) const               ; // Local
    inline void SetIsRoot(LO i, bool value=true) ; // Local
     
    inline const RCP<const Xpetra::Map<LO,GO> > GetMap() const ;

    /*! @brief Compute sizes of all the aggregates.

    - FIXME Is this dangerous, i.e., could the user change this?
    */
    Teuchos::ArrayRCP<LO> ComputeAggregateSizes() const
    ; //ComputeAggSizes

    /*! @brief Compute lookup table that provides DOFs belonging to a given aggregate.

    @param aggToRowMap aggToRowMap[i][j] is the jth local DOF in local aggregate i
    */
    void ComputeAggregateToRowMap(Teuchos::ArrayRCP<Teuchos::ArrayRCP<LO> > &aggToRowMap) const ; //AggregateToRowMap

    //! @name Overridden from Teuchos::Describable 
    //@{
     
    //! Return a simple one-line description of this object.
    std::string description() const ;
     
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const ;; //class Aggregates

    // Constructors to create aggregates.
    template <class LocalOrdinal ,
              class GlobalOrdinal,
              class Node         ,
              class LocalMatOps>
    Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Aggregates(const MueLu::Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & graph)
    ;

  } //namespace MueLu

#define MUELU_AGGREGATES_SHORT
#endif // MUELU_AGGREGATES_DECL_HPP
