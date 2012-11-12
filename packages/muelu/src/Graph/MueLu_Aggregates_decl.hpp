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

    /*! @brief Standard constructor for Aggregates structure
     *
     * Standard constructor of aggregates takes a Graph object as parameter.
     * Uses the graph.GetImportMap() to initialize the internal vector for mapping nodes to (local) aggregate ids as well as
     * th mapping of node to the owning processor id.
     *
     */
    Aggregates(const Graph & graph);

    /*! @brief Constructor for Aggregates structure
     *
     * This constructor takes a RCP pointer to a map which is used for the internal mappings of nodes to the (local) aggregate ids and the owning processor.
     *
     */
    Aggregates(const RCP<const Map> & map);

    /*! @brief Destructor
     *
     */
    virtual ~Aggregates() { }

    LO GetNumAggregates() const           { return nAggregates_;        } ///< returns the number of aggregates of the current processor. Note: could/should be renamed to GetNumLocalAggregates?
    void SetNumAggregates(LO nAggregates) { nAggregates_ = nAggregates; } ///< set number of local aggregates on current processor. This has to be done by the aggregation routines.
    RCP<LOVector> & GetVertex2AggIdNonConst()     { return vertex2AggId_;       } ///< provide access to the underlaying aggregation information. Vertex2AggId provides a mapping of the local node id to the local aggregate id on the current processor
    RCP<LOVector> & GetProcWinnerNonConst()       { return procWinner_;         } ///< provide access to the underlaying aggregation information. ProcWinner stores the processor id owning the given (local) node id.
    const RCP<LOVector> & GetVertex2AggId() const { return vertex2AggId_;       } ///< provide access to the underlaying aggregation information. Vertex2AggId provides a mapping of the local node id to the local aggregate id on the current processor
    const RCP<LOVector> & GetProcWinner() const   { return procWinner_;         } ///< provide access to the underlaying aggregation information. ProcWinner stores the processor id owning the given (local) node id.

    bool IsRoot(LO i) const               { return isRoot_[i];          } ///< returns true if node with given local node id is marked to be a root node
    void SetIsRoot(LO i, bool value=true) { isRoot_[i] = value;         } ///< set root node information. Used by aggregation methods only.

    const RCP<const Map> GetMap() const; ///< returns (overlapping) map of aggregate/node distribution

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
