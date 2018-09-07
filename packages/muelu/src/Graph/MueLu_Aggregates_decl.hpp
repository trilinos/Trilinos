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
#ifndef MUELU_AGGREGATES_DECL_HPP
#define MUELU_AGGREGATES_DECL_HPP

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Aggregates_fwd.hpp"

#include "MueLu_Graph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_IndexManager.hpp"

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

****************************************************************************/

namespace MueLu {

/*!
    @class Aggregates
    @brief Container class for aggregation information.

    @ingroup Aggregation

    Structure holding aggregate information. Right now, nAggregates, IsRoot,
    Vertex2AggId, procWinner are populated.  This allows us to look at a node
    and determine the aggregate to which it has been assigned and the id of the
    processor that owns this aggregate. It is not so easy to determine vertices
    within the kth aggregate or the size of the kth aggregate. Thus, it might be
    useful to have a secondary structure which would be a rectangular CrsGraph
    where rows (or vertices) correspond to aggregates and colunmns (or edges)
    correspond to nodes. While not strictly necessary, it might be convenient.
*/

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class Aggregates : public BaseClass {
#undef MUELU_AGGREGATES_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    /*! @brief Standard constructor for Aggregates structure
     *
     * Standard constructor of aggregates takes a Graph object as parameter.
     * Uses the graph.GetImportMap() to initialize the internal vector for mapping nodes to (local) aggregate ids as well as
     * the mapping of node to the owning processor id.
     *
     */
    Aggregates(const GraphBase & graph);

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

    /*! @brief Set number of local aggregates on current processor.

        This has to be done by the aggregation routines.
    */
    void SetNumAggregates(LO nAggregates) { nAggregates_ = nAggregates; }

    /*! @brief Get the index manager used by structured aggregation algorithms.

        This has to be done by the aggregation factory.
    */
    RCP<IndexManager>& GetIndexManager() { return geoData_; }

    /*! @brief Get the index manager used by structured aggregation algorithms.

        This has to be done by the aggregation factory.
    */
    void SetIndexManager(RCP<IndexManager> & geoData) { geoData_ = geoData; }

    //! @brief Record whether aggregates include DOFs from other processes.
    void AggregatesCrossProcessors(const bool &flag) {aggregatesIncludeGhosts_ = flag;};

    /*! @brief Return false if and only if no aggregates include DOFs from other processes.

        Used in construction of tentative prolongator to skip a communication phase.
    */
    bool AggregatesCrossProcessors() const {return aggregatesIncludeGhosts_;};

    /*! @brief Returns a nonconstant vector that maps local node IDs to local aggregates IDs.

        For local node ID i, the corresponding vector entry v[i] is the local aggregate id to which i belongs on the current processor.
    */
    RCP<LOMultiVector> & GetVertex2AggIdNonConst()     { return vertex2AggId_;       }

    /*! @brief Returns nonconsant vector that maps local node IDs to owning processor IDs.

        For local node ID i, the corresponding vector entry v[i] is the owning processor ID.
    */
    RCP<LOVector> & GetProcWinnerNonConst()       { return procWinner_;         }
    /*! @brief Returns constant vector that maps local node IDs to local aggregates IDs.

        For local node ID i, the corresponding vector entry v[i] is the local aggregate id to which i belongs on the current processor.
    */
    const RCP<LOMultiVector> & GetVertex2AggId() const { return vertex2AggId_;       }

    /*! @brief Returns constant vector that maps local node IDs to owning processor IDs.

        For local node ID i, the corresponding vector entry v[i] is the owning processor ID.
    */
    const RCP<LOVector> & GetProcWinner() const   { return procWinner_;         }

    //! Returns true if node with given local node id is marked to be a root node
    bool IsRoot(LO i) const               { return isRoot_[i];          }

    /*! @brief Set root node information.

    Used by aggregation methods only.
    */
    void SetIsRoot(LO i, bool value=true) { isRoot_[i] = value;         }

    const RCP<const Map> GetMap() const; ///< returns (overlapping) map of aggregate/node distribution

    /*! @brief Compute sizes of aggregates

      Returns the number of nodes in each aggregate in an array.
      If the aggregate sizes are not stored internally (which is the default), they are computed and returned.
      If the aggregate sizes have been stored internally, then they are *not* recomputed, but instead the
      stored sizes are returned.

      @param[in] forceRecompute if true, force recomputation of the aggregate sizes.
     */
    Teuchos::ArrayRCP<LO> ComputeAggregateSizes(bool forceRecompute = false) const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const;

  private:
    LO   nAggregates_;              ///< Number of aggregates on this processor

    /*! vertex2AggId[k] gives a local id corresponding to the aggregate to which
     * local id k has been assigned. While k is the local id on my processor (MyPID),
     * vertex2AggId[k] is the local id on the processor which actually owns the aggregate.
     */
    RCP<LOMultiVector> vertex2AggId_;

    /*!
     * If k is the local id on my processor (MyPID), the owning processor has the
     * id given by procWinner[k]
     */
    RCP<LOVector> procWinner_;

    /*! geoData stores an index manager object that is used to perform structured aggreation
     *  on a problem.
     */
    RCP<IndexManager> geoData_;

    Teuchos::ArrayRCP<bool> isRoot_;//< IsRoot[i] indicates whether vertex i is a root node.

    //! Set to false iff aggregates do not include any DOFs belong to other processes.
    bool aggregatesIncludeGhosts_;

    //! Array of sizes of each local aggregate.
    mutable Teuchos::ArrayRCP<LO> aggregateSizes_;

    //! Get global number of aggregates
    // This method is private because it is used only for printing and because with the current implementation, communication occurs each time this method is called.
    GO GetNumGlobalAggregates() const;
  };

} //namespace MueLu

#define MUELU_AGGREGATES_SHORT
#endif // MUELU_AGGREGATES_DECL_HPP
