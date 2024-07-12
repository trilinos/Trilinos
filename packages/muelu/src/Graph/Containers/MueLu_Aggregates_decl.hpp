// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATES_DECL_HPP
#define MUELU_AGGREGATES_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Kokkos_StaticCrsGraph.hpp>

#include "MueLu_Aggregates_fwd.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_BaseClass.hpp"

#include "MueLu_LWGraph_kokkos.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_IndexManager.hpp"
#include "MueLu_IndexManager_kokkos.hpp"

#define MUELU_UNAGGREGATED -1 /* indicates that a node is unassigned to  */
                              /* any aggregate.                          */

#define MUELU_UNASSIGNED -1 /* indicates a vertex is not yet claimed   */
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
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class Aggregates : public BaseClass {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using execution_space     = typename Node::execution_space;
  using node_type           = Node;
  using device_type         = typename Node::device_type;
  using range_type          = Kokkos::RangePolicy<local_ordinal_type, execution_space>;
  using LO_view             = Kokkos::View<local_ordinal_type*, device_type>;

  using aggregates_sizes_type = Kokkos::View<LocalOrdinal*, device_type>;

 private:
#undef MUELU_AGGREGATES_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  // Defining types that require the short names included above
  using local_graph_type = typename LWGraph_kokkos::local_graph_type;
  using colors_view_type = Kokkos::View<typename local_graph_type::entries_type::data_type,
                                        typename local_graph_type::device_type::memory_space>;

  /*! @brief Standard constructor for Aggregates structure
   *
   * Standard constructor of aggregates takes a Graph object as parameter.
   * Uses the graph.GetImportMap() to initialize the internal vector for mapping nodes to (local) aggregate ids as well as
   * the mapping of node to the owning processor id.
   *
   */
  Aggregates(const LWGraph& graph);

  /*! @brief Standard constructor for Aggregates structure
   *
   * Standard constructor of aggregates takes a LWGraph object as parameter.
   * Uses the graph.GetImportMap() to initialize the internal vector for mapping nodes to (local) aggregate ids as well as
   * the mapping of node to the owning processor id.
   *
   */
  Aggregates(LWGraph_kokkos graph);

  /*! @brief Constructor for Aggregates structure
   *
   * This constructor takes a RCP pointer to a map which is used for the internal mappings of nodes to the (local) aggregate ids and the owning processor.
   *
   */
  Aggregates(const RCP<const Map>& map);

  /*! @brief Destructor
   *
   */
  virtual ~Aggregates() {}

  //! @name Set/Get Methods for specific aggregation data
  //@{

  /*! @brief Get the index manager used by structured aggregation algorithms.
      This has to be done by the aggregation factory.
  */
  RCP<IndexManager_kokkos>& GetIndexManagerKokkos() { return geoDataKokkos_; }

  /*! @brief Set the index manager used by structured aggregation algorithms.
      This has to be done by the aggregation factory.
  */
  void SetIndexManagerKokkos(RCP<IndexManager_kokkos>& geoDataKokkos) { geoDataKokkos_ = geoDataKokkos; }

  /*! @brief Get the index manager used by various aggregation algorithms.
      This has to be done by the aggregation factory.
  */
  RCP<IndexManager>& GetIndexManager() { return geoData_; }

  /*! @brief Set the index manager used by various aggregation algorithms.
      This has to be done by the aggregation factory.
  */
  void SetIndexManager(RCP<IndexManager>& geoData) { geoData_ = geoData; }

  /*! @brief Get a distance 2 coloring of the underlying graph.
      The coloring is computed and set during Phase1 of aggregation.
  */
  colors_view_type& GetGraphColors() { return graphColors_; }

  /*! @brief Set a distance 2 coloring of the underlying graph.
      The coloring is computed and set during Phase1 of aggregation.
  */
  void SetGraphColors(colors_view_type graphColors) { graphColors_ = graphColors; }

  /*! @brief Get the number of colors needed by the distance 2 coloring.
   */
  LO GetGraphNumColors() { return graphNumColors_; }

  /*! @brief Set the number of colors needed by the distance 2 coloring.
   */
  void SetGraphNumColors(const LO graphNumColors) { graphNumColors_ = graphNumColors; }

  //@}

  /*! @brief Set number of local aggregates on current processor.

      This has to be done by the aggregation routines.
  */
  void SetNumAggregates(LO nAggregates) { numAggregates_ = nAggregates; }

  /*! @brief Set number of global aggregates on current processor.

      This has to be done by the aggregation routines.
  */
  void SetNumGlobalAggregates(GO nGlobalAggregates) { numGlobalAggregates_ = nGlobalAggregates; }

  ///< returns the number of aggregates of the current processor. Note: could/should be renamed to GetNumLocalAggregates?
  KOKKOS_INLINE_FUNCTION LO GetNumAggregates() const {
    return numAggregates_;
  }

  //! @brief Record whether aggregates include DOFs from other processes.
  KOKKOS_INLINE_FUNCTION void AggregatesCrossProcessors(const bool& flag) {
    aggregatesIncludeGhosts_ = flag;
  }

  /*! @brief Return false if and only if no aggregates include DOFs from other processes.

      Used in construction of tentative prolongator to skip a communication phase.
  */
  KOKKOS_INLINE_FUNCTION bool AggregatesCrossProcessors() const {
    return aggregatesIncludeGhosts_;
  }

  /*! @brief Returns a nonconstant vector that maps local node IDs to local aggregates IDs.

      For local node ID i, the corresponding vector entry v[i] is the local aggregate id to which i belongs on the current processor.
  */
  RCP<LOMultiVector>& GetVertex2AggIdNonConst() { return vertex2AggId_; }

  /*! @brief Returns nonconstant vector that maps local node IDs to owning processor IDs.

      For local node ID i, the corresponding vector entry v[i] is the owning processor ID.
  */
  RCP<LOVector>& GetProcWinnerNonConst() { return procWinner_; }
  /*! @brief Returns constant vector that maps local node IDs to local aggregates IDs.

      For local node ID i, the corresponding vector entry v[i] is the local aggregate id to which i belongs on the current processor.
  */
  const RCP<LOMultiVector>& GetVertex2AggId() const { return vertex2AggId_; }

  /*! @brief Returns constant vector that maps local node IDs to owning processor IDs.

      For local node ID i, the corresponding vector entry v[i] is the owning processor ID.
  */
  const RCP<LOVector>& GetProcWinner() const { return procWinner_; }

  //! Returns true if node with given local node id is marked to be a root node
  inline bool IsRoot(LO i) const { return isRoot_[i]; }

  /*! @brief Set root node information.

  Used by aggregation methods only.
  */
  inline void SetIsRoot(LO i, bool value = true) { isRoot_[i] = value; }

  const RCP<const Map> GetMap() const;  ///< returns (overlapping) map of aggregate/node distribution

  /*! @brief Compute sizes of aggregates

    Returns the number of nodes in each aggregate in an array.
    If the aggregate sizes are not stored internally (which is the default), they are computed and returned.
    If the aggregate sizes have been stored internally, then they are *not* recomputed, but instead the
    stored sizes are returned.

    @param[in] forceRecompute if true, force recomputation of the aggregate sizes.
   */
  typename aggregates_sizes_type::const_type ComputeAggregateSizes(bool forceRecompute = false) const;

  /*! @brief Compute sizes of aggregates

    Returns the number of nodes in each aggregate in an array.
    If the aggregate sizes are not stored internally (which is the default), they are computed and returned.
    If the aggregate sizes have been stored internally, then they are *not* recomputed, but instead the
    stored sizes are returned.

    @param[in] forceRecompute if true, force recomputation of the aggregate sizes.
   */
  Teuchos::ArrayRCP<LocalOrdinal> ComputeAggregateSizesArrayRCP(bool forceRecompute = false) const;

  local_graph_type GetGraph() const;

  /*! @brief Generates a compressed list of nodes in each aggregate, where
    the entries in aggNodes[aggPtr[i]] up to aggNodes[aggPtr[i+1]-1] contain the nodes in aggregate i.
    unaggregated contains the list of nodes which are, for whatever reason, not aggregated (e.g. Dirichlet)
   */
  void ComputeNodesInAggregate(LO_view& aggPtr, LO_view& aggNodes, LO_view& unaggregated) const;

  //! Get global number of aggregates
  //  If # of global aggregates is unknown, this method does coummunication and internally record the value
  GO GetNumGlobalAggregatesComputeIfNeeded();

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const;

 private:
  LO numAggregates_;        ///< Number of aggregates on this processor
  GO numGlobalAggregates_;  ///< Number of global aggregates

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
  RCP<IndexManager_kokkos> geoDataKokkos_;

  /*! geoData stores an index manager object that is used to perform structured aggreation
   *  on a problem.
   */
  RCP<IndexManager> geoData_;

  /*! graphColors_ stores a view that assigns a color to each node in the graph
   *  These colors are used to parallelize the aggregation process in UncoupledAggregation
   */
  colors_view_type graphColors_;

  /*! graphNumColors_ stores the number of colors that are needed to perform a distance 2
   *  coloring of the underlying graph.
   */
  LO graphNumColors_;

  //! An ArrayRCP of booleans specifying if a local entry is an aggregate root.
  Teuchos::ArrayRCP<bool> isRoot_;

  //! Set to false iff aggregates do not include any DOFs belong to other processes.
  bool aggregatesIncludeGhosts_;

  //! Array of sizes of each local aggregate.
  mutable aggregates_sizes_type aggregateSizes_;

  /*! aggragateSizesHost_ is a host copy of aggregate sizes, which
   * helps slightly reduce the cost of calling ComputeAggregateSizes
   * from different parts of MueLu that require such data on the host device.
   */
  mutable
      typename aggregates_sizes_type::HostMirror aggregateSizesHost_;

  //! Aggregates represented as Kokkos graph type
  mutable local_graph_type graph_;
};

}  // namespace MueLu

#define MUELU_AGGREGATES_SHORT
#endif  // MUELU_AGGREGATES_DECL_HPP
