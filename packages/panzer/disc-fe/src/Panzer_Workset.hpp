// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_WORKSET_HPP
#define PANZER_WORKSET_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <iostream>

#include "Panzer_Dimension.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_Dimension.hpp"

#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

  struct WorksetNeeds;

  struct LocalMeshPartition;

  class SubcellConnectivity;

  class OrientationsInterface;

  /**
   * \class WorksetOptions
   *
   * \brief Used to define options for lazy evaluation of BasisValues and IntegrationValues objects
   */
  struct
  WorksetOptions
  {
    /// Default constructor
    WorksetOptions():
      side_assembly_(false),
      align_side_points_(false)
    {

    }

    /// Build integration values for sides
    bool side_assembly_;

    /// If workset side integration values must align with another workset, there must be a unique order assigned
    bool align_side_points_;

    /// Must be set to apply orientations - if it is set, then orientations will be applied to basis values
    Teuchos::RCP<const OrientationsInterface> orientations_;
  };

  /** This is used within the workset to make edge based (DG like) assembly
    * an easier task. This basically allows separation of the workset abstraction
    * from how it is accessed.
    */
  class WorksetDetails {
  public:
    typedef PHX::MDField<double,Cell,NODE,Dim> CellCoordArray;

    //! Default constructor
    WorksetDetails();

    //! Constructs the workset details from a given chunk of the mesh
    void
    setup(const LocalMeshPartition & partition,
          const WorksetOptions & options);

    /// DEPRECATED - use: numCells()
    int num_cells;

    /// DEPRECATED - use: getSubcellDimension()
    int subcell_dim;

    // DEPRECATED - use: getLocalCellIDs()
    PHX::View<const int*> cell_local_ids_k;

    // DEPRECATED - use: getLocalCellIDs()
    std::vector<size_t> cell_local_ids;

    /// DEPRECATED - use: getCellNodes()
    CellCoordArray cell_node_coordinates;

    /// DEPRECATED - use: getElementBlock()
    std::string block_id;

    /// DEPRECATED - use: getSubcellIndex()
    int subcell_index; //! If workset corresponds to a sub cell, what is the index?

    //! Value correspondes to integration order.  Use the offest for indexing.
    //TEUCHOS_DEPRECATED
    Teuchos::RCP< std::vector<int> > ir_degrees;

    //TEUCHOS_DEPRECATED
    mutable std::vector<Teuchos::RCP<panzer::IntegrationValues2<double> > > int_rules;

    //! Value corresponds to basis type.  Use the offest for indexing.
    //TEUCHOS_DEPRECATED
    Teuchos::RCP< std::vector<std::string> > basis_names;

    //! Static basis function data, key is basis name, value is index in the static_bases vector
    //TEUCHOS_DEPRECATED
    mutable std::vector<Teuchos::RCP< panzer::BasisValues2<double> > > bases;

    /** Grab the face connectivity for this workset
     * DEPRECATED - use: workset.getSubcellConnectivity(workset.numDimensions()-1)
     */
    const panzer::SubcellConnectivity & getFaceConnectivity() const;

    /// Grab the integration rule for a given integration description (throws error if integration doesn't exist)
    const panzer::IntegrationRule & getIntegrationRule(const panzer::IntegrationDescriptor & description) const;

    /// Grab the pure basis (contains data layouts) for a given basis description (throws error if integration doesn't exist)
    const panzer::PureBasis & getBasis(const panzer::BasisDescriptor & description) const;

    /// Get the element block id
    const std::string &
    getElementBlock() const
    {return block_id;}

    /// Get the sideset id (returns "" if not a sideset)
    const std::string &
    getSideset() const
    {return sideset_;}

    /// Get the cell dimension for the mesh
    unsigned int
    numDimensions() const
    {return num_dimensions_;}

    /// Get the subcell index (returns -1 if not a subcell)
    int
    getSubcellIndex() const
    {return subcell_index;}

    /// Get the subcell dimension
    int
    getSubcellDimension() const
    {return subcell_dim;}

    /// Get the node coordinates for the cells
    CellCoordArray
    getCellNodes() const
    {return cell_node_coordinates;}

    /// Get the local cell IDs for the workset
    Kokkos::View<const int*,PHX::Device>
    getLocalCellIDs() const
    {return cell_local_ids_k;}

    /**
     * \brief Get the subcell connectivity for the workset topology
     *
     * Subcells are node/edge/face in association with dimension of mesh.
     *
     * For example:
     *     | Subcell Name | Subcell Dimension
     * ---------------------------------------
     * 1D: | Node         | 0
     *     | Edge         | 1
     * ---------------------------------------
     * 2D: | Node         | 0
     *     | Edge         | 1
     *     | Face         | 2
     * ---------------------------------------
     * 3D: | Node         | 0
     *     | Edge         | 1
     *     | Face         | 2
     *     | Cell         | 3
     * ---------------------------------------
     *
     * \throws If setup has not been called
     * \throws If requested subcell dimension larger than 'numDimensions()'
     *
     * \param[in] subcell_dimension Requested subcell_dimension
     *
     * \return SubcellConnectivity object containing connectivity information for subcell dimension
     *
     */
    const SubcellConnectivity &
    getSubcellConnectivity(const unsigned int subcell_dimension) const;

    /// Check if subcell connectivity exists for a given dimension
    bool
    hasSubcellConnectivity(const unsigned int subcell_dimension) const;

    /**
     * \brief Get the integration values for a given integration description
     *
     * \throws If setup has not been called
     *
     * \param[in] description Descriptor for integration scheme
     * \param[in] lazy_version Get an empty IntegrationValues2 object that will construct/allocate itself on demand (less memory - EXPERIMENTAL)
     *
     * \return Object containing integration values
     */
    const panzer::IntegrationValues2<double> &
    getIntegrationValues(const panzer::IntegrationDescriptor & description,
                         const bool lazy_version=false) const;

    /*
     * \brief Grab the basis values for a given basis description
     *
     * \note An arbitrary integration order is used - only use for getting basis coordinates
     *
     * \throws If setup has not been called
     *
     * \param[in] basis_description Description of requested basis
     * \param[in] lazy_version Get an empty BasisValues2 object that will construct/allocate itself on demand (less memory - EXPERIMENTAL)
     *
     * \return Object containing basis values
     */
    const panzer::BasisValues2<double> &
    getBasisValues(const panzer::BasisDescriptor & basis_description,
                   const bool lazy_version=false) const;

    /*
     * \brief Grab the basis values for a given basis description
     *
     * \throws If setup has not been called
     *
     * \todo This needs to be const, but one path for workset construction requires non-const
     *
     * \param[in] basis_description Description of requested basis
     * \param[in] integration_description Descriptor for integration scheme
     * \param[in] lazy_version Get an empty BasisValues2 object that will construct/allocate itself on demand (less memory - EXPERIMENTAL)
     *
     * \return Object containing basis values
     */
    panzer::BasisValues2<double> &
    getBasisValues(const panzer::BasisDescriptor & basis_description,
                   const panzer::IntegrationDescriptor & integration_description,
                   const bool lazy_version=false) const;

    /*
     * \brief Grab the basis values for a given basis description
     *
     * \throws If setup has not been called
     * \throws if point_descriptor has not been registered
     *
     * \param[in] basis_description Description of requested basis
     * \param[in] point_description Descriptor for points
     * \param[in] lazy_version Get an empty BasisValues2 object that will construct/allocate itself on demand (less memory - EXPERIMENTAL)
     *
     * \return Object containing basis values
     */
    const panzer::BasisValues2<double> &
    getBasisValues(const panzer::BasisDescriptor & basis_description,
                   const panzer::PointDescriptor & point_description,
                   const bool lazy_version=false) const;

    /**
     * \brief Grab the basis values for a given basis description and integration description (throws error if it doesn't exist)
     *
     * \todo This needs to be const, but one path for workset construction requires non-const
     */
    panzer::PointValues2<double> &
    getPointValues(const panzer::PointDescriptor & point_description) const;

    /// Number of total cells in workset (owned, ghost, and virtual)
    int numCells() const {return num_cells;}

    /// Number of cells owned by this workset
    int numOwnedCells() const {return num_owned_cells_;}

    /// Number of cells owned by a different workset
    int numGhostCells() const {return num_ghost_cells_;}

    /// Number of cells not owned by any workset - these are used for boundary conditions
    int numVirtualCells() const {return num_virtual_cells_;}

    /// Provides access to set numbers of cells (required for backwards compatibility)
    void setNumberOfCells(const int owned_cells,
                          const int ghost_cells,
                          const int virtual_cells);

  protected:

    bool setup_;

    int num_owned_cells_;
    int num_ghost_cells_;
    int num_virtual_cells_;

    int num_dimensions_;

    std::string sideset_;

    WorksetOptions options_;

    Teuchos::RCP<const shards::CellTopology> cell_topology_;

    // TODO: Const vs non-const is complicated here due to how point values are generated and orientations are applied
    // Unifying the construction method for worksets will help reduce the clutter here, but point values will almost always be non-const
    mutable std::map<size_t,Teuchos::RCP<const panzer::PureBasis > > _pure_basis_map;
    mutable std::map<size_t,Teuchos::RCP<const panzer::IntegrationRule > > _integration_rule_map;
    mutable std::map<size_t,Teuchos::RCP<const panzer::PointRule > > _point_rule_map;

    mutable std::map<size_t,Teuchos::RCP<const panzer::IntegrationValues2<double> > > integration_values_map_;
    mutable std::map<size_t,Teuchos::RCP<panzer::PointValues2<double> > > point_values_map_;
    mutable std::map<size_t,std::map<size_t,Teuchos::RCP<panzer::BasisValues2<double> > > > basis_integration_values_map_;
    mutable std::map<size_t,std::map<size_t,Teuchos::RCP<panzer::BasisValues2<double> > > > basis_point_values_map_;

    Teuchos::RCP<panzer::SubcellConnectivity> face_connectivity_;

  };

  /** This is the main workset object. Note that it inherits from WorksetDetails, this
    * is to maintain backwards compatibility in the use of the workset object. The addition
    * of a details vector supports things like DG based assembly.
    */
  class Workset : public WorksetDetails {
  public:
    //! Default constructor, identifier is a useless 0 by default
    Workset() : identifier_(0) {}

    //! Constructor that that requires a unique identifier
    Workset(std::size_t identifier) : identifier_(identifier) {}

    //! Set the unique identifier for this workset, this is not an index!
    void setIdentifier(std::size_t identifier) { identifier_ = identifier; }

    //! Get the unique identifier for this workset, this is not an index!
    std::size_t getIdentifier() const { return identifier_; }

    double alpha;
    double beta;
    double time;
    double step_size;
    double stage_number;
    std::vector<double> gather_seeds; // generic gather seeds
    bool evaluate_transient_terms;

    //! other contains details about the side-sharing elements on the other side
    //! of the interface. If Teuchos::nonnull(other), then Workset contains two
    //! WorksetDetails: itself, and other.
    Teuchos::RCP<WorksetDetails> other;

    //! op(0) return *this; op(1) returns *other.
    WorksetDetails& operator()(const int i) {
      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
      return i == 0 ? static_cast<WorksetDetails&>(*this) : *other;
    }
    //! const accessor.
    const WorksetDetails& operator()(const int i) const {
      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
      return i == 0 ? static_cast<const WorksetDetails&>(*this) : *other;
    }
    //! Convenience wrapper to operator() for pointer access.
    WorksetDetails& details(const int i) { return operator()(i); }
    const WorksetDetails& details(const int i) const { return operator()(i); }
    //! Return the number of WorksetDetails this Workset holds.
    size_t numDetails() const { return Teuchos::nonnull(other) ? 2 : 1; }

  private:
    std::size_t identifier_;
  };

  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w);

  /** This accessor may be used by an evaluator to abstract Details
    * Index. "Details Index" may be in the evaluator's constructor
    * ParameterList. If it is, then this accessor's constructor reads it and
    * checks for error. Regardless of whether it is, operator() then wraps
    * workset's access to WorksetDetails fields to provide the correct ones.
    */
  class WorksetDetailsAccessor {
  public:
    //! Default value is 0, which is backwards compatible.
    WorksetDetailsAccessor() : details_index_(0) {}
    //! An evaluator builder sets the details index.
    void setDetailsIndex(const int di) { details_index_ = di; }
    //! Get the details index. Generally, only a gather evaluator needs to know
    //! its details index.
    int getDetailsIndex() const { return details_index_; }
    //! Workset wrapper to extract the correct details. Example: wda(workset).bases[i].
    WorksetDetails& operator()(Workset& workset) const {
      return workset(details_index_);
    }
    //! const accessor.
    const WorksetDetails& operator()(const Workset& workset) const {
      return workset(details_index_);
    }
  private:
    int details_index_;
  };

} // namespace panzer

#endif
