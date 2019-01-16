// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER


#ifndef PANZER_WORKSET_HPP
#define PANZER_WORKSET_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <iostream>

#include "Panzer_Dimension.hpp"
#include "Panzer_BasisValues.hpp"
#include "Panzer_BasisIntegrationValues.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_Dimension.hpp"


#include "Panzer_SubcellConnectivity.hpp"

#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"



namespace panzer {

/// Used for construction of workset
struct WorksetNeeds;

/// Used for construction of workset
struct LocalMeshPartition;

class OrientationsInterface;

/// Used for defining subcell connectivity
class SubcellConnectivity;

template<typename T>
using BasisPointValues = BasisValues2<T>;

template<typename T>
using IntegrationValues = IntegrationValues2<T>;

template<typename T>
using PointValues = PointValues2<T>;

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

struct
WorksetDetails
{
  WorksetDetails() = default;
  virtual ~WorksetDetails() = default;
};

struct
WorksetFADDetails:
    public WorksetDetails
{

  WorksetFADDetails():
    alpha(0.), beta(0.)
  { }

  double alpha;
  double beta;
  std::vector<double> gather_seeds;
};

struct
WorksetStepperDetails:
    public WorksetDetails
{

  WorksetStepperDetails():
    evaluate_transient_terms(false),
    stage_number(0)
  { }

  bool evaluate_transient_terms;

  // TODO: why is this a double?
  double stage_number;
};

class
Workset
{
public:

  /// Default floating point type
  using Scalar = double;

  /// Default constructor
  Workset();

  /**
   * \brief Setup the workset for a given partition of the mesh
   *
   * \note If called again with a new partition, all internal information is cleared
   *
   * \param[in] partition Partition representing workset's chunk of mesh
   * \param[in] options Set of options for constructing workset values objects
   */
  void
  setup(const panzer::LocalMeshPartition & partition,
        const WorksetOptions & options);

  /**
   * \brief Setup the workset for a given interface on the mesh
   *
   * The other 'zones' for the workset can be accessed with 'operator()'
   * The number of zones is given by 'size'
   *
   * \param[in] partitions Partitions representing workset's zones of mesh
   * \param[in] options Set of options for constructing workset values objects
   */
  void
  setup(const std::vector<panzer::LocalMeshPartition> & partitions,
        const WorksetOptions & options);

  /**
   * \brief Array of cell indexes for indexing into mpi rank data level
   *
   * \throws If setup has not been called
   *
   * \return View of local cell indexes
   */
  const Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> &
  getLocalCellIDs() const;

  /**
   * \brief Access cell vertices for workset cells
   *
   * \throws If setup has not been called
   *
   * \note Intrepid requires this to be a DynRankView - there may be a way to wrap a View in a DynRankView
   *
   * \return View of local cell vertices
   */
//  Kokkos::View<const Scalar***, PHX::Device>
  const PHX::MDField<Scalar,Cell,NODE,Dim> &
  getCellVertices() const;

  /**
   * \brief Get the block identification string for workset
   *
   * \return String identifier for block
   */
  const std::string &
  getElementBlock() const
  {return element_block_;}

  /**
   * \brief Get the sideset identification string for workset
   *
   * \return String identifier for sideset - returns empty string if no sideset is assigned
   */
  const std::string &
  getSideset() const
  {return sideset_;}

  /**
   * \brief Get the subcell index for this workset
   *
   * \note If workset is not associated with a subcell, function will return -1
   *
   * \return Subcell index associated with this workset. Returns -1 for no subcell association.
   */
  int
  getSubcellIndex() const
  {return subcell_index_;}

  /**
   * \brief Get the subcell dimension for this workset
   *
   * \note If workset is not associated with a subcell, function will return -1
   *
   * \return Subcell index associated with this workset. Returns -1 for no subcell association.
   */
  int
  getSubcellDimension() const
  {return subcell_dimension_;}

  /// Number of dimensions in mesh
  unsigned int
  numDimensions() const
  {return num_dimensions_;}

  /**
   * \brief Check if subcell connectivity is available for a given dimension
   *
   * \return True if subcell connectivity exists
   */
  bool
  hasSubcellConnectivity(const unsigned int subcell_dimension) const;

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
   * \throws If requested subcell dimension is negative or larger than 'numDimensions()'
   *
   * \param[in] subcell_dimension Requested subcell_dimension
   *
   * \return SubcellConnectivity object containing connectivity information for subcell dimension
   *
   */
  const SubcellConnectivity &
  getSubcellConnectivity(const unsigned int subcell_dimension) const;

  /**
   * \brief Get the integration values for a given integration description
   *
   * \throws If setup has not been called
   *
   * \param[in] description Descriptor for integration scheme
   *
   * \return Object containing integration values
   */
  const panzer::IntegrationValues<Scalar> &
  getIntegrationValues(const panzer::IntegrationDescriptor & description) const;

  /*
   * \brief Grab the basis values for a given basis description
   *
   * \throws If setup has not been called
   *
   * \param[in] basis_description Description of requested basis
   *
   * \return Object containing basis values
   */
  const panzer::BasisValues<Scalar> &
  getBasisValues(const panzer::BasisDescriptor & basis_description) const;

  /**
   * \brief Get the point values for a given point description
   *
   * \throws If setup has not been called
   *
   * \note This returns a non-const object because we may need to directly set its values in an evaluator
   *
   * \param[in] point_description Description of requested points
   *
   * \return Object containing basis values evaluated using given integration scheme
   */
  panzer::PointValues<Scalar> &
  getPointValues(const panzer::PointDescriptor & point_description) const;

  /*
   * \brief Grab the basis integration values for a given basis description and integration description
   *
   * \throws If setup has not been called
   *
   * \param[in] basis_description Description of requested basis
   * \param[in] integration_description Description of requested integration scheme
   *
   * \return Object containing basis values evaluated using given integration scheme
   */
  const panzer::BasisIntegrationValues<Scalar> &
  getBasisIntegrationValues(const panzer::BasisDescriptor & basis_description,
                            const panzer::IntegrationDescriptor & integration_description) const;

  /*
   * \brief Grab the basis points values for a given basis description and points description
   *
   * \throws If setup has not been called
   *
   * \param[in] basis_description Description of requested basis
   * \param[in] points_description Description of requested points
   *
   * \return Object containing basis values evaluated using given integration scheme
   */
  const panzer::BasisPointValues<Scalar> &
  getBasisPointValues(const panzer::BasisDescriptor & basis_description,
                      const panzer::PointDescriptor & point_description) const;

  /// Get the total number of cells in workset
  int
  numCells() const
  {return num_cells_;}

  /// Number of cells owned by this workset
  int
  numOwnedCells() const
  {return num_owned_cells_;}

  /// Number of cells owned by other worksets
  int
  numGhostCells() const
  {return num_ghost_cells_;}

  /// Number of cells not owned by any workset - these are used for boundary conditions
  int
  numVirtualCells() const
  {return num_virtual_cells_;}

  /**
   * \brief Set the unique identifier for this workset
   *
   * \param[in] identifier Identifier for workset
   */
  void
  setIdentifier(const size_t identifier)
  {identifier_ = identifier;}

  /// Get the unique identifier for this workset
  size_t
  getIdentifier() const
  {return identifier_;}

  /**
   * \brief Check if a details object exists with a given name
   *
   * \param[in] name Name of details
   *
   * \return True if details has been registered
   */
  bool
  hasDetails(const std::string & name) const
  {return details_map_.find(name) != details_map_.end();}

  /**
   * \brief Get a details object stored in the workset
   *
   * \throws If details not found
   * \throws If type of details is incorrect
   *
   * \param[in] name Name of details
   *
   * \return Reference to details object
   */
  template<typename T>
  const T &
  getDetails(const std::string & name) const
  {
    const auto itr = details_map_.find(name);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == details_map_.end(), "Workset::getDetails : Details with name '"<<name<<"' does not exist");
    auto cast_ptr = Teuchos::rcp_dynamic_cast<const T>(itr->second, true);
    return *cast_ptr;
  }

  /**
   * \brief Get a details object stored in the workset
   *
   * \throws If details not found
   * \throws If type of details is incorrect
   *
   * \param[in] name Name of details
   *
   * \return Reference to details object
   */
  template<typename T>
  T &
  getDetails(const std::string & name)
  {
    auto itr = details_map_.find(name);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == details_map_.end(), "Workset::getDetails : Details with name '"<<name<<"' does not exist");
    auto cast_ptr = Teuchos::rcp_dynamic_cast<T>(itr->second, true);
    return *cast_ptr;
  }

  /**
   * \brief Set a details object in the workset
   *
   * \throws If details already exists and throw_if_exists is true
   *
   * \param[in] name Name of details
   * \param[in] details Details to set
   * \param[in] throw_if_exists If name is already associated with details, throw error
   */
  void
  setDetails(const std::string & name,
             const Teuchos::RCP<WorksetDetails> & details,
             const bool throw_if_exists=false);

  /// Get time
  double
  getTime() const
  {return time_;}

  /// Set time
  void
  setTime(const double time)
  {time_ = time;}

  /// Get dt
  double
  getTimeStepSize() const
  {return step_size_;}

  /// Set dt
  void
  setTimeStepSize(const double dt)
  {step_size_ = dt;}

  // **********************************************************************
  // MOST OF THE FOLLOWING HAS BEEN REPLACED OR MOVED
//
//  // Replaced with numCells()
//  index_t num_cells;
//
//  // Replaced with subcellDimension()
//  int subcell_dim;
//
//  // Replaced with getSubcellIndex()
//  int subcell_index;
//
//  // Replaced with getIntegrationRule(descriptor)
//  Teuchos::RCP< std::vector<int> > ir_degrees;
//
//  // Replaced with getIntegrationValues(descriptor)
//  std::vector<Teuchos::RCP<panzer::IntegrationValues2<double> > > int_rules;
//
//  // Replaced with getBasisValues(descriptor) and getBasisIntegrationValues(descriptor, descriptor)
//  std::vector<Teuchos::RCP< panzer::BasisValues2<double> > > bases;
//
//  // Replaced with getLocalCellIndexes()
//  Kokkos::View<const int*,PHX::Device> cell_local_ids_k;
//  std::vector<GO> cell_local_ids;
//
//  // Replaced with getCellVertices()
//  CellCoordArray cell_vertex_coordinates;
//
//  // Replaced with getElementBlock()
//  std::string block_id;
//
//  /// Replaced with getSubcellConnectivity(dimension)
//  const panzer::SubcellConnectivity & getFaceConnectivity() const;
//
//
//  // Replaced with getTime()
//  double time;
//
//  // Replaced with getTimeStepSize() -> getDT()?
//  double step_size;
//
//  double alpha;                       // Moved to WorksetFADDetails
//  double beta;                        // Moved to WorksetFADDetails
//  std::vector<double> gather_seeds;   // Moved to WorksetFADDetails
//
//  double stage_number;                // Moved to WorksetStepperDetails
//  bool evaluate_transient_terms;      // Moved to WorksetStepperDetails
//
//  // Remove completely? - not sure what this is used for
//  Teuchos::RCP< std::vector<std::string> > basis_names;

  /**
   * \brief Get Workset associated with index i
   *
   * \param[in] i Index of workset
   *
   * \return Workset associated with index i
   */
  Workset &
  operator()(const unsigned int i);

  /**
   * \brief Get Workset associated with index i
   *
   * \param[in] i Index of workset
   *
   * \return Workset associated with index i
   */
  const Workset &
  operator()(const unsigned int i) const;

  /**
   * \brief Get the number of worksets associated with this workset
   *
   * \return Number of worksets
   */
  unsigned int
  size() const;

  // **********************************************************************

protected:

  /// List of other worksets associated with this workset - e.g. boundary worksets
  std::vector<Teuchos::RCP<Workset>> others_;

  /// Current time
  double time_;

  /// I assume this is dt
  double step_size_;

  /// Unique identifier for this workset
  size_t identifier_;

  /// Number of total cells in workset
  int num_cells_;

  /// Number of cells owned by this workset
  int num_owned_cells_;

  /// Number of cells owned by other worksets
  int num_ghost_cells_;

  /// Number of cells not owned by worksets (normally used to represent boundary conditions)
  int num_virtual_cells_;

  /// Has workset been setup?
  bool setup_;

  /// Block ID string
  std::string element_block_;

  /// Sideset ID string
  std::string sideset_;

  // Mesh dimension
  unsigned int num_dimensions_;

  /// Subcell Dimension - only used if workset is associated with a single subcell
  int subcell_dimension_;

  /// Subcell Index - only used if workset is associated with a single subcell
  int subcell_index_;

  /// Local cell ids with respect to local (mpi rank) topology
  Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> local_cell_ids_;

  /// Cell geometry for local chunk of mesh
  PHX::MDField<Scalar,Cell,NODE,Dim> cell_vertices_;

  /// Cell topology for local chunk of mesh
  Teuchos::RCP<const shards::CellTopology> cell_topology_;

  /// Map of integration values
  mutable std::map<size_t,Teuchos::RCP<panzer::IntegrationValues<Scalar> > > integration_values_map_;

  /// Map of basis values
  mutable std::map<size_t,Teuchos::RCP<panzer::BasisValues<Scalar> > > basis_values_map_;

  /// Map of basis integration values
  mutable std::map<size_t,std::map<size_t,Teuchos::RCP<panzer::BasisIntegrationValues<Scalar> > > > basis_integration_values_map_;

  /// Map of basis integration values
  mutable std::map<size_t,std::map<size_t,Teuchos::RCP<panzer::BasisPointValues<Scalar> > > > basis_point_values_map_;

  /// Map of point values
  mutable std::map<size_t,Teuchos::RCP<panzer::PointValues<double> > > point_values_map_;

  /// Face connectivity
  // TODO: This needs to be expanded for nodal/edge connectivities
  Teuchos::RCP<panzer::SubcellConnectivity> subcell_connectivity_;

  /// List of workset details for auxiliary storage
  std::map<std::string, Teuchos::RCP<WorksetDetails>> details_map_;

  /// Options for generating values objects
  WorksetOptions options_;

};

/**
 * \brief Stream operator for workset
 *
 * \param[inout] os Stream object
 * \param[in] workset Workset to print to stream
 *
 * \return Reference to stream object
 */
std::ostream &
operator<<(std::ostream & os,
           const panzer::Workset & workset);








//
//  struct WorksetNeeds;
//
//  template<typename LO, typename GO>
//  struct LocalMeshPartition;
//
//  class SubcellConnectivity;
//
//  /** This is used within the workset to make edge based (DG like) assembly
//    * an easier task. This basically allows separation of the workset abstraction
//    * from how it is accessed.
//    */
//  class WorksetDetails {
//  public:
//    typedef PHX::MDField<double,Cell,NODE,Dim> CellCoordArray;
//
//    typedef std::size_t GO;
//    typedef int LO;
//
//    //! Default constructor
//    WorksetDetails()
//      : _num_owned_cells(-1)
//      , _num_ghost_cells(-1)
//      , _num_virtual_cells(-1)
//    { }
//
//    //! Constructs the workset details from a given chunk of the mesh
//    void setup(const panzer::LocalMeshPartition<int,panzer::Ordinal64> & partition, const panzer::WorksetNeeds & needs);
//
//    void setupNeeds(Teuchos::RCP<const shards::CellTopology> cell_topology,
//                    const Kokkos::View<double***,PHX::Device> & cell_vertices,
//                    const panzer::WorksetNeeds & needs);
//
//    Kokkos::View<const int*,PHX::Device> cell_local_ids_k;
//    std::vector<GO> cell_local_ids;
//    CellCoordArray cell_vertex_coordinates;
//    std::string block_id;
//
//    int subcell_index; //! If workset corresponds to a sub cell, what is the index?
//
//    //! Value correspondes to integration order.  Use the offest for indexing.
//    //TEUCHOS_DEPRECATED
//    Teuchos::RCP< std::vector<int> > ir_degrees;
//
//    //TEUCHOS_DEPRECATED
//    std::vector<Teuchos::RCP<panzer::IntegrationValues2<double> > > int_rules;
//
//    //! Value corresponds to basis type.  Use the offest for indexing.
//    //TEUCHOS_DEPRECATED
//    Teuchos::RCP< std::vector<std::string> > basis_names;
//
//    //! Static basis function data, key is basis name, value is index in the static_bases vector
//    //TEUCHOS_DEPRECATED
//    std::vector<Teuchos::RCP< panzer::BasisValues2<double> > > bases;
//
//    /// Grab the face connectivity for this workset
//    const panzer::SubcellConnectivity & getFaceConnectivity() const;
//
//    /// Grab the integration values for a given integration description (throws error if integration doesn't exist)
//    const panzer::IntegrationValues2<double> & getIntegrationValues(const panzer::IntegrationDescriptor & description) const;
//
//    /// Grab the integration rule (contains data layouts) for a given integration description (throws error if integration doesn't exist)
//    const panzer::IntegrationRule & getIntegrationRule(const panzer::IntegrationDescriptor & description) const;
//
//    /// Grab the basis values for a given basis description and integration description (throws error if it doesn't exist)
//    panzer::BasisValues2<double> & getBasisValues(const panzer::BasisDescriptor & basis_description,
//                                                  const panzer::IntegrationDescriptor & integration_description);
//
//    /// Grab the basis values for a given basis description and integration description (throws error if it doesn't exist)
//    const panzer::BasisValues2<double> & getBasisValues(const panzer::BasisDescriptor & basis_description,
//                                                        const panzer::IntegrationDescriptor & integration_description) const;
//
//    /// Grab the basis values for a given basis description and integration description (throws error if it doesn't exist)
//    const panzer::BasisValues2<double> & getBasisValues(const panzer::BasisDescriptor & basis_description,
//                                                        const panzer::PointDescriptor & point_description) const;
//
//    /// Grab the basis values for a given basis description and integration description (throws error if it doesn't exist)
//    const panzer::PointValues2<double> & getPointValues(const panzer::PointDescriptor & point_description) const;
//
//    /// Grab the pure basis (contains data layouts) for a given basis description (throws error if integration doesn't exist)
//    const panzer::PureBasis & getBasis(const panzer::BasisDescriptor & description) const;
//
//    /// Number of cells owned by this workset
//    int numOwnedCells() const {return _num_owned_cells;}
//
//    /// Number of cells owned by a different workset
//    int numGhostCells() const {return _num_ghost_cells;}
//
//    /// Number of cells not owned by any workset - these are used for boundary conditions
//    int numVirtualCells() const {return _num_virtual_cells;}
//
//    /// Provides access to set numbers of cells (required for backwards compatibility)
//    void setNumberOfCells(int o_cells,int g_cells,int v_cells)
//    {
//      _num_owned_cells = o_cells;
//      _num_ghost_cells = g_cells;
//      _num_virtual_cells = v_cells;
//    }
//
//  protected:
//
//    int _num_owned_cells;
//    int _num_ghost_cells;
//    int _num_virtual_cells;
//
//    std::map<size_t,Teuchos::RCP<const panzer::IntegrationRule > > _integration_rule_map;
//    std::map<size_t,Teuchos::RCP<const panzer::IntegrationValues2<double> > > _integrator_map;
//
//    std::map<size_t,Teuchos::RCP<const panzer::PureBasis > > _pure_basis_map;
//    std::map<size_t,std::map<size_t,Teuchos::RCP<panzer::BasisValues2<double> > > > _basis_map;
//
//    std::map<size_t,Teuchos::RCP<const panzer::PointRule > > _point_rule_map;
//    std::map<size_t,Teuchos::RCP<const panzer::PointValues2<double> > > _point_map;
//
//    Teuchos::RCP<panzer::SubcellConnectivity> _face_connectivity;
//
//  };
//
//  /** This is the main workset object. Not that it inherits from WorksetDetails, this
//    * is to maintain backwards compatibility in the use of the workset object. The addition
//    * of a details vector supports things like DG based assembly.
//    */
//  class Workset : public WorksetDetails {
//  public:
//    //! Default constructor, identifier is a useless 0 by default
//    Workset() : identifier_(0) {}
//
//    //! Constructor that that requires a unique identifier
//    Workset(std::size_t identifier) : identifier_(identifier) {}
//
//    //! Set the unique identifier for this workset, this is not an index!
//    void setIdentifier(std::size_t identifier) { identifier_ = identifier; }
//
//    //! Get the unique identifier for this workset, this is not an index!
//    std::size_t getIdentifier() const { return identifier_; }
//
//    index_t num_cells;
//    int subcell_dim; //! If workset corresponds to a sub cell, what is the dimension?
//
//    double alpha;
//    double beta;
//    double time;
//    double step_size;
//    double stage_number;
//    std::vector<double> gather_seeds; // generic gather seeds
//    bool evaluate_transient_terms;
//
//    //! other contains details about the side-sharing elements on the other side
//    //! of the interface. If Teuchos::nonnull(other), then Workset contains two
//    //! WorksetDetails: itself, and other.
//    Teuchos::RCP<WorksetDetails> other;
//
//    //! op(0) return *this; op(1) returns *other.
//    WorksetDetails& operator()(const int i) {
//      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
//      return i == 0 ? static_cast<WorksetDetails&>(*this) : *other;
//    }
//    //! const accessor.
//    const WorksetDetails& operator()(const int i) const {
//      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
//      return i == 0 ? static_cast<const WorksetDetails&>(*this) : *other;
//    }
//    //! Convenience wrapper to operator() for pointer access.
//    WorksetDetails& details(const int i) { return operator()(i); }
//    const WorksetDetails& details(const int i) const { return operator()(i); }
//    //! Return the number of WorksetDetails this Workset holds.
//    size_t numDetails() const { return Teuchos::nonnull(other) ? 2 : 1; }
//
//  private:
//    std::size_t identifier_;
//  };
//
//  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w);
//
//  /** This accessor may be used by an evaluator to abstract Details
//    * Index. "Details Index" may be in the evaluator's constructor
//    * ParameterList. If it is, then this accessor's constructor reads it and
//    * checks for error. Regardless of whether it is, operator() then wraps
//    * workset's access to WorksetDetails fields to provide the correct ones.
//    */
//  class WorksetDetailsAccessor {
//  public:
//    //! Default value is 0, which is backwards compatible.
//    WorksetDetailsAccessor() : details_index_(0) {}
//    //! An evaluator builder sets the details index.
//    void setDetailsIndex(const int di) { details_index_ = di; }
//    //! Get the details index. Generally, only a gather evaluator needs to know
//    //! its details index.
//    int getDetailsIndex() const { return details_index_; }
//    //! Workset wrapper to extract the correct details. Example: wda(workset).bases[i].
//    WorksetDetails& operator()(Workset& workset) const {
//      return workset(details_index_);
//    }
//    //! const accessor.
//    const WorksetDetails& operator()(const Workset& workset) const {
//      return workset(details_index_);
//    }
//  private:
//    int details_index_;
//  };

} // namespace panzer

#endif

