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

#include "Panzer_Workset.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_WorksetUtilities.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_PointGenerator.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_OrientationsInterface.hpp"

#include "Panzer_SubcellConnectivity.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

namespace {

void
applyBV2Orientations(BasisValues2<double> & basis_values,
                     const Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> & local_cell_ids,
                     const Teuchos::RCP<const OrientationsInterface> & orientations_interface)
{
  // This call exists because there is a middle man we have to go through.
  // orientations_interface is a 'local' object, not a workset object so we need to map local cells to workset cells
  if(orientations_interface.is_null())
    return;
  const int num_cells = local_cell_ids.extent_int(0);
  const auto & local_orientations = *orientations_interface->getOrientations();
  std::vector<Intrepid2::Orientation> workset_orientations(num_cells);
  for(int i=0; i<num_cells; ++i)
    workset_orientations[i] = local_orientations[local_cell_ids[i]];
  basis_values.applyOrientations(workset_orientations,num_cells);
}

}

Workset::
Workset():
  num_cells_(0),
  num_owned_cells_(0),
  num_ghost_cells_(0),
  num_virtual_cells_(0),
  setup_(false),
  num_dimensions_(0),
  subcell_index_(-1),
  subcell_dimension_(-1),
  time_(0.),
  identifier_(0),
  step_size_(0.)
{

}

const Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> &
Workset::
getLocalCellIDs() const
{
  TEUCHOS_ASSERT(setup_);
  return local_cell_ids_;
}

const PHX::MDField<Workset::Scalar,Cell,NODE,Dim> &
Workset::
getCellVertices() const
{
  TEUCHOS_ASSERT(setup_);
  return cell_vertices_;
}

void
Workset::
setup(const panzer::LocalMeshPartition & partition,
      const WorksetOptions & options)
{

  TEUCHOS_ASSERT(!setup_);

  options_ = options;

  // Initialize cells size
  num_owned_cells_ = partition.num_owned_cells;
  num_ghost_cells_ = partition.num_ghost_cells;
  num_virtual_cells_ = partition.num_virtual_cells;
  num_cells_ = num_owned_cells_ + num_ghost_cells_ + num_virtual_cells_;

  cell_topology_ = partition.cell_topology;
  subcell_dimension_ = partition.subcell_dimension;
  subcell_index_ = partition.subcell_index;

  // Mesh dimensions
  TEUCHOS_ASSERT(not cell_topology_.is_null());
  num_dimensions_ = cell_topology_->getDimension();

  // Get block id
  element_block_ = partition.element_block_name;
  sideset_ = partition.sideset_name;

  setup_ = true;

  // If there are no cells, then it would be bad to continue on
  if(num_cells_ == 0)
    return;

  // Allocate and fill the local cell indexes for this workset
  {
    Kokkos::View<panzer::LocalOrdinal*, PHX::Device> cell_ids = Kokkos::View<panzer::LocalOrdinal*, PHX::Device>("cell_ids",num_cells_);
    Kokkos::deep_copy(cell_ids, partition.local_cells);
    local_cell_ids_ = cell_ids;
  }

  // Allocate and fill the cell vertices
  {
    // Double check this
    TEUCHOS_ASSERT(partition.cell_vertices.Rank == 3);

    // Grab the size of the cell vertices array
    const size_t num_partition_cells = partition.cell_vertices.extent(0);
    const size_t num_vertices_per_cell = partition.cell_vertices.extent(1);
    const size_t num_dims_per_vertex = partition.cell_vertices.extent(2);

    // Make sure there isn't some strange problem going on
    TEUCHOS_ASSERT(num_partition_cells == size_t(num_cells_));
    TEUCHOS_ASSERT(num_vertices_per_cell > 0);
    TEUCHOS_ASSERT(num_dims_per_vertex > 0);

    // Allocate the worksets copy of the cell vertices
    MDFieldArrayFactory af("",true);
    cell_vertices_ = af.buildStaticArray<double, Cell, NODE, Dim>("cell vertices", num_partition_cells, num_vertices_per_cell, num_dims_per_vertex);

    // Copy vertices over
    const auto & partition_vertices = partition.cell_vertices;
    for(int i=0;i<num_cells_;++i)
      for(int j=0;j<num_vertices_per_cell;++j)
        for(int k=0;k<num_dims_per_vertex;++k)
          cell_vertices_(i,j,k) = partition_vertices(i,j,k);

  }

  // Add the connectivity
  if(partition.has_connectivity){
    auto face_connectivity = Teuchos::rcp(new FaceConnectivity);
    face_connectivity->setup(partition);
    subcell_connectivity_ = face_connectivity;
  }

}


void
Workset::
setup(const std::vector<panzer::LocalMeshPartition> & partitions,
      const WorksetOptions & options)
{
  TEUCHOS_ASSERT(!setup_);
  TEUCHOS_ASSERT(partitions.size() > 0);
  setup(partitions[0], options);
  for(unsigned int i=1; i<partitions.size(); ++i){
    others_.push_back(Teuchos::rcp(new Workset()));
    // We use this call to ensure align_side_integration_points_ is set to true
    others_.back()->setup(partitions[i], options);
  }
}

bool
Workset::
hasSubcellConnectivity(const unsigned int subcell_dimension) const
{
  TEUCHOS_ASSERT(setup_);
  return (subcell_dimension == (numDimensions() - 1)) and (not subcell_connectivity_.is_null());
}


const SubcellConnectivity &
Workset::
getSubcellConnectivity(const unsigned int subcell_dimension) const
{
  TEUCHOS_ASSERT(setup_);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(not hasSubcellConnectivity(subcell_dimension),
                              "Workset::getSubcellConnectivity : Requested subcell dimension "<<subcell_dimension<<" for a "<<num_dimensions_<<"D workset. This is not supported.");
  return *subcell_connectivity_;
}

const panzer::IntegrationValues<Workset::Scalar> &
Workset::
getIntegrationValues(const panzer::IntegrationDescriptor & description) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if exists
  const auto itr = integration_values_map_.find(description.getKey());
  if(itr != integration_values_map_.end())
    return *(itr->second);

  // Since it doesn't exist, we need to create it
  const unsigned int subcell_dimension = numDimensions()-1;
  int num_faces = -1;
  if(hasSubcellConnectivity(subcell_dimension))
    num_faces = getSubcellConnectivity(subcell_dimension).numSubcells();

  // For now, we need to make sure the descriptor lines up with the workset
  if(options_.side_assembly_){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(description.getSide() != getSubcellIndex(),
                                "Workset::getIntegrationValues : Attempted to build integration values for side '"<<description.getSide()<<"', but workset is constructed for side '"<<getSubcellIndex()<<"'");
  }

  auto ir = Teuchos::rcp(new IntegrationRule(description, cell_topology_, num_cells_, num_faces));

  auto iv = Teuchos::rcp(new IntegrationValues<Scalar>("",true));
  iv->setupArrays(ir);
  iv->evaluateValues(getCellVertices());

  if(description.getType() == IntegrationDescriptor::SURFACE)
    iv->makeOrderingUnique();

  if(options_.side_assembly_ and options_.align_side_points_)
    iv->makeOrderingUnique();

  integration_values_map_[description.getKey()] = iv;

  return *iv;

}

const panzer::BasisValues<Workset::Scalar> &
Workset::
getBasisValues(const panzer::BasisDescriptor & description) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if exists
  const auto itr = basis_values_map_.find(description.getKey());
  if(itr != basis_values_map_.end())
    return *(itr->second);

  // TODO: We currently overlap BasisIntegrationValues and BasisValues
  // To counter this we create a fake integration rule at this point to ensure the basis values exist

  IntegrationDescriptor id(2*description.getOrder(), IntegrationDescriptor::VOLUME);

  // We have to have the right integrator if this is a side workset
  if(options_.side_assembly_){
    TEUCHOS_ASSERT(getSubcellIndex() >= 0);
    id = IntegrationDescriptor(2*description.getOrder(), IntegrationDescriptor::SIDE, getSubcellIndex());
  }

  // Tragically we have to create this integration rule, this will be fixed when we no longer need
  // a point rule to generate the basis values
  const auto & iv = getIntegrationValues(id);

  auto bir = Teuchos::rcp(new BasisIRLayout(description.getType(), description.getOrder(), *iv.int_rule));

  auto bv = Teuchos::rcp(new BasisValues2<Scalar>("", true, false));
  bv->setupArrays(bir);
  bv->evaluateValues(iv.ip_coordinates,
                     iv.jac,
                     iv.jac_det,
                     iv.jac_inv,
                     iv.weighted_measure,
                     getCellVertices(),
                     true,
                     num_cells_);

  applyBV2Orientations(*bv,getLocalCellIDs(),options_.orientations_);

  basis_values_map_[description.getKey()] = bv;

  return *bv;

}

panzer::PointValues<Workset::Scalar> &
Workset::
getPointValues(const panzer::PointDescriptor & description) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if exists
  const auto itr = point_values_map_.find(description.getKey());
  if(itr != point_values_map_.end())
    return *(itr->second);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(not description.hasGenerator(),
                              "Point Descriptor of type '"<<description.getType()<<"' does not have associated generator.");

  auto pr = Teuchos::rcp(new PointRule(description, cell_topology_, num_cells_));

  auto pv = Teuchos::rcp(new PointValues2<Scalar>("",true));

  pv->setupArrays(pr);

  // Point values are not necessarily set at the workset level, but can be set by evaluators
  if(description.hasGenerator())
    if(description.getGenerator().hasPoints(*cell_topology_))
      pv->evaluateValues(cell_vertices_, description.getGenerator().getPoints(*cell_topology_),false, num_cells_);

  point_values_map_[description.getKey()] = pv;

  return *pv;

}

const panzer::BasisIntegrationValues<Workset::Scalar> &
Workset::
getBasisIntegrationValues(const panzer::BasisDescriptor & basis_description,
                          const panzer::IntegrationDescriptor & integration_description) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if exists
  const auto itr = basis_integration_values_map_.find(basis_description.getKey());
  if(itr != basis_integration_values_map_.end()){
    const auto & submap = itr->second;
    const auto itr2 = submap.find(integration_description.getKey());
    if(itr2 != submap.end())
      return *(itr2->second);

  }

  // Get the integration values for this description
  const auto & iv = getIntegrationValues(integration_description);

  auto bir = Teuchos::rcp(new BasisIRLayout(basis_description.getType(), basis_description.getOrder(), *iv.int_rule));

  auto biv = Teuchos::rcp(new BasisValues2<Scalar>("", true, true));
  biv->setupArrays(bir);
  if((integration_description.getType() == IntegrationDescriptor::CV_BOUNDARY) or
     (integration_description.getType() == IntegrationDescriptor::CV_SIDE) or
     (integration_description.getType() == IntegrationDescriptor::CV_VOLUME)){

    biv->evaluateValuesCV(iv.ref_ip_coordinates,
                          iv.jac,
                          iv.jac_det,
                          iv.jac_inv,
                          getCellVertices(),
                          true,
                          num_cells_);
  } else {

    if(integration_description.getType() == IntegrationDescriptor::VOLUME){

      // TODO: Eventually we will use the other call, however, that will be part of the BasisValues2 refactor
      // The reason we don't do it now is that there are small differences (machine precision) that break EMPIRE testing
      biv->evaluateValues(iv.cub_points,
                          iv.jac,
                          iv.jac_det,
                          iv.jac_inv,
                          iv.weighted_measure,
                          getCellVertices(),
                          true,
                          num_cells_);

    } else {

      biv->evaluateValues(iv.ref_ip_coordinates,
                          iv.jac,
                          iv.jac_det,
                          iv.jac_inv,
                          iv.weighted_measure,
                          getCellVertices(),
                          true,
                          num_cells_);
    }
  }

  applyBV2Orientations(*biv,getLocalCellIDs(),options_.orientations_);

  basis_integration_values_map_[basis_description.getKey()][integration_description.getKey()] = biv;

  return *biv;

}

const panzer::BasisPointValues<Workset::Scalar> &
Workset::
getBasisPointValues(const panzer::BasisDescriptor & basis_description,
                    const panzer::PointDescriptor & point_description) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if exists
  const auto itr = basis_point_values_map_.find(basis_description.getKey());
  if(itr != basis_point_values_map_.end()){
    const auto & submap = itr->second;
    const auto itr2 = submap.find(point_description.getKey());
    if(itr2 != submap.end())
      return *(itr2->second);

  }

  // Get the integration values for this description
  const auto & pv = getPointValues(point_description);

  auto bir = Teuchos::rcp(new BasisIRLayout(basis_description.getType(), basis_description.getOrder(), *pv.point_rule));

  auto bpv = Teuchos::rcp(new BasisValues2<Scalar>("", true, false));
  bpv->setupArrays(bir);
  bpv->evaluateValues(pv.coords_ref,
                      pv.jac,
                      pv.jac_det,
                      pv.jac_inv,
                      num_cells_);

  // TODO: We call this separetely due to how BasisValues2 is structured - needs to be streamlined
  bpv->evaluateBasisCoordinates(getCellVertices(),num_cells_);

  applyBV2Orientations(*bpv, getLocalCellIDs(), options_.orientations_);

  basis_point_values_map_[basis_description.getKey()][point_description.getKey()] = bpv;

  return *bpv;

}

void
Workset::
setDetails(const std::string & name,
           const Teuchos::RCP<WorksetDetails> & details,
           const bool throw_if_exists)
{
  if(throw_if_exists)
    TEUCHOS_TEST_FOR_EXCEPT_MSG(details_map_.find(name) != details_map_.end(), "Workset::setDetails : Details '"<<name<<"' already exists");
  details_map_[name] = details;
}


Workset &
Workset::
operator()(const unsigned int i)
{
  TEUCHOS_ASSERT(setup_);
  if(i == 0)
    return *this;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(i >= size(), "Workset::operator() : Requested workset "<<i<<". Only "<<size()<<" worksets are accessible.");
  return *others_[i-1];
}

const Workset &
Workset::
operator()(const unsigned int i) const
{
  TEUCHOS_ASSERT(setup_);
  if(i == 0)
    return *this;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(i >= size(), "Workset::operator() (const) : Requested workset "<<i<<". Only "<<size()<<" worksets are accessible.");
  return *others_[i-1];
}

unsigned int
Workset::
size() const
{
  if(not setup_) return 0;
  return 1 + others_.size();
}

std::ostream &
operator<<(std::ostream& os,
           const panzer::Workset& w)
{
  using std::endl;

  os << "Workset\n";
  os << "  Element Block: " << w.getElementBlock() << "\n";
  os << "  Number of Cells: " << w.numCells() << "\n";
  os << "  Local Cell IDs size: " << w.getLocalCellIDs().size() << "\n";
  os << "  Subcell Dimension: " << w.getSubcellDimension() << "\n";
  os << "  Subcell Index: " << w.getSubcellIndex() << "\n";

  return os;
}



































//  // Create the pure basis
//  for(const auto & basis_description : basis_descriptors){
//    // Create and store integration rule
//    Teuchos::RCP<panzer::PureBasis> basis = Teuchos::rcp(new panzer::PureBasis(basis_description, cell_topology_, num_cells));
//    basis_map_[basis_description.getKey()] = basis;
//  }
//
//  // Create the integration terms and basis-integrator pairs
//  for(const auto & integration_description : integration_descriptors){
//
//    int num_faces = -1;
//    if(integration_description.getType() == panzer::IntegrationRule::SURFACE){
//      num_faces = getSubcellConnectivity(numDimensions()-1).numSubcells();
//    }
//
//    // Create and store integration rule
//    Teuchos::RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(integration_description, cell_topology_, num_cells, num_faces));
//
//    // Add integration rule to map
//    integration_rule_map_[integration_description.getKey()] = ir;
//
//    // Create and store integration values
//    Teuchos::RCP<panzer::IntegrationValues<Scalar> > iv = Teuchos::rcp(new panzer::IntegrationValues<Scalar>());
//
//    // TODO: setup integration values
//    TEUCHOS_ASSERT(false);
//
//    // Add integration values to map
//    integration_values_map_[integration_description.getKey()] = iv;
//
//    // We need to generate a integration rule - basis pair for each basis
//    for(const auto & basis_description : basis_descriptors){
//
//      // Grab the basis that was pre-calculated
//      const Teuchos::RCP<const panzer::PureBasis> & basis = basis_map_[basis_description.getKey()];
//
//      // Create a basis ir layout for this pair of integrator and basis
//      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*ir));
//
//        // Create and store basis values
//      {
//        Teuchos::RCP<panzer::BasisValues<Scalar> > bv = Teuchos::rcp(new panzer::BasisValues<double>());
//
//        // TODO: Setub basis values
//        TEUCHOS_ASSERT(false);
//
//        // Add basis values to map
//        basis_values_map_[basis_description.getKey()] = bv;
//      }
//
//      // Setup basis integration values
//      {
//        Teuchos::RCP<panzer::BasisIntegrationValues<Scalar> > biv = Teuchos::rcp(new panzer::BasisIntegrationValues<double>());
//
//        // TODO: Setub basis integration values
//        TEUCHOS_ASSERT(false);
//
//        // Add basis integration values to map
//        basis_integration_values_map_[basis_description.getKey()][integration_description.getKey()] = biv;
//      }
//
//    }
//
//  }
//
//  // Create the point terms and basis-integrator pairs
//  for(const panzer::PointDescriptor & point_description : point_descriptors){
//
//    // get the generaotr, and build some points from a topology
//    auto points = point_description.getGenerator().getPoints(*cell_topology_);
//
//    // Create and store integration rule
//    Teuchos::RCP<panzer::PointRule> pr = Teuchos::rcp(new panzer::PointRule(point_description,cell_topology_,num_cells));
//    point_rule_map_[point_description.getKey()] = pr;
//
//    // Setup point values
//    Teuchos::RCP<panzer::PointValues<double> > pv = Teuchos::rcp(new panzer::PointValues<Scalar>());
//
//    // TODO: Setup point values
//    TEUCHOS_ASSERT(false);
//
//    // Add point values to map
//    point_values_map_[point_description.getKey()] = pv;
//
//    // We need to generate a integration rule - basis pair for each basis
//    for(const auto & basis_description : basis_descriptors){
//
//      // Grab the basis that was pre-calculated
//      const Teuchos::RCP<const panzer::PureBasis> & basis = basis_map_[basis_description.getKey()];
//
//      // Create a basis ir layout for this pair of integrator and basis
//      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*pr));
//
//      // Create and store basis values at points
//      Teuchos::RCP<panzer::BasisPointValues<Scalar> > bv = Teuchos::rcp(new panzer::BasisPointValues<Scalar>());
//
//      // TODO: Setup basis values
//      TEUCHOS_ASSERT(false);
//
//      // Add basis values to map
//      basis_point_values_map_[basis_description.getKey()][point_description.getKey()] = bv;
//    }
//
//  }



//void
//WorksetDetails::setup(const panzer::LocalMeshPartition<int,panzer::Ordinal64> & partition,
//                      const panzer::WorksetNeeds & needs)
//{
//
//
//  const size_t num_cells = partition.local_cells.extent(0);
//
//  num_owned_cells_ = partition.num_owned_cells;
//  num_ghost_cells_ = partition.num_ghost_cells;
//  num_virtual_cells_ = partition.num_virtual_cells;
//
//  subcell_index = -1;
//  block_id = partition.element_block_name;
//
//  Kokkos::View<int*, PHX::Device> cell_ids = Kokkos::View<int*, PHX::Device>("cell_ids",num_cells);
//  Kokkos::deep_copy(cell_ids, partition.local_cells);
//  cell_local_ids_k = cell_ids;
//
//  cell_local_ids.resize(num_cells,-1);
//  for(size_t cell=0;cell<num_cells;++cell){
//    const int local_cell = partition.local_cells(cell);
//    cell_local_ids[cell] = local_cell;
//  }
//
//  auto fc = Teuchos::rcp(new panzer::FaceConnectivity());
//  fc->setup(partition);
//  _face_connectivity = fc;
//
//  setupNeeds(partition.cell_topology, partition.cell_vertices, needs);
//}
//
//void WorksetDetails::setupNeeds(Teuchos::RCP<const shards::CellTopology> cell_topology,
//                                const Kokkos::View<double***,PHX::Device> & cell_vertices,
//                                const panzer::WorksetNeeds & needs)
//{
//
//  const size_t num_cells = cell_vertices.extent(0);
//  const size_t num_vertices_per_cell = cell_vertices.extent(1);
//  const size_t num_dims_per_vertex = cell_vertices.extent(2);
//
//  // Set cell vertices
//  {
//
//    MDFieldArrayFactory af("",true);
//
//    cell_vertex_coordinates = af.template buildStaticArray<double, Cell, NODE, Dim>("cell_vertices",num_cells, num_vertices_per_cell, num_dims_per_vertex);
//
//    for(size_t i=0;i<num_cells;++i)
//      for(size_t j=0;j<num_vertices_per_cell;++j)
//        for(size_t k=0;k<num_dims_per_vertex;++k)
//          cell_vertex_coordinates(i,j,k) = cell_vertices(i,j,k);
//
//  }
//
//  // DEPRECATED - makes sure deprecated arrays are filled with something - this will probably segfault or throw an error
//  panzer::populateValueArrays(num_cells, false, needs, *this);
//
//  const std::vector<panzer::BasisDescriptor> & basis_descriptors = needs.getBases();
//  const std::vector<panzer::IntegrationDescriptor> & integration_descriptors = needs.getIntegrators();
//  const std::vector<panzer::PointDescriptor> & point_descriptors = needs.getPoints();
//
//  // Create the pure basis
//  for(const panzer::BasisDescriptor & basis_description : basis_descriptors){
//    // Create and store integration rule
//    Teuchos::RCP<panzer::PureBasis> basis = Teuchos::rcp(new panzer::PureBasis(basis_description, cell_topology, num_cells));
//    _purebasis_map_[basis_description.getKey()] = basis;
//  }
//
//  // Create the integration terms and basis-integrator pairs
//  for(const panzer::IntegrationDescriptor & integration_description : integration_descriptors){
//
//    int num_faces = -1;
//    if(integration_description.getType() == panzer::IntegrationRule::SURFACE){
//      num_faces = getFaceConnectivity().numSubcells();
//    }
//
//    // Create and store integration rule
//    Teuchos::RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(integration_description, cell_topology, num_cells, num_faces));
//    integration_rule_map_[integration_description.getKey()] = ir;
//
//    // Create and store integration values
//    Teuchos::RCP<panzer::IntegrationValues2<double> > iv = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
//    iv->setupArrays(ir);
//    iv->evaluateValues(cell_vertex_coordinates,num_cells);
//    _integrator_map[integration_description.getKey()] = iv;
//
//    // We need to generate a integration rule - basis pair for each basis
//    for(const panzer::BasisDescriptor & basis_description : basis_descriptors){
//
//      // Grab the basis that was pre-calculated
//      const Teuchos::RCP<const panzer::PureBasis> & basis = _purebasis_map_[basis_description.getKey()];
//
//      // Create a basis ir layout for this pair of integrator and basis
//      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*ir));
//
//      // Create and store basis values
//      Teuchos::RCP<panzer::BasisValues2<double> > bv = Teuchos::rcp(new panzer::BasisValues2<double>("",true,true));
//      bv->setupArrays(b_layout);
//      if(ir->getType() == panzer::IntegrationDescriptor::SURFACE){
//        bv->evaluateValues(iv->ref_ip_coordinates,
//                           iv->jac,
//                           iv->jac_det,
//                           iv->jac_inv,
//                           iv->weighted_measure,
//                           cell_vertex_coordinates,
//                           true,
//                           num_cells);
//      } else if((ir->getType() == panzer::IntegrationDescriptor::CV_VOLUME)
//          or (ir->getType() == panzer::IntegrationDescriptor::CV_SIDE)
//          or (ir->getType() == panzer::IntegrationDescriptor::CV_BOUNDARY)){
//        bv->evaluateValuesCV(iv->ref_ip_coordinates,
//                             iv->jac,
//                             iv->jac_det,
//                             iv->jac_inv);
//      } else {
//        bv->evaluateValues(iv->cub_points,
//                           iv->jac,
//                           iv->jac_det,
//                           iv->jac_inv,
//                           iv->weighted_measure,
//                           cell_vertex_coordinates,
//                           true,
//                           num_cells);
//      }
//      basis_map_[basis_description.getKey()][integration_description.getKey()] = bv;
//    }
//
//  }
//
//  // Create the point terms and basis-integrator pairs
//  for(const panzer::PointDescriptor & point_description : point_descriptors){
//
//    // get the generaotr, and build some points from a topology
//    auto points = point_description.getGenerator().getPoints(*cell_topology);
//
//    // Create and store integration rule
//    Teuchos::RCP<panzer::PointRule> pr = Teuchos::rcp(new panzer::PointRule(point_description,
//                                                                            cell_topology,
//                                                                            num_cells));
//    point_rule_map_[point_description.getKey()] = pr;
//
//    // Create and store integration values
//    Teuchos::RCP<panzer::PointValues2<double> > pv = Teuchos::rcp(new panzer::PointValues2<double>("",true));
//    pv->setupArrays(pr);
//    pv->evaluateValues(cell_vertex_coordinates,points);
//
//    _point_map[point_description.getKey()] = pv;
//
//    // We need to generate a integration rule - basis pair for each basis
//    for(const panzer::BasisDescriptor & basis_description : basis_descriptors){
//
//      // Grab the basis that was pre-calculated
//      const Teuchos::RCP<const panzer::PureBasis> & basis = _purebasis_map_[basis_description.getKey()];
//
//      // Create a basis ir layout for this pair of integrator and basis
//      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*pr));
//
//      // Create and store basis values
//      Teuchos::RCP<panzer::BasisValues2<double> > bv = Teuchos::rcp(new panzer::BasisValues2<double>("",true,true));
//      bv->setupArrays(b_layout);
//
//      bv->evaluateValues(pv->coords_ref,
//                         pv->jac,
//                         pv->jac_det,
//                         pv->jac_inv);
//
//      basis_map_[basis_description.getKey()][point_description.getKey()] = bv;
//    }
//
//  }
//
//}
//
//const panzer::SubcellConnectivity &
//WorksetDetails::getFaceConnectivity() const
//{
//  TEUCHOS_ASSERT(_face_connectivity != Teuchos::null);
//  return *_face_connectivity;
//}
//
//const panzer::IntegrationValues2<double> &
//WorksetDetails::getIntegrationValues(const panzer::IntegrationDescriptor & description) const
//{
//  const auto itr = _integrator_map.find(description.getKey());
//  TEUCHOS_ASSERT(itr != _integrator_map.end());
//  return *(itr->second);
//}
//
//const panzer::IntegrationRule &
//WorksetDetails::getIntegrationRule(const panzer::IntegrationDescriptor & description) const
//{
//  const auto itr = integration_rule_map_.find(description.getKey());
//  TEUCHOS_ASSERT(itr != integration_rule_map_.end());
//  return *(itr->second);
//}
//
//panzer::BasisValues2<double> &
//WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description,
//                               const panzer::IntegrationDescriptor & integration_description)
//{
//  const auto itr = basis_map_.find(basis_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == basis_map_.end(),
//                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" "
//                              "of order " + std::to_string(basis_description.getOrder()));
//  const auto & integration_map = itr->second;
//  const auto itr2 = integration_map.find(integration_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == integration_map.end(),
//                              "Workset::getBasisValues: Can't find integration " + std::to_string(integration_description.getType()) + " "
//                              "of order " + std::to_string(integration_description.getOrder()));
//  return *(itr2->second);
//}
//
//const panzer::BasisValues2<double> &
//WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description,
//                               const panzer::PointDescriptor & point_description) const
//{
//  const auto itr = basis_map_.find(basis_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == basis_map_.end(),
//                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" "
//                              "of order " + std::to_string(basis_description.getOrder()));
//  const auto & point_map = itr->second;
//  const auto itr2 = point_map.find(point_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == point_map.end(),
//                              "Workset::getBasisValues: Can't find point values \"" + point_description.getType() + "\"");
//  return *(itr2->second);
//}
//
//const panzer::BasisValues2<double> &
//WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description,
//                               const panzer::IntegrationDescriptor & integration_description) const
//{
//  const auto itr = basis_map_.find(basis_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == basis_map_.end(),
//                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" "
//                              "of order " + std::to_string(basis_description.getOrder()));
//  const auto & integration_map = itr->second;
//  const auto itr2 = integration_map.find(integration_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == integration_map.end(),
//                              "Workset::getBasisValues: Can't find integration " + std::to_string(integration_description.getType()) + " "
//                              "of order " + std::to_string(integration_description.getOrder()));
//  return *(itr2->second);
//}
//
//const panzer::PointValues2<double> &
//WorksetDetails::getPointValues(const panzer::PointDescriptor & point_description) const
//{
//  const auto itr = _point_map.find(point_description.getKey());
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == _point_map.end(),
//                              "Workset::getPointValues: Can't find point values \"" + point_description.getType() + "\"");
//  return *(itr->second);
//}
//
//const panzer::PureBasis &
//WorksetDetails::getBasis(const panzer::BasisDescriptor & description) const
//{
//  const auto itr = _purebasis_map_.find(description.getKey());
//  TEUCHOS_ASSERT(itr != _purebasis_map_.end());
//  return *(itr->second);
//}


}
