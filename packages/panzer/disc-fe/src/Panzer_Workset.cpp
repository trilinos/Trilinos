// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Workset.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_PointGenerator.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_HashUtils.hpp"
#include "Panzer_ConvertNormalToRotationMatrix.hpp"

#include "Panzer_SubcellConnectivity.hpp"
#include "Panzer_OrientationsInterface.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

namespace {

void buildLocalOrientations(const int num_cells,
                            const Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> & local_cell_ids,
                            const Teuchos::RCP<const OrientationsInterface> & orientations_interface,
                            std::vector<Intrepid2::Orientation> & workset_orientations)
{
  // from a list of cells in the workset, extract the subset of orientations that correspond

  const auto & local_orientations = *orientations_interface->getOrientations();
  workset_orientations.resize(num_cells);

  // We can only apply orientations to owned and ghost cells - virtual cells are ignored (no orientations available)
  auto local_cell_ids_host = Kokkos::create_mirror_view(local_cell_ids);
  Kokkos::deep_copy(local_cell_ids_host, local_cell_ids);
  for(int i=0; i<num_cells; ++i)
    workset_orientations[i] = local_orientations[local_cell_ids_host[i]];
}

void
applyBV2Orientations(const int num_cells,
                     BasisValues2<double> & basis_values,
                     const Kokkos::View<const panzer::LocalOrdinal*,PHX::Device> & local_cell_ids,
                     const Teuchos::RCP<const OrientationsInterface> & orientations_interface)
{
  // This call exists because there is a middle man we have to go through.
  // orientations_interface is a 'local' object, not a workset object so we need to map local cells to workset cells

  // If the object doesn't exist, feel free to skip applying orientations, they aren't needed in some cases (e.g. DG/FV)
  if(orientations_interface.is_null())
    return;

  // Ignore this operation if it has already been applied
  if(basis_values.orientationsApplied())
    return;

  // pull out the subset of orientations required for this workset
  std::vector<Intrepid2::Orientation> workset_orientations(num_cells);
  buildLocalOrientations(num_cells,local_cell_ids,orientations_interface, workset_orientations);

  basis_values.applyOrientations(workset_orientations,num_cells);
}

}

WorksetDetails::
WorksetDetails()
  : num_cells(0)
  , subcell_dim(-1)
  , subcell_index(-1)
  , ir_degrees(new std::vector<int>())
  , basis_names(new std::vector<std::string>())
  , setup_(false)
  , num_owned_cells_(0)
  , num_ghost_cells_(0)
  , num_virtual_cells_(0)
  , num_dimensions_(-1)
{ }

void
WorksetDetails::
setup(const panzer::LocalMeshPartition & partition,
      const WorksetOptions & options)
{

  num_cells = partition.local_cells.extent(0);
  num_owned_cells_ = partition.num_owned_cells;
  num_ghost_cells_ = partition.num_ghstd_cells;
  num_virtual_cells_ = partition.num_virtual_cells;
  options_ = options;

  TEUCHOS_ASSERT(num_cells == num_owned_cells_ + num_ghost_cells_ + num_virtual_cells_);

  TEUCHOS_ASSERT(partition.cell_topology != Teuchos::null);
  cell_topology_ = partition.cell_topology;

  num_dimensions_ = cell_topology_->getDimension();
  subcell_dim = partition.subcell_dimension;
  subcell_index = partition.subcell_index;
  block_id = partition.element_block_name;
  sideset_ = partition.sideset_name;


  // Allocate and fill the local cell indexes for this workset
  {
    Kokkos::View<LocalOrdinal*, PHX::Device> cell_ids = Kokkos::View<LocalOrdinal*, PHX::Device>("cell_ids",num_cells);
    Kokkos::deep_copy(cell_ids, partition.local_cells);
    cell_local_ids_k = cell_ids;

    // DEPRECATED - only retain for backward compatability
    auto local_cells_h = Kokkos::create_mirror_view(partition.local_cells);
    Kokkos::deep_copy(local_cells_h, partition.local_cells);
    cell_local_ids.resize(num_cells,-1);
    for(int cell=0;cell<num_cells;++cell){
      const int local_cell = local_cells_h(cell);
      cell_local_ids[cell] = local_cell;
    }
  }

  // Allocate and fill the cell nodes
  {
    // Double check this
    TEUCHOS_ASSERT(partition.cell_nodes.rank == 3);

    // Grab the size of the cell node array
    const int num_partition_cells = partition.cell_nodes.extent(0);
    const int num_nodes_per_cell = partition.cell_nodes.extent(1);
    const int num_dims_per_node = partition.cell_nodes.extent(2);

    // Make sure there isn't some strange problem going on
    TEUCHOS_ASSERT(num_partition_cells == num_cells);
    TEUCHOS_ASSERT(num_nodes_per_cell > 0);
    TEUCHOS_ASSERT(num_dims_per_node > 0);

    // Allocate the worksets copy of the cell nodes 
    MDFieldArrayFactory af("",true);
    cell_node_coordinates = af.buildStaticArray<double, Cell, NODE, Dim>("cell nodes", num_partition_cells, num_nodes_per_cell, num_dims_per_node);

    // Copy nodes over
    const auto partition_nodes = partition.cell_nodes;
    auto cnc = cell_node_coordinates.get_view();
    Kokkos::parallel_for(num_cells, KOKKOS_LAMBDA (int i) {
      for(int j=0;j<num_nodes_per_cell;++j)
        for(int k=0;k<num_dims_per_node;++k)
          cnc(i,j,k) = partition_nodes(i,j,k);
      });
    Kokkos::fence();
  }

  // Add the subcell connectivity
  if(partition.has_connectivity){
    auto face_connectivity = Teuchos::rcp(new FaceConnectivity);
    face_connectivity->setup(partition);
    face_connectivity_ = face_connectivity;
  }
  // We have enough information to construct Basis/Point/Integration Values on the fly
  setup_ = true;

}

bool
WorksetDetails::
hasSubcellConnectivity(const unsigned int subcell_dimension) const
{
  TEUCHOS_ASSERT(setup_);
  return (subcell_dimension == (numDimensions() - 1)) and (not face_connectivity_.is_null());
}


const SubcellConnectivity &
WorksetDetails::
getSubcellConnectivity(const unsigned int subcell_dimension) const
{
  TEUCHOS_ASSERT(setup_);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(not hasSubcellConnectivity(subcell_dimension),
                              "Workset::getSubcellConnectivity : Requested subcell dimension "<<subcell_dimension<<" for a "<<num_dimensions_<<"D workset. This is not supported.");
  // Right now we have only one option
  return *face_connectivity_;
}

const panzer::SubcellConnectivity &
WorksetDetails::getFaceConnectivity() const
{
  TEUCHOS_ASSERT(face_connectivity_ != Teuchos::null);
  return *face_connectivity_;
}

const panzer::IntegrationValues2<double> &
WorksetDetails::
getIntegrationValues(const panzer::IntegrationDescriptor & description,
                     const bool lazy_version) const
{
  TEUCHOS_ASSERT(setup_);

  // We need unique keys for the lazy copy or else we get some weird behavior
  size_t key = description.getKey();
  if(lazy_version)
    panzer::hash_combine<int>(key, 123);

  // Check if exists
  const auto itr = integration_values_map_.find(key);
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
  auto ir = Teuchos::rcp(new IntegrationRule(description, cell_topology_, numCells(), num_faces));

  // Create the integration values object
  Teuchos::RCP<IntegrationValues2<double>> iv;
  if(lazy_version){
    iv = Teuchos::rcp(new IntegrationValues2<double>());

    iv->setup(ir,getCellNodes(),numCells());

    // Surface integration schemes need to properly "permute" their entries to line up the surface points between cells
    if(description.getType() == panzer::IntegrationDescriptor::SURFACE)
      iv->setupPermutations(face_connectivity_, numVirtualCells());

  } else {

    iv = Teuchos::rcp(new IntegrationValues2<double>("",true));
    iv->setupArrays(ir);
    iv->evaluateValues(getCellNodes(), numCells(), face_connectivity_, numVirtualCells());

    // This is an advanced feature that requires changes to the workset construction
    // Basically there needs to be a way to grab the side assembly for both "details" belonging to the same workset, which requires a refactor
    // Alternatively, we can do this using a face connectivity object, but the refactor is more important atm.
    TEUCHOS_ASSERT(not (options_.side_assembly_ and options_.align_side_points_));

  }

  integration_values_map_[key] = iv;
  ir_degrees->push_back(iv->int_rule->cubature_degree);
  int_rules.push_back(iv);

  return *iv;

}

const panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & description,
               const bool lazy_version) const
{
  TEUCHOS_ASSERT(setup_);

  // Check if one exists (can be of any integration order)
  const auto itr = basis_integration_values_map_.find(description.getKey());
  if(itr != basis_integration_values_map_.end()){
    for(const auto & pr : itr->second)
      return *pr.second;
  }

  // TODO: We currently overlap BasisIntegrationValues and BasisValues
  // To counter this we create a fake integration rule at this point to ensure the basis values exist

  IntegrationDescriptor id(2*description.getOrder(), IntegrationDescriptor::VOLUME);

  // We have to have the right integrator if this is a side workset
  if(options_.side_assembly_){
    TEUCHOS_ASSERT(getSubcellIndex() >= 0);
    id = IntegrationDescriptor(2*description.getOrder(), IntegrationDescriptor::SIDE, getSubcellIndex());
  }

  // Now just use the other call
  return getBasisValues(description, id, lazy_version);

}


panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & basis_description,
               const panzer::IntegrationDescriptor & integration_description,
               const bool lazy_version) const
{
  TEUCHOS_ASSERT(setup_);

  // We need unique keys for the lazy copy or else we get some weird behavior
  size_t basis_key = basis_description.getKey();
  if(lazy_version)
    panzer::hash_combine<int>(basis_key, 123);

  // We need unique keys for the lazy copy or else we get some weird behavior
  size_t integration_key = integration_description.getKey();
  if(lazy_version)
    panzer::hash_combine<int>(integration_key, 123);

  // Check if exists
  const auto itr = basis_integration_values_map_.find(basis_key);
  if(itr != basis_integration_values_map_.end()){
    const auto & submap = itr->second;
    const auto itr2 = submap.find(integration_key);
    if(itr2 != submap.end())
      return *(itr2->second);

  }

  // Get the integration values for this description
  const auto & iv = getIntegrationValues(integration_description,lazy_version);
  auto bir = Teuchos::rcp(new BasisIRLayout(basis_description.getType(), basis_description.getOrder(), *iv.int_rule));

  Teuchos::RCP<BasisValues2<double>> biv;

  if(lazy_version){

    // Initialized for lazy evaluation

    biv = Teuchos::rcp(new BasisValues2<double>());

    if(integration_description.getType() == IntegrationDescriptor::VOLUME)
      biv->setupUniform(bir, iv.getUniformCubaturePointsRef(false), iv.getJacobian(false), iv.getJacobianDeterminant(false), iv.getJacobianInverse(false));
    else
      biv->setup(bir, iv.getCubaturePointsRef(false), iv.getJacobian(false), iv.getJacobianDeterminant(false), iv.getJacobianInverse(false));

    // pull out the subset of orientations required for this workset
    std::vector<Intrepid2::Orientation> workset_orientations;
    buildLocalOrientations(numOwnedCells()+numGhostCells(),getLocalCellIDs(),options_.orientations_, workset_orientations);

    biv->setOrientations(workset_orientations, numOwnedCells()+numGhostCells());
    biv->setWeightedMeasure(iv.getWeightedMeasure(false));
    biv->setCellNodeCoordinates(cell_node_coordinates);

  } else {

    // Standard, fully allocated version of BasisValues2

    biv = Teuchos::rcp(new BasisValues2<double>("", true, true));
    biv->setupArrays(bir);
    if((integration_description.getType() == IntegrationDescriptor::CV_BOUNDARY) or
       (integration_description.getType() == IntegrationDescriptor::CV_SIDE) or
       (integration_description.getType() == IntegrationDescriptor::CV_VOLUME)){

      biv->evaluateValuesCV(iv.ref_ip_coordinates,
                            iv.jac,
                            iv.jac_det,
                            iv.jac_inv,
                            getCellNodes(),
                            true,
                            numCells());
    } else {

      if(integration_description.getType() == IntegrationDescriptor::VOLUME){

        // TODO: Eventually we will use the other call, however, that will be part of the BasisValues2 refactor
        // The reason we don't do it now is that there are small differences (machine precision) that break EMPIRE testing
        biv->evaluateValues(iv.cub_points,
                            iv.jac,
                            iv.jac_det,
                            iv.jac_inv,
                            iv.weighted_measure,
                            getCellNodes(),
                            true,
                            numCells());

      } else {

        biv->evaluateValues(iv.ref_ip_coordinates,
                            iv.jac,
                            iv.jac_det,
                            iv.jac_inv,
                            iv.weighted_measure,
                            getCellNodes(),
                            true,
                            numCells());
      }
    }

    applyBV2Orientations(numOwnedCells()+numGhostCells(),*biv,getLocalCellIDs(),options_.orientations_);

  }

  basis_integration_values_map_[basis_key][integration_key] = biv;
  bases.push_back(biv);
  basis_names->push_back(biv->basis_layout->name());

  return *biv;

}


panzer::PointValues2<double> &
WorksetDetails::
getPointValues(const panzer::PointDescriptor & description) const
{
  TEUCHOS_ASSERT(setup_);


  // Check if exists
  const auto itr = point_values_map_.find(description.getKey());
  if(itr != point_values_map_.end())
    return *(itr->second);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(not description.hasGenerator(),
                              "Point Descriptor of type '"<<description.getType()<<"' does not have associated generator.");

  auto pr = Teuchos::rcp(new PointRule(description, cell_topology_, numCells()));

  auto pv = Teuchos::rcp(new PointValues2<double>("",true));

  pv->setupArrays(pr);

  // Point values are not necessarily set at the workset level, but can be set by evaluators
  if(description.hasGenerator())
    if(description.getGenerator().hasPoints(*cell_topology_))
      pv->evaluateValues(getCellNodes(), description.getGenerator().getPoints(*cell_topology_),false, numCells());

  point_values_map_[description.getKey()] = pv;

  return *pv;

}

const panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & basis_description,
               const panzer::PointDescriptor & point_description,
               const bool lazy_version) const
{
  TEUCHOS_ASSERT(setup_);

  // We need unique keys for the lazy copy or else we get some weird behavior
  size_t basis_key = basis_description.getKey();
  if(lazy_version)
    panzer::hash_combine<int>(basis_key, 123);

  // Check if exists
  const auto itr = basis_point_values_map_.find(basis_key);
  if(itr != basis_point_values_map_.end()){
    const auto & submap = itr->second;
    const auto itr2 = submap.find(point_description.getKey());
    if(itr2 != submap.end())
      return *(itr2->second);

  }

  // Get the integration values for this description
  const auto & pv = getPointValues(point_description);

  auto bir = Teuchos::rcp(new BasisIRLayout(basis_description.getType(), basis_description.getOrder(), *pv.point_rule));

  Teuchos::RCP<BasisValues2<double>> bpv;

  if(lazy_version){

    // Initialized for lazy evaluation

    bpv = Teuchos::rcp(new BasisValues2<double>());

    bpv->setupUniform(bir, pv.coords_ref, pv.jac, pv.jac_det, pv.jac_inv);

    // pull out the subset of orientations required for this workset
    std::vector<Intrepid2::Orientation> workset_orientations;
    buildLocalOrientations(numOwnedCells()+numGhostCells(),getLocalCellIDs(),options_.orientations_, workset_orientations);

    bpv->setOrientations(workset_orientations, numOwnedCells()+numGhostCells());
    bpv->setCellNodeCoordinates(cell_node_coordinates);

  } else {

    // Standard fully allocated version

    bpv = Teuchos::rcp(new BasisValues2<double>("", true, false));
    bpv->setupArrays(bir);
    bpv->evaluateValues(pv.coords_ref,
                        pv.jac,
                        pv.jac_det,
                        pv.jac_inv,
                        numCells());

    // TODO: We call this separately due to how BasisValues2 is structured - needs to be streamlined
    bpv->evaluateBasisCoordinates(getCellNodes(),numCells());

    applyBV2Orientations(numOwnedCells()+numGhostCells(),*bpv, getLocalCellIDs(), options_.orientations_);

  }

  basis_point_values_map_[basis_key][point_description.getKey()] = bpv;

  return *bpv;

}

const panzer::IntegrationRule &
WorksetDetails::
getIntegrationRule(const panzer::IntegrationDescriptor & description) const
{
  const auto itr = _integration_rule_map.find(description.getKey());
  if(itr == _integration_rule_map.end()){

    // Must run setup or cell topology wont be set properly
    TEUCHOS_ASSERT(setup_);

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

    auto ir = Teuchos::rcp(new IntegrationRule(description, cell_topology_, numCells(), num_faces));

    _integration_rule_map[description.getKey()] = ir;

    return *ir;
  }
  return *(itr->second);
}

const panzer::PureBasis &
WorksetDetails::
getBasis(const panzer::BasisDescriptor & description) const
{
  const auto itr = _pure_basis_map.find(description.getKey());
  if(itr == _pure_basis_map.end()){

    // Must run setup or cell topology wont be set properly
    TEUCHOS_ASSERT(setup_);

    // Create and storethe pure basis
    Teuchos::RCP<panzer::PureBasis> basis = Teuchos::rcp(new panzer::PureBasis(description, cell_topology_, numCells()));
    _pure_basis_map[description.getKey()] = basis;
    return *basis;
  }
  return *(itr->second);
}


void
WorksetDetails::
setNumberOfCells(const int o_cells,
                 const int g_cells,
                 const int v_cells)
{
  num_owned_cells_ = o_cells;
  num_ghost_cells_ = g_cells;
  num_virtual_cells_ = v_cells;
  num_cells = o_cells + g_cells + v_cells;
}

std::ostream&
operator<<(std::ostream& os,
           const panzer::Workset& w)
{
  using std::endl;

  os << "Workset" << endl;
  os << "  block_id=" << w.getElementBlock() << endl;
  os << "  num_cells:" << w.num_cells << endl;
  os << "  num_owned_cells:" << w.numOwnedCells() << endl;
  os << "  num_ghost_cells:" << w.numGhostCells() << endl;
  os << "  num_virtual_cells:" << w.numVirtualCells() << endl;
  os << "  cell_local_ids (size=" << w.getLocalCellIDs().size() << ")" << endl;
  os << "  subcell_dim = " << w.getSubcellDimension() << endl;
  os << "  subcell_index = " << w.getSubcellIndex() << endl;

  os << "  ir_degrees: " << endl;
  for (std::vector<int>::const_iterator ir = w.ir_degrees->begin();
 ir != w.ir_degrees->end(); ++ir)
    os << "    " << *ir << std::endl;

  std::vector<int>::const_iterator ir = w.ir_degrees->begin();
  for (std::vector<Teuchos::RCP<panzer::IntegrationValues2<double> > >::const_iterator irv = w.int_rules.begin();
       irv != w.int_rules.end(); ++irv,++ir) {

    os << "  IR Values (Degree=" << *ir << "):" << endl;

    os << "    cub_points:" << endl;
    os << (*irv)->cub_points << endl;

    os << "    side_cub_points:" << endl;
    os << (*irv)->side_cub_points << endl;

    os << "    cub_weights:" << endl;
    os << (*irv)->cub_weights << endl;

    os << "    node_coordinates:" << endl;
    os << (*irv)->node_coordinates << endl;

    os << "    jac:" << endl;
    os << (*irv)->jac << endl;

    os << "    jac_inv:" << endl;
    os << (*irv)->jac_inv << endl;

    os << "    jac_det:" << endl;
    os << (*irv)->jac_det << endl;

    os << "    weighted_measure:" << endl;
    os << (*irv)->weighted_measure << endl;

    os << "    covarient:" << endl;
    os << (*irv)->covarient << endl;

    os << "    contravarient:" << endl;
    os << (*irv)->contravarient << endl;

    os << "    norm_contravarient:" << endl;
    os << (*irv)->norm_contravarient << endl;

    os << "    ip_coordinates:" << endl;
    os << (*irv)->ip_coordinates << endl;

    os << "    int_rule->getName():" << (*irv)->int_rule->getName() << endl;
  }


  os << "  basis_names: " << endl;
  for (std::vector<std::string>::const_iterator b = w.basis_names->begin();
 b != w.basis_names->end(); ++b)
    os << "    " << *b << std::endl;

  std::vector<std::string>::const_iterator b = w.basis_names->begin();

  for (std::vector<Teuchos::RCP< panzer::BasisValues2<double> > >::const_iterator bv = w.bases.begin(); bv != w.bases.end(); ++bv,++b) {

    os << "  Basis Values (basis_name=" << *b << "):" << endl;

/*
    os << "    basis_ref:" << endl;
    os << (*bv)->basis_ref << endl;

    os << "    basis:" << endl;
    os << (*bv)->basis_scalar << endl;

    os << "    grad_basis_ref:" << endl;
    os << (*bv)->grad_basis_ref << endl;

    os << "    grad_basis:" << endl;
    os << (*bv)->grad_basis << endl;

    os << "    curl_basis_ref:" << endl;
    os << (*bv)->curl_basis_ref_vector << endl;

    os << "    curl_basis:" << endl;
    os << (*bv)->curl_basis_vector << endl;

    os << "    basis_coordinates_ref:" << endl;
    os << (*bv)->basis_coordinates_ref << endl;

    os << "    basis_coordinates:" << endl;
    os << (*bv)->basis_coordinates << endl;
*/

    os << "    basis_layout->name():" << (*bv)->basis_layout->name() << endl;
  }



  return os;
}

}
