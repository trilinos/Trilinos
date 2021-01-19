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
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_PointGenerator.hpp"
#include "Panzer_PointValues2.hpp"

#include "Panzer_SubcellConnectivity.hpp"
#include "Panzer_OrientationsInterface.hpp"

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

  // If the object doesn't exist, feel free to skip applying orientations, they aren't needed in some cases (e.g. DG/FV)
  if(orientations_interface.is_null())
    return;
  const int num_cells = local_cell_ids.extent_int(0);
  const auto & local_orientations = *orientations_interface->getOrientations();
  std::vector<Intrepid2::Orientation> workset_orientations(num_cells);
  for(int i=0; i<num_cells; ++i)
    workset_orientations[i] = local_orientations[local_cell_ids[i]];
  basis_values.applyOrientations(workset_orientations,num_cells);
}

void
correctVirtualNormals(const Teuchos::RCP<panzer::IntegrationValues2<double> > iv,
                      const int num_real_cells,
                      const int num_virtual_cells,
                      const Teuchos::RCP<const shards::CellTopology> cell_topology,
                      const SubcellConnectivity & face_connectivity)
{

  if (iv->int_rule->getType() != panzer::IntegrationDescriptor::SURFACE)
    return;

  // IntegrationValues2 doesn't know anything about virtual cells, so it sets up incorrect normals for those.
  // What we want is for the adjoining face of the virtual cell to have normals that are the negated real cell's normals.
  // we correct the normals here:
  const int space_dim      = cell_topology->getDimension();
  const int faces_per_cell = cell_topology->getSubcellCount(space_dim-1);
  const int points          = iv->surface_normals.extent_int(1);
  const int points_per_face = points / faces_per_cell;
  for (int virtual_cell_ordinal=0; virtual_cell_ordinal<num_virtual_cells; virtual_cell_ordinal++)
  {
    const panzer::LocalOrdinal virtual_cell = virtual_cell_ordinal+num_real_cells;
    int virtual_local_face_id = -1; // the virtual cell face that adjoins the real cell
    int face_ordinal = -1;
    for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++)
    {
      face_ordinal = face_connectivity.subcellForCell(virtual_cell, local_face_id);
      if (face_ordinal >= 0)
      {
        virtual_local_face_id = local_face_id;
        break;
      }
    }
    if (face_ordinal >= 0)
    {
      const int first_cell_for_face = face_connectivity.cellForSubcell(face_ordinal, 0);
      const panzer::LocalOrdinal other_side = (first_cell_for_face == virtual_cell) ? 1 : 0;
      const panzer::LocalOrdinal real_cell = face_connectivity.cellForSubcell(face_ordinal,other_side);
      const panzer::LocalOrdinal face_in_real_cell = face_connectivity.localSubcellForSubcell(face_ordinal,other_side);
      TEUCHOS_ASSERT(real_cell < num_real_cells);
      for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++)
      {
        int virtual_cell_point = points_per_face * virtual_local_face_id + point_ordinal;
        int real_cell_point = points_per_face * face_in_real_cell + point_ordinal;
        // the following arrays will be used to produce/store the rotation matrix below
        double normal[3], transverse[3], binormal[3];
        for(int i=0;i<3;i++)
        {
          normal[i]=0.;
          transverse[i]=0.;
          binormal[i]=0.;
        }

        for (int d=0; d<space_dim; d++)
        {
          const auto n_d = iv->surface_normals(real_cell,real_cell_point,d);
          iv->surface_normals(virtual_cell,virtual_cell_point,d) = -n_d;
          normal[d] = -n_d;
        }

        panzer::IntegrationValues2<double>::convertNormalToRotationMatrix(normal,transverse,binormal);

        for(int dim=0; dim<3; ++dim){
          iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,0,dim) = normal[dim];
          iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,1,dim) = transverse[dim];
          iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,2,dim) = binormal[dim];
        }
      }
      // clear the other normals and rotation matrices for the virtual cell:
      for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++)
      {
        if (local_face_id == virtual_local_face_id) continue;
        for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++)
        {
          int point = local_face_id * points_per_face + point_ordinal;
          for (int dim=0; dim<space_dim; dim++)
          {
            iv->surface_normals(virtual_cell,point,dim) = 0.0;
          }
          for(int dim1=0; dim1<3; ++dim1)
          {
            for(int dim2=0; dim2<3; ++dim2)
            {
              iv->surface_rotation_matrices(virtual_cell,point,dim1,dim2) = 0;
            }
          }
        }
      }
    }
  }
}

}

WorksetDetails::
WorksetDetails()
  : num_cells(0)
  , subcell_dim(-1)
  , subcell_index(-1)
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
    cell_local_ids.resize(num_cells,-1);
    for(int cell=0;cell<num_cells;++cell){
      const int local_cell = partition.local_cells(cell);
      cell_local_ids[cell] = local_cell;
    }
  }

  // Allocate and fill the cell vertices
  {
    // Double check this
    TEUCHOS_ASSERT(partition.cell_vertices.Rank == 3);

    // Grab the size of the cell vertices array
    const int num_partition_cells = partition.cell_vertices.extent(0);
    const int num_vertices_per_cell = partition.cell_vertices.extent(1);
    const int num_dims_per_vertex = partition.cell_vertices.extent(2);

    // Make sure there isn't some strange problem going on
    TEUCHOS_ASSERT(num_partition_cells == num_cells);
    TEUCHOS_ASSERT(num_vertices_per_cell > 0);
    TEUCHOS_ASSERT(num_dims_per_vertex > 0);

    // Allocate the worksets copy of the cell vertices
    MDFieldArrayFactory af("",true);
    cell_vertex_coordinates = af.buildStaticArray<double, Cell, NODE, Dim>("cell vertices", num_partition_cells, num_vertices_per_cell, num_dims_per_vertex);

    // Copy vertices over
    const auto & partition_vertices = partition.cell_vertices;
    for(int i=0;i<num_cells;++i)
      for(int j=0;j<num_vertices_per_cell;++j)
        for(int k=0;k<num_dims_per_vertex;++k)
          cell_vertex_coordinates(i,j,k) = partition_vertices(i,j,k);

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

//void
//WorksetDetails::
//setupNeeds(Teuchos::RCP<const shards::CellTopology> cell_topology,
//           const Kokkos::View<double***,PHX::Device> & cell_vertices,
//           const panzer::WorksetNeeds & needs)
//{
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
//    _pure_basis_map[basis_description.getKey()] = basis;
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
//    _integration_rule_map[integration_description.getKey()] = ir;
//
//    // Create and store integration values
//    Teuchos::RCP<panzer::IntegrationValues2<double> > iv = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
//    iv->setupArrays(ir);
//    iv->evaluateValues(cell_vertex_coordinates,num_cells,face_connectivity_);
//    if (iv->int_rule->getType() == panzer::IntegrationDescriptor::SURFACE)
//    {
//
//      // IntegrationValues2 doesn't know anything about virtual cells, so it sets up incorrect normals for those.
//      // What we want is for the adjoining face of the virtual cell to have normals that are the negated real cell's normals.
//      // we correct the normals here:
//      auto num_real_cells = _num_owned_cells + _num_ghost_cells;
//      const int space_dim      = cell_topology->getDimension();
//      const int faces_per_cell = cell_topology->getSubcellCount(space_dim-1);
//      const int points          = iv->surface_normals.extent_int(1);
//      const int points_per_face = points / faces_per_cell;
//      for (int virtual_cell_ordinal=0; virtual_cell_ordinal<_num_virtual_cells; virtual_cell_ordinal++)
//      {
//        const panzer::LocalOrdinal virtual_cell = virtual_cell_ordinal+num_real_cells;
//        int virtual_local_face_id = -1; // the virtual cell face that adjoins the real cell
//        int face_ordinal = -1;
//        for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++)
//        {
//          face_ordinal = face_connectivity_->subcellForCell(virtual_cell, local_face_id);
//          if (face_ordinal >= 0)
//          {
//            virtual_local_face_id = local_face_id;
//            break;
//          }
//        }
//        if (face_ordinal >= 0)
//        {
//          const int first_cell_for_face = face_connectivity_->cellForSubcell(face_ordinal, 0);
//          const panzer::LocalOrdinal other_side = (first_cell_for_face == virtual_cell) ? 1 : 0;
//          const panzer::LocalOrdinal real_cell = face_connectivity_->cellForSubcell(face_ordinal,other_side);
//          const panzer::LocalOrdinal face_in_real_cell = face_connectivity_->localSubcellForSubcell(face_ordinal,other_side);
//          TEUCHOS_ASSERT(real_cell < num_real_cells);
//          for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++)
//          {
//            int virtual_cell_point = points_per_face * virtual_local_face_id + point_ordinal;
//            int real_cell_point = points_per_face * face_in_real_cell + point_ordinal;
//            // the following arrays will be used to produce/store the rotation matrix below
//            double normal[3], transverse[3], binormal[3];
//            for(int i=0;i<3;i++)
//            {
//              normal[i]=0.;
//              transverse[i]=0.;
//              binormal[i]=0.;
//            }
//
//            for (int d=0; d<space_dim; d++)
//            {
//              const auto n_d = iv->surface_normals(real_cell,real_cell_point,d);
//              iv->surface_normals(virtual_cell,virtual_cell_point,d) = -n_d;
//              normal[d] = -n_d;
//            }
//
//            panzer::IntegrationValues2<double>::convertNormalToRotationMatrix(normal,transverse,binormal);
//
//            for(int dim=0; dim<3; ++dim){
//              iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,0,dim) = normal[dim];
//              iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,1,dim) = transverse[dim];
//              iv->surface_rotation_matrices(virtual_cell,virtual_cell_point,2,dim) = binormal[dim];
//            }
//          }
//          // clear the other normals and rotation matrices for the virtual cell:
//          for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++)
//          {
//            if (local_face_id == virtual_local_face_id) continue;
//            for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++)
//            {
//              int point = local_face_id * points_per_face + point_ordinal;
//              for (int dim=0; dim<space_dim; dim++)
//              {
//                iv->surface_normals(virtual_cell,point,dim) = 0.0;
//              }
//              for(int dim1=0; dim1<3; ++dim1)
//              {
//                for(int dim2=0; dim2<3; ++dim2)
//                {
//                  iv->surface_rotation_matrices(virtual_cell,point,dim1,dim2) = 0;
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//    _integrator_map[integration_description.getKey()] = iv;
//
//    // We need to generate a integration rule - basis pair for each basis
//    for(const panzer::BasisDescriptor & basis_description : basis_descriptors){
//
//      // Grab the basis that was pre-calculated
//      const Teuchos::RCP<const panzer::PureBasis> & basis = _pure_basis_map[basis_description.getKey()];
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
//      _basis_map[basis_description.getKey()][integration_description.getKey()] = bv;
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
//    _point_rule_map[point_description.getKey()] = pr;
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
//      const Teuchos::RCP<const panzer::PureBasis> & basis = _pure_basis_map[basis_description.getKey()];
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
//      _basis_map[basis_description.getKey()][point_description.getKey()] = bv;
//    }
//
//  }
//
//}


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

  auto ir = Teuchos::rcp(new IntegrationRule(description, cell_topology_, numCells(), num_faces));

  auto iv = Teuchos::rcp(new IntegrationValues2<double>("",true));
  iv->setupArrays(ir);
  iv->evaluateValues(getCellVertices(), numCells(), face_connectivity_);

  correctVirtualNormals(iv, num_owned_cells_ + num_ghost_cells_, num_virtual_cells_, cell_topology_, *face_connectivity_);

  // This is an advanced feature that requires changes to the workset construction
  // Basically there needs to be a way to grab the side assembly for both "details" belonging to the same workset, which requires a refactor
  // Alternatively, we can do this using a face connectivity object, but the refactor is more important atm.
  TEUCHOS_ASSERT(not (options_.side_assembly_ and options_.align_side_points_));

  integration_values_map_[description.getKey()] = iv;

  return *iv;

}

const panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & description) const
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
  return getBasisValues(description, id);

}


panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & basis_description,
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

  auto biv = Teuchos::rcp(new BasisValues2<double>("", true, true));
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
                          getCellVertices(),
                          true,
                          numCells());

    } else {

      biv->evaluateValues(iv.ref_ip_coordinates,
                          iv.jac,
                          iv.jac_det,
                          iv.jac_inv,
                          iv.weighted_measure,
                          getCellVertices(),
                          true,
                          numCells());
    }
  }

  applyBV2Orientations(*biv,getLocalCellIDs(),options_.orientations_);

  basis_integration_values_map_[basis_description.getKey()][integration_description.getKey()] = biv;

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
      pv->evaluateValues(getCellVertices(), description.getGenerator().getPoints(*cell_topology_),false, numCells());

  point_values_map_[description.getKey()] = pv;

  return *pv;

}

const panzer::BasisValues2<double> &
WorksetDetails::
getBasisValues(const panzer::BasisDescriptor & basis_description,
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

  auto bpv = Teuchos::rcp(new BasisValues2<double>("", true, false));
  bpv->setupArrays(bir);
  bpv->evaluateValues(pv.coords_ref,
                      pv.jac,
                      pv.jac_det,
                      pv.jac_inv,
                      numCells());

  // TODO: We call this separetely due to how BasisValues2 is structured - needs to be streamlined
  bpv->evaluateBasisCoordinates(getCellVertices(),numCells());

  applyBV2Orientations(*bpv, getLocalCellIDs(), options_.orientations_);

  basis_point_values_map_[basis_description.getKey()][point_description.getKey()] = bpv;

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
