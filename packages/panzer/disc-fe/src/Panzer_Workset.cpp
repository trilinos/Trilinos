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

#include "Shards_CellTopology.hpp"

namespace panzer {

void
WorksetDetails::setup(const panzer::LocalMeshPartition<int,panzer::Ordinal64> & partition,
                      const panzer::WorksetNeeds & needs)
{


  const size_t num_cells = partition.local_cells.extent(0);

  _num_owned_cells = partition.num_owned_cells;
  _num_ghost_cells = partition.num_ghstd_cells;
  _num_virtual_cells = partition.num_virtual_cells;

  subcell_index = -1;
  block_id = partition.element_block_name;

  Kokkos::View<int*, PHX::Device> cell_ids = Kokkos::View<int*, PHX::Device>("cell_ids",num_cells);
  Kokkos::deep_copy(cell_ids, partition.local_cells);
  cell_local_ids_k = cell_ids;

  cell_local_ids.resize(num_cells,-1);
  for(size_t cell=0;cell<num_cells;++cell){
    const int local_cell = partition.local_cells(cell);
    cell_local_ids[cell] = local_cell;
  }

  auto fc = Teuchos::rcp(new panzer::FaceConnectivity());
  fc->setup(partition);
  _face_connectivity = fc;

  setupNeeds(partition.cell_topology, partition.cell_vertices, needs);
}

void WorksetDetails::setupNeeds(Teuchos::RCP<const shards::CellTopology> cell_topology,
                                const Kokkos::View<double***,PHX::Device> & cell_vertices,
                                const panzer::WorksetNeeds & needs)
{

  const size_t num_cells = cell_vertices.extent(0);
  const size_t num_vertices_per_cell = cell_vertices.extent(1);
  const size_t num_dims_per_vertex = cell_vertices.extent(2);

  // Set cell vertices
  {

    MDFieldArrayFactory af("",true);

    cell_vertex_coordinates = af.template buildStaticArray<double, Cell, NODE, Dim>("cell_vertices",num_cells, num_vertices_per_cell, num_dims_per_vertex);

    for(size_t i=0;i<num_cells;++i)
      for(size_t j=0;j<num_vertices_per_cell;++j)
        for(size_t k=0;k<num_dims_per_vertex;++k)
          cell_vertex_coordinates(i,j,k) = cell_vertices(i,j,k);

  }

  // DEPRECATED - makes sure deprecated arrays are filled with something - this will probably segfault or throw an error
  panzer::populateValueArrays(num_cells, false, needs, *this);

  const std::vector<panzer::BasisDescriptor> & basis_descriptors = needs.getBases();
  const std::vector<panzer::IntegrationDescriptor> & integration_descriptors = needs.getIntegrators();
  const std::vector<panzer::PointDescriptor> & point_descriptors = needs.getPoints();

  // Create the pure basis
  for(const panzer::BasisDescriptor & basis_description : basis_descriptors){
    // Create and store integration rule
    Teuchos::RCP<panzer::PureBasis> basis = Teuchos::rcp(new panzer::PureBasis(basis_description, cell_topology, num_cells));
    _pure_basis_map[basis_description.getKey()] = basis;
  }

  // Create the integration terms and basis-integrator pairs
  for(const panzer::IntegrationDescriptor & integration_description : integration_descriptors){

    int num_faces = -1;
    if(integration_description.getType() == panzer::IntegrationRule::SURFACE){
      num_faces = getFaceConnectivity().numSubcells();
    }

    // Create and store integration rule
    Teuchos::RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(integration_description, cell_topology, num_cells, num_faces));
    _integration_rule_map[integration_description.getKey()] = ir;

    // Create and store integration values
    Teuchos::RCP<panzer::IntegrationValues2<double> > iv = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
    iv->setupArrays(ir);
    iv->evaluateValues(cell_vertex_coordinates,num_cells);
    _integrator_map[integration_description.getKey()] = iv;

    // We need to generate a integration rule - basis pair for each basis
    for(const panzer::BasisDescriptor & basis_description : basis_descriptors){

      // Grab the basis that was pre-calculated
      const Teuchos::RCP<const panzer::PureBasis> & basis = _pure_basis_map[basis_description.getKey()];

      // Create a basis ir layout for this pair of integrator and basis
      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*ir));

      // Create and store basis values
      Teuchos::RCP<panzer::BasisValues2<double> > bv = Teuchos::rcp(new panzer::BasisValues2<double>("",true,true));
      bv->setupArrays(b_layout);
      if(ir->getType() == panzer::IntegrationDescriptor::SURFACE){
        bv->evaluateValues(iv->ref_ip_coordinates,
                           iv->jac,
                           iv->jac_det,
                           iv->jac_inv,
                           iv->weighted_measure,
                           cell_vertex_coordinates,
                           true,
                           num_cells);
      } else if((ir->getType() == panzer::IntegrationDescriptor::CV_VOLUME)
          or (ir->getType() == panzer::IntegrationDescriptor::CV_SIDE)
          or (ir->getType() == panzer::IntegrationDescriptor::CV_BOUNDARY)){
        bv->evaluateValuesCV(iv->ref_ip_coordinates,
                             iv->jac,
                             iv->jac_det,
                             iv->jac_inv);
      } else {
        bv->evaluateValues(iv->cub_points,
                           iv->jac,
                           iv->jac_det,
                           iv->jac_inv,
                           iv->weighted_measure,
                           cell_vertex_coordinates,
                           true,
                           num_cells);
      }
      _basis_map[basis_description.getKey()][integration_description.getKey()] = bv;
    }

  }

  // Create the point terms and basis-integrator pairs
  for(const panzer::PointDescriptor & point_description : point_descriptors){

    // get the generaotr, and build some points from a topology
    auto points = point_description.getGenerator().getPoints(*cell_topology);

    // Create and store integration rule
    Teuchos::RCP<panzer::PointRule> pr = Teuchos::rcp(new panzer::PointRule(point_description,
                                                                            cell_topology,
                                                                            num_cells));
    _point_rule_map[point_description.getKey()] = pr;

    // Create and store integration values
    Teuchos::RCP<panzer::PointValues2<double> > pv = Teuchos::rcp(new panzer::PointValues2<double>("",true));
    pv->setupArrays(pr);
    pv->evaluateValues(cell_vertex_coordinates,points);

    _point_map[point_description.getKey()] = pv;

    // We need to generate a integration rule - basis pair for each basis
    for(const panzer::BasisDescriptor & basis_description : basis_descriptors){

      // Grab the basis that was pre-calculated
      const Teuchos::RCP<const panzer::PureBasis> & basis = _pure_basis_map[basis_description.getKey()];

      // Create a basis ir layout for this pair of integrator and basis
      Teuchos::RCP<panzer::BasisIRLayout> b_layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*pr));

      // Create and store basis values
      Teuchos::RCP<panzer::BasisValues2<double> > bv = Teuchos::rcp(new panzer::BasisValues2<double>("",true,true));
      bv->setupArrays(b_layout);

      bv->evaluateValues(pv->coords_ref,
                         pv->jac,
                         pv->jac_det,
                         pv->jac_inv);

      _basis_map[basis_description.getKey()][point_description.getKey()] = bv;
    }

  }

}

const panzer::SubcellConnectivity &
WorksetDetails::getFaceConnectivity() const
{
  TEUCHOS_ASSERT(_face_connectivity != Teuchos::null);
  return *_face_connectivity;
}

const panzer::IntegrationValues2<double> &
WorksetDetails::getIntegrationValues(const panzer::IntegrationDescriptor & description) const
{
  const auto itr = _integrator_map.find(description.getKey());
  TEUCHOS_ASSERT(itr != _integrator_map.end());
  return *(itr->second);
}

const panzer::IntegrationRule &
WorksetDetails::getIntegrationRule(const panzer::IntegrationDescriptor & description) const
{
  const auto itr = _integration_rule_map.find(description.getKey());
  TEUCHOS_ASSERT(itr != _integration_rule_map.end());
  return *(itr->second);
}

panzer::BasisValues2<double> &
WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description, 
                               const panzer::IntegrationDescriptor & integration_description)
{
  const auto itr = _basis_map.find(basis_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == _basis_map.end(),
                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" " 
                              "of order " + std::to_string(basis_description.getOrder()));
  const auto & integration_map = itr->second;
  const auto itr2 = integration_map.find(integration_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == integration_map.end(),
                              "Workset::getBasisValues: Can't find integration " + std::to_string(integration_description.getType()) + " " 
                              "of order " + std::to_string(integration_description.getOrder()));
  return *(itr2->second);
}

const panzer::BasisValues2<double> &
WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description, 
                               const panzer::PointDescriptor & point_description) const
{
  const auto itr = _basis_map.find(basis_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == _basis_map.end(),
                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" " 
                              "of order " + std::to_string(basis_description.getOrder()));
  const auto & point_map = itr->second;
  const auto itr2 = point_map.find(point_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == point_map.end(),
                              "Workset::getBasisValues: Can't find point values \"" + point_description.getType() + "\""); 
  return *(itr2->second);
}

const panzer::BasisValues2<double> &
WorksetDetails::getBasisValues(const panzer::BasisDescriptor & basis_description, 
                               const panzer::IntegrationDescriptor & integration_description) const
{
  const auto itr = _basis_map.find(basis_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == _basis_map.end(),
                              "Workset::getBasisValues: Can't find basis \"" + basis_description.getType() + "\" " 
                              "of order " + std::to_string(basis_description.getOrder()));
  const auto & integration_map = itr->second;
  const auto itr2 = integration_map.find(integration_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr2 == integration_map.end(),
                              "Workset::getBasisValues: Can't find integration " + std::to_string(integration_description.getType()) + " " 
                              "of order " + std::to_string(integration_description.getOrder()));
  return *(itr2->second);
}

const panzer::PointValues2<double> &
WorksetDetails::getPointValues(const panzer::PointDescriptor & point_description) const
{
  const auto itr = _point_map.find(point_description.getKey());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == _point_map.end(),
                              "Workset::getPointValues: Can't find point values \"" + point_description.getType() + "\""); 
  return *(itr->second);
}

const panzer::PureBasis &
WorksetDetails::getBasis(const panzer::BasisDescriptor & description) const
{
  const auto itr = _pure_basis_map.find(description.getKey());
  TEUCHOS_ASSERT(itr != _pure_basis_map.end());
  return *(itr->second);
}

  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w)
  {
    using std::endl;

    os << "Workset" << endl;
    os << "  block_id=" << w.block_id << endl;
    os << "  num_cells:" << w.num_cells << endl;
    os << "  cell_local_ids (size=" << w.cell_local_ids.size() << ")" << endl;
    os << "  subcell_dim = " << w.subcell_dim << endl;
    os << "  subcell_index = " << w.subcell_index << endl;

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
