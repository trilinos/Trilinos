// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "STKElem.hpp"

namespace stk {
namespace transfer {

namespace STKElemUtil {

unsigned parametric(std::vector<double> &para_coords,
                    const double *to,
                    const mesh::Entity element,
                    const mesh::FieldBase &coords_field,
                    const mesh::BulkData& bulkData)
{
  const unsigned numCells = 1; // FIXME
  unsigned spatial_dimension = bulkData.mesh_meta_data().spatial_dimension();

  unsigned found_it = 0;

  MDArray input_phy_points (1,spatial_dimension);
  for (unsigned i=0; i<spatial_dimension; ++i) input_phy_points(0,i) = to[i];

  const mesh::Bucket & bucket = bulkData.bucket(element);
  const CellTopologyData * const bucket_cell_topo_data = mesh::get_cell_topology(bucket).getCellTopologyData();
  ThrowRequireMsg (bucket_cell_topo_data, __FILE__<<":"<<__LINE__<<" parametric::bogus topology");

  shards::CellTopology topo(bucket_cell_topo_data);
  const unsigned numNodes = topo.getNodeCount();

  mesh::Entity const* elem_node_rels = bulkData.begin_nodes(element);
  const unsigned num_nodes = bulkData.num_nodes(element);

  ThrowRequireMsg (topo.getDimension() == spatial_dimension,__FILE__<<":"<<__LINE__<<" Wrong spatical spatial_dimension"
    <<" for topology. Expected "<<spatial_dimension<<" found "<<topo.getDimension());
  ThrowRequireMsg (numNodes == num_nodes ,
    __FILE__<<":"<<__LINE__<<" Expected "<<numNodes<<" nodes but found "<<num_nodes);

  /// FIXME -- fill cellWorkset
  MDArray cellWorkset(numCells, numNodes, spatial_dimension);
  for (unsigned iCell = 0; iCell < numCells; iCell++) {
    for (unsigned iNode = 0; iNode < numNodes; iNode++) {
      const mesh::Entity node = elem_node_rels[iNode];
      const double * coords = static_cast<const double*>(stk::mesh::field_data(coords_field, node));
      for (unsigned iDim=0; iDim < spatial_dimension; iDim++) cellWorkset(iCell, iNode, iDim) = coords[iDim];
    }
  }

  MDArray parametric_coordinates(1,spatial_dimension);
  const unsigned cellOrd = 0;  // FIXME
  Intrepid::CellTools<double>::mapToReferenceFrame(parametric_coordinates,
                                                   input_phy_points,
                                                   cellWorkset,
                                                   topo,
                                                   cellOrd);
  MDArrayUInt inclusion_results(1);  // FIXME
  const double threshold = 1.e-4; // (INTREPID_THRESHOLD default = 10*double_eps ~ 20e-16)
  Intrepid::CellTools<double>::checkPointwiseInclusion(inclusion_results,
                                                       parametric_coordinates,
                                                       topo,
                                                       threshold);
  found_it = inclusion_results(0);
  if (found_it) {
    // for testing only
    if (0) {
      MDArray images(1, spatial_dimension );
      //Intrepid::CellTools<double>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);
      Intrepid::CellTools<double>::mapToPhysicalFrame(images,
                                                      parametric_coordinates,
                                                      cellWorkset,
                                                      topo,
                                                      cellOrd);
    }
  }
  para_coords.resize(spatial_dimension);
  for (unsigned i=0; i<spatial_dimension; ++i) para_coords[i] = parametric_coordinates(0,i);
  return found_it;
}

void parametric(std::vector<std::vector<double> > &val,
          const std::vector<double>               &para_coords,
          const mesh::Entity                       element,
          const std::vector<mesh::FieldBase*>     &values_field,
          const mesh::BulkData&                   bulkData)
{
  unsigned spatial_dimension = bulkData.mesh_meta_data().spatial_dimension();

  const unsigned numIntrp   = 1; // FIXME

  mesh::Entity const* elem_node_rels = bulkData.begin_nodes(element);
  const unsigned num_nodes = bulkData.num_nodes(element);

  const mesh::Bucket & elem_bucket = bulkData.bucket(element);
  const CellTopologyData * const bucket_cell_topo_data = mesh::get_cell_topology(elem_bucket).getCellTopologyData();
  ThrowRequireMsg (bucket_cell_topo_data, __FILE__<<":"<<__LINE__<<" parametric::bogus topology");

  const shards::CellTopology topo(bucket_cell_topo_data);
  const unsigned numNodes = topo.getNodeCount();

  ThrowRequireMsg (topo.getDimension() == spatial_dimension,__FILE__<<":"<<__LINE__<<" Wrong spatical spatial_dimension"
    <<" for topology. Expected "<<spatial_dimension<<" found "<<topo.getDimension());
  ThrowRequireMsg (numNodes == num_nodes ,
    __FILE__<<":"<<__LINE__<<" Expected "<<numNodes<<" nodes but found "<<num_nodes);


  MDArray refPoints (numIntrp,  spatial_dimension);
  MDArray refVals   (           num_nodes, numIntrp);


  for (unsigned intrp = 0; intrp < numIntrp; intrp++)
    for (unsigned i=0; i<spatial_dimension; ++i) refPoints(intrp,i) = para_coords[i];

  fill_ref_vals(refVals, refPoints, topo);

  for (unsigned intrp = 0; intrp < numIntrp; intrp++) {
    const unsigned num_values = values_field.size();
    val.resize(num_values);
    for (unsigned ival=0; ival < num_values; ++ival) {

      const mesh::FieldBase &field = *values_field[ival];
      const mesh::Bucket & node_bucket = bulkData.bucket(elem_node_rels[0]);
      const unsigned bytes = field_bytes_per_entity(field, node_bucket);
      const unsigned bytes_per_entry = field.data_traits().size_of;
      const unsigned num_entry = bytes/bytes_per_entry;

      ThrowRequireMsg (bytes == num_entry * bytes_per_entry,
         __FILE__<<":"<<__LINE__<<" Error:" <<"  bytes:" <<bytes<<"  num_entry:" <<num_entry
              <<"  bytes_per_entry:" <<bytes_per_entry);

      val[ival].resize(num_entry);
      MDArray basis     (num_entry,  num_nodes, numIntrp);
      MDArray dataVals  (num_entry,  num_nodes);
      MDArray fieldVals (num_entry,  numIntrp);
      fieldVals.initialize(0.0);

      // transfer reference basis values to physical frame values
      Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(basis, refVals);
      for (unsigned iNode = 0; iNode < num_nodes; iNode++) {
        const mesh::Entity node = elem_node_rels[iNode];
        const double * values = static_cast<const double*>(stk::mesh::field_data(field, node));
        for (unsigned e = 0; e < num_entry; ++e) dataVals(e, iNode) = values[e];
      }
      // evaluate function at specified points
      Intrepid::FunctionSpaceTools::evaluate<double>(fieldVals, dataVals, basis);
      for (unsigned e = 0; e < num_entry; ++e) val[ival][e] = fieldVals(e, intrp);
    }
  }
}

BasisTable setupBasisTable()
{
  BasisTable basisTable;
  basisTable[shards::getCellTopologyData<shards::Line<2> >()-> key]               .reset( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Triangle<3> >()-> key]           .reset( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Triangle<6> >()-> key]           .reset( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Quadrilateral<4> >()-> key]      .reset( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Quadrilateral<9> >()-> key]      .reset( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Hexahedron<8> >()-> key]         .reset( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Hexahedron<27> >()-> key]        .reset( new Intrepid::Basis_HGRAD_HEX_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Tetrahedron<4> >()-> key]        .reset( new Intrepid::Basis_HGRAD_TET_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Tetrahedron<10> >()-> key]       .reset( new Intrepid::Basis_HGRAD_TET_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Wedge<6> >()-> key]              .reset( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellTriangle<3> >()-> key]      .reset( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellTriangle<6> >()-> key]      .reset( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()-> key] .reset( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
  return basisTable;
}

} // namespace STKElemUtil

unsigned  STKElem::value_size(const EntityKey k, const unsigned i) const
{
  const mesh::Entity      elem = entity(k);
  const mesh::FieldBase &field = *m_values_field[i];

  mesh::Entity const* elem_node_rels = m_bulk_data.begin_nodes(elem);
  const mesh::Entity node = elem_node_rels[0];

  const mesh::Bucket    &bucket= m_bulk_data.bucket(node);

  const unsigned bytes = field_bytes_per_entity(field, bucket);
  const unsigned bytes_per_entry = field.data_traits().size_of;
  const unsigned num_entry = bytes/bytes_per_entry;

  ThrowRequireMsg (bytes == num_entry * bytes_per_entry,
    __FILE__<<":"<<__LINE__<<" Error:" <<"  bytes:" <<bytes<<"  num_entry:" <<num_entry
         <<"  bytes_per_entry:" <<bytes_per_entry);
  return  num_entry;
}

STKElem::EntityKeySet STKElem::entity_keys (
  const mesh::BulkData  &bulk_data,
  const       EntityVec &entities)
{
  EntityKeySet entity_keys;
  for (EntityVec::const_iterator e=entities.begin(); e!=entities.end(); ++e) {
    const mesh::EntityKey k = bulk_data.entity_key(*e);
    entity_keys.insert(k);
  }
  return entity_keys;
}

STKElem::STKElem(
          const            EntityVec          &entities,
          const   mesh::FieldBase             &coord,
          const std::vector<mesh::FieldBase*> &val) :
  m_bulk_data         (coord.get_mesh()),
  m_mesh_modified     (false),
  m_comm              (m_bulk_data.parallel()),
  m_entity_keys       (entity_keys(m_bulk_data, entities)),
  m_coordinates_field (coord),
  m_values_field      (val),
  m_entities_currently_ghosted(),
  m_spatial_dimension(m_bulk_data.mesh_meta_data().spatial_dimension())
{
  const std::string name = "Transfer Ghosting";
  m_bulk_data.modification_begin();
  m_transfer_entity_ghosting = &m_bulk_data.create_ghosting(name);
  m_bulk_data.modification_end();
}

void STKElem::bounding_boxes (std::vector< std::pair<Box,EntityProc> > &v) const {

  const int proc_id = parallel_machine_rank(m_comm);
  v.clear();

  for (EntityKeySet::const_iterator k=m_entity_keys.begin(); k!=m_entity_keys.end(); ++k) {
    const EntityKey id = *k;
    Point min_corner, max_corner;

    elem_coord_limits(min_corner, max_corner, id);
    v.push_back( std::make_pair(Box(min_corner, max_corner), EntityProc(id, proc_id)) );
  }
}

void STKElem::copy_entities(
  const EntityProcVec  &keys_to_copy,
  const std::string    &transfer_name)
{
  m_bulk_data.modification_begin();
  {
    mesh::EntityProcVec new_entities_to_copy(keys_to_copy.size());
    for (size_t i=0; i<keys_to_copy.size(); ++i) {
      // convert from EntityProc based on EntityKey to EntityProc based on raw Entity.
      const EntityProc key_proc = keys_to_copy[i];
      const EntityKey       key = key_proc.id();
      const unsigned       proc = key_proc.proc();
      const Entity            e = entity(key);
      const mesh::EntityProc ep( e, proc);
      new_entities_to_copy[i] = ep;
    }
    m_entities_currently_ghosted.insert(m_entities_currently_ghosted.end(),
                                        new_entities_to_copy.begin(),
                                        new_entities_to_copy.end());

    std::sort(m_entities_currently_ghosted.begin(), m_entities_currently_ghosted.end());
    mesh::EntityProcVec::iterator del = std::unique(m_entities_currently_ghosted.begin(), m_entities_currently_ghosted.end());
    m_entities_currently_ghosted.resize(std::distance(m_entities_currently_ghosted.begin(), del));
  }
  {
    m_bulk_data.change_ghosting(*m_transfer_entity_ghosting,
                                m_entities_currently_ghosted);

    std::vector<mesh::EntityKey> receive;
    std::vector<mesh::EntityProc> send;
    m_transfer_entity_ghosting->receive_list( receive );
    m_transfer_entity_ghosting->send_list( send );
  }
  m_mesh_modified = true;
  m_bulk_data.modification_end();
}

void STKElem::update_values ()
{
  std::vector<const mesh::FieldBase *> fields(m_values_field.begin(), m_values_field.end());
  if (m_mesh_modified) {
    // Copy coordinates to the newly ghosted nodes
    m_mesh_modified = false;
    fields.push_back(&m_coordinates_field);
  }
  mesh::communicate_field_data( *m_transfer_entity_ghosting , fields);
  mesh::copy_owned_to_shared  (  m_bulk_data, fields );
}

void STKElem::elem_coord_limits(Point &min_corner, Point &max_corner, const EntityKey k) const {


  for (unsigned j = 0; j < m_spatial_dimension; ++j ) {
    min_corner[j] =  std::numeric_limits<float>::max();
    max_corner[j] = -std::numeric_limits<float>::max();
  }

  const mesh::Entity elem = entity(k);
  //extract elem_node_relations
  mesh::Entity const* elem_node_rels = m_bulk_data.begin_nodes(elem);
  const unsigned num_nodes = m_bulk_data.num_nodes(elem);

  for ( unsigned ni = 0; ni < num_nodes; ++ni ) {
    const mesh::Entity node = elem_node_rels[ni];

    // pointers to real data
    const double * coords =
      static_cast<const double*>(stk::mesh::field_data(m_coordinates_field, node));
    // check max/min
    for ( unsigned j = 0; j < m_spatial_dimension; ++j ) {
      const float c  = coords[j];
      min_corner[j] = std::min(min_corner[j], c);
      max_corner[j] = std::max(max_corner[j], c);
    }
  }
}

double STKElem::parametric_coord(std::vector<double> &coords,
                                 const double *to,
                                 const EntityKey k ) const {
  const mesh::Entity element = entity(k);
  STKElemUtil::parametric( coords, to, element, m_coordinates_field, m_bulk_data);
  double dist = 0;
  for (unsigned i=0; i<m_spatial_dimension; ++i) {
    dist = std::max(dist, std::abs(coords[i]));
  }
  return dist;
}


} // namespace stk
} // namespace transfer
