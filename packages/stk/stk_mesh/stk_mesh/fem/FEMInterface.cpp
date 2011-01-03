#include <stdexcept>
#include <sstream>

#include <stk_mesh/fem/FEMInterface.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

namespace stk {
namespace mesh {
namespace fem {

void
set_fem_interface(
  MetaData &            meta_data,
  FEMInterface *        fem_interface)
{
  meta_data.declare_attribute_no_delete<FEMInterface>(fem_interface);
}


FEMInterface &
get_fem_interface(
  const MetaData &      meta_data)
{
  FEMInterface *fem_interface = const_cast<FEMInterface *>(meta_data.get_attribute<FEMInterface>());
  ThrowErrorMsgIf( !fem_interface, "FEMInterface plugin must be defined" );

  return *fem_interface;
}


FEMInterface &
get_fem_interface(
  const BulkData &      bulk_data)
{
  return get_fem_interface(bulk_data.mesh_meta_data());
}


FEMInterface &
get_fem_interface(
  const Part &        part)
{
  return get_fem_interface(part.mesh_meta_data());
}


FEMInterface &
get_fem_interface(
  const Bucket &      bucket)
{
  return get_fem_interface(bucket.mesh().mesh_meta_data());
}


FEMInterface &
get_fem_interface(
  const Entity &      entity)
{
  return get_fem_interface(entity.bucket().mesh().mesh_meta_data());
}


void
set_spatial_dimension(
  const MetaData &      meta_data, 
  size_t                spatial_dimension)
{
  FEMInterface &fem = get_fem_interface(meta_data);
  
  fem.set_spatial_dimension(spatial_dimension);
}


void
register_cell_topology(
  MetaData &            meta_data,
  const CellTopology    cell_topology,
  EntityRank            entity_rank)
{
  FEMInterface &fem = get_fem_interface(meta_data);

  fem.register_cell_topology(cell_topology, entity_rank);
}


void
set_cell_topology(
  Part &                part,
  CellTopology          cell_topology)
{
  FEMInterface &fem = get_fem_interface(part.mesh_meta_data());

  fem.set_cell_topology(part, cell_topology);
}


CellTopology
get_cell_topology(
  const Entity &        entity )
{
  return get_cell_topology( entity.bucket());
}


CellTopology
get_cell_topology(
  const Part &          part )
{
  FEMInterface &fem = get_fem_interface(part.mesh_meta_data());

  CellTopology cell_topology(fem.get_cell_topology( part ));

  if (!cell_topology.getCellTopologyData()) {
    const PartVector & supersets = part.supersets();

    for ( PartVector::const_iterator it = supersets.begin(); it != supersets.end() ; ++it ) {
      CellTopology top = fem.get_cell_topology( *(*it));

      if ( ! cell_topology.getCellTopologyData() ) {
        cell_topology = top ;
      }
      else {
        ThrowErrorMsgIf( top.getCellTopologyData() && top != cell_topology,
          "Cell topology is ambiguously defined. It is defined as " << cell_topology.getName() <<
          " on part " << part.name() << " and as " << top.getName() << " on its superset part " << (*it)->name() );
      }
    }
  }

  return cell_topology ;
}


CellTopology
get_cell_topology(
  const Bucket &                bucket)
{
  const BulkData   & bulk_data = bucket.mesh();
  const MetaData   & meta_data = bulk_data.mesh_meta_data();
  const PartVector & all_parts = meta_data.get_parts();

  FEMInterface &fem = get_fem_interface(meta_data);

  CellTopology cell_topology;

  const std::pair< const unsigned *, const unsigned * > supersets = bucket.superset_part_ordinals();

  if (supersets.first != supersets.second) {
    const Part *first_found_part = 0;

    for ( const unsigned * it = supersets.first ; it != supersets.second ; ++it ) {

      const Part & part = * all_parts[*it] ;

      if ( part.primary_entity_rank() == bucket.entity_rank() ) {

        CellTopology top = fem.get_cell_topology( part );

        if ( ! cell_topology.getCellTopologyData() ) {
          cell_topology = top ;

          if (!first_found_part)
            first_found_part = &part;
        }
        else {
          ThrowErrorMsgIf( top.getCellTopologyData() && top != cell_topology,
            "Cell topology is ambiguously defined. It is defined as " << cell_topology.getName() <<
            " on part " << first_found_part->name() << " and as " << top.getName() << " on its superset part " << part.name() );
        }
      }
    }
  }

  return cell_topology ;
}


Part &
get_part(
  const MetaData &      meta_data,
  const CellTopology    cell_topology)
{
  FEMInterface &fem = get_fem_interface(meta_data);

  return fem.get_part(cell_topology);
}


EntityRank
get_entity_rank(
  const MetaData &      meta_data,
  CellTopology          cell_topology)
{
  FEMInterface &fem = get_fem_interface(meta_data);

  return fem.get_entity_rank(cell_topology);
}


std::vector<std::string>
entity_rank_names(
  size_t                spatial_dimension )
{
  ThrowInvalidArgMsgIf( spatial_dimension < 1 || 3 < spatial_dimension,
                        "Invalid spatial dimension = " << spatial_dimension );

  std::vector< std::string > names ;

  names.reserve( spatial_dimension + 1 );

  names.push_back( std::string( "NODE" ) );

  if ( 1 < spatial_dimension ) { names.push_back( std::string("EDGE") ); }
  if ( 2 < spatial_dimension ) { names.push_back( std::string("FACE") ); }

  names.push_back( std::string("ELEMENT") );

  return names ;
}

} // namespace fem



Part &
declare_part(
  MetaData &            meta_data,
  const std::string &   name)
{
  return meta_data.declare_part(name);
}


Part &
declare_part(
  MetaData &            meta_data,
  const std::string &   name,
  EntityRank            entity_rank)
{
  return meta_data.declare_part(name, entity_rank);
}


Part &
declare_part(
  MetaData &            meta_data,
  const PartVector &    name)
{
  return meta_data.declare_part(name);
}


Part &
declare_part(
  MetaData &            meta_data,
  const std::string &   name,
  fem::CellTopology     cell_topology)
{
  fem::FEMInterface &fem_interface = fem::get_fem_interface(meta_data);

  // Get entity rank of cell topology from FEM interface
  EntityRank entity_rank = fem_interface.get_entity_rank(cell_topology);

  ThrowErrorMsgIf( entity_rank == fem::INVALID_RANK,
    "No entity rank has been defined for cell topology " << cell_topology.getName() <<
    "; FEMInterface::declare_cell_topology() defines these");

  Part & part = meta_data.declare_part(name, entity_rank);

  fem_interface.set_cell_topology(part, cell_topology);

  return part ;
}

} // namespace mesh
} // namespace stk
