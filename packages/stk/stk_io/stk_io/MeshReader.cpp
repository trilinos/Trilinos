/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <string.h>
#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_io/MeshReader.hpp>
#include <stk_io/IossBridge.hpp>

namespace stk {
namespace io {

typedef mesh::Field<double>                       ScalarField ;
typedef mesh::Field<double,mesh::Cartesian>       CartesianField ;
typedef mesh::Field<double,mesh::FullTensor>      FullTensorField ;
typedef mesh::Field<double,mesh::SymmetricTensor> SymmetricTensorField ;

//----------------------------------------------------------------------

bool MeshReaderFilter::accept_part(mesh::MetaData     & /* meta_data   */ ,
				   const std::string  & /* io_name     */,
				   mesh::EntityRank     /* entity_rank */ ,
				   const CellTopologyData * ) const
{ return true ; }


mesh::FieldBase *
MeshReaderFilter::map_attribute(
  mesh::Part        & part ,
  const std::string & io_name ,
  const std::vector< const shards::ArrayDimTag * > & tags ) const
{
  std::string attr_name ;

  if ( equal_case( io_name , std::string("distribution_factors") ) ) {
    attr_name = part.name();
    attr_name.append("_");
  }

  attr_name.append(io_name);

  mesh::MetaData & meta_data = part.mesh_meta_data();

  mesh::FieldBase * field =
    meta_data.declare_field_base( attr_name , mesh::data_traits<double>() ,
                                  unsigned(tags.size()) ,
				  & tags[0] , 1 );

  return field ;
}

mesh::FieldBase *
MeshReaderFilter::map_transient(
  mesh::Part        & part ,
  const std::string & io_name ,
  const std::vector< const shards::ArrayDimTag * > & tags ) const
{
  mesh::MetaData & meta_data = part.mesh_meta_data();

  mesh::FieldBase * field =
    meta_data.declare_field_base( io_name , mesh::data_traits<double>() ,
                                  unsigned(tags.size()) ,
				  & tags[0] , 1 );

  return field ;
}

MeshReaderFilter::~MeshReaderFilter() {}

MeshReader::~MeshReader() {}

//----------------------------------------------------------------------

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

class IossReader : public MeshReader {
public:

  IossReader( ParallelMachine     parallel ,
              mesh::MetaData    & meta_data ,
              mesh::VectorField & coord_field ,
              const std::string & file_format ,
              const std::string & file_name ,
              Ioss::Region      * ioss_region ,
              const MeshReaderFilter & );

  ~IossReader();


  void read_model( mesh::BulkData & ) const ;

private:
  IossReader(const IossReader& from); // do not implement
  IossReader& operator=(const IossReader& from); // do not implement

  void throw_error( const std::string & ) const ;

  //----------------------------------

  void read_meta_data_space( const MeshReaderFilter & );

  mesh::Part * read_meta_data_part( const MeshReaderFilter & ,
                                    const Ioss::GroupingEntity & ,
                                    mesh::EntityRank ,
                                    const CellTopologyData * const );

  void read_meta_data_block( const MeshReaderFilter & ,
                             const Ioss::EntityBlock & );

  void read_meta_data_side(  const MeshReaderFilter & ,
                             const Ioss::GroupingEntity & ,
                             mesh::EntityRank );

  void read_meta_data_node( const MeshReaderFilter     & ,
                            const Ioss::GroupingEntity & );

  //----------------------------------

  void read_bulk_data_block( mesh::BulkData & ,
                             mesh::Part     & ,
                             Ioss::GroupingEntity & ,
                             std::vector< mesh::Entity * > & ) const ;

  void read_bulk_data_side( mesh::BulkData & ,
                            mesh::Part     & ,
                            Ioss::EntityBlock & ,
                            std::vector< mesh::Entity * > & ) const ;

  void read_bulk_data_node( mesh::BulkData & ,
                            mesh::Part     & ,
                            Ioss::GroupingEntity & ,
                            std::vector< mesh::Entity * > & ) const ;

  void read_bulk_data_attr( const std::vector< mesh::Entity * > & ,
                            mesh::Part & ,
                            Ioss::GroupingEntity & ) const ;

  void read_bulk_data_coordinates( mesh::BulkData & mesh ) const ;

  void read_field( const std::vector< mesh::Entity *> & ,
                   const mesh::FieldBase      & ,
                   const Ioss::GroupingEntity & ,
                   const std::string          & ) const ;


  Ioss::Region * m_ioss_region ;
  unsigned       m_spatial_dim ;
};

//----------------------------------------------------------------------

void IossReader::throw_error( const std::string & msg ) const
{
  std::string ex_msg ;
  ex_msg = "stk::io::IossReader['" + m_file_format + "' '" + m_file_name + "'] FAILED: " + msg;
  throw std::runtime_error( ex_msg );
}

//----------------------------------------------------------------------

  void IossReader::read_meta_data_space( const MeshReaderFilter & /* filter */)
{
  // Exactly one node block to obtain the nodal coordinates dimension

  const Ioss::NodeBlockContainer& node_blocks =
    m_ioss_region->get_node_blocks();

  if ( node_blocks.size() != 1 ) {
    std::ostringstream msg ;
    msg << " 1 != ( node_blocks.size() = " << node_blocks.size() << " )" ;
    throw_error( msg.str() );
  }

  Ioss::NodeBlock & nb = * node_blocks[0];

  m_spatial_dim = nb.get_property("component_degree").get_int();

  stk::mesh::put_field( m_coord_field, mesh::Node,
                         m_meta_data.universal_part(), m_spatial_dim);
}

//----------------------------------------------------------------------

mesh::Part * IossReader::read_meta_data_part(
  const MeshReaderFilter     & filter ,
  const Ioss::GroupingEntity & io_part ,
  mesh::EntityRank             entity_rank ,
  const CellTopologyData * const   celltop )
{
  mesh::Part * part = NULL ;

  if ( filter.accept_part(m_meta_data,io_part.name(),entity_rank,celltop) ) {

    part = & m_meta_data.declare_part( io_part.name() , entity_rank );

    if ( celltop ) { mesh::set_cell_topology( *part , celltop ); }
  }

  return part ;
}

//----------------------------------------------------------------------

void IossReader::read_meta_data_block(
				      const MeshReaderFilter & filter ,
				      const Ioss::EntityBlock & block )
{
  const CellTopologyData * const celltop = stk::io::map_topology_ioss_to_cell( block.topology());

  mesh::Part * const part =
    read_meta_data_part( filter, block, mesh::Element, celltop );

  if (!part)
    return;

  m_parts.push_back( part );

  const int attribute_count =
    block.get_property("attribute_count").get_int() ;

  if ( 0 < attribute_count ) {
    std::ostringstream msg ;
    msg << "ERROR: Attributes are not currently handled in this routine";
    throw_error( msg.str() );
  }
}

//----------------------------------------------------------------------

void IossReader::read_meta_data_side(
				     const MeshReaderFilter     & filter ,
				     const Ioss::GroupingEntity & side ,
				     mesh::EntityRank             entity_rank )
{
  mesh::Part * part = read_meta_data_part( filter, side, entity_rank, NULL );

  if ( !part )
    return;

  m_parts.push_back( part );

  std::pair< stk::mesh::FieldBase * , std::string > attr ;
  attr.first  = NULL ;
  attr.second.assign("distribution_factors");

  const size_t block_count = side.block_count();

  bool has_distribution_factors = false ;

  std::vector< mesh::Part * > side_blocks( block_count );

  for ( size_t i = 0 ; i < block_count ; i++ ) {

    Ioss::EntityBlock & io_block = * side.get_block(i);

    const CellTopologyData * const celltop =
      stk::io::map_topology_ioss_to_cell( io_block.topology());

    side_blocks[i] =
      read_meta_data_part( filter, io_block, entity_rank , celltop );

    if ( side_blocks[i] ) {

      m_meta_data.declare_part_subset( *part , * side_blocks[i] );

      if (io_block.field_exists(attr.second)) {
	has_distribution_factors = true ;
      }
      else {
	side_blocks[i] = NULL ;
      }
    }
  }

  if ( has_distribution_factors ) {

    std::vector< const shards::ArrayDimTag * > tags(1);
    tags[0] = & mesh::ElementNode::tag();

    attr.first = filter.map_attribute( *part , attr.second , tags );

    if ( attr.first ) {

      m_attributes.push_back( attr );

      for ( unsigned i = 0 ; i < side_blocks.size() ; i++ ) {

	if ( side_blocks[i] ) {
	  mesh::Part & block = * side_blocks[i] ;

	  const CellTopologyData * const top = mesh::get_cell_topology( block );
	  if (top == NULL) {
	    std::ostringstream msg ;
	    msg << " INTERNAL_ERROR: Part " << block.name() << " returned NULL from get_cell_topology()";
	    throw std::runtime_error( msg.str() );
	  }

	  const unsigned stride[1] = { top->node_count };

	  m_meta_data.declare_field_restriction(
						*attr.first , block.primary_entity_rank() , block , stride );
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void IossReader::read_meta_data_node(
  const MeshReaderFilter     & filter ,
  const Ioss::GroupingEntity & io_part )
{
  std::pair< stk::mesh::FieldBase * , std::string > attr ;
  attr.first = NULL ;
  attr.second.assign( "distribution_factors" );

  stk::mesh::Part * part =
    read_meta_data_part( filter, io_part, mesh::Node, NULL );

  if ( part ) {

    m_parts.push_back( part );

    const bool has_distribution_factors = io_part.field_exists( attr.second );

    if ( has_distribution_factors ) {

      attr.first = filter.map_attribute( *part , attr.second ,
                                         std::vector<const shards::ArrayDimTag *>() );

      if ( attr.first ) {
        const unsigned stride[1] = { 1 };

        m_meta_data.declare_field_restriction(
         *attr.first, mesh::Node , *part , stride );

        m_attributes.push_back( attr );
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void IossReader::read_bulk_data_block(
  mesh::BulkData & bulk ,
  mesh::Part     & part ,
  Ioss::GroupingEntity & io_part ,
  std::vector< mesh::Entity * > & elements ) const
{
  const CellTopologyData* cell_topo = mesh::get_cell_topology(part);
  if (cell_topo == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
    throw std::runtime_error( msg.str() );
  }

  // Create elements:
  {
    std::vector<int> elem_ids ;
    std::vector<int> connectivity ;

    io_part.get_field_data("ids", elem_ids);
    io_part.get_field_data("connectivity", connectivity);

    const size_t num_elems = elem_ids.size();
    const int nodes_per_elem = cell_topo->node_count ;

    if ( elem_ids.size() * nodes_per_elem != connectivity.size() ) {
      std::ostringstream msg ;
      msg << "Erroneous connectivity information read for '" ;
      msg << io_part.name();
      msg << "' nodes_per_element( " << nodes_per_elem ;
      msg << " ) * element_count(  " << elem_ids.size();
      msg << " ) != element_node_count( " << connectivity.size();
      msg << " )" ;
      throw_error( msg.str() );
    }

    elements.resize( num_elems );

    for ( size_t i = 0 ; i < num_elems ; ++i ) {
      /// \todo REFACTOR cast from int to unsigned is unsafe and ugly.
      /// change function to take int[] argument.

      int *conn = &connectivity[i*nodes_per_elem];

      elements[i] = & mesh::declare_element( bulk, part, elem_ids[i], conn);
    }
  }
}

//----------------------------------------------------------------------

void IossReader::read_bulk_data_side(
  mesh::BulkData & bulk ,
  mesh::Part     & part ,
  Ioss::EntityBlock & io_part ,
  std::vector<mesh::Entity *> & sides ) const
{
  bool ok_size = true ;

  std::vector<int> side_ids ;
  std::vector<int> elem_side ;

  io_part.get_field_data("ids", side_ids);
  io_part.get_field_data("element_side", elem_side);

  const size_t num_sides = side_ids.size();

  ok_size = num_sides * 2 == elem_side.size();

  if ( ok_size ) {

    sides.resize( num_sides );

    for ( size_t i = 0 ; i < num_sides ; ++i ) {

      // Ioss is 1-based, mesh is 0-based.
      // Hence the '-1' in the following line.

      const int side_local_id  = elem_side[i*2+1] - 1 ;
      const int elem_global_id = elem_side[i*2] ;
      const int side_global_id = side_ids[i] ;

      mesh::Entity * const elem = bulk.get_entity(mesh::Element, elem_global_id);

      sides[i] = ! elem ? NULL
                 : & mesh::declare_element_side( bulk,
                                                 side_global_id,
                                                 *elem ,
                                                 side_local_id ,
                                                 & part );
    }
  }

  if ( ! ok_size ) {
    std::ostringstream msg ;
    msg << " Reading { io_part = '" << io_part.name()
        << "', io_size = " << num_sides
        << " } Data size error" ;
    throw_error( msg.str() );
  }
}

//----------------------------------------------------------------------

void IossReader::read_bulk_data_node(
  mesh::BulkData                & bulk ,
  mesh::Part                    & part ,
  Ioss::GroupingEntity          & io_part ,
  std::vector< mesh::Entity * > & entities ) const
{
  std::vector< mesh::Part * > add( 1 ); add[0] = & part ;

  std::vector<int> node_ids ;

  io_part.get_field_data("ids", node_ids);

  const size_t num_nodes = node_ids.size();

  entities.resize( num_nodes );

  for ( size_t i = 0 ; i < num_nodes ; ++i ) {
    entities[i] = bulk.get_entity(mesh::Node, node_ids[i]);

    if ( entities[i] ) { // Put the node into the part

      bulk.change_entity_parts( * entities[i] , add );
    }
  }
}

//----------------------------------------------------------------------
//  Copy data from input to field.
/// \todo  REFACTOR  What about permuted field data for faces and edges?

template< typename T >
bool copy_field_data( const std::vector< mesh::Entity * > & entities ,
                      const mesh::FieldBase               & field ,
                      const std::vector<double>           & data ,
                      const size_t                          ncomp )
{
  std::vector< mesh::Entity * >::const_iterator ie = entities.begin();
  std::vector< double >::const_iterator id = data.begin();

  for ( ; ie != entities.end() ; ++ie ) {

    if ( NULL != *ie ) {

      mesh::Entity & e = **ie ;

      const size_t n = mesh::field_data_size( field , e );

      if ( n != ncomp * sizeof(T) ) { return false ; }

      T * d = (T *) mesh::field_data( field , e );
      T * const d_end = d + ncomp ;

      while ( d != d_end ) { *d++ = (T) *id++ ; }
    }
  }
  return true ;
}


void IossReader::read_field( const std::vector< mesh::Entity *> & entities ,
                             const mesh::FieldBase      & field ,
                             const Ioss::GroupingEntity & io_part ,
                             const std::string          & io_name ) const
{
  const bool ok_exists = io_part.field_exists( io_name );
  bool ok_size = false ;
  bool ok_copy = false ;

  size_t entity_count = entities.size();
  size_t io_count     = 0 ;
  size_t io_data_size = 0 ;
  size_t io_ncomp     = 0 ;

  if ( ok_exists ) {

    Ioss::Field io_field = io_part.get_field(io_name);

    std::vector<double> io_data ;

    io_count     = io_part.get_field_data( io_name , io_data );
    io_data_size = io_data.size();
    io_ncomp     = io_field.transformed_storage()->component_count();

    ok_size = io_count     == entity_count &&
              io_data_size == entity_count * io_ncomp ;

    if ( ok_size ) {
      if ( field.type_is<double>() ) {
        ok_copy = copy_field_data<double>( entities, field, io_data, io_ncomp);
      }
      else if ( field.type_is<float>() ) {
        ok_copy = copy_field_data<float>( entities, field, io_data, io_ncomp);
      }
      else if ( field.type_is<int>() ) {
        ok_copy = copy_field_data<int>( entities, field, io_data, io_ncomp );
      }
      else if ( field.type_is<unsigned>() ) {
        ok_copy = copy_field_data<unsigned>(entities,field,io_data,io_ncomp );
      }
    }
  }

  if ( ! ok_exists || ! ok_size || ! ok_copy ) {
    std::ostringstream msg ;
    msg << " Reading { io_part = '" << io_part.name()
        << "', io_field = '" << io_name << "'"
        << "', io_size = " << io_data_size
        << ", io_components = " << io_ncomp
        << ", io_item_count = " << io_count
        << " } -> { " << field
        << " , entity_count = " << entity_count
        << " } : " ;

    if      ( ! ok_exists ) { msg << "No source data" ; }
    else if ( ! ok_size )   { msg << "Bad data sizes" ; }
    else                    { msg << "Failed data copy" ; }
    throw_error( msg.str() );
  }
}

//----------------------------------------------------------------------
// Read attributes:

void IossReader::read_bulk_data_attr(
  const std::vector< mesh::Entity * > & entities ,
  mesh::Part & part ,
  Ioss::GroupingEntity & io_part ) const
{
  typedef std::pair< stk::mesh::FieldBase * , std::string > field_map ;

  for ( std::vector< field_map >::const_iterator
        j = m_attributes.begin() ; j != m_attributes.end() ; ++j ) {
    const stk::mesh::EntityRank  type   =
      stk::mesh::fem_entity_rank( part.primary_entity_rank() );
    const stk::mesh::FieldBase & field  = * j->first ;
    const std::string          & ioname = j->second ;

    const stk::mesh::FieldBase::Restriction & r = field.restriction(type,part);

    if ( r.stride[0] ) {
      read_field( entities , field , io_part , ioname );
    }
  }
}


void IossReader::read_bulk_data_coordinates( mesh::BulkData & mesh ) const
{
  // Exactly one node block to obtain the nodal coordinates and ids:

  const Ioss::NodeBlockContainer& node_blocks =
    m_ioss_region->get_node_blocks();

  const bool ok_one_block = node_blocks.size() == 1 ;
  bool ok_data = true ;

  std::vector<int>    node_ids ;
  std::vector<double> node_coords ;

  if ( ok_one_block ) {

   Ioss::NodeBlock & nb = * node_blocks[0];

    nb.get_field_data("ids", node_ids);
    nb.get_field_data("mesh_model_coordinates", node_coords);

    const size_t num_nodes = node_ids.size();

    ok_data = num_nodes * m_spatial_dim == node_coords.size();

    if ( ok_data ) {

      for( unsigned i=0; i<num_nodes; ++i) {
        mesh::Entity * const node = mesh.get_entity( mesh::Node, node_ids[i] );

        if ( NULL != node ) {
          double* coord_data = mesh::field_data(m_coord_field, *node);
	  if (coord_data == NULL) {
	    std::ostringstream msg ;
	    msg << " INTERNAL_ERROR: Node " << node_ids[i]
		<< " does not have a coordinate data entry in read_bulk_data_coordinates()";
	    throw std::runtime_error( msg.str() );
	  }

          for(unsigned j=0; j<m_spatial_dim; ++j) {
            coord_data[j] = node_coords[i*m_spatial_dim+j];
          }
        }
      }
    }
  }

  if ( ! ok_one_block || ! ok_data ) {
    std::ostringstream msg ;
    msg << "stk::io::IossReader::read_model ERROR: " ;
    if ( ! ok_one_block ) {
      msg << "1 != " << node_blocks.size() << " = node_blocks.size()" ;
    }
    else {
     msg << " spatial_dim(" << m_spatial_dim
         << ") * num_nodes(" << node_ids.size()
         << ") != node_coord.size(" << node_coords.size()
         << ")" ;
    }
    throw_error( msg.str() );
  }
}

//----------------------------------------------------------------------
// Read elements, then sides, then nodes

void IossReader::read_model( mesh::BulkData & mesh ) const
{
  mesh.modification_begin();

  std::vector< mesh::Entity * > entities ;

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;

    if ( part.primary_entity_rank() == mesh::Element ) {

      Ioss::GroupingEntity & io_part =
        * m_ioss_region->get_entity( part.name() );

      read_bulk_data_block( mesh , part , io_part , entities );

      read_bulk_data_attr( entities , part , io_part );
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        ip = m_parts.begin() ; ip != m_parts.end() ; ++ip ) {

    mesh::Part & part = **ip ;

    const mesh::PartVector & subsets = part.subsets();

    if ( part.primary_entity_rank() == mesh::Face ) {

      for ( unsigned i = 0 ; i < subsets.size() ; ++i ) {
        mesh::Part & block = * subsets[i] ;

        Ioss::FaceBlock * const io_part =
          m_ioss_region->get_faceblock( block.name() );

        read_bulk_data_side( mesh , block , *io_part , entities );

        read_bulk_data_attr( entities , block , *io_part );
      }
    }
    else if ( part.primary_entity_rank() == mesh::Edge ) {

      for ( unsigned i = 0 ; i < subsets.size() ; ++i ) {
        mesh::Part & block = * subsets[i] ;

        Ioss::EdgeBlock * const io_part =
          m_ioss_region->get_edgeblock( block.name() );

        read_bulk_data_side( mesh , block , *io_part , entities );

        read_bulk_data_attr( entities , block , *io_part );
      }
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;

    if ( part.primary_entity_rank() == mesh::Node ) {

      Ioss::GroupingEntity & io_part =
        * m_ioss_region->get_entity( part.name() );

      read_bulk_data_node( mesh , part , io_part , entities );

      read_bulk_data_attr( entities , part , io_part );
    }
  }

  mesh.modification_end();

  read_bulk_data_coordinates( mesh );
}

//----------------------------------------------------------------------

IossReader::IossReader( ParallelMachine     parallel ,
                        mesh::MetaData    & meta_data ,
                        mesh::VectorField & coord_field ,
                        const std::string & file_format ,
                        const std::string & file_name ,
                        Ioss::Region      * ioss_region ,
                        const MeshReaderFilter & filter )
  : MeshReader( parallel ,
                meta_data ,
                coord_field ,
                file_format ,
                file_name ),
    m_ioss_region( ioss_region )
{
  read_meta_data_space( filter );

  {
    const Ioss::ElementBlockContainer& elem_blocks =
      m_ioss_region->get_element_blocks();

    for ( Ioss::ElementBlockContainer::const_iterator
          it = elem_blocks.begin(); it != elem_blocks.end(); ++it) {
      read_meta_data_block( filter , **it );
    }
  }

  {
    const Ioss::FaceSetContainer& face_sets =
      m_ioss_region->get_facesets();

    for ( Ioss::FaceSetContainer::const_iterator
          it = face_sets.begin(); it != face_sets.end(); ++it) {
      read_meta_data_side( filter , **it , mesh::Face );
    }
  }

  {
    const Ioss::EdgeSetContainer& edge_sets =
      m_ioss_region->get_edgesets();

    for ( Ioss::EdgeSetContainer::const_iterator
          it = edge_sets.begin(); it != edge_sets.end(); ++it) {
      read_meta_data_side( filter , **it , mesh::Edge );
    }
  }

  {
    const Ioss::NodeSetContainer& node_sets =
      m_ioss_region->get_nodesets();

    for( Ioss::NodeSetContainer::const_iterator
         it = node_sets.begin(); it != node_sets.end(); ++it) {

      read_meta_data_node( filter , **it );
    }
  }
}

IossReader::~IossReader()
{
  delete m_ioss_region ;
}

} // namespace <empty>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

MeshReader * new_reader(
  ParallelMachine          comm ,
  mesh::MetaData         & meta_data ,
  mesh::VectorField      & coord_field ,
  const MeshReaderFilter & filter ,
  const std::string      & file_format ,
  const std::string      & file_name )
{
  Ioss::Init::Initializer init_db ;

  Ioss::DatabaseIO* const iodb =
    Ioss::IOFactory::create(file_format, file_name, Ioss::READ_MODEL, comm);

  Ioss::Region * const ior =
    iodb ? new Ioss::Region( iodb , file_name ) : NULL ;

  if ( ! ior && iodb ) { delete iodb ; }

  IossReader * const reader = ! ior ? NULL :
    new IossReader( comm , meta_data , coord_field ,
                    file_format , file_name ,
                    ior , filter );

  return reader ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}//namespace io
}//namespace stk


