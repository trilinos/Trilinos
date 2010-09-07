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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_io/deprecated/MeshWriter.hpp>
#include <stk_io/IossBridge.hpp>

namespace stk {
namespace io {

MeshWriter::~MeshWriter() {}

//----------------------------------------------------------------------

typedef mesh::Field<double>                       ScalarField ;
typedef mesh::Field<double,mesh::Cartesian>       CartesianField ;
typedef mesh::Field<double,mesh::FullTensor>      FullTensorField ;
typedef mesh::Field<double,mesh::SymmetricTensor> SymmetricTensorField ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

  void get_io_field_type(const stk::mesh::FieldBase &field, int num_comp, std::pair<std::string, Ioss::Field::BasicType> *result)
  {
    const std::string invalid("invalid");
    const std::string scalar("scalar");
    const std::string vector_2d("vector_2d");
    const std::string vector_3d("vector_3d");
    const std::string quaternion_2d("quaternion_2d");
    const std::string quaternion_3d("quaternion_3d");
    const std::string full_tensor_36("full_tensor_36");
    const std::string full_tensor_32("full_tensor_32");
    const std::string full_tensor_22("full_tensor_22");
    const std::string full_tensor_16("full_tensor_16");
    const std::string full_tensor_12("full_tensor_12");
    const std::string sym_tensor_33("sym_tensor_33");
    const std::string sym_tensor_31("sym_tensor_31");
    const std::string sym_tensor_21("sym_tensor_21");
    const std::string sym_tensor_13("sym_tensor_13");
    const std::string sym_tensor_11("sym_tensor_11");
    const std::string sym_tensor_10("sym_tensor_10");
    const std::string asym_tensor_03("asym_tensor_03");
    const std::string asym_tensor_02("asym_tensor_02");
    const std::string asym_tensor_01("asym_tensor_01");
    const std::string matrix_22("matrix_22");
    const std::string matrix_33("matrix_33");

    const unsigned rank = field.rank();
    const shards::ArrayDimTag * const * const tags = field.dimension_tags();

    result->second = Ioss::Field::INVALID;

    if ( field.type_is<double>() ) {
      result->second = Ioss::Field::REAL;
    }
    else if ( field.type_is<int>() ) {
      result->second = Ioss::Field::INTEGER;
    }

    if ( 0 == rank ) {
      result->first = scalar ;
    }
    else if ( 1 == rank ) {
      if ( tags[0] == & stk::mesh::Cartesian::tag() && 2 == num_comp ) {
	result->first = vector_2d ;
      }
      else if ( tags[0] == & stk::mesh::Cartesian::tag() && 3 == num_comp ) {
	result->first = vector_3d ;
      }
      else if ( tags[0] == & stk::mesh::FullTensor::tag() && 9 == num_comp ) {
	result->first = full_tensor_36 ;
      }
      else if ( tags[0] == & stk::mesh::FullTensor::tag() && 5 == num_comp ) {
	result->first = full_tensor_32 ;
      }
      else if ( tags[0] == & stk::mesh::FullTensor::tag() && 4 == num_comp ) {
	result->first = full_tensor_22 ;
      }
      else if ( tags[0] == & stk::mesh::FullTensor::tag() && 3 == num_comp ) {
	result->first = full_tensor_12 ;
      }
      else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 6 == num_comp ) {
	result->first = sym_tensor_33 ;
      }
      else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 4 == num_comp ) {
	result->first = sym_tensor_31 ;
      }
      else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 3 == num_comp ) {
	result->first = sym_tensor_21 ;
      }
    }

    if ( result->first.empty() ) {
      std::ostringstream tmp ;
      tmp << "Real[" << num_comp << "]" ;
      result->first = tmp.str();
    }
  }

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

class IossWriter : public MeshWriter {
public:

  Ioss::Region * m_ioss_region ;
  int m_spatial_dim ;

  void write_model( mesh::BulkData & ) const ;

  void write_transients( mesh::BulkData & , double time ) const ;

  ~IossWriter();

  IossWriter(
    ParallelMachine     parallel ,
    mesh::MetaData    & meta_data ,
    mesh::VectorField & coord_field ,
    const std::vector< mesh::Part * > & parts ,
    const std::vector< std::pair<mesh::FieldBase*,std::string> > & attributes ,
    const std::vector< std::pair<mesh::FieldBase*,std::string> > & transients ,
    const std::string & file_format ,
    const std::string & file_name ,
    Ioss::Region      * ioss_region );

  //--------------------------------

  void throw_error( const std::string & ) const ;

  void write_field_info(
    Ioss::GroupingEntity  & io_part ,
    Ioss::Field::RoleType   iorole ,
    mesh::Part            & part ,
    mesh::EntityRank        type ,
    unsigned                count ) const ;

  void write_topology_info( mesh::BulkData & ) const ;
  void write_topology_bulk( mesh::BulkData & ) const ;

  void write_universal_node_block_info( mesh::BulkData & ) const ;
  void write_universal_node_block_bulk( mesh::BulkData & ) const ;

  void write_universal_node_block_transient( mesh::BulkData & ) const ;

  void write_element_block_info( mesh::BulkData & , mesh::Part & ) const;
  void write_element_block_bulk( mesh::BulkData & , mesh::Part & ) const;

  void write_node_set_info( mesh::BulkData & , mesh::Part & ) const;
  void write_node_set_bulk( mesh::BulkData & , mesh::Part & ) const;

  template< class IOSet , class IOBlock , mesh::EntityRank type >
  void write_side_block_info( mesh::BulkData    & bulk ,
                              mesh::Part        & part ,
                              IOSet             & iopart ) const ;

  void write_side_block_bulk( mesh::BulkData    & bulk ,
                              mesh::Part        & part ,
                              Ioss::GroupingEntity * const b ) const ;

private:
  IossWriter(const IossWriter& from); // do not implement
  IossWriter& operator=(const IossWriter& from); // do not implement
};

IossWriter::~IossWriter()
{
  delete m_ioss_region ;
}

void IossWriter::write_model( mesh::BulkData & bulk ) const
{
  write_topology_info( bulk );
  write_topology_bulk( bulk );
}

void IossWriter::write_transients( mesh::BulkData & bulk ,
                                   double time ) const
{
  m_ioss_region->begin_mode( Ioss::STATE_TRANSIENT );
  const int state = m_ioss_region->add_state( time );
  m_ioss_region->begin_state( state );

  write_universal_node_block_transient( bulk );

  m_ioss_region->end_state( state );
  m_ioss_region->end_mode( Ioss::STATE_TRANSIENT );
}

void IossWriter::throw_error( const std::string & msg ) const
{
  std::string ex_msg ;
  ex_msg.append( "stk::io::IossWriter[ '" );
  ex_msg.append( m_file_format );
  ex_msg.append( "' '" );
  ex_msg.append( m_file_name );
  ex_msg.append( "' ] FAILED: ");
  ex_msg.append( msg );
  throw std::runtime_error( ex_msg );
}

//----------------------------------------------------------------------

IossWriter::IossWriter(
  ParallelMachine     parallel ,
  mesh::MetaData    & meta_data ,
  mesh::VectorField & coord_field ,
  const std::vector< mesh::Part * > & parts ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & attributes ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & transients ,
  const std::string & file_format ,
  const std::string & file_name ,
  Ioss::Region      * ioss_region )
  : MeshWriter( parallel ,
                meta_data ,
                coord_field ,
                parts ,
                attributes ,
                transients ,
                file_format ,
                file_name ),
    m_ioss_region( ioss_region ),
    m_spatial_dim( 0 )
{
  const mesh::FieldBase::Restriction & res =
    m_coord_field.restriction( mesh::Node , m_meta_data.universal_part() );

  /** \todo REFACTOR  Need a clear way to query dimensions
   *                  from the field restriction.
   */
  m_spatial_dim = res.stride[0] ;
}

//----------------------------------------------------------------------

void IossWriter::write_field_info(
  Ioss::GroupingEntity  & io_part ,
  Ioss::Field::RoleType   iorole ,
  mesh::Part            & part ,
  mesh::EntityRank        type ,
  unsigned                count ) const
{
  const std::vector< std::pair<mesh::FieldBase*,std::string> > empty ;

  const std::vector< std::pair<mesh::FieldBase*,std::string> > & io_fields =
    iorole == Ioss::Field::ATTRIBUTE ? m_attributes : (
    iorole == Ioss::Field::TRANSIENT ? m_transients : empty );

  std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator i ;

  for ( i = io_fields.begin(); i != io_fields.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const std::string     & ioname = i->second ;

    const mesh::FieldBase::Restriction & res = field.restriction(type, part);

    if ( 0 < res.stride[0] ) {

      std::pair<std::string, Ioss::Field::BasicType> field_type;
      get_io_field_type( field, res.stride[0], &field_type );

      io_part.field_add(
        Ioss::Field( ioname, field_type.second, field_type.first,
                     iorole , count));
    }
  }
}


void write_field_data( const stk::mesh::FieldBase & field,
                       std::vector<stk::mesh::Entity*> & entities,
                       Ioss::GroupingEntity  & io_entity,
                       const std::string &io_fld_name)
{
  /// \todo REFACTOR Need some additional compatability checks between
  /// Ioss field and stk::mesh::Field; better error messages...

  const size_t entity_count = entities.size();
  const bool ok_io_field = io_entity.field_exists( io_fld_name );
  bool ok_field_data = true ;
  bool ok_io_write = true ;
  size_t io_entity_count = 0 ;

  if ( ok_io_field ) {

    Ioss::Field att_field = io_entity.get_field(io_fld_name);

    int field_component_count =
      att_field.transformed_storage()->component_count();

    std::vector<double> io_field(entity_count*field_component_count);

    for (size_t i=0; i < entity_count && ok_field_data ; ++i) {
      /// \todo REFACTOR Is there a way to get the data from a
      /// "FieldBase*" field as a double* without the cast?

      if (entities[i] != NULL) {

        double *fld_data = static_cast<double*>(stk::mesh::field_data(field, *entities[i]));

        ok_field_data = fld_data != NULL ;

        if ( ok_field_data ) {
          for(int j=0; j<field_component_count; ++j) {
            io_field[i*field_component_count+j] = fld_data[j];
          }
        }
      } else {
        for(int j=0; j<field_component_count; ++j) {
          io_field[i*field_component_count+j] = 0.0;
        }
      }
    }

    if ( ok_field_data ) {

      io_entity_count = io_entity.put_field_data(io_fld_name, io_field);

      ok_io_write = io_entity_count == entity_count ;
    }
  }

  if ( ! ok_io_field || ! ok_field_data || ! ok_io_write ) {
    std::ostringstream errmsg;

    errmsg << "stk::io::IossWriter::write_field_data ERROR writing from "
           << field << " to " << io_fld_name ;

    if ( ! ok_io_field ) {
      errmsg << " output field not defined" ;
    }
    else if ( ! ok_field_data ) {
      errmsg << " input field data is NULL on some mesh entities" ;
    }
    else {
      errmsg << " Field count mismatch: ioSystem has " << io_entity_count
             << " entries, but the stk:mesh system has " << entity_count
             << " entries." ;
    }
    throw std::runtime_error(errmsg.str());
  }
}

//----------------------------------------------------------------------


template< class IOSet , class IOBlock , mesh::EntityRank type >
void IossWriter::write_side_block_info(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ,
  IOSet             & iopart ) const
{
  mesh::Selector selector = ( m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part() ) & part;

  const int count = count_selected_entities(selector, bulk.buckets(type));

  const CellTopologyData *const ef_topology = mesh::get_cell_topology(part);
  if (ef_topology == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
    throw std::runtime_error( msg.str() );
  }

  const std::string io_topo = map_topology_cell_to_ioss(ef_topology, m_spatial_dim);

  IOBlock *efb =
    new IOBlock( iopart.get_database(), part.name(), io_topo, "unknown", count);

  if ( ! iopart.add( efb ) ) {
    throw_error( std::string("Ioss::Region::add FAILED") );
  }

  // Assume exactly one attribute and that attribute is the distribution factor

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;

    const mesh::FieldBase::Restriction & res = field.restriction(type, part);

    if ( 0 < res.stride[0] ) {

      int nodes_per_ef = ef_topology->node_count;

      // Requirement: node_per_ef == res.stride[0]
      // Requirement: i->second == "distribution_factors"

      std::ostringstream storage_type ;
      storage_type << "Real[" << nodes_per_ef << "]" ;
      efb->field_add( Ioss::Field( std::string("distribution_factors"),
                                   Ioss::Field::REAL, storage_type.str(),
                                   Ioss::Field::MESH, count) );
    }
  }
}

void IossWriter::write_side_block_bulk(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ,
  Ioss::GroupingEntity * const b ) const
{

  std::vector<mesh::Entity *> sides ;

  const mesh::MetaData & meta_data = part.mesh_meta_data();

  {
    mesh::Selector selector = ( m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part() ) & part;
    get_selected_entities( selector, bulk.buckets(part.primary_entity_rank()) ,
                                  sides );
  }

  const size_t num_sides = sides.size();

  std::vector<int> side_ids(num_sides);
  std::vector<int> elem_side_ids(num_sides*2);

  for(size_t i=0; i<num_sides; ++i) {

    const mesh::Entity & side = * sides[i] ;

    const mesh::PairIterRelation side_elem = side.relations( mesh::Element );

    // Which element to use?
    // Any locally owned element that has the "correct" orientation

    const size_t num_side_elem = side_elem.size();

    const mesh::Relation * rel = NULL ;

    for ( size_t j = 0 ; j < num_side_elem && ! rel ; ++j ) {
      rel = & side_elem[j] ;
      const mesh::Entity & elem = * rel->entity();

      if ( ! elem.bucket().member( meta_data.locally_owned_part() ) ||
           ! element_side_polarity( elem , side , rel->identifier() ) ) {
        rel = NULL ;
      }
    }

    if ( NULL == rel ) { // no suitable element found
      std::ostringstream oss;
      oss << "ERROR, no suitable element found";
      throw std::runtime_error(oss.str());
    }

    side_ids[i] = side.identifier();

    // Ioss is 1-based, mesh is 0-based.
    // Hence the '+1' in the following line.

    elem_side_ids[i*2]   = rel->entity()->identifier();
    elem_side_ids[i*2+1] = rel->identifier() + 1 ;
  }

  const size_t num_ids_written = b->put_field_data("ids", side_ids);
  const size_t num_side_written = b->put_field_data("element_side",elem_side_ids);

  if ( num_sides != num_ids_written || num_sides != num_side_written ) {
    std::ostringstream msg ;

    msg << "stk::io::IossWriter::write_side_block_bulk FAILED for " ;
    msg << part.name() ;
    msg << " in Ioss::GroupingEntity::put_field_data:" ;
    msg << " num_sides = " << num_sides ;
    msg << " , num_ids_written = " << num_ids_written ;
    msg << " , num_side_written = " << num_side_written ;
    throw std::runtime_error( msg.str() );
  }

  // Assume exactly one attribute and that attribute is the distribution factor

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;

    const mesh::FieldBase::Restriction & res =
      field.restriction(part.primary_entity_rank(), part);

    if ( 0 < res.stride[0] ) {
      write_field_data( field, sides, *b, std::string("distribution_factors"));
    }
  }
}

//----------------------------------------------------------------------

void IossWriter::write_topology_info( mesh::BulkData & bulk ) const
{
  /// \todo REFACTOR  The various element, node, edge, and face
  ///                 bulk-data counts are required as meta-data in the Ioss.
  ///                 However, these bulk-data counts are bulk-data

  /// \todo REFACTOR  Required some debugging to learn about
  ///                 the Ioss::Region state machine
  ///                 Ioss::STATE_DEFINE_MODEL
  ///                 Ioss::STATE_MODEL

  m_ioss_region->begin_mode( Ioss::STATE_DEFINE_MODEL );

  write_universal_node_block_info( bulk );

  /// \todo REFACTOR  There seems to be some ordering dependencies
  ///                 when defining parts.
  /// (GDS) The only restriction should be that the "ids" field for
  /// the nodes and elements must be defined before the ids are used
  /// in the output of another field.  For example, the element
  /// connectivity field uses node ids, so the node ids field must be
  /// defined to the Ioss prior to the connectivity field. There
  /// should be no ordering dependency of parts during the meta data...

  // Put IO parts into Ioss::Region

  for ( std::vector< mesh::Part* >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;
    if ( part.primary_entity_rank() == mesh::Element ) {
      write_element_block_info( bulk , part );
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;
    if ( part.primary_entity_rank() == mesh::Node ) {
      write_node_set_info( bulk , part );
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {

    mesh::Part & part = **i ;

    const stk::mesh::PartVector & subsets = part.subsets();

    if ( part.primary_entity_rank() == mesh::Face ) {

      Ioss::FaceSet * const ss =
        new Ioss::FaceSet( m_ioss_region->get_database(), part.name() );

      if ( ! m_ioss_region->add( ss ) ) {
        throw_error( std::string("Ioss::Region::add FAILED") );
      }

      if ( 0 < subsets.size() ) {
        for (size_t j = 0; j < subsets.size(); j++) {
          write_side_block_info<Ioss::FaceSet,Ioss::FaceBlock,mesh::Face>(
            bulk , *subsets[j] , *ss );
        }
      }
      else {
        write_side_block_info<Ioss::FaceSet,Ioss::FaceBlock,mesh::Face>(
          bulk , part , *ss );
      }
    }
    else if ( part.primary_entity_rank() == mesh::Edge ) {

      Ioss::EdgeSet * const ss =
        new Ioss::EdgeSet( m_ioss_region->get_database(), part.name() );

      if ( ! m_ioss_region->add( ss ) ) {
        throw_error( std::string("Ioss::Region::add FAILED") );
      }

      if ( 0 < subsets.size() ) {
        for (size_t j = 0; j < subsets.size(); j++) {
          write_side_block_info<Ioss::EdgeSet,Ioss::EdgeBlock,mesh::Edge>(
            bulk , *subsets[j] , *ss );
        }
      }
      else {
        write_side_block_info<Ioss::EdgeSet,Ioss::EdgeBlock,mesh::Edge>(
          bulk , part , *ss );
      }
    }
  }

  m_ioss_region->end_mode( Ioss::STATE_DEFINE_MODEL );

  m_ioss_region->begin_mode( Ioss::STATE_DEFINE_TRANSIENT );

  {
    mesh::Selector selector = m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part();
    const int num_nodes =
      count_selected_entities( selector, bulk.buckets(mesh::Node));

    Ioss::NodeBlock * const nb = m_ioss_region->get_node_blocks()[0];

    write_field_info( *nb, Ioss::Field::TRANSIENT ,
                      m_meta_data.universal_part(), mesh::Node, num_nodes );
  }
  m_ioss_region->end_mode( Ioss::STATE_DEFINE_TRANSIENT );
}

void IossWriter::write_topology_bulk( mesh::BulkData & bulk ) const
{
  m_ioss_region->begin_mode( Ioss::STATE_MODEL );

  write_universal_node_block_bulk( bulk );

  /// \todo REFACTOR  There seems to be some ordering dependencies
  ///                 when defining parts.
  /// (GDS) The only restriction should be that the "ids" field for
  /// the nodes and elements must be defined before the ids are used
  /// in the output of another field.  For example, the element
  /// connectivity field uses node ids, so the node ids field must be
  /// defined to the Ioss prior to the connectivity field. There
  /// should be no ordering dependency of parts during the meta data...

  // Put IO parts into Ioss::Region

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;
    if ( part.primary_entity_rank() == mesh::Element ) {
      write_element_block_bulk( bulk , part );
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {
    mesh::Part & part = **i ;
    if ( part.primary_entity_rank() == mesh::Node ) {
      write_node_set_bulk( bulk , part );
    }
  }

  for ( std::vector< mesh::Part * >::const_iterator
        i = m_parts.begin() ; i != m_parts.end() ; ++i ) {

    mesh::Part & part = **i ;

    const stk::mesh::PartVector & subsets = part.subsets();

    if ( part.primary_entity_rank() == mesh::Face ) {
      if ( 0 < subsets.size() ) {
        for ( size_t j = 0; j <  subsets.size(); j++ ) {
          mesh::Part & block = * subsets[j];
          Ioss::FaceBlock * efb = m_ioss_region->get_faceblock( block.name() );
          write_side_block_bulk( bulk, block , efb );
        }
      } else {
        Ioss::FaceBlock * efb = m_ioss_region->get_faceblock( part.name() );
        write_side_block_bulk( bulk, part, efb );
      }
    }
    else if ( part.primary_entity_rank() == mesh::Edge ) {
      if ( 0 < subsets.size() ) {
        for ( size_t j = 0; j <  subsets.size(); j++ ) {
          mesh::Part & block = * subsets[j];
          Ioss::EdgeBlock * efb = m_ioss_region->get_edgeblock( block.name() );
          write_side_block_bulk( bulk, block , efb );
        }
      } else {
        Ioss::EdgeBlock * efb = m_ioss_region->get_edgeblock( part.name() );
        write_side_block_bulk( bulk, part, efb );
      }
    }
  }

  m_ioss_region->end_mode( Ioss::STATE_MODEL );
}

//----------------------------------------------------------------------

void IossWriter::write_universal_node_block_info( mesh::BulkData & bulk ) const
{
  const std::string name("nodeblock_1");

  m_ioss_region->property_add(
    Ioss::Property( "spatial_dimension" , m_spatial_dim ) );

  mesh::Selector selector = m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part();
  const int num_nodes =
    count_selected_entities( selector, bulk.buckets(mesh::Node));

  Ioss::NodeBlock * const nb =
    new Ioss::NodeBlock( m_ioss_region->get_database(),
                         name, num_nodes, m_spatial_dim );

  if ( ! m_ioss_region->add( nb ) ) {
    throw_error( std::string("Ioss::Region::add FAILED") );
  }

  write_field_info( *nb, Ioss::Field::ATTRIBUTE ,
                    m_meta_data.universal_part(), mesh::Node, num_nodes );
}

void IossWriter::write_universal_node_block_bulk( mesh::BulkData & bulk ) const
{
  static const char method[] =
    "stk::io::IossWriter::write_universal_node_block_info" ;

  // Exactly one node block to obtain the nodal coordinates and ids:
  // Note that the "ids" field of the nodes needs to be written
  // before any other bulk data that uses node ids since it sets up
  // the global->local mapping of nodes for the output database.
  // Similarly for the element "ids" field related to bulk data
  // using element ids.

  const Ioss::NodeBlockContainer& node_blocks =
    m_ioss_region->get_node_blocks();

  if ( node_blocks.size() != 1 ) {
    std::ostringstream msg ;
    msg << method ;
    msg << ": ERROR node_blocks.size() = " << node_blocks.size();
    throw std::runtime_error( msg.str() );
  }

  Ioss::NodeBlock & nb = * node_blocks[0];

  std::vector<mesh::Entity *> nodes ;

  mesh::Selector selector = m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part();
  mesh::get_selected_entities( selector, bulk.buckets(mesh::Node), nodes );

  const size_t num_nodes = nodes.size();

  std::vector<int> node_ids(num_nodes);

  for (size_t i=0; i<num_nodes; ++i) {
    const mesh::Entity & node = * nodes[i] ;
    node_ids[i] = node.identifier();
  }

  size_t num_ids_written = nb.put_field_data("ids", node_ids);

  if ( num_nodes != num_ids_written ) {

    std::ostringstream msg ;
    msg << method << " FAILED in Ioss::NodeBlock::put_field_data:" ;
    msg << " num_nodes = " << num_nodes ;
    msg << " , num_ids_written = " << num_ids_written ;
    throw std::runtime_error( msg.str() );
  }

  write_field_data( m_coord_field, nodes, nb , "mesh_model_coordinates" );
}

void IossWriter::write_universal_node_block_transient(
  mesh::BulkData & bulk ) const
{
  bool no_nodes = true ;

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i =  m_transients.begin() ;
        i != m_transients.end() && no_nodes ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Node, m_meta_data.universal_part() );

    if ( 0 < res.stride[0] ) { no_nodes = false ; }
  }

  if ( no_nodes ) { return ; }

  Ioss::NodeBlock & nb = * m_ioss_region->get_node_blocks()[0];

  std::vector<mesh::Entity *> nodes ;

  mesh::Selector selector = m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part();
  mesh::get_selected_entities( selector, bulk.buckets(mesh::Node), nodes );

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_transients.begin(); i != m_transients.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Node, m_meta_data.universal_part() );

    if ( 0 < res.stride[0] ) {
      write_field_data( field, nodes, nb, i->second );
    }
  }
}

//----------------------------------------------------------------------

void IossWriter::write_element_block_info(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ) const
{

  const CellTopologyData * const cell_top = mesh::get_cell_topology( part );
  if (cell_top == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
    throw std::runtime_error( msg.str() );
  }

  mesh::Selector selector = m_meta_data.locally_owned_part() & part;
  const int num_elems =
    count_selected_entities( selector, bulk.buckets(mesh::Element));

  int attribute_count = 0;

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field = * i->first ;

    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Element, part);

    if ( 0 < res.stride[0] ) ++attribute_count ;
  }

  Ioss::ElementBlock *eb =
    new Ioss::ElementBlock( m_ioss_region->get_database() ,
                            part.name() ,
                            map_topology_cell_to_ioss( cell_top, m_spatial_dim ) ,
                            num_elems , attribute_count );

  if ( ! m_ioss_region->add( eb ) ) {
    throw_error( std::string("Ioss::Region::add FAILED") );
  }

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const std::string     & ioname = i->second ;

    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Element, part);

    if ( 0 < res.stride[0] ) {

      std::pair<std::string, Ioss::Field::BasicType> field_type;
      get_io_field_type( field, res.stride[0], &field_type );

        eb->field_add(Ioss::Field(ioname, field_type.second, field_type.first,
                                  Ioss::Field::ATTRIBUTE, num_elems));
    }
  }
}

void IossWriter::write_element_block_bulk(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ) const
{
  static const char method[] =
    "stk::io::IossWriter::write_element_block_bulk" ;

  const CellTopologyData* cell_topo = mesh::get_cell_topology(part);
  if (cell_topo == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
    throw std::runtime_error( msg.str() );
  }

  const int nodes_per_elem = cell_topo->node_count;

  std::vector<mesh::Entity *> elements ;

  mesh::Selector selector = part & m_meta_data.locally_owned_part();
  get_selected_entities( selector, bulk.buckets(mesh::Element), elements );

  const size_t num_elems = elements.size();

  std::cout << "elem-block = " << part.name()
            << " , num_elems = " << num_elems
            << " , nodes_per_elem = " << nodes_per_elem
            << std::endl;

  std::vector<int> elem_ids(num_elems);
  std::vector<int> connectivity(num_elems*nodes_per_elem);

  for ( size_t i = 0 ; i < num_elems ; ++i ) {

    elem_ids[i] = elements[i]->identifier();

    const mesh::PairIterRelation elem_nodes =
      elements[i]->relations( mesh::Node );

    for ( int j = 0 ; j < nodes_per_elem ; ++j ) {
      connectivity[ i * nodes_per_elem + j ] =
        elem_nodes[j].entity()->identifier();
    }
  }

  Ioss::ElementBlock * const eb = m_ioss_region->get_element_block(part.name());
  if (eb == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_element_block()";
    throw std::runtime_error( msg.str() );
  }

  const size_t num_ids_written = eb->put_field_data("ids", elem_ids);
  const size_t num_con_written = eb->put_field_data("connectivity", connectivity);

  if ( num_elems != num_ids_written || num_elems != num_con_written ) {
    std::ostringstream msg ;
    msg << method << " FAILED in Ioss::ElementBlock::put_field_data:"
        << std::endl ;
    msg << "  num_elems = " << num_elems << std::endl ;
    msg << "  num_ids_written = " << num_ids_written << std::endl ;
    msg << "  num_connectivity_written = " << num_con_written << std::endl ;
    throw std::runtime_error( msg.str() );
  }


  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const std::string     & ioname = i->second ;

    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Element, part);

    if ( 0 < res.stride[0] ) {
      write_field_data( field, elements, *eb, ioname);
    }
  }
}

//----------------------------------------------------------------------

void IossWriter::write_node_set_info(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ) const
{
  mesh::Selector selector = ( m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part() ) & part;

  const int num_nodes =
    count_selected_entities(selector, bulk.buckets(mesh::Node));

  Ioss::NodeSet * const ns =
    new Ioss::NodeSet( m_ioss_region->get_database(), part.name(), num_nodes );

  if ( ! m_ioss_region->add( ns ) ) {
    throw_error( std::string("Ioss::Region::add FAILED") );
  }

  // NOTE: A nodeset has the "distribution_factors" field added automatically.
}

void IossWriter::write_node_set_bulk(
  mesh::BulkData    & bulk ,
  mesh::Part        & part ) const
{
  std::vector<mesh::Entity *> nodes ;

  {
    mesh::Selector selector = ( m_meta_data.locally_owned_part() | m_meta_data.globally_shared_part() ) & part;

    get_selected_entities( selector, bulk.buckets(stk::mesh::Node), nodes );
  }

  const size_t num_nodes = nodes.size();

  std::vector< int > node_ids( num_nodes );

  for ( size_t i = 0 ; i < num_nodes ; ++i ) {
    node_ids[i] = nodes[i]->identifier();
  }

  Ioss::NodeSet * const ns = m_ioss_region->get_nodeset(part.name());
  if (ns == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_nodeset()";
    throw std::runtime_error( msg.str() );
  }

  const size_t num_ids_written = ns->put_field_data("ids", node_ids);

  if ( num_nodes != num_ids_written ) {
    std::ostringstream msg ;
    msg << " FAILED in Ioss::NodeSet::put_field_data:"
        << std::endl ;
    msg << "  num_nodes = " << num_nodes << std::endl ;
    msg << "  num_ids_written = " << num_ids_written << std::endl ;
    throw std::runtime_error( msg.str() );
  }

  for ( std::vector< std::pair<mesh::FieldBase*,std::string> >::const_iterator
        i = m_attributes.begin(); i != m_attributes.end() ; ++i ) {

    const mesh::FieldBase & field  = * i->first ;
    const std::string     & ioname = i->second ;

    const mesh::FieldBase::Restriction & res =
      field.restriction(mesh::Node, part);

    if ( 0 < res.stride[0] ) {
      write_field_data( field, nodes, *ns, ioname);
    }
  }
}

} // namespace <blank>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

MeshWriter * new_writer(
  ParallelMachine     comm ,
  mesh::MetaData    & meta_data ,
  mesh::VectorField & coord_field ,
  const std::vector< mesh::Part * >      & parts ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & attributes ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & transients ,
  const std::string & file_format ,
  const std::string & file_name )
{
  Ioss::Init::Initializer init_db;

  Ioss::DatabaseIO* const iodb =
    Ioss::IOFactory::create(file_format, file_name, Ioss::WRITE_RESULTS, comm);

  Ioss::Region * const ior =
    iodb ? new Ioss::Region( iodb , file_name ) : NULL ;

  if ( ! ior && iodb ) { delete iodb ; }

  IossWriter * writer = ! ior ? NULL :
    new IossWriter( comm , meta_data , coord_field ,
                parts , attributes , transients ,
                file_format , file_name , ior );

  return writer ;
}



} // namespace io
} // namespace stk

