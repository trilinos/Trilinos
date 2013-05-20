/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.               */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_MetaData_hpp
#define stk_mesh_MetaData_hpp

//----------------------------------------------------------------------

#include <iosfwd>

#include <stk_util/util/SameType.hpp>
#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/PropertyBase.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/CellTopology.hpp>
#include <stk_mesh/base/Trace.hpp>

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>
#include <stk_mesh/baseImpl/FieldRepository.hpp>

#include <stk_topology/topology.hpp>

namespace shards {
  class CellTopologyManagedData;
}

namespace stk {
namespace mesh {

class BulkData;

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print an entity key for this meta data */
std::ostream &
print_entity_key( std::ostream & os, const MetaData & meta_data, const EntityKey & key);

std::string
print_entity_key( const MetaData & meta_data, const EntityKey & key );

//----------------------------------------------------------------------
/** \brief  The manager of an integrated collection of
 *          \ref stk::mesh::Part "parts" and
 *          \ref stk::mesh::Field "fields".
 *
 *  Mesh meta data must be identical on all processors.
 *
 * The FEM features include the concept of spatial dimension with
 * entity ranks tied to the given spatial dimension, cell topology
 * mapping to parts along with induced cell topology membership
 * through part subsetting, and many additional invariants that are
 * enforced.
 *
 * Invariants for MetaData:
 * 1.  Each cell topology has one and only one root cell topology part.  The
 *     root cell topology part for a cell topology is a unique part that can be
 *     subsetted to induce cell topology on the subset part.
 *     -> Enforced by register_cell_topology is the only function that modifies
 *     the PartCellTopologyVector private data on MetaData.
 * 2.  Root cell topology parts cannot be subsets of parts with cell topologies
 *     -> Enforced by declare_part_subset
 * 3.  Incompatible cell topologies are prohibited.  I.e. parts with
 *     different cell topologies of the same rank (as the part) cannot be subsets
 *     of each other.
 *     -> Enforced by declare_part_subset
 */

// 02/10/11 MetaData Todo:
// * Implement get_cell_topology for Part.
// * Implement declare_part with cell topology
// Non-critical:
// * Implement stk::mesh::get namespace to include getters for MetaData,
//   BulkData, MetaData, BulkData, CellTopology from things like Part,
//   Bucket, Entity, etc.
// * Create impl class inside the handle classes to hold their parent pointer
//   and a friend to the getter above.
class MetaData {
public:

  /** \} */
  //------------------------------------
  /** \name  Meta data manager construction and destruction
   *  \{
   */

  /// CellTopologyPartEntityRankMap maps each Cell Topology to its root cell topology part and its associated rank
  typedef std::map<CellTopology, std::pair<Part *, EntityRank> > CellTopologyPartEntityRankMap;
  /// PartCellTopologyVector is a fast-lookup vector of size equal to the number of parts
  typedef std::vector<CellTopology> PartCellTopologyVector;

  enum EntityRankValue
  {
    NODE_RANK = stk::topology::NODE_RANK,
    EDGE_RANK = stk::topology::EDGE_RANK,
    FACE_RANK = stk::topology::FACE_RANK,
    ELEMENT_RANK = stk::topology::ELEMENT_RANK,
    CONSTRAINT_RANK = stk::topology::CONSTRAINT_RANK,
    INVALID_RANK = stk::topology::INVALID_RANK
  };

  inline static MetaData & get( const Part & part ) { return part.meta_data(); }
  inline static MetaData & get( const FieldBase & field ) { return field.meta_data(); }
  inline static MetaData & get( const PropertyBase & property ) { return property.meta_data(); }

  static MetaData & get( const BulkData & bulk_data );
  static MetaData & get( const Bucket & bucket );
  static MetaData & get( const Ghosting & ghost );

  /** \brief  Construct a meta data manager to own parts and fields.  */
  explicit MetaData(size_t spatial_dimension, const std::vector<std::string>& rank_names = std::vector<std::string>());

  /** \brief  Construct a meta data manager to own parts and fields.  */
  MetaData();

  /** \brief  Destroy the meta data manager and
   *          all of the parts and fields that it owns.
   */
  ~MetaData();

  //------------------------------------
  /** \name Predefined Parts
   *  \{
   */

  /** \brief  Universal subset for the problem domain.
   *          All other parts are a subset of the universal part.
   */
  Part & universal_part() const { return *m_universal_part; }

  /** \brief  Subset for the problem domain that is owned by the
   *          local process.  Ghost entities are not members of this part.
   */
  Part & locally_owned_part()  const { return *m_owns_part ; }

  /** \brief  Subset for the problem domain that is shared with another
   *          process.  Ghost entities are not members of this part.
   */
  Part & globally_shared_part() const { return *m_shares_part ; }

  /** \} */
  //------------------------------------
  /** \name  Declare and query parts
   *  \{
   */

  /** \brief  Get an existing part by its application-defined text name.
   *
   *  Return NULL if not present and required_by == NULL.
   *  If required and not present then throws an exception
   *  with the 'required_by' text.
   */
  /// \todo REFACTOR remove required_by argument
  Part * get_part( const std::string & p_name,
                   const char * required_by = NULL ) const ;

  /** \brief  Get an existing part by its ordinal */
  Part & get_part( unsigned ord ) const ;

  /** \brief  Query all parts of the mesh ordered by the parts' ordinal. */
  const PartVector & get_parts() const { return m_part_repo.get_all_parts(); }

  /** \brief  Query non-internal parts of the mesh ordered by the parts' ordinal. */
  const PartVector get_mesh_parts() const { return m_part_repo.get_mesh_parts(); }

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part will have member entities that are of equal or lesser rank.
   *  When an entity of equal rank becomes a member
   *  then all related entities of lesser rank also become members.
   */
  Part & declare_part( const std::string & p_name, EntityRank rank );

  /** \brief  Declare a part of the given name
   *          Redeclaration returns the previously declared part.
   *
   *  This part does not have an entity type rank.
   */
  Part & declare_part( const std::string & p_name);


  /** \brief  Declare a part with a given cell topology
   */
  Part &declare_part( const std::string &name, CellTopology cell_topology)
  {
    ThrowRequireMsg(is_initialized(),"MetaData::declare_part: initialize() must be called before this function");
    Part &root_part = get_cell_topology_root_part(cell_topology);
    EntityRank primary_entity_rank = root_part.primary_entity_rank();
    Part & part = declare_part(name, primary_entity_rank);
    declare_part_subset(root_part, part);
    return part;
  }

  /** \brief  Declare a part with a given cell topology
   */
  template< class Top >
  Part &declare_part(const std::string &name) {
    return declare_part(name, shards::getCellTopologyData<Top>());
  }

  /** \brief  Declare a superset-subset relationship between parts */
  void declare_part_subset( Part & superset , Part & subset );

    /** \brief  Declare an attribute on a part.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   * \todo REFACTOR  Should be using a shared pointer in this interface.
   *       declare_attribute( Part & , shared_ptr<const T> & );
   */
  template<class T>
  const T * declare_attribute_with_delete( Part & part, const T * attribute);
  template<class T>
  const T * declare_attribute_no_delete( Part & part, const T * attribute);
  template<class T>
  bool remove_attribute( Part & part, const T * attribute);

  /** \} */
  //------------------------------------
  /** \name  Entity-ranks
   *  \{
   */

  /** \brief initialize
   *
   */
  void initialize(size_t spatial_dimension, const std::vector<std::string> &rank_names = std::vector<std::string>());

  bool is_initialized() const
  { return !m_entity_rank_names.empty(); }

  EntityRank entity_rank( const std::string &name ) const;

  const std::vector<std::string> & entity_rank_names() const
    { return m_entity_rank_names ; }

  std::vector<std::string>::size_type entity_rank_count() const
    { return m_entity_rank_names.size(); }

  const std::string & entity_rank_name( EntityRank entity_rank ) const ;

  /** \brief Returns the side rank which changes depending on spatial dimension
   */
  EntityRank side_rank() const
  {
    switch (m_spatial_dimension)
    {
    case 1 : return stk::topology::NODE_RANK;
    case 2 : return stk::topology::EDGE_RANK;
    case 3 : return stk::topology::FACE_RANK;
    default: break;
    }
    return stk::topology::INVALID_RANK;
  }

  /**
   * Return true if rank is valid.
   */
  bool check_rank(EntityRank rank) const;

  /** \} */
  //------------------------------------
  /** \name  Declare and query fields
   *  \{
   */

  /** \brief  Get a field, return NULL if it does not exist.
   *
   *  \exception std::runtime_error
   *    If the field exits and the
   *    \ref stk::mesh::Field "field_type" does not match or
   *    if required_by != NULL and a field of that name is not found.
   */
  template< class field_type >
  field_type * get_field( const std::string & name ) const ;

  /**
   * \brief Get a field by name with unknown type, NULL if does not exist
   */
  FieldBase* get_field( const std::string& name ) const;

  /** \brief  Get all defined fields */
  const FieldVector & get_fields() const {
    return m_field_repo.get_fields() ;
  }

  /** \brief  Declare a field of the given
   *          \ref stk::mesh::Field "field_type", test name,
   *          and number of states.
   *
   *  A compatible redeclaration returns the previously declared field.
   *  \exception std::runtime_error  If a redeclaration is incompatible
   *
   *  See Field.hpp for a full discussion of Fields.
   */
  template< class field_type >
  field_type & declare_field( const std::string & name ,
                              unsigned number_of_states = 1 );

  /** \brief  Declare an attribute on a field.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   * \todo REFACTOR  Should be using a shared pointer in this interface.
   */
  template<class T>
  const T * declare_attribute_with_delete( FieldBase & field, const T * attribute);
  template<class T>
  const T * declare_attribute_no_delete( FieldBase & field, const T * attribute);

  /** \} */
  //------------------------------------

  template<class T>
  const T * get_attribute() const ;

  /** \brief  Declare an attribute on the meta data.
   *          Return the attribute of that type,
   *          which may be an already existing value.
   */
  template<class T>
  const T * declare_attribute_with_delete( const T * attribute);

  template<class T>
  const T * declare_attribute_no_delete( const T * attribute);

  template<class T>
  bool remove_attribute( const T * );

  //------------------------------------
  /** \} */
  /** \name  Declare and query properties associated with parts
   *  \{
   */

  /** \brief  Get a property, return NULL if it does not exist.
   *
   *  \exception std::runtime_error
   *    If the property exits and the
   *    \ref stk::mesh::Property "type" does not match or
   */
  template< typename DataType >
  Property<DataType> * get_property( const std::string & name ) const ;

  /** \brief  Get all defined properties */
  const std::vector< PropertyBase * > & get_properties() const
    { return m_properties ; }

  /** \brief  Declare a property of the given
   *          \ref stk::mesh::Property "type", name, and dimensions.
   *
   *  A compatible redeclaration returns the previously declared property.
   *  \exception std::runtime_error  If a redeclaration is incompatible
   */
  template< typename DataType >
  Property<DataType> & declare_property( const std::string & name ,
                                         unsigned size = 1 );

  /** \brief  Put a property on the given part */
  void put_property( PropertyBase & property, Part & part);

  /** get the spatial-dimension. */
  unsigned spatial_dimension() const { return m_spatial_dimension; }

  /** \brief  Commit the part and field declarations so that the
   *          meta data manager can be used to create
   *          \ref stk::mesh::BulkData "mesh bulk data".
   *
   *  Verifies consistency of the meta data and clean out redundant
   *  field data allocation rules.
   *  Once committed no further part or field declarations can be made.
   */
  void commit();

  /** \brief  Query if the meta data manager is committed */
  bool is_commit() const { return m_commit ; }

  /** \} */
  //------------------------------------

  /** \name  Field declaration with weak type information;
   *         direct use in application code is strongly discouraged.
   *  \{
   */

  /** \brief  Declare a field via runtime type information */
  FieldBase * declare_field_base(
    const std::string & arg_name,
    const DataTraits  & arg_traits ,
    unsigned            arg_rank ,
    const shards::ArrayDimTag * const * arg_dim_tags ,
    unsigned arg_num_states )
  {
    require_not_committed();

    return m_field_repo.declare_field(
                  arg_name, arg_traits, arg_rank, arg_dim_tags,
                  arg_num_states, this
                 );
  }

  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  EntityRank       arg_entity_rank ,
                                  const Part     & arg_part ,
                                  const unsigned * arg_stride ,
                                  const void*      arg_init_value = NULL );

  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  EntityRank       arg_entity_rank ,
                                  const Selector & arg_selector ,
                                  const unsigned * arg_stride ,
                                  const void*      arg_init_value = NULL );

  /** \brief This function is used to register new cell topologies and their associated ranks with MetaData.
   * Currently, several shards Cell Topologies are registered with appropriate ranks at initialization time.
   * See:  internal_declare_known_cell_topology_parts for the whole list.
   *
   * Note:  This function also creates the root cell topology part which is accessible from get_cell_topology_root_part
   */
  void register_cell_topology(const CellTopology cell_topology, EntityRank in_entity_rank);

  shards::CellTopology register_superelement_cell_topology(stk::topology t);

  /** \brief Return the root cell topology part associated with the given cell topology.
   * This Part is created in register_cell_topology
   */

  Part &get_cell_topology_root_part(const CellTopology cell_topology) const;

  /** \brief Return the cell topology associated with the given part.
   * The cell topology is set on a part through part subsetting with the root
   * cell topology part.
   */
  CellTopology get_cell_topology( const Part & part) const;

  CellTopology get_cell_topology( const std::string & topology_name) const;

  /** \brief Return the EntityRank that is associated with the given cell
   * topology.  In several cases, this rank is dependent on spatial
   * dimension.
   */
  EntityRank get_entity_rank(const CellTopology cell_topology) const;

  void dump_all_meta_info(std::ostream& out = std::cout) const;

  /** \} */
private:
  // Functions

  MetaData( const MetaData & );                ///< \brief  Not allowed
  MetaData & operator = ( const MetaData & );  ///< \brief  Not allowed

  Part & declare_internal_part( const std::string & p_name);

  Part & declare_internal_part( const std::string & p_name, EntityRank rank);

  void internal_declare_known_cell_topology_parts();

  void internal_declare_part_subset( Part & superset , Part & subset );

  void assign_cell_topology( Part & part, CellTopology topo);

  // Members

  bool   m_commit ;
  impl::PartRepository m_part_repo ;
  CSet   m_attributes ;

  Part * m_universal_part ;
  Part * m_owns_part ;
  Part * m_shares_part ;

  impl::FieldRepository        m_field_repo ;

  std::vector< PropertyBase* > m_properties ;
  std::vector< std::string >   m_entity_rank_names ;
  std::vector<shards::CellTopologyManagedData*> m_created_topologies;

  unsigned m_spatial_dimension;
  EntityRank m_side_rank;

  /// Used to store mapping between Cell Topologies and their associated root parts and specified ranks:
  CellTopologyPartEntityRankMap m_cellTopologyPartEntityRankMap;
  /// Fast-lookup vector that maps part ordinals to Cell Topologies.
  PartCellTopologyVector        m_partCellTopologyVector;

  /** \name  Invariants/preconditions for MetaData.
   * \{
   */
  void require_committed() const ;

  void require_not_committed() const ;

  void require_same_mesh_meta_data( const MetaData & rhs ) const ;

  void require_valid_entity_rank( EntityRank rank) const ;

  void require_not_relation_target( const Part * const part ) const ;
  /** \} */
  //------------------------------------

  Property<void> * get_property_base( const std::string & ,
                                      const std::type_info & ,
                                      unsigned = 0 ) const ;

  void clean_field_restrictions();
};

/** \brief  Verify that the meta data is identical on all processors */
void verify_parallel_consistency( const MetaData & , ParallelMachine );

bool is_cell_topology_root_part(const Part & part);

/** set a cell_topology on a part */
void set_cell_topology( Part &part, const CellTopology cell_topology);

/** set a cell_topology on a part */
template<class Topology>
inline void set_cell_topology(Part & part)
{
  stk::mesh::set_cell_topology(part, CellTopology(shards::getCellTopologyData<Topology>()));
}


/** Get the cell_topology off a bucket */
CellTopology get_cell_topology(const Bucket &bucket);

CellTopology get_cell_topology(Entity entity);

/** set a stk::topology on a part */
void set_topology(Part &part, stk::topology topology);

/** Get the cell_topology off an entity */
CellTopology get_cell_topology(const Entity entity);

/** get the stk::topology given a Shards Cell Topology */
stk::topology get_topology(CellTopology shards_topology, int spatial_dimension = 3);

/** Get the Shards Cell Topology given a stk::topology  */
CellTopology get_cell_topology(stk::topology topo);

/** Get default entity rank names */
const std::vector<std::string>& entity_rank_names();

/** \name  Declare field data allocation rules
 *  \{
 */

/** \brief  Declare a field to exist for a given entity type and Part.
 *
 * See Field.hpp for a full discussion of field restrictions.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Selector & selector ,
                        const void* init_value = NULL);

/** \brief Declare a field to exist for a given entity type and Part. The
 *         extra unsigned arguments specify the size of a dimension. So,
 *         put_field( field, rank, part, 3, 3 ) would create a 3x3 2D field.
 *         Fields of up to seven dimensions are supported.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Selector & selector ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 ,
                        const void* init_value = NULL);

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 ,
                        unsigned     n7 ,
                        const void* init_value = NULL);
/** \} */
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

namespace stk {
namespace mesh {

inline
Part & MetaData::get_part( unsigned ord ) const
{ return * m_part_repo.get_all_parts()[ord] ; }

template< class field_type >
inline
field_type * MetaData::get_field( const std::string & name ) const
{
  typedef FieldTraits< field_type > Traits ;

  const DataTraits & dt = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * tags[8] ;

  Traits::assign_tags( tags );

  FieldBase * const field =
    m_field_repo.get_field( "stk::mesh::MetaData::get_field" ,
                          name , dt , Traits::Rank , tags , 0 );
  if (field == NULL) {
    return static_cast<field_type*>(NULL);
  }
  else {
    return dynamic_cast< field_type * >( field );
  }
}

template< class field_type >
inline
field_type & MetaData::declare_field( const std::string & name ,
                                      unsigned number_of_states )
{
  typedef FieldTraits< field_type > Traits ;

  const DataTraits & traits = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * dim_tags[8] ;

  Traits::assign_tags( dim_tags );

  static const char* reserved_state_suffix[6] = {
    "_STKFS_OLD",
    "_STKFS_N",
    "_STKFS_NM1",
    "_STKFS_NM2",
    "_STKFS_NM3",
    "_STKFS_NM4"
  };

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len_name   = name.size();
    const int len_suffix = std::strlen( reserved_state_suffix[i] );
    const int offset     = len_name - len_suffix ;
    if ( 0 <= offset ) {
      const char * const name_suffix = name.c_str() + offset ;
      ThrowErrorMsgIf( equal_case( name_suffix , reserved_state_suffix[i] ),
          "For name = \"" << name_suffix <<
          "\" CANNOT HAVE THE RESERVED STATE SUFFIX \"" <<
          reserved_state_suffix[i] << "\"" );
    }
  }

  // Check that the field of this name has not already been declared

  field_type * f[ MaximumFieldStates ] ;

  f[0] = dynamic_cast<field_type*>(m_field_repo.get_field(
      "MetaData::declare_field" ,
      name ,
      traits ,
      Traits::Rank ,
      dim_tags ,
      number_of_states
      ));

  if ( NULL != f[0] ) {
    for ( unsigned i = 1 ; i < number_of_states ; ++i ) {
      f[i] = &f[0]->field_of_state(static_cast<FieldState>(i));
    }
  }
  else {
    // Field does not exist then create it

    std::string field_names[ MaximumFieldStates ];

    field_names[0] = name ;

    if ( 2 == number_of_states ) {
      field_names[1] = name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < number_of_states ; ++i ) {
        field_names[i] = name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < number_of_states ; ++i ) {

      f[i] = new field_type(
          this,
          m_field_repo.get_fields().size() ,
          field_names[i] ,
          traits ,
          Traits::Rank,
          dim_tags,
          number_of_states ,
          static_cast<FieldState>(i)
          );

      m_field_repo.add_field( f[i] );
    }

    for ( unsigned i = 0 ; i < number_of_states ; ++i ) {
      f[i]->m_impl.set_field_states( f );
    }
  }

  return *f[0] ;
}

template< class field_type >
inline
field_type & put_field(
  field_type & field ,
  EntityRank entity_rank ,
  const Part & part ,
  const void* init_value)
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field(
  field_type & field ,
  EntityRank entity_rank ,
  const Selector & selector ,
  const void* init_value)
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride );

  MetaData::get(field).declare_field_restriction( field, entity_rank, selector, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, selector, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, selector, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Selector &selector ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, selector, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        unsigned    n6 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5, n6 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 ,
                        unsigned    n4 ,
                        unsigned    n5 ,
                        unsigned    n6 ,
                        unsigned    n7 ,
                        const void* init_value )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5, n6, n7 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride, init_value);

  return field ;
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( const T * a )
{
  require_not_committed();
  return m_attributes.insert_with_delete( a );
}

template<class T>
inline
const T *
MetaData::get_attribute() const
{ return m_attributes.get<T>(); }

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( const T * attribute )
{
  require_not_committed();
  return m_attributes.insert_no_delete( attribute );
}

template<class T>
inline
bool
MetaData::remove_attribute( const T * a )
{
  return m_attributes.remove( a );
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( Part & part , const T * attribute )
{
  require_not_committed();
  return m_part_repo.declare_attribute_with_delete( part, attribute );
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( Part & part , const T * attribute )
{
  require_not_committed();
  return m_part_repo.declare_attribute_no_delete( part, attribute );
}

template<class T>
inline
bool
MetaData::remove_attribute( Part & part , const T * attribute )
{
  return m_part_repo.remove_attribute(part, attribute);
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( FieldBase & field , const T * attribute )
{
  require_not_committed();
  return m_field_repo.declare_attribute_with_delete(field, attribute);
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( FieldBase & field , const T * attribute )
{
  require_not_committed();
  return m_field_repo.declare_attribute_no_delete(field, attribute);
}

//----------------------------------------------------------------------


template< typename DataType >
inline
Property<DataType> *
MetaData::get_property( const std::string & name ) const
{
  Property<void> * const pv = get_property_base( name, typeid(DataType) );
  return pv ? pv->property<DataType>() : (Property<DataType>*) NULL ;
}

template< typename DataType >
inline
Property<DataType> &
MetaData::declare_property( const std::string & name , unsigned size )
{
  Property<void> * pv = get_property_base(name,typeid(DataType),size);
  Property<DataType> * prop = NULL ;

  if ( pv != NULL ) {
    prop = pv->property<DataType>();
  }
  else {
    if ( 1 == size ) {
      pv = prop = new Property<DataType>( *this , m_properties.size() , name );
    }
    else {
      pv = prop = new Property< std::vector<DataType> >(
                    *this , m_properties.size() , name , size );
    }
    m_properties.push_back( pv );
  }
  return *prop ;
}

inline
void MetaData::put_property( PropertyBase & property , Part & part )
{
  property.add_property( part.mesh_meta_data_ordinal() );
}

inline
bool MetaData::check_rank(EntityRank rank) const
{
  return rank < m_entity_rank_names.size();
}

inline
bool
is_auto_declared_part(const Part &part)
{
  const std::string &part_name = part.name();

  return !part_name.empty() && part_name[0] == '{';
}

} // namespace mesh
} // namespace stk

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_MetaData_hpp */
