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

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/PropertyBase.hpp>
#include <stk_mesh/base/EntityKey.hpp>

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>
#include <stk_mesh/baseImpl/FieldRepository.hpp>

namespace stk {
namespace mesh {


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
 */
class MetaData {
public:

  /** \} */
  //------------------------------------
  /** \name  Meta data manager construction and destruction
   *  \{
   */

  inline static MetaData & get( const Part & part ) { return part.meta_data(); }
  inline static MetaData & get( const FieldBase & field ) { return field.meta_data(); }
  inline static MetaData & get( const PropertyBase & property ) { return property.meta_data(); }

  static MetaData & get( const BulkData & bulk_data );
  static MetaData & get( const Bucket & bucket );
  static MetaData & get( const Entity & entity );
  static MetaData & get( const Ghosting & ghost );

  /** \brief  Construct a meta data manager to own parts and fields.  */
  explicit MetaData( const std::vector<std::string>& entity_rank_names );

  /** \brief  Construct a meta data manager to own parts and fields.  */
  MetaData();

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

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part will have member entities that are of equal or lesser rank.
   *  When an entity of equal rank becomes a member
   *  then all related entities of lesser rank also become members.
   */
  Part & declare_part( const std::string & p_name, EntityRank rank );

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part does not have an entity type rank.
   */
  Part & declare_part( const std::string & p_name);

  /** \brief  Declare a part as the defined-intersection
   *          of the given collection of parts.
   *
   *  The entity type rank will be the smallest entity type rank
   *  of the parts in the intersection, if any member has an
   *  entity type rank.
   */
  Part & declare_part( const PartVector & p_name);

  /** \brief  Declare a superset-subset relationship between parts */
  void declare_part_subset( Part & superset , Part & subset );

  /** \brief  Declare an entity-relationship between parts.
   *
   *  If \ref stk::mesh::Entity "entity" <b> e1 </b> is a member
   *  of <em> root_part </em> and there exists an
   *  \ref stk::mesh::Relation "entity relation"
   *  from <b> e1 </b> to <b> e2 </b> that satisfies the
   *  \ref stk_mesh_relations "relation stencil"
   *  then <b> e2 </b> must be a member of the <em> target_part </em>.
   */
  void declare_part_relation( Part & root_part ,
                              relation_stencil_ptr stencil ,
                              Part & target_part );

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

  /** \} */
  //------------------------------------
  /** \name  Entity-ranks
   *  \{
   */
  /** \brief entity-rank names
   *
   */
  void set_entity_rank_names(const std::vector<std::string> &entity_rank_names);

  EntityRank entity_rank( const std::string &name ) const;

  const std::vector<std::string> & entity_rank_names() const
    { return m_entity_rank_names ; }

  std::vector<std::string>::size_type entity_rank_count() const
    { return m_entity_rank_names.size(); }

  const std::string & entity_rank_name( EntityRank entity_rank ) const ;

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

  /** \brief Declare a field relation.
   *
   *  The pointer_field's scalar type must be a pointer to the
   *  scalar type of the reference_field.  The following
   *  derived field data relationship maintained.
   *
   *  Let   e_root -> Relation( e_target , ord , kind )
   *  Let   i = stencil( e_root.entity_rank() ,
   *                     e_target.entity_rank() , ord , kind )
   *  Let   Scalar ** ptr = field_data( pointer_field , e_root )
   *  then  ptr[i] = field_data( referenced_field , e_target )
   *
   *  This derived field data relationship is typically used
   *  to support fast access to field data on entities
   *  related to the root entity; e.g. field data associated with
   *  the nodes of an element.
   *
   *  See FieldRelation.hpp for a full discussion of field relations.
   */
  template< class PointerFieldType , class ReferencedFieldType >
  void declare_field_relation( PointerFieldType & pointer_field ,
                               relation_stencil_ptr stencil ,
                               ReferencedFieldType & referenced_field );

  /** \brief  Get all field relations */
  const std::vector<FieldRelation> & get_field_relations() const
    { return m_field_relations ; }

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

  /** \brief  Destroy the meta data manager and
   *          all of the parts and fields that it owns.
   */
  ~MetaData();

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
    unsigned arg_num_states );

  /** \brief  Declare a field restriction via runtime type information.
   */
  void declare_field_restriction( FieldBase      & arg_field ,
                                  EntityRank       arg_entity_rank ,
                                  const Part     & arg_part ,
                                  const unsigned * arg_stride );
  /** \} */
private:
  MetaData( const MetaData & );                ///< \brief  Not allowed
  MetaData & operator = ( const MetaData & );  ///< \brief  Not allowed

  bool   m_commit ;
  impl::PartRepository m_part_repo ;
  CSet   m_attributes ;

  Part * m_universal_part ;
  Part * m_owns_part ;
  Part * m_shares_part ;


  impl::FieldRepository        m_field_repo ;

  std::vector< FieldRelation > m_field_relations ;
  std::vector< PropertyBase* > m_properties ;
  std::vector< std::string >   m_entity_rank_names ;

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

  void internal_declare_field_relation( FieldBase & ,
                                        relation_stencil_ptr ,
                                        FieldBase & );

  void clean_field_restrictions();
};

/** \brief  Verify that the meta data is identical on all processors */
void verify_parallel_consistency( const MetaData & , ParallelMachine );

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
                        const Part & part );

/** \brief Declare a field to exist for a given entity type and Part. The
 *         extra unsigned arguments specify the size of a dimension. So,
 *         put_field( field, rank, part, 3, 3 ) would create a 3x3 2D field.
 *         Fields of up to seven dimensions are supported.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 );

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 );

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 );

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 );

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 );

template< class field_type >
field_type & put_field( field_type & field ,
                        EntityRank  entity_rank ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 );

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
                        unsigned     n7 );
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

  return static_cast< field_type * >( field );
}

template< class field_type >
inline
field_type & MetaData::declare_field( const std::string & name ,
                                      unsigned number_of_states )
{
  typedef FieldTraits< field_type > Traits ;

  const DataTraits & dt = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * tags[8] ;

  Traits::assign_tags( tags );

  return * static_cast< field_type * >(
    declare_field_base( name , dt , Traits::Rank , tags , number_of_states ) );
}

template< class field_type >
inline
field_type & put_field(
  field_type & field ,
  EntityRank entity_rank ,
  const Part & part )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        EntityRank entity_rank ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

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
                        unsigned    n4 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

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
                        unsigned    n5 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

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
                        unsigned    n6 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5, n6 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

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
                        unsigned    n7 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 , n4, n5, n6, n7 );

  MetaData::get(field).declare_field_restriction( field, entity_rank, part, stride);

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

template< class PointerFieldType , class ReferencedFieldType >
inline
void MetaData::declare_field_relation(
  PointerFieldType & pointer_field ,
  relation_stencil_ptr stencil ,
  ReferencedFieldType & referenced_field )
{
  typedef typename FieldTraits< PointerFieldType >::data_type pointer_type ;
  typedef typename FieldTraits< ReferencedFieldType >::data_type data_type ;

  StaticAssert< SameType< pointer_type , data_type * >::value >::ok();
  StaticAssert< FieldTraits< PointerFieldType >::Rank == 1 >::ok();

  internal_declare_field_relation( pointer_field , stencil , referenced_field );
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
