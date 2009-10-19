

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
#include <stk_mesh/base/Property.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print an entity key for this meta data */
std::ostream &
print_entity_key( std::ostream & , const MetaData & , const EntityKey & );

//----------------------------------------------------------------------
/** \brief  The manager of an integrated collection of
 *          \ref stk::mesh::Part "parts" and
 *          \ref stk::mesh::Field "fields".
 *
 *  Mesh meta data must be identical on all processors.
 */
class MetaData {
public:

  //------------------------------------
  /** \name Predefined Parts
   *  \{
   */

  /** \brief  Universal subset for the problem domain.
   *          All other parts are a subset of the universal part.
   */
  Part & universal_part() const { return *m_universal_part; }

  /** \brief  Subset for the problem domain that is used by the
   *          local processor.  Ghost entities are not members of this part.
   */
  Part & locally_used_part() const { return *m_uses_part ; }

  /** \brief  Subset for the problem domain that is owned by the
   *          local processor.  A subset of the locally_used_part.
   */
  Part & locally_owned_part()  const { return *m_owns_part ; }

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
  Part * get_part( const std::string & ,
                   const char * required_by = NULL ) const ;

  /** \brief  Get an existing part by its ordinal */
  Part & get_part( unsigned ) const ;

  /** \brief  Query all parts of the mesh ordered by the parts' ordinal. */
  const PartVector & get_parts() const { return m_universal_part->subsets(); }

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part will have member entities that are of equal or lesser rank.
   *  When an entity of equal rank becomes a member
   *  then all related entities of lesser rank also become members.
   */
  Part & declare_part( const std::string & , EntityType rank );

  /** \brief  Declare a part of the given name and entity rank
   *          Redeclaration returns the previously declared part.
   *
   *  This part does not have an entity type rank.
   */
  Part & declare_part( const std::string & );

  /** \brief  Declare a part as the defined-intersection
   *          of the given collection of parts.
   *
   *  The entity type rank will be the smallest entity type rank
   *  of the parts in the intersection, if any member has an
   *  entity type rank.
   */
  Part & declare_part( const PartVector & );

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
  const T * declare_attribute_with_delete( Part & , const T *);
  template<class T>
  const T * declare_attribute_no_delete( Part & , const T *);

  /** \} */
  //------------------------------------
  /** \name  Entity-types
   *  \{
   */
  /** \brief Query entity-type names
   *
   */
  const std::vector<std::string> & entity_type_names() const
    { return m_entity_type_names ; }

  std::vector<std::string>::size_type entity_type_count() const
    { return m_entity_type_names.size(); }

  const std::string & entity_type_name( unsigned ) const ;

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
  const std::vector< FieldBase * > & get_fields() const { return m_fields ; }

  /** \brief  Declare a field of the given
   *          \ref stk::mesh::Field "field_type", test name,
   *          and number of states.
   *
   *  A compatible redeclaration returns the previously declared field.
   *  \exception std::runtime_error  If a redeclaration is incompatible
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
  const T * declare_attribute_with_delete( FieldBase & , const T *);
  template<class T>
  const T * declare_attribute_no_delete( FieldBase & , const T *);

  /** \brief Declare a field relation.
   *
   *  The pointer_field's scalar type must be a pointer to the
   *  scalar type of the reference_field.  The following
   *  derived field data relationship maintained.
   *
   *  Let   e_root -> Relation( e_target , ord , kind )
   *  Let   i = stencil( e_root.entity_type() ,
   *                     e_target.entity_type() , ord , kind )
   *  Let   Scalar ** ptr = field_data( pointer_field , e_root )
   *  then  ptr[i] = field_data( referenced_field , e_target )
   *
   *  This derived field data relationship is typically used
   *  to support fast access to field data on entities
   *  related to the root entity; e.g. field data associated with
   *  the nodes of an element.
   */
  template< class PointerFieldType , class ReferencedFieldType >
  void declare_field_relation( PointerFieldType & pointer_field ,
                               relation_stencil_ptr stencil ,
                               ReferencedFieldType & referenced_field );

  /** \todo REFACTOR eliminate this method, it does not belong here. */
  void declare_field_lock_relation( FieldBase & pointer_field ,
                                    relation_stencil_ptr stencil );
  
  /** \brief  Get field relations */
  const std::vector<FieldRelation> & get_field_relations() const
    { return m_field_relations ; }

  /** \} */
  //------------------------------------
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
  void put_property( PropertyBase & , Part & );

  /** \} */
  //------------------------------------
  /** \name  Meta data manager construction and destruction
   *  \{
   */

  /** \brief  Construct a meta data manager to own parts and fields.  */
  explicit MetaData( const std::vector<std::string>& entity_type_names );

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

  /** \name  Frequently asserted conditions.
   * \{
   */
  void assert_committed( const char * ) const ;

  void assert_not_committed( const char * ) const ;

  void assert_same_mesh_meta_data( const char * , const MetaData & ) const ;

  void assert_entity_type( const char * , unsigned ) const ;

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
                                  unsigned         arg_entity_type ,
                                  const Part     & arg_part ,
                                  const unsigned * arg_stride );
  /** \} */
private:

  MetaData();                                  ///< \brief  Not allowed
  MetaData( const MetaData & );                ///< \brief  Not allowed
  MetaData & operator = ( const MetaData & );  ///< \brief  Not allowed

  bool   m_commit ;
  Part * m_universal_part ; /* Subset list contains all other parts */
  Part * m_uses_part ;
  Part * m_owns_part ;

  std::vector< FieldBase * >   m_fields ;
  std::vector< FieldRelation > m_field_relations ;
  std::vector< PropertyBase* > m_properties ;
  std::vector< std::string >   m_entity_type_names ;


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
 *
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part );

/** \brief  Declare a  field to exist
 *          for a given entity type and Part.
 */
template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
                        const Part & part ,
                        unsigned     n1 ,
                        unsigned     n2 ,
                        unsigned     n3 ,
                        unsigned     n4 ,
                        unsigned     n5 ,
                        unsigned     n6 );

template< class field_type >
field_type & put_field( field_type & field ,
                        unsigned  entity_type ,
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
{ return * m_universal_part->m_subsets[ord] ; }

template< class field_type >
inline
field_type * MetaData::get_field( const std::string & name ) const
{
  typedef FieldTraits< field_type > Traits ;

  const DataTraits & dt = data_traits< typename Traits::data_type >();

  const shards::ArrayDimTag * tags[8] ;

  Traits::assign_tags( tags );

  FieldBase * const field =
    stk::mesh::get_field( "stk::mesh::MetaData::get_field" ,
                          name , dt , Traits::Rank , tags , 0 , m_fields );

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
  unsigned entity_type ,
  const Part & part )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride );

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
                        const Part &part ,
                        unsigned    n1 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 );

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 );

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
                        const Part &part ,
                        unsigned    n1 ,
                        unsigned    n2 ,
                        unsigned    n3 )
{
  typedef FieldTraits< field_type > Traits ;
  typedef typename Traits::Helper   Helper ;

  unsigned stride[8] ;

  Helper::assign( stride , n1 , n2 , n3 );

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
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

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
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

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
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

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template< class field_type >
inline
field_type & put_field( field_type &field ,
                        unsigned entity_type ,
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

  field.mesh_meta_data().
    declare_field_restriction( field, entity_type, part, stride);

  return field ;
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( Part & p , const T * a )
{
  assert_not_committed( "stk::mesh::MetaData::declare_attribute_with_delete" );
  return p.m_attribute.template insert_with_delete<T>( a );
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( Part & p , const T * a )
{
  assert_not_committed( "stk::mesh::MetaData::declare_attribute_no_delete" );
  return p.m_attribute.template insert_no_delete<T>( a );
}

template<class T>
inline
const T *
MetaData::declare_attribute_with_delete( FieldBase & f , const T * a )
{
  assert_not_committed( "stk::mesh::MetaData::declare_attribute_with_delete" );
  return f.m_attribute.template insert_with_delete<T>( a );
}

template<class T>
inline
const T *
MetaData::declare_attribute_no_delete( FieldBase & f , const T * a )
{
  assert_not_committed( "stk::mesh::MetaData::declare_attribute_no_delete" );
  return f.m_attribute.template insert_no_delete<T>( a );
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
void MetaData::put_property( PropertyBase & prop , Part & part )
{
  prop.add_property( part.mesh_meta_data_ordinal() );
}

} // namespace mesh
} // namespace stk

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_MetaData_hpp */

