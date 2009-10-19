
#ifndef stk_mesh_Field_hpp
#define stk_mesh_Field_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <Shards_Array.hpp>

#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/util/CSet.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/DataTraits.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Enumeration of states for multi-state
 *          \ref stk::mesh::Field "fields".
 *
 *  A field may be declared to have field data for multiple states.
 *  -  Field states <b>StateNone</b>, <b>StateNew</b>, or <b>StateNP1</b>
 *     refer to the current or newest of a field.
 *  -  Field states <b>StateOld</b> or <b>StateN</b> refer to
 *     the previous state of a field with two or more states.
 *   - The remaining field states <b>StateNM1</b>, <b>StateNM2</b>,
 *     <b>StateNM3</b>, <b>StateNM4</b> refer to prior states
 *     N-1, N-2, N-3, and N-4 accordingly.
 */
enum FieldState {
  StateNone = 0,  ///< \brief State of a field with one state
  StateNew  = 0,  ///< \brief Newest state of a field with two states
  StateNP1  = 0,  ///< \brief Newest state of a field with three+ states
  StateOld  = 1,  ///< \brief Previous state of a field with two states
  StateN    = 1,  ///< \brief Previous state of a field with three+ states
  StateNM1  = 2,  ///< \brief Previous-1 state of a field with three+ states
  StateNM2  = 3,  ///< \brief Previous-2 state of a field with four+ states
  StateNM3  = 4,  ///< \brief Previous-3 state of a field with five+ states
  StateNM4  = 5   ///< \brief Previous-4 state of a field with six states
};

/** \brief Maximum number of states that a \ref stk::mesh::Field "field"
 *         can have.
 */
enum { MaximumFieldStates = 6 };

/** \brief  Return the string name of a field state. */
const char * field_state_name( FieldState );

/** \} */
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \ingroup stk_stk_mesh_module
 *  \brief  Field base class with an anonymous data type and
 *          anonymous multi-dimension.
 */
template<>
class Field< void , void , void , void , void , void , void , void > {
public:
  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & mesh_meta_data() const { return *m_mesh_meta_data ; }

  /** \brief  Internally generated ordinal of this field that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_mesh_meta_data_ordinal ; }

  /** \brief  Application-defined text name of this field */
  const std::string & name() const { return m_name ; }

  /** \brief  Query if the type is Type */
  template<class Type> bool type_is() const
  { return m_data_traits.type_info == typeid(Type) ; }

  /** \brief  Return the \ref stk::mesh::DataTraits "data traits"
   *          for this field's type 
   */
  const DataTraits & data_traits() const { return m_data_traits ; }

  /** \brief  Number of states of this field */
  unsigned number_of_states() const { return m_num_states ; }

  /** \brief  FieldState of this field */
  FieldState state() const { return m_this_state ; }

  /** \brief  Multi-dimensional array rank of this field,
   *          which is zero for a scalar field.
   */
  unsigned rank() const { return m_rank ; }

  /** \brief  Multi-dimensional
   *          \ref shards::ArrayDimTag "array dimension tags"
   *          of this field.
   */
  const shards::ArrayDimTag * const * dimension_tags() const
  { return m_dim_tags ; }

  /** \brief  Maximum field data allocation size declared for this
   *          field for the given entity type.
   */
  unsigned max_size( unsigned entity_type) const ;

  //----------------------------------------

  /** \brief  Query attribute that has been attached to this field */
  template<class A>
  const A * attribute() const { return m_attribute.template get<A>(); }

  /** \brief  Field restriction (i.e. field data allocation rule)
   *          for a field.  An internal class that should never need
   *          to be used <em> directly </em> within an application code.
   */
  struct Restriction {
    typedef shards::array_traits::int_t size_type ;

    /** \brief  Encoding of ( unsigned , Part ordinal ) */
    EntityKey key ;

    /** \brief  \ref shards::Array "Multi-dimensional array" stride */
    size_type stride[ MaximumFieldDimension ];

#ifndef DOXYGEN_COMPILE

    Restriction();
    Restriction( const Restriction & rhs );
    Restriction & operator = ( const Restriction & rhs );

    Restriction( unsigned t , unsigned );
    unsigned type()    const { return entity_type( key ); }
    unsigned ordinal() const { return entity_id( key ); }
#endif /* DOXYGEN_COMPILE */
  };

  /** \brief  A fields' restrictions are maintained in a std::vector */
  typedef std::vector<Restriction> RestrictionVector;
  
  /** \brief  Vector of field restriction which is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const RestrictionVector &restrictions() const ;

  /** \brief  Query a field restriction, result is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const Restriction & restriction( unsigned , const Part & ) const ;

  //----------------------------------------

private:

  /* \brief  A field is owned by a MetaData, as such only the owning
   *         MetaData can create, delete, or modify a field.
   *         The owner-modifies rule is enforced by all non-const
   *         methods being private and the MetaData be a friend.
   */
  friend class ::stk::mesh::MetaData ;
 
  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestMetaData ;

  template< typename Scalar ,
            class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 >
  friend class stk::mesh::Field ;

#ifndef DOXYGEN_COMPILE

  ~Field();

  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

  Field( MetaData * , unsigned arg_ordinal ,
         const std::string & ,
         const DataTraits & ,
         unsigned arg_number_of_states , FieldState arg_this_state );

  static FieldBase *
    declare_field( const std::string                 & arg_name ,
                   const DataTraits                  & arg_traits ,
                   unsigned                            arg_rank ,
                   const shards::ArrayDimTag * const * arg_dim_tags ,
                   unsigned                            arg_num_states ,
                   MetaData                          * arg_meta_data ,
                   std::vector<FieldBase*>           & arg_meta_data_fields );

  void insert_restriction( const char *       arg_method ,
                           EntityType           arg_entity_type ,
                           const Part       & arg_part ,
                           const unsigned   * arg_stride );

  void verify_and_clean_restrictions( const char       * arg_method ,
                                      const PartVector & arg_all_parts );

  RestrictionVector & restrictions();

  //----------------------------------

  const std::string  m_name ;                    ///< Name of the field
  CSet               m_attribute ;               ///< User's attributes
  const DataTraits & m_data_traits ;             ///< Data type traits
  MetaData * const   m_mesh_meta_data ;          ///< Owner of this field
  const unsigned     m_mesh_meta_data_ordinal ;  ///< Ordinal in the field set
  const unsigned     m_num_states ;              ///< Number of states
  const FieldState   m_this_state ;              ///< Field state of this field
  unsigned           m_rank ;                    ///< Number of dimensions
  RestrictionVector  m_dim_map ;                 ///< Only valid on StateNone
  Field *            m_field_states[ MaximumFieldStates ];
  const shards::ArrayDimTag * m_dim_tags[ MaximumFieldDimension ];

#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------
 
void print_field_type( std::ostream                      & arg_msg ,
                       unsigned                            arg_scalar_type ,
                       unsigned                            arg_rank ,
                       const shards::ArrayDimTag * const * arg_tags );

/** \brief  Get a field by name from a vector of fields.
 *          If found verify its type information,
 *          throw and exception if invalid.
 */
FieldBase * get_field(  
  const char                        * arg_method ,
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank , 
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ,
  const std::vector<FieldBase*>     & arg_meta_data_fields );

/** \brief  Print the field type, text name, and number of states. */
std::ostream & operator << ( std::ostream & , const FieldBase & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );

//----------------------------------------------------------------------

/** \ingroup stk_mesh_module
 *  \brief  Field with defined data type and multi-dimensions (if any)
 */
template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class Field : public FieldBase {
public:

  /** \brief  Query this field for a given field state. */
  const Field & field_of_state( FieldState state ) const
  { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

  /** \brief  Query this field for a given field state. */
  Field & field_of_state( FieldState state )
  { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

private:

#ifndef DOXYGEN_COMPILE

  ~Field();
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

#endif /* DOXYGEN_COMPILE */
};


//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
// Internal implementation details to follow.
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \brief  Scalar type and multi-dimensional array traits of a Field */
template< typename Scalar >
struct FieldTraits< Field<Scalar,void,void,void,void,void,void,void> >
{
public:
  typedef shards::array_traits::Helper<Scalar,shards::RankZero,
                                       void,void,void,void,void,void,void,void>
    Helper ;

  typedef Scalar data_type ; ///< \brief  Data type of the field's members
  typedef void   tag1 ;      ///< \brief  Array dimension tag
  typedef void   tag2 ;      ///< \brief  Array dimension tag
  typedef void   tag3 ;      ///< \brief  Array dimension tag
  typedef void   tag4 ;      ///< \brief  Array dimension tag
  typedef void   tag5 ;      ///< \brief  Array dimension tag
  typedef void   tag6 ;      ///< \brief  Array dimension tag
  typedef void   tag7 ;      ///< \brief  Array dimension tag

  /** \brief  Multidimensional array rank */
  enum { Rank = 0 };

  static void assign_tags( const shards::ArrayDimTag ** tags ) {}
};

/** \brief  Scalar type and multi-dimensional array traits of a Field */
template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct FieldTraits< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{
public:
  typedef shards::array_traits::Helper<Scalar,shards::FortranOrder,
                                       Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
    Helper ;
  
  typedef Scalar data_type ; ///< \brief  Data type of the field's members
  typedef Tag1   tag1 ;      ///< \brief  Array dimension tag
  typedef Tag2   tag2 ;      ///< \brief  Array dimension tag
  typedef Tag3   tag3 ;      ///< \brief  Array dimension tag
  typedef Tag4   tag4 ;      ///< \brief  Array dimension tag
  typedef Tag5   tag5 ;      ///< \brief  Array dimension tag
  typedef Tag6   tag6 ;      ///< \brief  Array dimension tag
  typedef Tag7   tag7 ;      ///< \brief  Array dimension tag

  /** \brief  Multidimensional array rank */
  enum { Rank = Helper::Rank };

  static void assign_tags( const shards::ArrayDimTag ** tags )
    { Helper::assign_tags( tags ); }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

inline
FieldBase::Restriction::Restriction()
  : key() { Copy<MaximumFieldDimension>( stride , size_type(0) ); }

inline
FieldBase::Restriction::Restriction( const FieldBase::Restriction & rhs )
  : key( rhs.key ) { Copy< MaximumFieldDimension >( stride , rhs.stride ); }

inline
FieldBase::Restriction &
FieldBase::Restriction::operator = ( const FieldBase::Restriction & rhs )
  {
    key = rhs.key ;
    Copy< MaximumFieldDimension >( stride , rhs.stride );
    return *this ;
  }

inline
FieldBase::Restriction::Restriction( unsigned t , unsigned ord )
  : key( EntityKey( t , ord ) )
    { Copy< MaximumFieldDimension >( stride , size_type(0) ); }

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relation_stencil
 *  \brief  A defined entity-relationship between a field of a pointer type
 *          and the field that it should point to.
 *          An internal class that should never need to be
 *          <em> directly </em> used within application code.
 *
 *  <b> Let </b>
 *  - <b> rel </b> be the \ref stk::mesh::Relation "relation" from
 *    \ref stk::mesh::Entity "entity" <b>e1</b> to
 *    \ref stk::mesh::Entity "entity" <b>e2</b>
 *
 *  <b> If </b>
 *  - \ref stk::mesh::Field "field" <b> m_root </b>
 *    has a pointer scalar type 'T *' AND
 *  - \ref stk::mesh::Field "field" <b> m_target </b>
 *    has a scalar type 'T' AND
 *  - \ref stk_mesh_field_data "field_data"( *m_root , e1 ) exits AND
 *  - \ref stk_mesh_field_data "field_data"( *m_target , e2 ) exits AND
 *  - \ref stk::mesh::Relation "relation" <b> rel </b> is in the domain of
 *    \ref stk_mesh_relation_stencil "relation stencil" <b> m_function </b>
 *
 *  <b> then </b>
 *  <PRE>
 *    index = (*m_function)( e1.entity_type() ,
 *                           e2.entity_type() ,
 *                           rel.identifier() ,
 *                           rel.kind() );
 *
 *    field_data(*m_root,e1)[index] == field_data(*m_target,e2)
 *  </PRE>
 */
struct FieldRelation {
  /** \brief  relation domain part */
  FieldBase          * m_root ;

  /** \brief  relation range part */
  FieldBase          * m_target ;

  /** \brief  \ref stk_mesh_relation_stencil "relation stencil" */
  relation_stencil_ptr m_function ;

#ifndef DOXYGEN_COMPILE

  FieldRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  FieldRelation( const FieldRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  FieldRelation & operator = ( const FieldRelation & rhs )
    {
      m_root = rhs.m_root ;
      m_target = rhs.m_target ;
      m_function = rhs.m_function ;
      return *this ;
    }

#endif /* DOXYGEN_COMPILE */
};

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Field_hpp */

