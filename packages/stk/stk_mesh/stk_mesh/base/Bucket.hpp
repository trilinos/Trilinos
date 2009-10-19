#ifndef stk_mesh_Bucket_hpp
#define stk_mesh_Bucket_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Part.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print the \ref stk::mesh::Part "part" names
 *          for which this bucket is a subset.
 */
std::ostream & operator << ( std::ostream & , const Bucket & );

/** \brief  Print the parts and entities of this bucket */
std::ostream &
print( std::ostream & , const std::string & indent , const Bucket & );

//----------------------------------------------------------------------
/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part"
 */
bool has_superset( const Bucket & ,  const Part & );

/** \brief  Is this bucket a subset of all of the given
 *          \ref stk::mesh::Part "parts"
 */
bool has_superset( const Bucket & , const PartVector & );

//----------------------------------------------------------------------
/** \brief  A random access iterator for a
 *  \ref stk::mesh::Bucket "bucket" that dereferences to a
 *  \ref stk::mesh::Entity "entity" reference.
 */
class BucketIterator : std::iterator<std::random_access_iterator_tag,Entity > {
private:
  const Bucket * m_bucket_ptr;
  size_t         m_current_entity;

  inline Entity & entity( const size_t ) const ;

  void throw_error(const char *err) const;

  template< class field_type >
    friend
    typename FieldTraits< field_type >::data_type *
    field_data( const field_type & f , const BucketIterator &i );

  template< class field_type > friend struct BucketArray ;

public:

  /** \brief Constructor 
    * \param bucket_ptr \ref stk::mesh::Bucket "bucket" pointer
    * \param offset int
    */
  template< typename intType >
  BucketIterator(const Bucket * const bucket_ptr, intType offset) {
    m_bucket_ptr = bucket_ptr;
    m_current_entity = offset;
  }

  /** \brief Default constructor */
  BucketIterator() { 
    m_bucket_ptr = NULL;
    m_current_entity = 0;
  }

  /** \brief Copy Constructor */
  BucketIterator(const BucketIterator &i) {
    m_bucket_ptr = i.m_bucket_ptr;
    m_current_entity = i.m_current_entity;
  }

  /** \brief Assignment operator */
  BucketIterator & operator=(const BucketIterator &i) {
    m_bucket_ptr = i.m_bucket_ptr;
    m_current_entity = i.m_current_entity;
    return *this;
  }

  /** \brief Dereference operator
    * \return \ref stk::mesh::Entity "entity" reference
   */
  inline Entity & operator*() const { return entity(0); }

  /** \brief Pointer operator
    * \return \ref stk::mesh::Entity "entity" pointer
   */
  inline Entity * operator->() const { return & entity(0); }

  /** \brief Pre increment */
  inline BucketIterator & operator++() {
    ++m_current_entity;
    return *this;
  }

  /** \brief Pre decrement */
  inline BucketIterator & operator--() {
    --m_current_entity;
    return *this;
  }

  /** \brief Post increment */
  inline BucketIterator operator++(int) {
    BucketIterator temp = *this;
    ++m_current_entity;
    return temp;
  }

  /** \brief Post decrement */
  inline BucketIterator operator--(int) {
    BucketIterator temp = *this;
    --m_current_entity;
    return temp;
  }

  /** \brief Less than */
  inline bool operator<(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator < given iterator from different bucket");
    return (m_current_entity < i.m_current_entity);
  }

  /** \brief Less than equal too */
  inline bool operator<=(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator <= given iterator from different bucket");
    return (m_current_entity <= i.m_current_entity);
  }

  /** \brief Greater than  */
  inline bool operator>(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator > given iterator from different bucket");
    return (m_current_entity > i.m_current_entity);
  }

  /** \brief Greater than equal too */
  inline bool operator>=(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator >= given iterator from different bucket");
    return (m_current_entity >= i.m_current_entity);
  }

  /** \brief Equal too */
  inline bool operator==(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator == given iterator from different bucket");
    return (m_current_entity == i.m_current_entity);
  }

  /** \brief Not equal */
  inline bool operator!=(const BucketIterator &i) const { 
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator != given iterator from different bucket");
    return (m_current_entity != i.m_current_entity);
  }

  inline BucketIterator & operator+=(int n) { 
    m_current_entity += n;
    return *this;
  }

  inline BucketIterator & operator-=(int n) { 
    m_current_entity -= n;
    return *this;
  }

  inline BucketIterator operator+(int n) const {
    return BucketIterator(m_bucket_ptr, m_current_entity + n);
  }

  inline BucketIterator operator-(int n) const {
    return BucketIterator(m_bucket_ptr, m_current_entity - n);
  }

  /** \brief Distance between iterators */
  inline ptrdiff_t operator-(const BucketIterator &i) const {
    if (m_bucket_ptr != i.m_bucket_ptr)
      throw_error("operator - given iterator from different bucket");
    return static_cast<ptrdiff_t>(m_current_entity - i.m_current_entity);
  }

  template< typename intType >
  inline Entity & operator[]( const intType & n ) const { return entity(n); }

}; // class BucketIterator

//----------------------------------------------------------------------
/** \brief  A container for the \ref stk_mesh_field_data "field data"
 *          of a homogeneous collection of 
 *          \ref stk::mesh::Entity "entities".
 *
 *  The entities are homogeneous in that they are of the same entity type
 *  and are members of the same of parts.
 */
class Bucket {
private:
  friend class BulkData ;

  struct DataMap {
    typedef FieldBase::Restriction::size_type size_type ;

    const size_type * m_stride ;
    size_type         m_base ;
    size_type         m_size ;
  };

  BulkData             & m_mesh ;        // Where this bucket resides
  const unsigned         m_entity_type ; // Type of entities for this bucket
  const unsigned * const m_key ;         // Unique key in the bulk data
  const size_t           m_alloc_size ;  // Allocation size of this bucket
  const size_t           m_capacity ;    // Capacity for entities
  size_t                 m_size ;        // Number of entities
  Bucket               * m_bucket ;
  DataMap        * const m_field_map ;   // Field value data map, shared
  Entity        ** const m_entities ;    // Array of entity pointers,
                                         // begining of field value memory.
public:

  //--------------------------------
  // Container-like types and methods:

  typedef BucketIterator iterator ;

  /** \brief Beginning of the bucket */
  inline iterator begin() const { return iterator(this,0); }

  /** \brief End of the bucket */
  inline iterator end() const { return iterator(this,m_size); }

  /** \brief  Number of entities associated with this bucket */
  size_t size() const { return m_size ; }

  /** \brief  Capacity of this bucket */
  size_t capacity() const { return m_capacity ; }

  /** \brief  Query the i^th entity */
  Entity & operator[] ( size_t i ) const { return *(m_entities[i]) ; }

  //--------------------------------
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & mesh() const { return m_mesh ; }

  /** \brief  Type of entities in this bucket */
  unsigned entity_type() const { return m_entity_type ; }

  /** \brief  This bucket is a subset of these \ref stk::mesh::Part "parts" */
  void supersets( PartVector & ) const ;

  //--------------------------------
  /** \brief  Bucket is a subset of the given part */
  bool member( const Part & ) const ;

  /** \brief  Bucket is a subset of all of the given parts */
  bool member_all( const std::vector<Part*> & ) const ;

  /** \brief  Bucket is a subset of any of the given parts */
  bool member_any( const std::vector<Part*> & ) const ;

  //--------------------------------
  /** Query bucket's supersets' ordinals. */
  std::pair<const unsigned *, const unsigned *>
    superset_part_ordinals() const
    {
      return std::pair<const unsigned *, const unsigned *>
             ( m_key + 1 , m_key + m_key[0] );
    }

#ifndef DOXYGEN_COMPILE
  const unsigned * key() const { return m_key ; }
#endif /* DOXYGEN_COMPILE */

  /** \brief  The allocation size, in bytes, of this bucket */
  unsigned allocation_size() const { return m_alloc_size ; }

private:

  ~Bucket();
  Bucket();
  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket( BulkData & ,
          unsigned          arg_entity_type ,
          const unsigned  * arg_key ,
          size_t            arg_alloc_size ,
          size_t            arg_capacity ,
          Bucket::DataMap * arg_field_map , 
          Entity         ** arg_entity_array );

  void update_state();

  static void copy_fields( Bucket & k_dst , unsigned i_dst ,
                           Bucket & k_src , unsigned i_src );

  static void zero_fields( Bucket & k_dst , unsigned i_dst );

  static void destroy_bucket( Bucket * );

  static Bucket * declare_nil_bucket( BulkData & , unsigned );

  static Bucket *
    declare_bucket( BulkData & ,
                    const unsigned entity_type ,
                    const unsigned part_count ,
                    const unsigned part_ord[] ,
                    const unsigned bucket_capacity ,
                    const std::vector< FieldBase * > & field_set ,
                          std::vector<Bucket*>       & bucket_set );


  friend
  unsigned field_data_size( const FieldBase & f , const Bucket & k );

  template< class field_type >
    friend
    typename FieldTraits< field_type >::data_type *
    field_data( const field_type & f , const Bucket::iterator &i );

  template< class field_type >
  friend
  typename FieldTraits< field_type >::data_type *
  field_data( const field_type & f , const Entity & e );

  template< class field_type > friend struct EntityArray ;
  template< class field_type > friend struct BucketArray ;
};

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

inline Entity & BucketIterator::entity( const size_t i ) const
{
  if ( ! m_bucket_ptr ) { throw_error("is NULL"); }
  return (*m_bucket_ptr)[ m_current_entity + i ] ;
}

}
}

#endif

