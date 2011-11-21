/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Bucket_hpp
#define stk_mesh_Bucket_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>
#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/baseImpl/BucketImpl.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {
class BucketRepository;
} // namespace impl


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

// The part count and parts are equal
bool bucket_part_equal( const unsigned * lhs , const unsigned * rhs );

//----------------------------------------------------------------------
/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part"
 */
bool has_superset( const Bucket & ,  const Part & p );

/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part" by partID
 */
bool has_superset( const Bucket & ,  const unsigned & ordinal );

/** \brief  Is this bucket a subset of all of the given
 *          \ref stk::mesh::Part "parts"
 */
bool has_superset( const Bucket & , const PartVector & );


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
  friend class impl::BucketRepository ;
  friend class impl::BucketImpl ;

  impl::BucketImpl       m_bucketImpl ;

public:

  //--------------------------------
  // Container-like types and methods:

  typedef boost::indirect_iterator<Entity*const*,Entity&> iterator ;

  /** \brief Beginning of the bucket */
  inline iterator begin() const { return iterator(m_bucketImpl.begin()); }

  /** \brief End of the bucket */
  inline iterator end() const { return iterator(m_bucketImpl.end()); }

  /** \brief  Number of entities associated with this bucket */
  size_t size() const { return m_bucketImpl.size() ; }

  /** \brief  Capacity of this bucket */
  size_t capacity() const { return m_bucketImpl.capacity() ; }

  /** \brief  Query the i^th entity */
  Entity & operator[] ( size_t i ) const { return m_bucketImpl[i] ; }

  /** \brief  Query the size of this field data specified by FieldBase */
  unsigned field_data_size(const FieldBase & field) const
  { return m_bucketImpl.field_data_size(field); }

  /** \brief  Query the stride of this field data specified by FieldBase */
  const FieldBase::Restriction::size_type * field_data_stride( const FieldBase & field ) const
  { return m_bucketImpl.field_data_stride(field); }

  /** \brief  Query the location of this field data specified by FieldBase and Entity */
  unsigned char * field_data_location( const FieldBase & field, const Entity & entity ) const
  { return m_bucketImpl.field_data_location(field,entity); }

  /** \brief  Query the location of this field data specified by FieldBase */
  unsigned char * field_data_location( const FieldBase & field ) const
  { return m_bucketImpl.field_data_location(field); }

  /** \brief  Query the location of this field data specified by FieldBase and Entity */
  template< class field_type >
  typename FieldTraits< field_type >::data_type *
  field_data( const field_type & field , const Entity & entity ) const
  { return m_bucketImpl.field_data(field,entity.bucket_ordinal()); }

  //--------------------------------
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & mesh() const { return m_bucketImpl.mesh(); }

  /** \brief  Type of entities in this bucket */
  unsigned entity_rank() const { return m_bucketImpl.entity_rank(); }

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
    superset_part_ordinals() const { return m_bucketImpl.superset_part_ordinals() ; }

  /** \brief Equivalent buckets have the same parts
   */
  bool equivalent( const Bucket& b ) const {
    return m_bucketImpl.equivalent(b.m_bucketImpl);
  }

#ifndef DOXYGEN_COMPILE
  const unsigned * key() const { return m_bucketImpl.key() ; }
#endif /* DOXYGEN_COMPILE */

  /** \brief  The allocation size, in bytes, of this bucket */
  unsigned allocation_size() const { return m_bucketImpl.allocation_size() ; }

  /** \brief  A method to assist in unit testing - accesses private data as necessary. */
  bool assert_correct() const;

private:
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & bulk_data() const { return m_bucketImpl.mesh(); }

  ~Bucket();
  Bucket();
  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket( BulkData & arg_mesh ,
          EntityRank arg_entity_rank,
          const std::vector<unsigned> & arg_key,
          size_t arg_capacity
        );


  friend class ::stk::mesh::BulkData;
};



struct BucketLess {
  bool operator()( const Bucket * lhs_bucket , const unsigned * rhs ) const ;
  bool operator()( const unsigned * lhs , const Bucket * rhs_bucket ) const ;
};


inline
std::vector<Bucket*>::iterator
lower_bound( std::vector<Bucket*> & v , const unsigned * key )
{ return std::lower_bound( v.begin() , v.end() , key , BucketLess() ); }



/** \} */

struct To_Ptr : std::unary_function<Entity&, Entity*>
{
  Entity* operator()(Entity& entity) const
  {
    return &entity;
  }
};

// Sometimes, we want a bucket-iterator to dereference to an Entity*
typedef boost::transform_iterator<To_Ptr, Bucket::iterator> BucketPtrIterator;

typedef Bucket::iterator BucketIterator;

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

