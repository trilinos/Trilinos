// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef stk_mesh_base_FieldBase_hpp
#define stk_mesh_base_FieldBase_hpp

#include <stddef.h>                     // for NULL
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, Bucket::size_type
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldRestriction.hpp>  // for FieldRestriction, etc
#include <stk_mesh/base/FieldState.hpp>  // for FieldState
#include <stk_mesh/base/Types.hpp>      // for FieldTraits, EntityRank, etc
#include <stk_mesh/base/NgpFieldBase.hpp>
#include <stk_mesh/base/FieldSyncDebugging.hpp>
#include <stk_mesh/base/StkFieldSyncDebugger.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>  // for FieldBaseImpl
#include <string>                       // for string
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include "Teuchos_any.hpp"
#include <type_traits>

namespace shards { class ArrayDimTag; }

namespace stk {
namespace mesh {

class BulkData;
class MetaData;
class UnitTestFieldImpl;
class DataTraits;
class FieldBase;
template<typename T, template <typename> class NgpDebugger> class DeviceField;

namespace impl {
class FieldRepository;
NgpFieldBase* get_ngp_field(const FieldBase & field);
void set_ngp_field(const FieldBase & stkField, NgpFieldBase * ngpField);
}

struct FieldMetaData
{
  unsigned char* m_data;
  int m_bytes_per_entity;   // num bytes per entity, 0 means bucket does not have this field
};

typedef std::vector<FieldMetaData> FieldMetaDataVector;

//----------------------------------------------------------------------
/** \ingroup stk_stk_mesh_module
 *  \brief  Field base class with an anonymous data type and
 *          anonymous multi-dimension.
 *
 * This class is the base class for all Fields. It defines the member-function
 * API for Field objects. Note that the entire Field API is much broader than
 * what's defined in this class. See Field.hpp for more.
 */
// Implementation details:
//   Simply wraps a FieldBaseImpl object which is kept as a member, all calls are
//   inlined and passed through to the impl object. This design is analogous
//   to the "P-impl" pattern.
class FieldBase
{
public:
  FieldBase() = delete;
  FieldBase(const FieldBase &) = delete;
  FieldBase & operator=(const FieldBase &) = delete;

   /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & mesh_meta_data() const { return m_impl.meta_data(); }
  /** \brief  Internally generated ordinal of this field that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_impl.ordinal(); }

  /** \brief  Application-defined text name of this field */
  const std::string & name() const { return m_impl.name() ; }

  /** \brief  Query if the type is Type */
  template<class Type> bool type_is() const
  { return m_impl.type_is<Type>(); }

  /** \brief  Return the \ref stk::mesh::DataTraits "data traits"
   *          for this field's type
   */
  const DataTraits & data_traits() const { return m_impl.data_traits() ; }

  /** \brief  Number of states of this field */
  unsigned number_of_states() const { return m_impl.number_of_states() ; }

  /** \brief  FieldState of this field */
  FieldState state() const { return m_impl.state() ; }

  /** \brief  Multi-dimensional array rank of this field,
   *          which is zero for a scalar field.
   */
  unsigned field_array_rank() const { return m_impl.field_array_rank(); }

  EntityRank entity_rank() const { return m_impl.entity_rank(); }

  /** \brief  Multi-dimensional
   *          \ref shards::ArrayDimTag "array dimension tags"
   *          of this field.
   */
  const shards::ArrayDimTag * const * dimension_tags() const
  { return m_impl.dimension_tags() ; }

  /** \brief  Maximum field data allocation size declared for this
   *          field for the given entity rank.
   */
  unsigned max_size( EntityRank ent_rank) const {
    return m_impl.max_size( ent_rank );
  }

  //----------------------------------------

  /** \brief  Query attribute that has been attached to this field */
  template<class A>
  const A * attribute() const { return m_impl.attribute<A>(); }

  typedef FieldRestriction Restriction;

  /** \brief  A fields' restrictions are maintained in a std::vector */
  typedef FieldRestrictionVector RestrictionVector;

  /** \brief  Vector of field restriction which is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const RestrictionVector &restrictions() const {
    return m_impl.restrictions();
  }

  //----------------------------------------

  FieldBase * field_state( FieldState fstate) const {
    return m_impl.field_state(fstate);
  }

  bool is_state_valid(FieldState fstate) const {
    return field_state(fstate) != nullptr;
  }

  const void* get_initial_value() const { return m_impl.get_initial_value(); }

  void* get_initial_value() { return m_impl.get_initial_value(); }

  unsigned get_initial_value_num_bytes() const { return m_impl.get_initial_value_num_bytes(); }

  virtual ~FieldBase() {}

  virtual std::ostream& print_data(std::ostream& out, void* data, unsigned size_per_entity) const { return out; }

  stk::mesh::BulkData& get_mesh() const
  { return *m_mesh; }

  // This should only be accessed by the stk::mesh::BulkData class
  inline const FieldMetaDataVector& get_meta_data_for_field() const {
    return m_field_meta_data;
  }

  // This should only be accessed by the stk::mesh::BulkData class
  inline FieldMetaDataVector& get_meta_data_for_field() {
    return m_field_meta_data;
  }

  unsigned length(const stk::mesh::Part& part) const;

  template<typename PARTVECTOR>
  bool defined_on_any(const PARTVECTOR& parts) const
  {
    bool defined_on_any_part = false;
    size_t i = 0;
    while(!defined_on_any_part && i < parts.size()) {
      defined_on_any_part = defined_on_any_part || defined_on(*parts[i]);
      ++i;
    }

    return defined_on_any_part;
  }

  bool defined_on(const stk::mesh::Part& part) const;

  void modify_on_host() const { m_impl.modify_on_host(); }
  void modify_on_device() const { m_impl.modify_on_device(); }
  void modify_on_host(const Selector& s) const { m_impl.modify_on_host(s); }
  void modify_on_device(const Selector& s) const { m_impl.modify_on_device(s); }
  void sync_to_host() const { m_impl.sync_to_host(); }
  void sync_to_device() const { m_impl.sync_to_device(); }
  void clear_sync_state() const { m_impl.clear_sync_state(); }
  void clear_host_sync_state() const { m_impl.clear_host_sync_state(); }
  void clear_device_sync_state() const { m_impl.clear_device_sync_state(); }

  unsigned synchronized_count() const
  {
    if (get_ngp_field() != nullptr) {
      return get_ngp_field()->synchronized_count();
    }
    return get_mesh().synchronized_count();
  }

  size_t num_syncs_to_host() const { return m_impl.num_syncs_to_host(); }
  size_t num_syncs_to_device() const { return m_impl.num_syncs_to_device(); }
  bool has_ngp_field() const { return m_impl.get_ngp_field() != nullptr; }

  template <typename StkDebugger>
  void make_field_sync_debugger() const {
    m_stkFieldSyncDebugger = Teuchos::any(StkDebugger(this));
  }

  template <typename StkDebugger>
  StkDebugger & get_field_sync_debugger() const {
    return Teuchos::any_cast<StkDebugger>(m_stkFieldSyncDebugger);
  }

  void rotate_multistate_data();


private:

  //  Associate this field with a bulk data.
  //  Note, a field can be associated with one and only one bulk data object
  void set_mesh(stk::mesh::BulkData* bulk);

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & meta_data() const { return m_impl.meta_data(); }

  NgpFieldBase * get_ngp_field() const { return m_impl.get_ngp_field(); }
  void set_ngp_field(NgpFieldBase * ngpField) const { m_impl.set_ngp_field(ngpField); }

  void increment_num_syncs_to_host() const { m_impl.increment_num_syncs_to_host(); }
  void increment_num_syncs_to_device() const { m_impl.increment_num_syncs_to_device(); }

  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::BulkData;
  friend class ::stk::mesh::impl::FieldRepository;
  friend class ::stk::mesh::impl::FieldBaseImpl;

  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestFieldImpl ;

  friend NgpFieldBase* impl::get_ngp_field(const FieldBase & stkField);
  friend void impl::set_ngp_field(const FieldBase & stkField, NgpFieldBase * ngpField);

  template <typename T, template <typename> class NgpDebugger> friend class HostField;
  template <typename T, template <typename> class NgpDebugger> friend class DeviceField;

  FieldMetaDataVector m_field_meta_data;

protected:
  FieldBase(MetaData                   * arg_mesh_meta_data,
            unsigned                     arg_ordinal,
            const std::string          & arg_name,
            const DataTraits           & arg_traits,
            unsigned                     arg_rank,
            const shards::ArrayDimTag  * const * arg_dim_tags,
            unsigned                     arg_number_of_states,
            FieldState                   arg_this_state)
    : m_mesh(nullptr),
      m_impl(arg_mesh_meta_data,
             stk::topology::INVALID_RANK,
             arg_ordinal,
             arg_name,
             arg_traits,
             arg_rank,
             arg_dim_tags,
             arg_number_of_states,
             arg_this_state)
  {}

  FieldBase(MetaData                   * arg_mesh_meta_data,
            stk::topology::rank_t        arg_entity_rank,
            unsigned                     arg_ordinal,
            const std::string          & arg_name,
            const DataTraits           & arg_traits,
            unsigned                     arg_rank,
            const shards::ArrayDimTag  * const * arg_dim_tags,
            unsigned                     arg_number_of_states,
            FieldState                   arg_this_state)
    : m_mesh(nullptr),
      m_impl(arg_mesh_meta_data,
             arg_entity_rank,
             arg_ordinal,
             arg_name,
             arg_traits,
             arg_rank,
             arg_dim_tags,
             arg_number_of_states,
             arg_this_state)
  {}

private:
  stk::mesh::BulkData* m_mesh;
  impl::FieldBaseImpl  m_impl;
  mutable Teuchos::any m_stkFieldSyncDebugger;
};

/** \brief  Print the field type, text name, and number of states. */
std::ostream & operator << ( std::ostream & , const FieldBase & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const ,
                      const FieldBase & );


std::ostream & print_restrictions( std::ostream & ,
                                   const char * const ,
                                   const FieldBase & );

namespace impl {
inline NgpFieldBase* get_ngp_field(const FieldBase & stkField) {
  return stkField.get_ngp_field();
}

inline void set_ngp_field(const FieldBase & stkField, NgpFieldBase * ngpField) {
  stkField.set_ngp_field(ngpField);
}
}

//
//  Field free access methods
//

inline unsigned field_bytes_per_entity(const FieldBase& f, const Bucket& b) {
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&f.get_mesh() == &b.mesh());
  return f.get_meta_data_for_field()[b.bucket_id()].m_bytes_per_entity;
}

inline unsigned field_bytes_per_entity(const FieldBase& f, unsigned bucket_id) {
  ThrowAssert(bucket_id < f.get_meta_data_for_field().size());
  return f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity;
}

inline unsigned field_bytes_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return field_bytes_per_entity(f, bulk.bucket(e));
}

inline bool is_matching_rank(const FieldBase& f, const Bucket& b) {
  ThrowAssert(&f.get_mesh() == &b.mesh());
  return(b.entity_rank() == f.entity_rank());
}

inline bool is_matching_rank(const FieldBase& f, Entity e) {
  return is_matching_rank(f, f.get_mesh().bucket(e));
}

inline bool is_matching_rank(const FieldBase& f, EntityRank rank) {
  return f.entity_rank() == rank;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, const Bucket& b) {
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&f.get_mesh() == &b.mesh());
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  return f.get_meta_data_for_field()[b.bucket_id()].m_bytes_per_entity/bytes_per_scalar;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, unsigned bucket_id) {
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  return f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity/bytes_per_scalar;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, Entity e) {
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  BulkData& bulk(f.get_mesh());
  ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return field_bytes_per_entity(f, bulk.bucket(e))/bytes_per_scalar;
}

inline bool field_is_allocated_for_bucket(const FieldBase& f, const Bucket& b) {
  ThrowAssert(&b.mesh() == &f.get_mesh());
  //return true if field and bucket have the same rank and the field is associated with the bucket
  ThrowAssert(f.get_meta_data_for_field().size() > b.bucket_id());
  return (is_matching_rank(f, b) && 0 != f.get_meta_data_for_field()[b.bucket_id()].m_bytes_per_entity);
}

inline size_t get_total_ngp_field_allocation_bytes(const FieldBase & f) {
  const Selector selector = stk::mesh::selectField(f);
  const BucketVector & buckets = f.get_mesh().get_buckets(f.entity_rank(), selector);
  const size_t numBuckets = buckets.size();
  const size_t bucketCapacity = (buckets.empty()) ? 0 : buckets[0]->capacity();
  const size_t numPerEntity = f.max_size(f.entity_rank());
  const size_t bytesPerScalar = f.data_traits().size_of;

  return numBuckets * bucketCapacity * numPerEntity * bytesPerScalar;
}

struct FieldBasePtrLess {
  bool operator()(const stk::mesh::FieldBase* lhs, const stk::mesh::FieldBase* rhs) const
  {
    return lhs->mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal();
  }
};


//
//  Optimized field data access, here the size of the field data is passed in rather than looked up.
//  This accessor can be used if the field is known to exist everywhere and known to have the same
//  size everywhere.
//

template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, EmptyStkFieldSyncDebugger>::value, void>::type
debug_stale_access_entity_check(const FieldBase & stkField, const Entity & entity,
                                const char * fileName, int lineNumber)
{
}

template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, StkFieldSyncDebugger>::value, void>::type
debug_stale_access_entity_check(const FieldBase & stkField, const Entity & entity,
                                const char * fileName, int lineNumber)
{
  if (stkField.has_ngp_field()) {
    stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_entity_check(entity, fileName, lineNumber);
  }
}


template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, EmptyStkFieldSyncDebugger>::value, void>::type
debug_stale_access_entity_check(const FieldBase & stkField, const unsigned bucketId, Bucket::size_type bucketOrd,
                                const char * fileName, int lineNumber)
{
}

template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, StkFieldSyncDebugger>::value, void>::type
debug_stale_access_entity_check(const FieldBase & stkField, const unsigned bucketId, Bucket::size_type bucketOrd,
                                const char * fileName, int lineNumber)
{
  if (stkField.has_ngp_field()) {
    stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_entity_check(bucketId, bucketOrd, fileName, lineNumber);
  }
}


template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, EmptyStkFieldSyncDebugger>::value, void>::type
debug_stale_access_bucket_check(const FieldBase & stkField, const Bucket & bucket,
                                const char * fileName, int lineNumber)
{
}

template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, StkFieldSyncDebugger>::value, void>::type
debug_stale_access_bucket_check(const FieldBase & stkField, const Bucket & bucket,
                                const char * fileName, int lineNumber)
{
  if (stkField.has_ngp_field()) {
    stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_bucket_check(bucket, fileName, lineNumber);
  }
}


template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, EmptyStkFieldSyncDebugger>::value, void>::type
debug_stale_access_bucket_check(const FieldBase & stkField, const unsigned & bucketId,
                                const char * fileName, int lineNumber)
{
}

template <typename StkDebugger>
typename std::enable_if<std::is_same<StkDebugger, StkFieldSyncDebugger>::value, void>::type
debug_stale_access_bucket_check(const FieldBase & stkField, const unsigned & bucketId,
                                const char * fileName, int lineNumber)
{
  if (stkField.has_ngp_field()) {
    stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_bucket_check(bucketId, fileName, lineNumber);
  }
}


template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord, const int knownSize,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssertMsg(f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity >= knownSize,
                 "field name= " << f.name() << "knownSize= " << knownSize << ", m_bytes_per_entity= "
                 << f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity);
  ThrowAssert(f.get_meta_data_for_field()[bucket_id].m_data != NULL);

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, bucket_ord, fileName, lineNumber);

  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data +
                                                                       knownSize * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  debug_stale_access_bucket_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, fileName, lineNumber);

  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, bucket_ord, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[bucket_id];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b, Bucket::size_type bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&f.get_mesh() == &b.mesh());

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), b[bucket_ord], fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data +
                                                                       field_meta_data.m_bytes_per_entity * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&b.mesh() == &f.get_mesh());

  debug_stale_access_bucket_check<StkDebugger>(static_cast<const FieldBase&>(f), b, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, Entity e,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  const MeshIndex& mi = f.get_mesh().mesh_index(e);
  ThrowAssertMsg(f.entity_rank() == mi.bucket->entity_rank(),
                 "field_data called with " << f.entity_rank() << " field (" << f.name() << ") and different-rank entity "
                 << f.get_mesh().entity_key(e) << ". The rank of the field and entity must match.");
  ThrowAssert(&f.get_mesh() == &mi.bucket->mesh());

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), e, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[mi.bucket->bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data +
                                                                       field_meta_data.m_bytes_per_entity * mi.bucket_ordinal);
}

} //namespace mesh
} //namespace stk

#endif //stk_mesh_base_FieldBase_hpp
