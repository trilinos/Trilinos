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
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>  // for FieldBaseImpl
#include <string>                       // for string
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc

namespace shards { class ArrayDimTag; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class UnitTestFieldImpl; } }
namespace stk { namespace mesh { namespace impl { class FieldRepository; } } }
namespace stk { namespace mesh { class DataTraits; } }
namespace stk { namespace mesh { template<typename T> class DeviceField; } }

namespace stk {
namespace mesh {

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
  void sync_to_host() const { m_impl.sync_to_host(); }
  void sync_to_device() const { m_impl.sync_to_device(); }
  void clear_sync_state() const { m_impl.clear_sync_state(); }

  unsigned synchronized_count() const
  {
    if (get_ngp_field() != nullptr) {
      return get_ngp_field()->synchronized_count();
    }
    return get_mesh().synchronized_count();
  }

  size_t num_syncs_to_host() const { return m_impl.num_syncs_to_host(); }
  size_t num_syncs_to_device() const { return m_impl.num_syncs_to_device(); }

  void rotate_multistate_data();

#ifdef STK_DEBUG_FIELD_SYNC
  unsigned get_num_components(const FastMeshIndex & index) const;
  bool has_device_field_debug_data() const { return ((get_debug_ngp_field() != nullptr) && (m_debugFieldLastModification.extent(0) > 0)); }
  void set_last_modification_view(const LastFieldModLocationType & lastModView) const;
  bool last_accessed_entity_value_has_changed(const FastMeshIndex & index, unsigned component) const;
  void set_last_modification(const FastMeshIndex & index, unsigned component, LastModLocation location) const;
  void detect_host_field_entity_modification() const;
  bool last_entity_modification_not_on_host(const FastMeshIndex & index, unsigned component) const;
  void store_last_entity_access_location(const FastMeshIndex & index) const;
  void check_stale_field_entity_access(const FastMeshIndex & index, const char * fileName, int lineNumber) const;

  void store_last_bucket_access_location(const stk::mesh::Bucket & bucket) const;
  void detect_host_field_bucket_modification() const;
  void check_stale_field_bucket_access(const stk::mesh::Bucket & bucket, const char * fileName, int lineNumber) const;

  void reset_debug_state();
  NgpFieldBase * get_debug_ngp_field() const { return m_impl.get_ngp_field(); }
#endif

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

  template <typename T> friend NgpField<T> & get_updated_ngp_field(const FieldBase & stkField);
  template <typename T> friend class DeviceField;

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

#ifdef STK_DEBUG_FIELD_SYNC
  mutable LastFieldModLocationType m_debugFieldLastModification;
  mutable FastMeshIndex            m_lastFieldEntityLocation = {INVALID_BUCKET_ID, 0};
  mutable Bucket *                 m_lastFieldBucketLocation = nullptr;
  mutable std::vector<uint8_t>     m_lastFieldValue;
#endif
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

#ifdef STK_DEBUG_FIELD_SYNC

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord, const int knownSize,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssertMsg(f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity >= knownSize, "field name= " << f.name() << "knownSize= " << knownSize << ", m_bytes_per_entity= " << f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity);
  ThrowAssert(f.get_meta_data_for_field()[bucket_id].m_data != NULL);

  if (f.has_device_field_debug_data()) {
    stk::mesh::FastMeshIndex index {bucket_id, static_cast<unsigned>(bucket_ord)};
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_entity_modification();
    f.check_stale_field_entity_access(index, fileName, lineNumber);
    f.store_last_entity_access_location(index);
  }

  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data + knownSize * bucket_ord);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  if (f.has_device_field_debug_data()) {
    stk::mesh::Bucket & bucket = *f.get_mesh().buckets(f.entity_rank())[bucket_id];
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_bucket_modification();
    f.check_stale_field_bucket_access(bucket, fileName, lineNumber);
    f.store_last_bucket_access_location(bucket);
  }

  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  if (f.has_device_field_debug_data()) {
    stk::mesh::FastMeshIndex index {bucket_id, static_cast<unsigned>(bucket_ord)};
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_entity_modification();
    f.check_stale_field_entity_access(index, fileName, lineNumber);
    f.store_last_entity_access_location(index);
  }

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[bucket_id];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * bucket_ord);
}


template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b, Bucket::size_type bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&f.get_mesh() == &b.mesh());

  if (f.has_device_field_debug_data()) {
    stk::mesh::FastMeshIndex index {b.bucket_id(), static_cast<unsigned>(bucket_ord)};
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_entity_modification();
    f.check_stale_field_entity_access(index, fileName, lineNumber);
    f.store_last_entity_access_location(index);
  }

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * bucket_ord);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&b.mesh() == &f.get_mesh());

  if (f.has_device_field_debug_data()) {
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_bucket_modification();
    f.check_stale_field_bucket_access(b, fileName, lineNumber);
    f.store_last_bucket_access_location(b);
  }

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
ngp_debug_field_data(const FieldType & f, const Bucket& b)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&b.mesh() == &f.get_mesh());
  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, Entity e,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  const MeshIndex& mi = f.get_mesh().mesh_index(e);
  ThrowAssert(f.entity_rank() == mi.bucket->entity_rank());
  ThrowAssert(&f.get_mesh() == &mi.bucket->mesh());

  if (f.has_device_field_debug_data()) {
    stk::mesh::FastMeshIndex index {mi.bucket->bucket_id(), mi.bucket_ordinal};
    f.get_debug_ngp_field()->detect_device_field_modification();
    f.detect_host_field_entity_modification();
    f.check_stale_field_entity_access(index, fileName, lineNumber);
    f.store_last_entity_access_location(index);
  }

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[mi.bucket->bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * mi.bucket_ordinal);
}

#else

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord, const int knownSize)
{
  ThrowAssertMsg(f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity >= knownSize, "field name= " << f.name() << "knownSize= " << knownSize << ", m_bytes_per_entity= " << f.get_meta_data_for_field()[bucket_id].m_bytes_per_entity);
  ThrowAssert(f.get_meta_data_for_field()[bucket_id].m_data != NULL);

  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data + knownSize * bucket_ord);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id)
{
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(f.get_meta_data_for_field()[bucket_id].m_data);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const unsigned bucket_id, Bucket::size_type bucket_ord)
{
  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[bucket_id];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * bucket_ord);
}


template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b, Bucket::size_type bucket_ord)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&f.get_mesh() == &b.mesh());
  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * bucket_ord);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, const Bucket& b)
{
  ThrowAssert(f.entity_rank() == b.entity_rank());
  ThrowAssert(&b.mesh() == &f.get_mesh());
  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data);
}

template<class FieldType>
inline
typename FieldTraits<FieldType>::data_type*
field_data(const FieldType & f, Entity e)
{
  const MeshIndex& mi           = f.get_mesh().mesh_index(e);

  ThrowAssert(f.entity_rank() == mi.bucket->entity_rank());
  ThrowAssert(&f.get_mesh() == &mi.bucket->mesh());
  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[mi.bucket->bucket_id()];
  return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_bytes_per_entity * mi.bucket_ordinal);
}

#endif

} //namespace mesh
} //namespace stk

#endif //stk_mesh_base_FieldBase_hpp
