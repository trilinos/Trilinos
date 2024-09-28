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

#include <any>
#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldRestriction.hpp>  // for FieldRestriction, etc
#include <stk_mesh/base/FieldState.hpp>  // for FieldState
#include <stk_mesh/base/Types.hpp>      // for EntityRank, etc
#include <stk_mesh/base/NgpFieldBase.hpp>
#include <stk_mesh/base/FieldSyncDebugging.hpp>
#include <stk_mesh/base/StkFieldSyncDebugger.hpp>
#include <string>                       // for string
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include <stk_util/util/SimpleArrayOps.hpp>  // for Copy
#include <stk_util/util/CSet.hpp>
#include <type_traits>

namespace stk::mesh 
{

class BulkData;
class MetaData;
class UnitTestFieldImpl;
class FieldBase;
template<typename T, template <typename> class NgpDebugger> class DeviceField;

namespace impl {
class FieldRepository;
NgpFieldBase* get_ngp_field(const FieldBase & field);
void set_ngp_field(const FieldBase & stkField, NgpFieldBase * ngpField);
stk::CSet & get_attributes(FieldBase & field);
}

struct FieldMetaData
{
  unsigned char* m_data = nullptr;
  int m_bytesPerEntity = 0;
  int m_firstDimension = 0;
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
class FieldBase
{
public:
  using value_type = void;

  FieldBase() = delete;
  FieldBase(const FieldBase &) = delete;
  FieldBase & operator=(const FieldBase &) = delete;

  virtual FieldBase * clone(stk::mesh::impl::FieldRepository & fieldRepo) const = 0;

   /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & mesh_meta_data() const { return *m_meta_data; }

  /** \brief  Internally generated ordinal of this field that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_ordinal; }

  /** \brief  Application-defined text name of this field */
  const std::string & name() const { return m_name ; }

  /** \brief  Query if the type is Type */
  template<class Type> bool type_is() const
  { return m_data_traits.type_info == typeid(Type); }

  /** \brief  Return the \ref stk::mesh::DataTraits "data traits"
   *          for this field's type
   */
  const DataTraits & data_traits() const { return m_data_traits ; }

  /** \brief  Number of states of this field */
  unsigned number_of_states() const { return m_num_states ; }

  /** \brief  FieldState of this field */
  FieldState state() const { return m_this_state; }

  EntityRank entity_rank() const { return m_entity_rank; }

  /** \brief  Maximum field data allocation size declared for this
   *          field for the given entity rank.
   */
  unsigned max_size() const;
  unsigned max_extent(unsigned dimension) const;

  //----------------------------------------

  /** \brief  Query attribute that has been attached to this field */
  template<class A>
  const A * attribute() const { return m_attribute.template get<A>(); }

  typedef FieldRestriction Restriction;

  /** \brief  A fields' restrictions are maintained in a std::vector */
  typedef FieldRestrictionVector RestrictionVector;

  /** \brief  Vector of field restriction which is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const FieldRestrictionVector &restrictions() const { return m_field_states[0]->m_restrictions; }
  FieldRestrictionVector & restrictions() { return m_field_states[0]->m_restrictions; }

  //----------------------------------------

  FieldBase * field_state( FieldState fstate) const {
    return m_field_states[fstate];
  }

  bool is_state_valid(FieldState fstate) const {
    return field_state(fstate) != nullptr;
  }

  const void* get_initial_value() const
  {
    return m_field_states[0]->m_initial_value;
  }
  
  void* get_initial_value() {
    return m_field_states[0]->m_initial_value;
  }
  
  unsigned get_initial_value_num_bytes() const {
    return m_field_states[0]->m_initial_value_num_bytes;
  }

  virtual ~FieldBase() {
    delete m_ngpField;

    if (state() == StateNone) {
      void*& init_val = m_initial_value;

      delete [] reinterpret_cast<char*>(init_val);
      init_val = nullptr;
    }
  }

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

  void modify_on_host() const;
  void modify_on_device() const;
  void modify_on_host(const Selector& s) const;
  void modify_on_device(const Selector& s) const;
  bool need_sync_to_host() const;
  bool need_sync_to_device() const;
  void sync_to_host() const;
  void sync_to_host(const stk::ngp::ExecSpace& newExecSpace) const;
  void sync_to_device() const;
  void sync_to_device(const stk::ngp::ExecSpace& newExecSpace) const;
  void clear_sync_state() const;
  void clear_host_sync_state() const;
  void clear_device_sync_state() const;

  unsigned synchronized_count() const
  {
    if (get_ngp_field() != nullptr) {
      return get_ngp_field()->synchronized_count();
    }
    return get_mesh().synchronized_count();
  }

  size_t num_syncs_to_host() const;
  size_t num_syncs_to_device() const;
  bool has_ngp_field() const { return get_ngp_field() != nullptr; }
  void fence() const;

  template <typename StkDebugger>
  void make_field_sync_debugger() const {
    if (!(m_stkFieldSyncDebugger.has_value())) {
      m_stkFieldSyncDebugger = std::any(StkDebugger(this));
    }
  }

  template <typename StkDebugger>
  StkDebugger & get_field_sync_debugger() const {
    return std::any_cast<StkDebugger&>(m_stkFieldSyncDebugger);
  }

  void rotate_multistate_data(bool rotateNgpFieldViews = false);

 private:
  stk::ngp::ExecSpace& get_execution_space() const {
    return m_execSpace;
  }

  void set_execution_space(const stk::ngp::ExecSpace& executionSpace) const {
    m_execSpace = executionSpace;
  }

  void set_execution_space(stk::ngp::ExecSpace&& executionSpace) const {
    m_execSpace = std::move(executionSpace);
  }

  void reset_execution_space() const {
    m_execSpace = Kokkos::DefaultExecutionSpace();
  }

  CSet & get_attributes() { return m_attribute; }

  template<class A>
    const A * declare_attribute_no_delete(const A * a) {
      return m_attribute.template insert_no_delete<A>(a);
    }   

  template<class A>
    const A * declare_attribute_with_delete(const A * a) {
      return m_attribute.template insert_with_delete<A>(a);
    }   

  template<class A>
    bool remove_attribute(const A * a) {
      return m_attribute.template remove<A>(a);
    }

  template<typename FieldType>
  void set_field_states(FieldType ** field_states)
  {
    for (unsigned i = 0; i < m_num_states; ++i) {
      m_field_states[i] = field_states[i];
    }
  }


  //  Associate this field with a bulk data.
  //  Note, a field can be associated with one and only one bulk data object
  void set_mesh(stk::mesh::BulkData* bulk);

  NgpFieldBase * get_ngp_field() const { return m_ngpField; }
  void set_ngp_field(NgpFieldBase * ngpField) const;

  void increment_num_syncs_to_host() const;
  void increment_num_syncs_to_device() const;

  void set_initial_value(const void* new_initial_value, unsigned num_scalars, unsigned num_bytes);

  void insert_restriction(const char     * arg_method ,
                          const Part       & arg_part ,
                          const unsigned     arg_num_scalars_per_entity ,
                          const unsigned     arg_first_dimension ,
                          const void*        arg_init_value = nullptr);

  void insert_restriction(const char     * arg_method ,
                          const Selector   & arg_selector ,
                          const unsigned     arg_num_scalars_per_entity ,
                          const unsigned     arg_first_dimension ,
                          const void*        arg_init_value = nullptr);

  void verify_and_clean_restrictions( const Part& superset, const Part& subset );

  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::BulkData;
  friend class ::stk::mesh::impl::FieldRepository;
  friend CSet & impl::get_attributes(stk::mesh::FieldBase & field);

  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestFieldImpl ;

  friend NgpFieldBase* impl::get_ngp_field(const FieldBase & stkField);
  friend void impl::set_ngp_field(const FieldBase & stkField, NgpFieldBase * ngpField);

  template <typename T, template <typename> class NgpDebugger> friend class HostField;
  template <typename T, template <typename> class NgpDebugger> friend class DeviceField;
  template <typename Scalar> friend class Field;

protected:
  FieldBase(MetaData                   * arg_mesh_meta_data,
            stk::topology::rank_t        arg_entity_rank,
            unsigned                     arg_ordinal,
            const std::string          & arg_name,
            const DataTraits           & arg_traits,
            unsigned                     arg_number_of_states,
            FieldState                   arg_this_state)
    : m_mesh(nullptr),
      m_entity_rank(arg_entity_rank),
      m_name(arg_name),
      m_num_states(arg_number_of_states),
      m_restrictions(),
      m_initial_value(nullptr),
      m_initial_value_num_bytes(0),
      m_data_traits( arg_traits ),
      m_meta_data(arg_mesh_meta_data),
      m_ordinal(arg_ordinal),
      m_this_state(arg_this_state),
      m_ngpField(nullptr),
      m_numSyncsToHost(0),
      m_numSyncsToDevice(0),
      m_modifiedOnHost(false),
      m_modifiedOnDevice(false),
      m_execSpace(Kokkos::DefaultExecutionSpace())
  {
    FieldBase * const pzero = nullptr ;
    Copy<MaximumFieldStates>(    m_field_states , pzero );
  }

private:
  FieldMetaDataVector m_field_meta_data;

  stk::mesh::BulkData* m_mesh;
  EntityRank m_entity_rank;
  const std::string m_name;
  const unsigned m_num_states;
  FieldBase                  * m_field_states[ MaximumFieldStates ];
  FieldRestrictionVector       m_restrictions;
  void*                        m_initial_value;
  unsigned                     m_initial_value_num_bytes;
  CSet                         m_attribute;
  const DataTraits           & m_data_traits;
  MetaData             * const m_meta_data;
  const unsigned               m_ordinal;
  const FieldState             m_this_state;
  mutable NgpFieldBase       * m_ngpField;
  mutable size_t               m_numSyncsToHost;
  mutable size_t               m_numSyncsToDevice;
  mutable bool                 m_modifiedOnHost;
  mutable bool                 m_modifiedOnDevice;
  mutable stk::ngp::ExecSpace  m_execSpace;
  mutable std::any m_stkFieldSyncDebugger;
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
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  return f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity;
}

inline unsigned field_bytes_per_entity(const FieldBase& f, unsigned bucket_id) {
  STK_ThrowAssert(bucket_id < f.get_meta_data_for_field().size());
  return f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity;
}

inline unsigned field_bytes_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return field_bytes_per_entity(f, bulk.bucket(e));
}

inline bool is_matching_rank(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  return(b.entity_rank() == f.entity_rank());
}

inline bool is_matching_rank(const FieldBase& f, Entity e) {
  return is_matching_rank(f, f.get_mesh().bucket(e));
}

inline bool is_matching_rank(const FieldBase& f, EntityRank rank) {
  return f.entity_rank() == rank;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  return f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity/bytes_per_scalar;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, unsigned bucket_id) {
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  return f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity/bytes_per_scalar;
}

inline unsigned field_scalars_per_entity(const FieldBase& f, Entity e) {
  const unsigned bytes_per_scalar = f.data_traits().size_of;
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return field_bytes_per_entity(f, bulk.bucket(e))/bytes_per_scalar;
}


inline unsigned field_extent0_per_entity(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  return f.get_meta_data_for_field()[b.bucket_id()].m_firstDimension;
}

inline unsigned field_extent0_per_entity(const FieldBase& f, unsigned bucket_id) {
  return f.get_meta_data_for_field()[bucket_id].m_firstDimension;
}

inline unsigned field_extent0_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_firstDimension;
}


inline unsigned field_extent1_per_entity(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  const unsigned numBytesPerScalar = f.data_traits().size_of;
  const unsigned numBytesPerEntity = f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity;
  const unsigned firstDimension = f.get_meta_data_for_field()[b.bucket_id()].m_firstDimension;
  return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
}

inline unsigned field_extent1_per_entity(const FieldBase& f, unsigned bucket_id) {
  const unsigned numBytesPerScalar = f.data_traits().size_of;
  const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity;
  const unsigned firstDimension = f.get_meta_data_for_field()[bucket_id].m_firstDimension;
  return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
}

inline unsigned field_extent1_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  const unsigned numBytesPerScalar = f.data_traits().size_of;
  const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_bytesPerEntity;
  const unsigned firstDimension = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_firstDimension;
  return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
}


inline unsigned field_extent_per_entity(const FieldBase& f, unsigned dimension, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  if (dimension == 0) {
    return f.get_meta_data_for_field()[b.bucket_id()].m_firstDimension;
  }
  else if (dimension == 1) {
    const unsigned numBytesPerScalar = f.data_traits().size_of;
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity;
    const unsigned firstDimension = f.get_meta_data_for_field()[b.bucket_id()].m_firstDimension;
    return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
  }
  else {
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity;
    return (numBytesPerEntity != 0) ? 1 : 0;
  }
}

inline unsigned field_extent_per_entity(const FieldBase& f, unsigned dimension, unsigned bucket_id) {
  if (dimension == 0) {
    return f.get_meta_data_for_field()[bucket_id].m_firstDimension;
  }
  else if (dimension == 1) {
    const unsigned numBytesPerScalar = f.data_traits().size_of;
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity;
    const unsigned firstDimension = f.get_meta_data_for_field()[bucket_id].m_firstDimension;
    return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
  }
  else {
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity;
    return (numBytesPerEntity != 0) ? 1 : 0;
  }
}

inline unsigned field_extent_per_entity(const FieldBase& f, unsigned dimension, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  if (dimension == 0) {
    return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_firstDimension;
  }
  else if (dimension == 1) {
    const unsigned numBytesPerScalar = f.data_traits().size_of;
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_bytesPerEntity;
    const unsigned firstDimension = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_firstDimension;
    return (numBytesPerEntity != 0) ? numBytesPerEntity/(numBytesPerScalar*firstDimension) : 0;
  }
  else {
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_bytesPerEntity;
    return (numBytesPerEntity != 0) ? 1 : 0;
  }
}


inline bool field_is_allocated_for_bucket(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(&b.mesh() == &f.get_mesh());
  //return true if field and bucket have the same rank and the field is associated with the bucket
  STK_ThrowAssert(f.get_meta_data_for_field().size() > b.bucket_id());
  return (is_matching_rank(f, b) && 0 != f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity);
}

inline size_t get_total_ngp_field_allocation_bytes(const FieldBase & f) {
  const Selector selector = stk::mesh::selectField(f);
  const BucketVector & buckets = f.get_mesh().get_buckets(f.entity_rank(), selector);
  const size_t numBuckets = buckets.size();
  const size_t bucketCapacity = (buckets.empty()) ? 0 : f.get_mesh().get_maximum_bucket_capacity();
  const size_t numPerEntity = f.max_size();
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
void debug_stale_access_entity_check(const FieldBase& stkField, const Entity& entity, const char* fileName, int lineNumber)
{
  if constexpr (std::is_same_v<StkDebugger, StkFieldSyncDebugger>)
  {
    if (stkField.has_ngp_field()) {
      stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_entity_check(entity, fileName, lineNumber);
    }
  }
}

template <typename StkDebugger>
void debug_stale_access_entity_check(const FieldBase& stkField, const unsigned bucketId, unsigned bucketOrd, const char* fileName, int lineNumber)
{
  if constexpr (std::is_same_v<StkDebugger, StkFieldSyncDebugger>)
  {
    if (stkField.has_ngp_field()) {
      stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_entity_check(bucketId, bucketOrd, fileName, lineNumber);
    }
  }
}

template <typename StkDebugger>
void debug_stale_access_bucket_check(const FieldBase& stkField, const Bucket& bucket, const char* fileName, int lineNumber)
{
  if constexpr (std::is_same_v<StkDebugger, StkFieldSyncDebugger>)
  {
    if (stkField.has_ngp_field()) {
      stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_bucket_check(bucket, fileName, lineNumber);
    }
  }
}

template <typename StkDebugger>
void debug_stale_access_bucket_check(const FieldBase& stkField, const unsigned& bucketId, const char* fileName, int lineNumber)
{
  if constexpr (std::is_same_v<StkDebugger, StkFieldSyncDebugger>)
  {
    if (stkField.has_ngp_field()) {
      stkField.get_field_sync_debugger<StkDebugger>().host_stale_access_bucket_check(bucketId, fileName, lineNumber);
    }
  }
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id, unsigned bucket_ord, const int knownSize,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  STK_ThrowAssertMsg(f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity >= knownSize,
                 "field name= " << f.name() << "knownSize= " << knownSize << ", m_bytesPerEntity= "
                 << f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity);
  STK_ThrowAssert(f.get_meta_data_for_field()[bucket_id].m_data != nullptr);

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, bucket_ord, fileName, lineNumber);

  return reinterpret_cast<typename FieldType::value_type*>(f.get_meta_data_for_field()[bucket_id].m_data +
                                                           knownSize * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  debug_stale_access_bucket_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, fileName, lineNumber);

  return reinterpret_cast<typename FieldType::value_type*>(f.get_meta_data_for_field()[bucket_id].m_data);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id, unsigned bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), bucket_id, bucket_ord, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[bucket_id];
  return reinterpret_cast<typename FieldType::value_type*>(field_meta_data.m_data +
                                                           field_meta_data.m_bytesPerEntity * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const Bucket& b, unsigned bucket_ord,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), b[bucket_ord], fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(field_meta_data.m_data +
                                                           field_meta_data.m_bytesPerEntity * bucket_ord);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const Bucket& b,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&b.mesh() == &f.get_mesh());

  debug_stale_access_bucket_check<StkDebugger>(static_cast<const FieldBase&>(f), b, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(field_meta_data.m_data);
}

template<class FieldType, typename StkDebugger = DefaultStkFieldSyncDebugger>
inline
typename FieldType::value_type*
field_data(const FieldType & f, Entity e,
           DummyOverload dummyArg = DummyOverload(), const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER)
{
  const MeshIndex& mi = f.get_mesh().mesh_index(e);
  STK_ThrowAssertMsg(f.entity_rank() == mi.bucket->entity_rank(),
                 "field_data called with " << f.entity_rank() << " field (" << f.name() << ") and different-rank entity "
                 << f.get_mesh().entity_key(e) << ". The rank of the field and entity must match.");
  STK_ThrowAssert(&f.get_mesh() == &mi.bucket->mesh());

  debug_stale_access_entity_check<StkDebugger>(static_cast<const FieldBase&>(f), e, fileName, lineNumber);

  const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[mi.bucket->bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(field_meta_data.m_data +
                                                           field_meta_data.m_bytesPerEntity * mi.bucket_ordinal);
}

} //namespace stk::mesh

#endif //stk_mesh_base_FieldBase_hpp
