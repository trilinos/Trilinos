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

#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldRestriction.hpp>  // for FieldRestriction, etc
#include <stk_mesh/base/FieldState.hpp>  // for FieldState
#include <stk_mesh/base/Types.hpp>      // for EntityRank, etc
#include <stk_mesh/base/NgpFieldBase.hpp>
#include <stk_mesh/base/ConstFieldData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/ConstFieldBytes.hpp>
#include <stk_mesh/base/FieldBytes.hpp>
#include <stk_mesh/base/FieldDataBase.hpp>
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include <stk_util/util/SimpleArrayOps.hpp>  // for Copy
#include <stk_util/util/CSet.hpp>
#include <string>                       // for string
#include <type_traits>

namespace stk::mesh 
{

class BulkData;
class MetaData;
class UnitTestFieldImpl;
class FieldBase;
template<typename T, typename NgpMemSpace> class DeviceField;

template <typename T, typename NgpMemSpace = NgpMeshDefaultMemSpace>
NgpField<T, NgpMemSpace>& get_updated_ngp_field_async(const FieldBase&, const stk::ngp::ExecSpace&);

template <typename T, typename NgpMemSpace = NgpMeshDefaultMemSpace>
NgpField<T, NgpMemSpace>& get_updated_ngp_field_async(const FieldBase&, stk::ngp::ExecSpace&&);

namespace impl {
class FieldRepository;
FieldDataBase* get_device_data(const FieldBase& field);
NgpFieldBase* get_ngp_field(const FieldBase& field);
void set_ngp_field(const FieldBase& stkField, NgpFieldBase* ngpField);
stk::CSet & get_attributes(FieldBase& field);
}

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

  Layout host_data_layout() const { return m_hostDataLayout; }
  Layout device_data_layout() const { return m_deviceDataLayout; }

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
    delete m_hostFieldData;
    delete m_deviceFieldData;
    delete m_ngpField;

    if (state() == StateNone) {
      void*& init_val = m_initial_value;

      delete [] reinterpret_cast<char*>(init_val);
      init_val = nullptr;
    }
  }

  virtual std::ostream& print_data(std::ostream& out, void* /*data*/, unsigned /*size_per_entity*/) const { return out; }

  stk::mesh::BulkData& get_mesh() const
  { return *m_mesh; }

  // Access the main host-side FieldMetaData array
  inline const FieldMetaDataArrayType& get_internal_field_meta_data() const {
    return static_cast<FieldBytes<stk::ngp::HostMemSpace>&>(*m_hostFieldData).m_fieldMetaData;
  }

  inline FieldMetaDataArrayType& get_internal_field_meta_data() {
    auto& hostFieldBytes = static_cast<FieldBytes<stk::ngp::HostMemSpace>&>(*m_hostFieldData);
    hostFieldBytes.modify_field_meta_data();

    return hostFieldBytes.m_fieldMetaData;
  }

  // Access the cached copy of the host-side FieldMetaData array
  inline const FieldMetaDataArrayType& get_meta_data_for_field() const {
    return m_cachedFieldMetaData;
  }

  inline FieldMetaDataArrayType& get_meta_data_for_field() {
    return m_cachedFieldMetaData;
  }

  // Access the cached copy of the host-side FieldMetaData array
  inline void update_cached_field_meta_data() {
    m_cachedFieldMetaData = static_cast<FieldBytes<stk::ngp::HostMemSpace>&>(*m_hostFieldData).m_fieldMetaData;
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
  
  bool defined_on(const stk::mesh::Entity& entity) const;

  void modify_on_host() const;
  void modify_on_device() const;
  void modify_on_host(const Selector& s) const;
  void modify_on_device(const Selector& s) const;
  bool need_sync_to_host() const;
  bool need_sync_to_device() const;
  void sync_to_host() const;
  void sync_to_host(const stk::ngp::ExecSpace& execSpace) const;
  void sync_to_device() const;
  void sync_to_device(const stk::ngp::ExecSpace& execSpace) const;
  void clear_sync_state() const;
  void clear_host_sync_state() const;
  void clear_device_sync_state() const;
  void fence() const;
  void fence(const stk::ngp::ExecSpace& execSpace) const;

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
  bool has_device_data() const { return m_deviceFieldData != nullptr; }

  void rotate_multistate_data(bool rotateNgpFieldViews = false);

protected:
  void update_host_field_data() const
  {
    if (m_hostFieldData->needs_update()) {
      m_hostFieldData->update(stk::ngp::ExecSpace(), host_data_layout());
    }
  }

  template <typename T, typename MemSpace, Layout DataLayout>
  void update_or_create_device_field_data() const
  {
    static_assert(Kokkos::SpaceAccessibility<stk::ngp::ExecSpace, MemSpace>::accessible);

    if (m_deviceFieldData != nullptr) {
      if (m_deviceFieldData->needs_update()) {
        m_deviceFieldData->update(stk::ngp::ExecSpace(), host_data_layout());
        increment_num_syncs_to_device();
      }
    }
    else {
      m_deviceFieldData =
          new FieldData<T, MemSpace, DataLayout>(static_cast<FieldBytes<stk::ngp::HostMemSpace>*>(m_hostFieldData));

      m_deviceFieldData->set_mesh(&get_mesh());
      if (m_deviceFieldData->needs_update()) {
        m_deviceFieldData->update(stk::ngp::ExecSpace(), host_data_layout());
        increment_num_syncs_to_device();
      }
      clear_host_sync_state();
    }
  }

  template <typename T, typename MemSpace, Layout DataLayout, typename ExecSpace>
  void async_update_or_create_device_field_data(const ExecSpace& execSpace) const
  {
    static_assert(Kokkos::SpaceAccessibility<stk::ngp::ExecSpace, MemSpace>::accessible);

    if (m_deviceFieldData != nullptr) {
      if (m_deviceFieldData->needs_update()) {
        m_deviceFieldData->update(execSpace, host_data_layout());
        increment_num_syncs_to_device();
      }
    }
    else {
      m_deviceFieldData =
          new FieldData<T, MemSpace, DataLayout>(static_cast<FieldBytes<stk::ngp::HostMemSpace>*>(m_hostFieldData));

      m_deviceFieldData->set_mesh(&get_mesh());
      if (m_deviceFieldData->needs_update()) {
        m_deviceFieldData->update(execSpace, host_data_layout());
        increment_num_syncs_to_device();
      }
      clear_host_sync_state();
    }
  }

  template <typename FieldDataType>
  void check_field_data_cast(FieldDataType fieldDataPtr) const {
    STK_ThrowRequireMsg(fieldDataPtr != nullptr, "Called FieldBase::data() with inconsistent template parameters for "
                        "Field '" << name() << "'.  The datatype and layout template parameters must match the "
                        "original Field registration.");
  }

  template <typename T, typename MemSpace, Layout DataLayout>
  FieldData<T, MemSpace, DataLayout>& field_data_handle(const FieldBase& fieldBase) const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      update_host_field_data();
      auto ptr = dynamic_cast<FieldData<T, stk::ngp::HostMemSpace, DataLayout>*>(fieldBase.m_hostFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
    else {
      update_or_create_device_field_data<T, MemSpace, DataLayout>();
      auto ptr = dynamic_cast<FieldData<T, MemSpace, DataLayout>*>(fieldBase.m_deviceFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
  }

  template <typename T, typename MemSpace, Layout DataLayout>
  ConstFieldData<T, MemSpace, DataLayout>& const_field_data_handle(const FieldBase& fieldBase) const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      update_host_field_data();
      auto ptr = dynamic_cast<ConstFieldData<T, stk::ngp::HostMemSpace, DataLayout>*>(fieldBase.m_hostFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
    else {
      update_or_create_device_field_data<T, MemSpace, DataLayout>();
      auto ptr = dynamic_cast<ConstFieldData<T, MemSpace, DataLayout>*>(fieldBase.m_deviceFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
  }


  template <typename T, typename MemSpace, Layout DataLayout, typename ExecSpace>
  FieldData<T, MemSpace, DataLayout>& async_field_data_handle(const FieldBase& fieldBase,
                                                              const ExecSpace& execSpace) const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      update_host_field_data();
      auto ptr = dynamic_cast<FieldData<T, stk::ngp::HostMemSpace, DataLayout>*>(fieldBase.m_hostFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
    else {
      async_update_or_create_device_field_data<T, MemSpace, DataLayout, ExecSpace>(execSpace);
      auto ptr = dynamic_cast<FieldData<T, MemSpace, DataLayout>*>(fieldBase.m_deviceFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
  }

  template <typename T, typename MemSpace, Layout DataLayout, typename ExecSpace>
  ConstFieldData<T, MemSpace, DataLayout>& async_const_field_data_handle(const FieldBase& fieldBase,
                                                                         const ExecSpace& execSpace) const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      update_host_field_data();
      auto ptr = dynamic_cast<ConstFieldData<T, stk::ngp::HostMemSpace, DataLayout>*>(fieldBase.m_hostFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
    else {
      async_update_or_create_device_field_data<T, MemSpace, DataLayout, ExecSpace>(execSpace);
      auto ptr = dynamic_cast<ConstFieldData<T, MemSpace, DataLayout>*>(fieldBase.m_deviceFieldData);
      check_field_data_cast(ptr);
      return *ptr;
    }
  }


  // Can't specialize a class member function -- only a class or free function
  template <typename T, FieldAccessTag FieldAccess, typename MemSpace, Layout DataLayout>
  struct FieldDataHelper {
    using FieldDataType = FieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase) {
      return fieldBase.field_data_handle<T, MemSpace, DataLayout>(fieldBase);
    }
  };

  template <typename T, typename MemSpace, Layout DataLayout>
  struct FieldDataHelper<T, ReadOnly, MemSpace, DataLayout> {
    using FieldDataType = ConstFieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase) {
      return fieldBase.const_field_data_handle<T, MemSpace, DataLayout>(fieldBase);
    }
  };

  template <typename T, typename MemSpace, Layout DataLayout>
  struct FieldDataHelper<T, ConstUnsynchronized, MemSpace, DataLayout> {
    using FieldDataType = ConstFieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase) {
      return fieldBase.const_field_data_handle<T, MemSpace, DataLayout>(fieldBase);
    }
  };

  template <typename T, FieldAccessTag FieldAccess, typename MemSpace>
  struct AutoLayoutFieldDataBuilder {
    using FieldDataType = typename FieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto>::FieldDataType;
    static FieldDataType build_field_data(const FieldBase& fieldBase) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        // Build from a properly-cast FieldData so that we can update it properly first while still taking
        // advantage of copy elision by directly returning the newly-constructed FieldData object.
        if (fieldBase.host_data_layout() == Layout::Right) {
          auto fieldDataRight = FieldDataHelper<T, FieldAccess, MemSpace, Layout::Right>::field_data(fieldBase);
          return typename FieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto>::FieldDataType(fieldDataRight,
                                                                                                 FieldAccess);
        }
        else if (fieldBase.host_data_layout() == Layout::Left) {
          auto fieldDataLeft = FieldDataHelper<T, FieldAccess, MemSpace, Layout::Left>::field_data(fieldBase);
          return typename FieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto>::FieldDataType(fieldDataLeft,
                                                                                                 FieldAccess);
        }
        else {
          STK_ThrowErrorMsg("Host data layout of " << fieldBase.host_data_layout() << " is unsupported.  It must be "
                            "either Layout::Right or Layout::Left.");
          return typename FieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto>::FieldDataType(); // Keep compiler happy
        }
      }
      else {
        STK_ThrowErrorMsg("Layout::Auto access to Field data is only available on the host.");
        return typename FieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto>::FieldDataType(); // Keep compiler happy
      }
    }
  };


  template <typename T, FieldAccessTag FieldAccess, typename MemSpace, Layout DataLayout, typename ExecSpace>
  struct AsyncFieldDataHelper {
    using FieldDataType = FieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase, const ExecSpace& execSpace) {
      return fieldBase.async_field_data_handle<T, MemSpace, DataLayout, ExecSpace>(fieldBase, execSpace);
    }
  };

  template <typename T, typename MemSpace, Layout DataLayout, typename ExecSpace>
  struct AsyncFieldDataHelper<T, ReadOnly, MemSpace, DataLayout, ExecSpace> {
    using FieldDataType = ConstFieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase, const ExecSpace& execSpace) {
      return fieldBase.async_const_field_data_handle<T, MemSpace, DataLayout, ExecSpace>(fieldBase, execSpace);
    }
  };

  template <typename T, typename MemSpace, Layout DataLayout, typename ExecSpace>
  struct AsyncFieldDataHelper<T, ConstUnsynchronized, MemSpace, DataLayout, ExecSpace> {
    using FieldDataType = ConstFieldData<T, MemSpace, DataLayout>;
    static FieldDataType& field_data(const FieldBase& fieldBase, const ExecSpace& execSpace) {
      return fieldBase.async_const_field_data_handle<T, MemSpace, DataLayout, ExecSpace>(fieldBase, execSpace);
    }
  };

  template <typename T, FieldAccessTag FieldAccess, typename MemSpace, typename ExecSpace>
  struct AsyncAutoLayoutFieldDataBuilder {
    using FieldDataType =
        typename AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto, ExecSpace>::FieldDataType;
    static FieldDataType build_field_data(const FieldBase& fieldBase, const ExecSpace& execSpace) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        // Build from a properly-cast FieldData so that we can update it properly first while still taking
        // advantage of copy elision by directly returning the newly-constructed FieldData object.
        if (fieldBase.host_data_layout() == Layout::Right) {
          auto fieldDataRight =
              AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Right, ExecSpace>::field_data(fieldBase,
                                                                                                   execSpace);
          return typename
              AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto, ExecSpace>::FieldDataType(fieldDataRight,
                                                                                                     FieldAccess);
        }
        else if (fieldBase.host_data_layout() == Layout::Left) {
          auto fieldDataLeft =
              AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Left, ExecSpace>::field_data(fieldBase,
                                                                                                  execSpace);
          return typename
              AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto, ExecSpace>::FieldDataType(fieldDataLeft,
                                                                                                     FieldAccess);
        }
        else {
          STK_ThrowErrorMsg("Host data layout of " << fieldBase.host_data_layout() << " is unsupported.  It must be "
                            "either Layout::Right or Layout::Left.");
          return typename AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto, ExecSpace>::FieldDataType(); // Keep compiler happy
        }
      }
      else {
        STK_ThrowErrorMsg("Layout::Auto access to Field data is only available on the host.");
        return typename AsyncFieldDataHelper<T, FieldAccess, MemSpace, Layout::Auto, ExecSpace>::FieldDataType(); // Keep compiler happy
      }
    }
  };

  void print_host_access_error(FieldAccessTag hostAccessTag, FieldAccessTag deviceAccessTag) const
  {
    STK_ThrowErrorMsg("Trying to access host-side FieldData with access tag " << hostAccessTag <<
                      " while device-side FieldData with access tag " << deviceAccessTag <<
                      " is still in use for Field '" << name() << "'.");
  }

  void print_device_access_error(FieldAccessTag deviceAccessTag, FieldAccessTag hostAccessTag) const
  {
    STK_ThrowErrorMsg("Trying to access device-side FieldData with access tag " << deviceAccessTag <<
                      " while host-side FieldData with access tag " << hostAccessTag <<
                      " is still in use for Field '" << name() << "'.");
  }

  template <FieldAccessTag FieldAccess, typename MemSpace>
  void check_lifetimes() const
  {
#ifdef STK_USE_DEVICE_MESH
    if constexpr (FieldAccess == ReadOnly) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        if (has_device_data()) {
          if (m_deviceFieldData->has_copies(ReadWrite) || m_deviceFieldData->has_copies(OverwriteAll)) {
            print_host_access_error(FieldAccess, m_deviceFieldData->access_tag());
          }
        }
      }
      else {
        if (m_hostFieldData->has_copies(ReadWrite) || m_hostFieldData->has_copies(OverwriteAll)) {
          print_device_access_error(FieldAccess, m_hostFieldData->access_tag());
        }
      }
    }
    if constexpr (FieldAccess == ReadWrite || FieldAccess == OverwriteAll) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        if (has_device_data()) {
          if (m_deviceFieldData->has_copies(ReadWrite) || m_deviceFieldData->has_copies(ReadOnly) ||
              m_deviceFieldData->has_copies(OverwriteAll)) {
            print_host_access_error(FieldAccess, m_deviceFieldData->access_tag());
          }
        }
      }
      else {
        if (m_hostFieldData->has_copies(ReadWrite) || m_hostFieldData->has_copies(ReadOnly) ||
            m_hostFieldData->has_copies(OverwriteAll)) {
          print_device_access_error(FieldAccess, m_hostFieldData->access_tag());
        }
      }
    }
#endif
  }

public:

  // This function is intended to be used when you have a persistent copy of a FieldData object
  // that was acquired with the stk::mesh::Unsynchronized or stk::mesh::ConstUnsynchronized
  // access tag.  The unsynchronized variants will not do any automatic data movement between
  // memory spaces, so this synchronize() function is your opportunity to force the same
  // operations.  Call this right before using your persistent FieldData copy, with the
  // FieldAccessTag and MemorySpace matching what you are about to do with your data.  Your
  // persistent copy will then be updated appropriately before use.
  //
  // This function is completely unnecessary for the normal workflow of reacquiring your
  // FieldData instance with a FieldBase::data() call before each use.

  template <FieldAccessTag FieldAccess = ReadWrite,
            typename MemSpace = stk::ngp::HostMemSpace>
  void synchronize() const
  {
#ifdef STK_USE_DEVICE_MESH
    if constexpr (FieldAccess == ReadWrite || FieldAccess == ReadOnly) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        sync_to_host();
      }
      else {
        sync_to_device();
      }
    }

    if constexpr (FieldAccess == ReadWrite) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        modify_on_host();
      }
      else {
        modify_on_device();
      }
    }
    else if constexpr (FieldAccess == OverwriteAll) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        clear_device_sync_state();
        modify_on_host();
      }
      else {
        clear_host_sync_state();
        modify_on_device();
      }
    }
#endif
  }

  // This function is the same as described above, except that it takes a Kokkos execution space
  // argument that will be used to run any synchronization operations that may need to take
  // place.  This is intended for asynchronous execution.

  template <FieldAccessTag FieldAccess = ReadWrite,
            typename MemSpace = stk::ngp::HostMemSpace,
            typename ExecSpace = stk::ngp::ExecSpace>
  void synchronize([[maybe_unused]] const ExecSpace& execSpace) const
  {
#ifdef STK_USE_DEVICE_MESH
    if constexpr (FieldAccess == ReadWrite || FieldAccess == ReadOnly) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        sync_to_host(execSpace);
      }
      else {
        sync_to_device(execSpace);
      }
    }

    if constexpr (FieldAccess == ReadWrite) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        modify_on_host();
      }
      else {
        modify_on_device();
      }
    }
    else if constexpr (FieldAccess == OverwriteAll) {
      if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
        clear_device_sync_state();
        modify_on_host();
      }
      else {
        clear_host_sync_state();
        modify_on_device();
      }
    }
#endif
  }

  // Acquire a temporary copy of a FieldData object to access your data in the desired memory space.
  // It is best to acquire it right before using it in a computational algorithm, and then let it go
  // out of scope and be destroyed when you are finished with it.  The lifetime of this object cannot
  // overlap with one for the same Field in the opposite memory space.  Lifetimes can overlap in the
  // same memory space, but you must remember to always reacquire this object after a mesh modification.
  // Debug builds will trap these errors.  The four template parameters are:
  //
  //   T : Mandatory parameter with the Field datatype.  This must match the type that the Field
  //     was originally registered with, or this function will throw.  Use a templated Field
  //     instance to omit this parameter.
  //
  //   FieldAccessTag : Optional tag indicating how you will access the data.  Options are:

  //     - stk::mesh::ReadWrite     : Sync data to memory space and mark modified; Allow modification [default]
  //     - stk::mesh::ReadOnly      : Sync data to memory space and do not mark modified; Disallow modification
  //     - stk::mesh::OverwriteAll  : Do not sync data and mark modified; Allow modification

  //     - stk::mesh::Unsynchronized       : Do not sync data and do not mark modified; Allow modification
  //     - stk::mesh::ConstUnsynchronized  : Do not sync data and do not mark modified; Disallow modification
  //
  //     This will control automatic synchronization between memory spaces so that you are guaranteed
  //     to be using up-to-date data wherever accessed.  The Unsynchronized variants are intended for
  //     users who must store persistent copies of their FieldData objects.  You must call the above
  //     synchronize() function before running your algorithm that uses your persistent copy, with the
  //     same access tag and memory space that you would otherwise have used, to get the data movement
  //     correct.  Do not use the Unsynchronized access tags for normal workflows.
  //
  //   MemSpace : Optional Kokkos memory space of the data that you want to access.  It can be either
  //     a Kokkos host space or a device space.  You can use the aliases "stk::ngp::HostMemSpace" and
  //     "stk::ngp::MemSpace" as convenient shortcuts.  The HostMemSpace alias is always the host space
  //     and the MemSpace alias is the default device space in a device build or the host space in
  //     a host build.  The default is "stk::ngp::HostMemSpace".
  //
  //   Layout : Optional data layout that must match the layout that the Field was originally
  //     registered with.  This can be values of stk::mesh::Layout::Left or stk::mesh::Layout::Right.
  //     If you did nothing special for field registration, this will always default to the correct
  //     value.  If you overrode the host-side data layout at Field registration, you need to
  //     manually provide the same value here.  This will always be Layout::Left on device and
  //     will default to Layout::Right on host in a "normal" build and Layout::Left on host in
  //     a STK_UNIFIED_MEMORY build.  You can use a templated Field instance to always omit this
  //     tag, regardless of customization.
  //
  // Some sample usage for a FieldBase instance of a Field<double>:
  //
  //   auto fieldData = myField.data<double>();                       <-- Read-write access to host data
  //   auto fieldData = myField.data<double, stk::mesh::ReadOnly>();  <-- Read-only access to host access
  //   auto fieldData = myField.data<double, stk::mesh::ReadWrite, stk::ngp::MemSpace>(); <-- Read-write access to device data
  //   auto fieldData = myField.data<double, stk::mesh::ReadOnly, stk::ngp::MemSpace>();  <-- Read-only access to device data

  template <typename T,
            FieldAccessTag FieldAccess = ReadWrite,
            typename MemSpace = stk::ngp::HostMemSpace,
            Layout DataLayout = DefaultLayoutSelector<MemSpace>::layout>
  typename std::enable_if_t<is_field_datatype_v<T>,
                            typename FieldDataHelper<T, FieldAccess, MemSpace, DataLayout>::FieldDataType>
  data() const
  {
    check_lifetimes<FieldAccess, MemSpace>();
    synchronize<FieldAccess, MemSpace>();

    if constexpr (DataLayout == Layout::Auto) {
      return AutoLayoutFieldDataBuilder<T, FieldAccess, MemSpace>::build_field_data(*this);
    }
    else {
      return typename FieldDataHelper<T, FieldAccess, MemSpace, DataLayout>::FieldDataType(
            FieldDataHelper<T, FieldAccess, MemSpace, DataLayout>::field_data(*this), FieldAccess);
    }
  }


  // This is the same as described above, except you can pass in a custom execution space argument
  // that will be used to run any data syncing or updating after a mesh modification that may
  // be necessary.  This is for asynchronous execution.

  template <typename T,
            FieldAccessTag FieldAccess = ReadWrite,
            typename MemSpace = stk::ngp::HostMemSpace,
            Layout DataLayout = DefaultLayoutSelector<MemSpace>::layout,
            typename ExecSpace = stk::ngp::ExecSpace>
  typename std::enable_if_t<is_field_datatype_v<T>,
                            typename FieldDataHelper<T, FieldAccess, MemSpace, DataLayout>::FieldDataType>
  data(const ExecSpace& execSpace) const
  {
    check_lifetimes<FieldAccess, MemSpace>();
    synchronize<FieldAccess, MemSpace>(execSpace);

    if constexpr (DataLayout == Layout::Auto) {
      return AsyncAutoLayoutFieldDataBuilder<T, FieldAccess, MemSpace, ExecSpace>::build_field_data(*this, execSpace);
    }
    else {
      return typename AsyncFieldDataHelper<T, FieldAccess, MemSpace, DataLayout, ExecSpace>::FieldDataType(
            AsyncFieldDataHelper<T, FieldAccess, MemSpace, DataLayout, ExecSpace>::field_data(*this, execSpace),
            FieldAccess);
    }
  }


  // Acquire a reference to a FieldBytes object in the desired memory space.  This is intended for
  // low-level access to the raw bytes behind the data when the datatype and layout of the Field
  // are not known.  There is no automatic data synchronization or marking performed, and very
  // few consistency checks and no bounds-checking performed.  This is intended for internal
  // STK Mesh use only.

  template <typename MemSpace = stk::ngp::HostMemSpace>
  FieldBytes<MemSpace>& bytes() const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      return static_cast<FieldBytes<stk::ngp::HostMemSpace>&>(*m_hostFieldData);
    }
    else {
      STK_ThrowRequireMsg(has_device_data(),
                          "You must fully-construct the device FieldData object through a FieldBase::data() call "
                          "before requesting a subset of it through FieldBase::bytes()");
      return dynamic_cast<FieldBytes<MemSpace>&>(*m_deviceFieldData);
    }
  }

  template <typename MemSpace = stk::ngp::HostMemSpace>
  ConstFieldBytes<MemSpace>& const_bytes() const
  {
    if constexpr (std::is_same_v<MemSpace, stk::ngp::HostMemSpace>) {
      return static_cast<ConstFieldBytes<stk::ngp::HostMemSpace>&>(*m_hostFieldData);
    }
    else {
      STK_ThrowRequireMsg(has_device_data(),
                          "You must fully-construct the device FieldData object through a FieldBase::data() call "
                          "before requesting a subset of it through FieldBase::bytes()");
      return dynamic_cast<ConstFieldBytes<MemSpace>&>(*m_deviceFieldData);
    }
  }

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

  FieldDataBase* get_device_data() const { return m_deviceFieldData; }
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

  friend FieldDataBase* impl::get_device_data(const FieldBase& stkField);
  friend NgpFieldBase* impl::get_ngp_field(const FieldBase& stkField);
  friend void impl::set_ngp_field(const FieldBase& stkField, NgpFieldBase* ngpField);

  template <typename T, typename NgpMemSpace> friend class HostField;
  template <typename T, typename NgpMemSpace> friend class DeviceField;
  template <typename Scalar, Layout HostLayout> friend class Field;
  template <typename T, typename NgpMemSpace> friend
      NgpField<T, NgpMemSpace>& get_updated_ngp_field_async(const FieldBase&, const stk::ngp::ExecSpace&);
  template <typename T, typename NgpMemSpace> friend
      NgpField<T, NgpMemSpace>& get_updated_ngp_field_async(const FieldBase&, stk::ngp::ExecSpace&&);

protected:
  FieldBase(MetaData* arg_mesh_meta_data,
            stk::topology::rank_t arg_entity_rank,
            unsigned arg_ordinal,
            const std::string& arg_name,
            const DataTraits& arg_traits,
            unsigned arg_number_of_states,
            FieldState arg_this_state,
            Layout hostDataLayout,
            Layout deviceDataLayout,
            FieldDataBase* hostFieldData)
    : m_mesh(nullptr),
      m_hostFieldData(hostFieldData),
      m_deviceFieldData(nullptr),
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
      m_hostDataLayout(hostDataLayout),
      m_deviceDataLayout(deviceDataLayout),
      m_execSpace(Kokkos::DefaultExecutionSpace())
  {
    FieldBase * const pzero = nullptr ;
    Copy<MaximumFieldStates>(m_field_states, pzero);
  }

private:
  stk::mesh::BulkData* m_mesh;

  FieldDataBase* m_hostFieldData;
  mutable FieldDataBase* m_deviceFieldData;
  FieldMetaDataArrayType m_cachedFieldMetaData;

  EntityRank m_entity_rank;
  const std::string m_name;
  const unsigned m_num_states;
  FieldBase* m_field_states[ MaximumFieldStates ];
  FieldRestrictionVector m_restrictions;
  void* m_initial_value;
  unsigned m_initial_value_num_bytes;
  CSet m_attribute;
  const DataTraits& m_data_traits;
  MetaData* const m_meta_data;
  const unsigned m_ordinal;
  const FieldState m_this_state;
  mutable NgpFieldBase* m_ngpField;
  mutable size_t m_numSyncsToHost;
  mutable size_t m_numSyncsToDevice;
  mutable bool m_modifiedOnHost;
  mutable bool m_modifiedOnDevice;
  Layout m_hostDataLayout;
  Layout m_deviceDataLayout;
  mutable stk::ngp::ExecSpace m_execSpace;
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
inline FieldDataBase* get_device_data(const FieldBase& stkField) {
  return stkField.get_device_data();
}

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
  STK_ThrowAssert(bucket_id < f.get_meta_data_for_field().extent(0));
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
  return f.get_meta_data_for_field()[b.bucket_id()].m_numComponentsPerEntity;
}

inline unsigned field_extent0_per_entity(const FieldBase& f, unsigned bucket_id) {
  return f.get_meta_data_for_field()[bucket_id].m_numComponentsPerEntity;
}

inline unsigned field_extent0_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_numComponentsPerEntity;
}


inline unsigned field_extent1_per_entity(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  return f.get_meta_data_for_field()[b.bucket_id()].m_numCopiesPerEntity;
}

inline unsigned field_extent1_per_entity(const FieldBase& f, unsigned bucket_id) {
  return f.get_meta_data_for_field()[bucket_id].m_numCopiesPerEntity;
}

inline unsigned field_extent1_per_entity(const FieldBase& f, Entity e) {
  BulkData& bulk(f.get_mesh());
  STK_ThrowAssert(f.entity_rank() == bulk.entity_rank(e));
  return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_numCopiesPerEntity;
}


inline unsigned field_extent_per_entity(const FieldBase& f, unsigned dimension, const Bucket& b) {
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());
  if (dimension == 0) {
    return f.get_meta_data_for_field()[b.bucket_id()].m_numComponentsPerEntity;
  }
  else if (dimension == 1) {
    return f.get_meta_data_for_field()[b.bucket_id()].m_numCopiesPerEntity;
  }
  else {
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[b.bucket_id()].m_bytesPerEntity;
    return (numBytesPerEntity != 0) ? 1 : 0;
  }
}

inline unsigned field_extent_per_entity(const FieldBase& f, unsigned dimension, unsigned bucket_id) {
  if (dimension == 0) {
    return f.get_meta_data_for_field()[bucket_id].m_numComponentsPerEntity;
  }
  else if (dimension == 1) {
    return f.get_meta_data_for_field()[bucket_id].m_numCopiesPerEntity;
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
    return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_numComponentsPerEntity;
  }
  else if (dimension == 1) {
    return f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_numCopiesPerEntity;
  }
  else {
    const unsigned numBytesPerEntity = f.get_meta_data_for_field()[bulk.bucket(e).bucket_id()].m_bytesPerEntity;
    return (numBytesPerEntity != 0) ? 1 : 0;
  }
}


inline bool field_is_allocated_for_bucket(const FieldBase& f, const Bucket& b) {
  STK_ThrowAssert(&b.mesh() == &f.get_mesh());
  //return true if field and bucket have the same rank and the field is associated with the bucket
  STK_ThrowAssert(f.get_meta_data_for_field().extent(0) > b.bucket_id());
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

inline bool FieldBase::defined_on(const stk::mesh::Entity& entity) const
{ 
  return field_bytes_per_entity(*this, entity) > 0;
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

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id, unsigned bucket_ord, const int knownSize)
{
  STK_ThrowAssertMsg(f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity >= knownSize,
                 "field name= " << f.name() << "knownSize= " << knownSize << ", m_bytesPerEntity= "
                 << f.get_meta_data_for_field()[bucket_id].m_bytesPerEntity);
  STK_ThrowAssert(f.get_meta_data_for_field()[bucket_id].m_data != nullptr);

  return reinterpret_cast<typename FieldType::value_type*>(f.get_meta_data_for_field()[bucket_id].m_data +
                                                           knownSize * bucket_ord);
}

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id)
{
  return reinterpret_cast<typename FieldType::value_type*>(f.get_meta_data_for_field()[bucket_id].m_data);
}

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const unsigned bucket_id, unsigned bucket_ord)
{
  const FieldMetaData& fieldMetaData = f.get_meta_data_for_field()[bucket_id];
  return reinterpret_cast<typename FieldType::value_type*>(fieldMetaData.m_data +
                                                           fieldMetaData.m_bytesPerEntity * bucket_ord);
}

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const Bucket& b, unsigned bucket_ord)
{
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&f.get_mesh() == &b.mesh());

  const FieldMetaData& fieldMetaData = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(fieldMetaData.m_data +
                                                           fieldMetaData.m_bytesPerEntity * bucket_ord);
}

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, const Bucket& b)
{
  STK_ThrowAssert(f.entity_rank() == b.entity_rank());
  STK_ThrowAssert(&b.mesh() == &f.get_mesh());

  const FieldMetaData& fieldMetaData = f.get_meta_data_for_field()[b.bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(fieldMetaData.m_data);
}

template<class FieldType>
inline
typename FieldType::value_type*
field_data(const FieldType & f, Entity e)
{
  const MeshIndex& mi = f.get_mesh().mesh_index(e);
  STK_ThrowAssertMsg(f.entity_rank() == mi.bucket->entity_rank(),
                 "field_data called with " << f.entity_rank() << " field (" << f.name() << ") and different-rank entity "
                 << f.get_mesh().entity_key(e) << ". The rank of the field and entity must match.");
  STK_ThrowAssert(&f.get_mesh() == &mi.bucket->mesh());

  const FieldMetaData& fieldMetaData = f.get_meta_data_for_field()[mi.bucket->bucket_id()];
  return reinterpret_cast<typename FieldType::value_type*>(fieldMetaData.m_data +
                                                           fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal);
}

} //namespace stk::mesh

#endif //stk_mesh_base_FieldBase_hpp
