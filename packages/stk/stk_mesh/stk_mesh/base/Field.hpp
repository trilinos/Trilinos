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

#ifndef stk_mesh_Field_hpp
#define stk_mesh_Field_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/baseImpl/FieldRepository.hpp>
#include <stk_util/util/string_case_compare.hpp>  // for equal_case
#include <iomanip>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \ingroup stk_mesh_module
 *  \brief  Field with defined data type and multi-dimensions (if any)
 *
 * Conceptually, a field describes computational data that is associated with
 * mesh entities. The data is not related to connectivity or mesh structure; rather,
 * it is used to track some property (like velocity) that is relevant to the
 * problem being solved.
 *
 * A field specification has three components
 *   - Field datatype: Defines the type of this field's data.
 *     - Examples:
 *       - A double field-type:
 *           stk::mesh::Field<double>
 *       - An integer field-type:
 *           stk::mesh::Field<int>
 *
 *   - Field declaration: Defines/creates a specific field of a specified type. The API for
 *     declaring fields is in MetaData. Declaration of a field requires a field datatype,
 *     field name, and number of states. The number of states defines the "memory" of the Field;
 *     if number-of-states is N, the last N values of the field are stored, and the user
 *     can advance/rotate the states by calling BulkData::update_field_data_states(). FieldState.hpp
 *     contains an enum that clients should use to refer to the various states.
 *     - Example:
 *       using CoordFieldType stk::mesh::Field<double> CoordFieldType;
 *       CoordFieldType& coord_field = meta.declare_field<double>("<name>", <num_states>);
 *
 *   - Field restrictions: Defines the set of entities that have a field and the memory layout of
 *     the data (e.g. a vector-3).  FieldRestrictions are applied to Parts and entities pick
 *     up the field by being placed in the associated part. Also, FieldRestrictions allow the same
 *     Field to be applied to different ranks of entities. E.g., you could apply a velocity Field
 *     to both Faces and Edges through two different FieldRestrictions for each of the different ranks.
 *     Also, FieldRestrictions can apply a nodal field to all nodes in an element ranked part. This
 *     depends on the concept of induced part membership since you can't put a node directly in an
 *     element ranked part, but any nodes reachable through downward relations from elements that are
 *     in the part will be induced into the part and therefore pick up the field. This allows you to
 *     specify special nodal fields for different element types.
 *
 *     The public API for declaring field restrictions is in MetaData.hpp;
 *     we strongly recommend using the put_field_on_mesh free functions to create your field restrictions. A
 *     field may have many restrictions, but the restrictions have to be compatible with each other;
 *     restrictions can overlap if an entity is a member of two parts A,B and you declare restrictions
 *     for a field for both A and B. If such a situation arises, the dimensionality of both restrictions
 *     must be identical.
 *
 *     If you care about field component subscripting for output to a file, then call
 *     stk::io::set_field_output_type() from IossBridge.hpp.  This will give you output subscripts of, say,
 *     [_x, _y, _z] for the "Vector_3D" type instead of the default [_1, _2, _3].
 *
 *     - Examples:
 *       - Put a scalar field on all nodes in a part
 *           using ScalarFieldType = stk::mesh::Field<double>;
 *           ScalarFieldType& field = meta.declare_field<double>("<name>", <num_states>);
 *           stk::mesh::put_field_on_mesh(field, stk::topology::NODE_RANK, <part>, <initialValue>);
 *
 *       - Put a vector (of size 3) of doubles on all nodes in a part
 *           using CoordFieldType = stk::mesh::Field<double>;
 *           CoordFieldType& field = meta.declare_field<double>("<name>", <num_states>);
 *           stk::mesh::put_field_on_mesh(field, stk::topology::NODE_RANK, <part>, 3, <initialValues>);
 *           stk::io::set_field_output_type(field, "Vector_3D");
 *
 *       - Put 2 copies of a vector (of size 3) of doubles on all nodes in a part
 *           using CoordFieldType = stk::mesh::Field<double>;
 *           CoordFieldType& field = meta.declare_field<double>("<name>", <num_states>);
 *           stk::mesh::put_field_on_mesh(field, stk::topology::NODE_RANK, <part>, 3, 2, <initialValues>);
 *           stk::io::set_field_output_type(field, "Vector_3D");
 *
 *       - Put a tensor (of size 3x3) of doubles on all nodes in a part
 *           using TensorFieldType = stk::mesh::Field<double>;
 *           TensorFieldType& field = meta.declare_field<double>("<name>", <num_states>);
 *           stk::mesh::put_field_on_mesh(field, stk::topology::NODE_RANK, <part>, 9, <initialValues>);
 *           stk::io::set_field_output_type(field, "Full_Tensor_36");
 *
 * Items of interest
 *   - Accessing field data: see FieldData.hpp
 *
 * Field-related API
 *   - stk_mesh/base/MetaData - Methods for...
 *     - declaring fields
 *     - declaring field restrictions
 *     - declaring field attributes
 *     - declaring field relations
 *     - querying fields by name
 *   - stk_mesh/base/BulkData - Methods for advancing field state
 *   - stk_mesh/base/FieldData - Defines API for accessing field data
 *   - stk_mesh/base/FieldBase - Defines member-function API for Field objects
 *   - stk_mesh/base/Field - Place to put Field template, not much API here
 *   - stk_mesh/base/FieldParallel - Free functions for some parallel operations
 *   - stk_mesh/base/FieldRestriction - Defines FieldRestriction class
 *   - stk_mesh/base/FieldState - Defines enums to refer to field states
 *   - Type-related in ArrayDim, DataTraits
 *
 * - TODO Describe relationship with Buckets
 */

template <typename T, Layout HostLayout>
class Field : public FieldBase {
public:
  using value_type = T;
  static constexpr Layout host_layout = HostLayout;
  static constexpr Layout device_layout = DefaultDeviceLayout;

  Field(MetaData* meta,
        stk::topology::rank_t entityRank,
        unsigned ordinal,
        unsigned rankedOrdinal,
        const std::string& name,
        const DataTraits& dataTraits,
        unsigned numberOfStates,
        FieldState thisState)
    : FieldBase(meta, entityRank, ordinal, rankedOrdinal, name, dataTraits, numberOfStates, thisState,
                host_layout, device_layout,
                new FieldData<T, stk::ngp::HostSpace, HostLayout>(entityRank, ordinal, name, dataTraits))
  {
  }

  /** \brief  Query this field for a given field state. */
  Field & field_of_state( FieldState input_state ) const {
#ifndef NDEBUG
    return dynamic_cast<Field &>( * FieldBase::field_state(input_state) );
#else // NDEBUG
    return static_cast<Field &>( * FieldBase::field_state(input_state) );
#endif
  }

  virtual ~Field() override = default;

  virtual std::ostream& print_data(std::ostream& out, const MeshIndex& mi) const override
  {
    // Don't trigger any syncs or updates, since this can be called for debugging inside mesh modification
    auto fieldData = data<ConstUnsynchronized>();
    auto fieldValues = fieldData.entity_values(mi);

    auto previousPrecision = out.precision();
    constexpr auto thisPrecision = std::numeric_limits<T>::digits10;

    out << "{";
    for (stk::mesh::CopyIdx copy : fieldValues.copies()) {
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        if constexpr (sizeof(T) == 1) {  // std::cout of char types insist on printing character and not numerical value
          if constexpr (std::is_signed_v<T>) {
            out << std::setprecision(thisPrecision) << static_cast<int>(fieldValues(copy, component)) << " ";
          }
          else {
            out << std::setprecision(thisPrecision) << static_cast<unsigned>(fieldValues(copy, component)) << " ";
          }
        }
        else {
          out << std::setprecision(thisPrecision) << fieldValues(copy, component) << " ";
        }
      }
    }
    out << "}";

    out << std::setprecision(previousPrecision);
    return out;
  }

  virtual FieldBase * clone(stk::mesh::impl::FieldRepository & fieldRepo, const std::string fieldName = "") const override
  {
    FieldBase * f[MaximumFieldStates] {nullptr};

    static const char* reserved_state_suffix[6] = {
      "_STKFS_OLD",
      "_STKFS_N",
      "_STKFS_NM1",
      "_STKFS_NM2",
      "_STKFS_NM3",
      "_STKFS_NM4"
    };

    std::string clonedFieldName = (fieldName == "") ? name() : fieldName;

    for (unsigned i = 0 ; i < 6 ; ++i) {
      const int len_name   = clonedFieldName.size();
      const int len_suffix = std::strlen(reserved_state_suffix[i]);
      const int offset     = len_name - len_suffix ;
      if ( 0 <= offset ) {
        [[maybe_unused]] const char * const name_suffix = clonedFieldName.c_str() + offset;
        STK_ThrowErrorMsgIf(equal_case(name_suffix , reserved_state_suffix[i]),
                        "For name = \"" << name_suffix << "\" CANNOT HAVE THE RESERVED STATE SUFFIX \"" <<
                        reserved_state_suffix[i] << "\"");
      }
    }

    std::string fieldNames[MaximumFieldStates];

    fieldNames[0] = clonedFieldName;

    if (number_of_states() == 2) {
      fieldNames[1] = clonedFieldName;
      fieldNames[1].append(reserved_state_suffix[0]);
    }
    else {
      for (unsigned i = 1; i < number_of_states(); ++i) {
        fieldNames[i] = clonedFieldName;
        fieldNames[i].append(reserved_state_suffix[i]);
      }
    }

    for (unsigned i = 0; i < number_of_states(); ++i) {
      f[i] = new Field(&fieldRepo.mesh_meta_data(),
                       entity_rank(),
                       fieldRepo.get_fields().size(),
                       fieldRepo.get_fields(entity_rank()).size(),
                       fieldNames[i],
                       data_traits(),
                       number_of_states(),
                       static_cast<FieldState>(i));

      fieldRepo.add_field(f[i]);
    }

    const InitValsViewType& initVals = get_initial_value_bytes();
    f[0]->m_initial_value = FieldBase::InitValsViewType("Init-Vals-"+name(), initVals.extent(0));

    for(unsigned j=0; j<initVals.extent(0); ++j) {
      f[0]->m_initial_value(j) = initVals(j);
    }

    for (unsigned i = 0 ; i < number_of_states() ; ++i) {
      f[i]->set_field_states(f);
    }

    return f[0];
  }

  template <typename Space, typename Enable = void>
  struct LayoutSelector {
    static constexpr Layout layout = device_layout;
  };

  template <typename Enable>
  struct LayoutSelector<stk::ngp::HostSpace, Enable> {
    static constexpr Layout layout = host_layout;
  };


  // Acquire a temporary copy of a FieldData object to access your data in the desired memory space.
  // It is best to acquire it right before using it in a computational algorithm, and then let it go
  // out of scope and be destroyed when you are finished with it.  The lifetime of this object cannot
  // overlap with one for the same Field in the opposite memory space.  Lifetimes can overlap in the
  // same memory space, but you must remember to always reacquire this object after a mesh modification.
  // Debug builds will trap these errors.  The two template parameters are:
  //
  //   FieldAccessTag : Optional tag indicating how you will access the data.  Options are:

  //     - stk::mesh::ReadOnly      : Sync data to memory space and do not mark modified; Disallow modification [default]
  //     - stk::mesh::ReadWrite     : Sync data to memory space and mark modified; Allow modification
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
  //   Space : Optional struct that defines the Kokkos memory space and execution space that you want
  //     to access.  It can be either a stk::ngp::HostSpace or an arbitrary device-space struct of your
  //     choosing.  STK provides a stk::ngp::DeviceSpace struct for you if you want to use it on device.
  //     The default template parameter is stk::ngp::HostSpace.
  //
  // Some sample usage for a Field<double> instance:
  //
  //   auto fieldData = myField.data();                       <-- Read-only access to host data
  //   auto fieldData = myField.data<stk::mesh::ReadWrite>(); <-- Read-write access to host data
  //   auto fieldData = myField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  <-- Read-only access to device data
  //   auto fieldData = myField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(); <-- Read-write access to device data

  template <FieldAccessTag FieldAccess = ReadOnly,
            typename Space = stk::ngp::HostSpace>
  typename FieldDataHelper<T, FieldAccess, Space, LayoutSelector<Space>::layout>::FieldDataType
  data() const
  {
    return FieldBase::data<T, FieldAccess, Space, LayoutSelector<Space>::layout>();
  }


  // This is the same as described above, except you can pass in a custom execution space argument
  // that will be used to run any data syncing or updating after a mesh modification that may
  // be necessary.  This is intended for asynchronous execution.

  template <FieldAccessTag FieldAccess = ReadOnly,
            typename Space = stk::ngp::HostSpace,
            typename ExecSpace = stk::ngp::ExecSpace>
  typename AsyncFieldDataHelper<T, FieldAccess, Space, LayoutSelector<Space>::layout, ExecSpace>::FieldDataType
  data(const ExecSpace& execSpace) const
  {
    return FieldBase::data<T, FieldAccess, Space, LayoutSelector<Space>::layout, ExecSpace>(execSpace);
  }

private:

#ifndef DOXYGEN_COMPILE
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );
#endif /* DOXYGEN_COMPILE */

};

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_Field_hpp */

