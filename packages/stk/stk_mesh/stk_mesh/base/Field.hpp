// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/FieldTraits.hpp>

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
 *   - Field type: Defines the type of this field's data. This can be scalar or multidimensional.
 *     If multidimensional (dimensions >= 1), a specification of each dimension will need
 *     to be provided in the form of a shards::ArrayDimTag. Many common ArrayDimTags are defined
 *     in stk_mesh/base/CoordinateSystems.hpp, but apps are free to define their own as well.
 *     - Examples:
 *       - A scalar double field-type:
 *           stk::mesh::Field<double>
 *       - A vector field type using the Cartesian array-dim-tag from CoordinateSystems:
 *           stk::mesh::Field<double, stk::mesh::Cartesian>
 *
 *   - Field declaration: Defines/creates a specific field of a specified type. The API for
 *     declaring fields is in MetaData. Declaration of a field requires a field-type,
 *     field name, and number of states. The number of states defines the "memory" of the Field;
 *     if number-of-states is N, the last N values of the field are stored, and the user
 *     can advance/rotate the states by calling BulkData::update_field_data_states(). FieldState.hpp
 *     contains an enum that clients should use to refer to the various states.
 *     - Examples:
 *       typedef stk::mesh::Field<double, stk::mesh::Cartesian>  CoordFieldType;
 *       CoordFieldType& coord_field = meta.declare_field<CoordFieldType>("<name>", <num_states>);
 *
 *   - Field restrictions: Defines the set of entities that have a field and the dimensions of the
 *     field (maximum of 7 dimensions).  FieldRestrictions are applied to Parts and entities pick
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
 *     we strongly recommend using the put_field free functions to create your field restrictions. A
 *     field may have many restrictions, but the restrictions have to be compatible with each other;
 *     restrictions can overlap if an entity is a member of two parts A,B and you declare restrictions
 *     for a field for both A and B. If such a situation arises, the dimensionality of both restrictions
 *     must be identical.
 *     - Examples:
 *       - Put a scalar field on all nodes in a part
 *           typedef stk::mesh::Field<double>  ScalarFieldType;
 *           ScalarFieldType& field = meta.declare_field<ScalarFieldType>("<name>", <num_states>);
 *           stk::mesh::put_field(field, stk::topology::NODE_RANK, <part>);
 *
 *       - Put a 1d (of size 3) field of doubles on all nodes in a part
 *           typedef stk::mesh::Field<double, stk::mesh::Cartesian3d>  CoordFieldType;
 *           CoordFieldType& field = meta.declare_field<CoordFieldType>("<name>", <num_states>);
 *           stk::mesh::put_field(field, stk::topology::NODE_RANK, <part>, Cartesian3d::Size);
 *
 *       - Put a 2d (of sizes 3 and 3) field of doubles on all nodes in a part
 *           typedef stk::mesh::Field<double, stk::mesh::Cartesian3d, stk::mesh::Cartesian3d> MultiDimFieldType;
 *           MultiFieldType& field = meta.declare_field<MultiFieldType>("<name>", <num_states>);
 *           stk::mesh::put_field(field, stk::topology::NODE_RANK, <part>, Cartesian3d::Size, Cartesian3d::Size);
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
 *   - stk_mesh/base/FieldTraits - Defines API for querying field-types
 *   - stk_mesh/base/CoordinateSystems - Defines common ArrayDimTags (used for defining Field types)
 *   - stk_mesh/base/TopologyDimensions - Contains useful Field types / ArrayDimTags for setting up field relations
 *   - Type-related in ArrayDim, DataTraits
 *
 * - TODO Describe relationship with Buckets
 */
// Implementation Details:
// The template arguments below describe the field type. Scalar is the scalar
// type of data contained by the field. The TagN describe each dimension of the
// Field, these are expected to be ArrayDimTags. Unused dimensions can be ignored.
template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class Field : public FieldBase {
public:
  Field(
       MetaData                   * arg_mesh_meta_data ,
       unsigned                     arg_ordinal ,
       const std::string          & arg_name ,
       const DataTraits           & arg_traits ,
       unsigned                     arg_rank,
       const shards::ArrayDimTag  * const * arg_dim_tags,
       unsigned                     arg_number_of_states ,
       FieldState                   arg_this_state
       )
    : FieldBase(arg_mesh_meta_data,
        arg_ordinal,
        arg_name,
        arg_traits,
        arg_rank,
        arg_dim_tags,
        arg_number_of_states,
        arg_this_state
        )
  {}

  Field(
       MetaData                   * arg_mesh_meta_data ,
       stk::topology::rank_t        arg_entity_rank ,
       unsigned                     arg_ordinal ,
       const std::string          & arg_name ,
       const DataTraits           & arg_traits ,
       unsigned                     arg_rank,
       const shards::ArrayDimTag  * const * arg_dim_tags,
       unsigned                     arg_number_of_states ,
       FieldState                   arg_this_state
       )
    : FieldBase(arg_mesh_meta_data,
        arg_entity_rank,
        arg_ordinal,
        arg_name,
        arg_traits,
        arg_rank,
        arg_dim_tags,
        arg_number_of_states,
        arg_this_state
        )
  {}

  /** \brief  Query this field for a given field state. */
  Field & field_of_state( FieldState input_state ) const {
#ifndef NDEBUG
    return dynamic_cast<Field &>( * FieldBase::field_state(input_state) );
#else // NDEBUG
    return static_cast<Field &>( * FieldBase::field_state(input_state) );
#endif
  }

  virtual ~Field(){}

  virtual std::ostream& print_data(std::ostream& out, void* data, unsigned size_per_entity) const
  {
    const unsigned num_scalar_values = size_per_entity / sizeof(Scalar);
    Scalar* casted_data = reinterpret_cast<Scalar*>(data);

    out << "{";
    for (unsigned i = 0; i < num_scalar_values; ++i) {
      out << casted_data[i] << " ";
    }
    out << "}";

    return out;
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

