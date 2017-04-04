// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#ifndef IOSS_Ioss_Field_h
#define IOSS_Ioss_Field_h

#include <Ioss_CodeTypes.h>
#include <stddef.h> // for size_t
#include <string>   // for string
#include <vector>   // for vector

namespace Ioss {

  class GroupingEntity;
  class Transform;
  class VariableType;

  /** \brief Holds metadata for bulk data associated with a GroupingEntity.
   */
  class Field
  {
  public:
    /** \brief The basic data type held in the field.
     */
    enum BasicType {
      INVALID = -1,
      REAL    = 1,
      DOUBLE  = 1,
      INTEGER = 4,
      INT32   = 4,
      INT64   = 8,
      COMPLEX,
      STRING,
      CHARACTER
    };

    static Ioss::Field::BasicType get_field_type(char /*dummy*/) { return CHARACTER; }
    static Ioss::Field::BasicType get_field_type(double /*dummy*/) { return DOUBLE; }
    static Ioss::Field::BasicType get_field_type(int /*dummy*/) { return INTEGER; }
    static Ioss::Field::BasicType get_field_type(int64_t /*dummy*/) { return INT64; }
    static Ioss::Field::BasicType get_field_type(Complex /*dummy*/) { return COMPLEX; }
    static Ioss::Field::BasicType get_field_type(std::string /*dummy*/) { return STRING; }

    /* \brief Categorizes the type of information held in the field.
     */
    enum RoleType {
      INTERNAL,
      MESH,      /**< A field which is used to define the basic geometry
                      or topology of the model and is not normally transient
                      in nature. Examples would be element connectivity or
                      nodal coordinates. */
      ATTRIBUTE, /**< A field which is used to define an attribute on an
                      EntityBlock derived class. Examples would be thickness
                      of the elements in a shell element block or the radius
                      of particles in a particle element block. */
      COMMUNICATION,
      INFORMATION,
      REDUCTION, /**< A field which typically summarizes some transient data
                      about an entity. The size of this field is typically not
                      proportional to the number of entities in a GroupingEntity.
                      An example would be average displacement over a group of
                      nodes or the kinetic energy of a model. This data is also
                      transient. */
      TRANSIENT  /**< A field which is typically calculated at multiple steps
                      or times in an analysis. These are typically "results"
                      data. Examples would be nodal displacement or element
                      stress. */
    };

    Field();

    // Create a field named 'name' that contains values of type 'type'
    // in a storage format of type 'storage'.  There are 'value_count'
    // items in the field.
    Field(std::string name, const BasicType type, const std::string &storage, const RoleType role,
          size_t value_count, size_t index = 0);

    Field(std::string name, const BasicType type, const std::string &storage, int copies,
          const RoleType role, size_t value_count, size_t index = 0);

    Field(std::string name, const BasicType type, const VariableType *storage, const RoleType role,
          size_t value_count, size_t index = 0);

    // Create a field from another field.
    Field(const Field & /*from*/);
    Field &operator=(const Field & /*from*/);

    // Compare two fields (used for STL container)
    bool operator<(const Field &other) const;

    ~Field();

    bool is_valid() const { return type_ != INVALID; }
    bool is_invalid() const { return type_ == INVALID; }

    const std::string &get_name() const { return name_; }

    /** \brief Get the basic data type of the data held in the field.
     *
     * \returns the basic data type of the data held in the field.
     */
    BasicType get_type() const { return type_; }

    const VariableType *raw_storage() const { return rawStorage_; }
    const VariableType *transformed_storage() const { return transStorage_; }

    size_t raw_count() const { return rawCount_; }           // Number of items in field
    size_t transformed_count() const { return transCount_; } // Number of items in field

    size_t get_size() const; // data size (in bytes) required to hold entire field

    /** \brief Get the role (MESH, ATTRIBUTE, TRANSIENT, REDUCTION, etc.) of the data in the field.
     *
     * \returns The RoleType of the data in the field.
     */
    RoleType get_role() const { return role_; }

    size_t get_index() const { return index_; }
    void set_index(size_t index) const { index_ = index; }

    void reset_count(size_t new_count);  // new number of items in field
    void reset_type(BasicType new_type); // new type of items in field.

    // Verify that data_size is valid.  Return value is the maximum
    // number of entities to get ('RawCount')
    // Throws runtime error if data_size too small.
    size_t verify(size_t data_size) const;

    // Verify that the type 'the_type' matches the field's type.
    // throws exception if the types don't match.
    void check_type(BasicType the_type) const;

    bool is_type(BasicType the_type) const { return the_type == type_; }

    bool add_transform(Transform *my_transform);
    bool transform(void *data);
    bool has_transform() const { return !transforms_.empty(); }

  private:
    std::string name_;

    size_t         rawCount_;   // Count of items in field before transformation
    size_t         transCount_; // Count of items in field after transformed
    size_t         size_;       // maximum data size (in bytes) required to hold entire field
    mutable size_t index_; // Optional flag that can be used by a client to indicate an ordering.
                           // Unused by field itself.
    BasicType type_;
    RoleType  role_;

    const VariableType *rawStorage_;   // Storage type of raw field
    const VariableType *transStorage_; // Storage type after transformation

    std::vector<Transform *> transforms_;
  };
}
#endif
