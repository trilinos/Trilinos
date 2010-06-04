/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_Field_h
#define SIERRA_Ioss_Field_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <vector>

namespace Ioss {

  class VariableType;
  class GroupingEntity;
  class Transform;

  class Field {
  public:
    enum BasicType {INVALID = -1, REAL, INTEGER, COMPLEX, STRING};
    enum RoleType {INTERNAL, MESH, ATTRIBUTE, COMMUNICATION, INFORMATION, REDUCTION, TRANSIENT};

    Field();

    // Create a field named 'name' that contains values of type 'type'
    // in a storage format of type 'storage'.  There are 'value_count'
    // items in the field.
    Field(const std::string &name, const BasicType type,
	  const std::string &storage,
	  const RoleType role, size_t value_count);
    
    Field(const std::string &name, const BasicType type,
	  const std::string &storage, int copies,
	  const RoleType role, size_t value_count);

    Field(const std::string &name, const BasicType type,
	  const VariableType *storage,
	  const RoleType role, size_t value_count);
    
    // Create a field from another field.
    Field(const Field&);
    Field& operator=(const Field&);

    // Compare two fields (used for STL container)
    bool operator<(const Field& other) const;

    ~Field();

    bool is_valid()    const {return type_ != INVALID;}
    bool is_invalid()  const {return type_ == INVALID;}

    const std::string &get_name()  const {return name_;}
    BasicType     get_type() const {return type_;}

    const VariableType *raw_storage() const {return rawStorage_;}
    const VariableType *transformed_storage() const {return transStorage_;}

    size_t raw_count() const {return rawCount_;} // Number of items in field
    size_t transformed_count() const {return transCount_;} // Number of items in field

    size_t get_size() const; // data size (in bytes) required to hold entire field

    RoleType      get_role() const {return role_;}

    void reset_count(size_t new_count); // new number of items in field
    void reset_type(BasicType new_type); // new type of items in field.

    // Verify that data_size is valid.  Return value is the maximum
    // number of entities to get ('RawCount')
    // Throws runtime error if data_size too small.
    size_t verify(size_t data_size) const;

    // Verify that the type 'the_type' matches the field's type.
    // throws exception if the types don't match.
    void check_type(BasicType the_type) const;
    
    bool add_transform(Transform* transform);
    bool transform(void *data);

  private:
    std::string  name_;

    size_t    rawCount_; // Count of items in field before transformation
    size_t    transCount_; // Count of items in field after transformed
    BasicType     type_;
    RoleType      role_;

    const VariableType* rawStorage_;  // Storage type of raw field
    const VariableType* transStorage_; // Storage type after transformation

    size_t size_; // maximum data size (in bytes) required to hold entire field
    std::vector<Transform*>  transforms_;
  };
}
#endif
