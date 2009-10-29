#ifndef stk_linsys_ImplDetails_hpp
#define stk_linsys_ImplDetails_hpp

#include <stk_linsys/FieldIdMap.hpp>

#include <limits>
#include <map>
#include <stk_mesh/base/Field.hpp>

namespace stk {

/** Linear-System Assembly
*/
namespace linsys {

/** Implementation Details -- not generally of interest for the public API
*/
namespace impl {

/** Given a map and a Field, return the int id that the field is mapped to.
  The field will be added to the map if not already present.
*/
int map_field_to_int(FieldIdMap& field_id_map,
                     const stk::mesh::FieldBase& field);

/** Given a map and a Field, return the int id that the field is mapped to.
  If the field is not found in the map, an exception is thrown.
*/
int
query_field_to_int_mapping(const FieldIdMap& field_id_map,
                           const stk::mesh::FieldBase& field);

/** Given an integer field_id, return a reference to the corresponding field.
 Throw an exception if field_id not found.
*/
const stk::mesh::FieldBase* get_field(const FieldIdMap& field_id_map,
                                      int field_id);

/** Given an EntityId, return the value as an int.
  Throws an exception if id is too large to represent as an int.
*/
int entityid_to_int(stk::mesh::EntityId id);

/** Given an EntityType, return the value as an int.
  Throws an exception if id is too large to represent as an int.
*/
int entitytype_to_int(stk::mesh::EntityType entity_type);

/** Determine whether 'id' can be converted to an int.
 * If so, do nothing. If 'id' is too large to be represented
 * as an int, throw an exception (std::runtime_error).
 */
template<typename T>
void verify_convertible_to_int(T id, const char* caller)
{
  if (sizeof(T) <= sizeof(int)) return;

  T intmax = std::numeric_limits<int>::max();
  if (intmax < id) {
    std::ostringstream msg;
    msg << caller << " ERROR, id " << id << " is too large to convert to int.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }
}

}//namespace impl
}//namespace linsys
}//namespace stk

#endif

