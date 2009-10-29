
#include <limits>
#include <sstream>
#include <stdexcept>

#include <stk_linsys/ImplDetails.hpp>

namespace stk {
namespace linsys {
namespace impl {

int
map_field_to_int(FieldIdMap& field_id_map,
                 const stk::mesh::FieldBase& field)
{
  FieldIdMap::iterator iter = field_id_map.find(&field);

  if (iter == field_id_map.end()) {
    iter = field_id_map.insert(iter, std::make_pair(&field,field_id_map.size()));
  }

  return iter->second;
}

int
query_field_to_int_mapping(const FieldIdMap& field_id_map,
                           const stk::mesh::FieldBase& field)
{
  FieldIdMap::const_iterator iter = field_id_map.find(&field);

  if (iter == field_id_map.end()) {
    std::ostringstream msg;
    msg << "stk::linsys::query_field_to_int_mapping ERROR: "
        << " field with name '"<<field.name()<<"' not found in field-to-int map.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  return iter->second;
}

const stk::mesh::FieldBase*
get_field(const FieldIdMap& field_id_map,
          int field_id)
{
  FieldIdMap::const_iterator
    iter = field_id_map.begin(), iter_end = field_id_map.end();

  while(iter!=iter_end && iter->second != field_id) ++iter;

  if (iter == iter_end) {
    std::ostringstream msg;
    msg << "stk::linsys::get_dof ERROR: "
     << "field_id ("<<field_id<<") returned from fei query is not mapped to a stk::mesh::Field.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  return iter->first;
}

int entitytype_to_int(stk::mesh::EntityType entity_type)
{
  int int_type = static_cast<int>(entity_type);

  verify_convertible_to_int(entity_type, "stk::linsys::entitytype_to_int");

  return int_type;
}

int entityid_to_int(stk::mesh::EntityId id)
{
  int int_id = static_cast<int>(id);

  verify_convertible_to_int(id, "stk::linsys::entityid_to_int");

  return int_id;
}

}//namespace impl
}//namespace linsys
}//namespace stk

