/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <limits>
#include <sstream>
#include <stdexcept>

#include <stk_linsys/ImplDetails.hpp>

namespace stk_classic {
namespace linsys {
namespace impl {

int
map_field_to_int(FieldIdMap& field_id_map,
                 const stk_classic::mesh::FieldBase& field)
{
  FieldIdMap::iterator iter = field_id_map.find(&field);

  if (iter == field_id_map.end()) {
    iter = field_id_map.insert(iter, std::make_pair(&field,field_id_map.size()));
  }

  return iter->second;
}

int
query_field_to_int_mapping(const FieldIdMap& field_id_map,
                           const stk_classic::mesh::FieldBase& field)
{
  FieldIdMap::const_iterator iter = field_id_map.find(&field);

  if (iter == field_id_map.end()) {
    std::ostringstream msg;
    msg << "stk_classic::linsys::query_field_to_int_mapping ERROR: "
        << " field with name '"<<field.name()<<"' not found in field-to-int map.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  return iter->second;
}

const stk_classic::mesh::FieldBase*
get_field(const FieldIdMap& field_id_map,
          int field_id)
{
  FieldIdMap::const_iterator
    iter = field_id_map.begin(), iter_end = field_id_map.end();

  while(iter!=iter_end && iter->second != field_id) ++iter;

  if (iter == iter_end) {
    std::ostringstream msg;
    msg << "stk_classic::linsys::get_dof ERROR: "
     << "field_id ("<<field_id<<") returned from fei query is not mapped to a stk_classic::mesh::Field.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  return iter->first;
}

int entitytype_to_int(stk_classic::mesh::EntityRank entity_rank)
{
  int int_type = static_cast<int>(entity_rank);

  verify_convertible_to_int(entity_rank, "stk_classic::linsys::entitytype_to_int");

  return int_type;
}

int entityid_to_int(stk_classic::mesh::EntityId id)
{
  int int_id = static_cast<int>(id);

  verify_convertible_to_int(id, "stk_classic::linsys::entityid_to_int");

  return int_id;
}

}//namespace impl
}//namespace linsys
}//namespace stk_classic

