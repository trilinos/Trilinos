#ifndef SIERRA_SIERRA_MESH_ENTITY_DESCRIPTOR_HPP
#define SIERRA_SIERRA_MESH_ENTITY_DESCRIPTOR_HPP

#include <sierra/util/descriptor.hpp>
#include <sierra/util/descriptor_manager.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct entity_descriptor_type {
  typedef size_t                          base_type;
  typedef sierra::util::parameters_by_value  parameter_catagory;

  static base_type invalid() { return static_cast<base_type>(-1); }
  //needed for descriptor_manager_traits
  static base_type first_valid()   { return static_cast<base_type>(0); }
  static base_type last_valid()    { return static_cast<base_type>(-2); }
};

typedef sierra::util::descriptor<details::entity_descriptor_type> entity_descriptor;
typedef sierra::util::descriptor_manager<details::entity_descriptor_type> entity_descriptor_manager;

} // details
} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_DETAILS_ENTITY_DESCRIPTOR_HPP
