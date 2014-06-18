#ifndef SIERRA_SIERRA_MESH_BUCKET_KEY_HPP
#define SIERRA_SIERRA_MESH_BUCKET_KEY_HPP

#include <sierra/util/descriptor.hpp>
#include <sierra/util/descriptor_manager.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct bucket_key_type {
  typedef size_t                          base_type;
  typedef sierra::util::parameters_by_value  parameter_catagory;

  static size_t invalid() { return static_cast<size_t>(-1); }
  //needed for descriptor_manager_traits
  static base_type first_valid()   { return static_cast<base_type>(0); }
  static base_type last_valid()    { return static_cast<base_type>(-2); }
};

typedef sierra::util::descriptor<details::bucket_key_type> bucket_key;
typedef sierra::util::descriptor_manager<details::bucket_key_type> bucket_key_manager;

} // namespace details
} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_DETAILS_BUCKET_KEY_HPP
