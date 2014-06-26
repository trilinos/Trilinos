#ifndef STK_SIERRA_MESH_DETAILS_RELATION_DESCRIPTOR_HPP
#define STK_SIERRA_MESH_DETAILS_RELATION_DESCRIPTOR_HPP

#include <sierra/util/descriptor.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct relation_descriptor_type {
  typedef size_t       base_type;
  typedef sierra::util::parameters_by_value  parameter_catagory;

  static size_t invalid() { return static_cast<size_t>(-1); }
};

typedef sierra::util::descriptor<details::relation_descriptor_type> relation_descriptor;

} // details
} // mesh
} // sierra

#endif // STK_SIERRA_MESH_DETAILS_RELATION_DESCRIPTOR_HPP
