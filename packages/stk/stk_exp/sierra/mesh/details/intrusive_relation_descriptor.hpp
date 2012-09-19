#ifndef SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_DESCRIPTOR_HPP
#define SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_DESCRIPTOR_HPP

#include <sierra/util/descriptor.hpp>

namespace sierra {
namespace mesh {
namespace details {

class intrusive_relation;

struct intrusive_relation_descriptor_type {
  typedef const intrusive_relation *       base_type;
  typedef sierra::util::parameters_by_value  parameter_catagory;

  static size_t invalid() { return (size_t)NULL; }
};

typedef sierra::util::descriptor<details::intrusive_relation_descriptor_type> intrusive_relation_descriptor;

} // details
} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_RELATION_DESCRIPTOR_HPP
