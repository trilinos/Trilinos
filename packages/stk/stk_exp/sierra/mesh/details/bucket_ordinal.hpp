#ifndef SIERRA_SIERRA_MESH_BUCKET_ORDINAL_HPP
#define SIERRA_SIERRA_MESH_BUCKET_ORDINAL_HPP

#include <sierra/util/descriptor.hpp>

namespace sierra {
namespace mesh {
namespace details {

struct bucket_ordinal_type {
  typedef size_t                          base_type;
  typedef sierra::util::parameters_by_value  parameter_catagory;

  static size_t invalid() { return static_cast<size_t>(-1); }
};

typedef sierra::util::descriptor<details::bucket_ordinal_type> bucket_ordinal;

} // namespace details
} // mesh
} // sierra


#endif // SIERRA_SIERRA_MESH_DETAILS_BUCKET_ORDINAL_HPP

