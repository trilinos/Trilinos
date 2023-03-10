#ifndef GEO_CLASSIFICATION_H
#define GEO_CLASSIFICATION_H

#include <ostream>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct GeoClassification
{
    explicit GeoClassification(const int_least8_t dim_ = 2, const int id_ = 0)
      : dim(dim_)
      , id(id_)
    {}

    int_least8_t dim;
    int id;
};

inline std::ostream& operator<<(std::ostream& os, const GeoClassification& g)
{
  os << "dim = " << static_cast<int>(g.dim) << ", id = " << g.id;
  return os;
}

inline bool operator==(const GeoClassification& lhs, const GeoClassification& rhs)
{
  return lhs.dim == rhs.dim && lhs.id == rhs.id;
}

inline bool operator!=(const GeoClassification& lhs, const GeoClassification& rhs)
{
  return !(lhs == rhs);
}

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
