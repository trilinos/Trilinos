#ifndef ABSTRACT_CDT_INTERFACE
#define ABSTRACT_CDT_INTERFACE

#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class AbstractCDTInterface
{
  public:
    virtual ~AbstractCDTInterface() = default;

  private:
    virtual void triangulate(const utils::impl::Projection& proj) = 0;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif