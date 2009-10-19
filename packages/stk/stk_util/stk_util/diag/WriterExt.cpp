#include <stk_util/diag/WriterExt.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/environment/Demangle.hpp>

namespace stk {
namespace diag {

Writer &
operator<<(
  Writer &                      dout,
  const std::type_info &        t)
{
  if (dout.shouldPrint())
    dout << stk::demangle(t.name());
  return dout;
}

} // namespace diag
} // namespace stk
