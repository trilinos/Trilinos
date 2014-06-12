#include <stk_expreval/Constants.hpp>

namespace stk_classic {
namespace expreval {

ConstantMap &
getConstantMap()
{
  static ConstantMap s_constantMap;

  if (s_constantMap.empty()) {
    s_constantMap["E"]     = s_e;
    s_constantMap["PI"]    = s_pi;
    s_constantMap["FALSE"] = s_false;
    s_constantMap["TRUE"]  = s_true;
  }

  return s_constantMap;
}

} // namespace expreval
} // namespace stk_classic
