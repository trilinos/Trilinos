#include <stk_expreval/Constants.hpp>

namespace stk {
namespace expreval {

// const double s_false	= 0.0;
// const double s_true	= 1.0;
// const double s_e	= 2.7182818284590452354;
// const double s_pi	= 3.14159265358979323846;


ConstantMap &
getConstantMap()
{
  static ConstantMap s_constantMap;

  if (s_constantMap.empty()) {
    s_constantMap["E"] = s_e;
    s_constantMap["PI"] = s_pi;
    s_constantMap["FALSE"] = s_false;
    s_constantMap["TRUE"] = s_true;
  }

  return s_constantMap;
}

} // namespace expreval
} // namespace stk
