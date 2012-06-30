#ifndef stk_expreval_Constants_hpp
#define stk_expreval_Constants_hpp

#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>

#include <stk_util/util/string_case_compare.hpp>

namespace stk {
namespace expreval {

/**
 * @brief Typedef <b>ConstantMap</b> maps a constant name to a double constant.
 * The mapping is case insensitive.
 */
typedef std::map<std::string, double, LessCase> ConstantMap;

const double s_false	= 0.0;
const double s_true	= 1.0;
const double s_e	= 2.7182818284590452354;
const double s_pi	= 3.14159265358979323846;

/**
 * @brief Member function <b>getConstantMap</b> returns s reference to the defined
 * constants.
 *
 * @return			a <b>ConstantMap</b> reference to the defined
 *				constants.
 */
ConstantMap &getConstantMap();

} // namespace expreval
} // namespace stk

#endif // stk_expreval_Constants_hpp
