#ifndef stk_util_environment_CPUTime_hpp
#define stk_util_environment_CPUTime_hpp

namespace stk {

/**
 * @brief Member function <b>cpu_time</b> returns the accumulated cpu time for the process as a
 * double precision value in seconds to "millisecond" accuracy.
 *
 * @return              a <b>double</b> value of the accumulated process CPU runtime.
 */
double cpu_time();

} // namespace stk

#endif // stk_util_environment_CPUTime_hpp
