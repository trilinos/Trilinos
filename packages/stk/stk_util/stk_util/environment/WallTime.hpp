#ifndef stk_util_environment_WallTime_hpp
#define stk_util_environment_WallTime_hpp

namespace stk {

/**
 * @brief Member function <b>wall_time</b> returns the epoch as a double precision value in seconds
 * to "millisecond" accuracy.
 *
 * @return              a <b>double</b> value of the seconds since the 1/1/1970 epoch.
 */
double wall_time();

double wall_dtime(double &t);

} // namespace stk

#endif // stk_util_environment_WallTime_hpp
