/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_environment_WallTime_hpp
#define stk_util_environment_WallTime_hpp

namespace stk_classic {

/**
 * @brief Member function <b>wall_time</b> returns the epoch as a double precision value in seconds
 * to "millisecond" accuracy.
 *
 * @return              a <b>double</b> value of the seconds since the 1/1/1970 epoch.
 */
double wall_time();

double wall_dtime(double &t);

} // namespace stk_classic

#endif // stk_util_environment_WallTime_hpp
