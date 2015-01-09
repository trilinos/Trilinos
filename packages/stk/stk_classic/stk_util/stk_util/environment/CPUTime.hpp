/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_environment_CPUTime_hpp
#define stk_util_environment_CPUTime_hpp

namespace stk_classic {

/**
 * @brief Member function <b>cpu_time</b> returns the accumulated cpu time for the process as a
 * double precision value in seconds to "millisecond" accuracy.
 *
 * @return              a <b>double</b> value of the accumulated process CPU runtime.
 */
double cpu_time();

} // namespace stk_classic

#endif // stk_util_environment_CPUTime_hpp
