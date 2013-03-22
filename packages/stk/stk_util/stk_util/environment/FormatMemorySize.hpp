/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_environment_FormatMemorySize_hpp
#define stk_util_environment_FormatMemorySize_hpp

#include <string>

namespace stk {

typedef size_t MemorySize;

std::string formatMemorySize(double time);
std::string formatMemorySize(MemorySize time);

} // namespace stk

#endif // stk_util_environment_FormatMemorySize_hpp
