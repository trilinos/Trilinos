/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_StaticAssert_hpp
#define STK_UTIL_UTIL_StaticAssert_hpp

namespace stk {

//----------------------------------------------------------------------
/** \brief  Compile-time assertion
 *  \ingroup util_module
 *  If the compile-time <b> expression </b> is true then defines
 *  - <b> enum { OK = true }; </b>
 *  - <b> static bool ok() { return true ; } </b>
 *
 * \todo REFACTOR Does StaticAssert belong here? Anywhere?
 */
template<bool expression> struct StaticAssert {};

template<> struct StaticAssert<true> {
  enum { OK = true };
  static bool ok() { return true ; }
};

//----------------------------------------------------------------------

} //namespace stk

#endif // STK_UTIL_UTIL_StaticAssert_hpp

