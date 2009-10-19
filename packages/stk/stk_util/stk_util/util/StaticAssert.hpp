#ifndef stk_util_util_StaticAssert_hpp
#define stk_util_util_StaticAssert_hpp

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

#endif

