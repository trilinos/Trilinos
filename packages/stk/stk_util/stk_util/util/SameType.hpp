#ifndef stk_util_util_SameType_hpp
#define stk_util_util_SameType_hpp

namespace stk {

//----------------------------------------------------------------------
/** \class  SameType
 *  \brief  Member <b> enum { value = ... }; </b>
 *          is true if <b> T1 </b> and <b> T2 </b> are the same type.
 *  \ingroup typelist_module
 */
template<typename T1, typename T2>
struct SameType
{ enum { value = false }; };
 
template<typename T>
struct SameType<T,T>
{ enum { value = true }; };

//----------------------------------------------------------------------

} //namespace stk

#endif

