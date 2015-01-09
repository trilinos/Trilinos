/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_util_SameType_hpp
#define stk_util_util_SameType_hpp

namespace stk_classic {

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

} //namespace stk_classic

#endif

