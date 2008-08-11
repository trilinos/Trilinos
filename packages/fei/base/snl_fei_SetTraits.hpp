#ifndef _snl_fei_SetTraits_hpp_
#define _snl_fei_SetTraits_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

namespace snl_fei {

/** Define traits 'set' objects. For now, only 1 trait: inserting an
   item into a set, in cases where we don't care about the return-value.
   In general, we just call the set's insert method.
   But this trait will be specialized (in other headers) to call the
   'insert2' method for special fei set objects.
*/
template<typename SET_TYPE>
struct SetTraits {
  /** insert item in specified set */
  static void insert(SET_TYPE* set_obj, typename SET_TYPE::key_type item)
  { set_obj->insert(item); }
};

}//namespace snl_fei

#endif

