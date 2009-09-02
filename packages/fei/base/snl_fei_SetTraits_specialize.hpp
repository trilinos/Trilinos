#ifndef _snl_fei_SetTraits_specialize_hpp_
#define _snl_fei_SetTraits_specialize_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_SetTraits.hpp>

namespace snl_fei {

/** specialization for fei::ctg_set<int> */
template<>
struct SetTraits<fei::ctg_set<int> > {
  /** insert item in specified set */
  static void insert(fei::ctg_set<int>* set_obj, int item)
  { set_obj->insert2(item); }
};

}//namespace snl_fei
#endif

