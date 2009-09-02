#ifndef _snl_fei_MapTraits_specialize_hpp_
#define _snl_fei_MapTraits_specialize_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_MapTraits.hpp>
#include <snl_fei_MapContig.hpp>

namespace snl_fei {
#if defined(FEI_HASH_MAP) && defined(FEI_HASH_SET)
template<>
struct MapTraits<FEI_HASH_MAP<int,FEI_HASH_SET<int>*> > {
  static FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,FEI_HASH_SET<int>*>& map_obj,
                int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,FEI_HASH_SET<int>*>& map_obj,
                     FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::value_type& val)
  { map_obj.insert(val); }
};
#endif

#if defined(FEI_HASH_MAP)

template<>
struct MapTraits<FEI_HASH_MAP<int,fei::ctg_set<int>*> >
{
  static FEI_HASH_MAP<int,fei::ctg_set<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,fei::ctg_set<int>*>& map_obj, int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,fei::ctg_set<int>*>& map_obj,
                     FEI_HASH_MAP<int,fei::ctg_set<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,fei::ctg_set<int>*>::value_type& val)
 { map_obj.insert(val); }};

template<>
struct MapTraits<FEI_HASH_MAP<int,std::set<int>*> > {
  static FEI_HASH_MAP<int,std::set<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,std::set<int>*>& map_obj,
                int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,std::set<int>*>& map_obj,
                     FEI_HASH_MAP<int,std::set<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,std::set<int>*>::value_type& val)
  { map_obj.insert(val); }
};
#endif

}//namespace snl_fei
#endif

