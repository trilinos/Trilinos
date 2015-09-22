/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_set_dof_mappings_hpp_
#define _fei_set_dof_mappings_hpp_

#include <fei_macros.hpp>

#include <fei_DofMapper.hpp>

namespace fei {

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
void set_dof_mappings(GlobalOrdinal first_index,
                      fei::DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>& dof_mapper)
{
  typedef typename fei::DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::DofMap DofMap;
  typedef typename fei::DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::IdxMap IdxMap;

  typename DofMap::iterator
    d_iter = dof_mapper.begin_dof(), d_end = dof_mapper.end_dof();

  IdxMap& idxmap = dof_mapper.get_idx_dof_map();

  GlobalOrdinal index = first_index;
  for(; d_iter != d_end; ++d_iter) {
    LocalOrdinal fieldsize = dof_mapper.getFieldSize(d_iter->first.field());
    d_iter->second = index;
    idxmap.insert(std::make_pair(index, &(d_iter->first)));
    index += fieldsize;
  }
  dof_mapper.set_maps_are_valid(true);
}

}//namespace fei

#endif

