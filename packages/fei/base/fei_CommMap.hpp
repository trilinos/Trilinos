/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_CommMap_hpp_
#define _fei_CommMap_hpp_

#include <fei_macros.hpp>
#include <fei_ArrayUtils.hpp>
#include <set>
#include <map>

namespace fei {

/** Map that maps processors to associated vectors of items to be sent or recv'd. */
template<typename T>
struct CommMap {
  typedef std::map<int,std::vector<T> > Type;
};

/** Given a proc and an array of items, add the mapping
       proc -> items
   to the given comm_map.
   Optionally ensure that the comm_map's vector of items for proc remains
   sorted and unique.
*/
template<typename T>
void addItemsToCommMap(int proc, size_t numItems, const T* items,
                       typename CommMap<T>::Type& comm_map,
                       bool keep_sorted_and_unique = true)
{
  typename CommMap<T>::Type::iterator iter = comm_map.find(proc);
  if (iter == comm_map.end()) {
    iter = comm_map.insert(std::make_pair(proc,std::vector<T>())).first;
  }

  std::vector<T>& comm_items = iter->second;

  if (keep_sorted_and_unique) {
    for(size_t i=0; i<numItems; ++i) {
      fei::sortedListInsert(items[i], comm_items);
    }
  }
  else {
    for(size_t i=0; i<numItems; ++i) {
      comm_items.push_back(items[i]);
    }
  }
}

} //namespace fei

#endif // _fei_CommMap_hpp_

