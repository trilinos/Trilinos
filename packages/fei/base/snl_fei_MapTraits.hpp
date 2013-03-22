#ifndef _snl_fei_MapTraits_hpp_
#define _snl_fei_MapTraits_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

namespace snl_fei {

/** Define map traits. For now, only a trait for map.lower_bound and
    map.insert. Generally we simply call the corresponding methods on
    std::map, but these traits will be specialized for map classes that
    may not always have directly equivalent methods (such as hash_map, etc.).
    (hash_map doesn't have lower_bound, and hash_map's insert doesn't allow
    the use of an iterator argument to hint where the object should be
    inserted.)
*/
template<typename MAP_TYPE>
struct MapTraits {
  /** lower_bound */
  static typename MAP_TYPE::iterator
    lower_bound(MAP_TYPE& map_obj,
                typename MAP_TYPE::key_type item)
  { return( map_obj.lower_bound(item) ); }

  /** insert a value using iterator for position hint */
  static void insert(MAP_TYPE& map_obj,
                     typename MAP_TYPE::iterator& pos,
                     typename MAP_TYPE::value_type& val)
  { map_obj.insert(pos, val); }
};

}//namespace snl_fei
#endif

