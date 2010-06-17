/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_Map_h
#define IOSS_Ioss_Map_h

#include <Ioss_CodeTypes.h>
#include <vector>
#include <string>

namespace Ioss {
  typedef std::vector<int> MapContainer;
  // map global to local ids
  typedef std::pair<int,int> IdPair;
  typedef std::vector<IdPair> ReverseMapContainer;

  typedef ReverseMapContainer::value_type RMCValuePair;

  // Class to support storing global/local element id map in sorted vector...
  class IdPairCompare
    {
    public:
      IdPairCompare() {}
      bool operator()(const IdPair& lhs, const IdPair &rhs) const
	{ return key_less(lhs.first, rhs.first); }
      bool operator()(const IdPair& lhs, const IdPair::first_type &k) const
	{ return key_less(lhs.first, k); }
      bool operator()(const IdPair::first_type& k, const IdPair &rhs) const
	{ return key_less(k, rhs.first); }
      // Assignment operator
      // Copy constructor
    private:
      bool key_less(const IdPair::first_type &k1, const IdPair::first_type &k2) const
	{ return k1 < k2; }
    };

  class IdPairEqual
    {
    public:
      IdPairEqual() {}
      bool operator()(const IdPair& lhs, const IdPair &rhs) const
	{ return key_equal(lhs.first, rhs.first); }
      bool operator()(const IdPair& lhs, const IdPair::first_type &k) const
	{ return key_equal(lhs.first, k); }
      bool operator()(const IdPair::first_type& k, const IdPair &rhs) const
	{ return key_equal(k, rhs.first); }
      // Assignment operator
      // Copy constructor
    private:
      bool key_equal(const IdPair::first_type &k1, const IdPair::first_type &k2) const
	{ return k1 == k2; }
    };

  class Map {
  public:
    Map(); // Default constructor

    int global_to_local(int global, bool must_exist = true) const;
    int local_to_global(int local) const;

    static void build_reverse_map(ReverseMapContainer *Map, const int *ids,
				  int num_to_get, int offset,
				  const std::string& type, int processor);

    // Determines whether the input map is sequential (map[i] == i)
    // Assumes that map is '1-based', size stored in [0]
    static bool is_sequential(const std::vector<int>& the_map);

    static void verify_no_duplicate_ids(ReverseMapContainer &reverse_map,
					const std::string &type, int processor);

    int  entityCount;
    bool sequentialG2L; // true if reverse node map is sequential (local==global)
    bool entityReordered;
    MapContainer forwardMap;
    ReverseMapContainer reverseMap;
    MapContainer reorderMap;
  private:
    Map(const Map& from); // do not implement
    Map& operator=(const Map& from); // do not implement
  };
}

#endif // IOSS_Ioss_Map_h
