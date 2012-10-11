// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef IOSS_Ioss_Map_h
#define IOSS_Ioss_Map_h

#include <Ioss_CodeTypes.h>
#include <vector>
#include <string>
#include <stdint.h>

namespace Ioss {
  typedef std::vector<int64_t> MapContainer;
  // map global to local ids
  typedef std::pair<int64_t,int64_t> IdPair;
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

    int64_t global_to_local(int64_t global, bool must_exist = true) const;
    int64_t local_to_global(int64_t local) const;

    void build_reverse_map(int processor);

    static void build_reverse_map(ReverseMapContainer *Map, const int64_t *ids,
				  int64_t num_to_get, int64_t offset, int processor);

    static void build_reverse_map(ReverseMapContainer *Map, const int *ids,
				  int num_to_get, int offset, int processor);

    void build_reorder_map(int64_t start, int64_t count);

    // Determines whether the input map is sequential (map[i] == i)
    // Assumes that map is '1-based', size stored in [0]
    static bool is_sequential(const MapContainer& the_map);

    static void verify_no_duplicate_ids(ReverseMapContainer &reverse_map,
					int processor);

    int64_t  entityCount;
    bool sequentialG2L; // true if reverse node map is sequential (local==global)
    bool entityReordered;

    MapContainer        map;
    MapContainer        reorder;
    ReverseMapContainer reverse;

  private:
    Map(const Map& from); // do not implement
    Map& operator=(const Map& from); // do not implement
  };
}

#endif // IOSS_Ioss_Map_h
