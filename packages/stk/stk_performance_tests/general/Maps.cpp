// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
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
// 

#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>

#include <gtest/gtest.h>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>

#include <stk_util/environment/WallTime.hpp>

// Only run tests for GCC 4 for now
// Altix preprocessor variables are bizzare: defines __GNUC__, does not define __ICC even though it's intel
#if __GNUC__ == 4 && !defined __itanium__ && !defined __ICC && !defined __PATHSCALE__ && !defined __PGI

#define FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_STD 1
#if defined(__PGI) || defined(__PATHSCALE__)
  #define FLAG_GOOGLE_SPARSE_HASH_MAP 0
  #define FLAG_GOOGLE_DENSE_HASH_MAP 0
#else
  #define FLAG_GOOGLE_SPARSE_HASH_MAP 1
  #define FLAG_GOOGLE_DENSE_HASH_MAP 1
#endif
#define FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1 1
#define FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_BOOST 1
#define FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TEUCHOS_HASHTABLE 0

#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
  #include <tr1/unordered_map>
#endif

#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_BOOST
  #include <boost/unordered_map.hpp>
  #include <algorithm>
#endif

#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TEUCHOS_HASHTABLE
  #include <ostream>
  #include <Teuchos_Hashtable.hpp>
  using Teuchos::Hashtable;
#endif

#if FLAG_GOOGLE_SPARSE_HASH_MAP
  #include <stk_util/util/config_google.h>
  #include HASH_MAP_H
  #include <functional>
  #include <stk_util/util/sparse_hash_map>
  using google::sparse_hash_map;
#endif

#if FLAG_GOOGLE_DENSE_HASH_MAP
  #include <stk_util/util/config_google.h>
  #include HASH_MAP_H
  #include <functional>
  #include <stk_util/util/densehashtable.h>
  #include <stk_util/util/dense_hash_map>
  using google::dense_hash_map;
#endif

using stk::mesh::EntityKey;
using stk::mesh::Entity;

namespace {

const int NUM_KEYS = 10000000;

// Generic performance test
template<class MapType>
void time_map(MapType& testmap, int iters, const std::string& mapname)
{
  std::cout << "For " << mapname << std::endl;

  // Populate map
  double start_time = stk::wall_time();
  for (int i = 0; i < iters ; ++i) {
    testmap[EntityKey(stk::topology::EDGE_RANK,i)] = i;
  }
  double timing0 = stk::wall_dtime(start_time);
  std::cout << mapname << "\ttiming map growth: \t"  << timing0 << " s" << std::endl;

  // Time item lookup - sequential
  start_time = stk::wall_time();
  int counter_to_fool_optimizer = 0;
  for (int i = 0; i < iters; ++i) {
    counter_to_fool_optimizer += testmap[EntityKey(stk::topology::EDGE_RANK,i)];
  }
  double timing = stk::wall_dtime(start_time);
  std::cout << mapname << "\tsequential timing map lookups: \t"  << timing << " s" << std::endl;

  //
  // Time item lookup - random
  //

  // Create list of random ids; pulls the cost of rand and modulo out of the timed section
  int* num_array = new int[iters];
  for (int i=0; i < iters; i++ ) {
    num_array[i] = std::rand() % NUM_KEYS;
  }

  // Do random key lookups
  start_time = stk::wall_time();
  counter_to_fool_optimizer = 0;
  for (int i = 0; i < iters; ++i) {
    counter_to_fool_optimizer += testmap[EntityKey(stk::topology::EDGE_RANK,num_array[i])];
  }
  double timing2 = stk::wall_dtime(start_time);
  delete[] num_array;

  std::cout << mapname << "\trandom timing map lookups: \t"  << timing2 << " s" << std::endl;

  // Code to fool optimizer
  if (counter_to_fool_optimizer == 666) std::cout << "";
}

//
// Hash functions
//

struct stk_entity_rep_hash : public std::unary_function< EntityKey, std::size_t >
{
  inline std::size_t
  operator()(const EntityKey& x) const
  {
    return (std::size_t)(x);
  }
};

} // empty namespace

#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos {

// Teuchos expects the developer to define Teuchos::hashCode
template<>
int hashCode(const EntityKey& x)
{
  return static_cast<int>(x);
}

}
#endif

namespace {

#if FLAG_GOOGLE_SPARSE_HASH_MAP || FLAG_GOOGLE_DENSE_HASH_MAP
struct google_stk_entity_rep_hash : public google::is_integral< EntityKey >
{
  inline std::size_t
  operator()(const EntityKey& x) const
  {
    return (std::size_t)(x);
  }
};
#endif

//
// Performance tests
//

//First map type: std
#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_STD
TEST( map_tests, stdmap)
{
  typedef std::map<EntityKey, int> EntityMap;
  EntityMap testmap;

  time_map(testmap, NUM_KEYS, "std map");
}
#endif //end of std map

//Second map type: tr1
#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
TEST( map_tests, tr1map)
{
  typedef std::tr1::unordered_map<EntityKey, int, stk_entity_rep_hash,
                                  std::equal_to<EntityKey>, std::allocator<std::pair<EntityKey const, int> > > EntityMap;

  EntityMap testmap;

  time_map(testmap, NUM_KEYS, "tr1 map");
}
#endif //end of tr1 map

//Third map type: boost
#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_BOOST
TEST( map_tests, boostmap)
{
  typedef boost::unordered_map< EntityKey, int, stk_entity_rep_hash > EntityMap;

  EntityMap testmap;

  time_map(testmap, NUM_KEYS, "boost map");
}
#endif //end of boost map

//Fourth map type: Teuchos
//NB: include files taken from stk_percept and adapted
#if FLAG_STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TEUCHOS_HASHTABLE
TEST( map_tests, teuchosmap)
{
  typedef Teuchos::Hashtable<EntityKey, int > EntityMap;
  EntityMap testmap;
  EntityMap testmap2;

  time_map(testmap, NUM_KEYS, "teuchos map");
  time_map2(testmap2, NUM_KEYS, "teuchos map");
}
#endif // end of teuchos map

//Sixth map type: google sparse
//From http://code.google.com/p/google-sparsehash
#if FLAG_GOOGLE_SPARSE_HASH_MAP
TEST( map_tests, googlesparsemap)
{
  typedef google::sparse_hash_map< EntityKey, int, stk_entity_rep_hash > EntityMap;

  EntityMap testmap;

  time_map(testmap, NUM_KEYS, "google sparse map");
}
#endif //end of sparse_hash_map

//Seventh map type: google dense
//From http://code.google.com/p/google-sparsehash
#if FLAG_GOOGLE_DENSE_HASH_MAP

template<class MapType> inline void SET_EMPTY_KEY(MapType&, EntityKey /*key*/) {}

template<class K, class V, class H, class E, class A>
inline void SET_EMPTY_KEY(dense_hash_map<K,V,H,E,A>& m, EntityKey key)
{
  m.set_empty_key(key);
}

TEST( map_tests, googledensemap)
{
  typedef google::dense_hash_map< EntityKey, int, stk_entity_rep_hash > EntityMap;

  EntityMap testmap;

  SET_EMPTY_KEY(testmap, EntityKey(stk::topology::NODE_RANK,0));

  time_map(testmap, NUM_KEYS, "google dense map");
}
#endif //end of dense_hash_map

}//namespace <anonymous>

#endif // GCC 4

