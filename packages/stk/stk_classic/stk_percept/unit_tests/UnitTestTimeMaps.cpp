/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#define DO_TEST_TIME_MAPS 0

#define NOMALLOC_ARRAY_CHECK_SIZES
#include <stk_percept/NoMallocArray.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <boost/unordered_map.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>
#include <map>

//#include <stk_percept/sparsehash/include/google/sparse_hash_map>
#define USE_SPARSEHASH 0
#if USE_SPARSEHASH
#include <google/sparse_hash_map>
#include <google/dense_hash_map>
#endif

namespace stk {
namespace percept {
namespace unit_tests {

#define EXTRA_PRINT 0


//=============================================================================
//=============================================================================
//=============================================================================
STKUNIT_UNIT_TEST(unit_tests_percept, noMallocArray)
{
  NoMallocArray<unsigned, 4u> nma;
  STKUNIT_EXPECT_TRUE(nma.max_size() == 4u);
  STKUNIT_EXPECT_TRUE(nma.size() == 0u);
  STKUNIT_EXPECT_TRUE(nma.max_size() == nma.max_capacity());

  bool got_it = false;
  try {
    NoMallocArray<unsigned, 4> nma1(5, 0u);
  }
  catch (std::runtime_error& e)
  {
    got_it = true;
  }
  if (!got_it)
  {
    STKUNIT_EXPECT_TRUE(false);
  }

  NoMallocArray<unsigned, 4> nma2(4u, 0u);
  std::cout << "nma.size()= " << nma2.size() << std::endl;
  STKUNIT_EXPECT_TRUE(nma2.size() == 4u);

  nma2.resize(0u);
  STKUNIT_EXPECT_TRUE(nma2.size() == 0u);
  nma2.insert(0u);
  nma2.insert(1u);
  nma2.insert(2u);
  nma2.insert(3u);

  STKUNIT_EXPECT_TRUE(nma2.size() == 4u);
  STKUNIT_EXPECT_TRUE(nma2[1] == 1);
  STKUNIT_EXPECT_TRUE(*nma2.begin() == 0);
  STKUNIT_EXPECT_TRUE(*(nma2.end()-1) == 3);
}

//=============================================================================
//=============================================================================
//=============================================================================

//template<class K, class V, class MAP >
template<class MAP >
static void setupMap(MAP& map, unsigned N)
{
  for (unsigned i = 0; i < N; i++)
  {
    map[i] = new unsigned(i);
  }
}


template<class MAP, class ITER >
static unsigned *find1(MAP& map,  ITER& i, unsigned key)
{
  i = map.find( key );
  return i != map.end() ? i->second : 0 ;
}

template<class MAP, class ITER >
static unsigned *find2(MAP& map,  ITER& i, unsigned key)
{
  return map[key];
}

template<class MAP, class ITER >
struct FindMapItem1
{
  inline unsigned *operator()(MAP& map,  ITER& i, unsigned key)
  {
    return map[key];
  }
};

template<class MAP, class ITER >
struct FindMapItem2
{
  inline unsigned *operator()(MAP& map,  ITER& i, unsigned key)
  {
    i = map.find( key );
    return i != map.end() ? i->second : 0 ;
    //return i->second;
  }
};

template<class MAP, class ITER, class FUNC >
static double dot1(MAP& map,  ITER& it, unsigned N, unsigned niter, FUNC& fm)
{
  //MAP ::iterator ii = 0;
  double sum=0.0;
  unsigned ll=0u;
  for (unsigned iter = 0; iter < niter; iter++)
  {
    //std::cout << "iter= " << iter << std::endl;
    for (unsigned i = 0; i < N; i++)
    {
      //unsigned *pm = find1(map, it, i);
      unsigned *pm = fm(map, it, i);
      ll+= *pm;
      sum += (double)( *pm)/(double(N));
      map[i] = pm;
    }
  }
  return sum;
}


template<class MAP, class ITER, class FUNC >
static void doTest(MAP& map,  ITER& it, unsigned N, unsigned niter, FUNC& fm, std::string msg)
{
  EXCEPTWATCH;

  double t0 =  stk::cpu_time();
  setupMap(map, N);
  double t1 =  stk::cpu_time();
  std::cout << "maptest:   setup time  = " << (t1-t0)/60. << " [min] for " << msg << std::endl;

  double t2s =  stk::cpu_time();
  double dd= dot1(map, it, N, niter, fm);
  double t2e =  stk::cpu_time();

  std::cout << "maptest:  lookup time  = " << (t2e-t2s)/60. << " [min] for " << msg << " dd= " << dd << std::endl;
}


STKUNIT_UNIT_TEST(time_maps, compare_different_maps)
{
#if DO_TEST_TIME_MAPS

  EXCEPTWATCH;

  unsigned N = 10000000; // 10M
  //unsigned N = 100000; // 100K
  //unsigned N = 10; // 10M
  unsigned niter = 10;

  unsigned init_capacity = N;

  typedef std::map<unsigned, unsigned *> std_map_type;
  std_map_type std_map1;
  std_map_type std_map2;
  //std::map<unsigned, unsigned*>::iterator std_map_it;
  std_map_type::iterator std_map_it1;
  std_map_type::iterator std_map_it2;

  typedef boost::unordered_map<unsigned, unsigned *> boost_map_type;
  boost_map_type boost_map1(init_capacity);
  boost_map_type boost_map2(init_capacity);
  boost_map_type::iterator boost_map_it1;
  boost_map_type::iterator boost_map_it2;

  //std::map<unsigned, unsigned *> std_map;

#if USE_SPARSEHASH
  typedef google::dense_hash_map<unsigned, unsigned *> google_dense_map_type;
  google_dense_map_type google_dense_map1(init_capacity);
  google_dense_map1.set_empty_key(2*N);
  google_dense_map_type::iterator google_dense_map_it1;

  {
    FindMapItem1< google_dense_map_type, google_dense_map_type::iterator > fm1;
    if (1) doTest(google_dense_map1, google_dense_map_it1, N, niter, fm1, "google_dense_map, map[key]");

  }
#endif

#if USE_SPARSEHASH
  typedef google::sparse_hash_map<unsigned, unsigned *> google_sparse_map_type;
  google_sparse_map_type google_sparse_map1(init_capacity);
  google_sparse_map_type google_sparse_map2(init_capacity);
  google_sparse_map_type::iterator google_sparse_map_it1;
  google_sparse_map_type::iterator google_sparse_map_it2;

  {
    FindMapItem1< google_sparse_map_type, google_sparse_map_type::iterator > fm1;
    if (1) doTest(google_sparse_map1, google_sparse_map_it1, N, niter, fm1, "google_sparse_map, map[key]");

  }
#endif

  {
    FindMapItem1< boost_map_type, boost_map_type::iterator > fm1;
    FindMapItem2< boost_map_type, boost_map_type::iterator > fm2;
    if (1) doTest(boost_map1, boost_map_it1, N, niter, fm1, "boost_map, map[key]");
    if (1) doTest(boost_map2, boost_map_it2, N, niter, fm2, "boost_map, find(key)");
  }

  {
    FindMapItem1< std_map_type, std_map_type::iterator > fm1;
    FindMapItem2< std_map_type, std_map_type::iterator > fm2;
    if (1) doTest(std_map1, std_map_it1, N, niter, fm1, "std_map, map[key]");
    if (1) doTest(std_map2, std_map_it2, N, niter, fm2, "std_map, find(key)");
  }

  //doTest(boost_map, N, niter);

#endif

}

}
}
}
