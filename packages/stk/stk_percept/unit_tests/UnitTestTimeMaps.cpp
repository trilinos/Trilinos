/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#define DO_TEST_TIME_MAPS 0

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <mpi.h>
#include <math.h>


#include <boost/unordered_map.hpp>

#include <map>

#include <tr1/unordered_map>


#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>

namespace stk
{
  namespace percept
  {
    namespace unit_tests
    {

#define EXTRA_PRINT 0

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

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
        
        double t0 =  stk::wall_time(); 
        setupMap(map, N);
        double t1 =  stk::wall_time(); 
        std::cout << "maptest: setup time = " << (t1-t0)/60. << " [min]" << std::endl;

        double t2s =  stk::wall_time(); 
        double dd= dot1(map, it, N, niter, fm);
        double t2e =  stk::wall_time(); 

        std::cout << "maptest: dotest time, " << msg << " = " << (t2e-t2s)/60. << " [min] dd= " << dd << std::endl;
      }

      TEST(time_maps, time_boost)
      {
        EXCEPTWATCH;
        if (!DO_TEST_TIME_MAPS) return;
        typedef std::map<unsigned, unsigned *> std_map_type;
        std_map_type std_map1;
        std_map_type std_map2;
        //std::map<unsigned, unsigned*>::iterator std_map_it;
        std_map_type::iterator std_map_it1;
        std_map_type::iterator std_map_it2;

        
        typedef boost::unordered_map<unsigned, unsigned *> boost_map_type;
        boost_map_type boost_map1;
        boost_map_type boost_map2;
        boost_map_type::iterator boost_map_it1;
        boost_map_type::iterator boost_map_it2;

        //std::map<unsigned, unsigned *> std_map;

        unsigned N = 10000000; // 10M
        //unsigned N = 100000; // 100K
        //unsigned N = 10; // 10M
        unsigned niter = 10;

        {
          FindMapItem1< std_map_type, std_map_type::iterator > fm1;
          FindMapItem2< std_map_type, std_map_type::iterator > fm2;
          if (1) doTest(std_map1, std_map_it1, N, niter, fm1, "std_map, map[key]");
          if (1) doTest(std_map2, std_map_it2, N, niter, fm2, "std_map, find(key)");
        }

        {
          FindMapItem1< boost_map_type, boost_map_type::iterator > fm1;
          FindMapItem2< boost_map_type, boost_map_type::iterator > fm2;
          if (1) doTest(boost_map1, boost_map_it1, N, niter, fm1, "boost_map, map[key]");
          if (1) doTest(boost_map2, boost_map_it2, N, niter, fm2, "boost_map, find(key)");
        }

        //doTest(boost_map, N, niter);

      }

    }
  }
}
