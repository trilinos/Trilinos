#ifndef stk_percept_UnitTestFixture_hpp
#define stk_percept_UnitTestFixture_hpp

#include <string>
#include <iostream>
#include <cmath>

#include <stk_util/parallel/Parallel.hpp>


#include <stk_util/diag/Writer.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>


namespace stk { 
  namespace percept { 

#define STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(expected, actual, tol)                  \
      {                                                                 \
        if (std::fabs(expected-actual) > 0.5*(std::fabs(expected)+std::fabs(actual))*tol) \
          {                                                             \
            std::cout << "std::fabs(expected-actual) = "                \
                      << std::fabs(expected-actual)                     \
                      << " 0.5*(std::fabs(expected)+std::fabs(actual))*tol = " \
                      << 0.5*(std::fabs(expected)+std::fabs(actual))*tol \
                      << " tol= " << tol << std::endl;                  \
            STKUNIT_EXPECT_DOUBLE_EQ(expected, actual);                         \
          }                                                             \
      }

#define STKUNIT_EXPECT_DOUBLE_EQ_APPROX(expected, actual) STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(expected, actual, 1.e-6)


#define TIME_IT(expr,total_time)                \
      {                                         \
        double tstart =  stk::wall_time();      \
        {                                       \
          expr                                  \
            }                                   \
        total_time = stk::wall_dtime(tstart);   \
      }


      //using namespace use_case;
      enum { 
        LOG_NORM              = 2*LOG_APPLICATION,
        LOG_GEOMETRY_VERIFIER = 2*LOG_NORM,
        LOG_MESH_COLORER      = 2*LOG_GEOMETRY_VERIFIER
        
      };


#if 0
      class Fixture : public RunEnvironment
      {
      public:
        static void setup_read_mesh_create_coordsMag_field(PerceptMesh& meshUtil, bool create_field, 
                                                           const std::string file="./cube_hex8.e", bool print=false);

      };
#endif
      typedef RunEnvironment Fixture;

  }
}

#endif
