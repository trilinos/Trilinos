// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_UnitTestFixture_hpp
#define percept_UnitTestFixture_hpp

#include <string>
#include <iostream>
#include <cmath>

#include <stk_util/parallel/Parallel.hpp>


#include <stk_util/util/Writer.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/RunEnvironment.hpp>


  namespace percept {

#define EXPECT_DOUBLE_EQ_APPROX_TOL(expected, actual, tol)                  \
      {                                                                 \
        if (std::fabs(expected-actual) > 0.5*(std::fabs(expected)+std::fabs(actual))*tol) \
          {                                                             \
            std::cout << "std::fabs(expected-actual) = "                \
                      << std::fabs(expected-actual)                     \
                      << " 0.5*(std::fabs(expected)+std::fabs(actual))*tol = " \
                      << 0.5*(std::fabs(expected)+std::fabs(actual))*tol \
                      << " tol= " << tol << std::endl;                  \
            EXPECT_DOUBLE_EQ(expected, actual);                         \
          }                                                             \
      }

#define EXPECT_DOUBLE_EQ_APPROX(expected, actual) EXPECT_DOUBLE_EQ_APPROX_TOL(expected, actual, 1.e-6)


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


      typedef RunEnvironment Fixture;

  }

#endif
