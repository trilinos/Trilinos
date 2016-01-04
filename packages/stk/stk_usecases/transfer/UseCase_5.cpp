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

#include <Intrepid_FieldContainer.hpp>
#include <boost/shared_ptr.hpp>



#include <stk_mesh/base/Comm.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include "MDMesh.hpp"
#include "LinearInterpolate.hpp"
#include <stk_transfer/GeometricTransfer.hpp>



bool use_case_5_driver(stk::ParallelMachine  comm)
{
  stk::diag::Timer timer("Transfer Use Case 5",
                          use_case::TIMER_TRANSFER,
                          use_case::timer());
  stk::diag::Timer timer_point_to_point(" Point To Point", timer);
  use_case::timerSet().setEnabledTimerMask(use_case::TIMER_ALL);

  bool status = true;

  const double TOLERANCE = 0.000001;
  const double  rand_max = RAND_MAX;
  enum {           DIM = 3  };
  enum { FROMNUMPOINTS = 100  };
  enum {   TONUMPOINTS = 100  };

  typedef Intrepid::FieldContainer<double>  MDArray;

  MDArray FromPoints (FROMNUMPOINTS,DIM),
          ToPoints   (  TONUMPOINTS,DIM),
          FromValues (FROMNUMPOINTS,  2),
          ToValues   (  TONUMPOINTS,  2);
  for (unsigned i=0 ; i<FROMNUMPOINTS; ++i) {
    double l=0, q=0;
    for (unsigned j=0 ; j<DIM; ++j) {
      FromPoints(i,j) = rand()/rand_max;
      l +=   FromPoints(i,j);
      q += j*FromPoints(i,j);
    }
    FromValues(i,0) = l;
    FromValues(i,1) = q;
  }
  for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
    for (unsigned j=0 ; j<DIM; ++j) {
      ToPoints(i,j) = rand()/rand_max;
    }
  }

  const double initial_radius   = .05;
  boost::shared_ptr<stk::transfer::MDMesh >
    transfer_domain_mesh (new stk::transfer::MDMesh(FromValues, FromPoints, initial_radius, comm));
  boost::shared_ptr<stk::transfer::MDMesh >
    transfer_range_mesh  (new stk::transfer::MDMesh(  ToValues, ToPoints,   initial_radius, comm));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::MDMesh,
      class stk::transfer::MDMesh
    >
  >
  transfer(transfer_domain_mesh, transfer_range_mesh, "STK Transfer test Use case 5");
  {
    stk::diag::TimeBlock __timer_point_to_point(timer_point_to_point);
    try {
      transfer.initialize();
      transfer.apply();
    } catch (std::exception &e) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an std::exception with what string:"
                <<e.what()
                <<"      rethrowing....."
                <<std::endl;
      status = status && false;
    } catch (...) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an exception, rethrowing..."
                <<std::endl;
      status = status && false;
    }
  }

  if (status) {
    bool success = true;
    for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
      double check_l = 0;
      double check_q = 0;
      for (unsigned j=0 ; j<DIM; ++j) {
        check_l +=   ToPoints(i,j);
        check_q += j*ToPoints(i,j);
      }
      if (TOLERANCE < fabs(check_l-ToValues(i,0)) ) {
        std::cout <<__FILE__<<":"<<__LINE__
                  <<" ToValues:"<<ToValues(i,0)
                  <<" check:"<<check_l
                  <<" error:"<<fabs(check_l-ToValues(i,0))
                  <<std::endl;
        success = false;
      }
      if (TOLERANCE < fabs(check_q-ToValues(i,1)) ) {
        std::cout <<__FILE__<<":"<<__LINE__
                  <<" ToValues:"<<ToValues(i,1)
                  <<" check:"<<check_q
                  <<" error:"<<fabs(check_q-ToValues(i,1))
                  <<std::endl;
        success = false;
      }
    }
    status = status && success;
  }
  timer.stop();

  const bool collective_result = use_case::print_status(comm, status);
  return collective_result;
}
