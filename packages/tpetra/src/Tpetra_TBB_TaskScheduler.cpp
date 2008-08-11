// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Tpetra_TBB_TaskScheduler.hpp"

#ifdef HAVE_TPETRA_TBB

#include <tbb/task_scheduler_init.h>

namespace Tpetra {

static tbb::task_scheduler_init* tpetra_tbb_task_scheduler = NULL;

//--------------------------------
bool task_scheduler_is_alive()
{
  return( tpetra_tbb_task_scheduler != NULL );
}

//--------------------------------
void task_scheduler(TASK_SCHEDULER_ACTION action, int num_threads)
{
  if (num_threads == 0) {
    num_threads = tbb::task_scheduler_init::automatic;
  }

  static int num_start_requests = 0;

  static int num_stop_requests = 0;

  if (action == START_TASK_SCHEDULER) {
    if (tpetra_tbb_task_scheduler == NULL) {
      tpetra_tbb_task_scheduler = new tbb::task_scheduler_init(tbb::task_scheduler_init::deferred);
      num_start_requests = 0;
      num_stop_requests = 0;
    }

    if (num_start_requests == 0) {
      std::cout << "!initializing task_scheduler!"<<std::endl;
      tpetra_tbb_task_scheduler->initialize(num_threads);
    }

    ++num_start_requests;
    std::cout << "num_start_requests="<<num_start_requests<<std::endl;
  }

  if (action == STOP_TASK_SCHEDULER) {
    if (num_stop_requests == num_start_requests-1) {
      std::cout << "!destroying task_scheduler!"<<std::endl;
      delete tpetra_tbb_task_scheduler;
      tpetra_tbb_task_scheduler = NULL;
    }

    ++num_stop_requests;
    std::cout << "num_stop_requests="<<num_stop_requests<<std::endl;
  }
}

}//namespace Tpetra

#endif

