// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"

#include "krp.hpp"
#include "papi.h"
#include "Panzer_PAPI_Counter.hpp"
#include <string>

namespace panzer {

  TEUCHOS_UNIT_TEST(papi, PAPICounter)
  {
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    panzer::PAPICounter counter("Panzer Jacobian",rank,comm);

    counter.start();

    double a = 0.0;
    for (int i=0; i < 1000; ++i)
      a += static_cast<double>(i);

    counter.stop();

    counter.report(std::cout);
    
    

  }

  /*
  TEUCHOS_UNIT_TEST(papi, NestedPAPICounter)
  {
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    panzer::PAPICounter outer_counter("Outer Counter",rank,comm);
    panzer::PAPICounter inner_counter_1("Inner Counter 1",rank,comm);
    panzer::PAPICounter inner_counter_2("Inner Counter 2",rank,comm);

    outer_counter.start();
    for (int i=0; i < 1000; ++i) {

      {
	inner_counter_1.start();
	
	double a = 0;
	for (int i=0; i < 1000000; ++i)
	  a += static_cast<double>(i);
	
	inner_counter_1.stop();
      }
      
      {
	inner_counter_2.start();
	double a = 0;
	for (int i=0; i < 1000000; ++i)
	a += static_cast<double>(i);
	inner_counter_2.stop();
      }

    }
    outer_counter.stop();

    outer_counter.report(std::cout);
    
//     panzer::PAPICounter inner_counter_1("Inner Counter 1",rank,comm);
//     panzer::PAPICounter inner_counter_2("Inner Counter 2",rank,comm);
    
    inner_counter_1.report(std::cout);
    inner_counter_2.report(std::cout);

  }
  */

}
