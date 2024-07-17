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

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"

#include "papi.h"
#include "Teuchos_DefaultComm.hpp"
#include "Panzer_PAPI_Counter2.hpp"
#include <string>

namespace panzer {


  TEUCHOS_UNIT_TEST(papi, NestedPAPICounter)
  {
    
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    panzer::PAPICounter2::addEventCounter(PAPI_TOT_INS);
    panzer::PAPICounter2::addEventCounter(PAPI_FP_OPS);
    panzer::PAPICounter2::addEventCounter(PAPI_L2_DCM);

    {
      panzer::PAPICounter2 outer_counter("Outer Counter");

      //    outer_counter.start();
      for (int i=0; i < 10; ++i) {
	
	{
	  panzer::PAPICounter2 inner_counter("Inner Loop 1");
	  
	  double a = 0;
	  for (int i=0; i < 100; ++i)
	    a += static_cast<double>(i);
	  
	}
	
	{
	  panzer::PAPICounter2 inner_counter("Inner Loop 2");
	  double a = 0;
	  for (int i=0; i < 100; ++i)
	    a += static_cast<double>(i);
	}
	
      }
      //outer_counter.stop();
      

    }  // dtor stops outer counter and records results
    Teuchos::RCP<const Teuchos::Comm<int> > t_comm = Teuchos::DefaultComm<int>::getComm();
    panzer::PAPICounter2::report(std::cout,*t_comm);
      
  }

}
