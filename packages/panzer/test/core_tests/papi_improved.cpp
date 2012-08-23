// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

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
