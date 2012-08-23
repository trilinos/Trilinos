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

#include "krp.hpp"
#include "papi.h"
#include "Panzer_PAPI_Counter.hpp"
#include <string>

namespace panzer {

  TEUCHOS_UNIT_TEST(papi, krp_c_fcn_calls)
  {
    int iam;
    MPI_Comm_rank(MPI_COMM_WORLD, &iam);
    int hw_counters;
    long long int rcy=0;
    long long int rus=0;
    long long int ucy=0;
    long long int uus=0;
    
    MPI_Comm comm = MPI_COMM_WORLD;
    long int rt_rus=0;
    long int rt_ins=0;
    long int rt_fp=0;
    long int rt_dcm=0; 
    
    std::string name = "TESTING RAW INTERFACE\0";

    panzer::krp_init_(&iam,&hw_counters,&rcy,&rus,&ucy,&uus);
    
    double a = 0.0;
    for (int i=0; i < 1000; ++i)
      a += static_cast<double>(i);

    out << "a = "  << a << std::endl;

    //panzer::krp_init_sum_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,&rt_rus,&rt_ins,&rt_fp,&rt_dcm);

    //panzer::krp_rpt_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,const_cast<char*>(name.c_str()));

    //panzer::krp_rpt_init_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,const_cast<char*>(name.c_str()));

    panzer::krp_rpt_init_sum_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,&rt_rus,&rt_ins,&rt_fp,&rt_dcm,const_cast<char*>(name.c_str()));

    PAPI_reset(hw_counters);

  }

}
