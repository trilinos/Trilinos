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
