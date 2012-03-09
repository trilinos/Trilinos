
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"

#include "krp.hpp"

#include <string>

namespace panzer {

  TEUCHOS_UNIT_TEST(papi, krp_c_fcn_calls)
  {
    int iam;
    MPI_Comm_rank(MPI_COMM_WORLD, &iam);
    int hw_counters;
    long long int rcy;
    long long int rus;
    long long int ucy;
    long long int uus;
    
    MPI_Comm comm = MPI_COMM_WORLD;
    long int rt_rus;
    long int rt_ins;
    long int rt_fp;
    long int rt_dcm; 
    
    std::string name = "ROGER1\0";

    panzer::krp_init_(&iam,&hw_counters,&rcy,&rus,&ucy,&uus);
    
    double a = 0;
    for (int i=0; i < 1000000; ++i)
      a += static_cast<double>(i);

    panzer::krp_init_sum_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,&rt_rus,&rt_ins,&rt_fp,&rt_dcm);

    panzer::krp_rpt_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,const_cast<char*>(name.c_str()));

    panzer::krp_rpt_init_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,const_cast<char*>(name.c_str()));

    panzer::krp_rpt_init_sum_(&iam,comm,&hw_counters,&rcy,&rus,&ucy,&uus,&rt_rus,&rt_ins,&rt_fp,&rt_dcm,const_cast<char*>(name.c_str()));

  }

}
