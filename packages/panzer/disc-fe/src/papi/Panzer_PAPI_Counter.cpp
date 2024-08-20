// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_PAPI_Counter.hpp"

#include "krp.hpp"
namespace panzer {

  std::map<std::string,PAPICounter::InternalCounter> PAPICounter::m_counters;

  PAPICounter::PAPICounter(const std::string name, const int my_rank, MPI_Comm comm)
    : m_name(name), m_rank(my_rank), m_comm(comm)
  {

  }

  void PAPICounter::start()
  {
    InternalCounter& c = m_counters[m_name];

    panzer::krp_init_(&m_rank,&c.hw_counters,&c.rcy,&c.rus,&c.ucy,&c.uus);
  }

  void PAPICounter::stop()
  {
    InternalCounter& c = m_counters[m_name];
    
    //panzer::krp_rpt_init_sum_(&m_rank,m_comm,&c.hw_counters,&c.rcy,&c.rus,&c.ucy,&c.uus,&c.rt_rus,&c.rt_ins,&c.rt_fp,&c.rt_dcm,&c.uus);
    panzer::krp_rpt_init_sum_(&m_rank,m_comm,&c.hw_counters,&c.rcy,&c.rus,&c.ucy,&c.uus,&c.rt_rus,&c.rt_ins,&c.rt_fp,&c.rt_dcm,const_cast<char*>(m_name.c_str()));
  }

  void PAPICounter::report(std::ostream& os)
  {
//     InternalCounter& c = m_counters[m_name];
    
//     panzer::krp_rpt_(&m_rank,m_comm,&c.hw_counters,&c.rcy,&c.rus,&c.ucy,&c.uus,const_cast<char*>(m_name.c_str()));
  }

}
