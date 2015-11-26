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
