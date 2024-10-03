// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PAPI_COUNTER_HPP
#define PANZER_PAPI_COUNTER_HPP

#include <string>
#include <map>
#include <mpi.h>

namespace panzer {

  class PAPICounter {

    struct InternalCounter {
      int hw_counters;
      long long int rcy;
      long long int rus;
      long long int ucy;
      long long int uus;
      
      long int rt_rus;
      long int rt_ins;
      long int rt_fp;
      long int rt_dcm;

      InternalCounter() : 
	hw_counters(0),
	rcy(0),
	rus(0),
	ucy(0),
	uus(0),
	rt_rus(0),
	rt_ins(0),
	rt_fp(0),
	rt_dcm(0)
      { }

    };
    
  public:
    
    PAPICounter(const std::string, const int my_rank, MPI_Comm comm);

    void start();

    void stop();

    void report(std::ostream& os);

  private:

    //! dangerous in a multithreaded world!
    static std::map<std::string,InternalCounter> m_counters;
    std::string m_name;      
    int m_rank;
    MPI_Comm m_comm;

  };

}

#endif
