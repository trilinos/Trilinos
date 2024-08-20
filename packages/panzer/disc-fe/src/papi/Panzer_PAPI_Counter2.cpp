// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_PAPI_Counter2.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <algorithm>
#include <cstring>

namespace panzer {

  int PAPICounter2::m_event_set = PAPI_NULL;
  std::vector<int> PAPICounter2::m_events;
  bool PAPICounter2::m_is_initialized = false;
  std::map<std::string,PAPICounter2::InternalCounter2> PAPICounter2::m_counters;

  /** initialize papi library */

  PAPICounter2::PAPICounter2(const std::string counter_name) :
    m_name(counter_name)
  {
    if (!m_is_initialized) {
      TEUCHOS_ASSERT( PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT );
 
      TEUCHOS_ASSERT( PAPI_create_eventset(&m_event_set) == PAPI_OK );
      
      for (std::vector<int>::const_iterator event = m_events.begin();
	   event != m_events.end(); ++event) {
	TEUCHOS_ASSERT( PAPI_add_event(m_event_set,*event) == PAPI_OK );
      }

      TEUCHOS_ASSERT( PAPI_start(m_event_set) == PAPI_OK );
      m_is_initialized = true;
    }
    
    // initialize the specific timer first time in
    std::map<std::string,InternalCounter2>::const_iterator counter;
    counter = m_counters.find(m_name);
    if (counter == m_counters.end()) {
      InternalCounter2& c = m_counters[m_name];
      c.accumulated_time = 0;
      c.start_counters.resize(m_events.size());
      c.stop_counters.resize(m_events.size());
      c.accumulated_counters.resize(m_events.size());
      c.num_calls = 0;
    }

    // mark start time
    InternalCounter2& c= m_counters[m_name];
    TEUCHOS_ASSERT( PAPI_read(m_event_set, &(c.start_counters[0])) == PAPI_OK );
    c.start_time = PAPI_get_real_usec();
    c.num_calls +=1;
  }

  PAPICounter2::~PAPICounter2()
  {
    // accumulate totals
    InternalCounter2& c= m_counters[m_name];

    TEUCHOS_ASSERT( PAPI_read(m_event_set, &(c.stop_counters[0])) == PAPI_OK );
    c.accumulated_time += PAPI_get_real_usec() - c.start_time;

    std::vector<long_long>::iterator accum = c.accumulated_counters.begin();
    std::vector<long_long>::const_iterator start = c.start_counters.begin();
    std::vector<long_long>::const_iterator stop = c.stop_counters.begin();
    for (; accum !=  c.accumulated_counters.end(); ++accum,++start,++stop)
      *accum += *stop - *start;
    
  }

  void PAPICounter2::addEventCounter(const int event)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(m_is_initialized,
			       std::logic_error,
			       "Error - cannot add event after PAPICounter is initialized!");

    m_events.push_back(event);
  }

  void PAPICounter2::startCounters()
  {
    TEUCHOS_ASSERT(PAPI_start(m_event_set) == PAPI_OK);
  }

  void PAPICounter2::stopCounters()
  {
    //TEUCHOS_ASSERT(PAPI_stop(m_event_set) == PAPI_OK);
  }

  void PAPICounter2::report(std::ostream& os, const Teuchos::Comm<int>& comm)
  {

    os << std::endl;
    os << "************************************************************" << std::endl;
    os << "* PAPI Counter Report (over all processes) " << std::endl;
    os << "************************************************************" << std::endl;

    for (std::map<std::string,InternalCounter2>::const_iterator timer = m_counters.begin();
	 timer != m_counters.end(); ++timer) {

      // Communicate totals across processes
      
      const std::vector<long long int>& accum = timer->second.accumulated_counters;
      std::vector<long long int> global_min(accum.size(),0);
      std::vector<long long int> global_max(accum.size(),0);
      std::vector<long long int> global_sum(accum.size(),0);
      std::vector<long long int> global_avg(accum.size(),0);
      long long int average_time = 0;

      Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, static_cast<int>(accum.size()), &accum[0], &global_min[0]);
      Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, static_cast<int>(accum.size()), &accum[0], &global_max[0]);
      Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, static_cast<int>(accum.size()), &accum[0], &global_sum[0]);
      Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, static_cast<int>(accum.size()), &accum[0], &global_avg[0]);

      for (std::vector<long long int>::iterator i = global_avg.begin();
	   i != global_avg.end(); ++i)
	(*i) = *i / Teuchos::as<long long int>(comm.getSize());

      Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &(timer->second.accumulated_time), &average_time);
      average_time /= Teuchos::as<long long int>(comm.getSize());

      os << timer->first<< ": Average Process Time (seconds) = " 
	 << timer->second.accumulated_time / 1.0e6 << std::endl;
      os << timer->first<< ": Number of Calls = " << timer->second.num_calls << std::endl;
      
      int i=0;
      for (std::vector<long_long>::const_iterator event=timer->second.accumulated_counters.begin();
	   event != timer->second.accumulated_counters.end(); ++event,++i) {
	char event_name[PAPI_MAX_STR_LEN];
	TEUCHOS_ASSERT( PAPI_event_code_to_name(m_events[i],event_name) == PAPI_OK);
	std::string string_event_name(event_name);
	os << timer->first << ": " << string_event_name << " = " 
	   << "min:" << global_min[i] 
	   << ", max:" << global_max[i]
	   << ", total:" << global_sum[i]
	   << ", avg:" << global_avg[i]
	   << std::endl;
      }
      
    }

    os << "************************************************************" << std::endl;
  }

}
