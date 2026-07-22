// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PAPI_COUNTER_2_HPP
#define PANZER_PAPI_COUNTER_2_HPP

#include <string>
#include <map>
#include <vector>
#include "papi.h"
#include "Teuchos_Comm.hpp"
#include "Teuchos_SerializationTraits.hpp"


// namespace Teuchos {
//   template<typename Ordinal>
//   class SerializationTraits<Ordinal,long_long>
//     : public DirectSerializationTraits<Ordinal,long_long>
//   {};
// }

namespace panzer {


  /** \brief Interface to papi counters

      Mimics the Teuchos::TimeMonitor functionality
  */ 
  class PAPICounter2 {

  public:
    
    /** Timer starts at construction.  Stops when destructor called. */
    PAPICounter2(const std::string);

    /** Stops timer */
    ~PAPICounter2();

    /** Add PAPI events. Can only be called before the initializePAPI() method. */
    static void addEventCounter(const int event);

    static void report(std::ostream& os, const Teuchos::Comm<int>& comm);

    static void startCounters();

    static void stopCounters();

  private:

    struct InternalCounter2 {
      long_long start_time;
      long_long accumulated_time;

      std::vector<long_long> start_counters;
      std::vector<long_long> stop_counters;
      std::vector<long_long> accumulated_counters;

      long_long num_calls;
    };

    //! PAPI event set
    static int m_event_set;
    //! papi event index
    static std::vector<int> m_events;
    //! true if the static members have been intitialized
    static bool m_is_initialized;
    //! maps the counter name to the data object
    static std::map<std::string,InternalCounter2> m_counters;
    //! name of this counter
    std::string m_name;

  };

}

#endif
