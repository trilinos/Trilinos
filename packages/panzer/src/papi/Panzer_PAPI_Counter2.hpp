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
