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

#ifndef PANZER_PAPI_COUNTER_HPP
#define PANZER_PAPI_COUNTER_HPP

#include <string>
#include <map>
#include "mpi.h"

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
