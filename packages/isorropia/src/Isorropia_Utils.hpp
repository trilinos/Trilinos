//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_Utils_hpp_
#define _Isorropia_Utils_hpp_

#include <Isorropia_ConfigDefs.hpp>

namespace Isorropia {

/** Utils is the namespace within the Isorropia namespace that defines
    general definitions of utilities that may be of use to any specific 
    partitioner.
*/

namespace Utils {

/** Internal Isorropia implementation utility.
  Given a vector that specifies all processors' old or current
  offsets into a global element list, and another vector that
  specifies all processors' new or desired offsets into a global
  element list, fill a send_info and recv_info vector with data that
  can be unpacked/interpreted as follows.

  while(i<send_info.size()) {
    send_info[i]   == proc to send data to
    send_info[i+1] == starting offset of local elements to be sent
    send_info[i+2] == number of elements to send
    i += 3;
  }

  while(i<recv_info.size()) {
    recv_info[i]   == proc to recv from
    recv_info[i+1] == offset at which incoming elements will be stored
    recv_info[i+2] == number of incoming elements
  }
*/
void create_comm_plan(int myPID,
                      const std::vector<int>& all_proc_old_offsets,
                      const std::vector<int>& all_proc_new_offsets,
                      std::vector<int>& send_info,
                      std::vector<int>& recv_info);

/** Return CPU time. To measure an elapsed time, take the difference
    between two returned values.
*/
double cpu_time();
}//namespace Utils
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

