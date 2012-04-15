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

#include <Isorropia_Utils.hpp>
#include <Isorropia_Exception.hpp>

namespace Isorropia {

namespace Utils {

void create_comm_plan(int myPID,
                      const std::vector<int>& all_proc_old_offsets,
                      const std::vector<int>& all_proc_new_offsets,
                      std::vector<int>& send_info,
                      std::vector<int>& recv_info)
{
  send_info.clear(); recv_info.clear();

  //The map elements that this processor currently holds are
  //given by the range:
  //      [all_proc_old_offsets[myPID]...all_proc_old_offsets[myPID+1]-1]
  //and the map elements that this processor needs to hold are
  //      [all_proc_new_offsets[myPID]...all_proc_new_offsets[myPID+1]-1]
  //
  //So now we need to figure out which elements we need to send/recv
  //to/from neighboring processors.

  int my_old_start = all_proc_old_offsets[myPID];
  int my_old_end = all_proc_old_offsets[myPID+1]-1;

  int p, numProcs = all_proc_old_offsets.size()-1;

  for(p=0; p<numProcs; ++p) {
    if (p==myPID) continue;

    int start_p = all_proc_new_offsets[p];
    int end_p = all_proc_new_offsets[p+1]-1;

    if (start_p <= my_old_start && my_old_start <= end_p) {
      //if my_old_start contained in [start_p .. end_p]

      //we'll send to proc p
      send_info.push_back(p);

      //we'll send elements starting from this position
      int send_offset = 0;
      send_info.push_back(send_offset);

      int num_send = end_p <= my_old_end ?
         end_p-my_old_start + 1 : my_old_end-my_old_start + 1;
      send_info.push_back(num_send);
    }
    else if (start_p <= my_old_end && my_old_end <= end_p) {
      //if my_old_end contained in [start_p .. end_p]

      //we'll send to proc p
      send_info.push_back(p);

      //we'll send elements starting from this position
      int send_offset = start_p - my_old_start;
      send_info.push_back(send_offset);

      int num_send = my_old_end - start_p + 1;
      send_info.push_back(num_send);
    }
    else if (my_old_start <= start_p && start_p <= my_old_end) {
      //if start_p contained in [my_old_start .. my_old_end]
      //we'll send to proc p
      send_info.push_back(p);
      int send_offset = start_p - my_old_start;
      send_info.push_back(send_offset);

      int num_send = my_old_end <= end_p ?
        my_old_end - start_p + 1 : end_p - start_p + 1;
      send_info.push_back(num_send);
    }
  }

  int my_new_start = all_proc_new_offsets[myPID];
  int my_new_end = all_proc_new_offsets[myPID+1]-1;
  for(p=0; p<numProcs; ++p) {
    if (p==myPID) continue;

    int old_start_p = all_proc_old_offsets[p];
    int old_end_p = all_proc_old_offsets[p+1]-1;

    if (old_start_p <= my_new_start && my_new_start <= old_end_p) {
      //if my_new_start contained in [old_start_p .. old_end_p]

      //we'll recv from proc p
      recv_info.push_back(p);

      //recv'd data will be placed at this offset in our new element list:
      int recv_offset = 0;
      recv_info.push_back(recv_offset);

      int num_recv = old_end_p <= my_new_end ?
         old_end_p - my_new_start + 1 : my_new_end - my_new_start + 1;
      recv_info.push_back(num_recv);
    }
    else if (old_start_p <= my_new_end && my_new_end <= old_end_p) {
      //if my_new_end contained in [old_start_p .. old_end_p]

      //we'll recv from proc p
      recv_info.push_back(p);

      //recv'd data will be placed at this offset in our new element list:
      int recv_offset = old_start_p - my_new_start;
      recv_info.push_back(recv_offset);

      int num_recv = my_new_end - old_start_p + 1;
      recv_info.push_back(num_recv);
    }
    else if (my_new_start <= old_start_p && old_start_p <= my_new_end) {
      //if old_start_p contained in [my_new_start .. my_new_end]

      //we'll recv from proc p
      recv_info.push_back(p);

      //recv'd data will be placed at this offset in our new element list:
      int recv_offset = old_start_p - my_new_start;
      recv_info.push_back(recv_offset);

      int num_recv = old_end_p <= my_new_end ?
         old_end_p - old_start_p + 1 : my_new_end - old_start_p + 1;
      recv_info.push_back(num_recv);
    }
  }
}

//----------------------------------------------------------------------------
double cpu_time()
{
  double cpu_seconds = 0.0;

#ifdef HAVE_TIME_H
  cpu_seconds = (1.0*clock())/CLOCKS_PER_SEC;
#endif

  return(cpu_seconds);
}

}//namespace Utils
}//namespace Isorropia

