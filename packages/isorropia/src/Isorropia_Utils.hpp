//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
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

