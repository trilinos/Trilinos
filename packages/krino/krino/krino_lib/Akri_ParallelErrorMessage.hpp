// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_PARALLELERRORMESSAGE_H_
#define KRINO_INCLUDE_AKRI_PARALLELERRORMESSAGE_H_

#include <stk_util/parallel/Parallel.hpp>
#include <sstream>
#include <string>

namespace krino
{

class ParallelErrorMessage
{
public:
  ParallelErrorMessage(const stk::Parallel & c);

  template <typename T> ParallelErrorMessage & operator<<(const T & t)
  {
    local_message << t;
    return *this;
  }

  // This will gather the messages from all ranks to rank 0 where it can be output.
  // and additionally return a parallel-consistent bool for whether or not any errors
  // were logged.
  std::pair<bool, std::string> gather_message() const;

  bool have_local_error() const { return !local_message.str().empty(); }

private:
  const stk::Parallel comm;
  std::ostringstream local_message;
};
}

#endif /* KRINO_INCLUDE_AKRI_PARALLELERRORMESSAGE_H_ */
