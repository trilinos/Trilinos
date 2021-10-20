// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ParallelErrorMessage.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino
{

ParallelErrorMessage::ParallelErrorMessage(const stk::Parallel & c) : comm(c) {}

std::pair<bool, std::string> ParallelErrorMessage::gather_message() const
{
  const std::string & local_str = local_message.str();

  std::ostringstream result;
  stk::all_write_string(comm.parallel(), result, local_str);

  int err = result.str().size();
  int global_err;
  stk::all_reduce_sum(comm.parallel(), &err, &global_err, 1);

  return make_pair(static_cast<bool>(global_err), result.str());
}
}
