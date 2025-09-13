// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ParallelErrorMessage.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino
{

ParallelErrorMessage::ParallelErrorMessage(const stk::Parallel & c, const unsigned precision) : comm(c) {local_message.precision(precision);}

std::pair<bool, std::string> gather_messages(const stk::ParallelMachine comm, const std::string & localMsg)
{
  std::ostringstream result;
  const bool globalErr = stk::is_true_on_any_proc(comm, !localMsg.empty());
  if (globalErr)
    stk::all_write_string(comm, result, localMsg);

  return make_pair(globalErr, result.str());
}

std::pair<bool, std::string> ParallelErrorMessage::gather_message() const
{
  return gather_messages(comm.parallel(), local_message.str());
}

void RequireEmptyErrorMsg(const stk::ParallelMachine comm, const std::string & localMsg, const std::string & errorHeaderMsg)
{
  const auto & [error, msg] = gather_messages(comm, localMsg);
  STK_ThrowRequireMsg(!error, errorHeaderMsg << "\n" << msg);
}

}
