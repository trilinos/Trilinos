// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/Null_Streambuf.hpp>
#include <iostream>
#include <sstream>
#include <memory>

namespace stk {

struct OutputStreams
{
  OutputStreams(ParallelMachine comm);
  ~OutputStreams();
  static OutputStreams& instance(ParallelMachine comm = MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD

  null_streambuf m_nullBuf;
  std::ostream   m_outputNull;
  std::ostream*  m_outputP0;
  std::ostringstream m_output;
  ParallelMachine m_comm;
};

static std::unique_ptr<OutputStreams> s_staticOutputStreams;

OutputStreams::OutputStreams(ParallelMachine comm)
: m_nullBuf(),
  m_outputNull(&m_nullBuf),
  m_outputP0(&std::cout),
  m_output(),
  m_comm(comm)
{
  if (stk::parallel_machine_rank(m_comm) != 0) {
    m_outputP0 = &m_outputNull;
  }
}

OutputStreams::~OutputStreams()
{
}

OutputStreams& OutputStreams::instance(ParallelMachine comm)
{
  if (s_staticOutputStreams == nullptr) {
    s_staticOutputStreams = std::make_unique<OutputStreams>(comm);
  }
  return *s_staticOutputStreams;
}

std::ostream& outputP0()
{
  return *OutputStreams::instance().m_outputP0;
}

std::ostream& outputNull()
{
  return OutputStreams::instance().m_outputNull;
}

std::ostream& output()
{
  return OutputStreams::instance().m_output;
}

void output_flush()
{
  std::ostringstream& oss = OutputStreams::instance().m_output;
  stk::all_write_string(OutputStreams::instance().m_comm, outputP0(), oss.str());
  oss.str("");
}

void set_outputP0(std::ostream* ostreamPtr, ParallelMachine comm)
{
  if (comm != OutputStreams::instance().m_comm) {
    reset_default_output_streams(comm);
  }
  if (stk::parallel_machine_rank(comm) == 0) {
    OutputStreams::instance().m_outputP0 = ostreamPtr;
  }
}

void reset_default_output_streams(ParallelMachine comm)
{
  s_staticOutputStreams = std::make_unique<OutputStreams>(comm);
}

}

