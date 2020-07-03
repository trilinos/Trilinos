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

#include "FileValidator.hpp"

#include <fstream>

#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace balance {

FileValidator::FileValidator(MPI_Comm comm)
  : m_comm(comm),
  m_isSerial(stk::parallel_machine_size(m_comm) == 1)
{ }

void FileValidator::require_file_exists(const std::string& filename) const
{
  ThrowRequireMsg(does_file_exist(filename), "Input file does not exist.\n");
}

bool FileValidator::serial_input_equals_output(const std::string& infile, const std::string& outfile) const
{
  return m_isSerial && (trim_filename(infile) == trim_filename(outfile));
}

bool FileValidator::does_file_exist(const std::string& filename) const
{
  bool exists = true;
  if (!std::ifstream(filename)) {
    exists = false;
  }
  return exists;
}

std::string FileValidator::trim_filename(const std::string& filename) const
{
  return (filename.substr(0,2) == "./") ? filename.substr(2) : filename;
}

} }
