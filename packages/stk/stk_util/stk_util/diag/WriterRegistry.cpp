// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_util/diag/Option.hpp>     // for OptionMaskParser
#include "stk_util/util/Writer.hpp"     // for WriterThrowSafe
#include "stk_util/util/Writer_fwd.hpp"  // for Writer
#include "stk_util/util/string_case_compare.hpp"  // for LessCase

namespace sierra {
namespace Diag {

WriterThrowSafe::WriterThrowSafe()
{
  for (Diag::WriterRegistry::iterator it = Diag::getWriterRegistry().begin(); it != Diag::getWriterRegistry().end(); ++it)
    m_writerVector.push_back(new stk::diag::WriterThrowSafe(*(*it).second.first));
}


WriterThrowSafe::~WriterThrowSafe()
{
  for (std::vector<stk::diag::WriterThrowSafe *>::iterator it = m_writerVector.begin(); it != m_writerVector.end(); ++it)
    delete (*it);
}


WriterRegistry::WriterRegistry()
{}


WriterRegistry::~WriterRegistry()
{}


WriterRegistry &
getWriterRegistry()
{
  static WriterRegistry s_writerRegistry;

  return s_writerRegistry;
}


void
registerWriter(
  const std::string &	name,
  Writer &		diag_writer,
  OptionMaskParser &    option_parser)
{
  getWriterRegistry().insert(std::make_pair(name, std::make_pair(&diag_writer, &option_parser)));
}


void
unregisterWriter(
  const std::string &	name,
  Writer &		writer)
{
  WriterRegistry::iterator it = getWriterRegistry().find(name);
  if (it != getWriterRegistry().end() && (*it).second.first == &writer)
    getWriterRegistry().erase(it);
}

} // namespace Diag
} // namespace sierra
