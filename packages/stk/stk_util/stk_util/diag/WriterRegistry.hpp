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

#ifndef STK_UTIL_DIAG_WriterRegistry_h
#define STK_UTIL_DIAG_WriterRegistry_h

#include "stk_util/diag/Option.hpp"               // for OptionMaskParser
#include "stk_util/util/Writer_fwd.hpp"           // for Writer
#include "stk_util/util/string_case_compare.hpp"  // for LessCase
#include <map>                                    // for map
#include <string>                                 // for string
#include <utility>                                // for pair
#include <vector>                                 // for vector

namespace stk { namespace diag { class Writer; } }
namespace stk { namespace diag { class WriterThrowSafe; } }

namespace sierra {
namespace Diag {

///
/// @addtogroup DiagWriterDetail
/// @{
///

/**
 * @brief Typedef <b>WriterRegistry</b> is a mapping from name to diagnostic writer.
 */
class WriterRegistry : public std::map<std::string, std::pair<stk::diag::Writer *, OptionMaskParser *>, stk::LessCase>
{
public:  
  WriterRegistry();
  
  ~WriterRegistry();
};

class WriterThrowSafe 
{
public:
  WriterThrowSafe();

  ~WriterThrowSafe();

private:
  std::vector<stk::diag::WriterThrowSafe *>     m_writerVector;
};
  
/**
 * @brief Function <b>getWriterRegistry</b> returns a reference to the diagnostic
 * writer registry.
 *
 * @return		a <b>WriterRegistry</b> reference to the diagnostic writer
 *			registry.
 */
WriterRegistry &getWriterRegistry();

/**
 * @brief Function <b>registerWriter</b> registers a diagnostic writer with the
 * diagnostic writer registry.
 *
 * @param name		a <b>std::string</b> const reference to the name to use for the
 *			diagnostic writer.
 *
 * @param diag_writer	a <b>Writer</b> reference to the diagnostic writer.
 *
 */
void registerWriter(const std::string &name, Writer &diag_writer, OptionMaskParser &option_parser);

///
/// @}
///

} // namespace Diag
} // namespace sierra


#endif // STK_UTIL_DIAG_WriterRegistry_h

