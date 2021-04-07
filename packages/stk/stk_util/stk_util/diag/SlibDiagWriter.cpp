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

#include "stk_util/diag/SlibDiagWriter.hpp"
#include "stk_util/diag/WriterParser.hpp"      // for WriterParser
#include "stk_util/diag/WriterRegistry.hpp"    // for registerWriter
#include "stk_util/environment/OutputLog.hpp"  // for dwout
#include "stk_util/util/Bootstrap.hpp"         // for Bootstrap
#include "stk_util/util/Writer.hpp"            // for Writer
#include "stk_util/util/Writer_fwd.hpp"        // for PrintMask, LOG_GLOBAL_VARIABLE, LOG_MEMORY
#include <cstdlib>                             // for getenv
#include <iosfwd>                              // for ostream



namespace sierra {
namespace Slib {

/**
 * Class <code>Parser</code> defines a mapping between strings and bit masks.
 *
 * After populating a Parser object, Diag::WriterParser::setPrintMask() will parse the
 * input string and return the corresponding print mask.
 *
 */
class DiagWriterParser : public Diag::WriterParser
{
public:
  /**
   * Creates a new <code>DiagWriterParser</code> instance.
   *
   */
  DiagWriterParser() {
    mask("resources", (Diag::PrintMask) (LOG_RESOURCE), "Display resource assignments");
    mask("plugins", (Diag::PrintMask) (LOG_PLUGIN), "Display plugin information");
    mask("global-variables", (Diag::PrintMask) (LOG_GLOBAL_VARIABLE), "Display global variable operations");
    mask("memory", (Diag::PrintMask) (LOG_MEMORY), "Display platform specific memory usage information");
  }
};

DiagWriterParser &
theDiagWriterParser()
{
  static DiagWriterParser parser;

  return parser;
}

stk::diag::Writer &
theDiagWriter()
{
  static stk::diag::Writer s_diagWriter(sierra::dwout().rdbuf(), theDiagWriterParser().parse(std::getenv("SIERRA_SLIBOUT")));
  return s_diagWriter;
}


namespace {

void bootstrap()
{
  Diag::registerWriter("slibout", slibout, theDiagWriterParser());
}

stk::Bootstrap x(&bootstrap);

} // namespace <unnamed>

} // namespace Slib
} // namespace sierra
