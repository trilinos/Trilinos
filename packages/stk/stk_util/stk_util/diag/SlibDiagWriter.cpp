/*--------------------------------------------------------------------*/
/*    Copyright 2004 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/util/Bootstrap.hpp>

#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/diag/WriterRegistry.hpp>

#include <stk_util/diag/SlibDiagWriter.hpp>

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
    /* %TRACE[NONE]% */  /* %TRACE% */

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
  /* %TRACE[NONE]% */  /* %TRACE% */
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
