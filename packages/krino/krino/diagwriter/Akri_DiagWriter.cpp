// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_util/util/Bootstrap.hpp>

#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/diag/WriterRegistry.hpp>

#include <Akri_DiagWriter.hpp>

namespace krino {

DiagWriterParser &
theDiagWriterParser()
{
  /* %TRACE% */  /* %TRACE% */
  static DiagWriterParser parser;
  return parser;
}


stk::diag::Writer &
theDiagWriter()
{
  /* %TRACE[NONE]% */  /* %TRACE% */
  static stk::diag::Writer s_diagWriter(sierra::dwout().rdbuf(), theDiagWriterParser().parse(std::getenv("KRINOLOG")));

  return s_diagWriter;
}

DiagWriterParser::DiagWriterParser()
    : stk::diag::WriterParser()
{
  /* %TRACE% */  /* %TRACE% */
  mask("debug",      (unsigned long) (LOG_DEBUG),      "Display debug diagnostic information");
  mask("subelement", (unsigned long) (LOG_SUBELEMENT), "Display subelement decomposition diagnostic information");
  mask("facets",     (unsigned long) (LOG_FACETS),     "Output exodus file with facets data");
  mask("parts",      (unsigned long) (LOG_PARTS),      "Display CDFEM parts diagnostic information");
}

namespace
{

void bootstrap() { sierra::Diag::registerWriter("krinolog", krinolog, theDiagWriterParser()); }

stk::Bootstrap x(&bootstrap);
} // namespace <unnamed>

} // namespace krino
