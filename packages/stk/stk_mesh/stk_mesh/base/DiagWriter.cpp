// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stk_mesh/base/DiagWriter.hpp>
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Types.hpp>      // for EntityProc, EntityState
#include <string>                       // for string
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg
#include "stk_mesh/baseImpl/Partition.hpp"

#ifdef STK_MESH_TRACE_ENABLED
#include <stk_util/util/Bootstrap.hpp>
#endif

namespace stk {
namespace mesh {

#ifdef STK_MESH_TRACE_ENABLED

namespace {

static stk::diag::Writer* s_diagWriter = NULL;

}

void initDiagWriter(std::ostream& stream)
{
  s_diagWriter = new stk::diag::Writer(stream.rdbuf(),
                                       theDiagWriterParser().parse(std::getenv("MESHLOG")));
}

stk::diag::Writer & theDiagWriter()
{
  ThrowRequireMsg(s_diagWriter != NULL, "Please call initDiagWwriter before theDiagWriter");
  return *s_diagWriter;
}

DiagWriterParser & theDiagWriterParser()
{
  static DiagWriterParser parser;

  return parser;
}

DiagWriterParser::DiagWriterParser()
  : stk::diag::WriterParser()
{
  mask("entity",       static_cast<unsigned long>(LOG_ENTITY),       "Display entity diagnostic information");
  mask("bucket",       static_cast<unsigned long>(LOG_BUCKET),       "Display bucket diagnostic information");
  mask("part",         static_cast<unsigned long>(LOG_PART),         "Display part diagnostic information");
  mask("field",        static_cast<unsigned long>(LOG_FIELD),        "Display field diagnostic information");
  mask("partition",    static_cast<unsigned long>(LOG_PARTITION),    "Display partition diagnostic information");
  mask("connectivity", static_cast<unsigned long>(LOG_CONNECTIVITY), "Display connectivity diagnostic information");
}

namespace {

void bootstrap()
{
//  diag::registerWriter("meshlog", meshlog, theDiagWriterParser());
}

stk::Bootstrap x(&bootstrap);

} // namespace <unnamed>

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Part& part)
{
  return writer << "Part[" << part.name() << ", " << part.mesh_meta_data_ordinal() << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityKey& key)
{
  return writer << "Entity[rank:" << key.rank() << ", id:" << key.id() << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityProc& entity_proc)
{
  return writer << "EntityProc[entity:" << entity_proc.first.local_offset() << ", proc: " << entity_proc.second << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Bucket& bucket)
{
  std::ostringstream out;
  out << bucket;
  return writer << out.str();
}

namespace impl {

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Partition& partition)
{
  std::ostringstream out;
  out << partition;
  return writer << out.str();
}

}

#endif

} // namespace mesh
} // namespace stk

int dummy_DiagWriter()
{
  // This function is present just to put a symbol in the object
  // file and eliminate a "empty object file" warning on the mac...
  return 1;
}
