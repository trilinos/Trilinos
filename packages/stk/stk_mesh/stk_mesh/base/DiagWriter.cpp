/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2011 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_mesh/base/DiagWriter.hpp>
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Types.hpp>      // for EntityProc, EntityState
#include <string>                       // for string
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_util/diag/Writer.hpp"     // for operator<<, Writer
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg


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


#endif

std::string log_to_str(EntityState log)
{
  if (log == 0) {
    return "Not changed";
  }
  else if (log == 1) {
    return "Created";
  }
  else if (log == 2) {
    return "Modified";
  }
  else if (log == 3) {
    return "Marked deleted";
  }
  else {
    ThrowRequireMsg(false, "Unknown log " << log);
  }
  return "";
}

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

} // namespace mesh
} // namespace stk

int dummy_DiagWriter()
{
  // This function is present just to put a symbol in the object
  // file and eliminate a "empty object file" warning on the mac...
  return 1;
}
