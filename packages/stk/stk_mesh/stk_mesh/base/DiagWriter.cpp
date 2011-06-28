/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2011 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifdef STK_MESH_TRACE_ENABLED

#include <stk_util/util/Bootstrap.hpp>

#include <stk_mesh/base/DiagWriter.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

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
  mask("entity", (unsigned long) (LOG_ENTITY), "Display entity diagnostic information");
  mask("bucket", (unsigned long) (LOG_BUCKET), "Display bucket diagnostic information");
  mask("part",   (unsigned long) (LOG_PART),   "Display bucket diagnostic information");
  mask("field",  (unsigned long) (LOG_FIELD),  "Display bucket diagnostic information");
}

namespace {

void bootstrap()
{
//  diag::registerWriter("meshlog", meshlog, theDiagWriterParser());
}

std::string log_to_str(EntityModificationLog log)
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

stk::Bootstrap x(&bootstrap);

} // namespace <unnamed>

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Part& part)
{
  return writer << "Part[" << part.name() << ", " << part.mesh_meta_data_ordinal() << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Entity& entity)
{
  // Get bucket of entity
  Bucket* bucket = NULL;
  try {
    bucket = &(entity.bucket());
  }
  catch (...) {} // leave bucket as NULL if it's not found

  std::string ownership_info = "unregistered";
  std::string entity_key_str;
  EntityKey key = entity.key();
  if (bucket) {
    MetaData& meta_data = MetaData::get(*bucket);
    Part &   owned  = meta_data.locally_owned_part();
    Part &   shared = meta_data.globally_shared_part();
    if (bucket->member(owned)) {
      ownership_info = "owned";
    }
    else if (bucket->member(shared)) {
      ownership_info = "shared";
    }
    else if (bucket->size() == 0) {
      ownership_info = "marked deleted";
    }
    else {
      ownership_info = "ghosted";
    }
    entity_key_str = print_entity_key(meta_data, key);
  }
  else {
    std::ostringstream out;
    out << "(rank:" << key.rank() << ",id:" << key.id() << ")";
    entity_key_str = out.str();
  }

  writer << "Entity[key:" << entity_key_str <<
                 ", ownership:" << ownership_info <<
                 ", log:" << log_to_str(entity.log_query()) <<
                 ", owner:" << entity.owner_rank();

  // print comm info
  writer << ", COMM: ";
  PairIterEntityComm comm_itr = entity.comm();
  for ( ; !comm_itr.empty(); ++comm_itr ) {
    writer << "(ghost:" << comm_itr->ghost_id << ", proc:" << comm_itr->proc << ") ";
  }
  return writer << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityKey& key)
{
  return writer << "Entity[rank:" << key.rank() << ", id:" << key.id() << "]";
}

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityProc& entity_proc)
{
  return writer << "EntityProc[entity:" << *entity_proc.first << ", proc: " << entity_proc.second << "]";
}

} // namespace mesh
} // namespace stk

#endif // STK_MESH_TRACE_ENABLED
