#ifndef stk_mesh_DiagWriter_h
#define stk_mesh_DiagWriter_h

#ifdef STK_MESH_TRACE_ENABLED

#include <stk_util/diag/Trace.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterOStream.hpp>
#include <stk_util/diag/WriterParser.hpp>

#include <stk_mesh/base/DiagWriter_fwd.hpp>
#include <stk_mesh/base/Part.hpp>

namespace stk {
namespace mesh {

class Part;
class Entity;
union EntityKey;

stk::diag::Writer &theDiagWriter();
#define meshlog stk::mesh::theDiagWriter()

class DiagWriterParser : public diag::WriterParser
{
public:
  DiagWriterParser();
};

DiagWriterParser &theDiagWriterParser();

typedef diag::Tracespec Tracespec;
typedef diag::Traceback Traceback;

/**
 * Defines the Trace objects that will be used in stk... this class
 * mostly just ensures that the correct DiagWriter is used for
 * tracing.
 */
class Trace : public diag::Trace
{
public:
  explicit Trace(const char *message)
    : diag::Trace(meshlog, message)
  {}

  Trace(const char *message, int print_mask)
    : diag::Trace(meshlog, message, print_mask)
  {}

  Trace(const char *message, int print_mask, bool do_trace)
    : diag::Trace(meshlog, message, print_mask, do_trace)
  {}
};

// If Writer does not know how to output an object you want to trace, you
// can address that here by defining an operator<< for that object. Note
// that Writer handles vectors and pointers automatically.

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Part& part);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Entity& entity);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityKey& key);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityProc& entity_proc);

} // namespace mesh
} // namespace stk

#endif

#endif // stk_mesh_DiagWriter_h
