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

#ifndef stk_mesh_DiagWriter_h
#define stk_mesh_DiagWriter_h

#include <sstream>

#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey

#ifdef STK_MESH_TRACE_ENABLED
#include <stk_mesh/base/DiagWriter_fwd.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/util/Trace.hpp>
#include <stk_util/util/WriterParser.hpp>
#endif

namespace stk { namespace diag { class Writer; } }
namespace stk { namespace mesh { class Part; } }


// Note, this classes/functions in this header are for internal use only.
// The API for tracing is defined in Trace.hpp

namespace stk {
namespace mesh {

namespace impl {
class Partition;
}

#ifdef STK_MESH_TRACE_ENABLED

// Must be called before theDiagWriter/meshlog
void initDiagWriter(std::ostream& stream);

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

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityKey& key);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const EntityProc& entity_proc);

stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Bucket& bucket);

namespace impl {
stk::diag::Writer& operator<<(stk::diag::Writer& writer, const Partition& partition);
}

#endif // STKMESH_TRACE_ENABLED

} // namespace mesh
} // namespace stk

#endif // stk_mesh_DiagWriter_h
