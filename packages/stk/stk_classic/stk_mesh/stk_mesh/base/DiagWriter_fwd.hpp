#ifndef stk_mesh_DiagWriter_fwd_h
#define stk_mesh_DiagWriter_fwd_h

#include <stk_util/diag/Writer_fwd.hpp>

namespace stk {
namespace mesh {

enum LogMask {
  LOG_ALWAYS            = stk::LOG_ALWAYS,
  LOG_TRACE             = stk::LOG_TRACE,
  LOG_TRACE_STATS       = stk::LOG_TRACE_STATS,
  LOG_TRACE_SUB_CALLS   = stk::LOG_TRACE_SUB_CALLS,
  LOG_MEMBERS           = stk::LOG_MEMBERS,

  LOG_STREAM_COMMON     = LOG_TRACE | LOG_TRACE_STATS,
  LOG_BUCKET            = 0x00000100,
  LOG_ENTITY            = 0x00000200,
  LOG_PART              = 0x00000300,
  LOG_FIELD             = 0x00000400,

  LOG_MOD_END           = 0x00000800
};

} // namespace mesh
} // namespace stk

#endif // stk_mesh_DiagWriter_fwd_h
