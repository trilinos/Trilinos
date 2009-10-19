#ifndef STK_UTIL_DIAG_WRITER_FWD_HPP
#define STK_UTIL_DIAG_WRITER_FWD_HPP

namespace stk {

///
/// @addtogroup diag_writer_detail
/// @{
///

/**
 * @brief Enumeration <code>LogMask</code> enumerates global bit assignments.
 *
 */
enum LogMask {
  LOG_ALWAYS		= 0x00000000,
  LOG_TRACE		= 0x00000001,
  LOG_TRACE_STATS	= 0x00000002,
  LOG_TRACE_SUB_CALLS	= 0x00000004,
  LOG_MEMBERS		= 0x00000008,

  LOG_STREAM_COMMON	= LOG_TRACE | LOG_TRACE_STATS
};

///
/// @}
///

namespace diag {

class Writer;

typedef unsigned long PrintMask;

} // namespace diag
} // namespace stk

#endif // STK_UTIL_DIAG_WRITER_FWD_HPP
