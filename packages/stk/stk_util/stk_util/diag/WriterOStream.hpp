#ifndef STK_UTIL_DIAG_WRITEROSTREAM_HPP
#define STK_UTIL_DIAG_WRITEROSTREAM_HPP

#include <stk_util/diag/Writer.hpp>

namespace stk {
namespace diag {

///
/// @addtogroup diag_writer_detail
/// @{
///

/**
 * @brief Function <code>operator<<</code> is the catch all std::ostream output put-to operator to
 * Writer put-to operator.  When using this, if you attempt to put and object that has no put-to
 * operator to std::ostream, expect to get a list of all opt-to operator defined for the
 * std::ostream.
 *
 * @param dout          a <b>Writer</b> reference to the writer to put to.
 *
 * @param t             a <b>T</b> const reference to the object to put.
 * 
 */
template <class T>
Writer &operator<<(Writer &dout, const T &t) {
  if (dout.shouldPrint())
    dout.getStream() << t;
  
  return dout;
}

///
/// @}
///

} // namespace diag
} // namespace stk

#endif // STK_UTIL_DIAG_WRITEROSTREAM_HPP
