/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_MANIP_HPP
#define STK_UTIL_DIAG_MANIP_HPP

#include <iomanip>

#include <stk_util/diag/Writer_fwd.hpp>
#include <stk_util/diag/Writer.hpp>

namespace stk {
namespace diag {

///
/// @addtogroup diag_writer_detail
/// @{
///

/**
 * @brief Class <b>_setw</b> is the width manipulator.
 *
 */
struct _setw
{
  _setw(int width)
    : m_width(width)
  {}

  int           m_width;
};

/**
 * @brief Function <code>setw</code> sets the width for the next field as a manipulator.
 *
 * @param width         a <b>int</b> value for the width of the next field.
 * 
 */
inline _setw setw(int width) {
  return _setw(width);
}

Writer &operator<<(Writer &dout, _setw set_width);


/**
 * @brief Class <b>_setprecision</b> is the precision manipulator.
 *
 */
struct _setprecision
{
  _setprecision(int precision)
    : m_precision(precision)
  {}

  int           m_precision;
};

/**
 * @brief Function <code>setprecision</code> sets the numeric precision as a manipulator.
 *
 * @param precision     a <b>int</b> value of the precision.
 * 
 */
inline _setprecision setprecision(int precision) {
  return _setprecision(precision);
}

Writer &operator<<(Writer &dout, _setprecision set_precision);


/**
 * @brief Class <b>_setfill</b> is the fill character manipulator.
 *
 */
struct _setfill
{
  _setfill(char fill)
    : m_fill(fill)
  {}

  char          m_fill;
};

/**
 * @brief Function <code>setfill</code> sets the fill character as a manipulator.
 *
 * @param fill          a <b>char</b> value of the fill character.
 * 
 */
inline _setfill setfill(char fill) {
  return _setfill(fill);
}

Writer &operator<<(Writer &dout, _setfill set_fill);


/**
 * @brief Class <b>_setiosflags</b> is the flags manipulator.
 *
 */
struct _setiosflags
{
  _setiosflags(std::ios_base::fmtflags flags)
    : m_flags(flags)
  {}

  std::ios_base::fmtflags  m_flags;
};

/**
 * @brief Function <code>setiosflags</code> sets the ios flags as a manipulator.
 *
 * @param flags         a <b>std::ios_base::fmtflags</b> value of the flags.
 * 
 */
inline _setiosflags setiosflags(std::ios_base::fmtflags flags) {
  return _setiosflags(flags);
}

Writer &operator<<(Writer &dout, _setiosflags set_flags);


/**
 * @brief Class <b>_resetiosflags</b> is the reset ios flags reset manipulator.
 *
 */
struct _resetiosflags
{
  _resetiosflags(std::ios_base::fmtflags flags)
    : m_flags(flags)
  {}

  std::ios_base::fmtflags  m_flags;
};

/**
 * @brief Function <code>resetiosflags</code> clears the ios flags as a manipulator.
 *
 * @param flags         a <b>std::ios_base::fmtflags</b> value of the flags.
 * 
 */
inline _resetiosflags resetiosflags(std::ios_base::fmtflags flags) {
  return _resetiosflags(flags);
}

Writer &operator<<(Writer &dout, _resetiosflags reset_flags);

Writer &fixed(Writer &dout);
Writer &scientific(Writer &dout);
Writer &dec(Writer &dout);
Writer &hex(Writer &dout);
Writer &oct(Writer &dout);

///
/// @}
///

} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

using stk::diag::setw;
using stk::diag::setfill;
using stk::diag::setprecision;

} // namespace Diag
} // namespace sierra

#endif // STK_UTIL_DIAG_MANIP_HPP
