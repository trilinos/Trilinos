// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H
#define TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/** \brief <tt>basic_ostream<></tt> subclass that does nothing but discard output.
 *
 * \ingroup teuchos_outputting_grp
 *
 * Use the class anytime you must pass an <tt>basic_ostream<></tt> object
 * but don't want the output for any reason.
 *
 * This subclass just sets the stream buffer to NULL and that is all you need to do!
 */
template<typename _CharT, typename _Traits>
class basic_oblackholestream
	: virtual public std::basic_ostream<_CharT, _Traits>
{
public:
  /** \brief . */
	explicit basic_oblackholestream() : std::basic_ostream<_CharT, _Traits>(NULL) {}
}; // end class basic_oblackholestream

} // end namespace Teuchos

#endif // TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H
