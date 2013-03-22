// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
