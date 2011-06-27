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

#ifndef TEUCHOS_ABSTRACT_FACTORY_HPP
#define TEUCHOS_ABSTRACT_FACTORY_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos {

/** \brief Simple, universal "Abstract Factory" interface for the
 * dynamic creation of objects.
 *
 * While <tt>RCP</tt> provides for specialized deallocation
 * policies it does not abstract, in any way, how an object is first
 * allocated.  The most general way to abstract how an object is
 * allocated is to use an "Abstract Factory".  This base class defines
 * the most basic "Abstract Factory" interface and defines only one
 * virtual function, <tt>create()</tt> that returns a
 * <tt>RCP</tt>-wrapped object.
 */
template<class T>
class AbstractFactory {
public:

#ifndef DOXYGEN_COMPILE
	/** \brief . */
	typedef Teuchos::RCP<T>   obj_ptr_t;
#endif

	/** \brief . */
	virtual ~AbstractFactory() {}

	/** \brief Create an object of type T returned as a smart reference
	 * counting pointer object.
	 */
	virtual obj_ptr_t create() const = 0;

}; // class AbstractFactory

} // end Teuchos

#endif // TEUCHOS_ABSTRACT_FACTORY_HPP
