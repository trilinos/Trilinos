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

#ifndef TEUCHOS_REDUCTION_OP_HPP
#define TEUCHOS_REDUCTION_OP_HPP

#include "Teuchos_Describable.hpp"


namespace Teuchos {


/** \brief Base interface class for user-defined reduction operations for
 * objects that use value semantics.
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal, typename T>
class ValueTypeReductionOp : public Describable {
public:
  /** \brief . */
  virtual void reduce(
    const Ordinal count,
    const T inBuffer[],
    T inoutBuffer[]
    ) const = 0;
};


/** \brief Base interface class for user-defined reduction operations for
 * objects that use reference semantics.
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal, typename T>
class ReferenceTypeReductionOp : public Describable {
public:
  /** \brief . */
  virtual void reduce(
    const Ordinal count,
    const T*const inBuffer[],
    T*const inoutBuffer[]
    ) const = 0;
};


} // namespace Teuchos


#endif // TEUCHOS_REDUCTION_OP_HPP
