// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
#define THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {


/** \brief Base interface for a strategy object for randomizing a
 * multi-vector.
 *
 * This object is *not* stateless in its use!  Every time it generates a new
 * random multi-vector its behavior changes.
 *
 * A single <tt>MultiVectorRandomizerBase</tt> object may be compatible with
 * many different types of concrete vector space implementations or may
 * compatible with only a specific instantiation of a concrete vector space
 * subclass.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class MultiVectorRandomizerBase {
public:

  /** \brief . */
  virtual ~MultiVectorRandomizerBase() {}

  /** \brief Determines if <tt>*this</tt> is compatible with multi-vectors
   * from the <tt>VectorSpace</tt> <tt>space</tt>.
   */
  virtual bool isCompatible( const VectorSpaceBase<Scalar> &space ) const = 0;

  /** \brief Randomize a "compatible" multi-vector.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>mv!=NULL</tt>
   * <li><tt>this->isCompatible(*mv->range()) == true</tt>
   * </ul>
   */
  void randomize(const Ptr<MultiVectorBase<Scalar> > &mv)
    { randomizeImpl(mv); }

private:

  /** \brief . */
  virtual void randomizeImpl(const Ptr<MultiVectorBase<Scalar> > &mv) = 0;
  
};


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
