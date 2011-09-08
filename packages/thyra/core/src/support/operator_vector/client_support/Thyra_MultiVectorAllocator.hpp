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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_ALLOCATOR_HPP
#define THYRA_MULTI_VECTOR_ALLOCATOR_HPP

#include "Thyra_VectorSpaceBase.hpp"
#include "Teuchos_TestForException.hpp"

namespace Thyra {

/** \brief Allocator class to be used with <tt>Teuchos::AbstractFactoryStd</tt> to create
 * <tt>MultiVectorBase</tt> objects of a given size.
 */
template<class Scalar>
class MultiVectorAllocator {
public:
  /** \brief . */
  MultiVectorAllocator() : numMembers_(0) {}
  /** \brief . */
  typedef Teuchos::RCP<MultiVectorBase<Scalar> >  ptr_t;         // required!
  /** \brief . */
  MultiVectorAllocator( const Teuchos::RCP<const VectorSpaceBase<Scalar> > &vs, int numMembers )
    : vs_(vs), numMembers_(numMembers)
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( vs.get()==NULL, std::logic_error, "Error!" );
#endif			
    }
  /** \brief . */
  const ptr_t allocate() const { return vs_->createMembers(numMembers_); }  // required!
private:
  Teuchos::RCP<const VectorSpaceBase<Scalar> >      vs_;
  int                                                       numMembers_;
};

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_ALLOCATOR_HPP
