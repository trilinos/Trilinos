// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
