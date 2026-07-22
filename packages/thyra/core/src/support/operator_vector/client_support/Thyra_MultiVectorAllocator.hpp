// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_ALLOCATOR_HPP
#define THYRA_MULTI_VECTOR_ALLOCATOR_HPP

#include "Thyra_VectorSpaceBase.hpp"
#include "Teuchos_Assert.hpp"

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
      TEUCHOS_TEST_FOR_EXCEPTION( vs.get()==NULL, std::logic_error, "Error!" );
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
