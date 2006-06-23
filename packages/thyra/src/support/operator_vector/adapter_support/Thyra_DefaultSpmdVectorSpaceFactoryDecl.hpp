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

#ifndef THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP
#define THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP

#include "Thyra_VectorSpaceFactoryBaseDecl.hpp"

namespace Thyra {

/** \brief Concrete implementation of a vector-space factory for a
 * locally-replicated distributed <tt>DefaultSpmdVectorSpace</tt> objects.
 *
 * This will create either serial independent (<tt>comm.get()==NULL</tt>) or
 * locally replicated (<tt>comm.get()!=NULL</tt>) vector space objects (see
 * <tt>createVecSpc()</tt>).  The primary motivation for this subclass is to
 * create locally replicated vector spaces for the domain space of
 * <tt>DefaultSpmdMultiVector</tt>.  In addition, an object of this type is
 * also returned from <tt>SpmdVectorSpaceDefaultBase::smallVecSpcFtcy()</tt>.
 *
 * Note that the default constructor is not allowed to avoid mistakes in using
 * this class.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class DefaultSpmdVectorSpaceFactory : public VectorSpaceFactoryBase<Scalar> {
public:

  /** \brief Construct with an <tt>Spmd_Comm</tt> object.
   *
   * \param  comm
   *           [in] The  communicator.  This object must be maintained
   *           by the client the entire time that <tt>this</tt> is in use.
   * 
   * Postconditions:<ul>
   * <li><tt>this->getComm().get() == comm.get()</tt>
   * </ul>
   */
  DefaultSpmdVectorSpaceFactory(
    const Teuchos::RefCountPtr<const Teuchos::Comm<Index> > &comm
    );

  /** \brief Return the Spmd communicator. */
  Teuchos::RefCountPtr<const Teuchos::Comm<Index> > getComm() const;

  /** @name Overridden from VectorSpaceFactoryBase */
  //@{

  /** \brief Create a new locally-replicated <tt>DefaultSpmdVectorSpace</tt>
   * object given its dimension!
   *
   * \param  dim
   *           [in] The dimension of the (locally replicated) vector space to create.
   *
   * This function returns:
   
   \code
   return Teuchos::rcp(new DefaultSpmdVectorSpace(this->getComm(),dim,dim))</tt>
   \endcode
   
   * and therefore <tt>return->dim()==dim</tt> and this implementation fully
   * satisfies the specification of
   * <tt>VectorSpaceFactoryBase::createVecSpc()</tt>.
   */
   Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > createVecSpc(int dim) const;

  //@}

private:

  Teuchos::RefCountPtr<const Teuchos::Comm<Index> >  comm_;

  DefaultSpmdVectorSpaceFactory(); // Not defined and not to be called!
    
}; // end class DefaultSpmdVectorSpaceFactory

// ///////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<const Teuchos::Comm<Index> >
DefaultSpmdVectorSpaceFactory<Scalar>::getComm() const
{
  return comm_;
}

} // end namespace Thyra

#endif  // THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP
