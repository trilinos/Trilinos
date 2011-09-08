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

#ifndef THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP
#define THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP

#include "Thyra_VectorSpaceFactoryBase.hpp"

namespace Teuchos { template<typename Ordinal> class Comm; }

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

  /** \brief Return the Spmd communicator. */
  Teuchos::RCP<const Teuchos::Comm<Ordinal> > getComm() const;

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
   return defaultSpmdVectorSpace(this->getComm(), dim, dim)</tt>
   \endcode
   
   * and therefore <tt>return->dim()==dim</tt> and this implementation fully
   * satisfies the specification of
   * <tt>VectorSpaceFactoryBase::createVecSpc()</tt>.
   */
   Teuchos::RCP<const VectorSpaceBase<Scalar> > createVecSpc(int dim) const;

  //@}

private:

  Teuchos::RCP<const Teuchos::Comm<Ordinal> >  comm_;

public:

  /** \brief Depreciated . */
  DefaultSpmdVectorSpaceFactory(
    const Teuchos::RCP<const Teuchos::Comm<Ordinal> > &comm = Teuchos::null
    );
    
}; // end class DefaultSpmdVectorSpaceFactory


/** \brief Construct with a <tt>Teuchos::Comm</tt> object.
 *
 * \param comm [in] The communicator.  This object must be maintained by the
 * client the entire time that <tt>this</tt> is in use.  It is allowed for
 * <tt>comm.get()==NULL</tt>.
 * 
 * Postconditions:<ul>
 * <li><tt>returnVal->getComm().get() == comm.get()</tt>
 * </ul>
 *
 * \relates DefaultSpmdVectorSpaceFactory
 */
template<class Scalar>
RCP<DefaultSpmdVectorSpaceFactory<Scalar> >
defaultSpmdVectorSpaceFactory(
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> > &comm = Teuchos::null
  )
{
  return Teuchos::rcp(new DefaultSpmdVectorSpaceFactory<Scalar>(comm));
}


// ///////////////////////////
// Inline members


template<class Scalar>
inline
Teuchos::RCP<const Teuchos::Comm<Ordinal> >
DefaultSpmdVectorSpaceFactory<Scalar>::getComm() const
{
  return comm_;
}

} // end namespace Thyra

#endif  // THYRA_SPMD_VECTOR_SPACE_FACTORY_STD_DECL_HPP
