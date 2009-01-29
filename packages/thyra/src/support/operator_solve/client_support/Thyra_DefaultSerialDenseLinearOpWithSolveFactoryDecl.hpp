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

#ifndef THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP


#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Thyra {


/** \brief Concreate LinearOpWithSolveFactoryBase subclass that creates
 * DefaultSerialDenseLinearOpWithSolve objects that use LAPACK.
 *
 * This class will work with any serial MultiVectorBase object for which there
 * are BLAS and LAPACK wrappers in Teuchos for.
 */
template<class Scalar>
class DefaultSerialDenseLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar>,
    virtual protected Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from LinearOpWithSolveFactoyBase */
  //@{
  
  /** \brief returns false. */
  virtual bool acceptsPreconditionerFactory() const;

  /** \brief Throws exception. */
  virtual void setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
    const std::string &precFactoryName
    );

  /** \brief Returns null . */
  virtual RCP<PreconditionerFactoryBase<Scalar> >
  getPreconditionerFactory() const;

  /** \brief Throws exception. */
  virtual void unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
    std::string *precFactoryName
    );

  /** \brief . */
  virtual bool isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc
    ) const;

  /** \brief Returns a DefaultSerialDenseLinearOpWithSolve object . */
  virtual RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  /** \brief . */
  virtual void initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  virtual void initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op
    ) const;

  /** \brief . */
  virtual void uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;
 
  /** \brief . */
  virtual bool supportsPreconditionerInputType(
    const EPreconditionerInputType precOpType
    ) const;

  /** \brief . */
  virtual void initializePreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  virtual void initializeApproxPreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  //@}

};


/** \brief Nonmember constructor.
 *
 * \releates DefaultSerialDenseLinearOpWithSolveFactory
 */
template<class Scalar>
RCP<DefaultSerialDenseLinearOpWithSolveFactory<Scalar> >
defaultSerialDenseLinearOpWithSolveFactory()
{
  return Teuchos::rcp(new DefaultSerialDenseLinearOpWithSolveFactory<Scalar>);
}


} // namespace Thyra


#endif // THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
