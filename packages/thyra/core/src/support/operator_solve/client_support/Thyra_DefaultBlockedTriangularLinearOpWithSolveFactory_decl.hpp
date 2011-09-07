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

#ifndef THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpSourceBase.hpp"


namespace Thyra {


/** \brief Implicit subclass that takes a blocked triangular LOWB object and
 * turns it into a LOWSB object.
 *
 * This class takes any upper or lower triangular
 * <tt>PhysicallyBlockedLinearOpBase</tt> object and compatible
 * <tt>LinearOpWithSolveFactoryBase</tt> object(s) and creates a LOWSB version
 * by creating LOWSB objects along the diagonal.
 *
 *
 * For example, consider the lower block triangular linear operator:

 \verbatim

       [ M(0,0)                   ]
   M = [ M(1,0)   M(1,1)          ]
       [ M(2,0)   M(2,1)   M(2,2) ]  

 \endverbatim
 
 * This class object will then create a new LOWSB object (of type
 * <tt>DefaultBlockedTriangularLinearOpWithSolve</tt>) that looks like:

 \verbatim

          [ invM(0,0)                       ]
   invM = [ M(1,0)     invM(1,1)            ]
          [ M(2,0)     M(2,1)     invM(2,2) ]  

 \endverbatim

 * where <tt>invM(k,k)</tt> are LOWSB objects created from the LOB objects
 * <tt>M(k,k)</tt> given a LOWSFB object.
 *
 * This class is not very compliciated, see the function
 * <tt>initializeOp()</tt> see what this class actually does!
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class DefaultBlockedTriangularLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar>
{
public:

  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{
  
  /** \brief Create given a single non-const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  DefaultBlockedTriangularLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
    );

  
  /** \brief Create given a single const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  DefaultBlockedTriangularLinearOpWithSolveFactory(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf
    );

  // 2007/10/02: rabartl: Add versions of constructor that accept an array of
  // LOWSFB objects when needed.  This will be needed for multi-physics
  // problems for instance!

  /** \brief . */
  RCP<LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF();

  /** \brief . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF() const;

  //@}

  /** \name Overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
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

  /** \brief . */
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

protected:

  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  /** \brief . */
  void informUpdatedVerbosityState() const;

  //@}

private:

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveFactoryBase<Scalar> > LOWSF_t;
  
  LOWSF_t lowsf_;

  // Not defined and not to be called
  DefaultBlockedTriangularLinearOpWithSolveFactory();

};


/** \brief Nonmember constructor.
 *
 * \releates DefaultBlockedTriangularLinearOpWithSolveFactory
 */
template<class Scalar>
RCP<DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar> >
defaultBlockedTriangularLinearOpWithSolveFactory(
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
  return Teuchos::rcp(
    new DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>(lowsf)
    );
}


/** \brief Nonmember constructor.
 *
 * \releates DefaultBlockedTriangularLinearOpWithSolveFactory
 */
template<class Scalar>
RCP<DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar> >
defaultBlockedTriangularLinearOpWithSolveFactory(
  const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
  return Teuchos::rcp(
    new DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>(lowsf)
    );
}


} // namespace Thyra


#endif // THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
