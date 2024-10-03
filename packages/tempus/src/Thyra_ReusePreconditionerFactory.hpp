//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_ReusePreconditionerFactory_hpp
#define Thyra_ReusePreconditionerFactory_hpp

#include "Thyra_PreconditionerFactoryBase.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerFactoryBase</tt> subclass that
 * just returns an already created/initialized preconditioner object.
 */
template <class Scalar>
class ReusePreconditionerFactory
  : virtual public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  ReusePreconditionerFactory() {}

  void initialize(const RCP<PreconditionerBase<Scalar> > &prec)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
#endif
    prec_ = prec;
  }

  RCP<PreconditionerBase<Scalar> > getNonconstPreconditioner() { return prec_; }

  RCP<const PreconditionerBase<Scalar> > getPreconditioner() const
  {
    return prec_;
  }

  void uninitialize() { prec_ = Teuchos::null; }

  /** \name Overridden from Teuchos::Describable. */
  //@{

  std::string description() const
  {
    std::ostringstream oss;
    oss << this->Teuchos::Describable::description() << "{"
        << "prec=";
    if (!is_null(prec_))
      oss << prec_->description();
    else
      oss << "NULL";
    oss << "}";
    return oss.str();
  }

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions)
   */
  //@{

  void setParameterList(RCP<ParameterList> const & /* paramList */) {}

  RCP<ParameterList> getNonconstParameterList() { return Teuchos::null; }

  RCP<ParameterList> unsetParameterList() { return Teuchos::null; }

  RCP<const ParameterList> getParameterList() const { return Teuchos::null; }

  RCP<const ParameterList> getValidParameters() const
  {
    return rcp(new ParameterList);
  }

  //@}

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  bool isCompatible(const LinearOpSourceBase<Scalar> & /* fwdOpSrc */) const
  {
    return false;
  }

  RCP<PreconditionerBase<Scalar> > createPrec() const { return prec_; }

  void initializePrec(
      const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
      PreconditionerBase<Scalar> * /* precOp */,
      const ESupportSolveUse /* supportSolveUse */ =
          SUPPORT_SOLVE_UNSPECIFIED) const
  {
  }

  void uninitializePrec(PreconditionerBase<Scalar> * /* precOp */,
                        RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
                        ESupportSolveUse *supportSolveUse                = NULL) const
  {
  }

  //@}

 private:
  // //////////////////////////////
  // Private data members

  RCP<PreconditionerBase<Scalar> > prec_;
};

/** \brief Nonmember constructor function.
 *
 * \relates ReusePreconditionerFactory
 */
template <class Scalar>
RCP<ReusePreconditionerFactory<Scalar> > reusePreconditionerFactory()
{
  return Teuchos::rcp(new ReusePreconditionerFactory<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates ReusePreconditionerFactory
 */
template <class Scalar>
RCP<ReusePreconditionerFactory<Scalar> > reusePreconditionerFactory(
    const RCP<PreconditionerBase<Scalar> > &prec)
{
  RCP<ReusePreconditionerFactory<Scalar> > fac =
      Teuchos::rcp(new ReusePreconditionerFactory<Scalar>());
  fac->initialize(prec);
  return fac;
}

}  // end namespace Thyra

#endif
