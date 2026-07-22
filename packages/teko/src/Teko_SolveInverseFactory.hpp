// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_SolveInverseFactory_hpp__
#define __Teko_SolveInverseFactory_hpp__

#include "Teko_InverseFactory.hpp"

namespace Teko {

class SolveInverseFactory : public InverseFactory {
 public:
  //! \name Constructors
  //@{

  /** \brief Constructor that takes a Thyra solve factory and
   *        makes it look like an InverseFactory
   *
   * Constructor that takes a Thyra solve factory and
   * makes it look like an InverseFactory.
   *
   * \param[in] lowsFactory Thyra LineaerOpWithSolveFactoryBase used for building
   *                        the inverse.
   */
  SolveInverseFactory(
      const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >& lowsFactory);

  //! Copy constructor
  SolveInverseFactory(const SolveInverseFactory& siFactory);
  //@}

  virtual ~SolveInverseFactory() {}

  /** \brief Build an inverse operator
   *
   * Build the inverse operator using this factory.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp) const;

  /** \brief Build a preconditioned inverse operator
   *
   * Build the inverse operator using this factory and a user specified
   * preconditioning operator. The default behavior is to call buildInverse
   * ignoring the preconditioner.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   * \param[in] precOp Preconditioning operator
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp, const LinearOp& precOp) const;

  /** \brief Pass in an already constructed inverse operator. Update
   *        the inverse operator based on the new source operator.
   *
   * Pass in an already constructed inverse operator. Update
   * the inverse operator based on the new source operator.
   *
   * \param[in]     source Source operator to be inverted.
   * \param[in,out] dest   Pre constructed inverse operator to be
   *                        rebuilt using the <code>source</code>
   *                        object.
   */
  virtual void rebuildInverse(const LinearOp& source, InverseLinearOp& dest) const;

  /** \brief Pass in an already constructed inverse operator. Update
   *        the inverse operator based on the new source operator.
   *
   * Pass in an already constructed inverse operator. Update
   * the inverse operator based on the new source operator.
   *
   * \param[in]     source Source operator to be inverted.
   * \param[in]     precOp Preconditioning operator
   * \param[in,out] dest   Pre constructed inverse operator to be
   *                        rebuilt using the <code>source</code>
   *                        object.
   */
  virtual void rebuildInverse(const LinearOp& source, const LinearOp& precOp,
                              InverseLinearOp& dest) const;

  /** \brief A function that permits inspection of the parameters used to create
   *        this object.
   *
   * A function that permits inspection of the parameters used to create this
   * object. Useful for determining defaults and settings used.
   *
   * \returns A list used to parameterize this object.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  //! Accessor primarily for testing purposes
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > getLowsFactory() const {
    return lowsFactory_;
  }

  /** Return a string that describes this factory */
  virtual std::string toString() const { return lowsFactory_->description(); }

 protected:
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory_;

 private:
  // hide me!
  SolveInverseFactory();
};

}  // end namespace Teko

#endif
