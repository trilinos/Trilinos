// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_PreconditionerInverseFactory_hpp__
#define __Teko_PreconditionerInverseFactory_hpp__

#include "Teko_InverseFactory.hpp"

namespace Teko {

class PreconditionerInverseFactory : public InverseFactory {
 public:
  using InverseFactory::buildInverse;
  using InverseFactory::rebuildInverse;

  //! \name Constructors
  //@{

  /** \brief Constructor that takes a Thyra solve factory and
   *        makes it look like an InverseFactory
   *
   * Constructor that takes a Thyra solve factory and
   * makes it look like an InverseFactory.
   *
   * \param[in] precFactory Thyra PreconditionerFactoryBase used for building
   *                        the inverse.
   */
  PreconditionerInverseFactory(
      const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> >& precFactory,
      const Teuchos::RCP<Teko::RequestHandler>& rh);

  /** \brief Constructor that takes a Thyra solve factory and
   *        makes it look like an InverseFactory. This constructor
   *        also permits the passing of an "Extra Parameters" parameter
   *        list.
   *
   * Constructor that takes a Thyra solve factory and
   * makes it look like an InverseFactory.  This constructor
   * also permits the passing of an "Extra Parameters" parameter
   * list to be used and updated through the "RequestedParameters" function.
   *
   * \param[in] precFactory Thyra PreconditionerFactoryBase used for building
   *                        the inverse.
   * \param[in] xtraParam Parameter list containing extra parameters.
   */
  PreconditionerInverseFactory(
      const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> >& precFactory,
      const Teuchos::RCP<const Teuchos::ParameterList>& xtraParam,
      const Teuchos::RCP<Teko::RequestHandler>& rh);

  //! Copy constructor
  PreconditionerInverseFactory(const PreconditionerInverseFactory& pFactory);
  //@}

  virtual ~PreconditionerInverseFactory() {}

  /** \brief Build an inverse operator
   *
   * Build the inverse operator using this factory. This returns
   * a linear operator that wraps a Thyra::PreconditionerBase object.
   * This PreconditionerBase object will be utilized when
   * <code>rebuildInverse</code> is called.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp) const;

  /** \brief Build an inverse operator and make sure it aware of some parents state
   *        This functionality is only useful for Teko::PreconditionerFactory inverses.
   *
   * Build an inverse operator and make sure it aware of some parents state
   * This functionality is only useful for Teko::PreconditionerFactory inverses.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   * \param[in] parentState Current state object to be used. Only useful for preconditioners.
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp,
                                       const PreconditionerState& parentState) const;

  /** \brief Pass in an already constructed inverse operator. Update
   *        the inverse operator based on the new source operator.
   *
   * Pass in an already constructed inverse operator. Update
   * the inverse operator based on the new source operator. This
   * method assumes the <code>dest</code> object also contains
   * the associated PreconditionerBase object as "prec" as extra
   * data in the RCP.
   *
   * \param[in]     source Source operator to be inverted.
   * \param[in,out] dest   Pre constructed inverse operator to be
   *                        rebuilt using the <code>source</code>
   *                        object.
   */
  virtual void rebuildInverse(const LinearOp& source, InverseLinearOp& dest) const;

  /** \brief A function that permits inspection of the parameters used to create
   *        this object.
   *
   * A function that permits inspection of the parameters used to create this
   * object. Useful for determining defaults and settings used.
   *
   * \returns A list used to parameterize this object.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  /** \brief Request the additional parameters this preconditioner factory
   *        needs.
   *
   * Request the additonal parameters needed by this preconditioner factory.
   * The parameter list will have a set of fields that can be filled with
   * the requested values. These fields include all requirements, even those
   * of the sub-solvers if there are any.  Once correctly filled the object
   * can be updated by calling the updateRequestedParameters with the filled
   * parameter list.
   *
   * \returns A parameter list with the requested parameters.
   *
   * \note The default implementation returns Teuchos::null.
   */
  virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

  /** \brief Update this object with the fields from a parameter list.
   *
   * Update the requested fields using a parameter list. This method is
   * expected to pair with the getRequestedParameters method (i.e. the fields
   * requested are going to be update using this method).
   *
   * \param[in] pl Parameter list containing the requested parameters.
   *
   * \returns If the method succeeded (found all its required parameters) this
   *          method returns true, otherwise it returns false.
   *
   * \note The default implementation returns true (it does nothing!).
   */
  virtual bool updateRequestedParameters(const Teuchos::ParameterList& pl);

  /** Return a string that describes this factory */
  virtual std::string toString() const { return precFactory_->description(); }

  /** Get the preconditioner factroy */
  Teuchos::RCP<const Thyra::PreconditionerFactoryBase<double> > getPrecFactory() const {
    return precFactory_;
  }

  /** Get the preconditioner factroy */
  Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > getPrecFactory() { return precFactory_; }

  /** This process the extra parameters passed in through the constructor.
   * Including the preRequest call and request calls to the request handler.
   * Another option would be to move the request call to the inverse construction
   * function.  However that breaks the "const" nature of the function and requires
   * the precFactory_ member to be mutable.  This is OK but not ideal. Note the user
   * should not call this directly because its called from the InverseLibrary::getInverseFactory
   * function.
   */
  void setupParameterListFromRequestHandler();

 protected:
  Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > precFactory_;
  Teuchos::RCP<Teuchos::ParameterList> extraParams_;

 private:
  // hide me!
  PreconditionerInverseFactory();
};

}  // namespace Teko

#endif
