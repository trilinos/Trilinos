// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_PreconditionerFactory_hpp__
#define __Teko_PreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_CloneFactory.hpp"
#include "Teko_PreconditionerState.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestHandlerContainer.hpp"

namespace Teko {

using Thyra::DefaultPreconditioner;
using Thyra::LinearOpBase;

/** \brief Abstract class which block preconditioner factories in Teko
 *        should be based on.
 *
 * Abstract class which block preconditioner factories in Teko should
 * be based on. All that is needed is the implementation of
 * "buildPreconditionerOperator". This class also uses the
 * <code>RequestHandler</code> concrete interface. This is useful for
 * extracting information from the user in an unobtrusive and modular
 * way.
 */
class PreconditionerFactory : public virtual Thyra::PreconditionerFactoryBase<double>,
                              public RequestHandlerContainer {
 public:
  /** \brief Function that is called to build the preconditioner
   *        for the linear operator that is passed in.
   *
   * This function builds a preconditioner based on the passed
   * in LinearOp.
   *
   * \param[in] lo    Source linear operator that is to be preconditioned.
   * \param[in] state An object associated with this operator to store
   *                  the preconditioner state.
   *
   * \returns The preconditioner as a linear operator (i.e. to perform
   *           a matrix-vector operation simply call "apply").
   */
  virtual LinearOp buildPreconditionerOperator(LinearOp &lo, PreconditionerState &state) const = 0;

  /** \brief Function that permits the construction of an arbitrary
   *        PreconditionerState object.
   *
   * Function that permits the construction of an arbitrary
   * PreconditionerState object. If the basic state object,
   * which takes a parameter list, is sufficient the default behavior
   * does precisely what is needed. Otherwise, an author of a
   * PreconditionerFactory would need to reimplement this method to
   * return a new state object.
   *
   * \returns A state object associated with this factory.
   */
  virtual Teuchos::RCP<PreconditionerState> buildPreconditionerState() const {
    return Teuchos::rcp(new PreconditionerState());
  }

  //! @name Methods for construction from a parameter list entry
  //@{

  /** \brief This function builds the internals of the preconditioner factory
   *        from a parameter list.
   *
   * This function builds the internals of the preconditioner factory
   * from a parameter list. Furthermore, it allows a preconditioner factory
   * developer to easily add a factory to the build system. This function
   * is required for building a preconditioner from a parameter list.
   *
   * \param[in] settings Parmaeter list to use as the internal settings
   *
   * \note The default implementation does nothing.
   */
  virtual void initializeFromParameterList(const Teuchos::ParameterList & /* settings */) {}

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
  virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const {
    return Teuchos::null;
  }

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
  virtual bool updateRequestedParameters(const Teuchos::ParameterList & /* pl */) { return true; }

  //@}

  //! Set the inverse library used by this preconditioner factory
  void setInverseLibrary(const Teuchos::RCP<const InverseLibrary> &il);

  //! Get the inverse library used by this preconditioner factory
  Teuchos::RCP<const InverseLibrary> getInverseLibrary() const;

  //! @name Methods inherited from Thyra::PreconditionerFactoryBase
  //@{

  //! is this operator compatiable with the preconditioner factory?
  bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const;

  //! create an instance of the preconditioner
  Teuchos::RCP<Thyra::PreconditionerBase<double> > createPrec() const;

  /** \brief initialize a newly created preconditioner object
   *
   * Initialize a newly created preconditioner object. For use with
   * nonlinear solvers.
   *
   * \param[in] fwdOpSrc Forward operator to be preconditioned
   * \param[in] solnVec Vector associated with this linear operator.
   * \param[in,out] precOp Return location for the preconditioner
   * \param[in] supportSolveUse Thyra information (?)
   */
  void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
                      const Teuchos::RCP<const Thyra::MultiVectorBase<double> > &solnVec,
                      Thyra::PreconditionerBase<double> *precOp,
                      const Thyra::ESupportSolveUse supportSolveUse) const;

  //! initialize a newly created preconditioner object
  void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
                      Thyra::PreconditionerBase<double> *precOp,
                      const Thyra::ESupportSolveUse supportSolveUse) const;

  //! wipe clean a already initialized preconditioner object
  void uninitializePrec(Thyra::PreconditionerBase<double> *prec,
                        Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > *fwdOpSrc,
                        Thyra::ESupportSolveUse *supportSolveUse) const;
  //@}

  //! @name Methods inherited from Teuchos::ParameterListAcceptor
  //@{

  //! Set parameters from a parameter list and return with default values.
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> &paramList);

  //! Get the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  //! Unset the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  //@}

  //! Set the request handler with pointers to the appropriate callbacks
  void setRequestHandler(const Teuchos::RCP<RequestHandler> &rh) { callbackHandler_ = rh; }

  //! Get the request handler with pointers to the appropriate callbacks
  Teuchos::RCP<RequestHandler> getRequestHandler() const { return callbackHandler_; }

 protected:
  //! for ParameterListAcceptor
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  //! For handling requests and send requests back to the user
  Teuchos::RCP<RequestHandler> callbackHandler_;

 private:
  //! Inverse library to be used by this factory
  Teuchos::RCP<const InverseLibrary> inverseLibrary_;

  //! If supported, set the request handler in this operator
  static void setOpRequestHandler(const RequestHandlerContainer &rhc, const LinearOp &op);

 public:
  /** \brief Builder function for creating preconditioner factories (yes
   *        this is a factory factory).
   *
   * Builder function for creating preconditioner factories (yes
   * this is a factory factory).
   *
   * \param[in] name     String name of factory to build
   * \param[in] settings Parameter list describing the parameters for the
   *                     factory to build
   * \param[in] invLib   Inverse library for the factory to use.
   *
   * \returns If the name is associated with a preconditioner
   *          a pointer is returned, otherwise Teuchos::null is returned.
   */
  static Teuchos::RCP<PreconditionerFactory> buildPreconditionerFactory(
      const std::string &name, const Teuchos::ParameterList &settings,
      const Teuchos::RCP<const InverseLibrary> &invLib = Teuchos::null);

  /** \brief Add a preconditioner factory to the builder. This is done using the
   *        clone pattern.
   *
   * Add a preconditioner factory to the builder. This is done using the
   * clone pattern. If your class does not support the Cloneable interface then
   * you can use the AutoClone class to construct your object.
   *
   * \note If this method is called twice with the same string, the latter clone pointer
   *       will be used.
   *
   * \param[in] name String to associate with this object
   * \param[in] clone Pointer to Cloneable object
   */
  static void addPreconditionerFactory(const std::string &name,
                                       const Teuchos::RCP<Cloneable> &clone);

  /** \brief Get the names of the block preconditioner factories
   */
  static void getPreconditionerFactoryNames(std::vector<std::string> &names);

 private:
  //! for creating the preconditioner factories objects
  static CloneFactory<PreconditionerFactory> precFactoryBuilder_;

  //! This is where the default objects are put into the precFactoryBuilder_
  static void initializePrecFactoryBuilder();
};

}  // end namespace Teko

#endif
