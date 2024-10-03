// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_IterativePreconditionerFactory_hpp__
#define __Teko_IterativePreconditionerFactory_hpp__

// Teko includes
#include "Teko_PreconditionerFactory.hpp"

namespace Teko {

/** \brief A class which applies a preconditioner repeatedly.
  *        The inherit assumption is that the preconditioner corresponds
  *        to a residual correction.
  *
  * A class which applies a preconditioner repeatedly.
  * The inherit assumption is that the preconditioner corresponds
  * to a residual correction. For a linear operator \f$A\f$ preconditioned
  * by \f$P^{-1}\f$ the initial residual correction is
  *
    \f$
    e_0 = P^{-1} r_0.
    \f$
  *
  * Now because the preconditioner approximates \f$A^{-1}\f$ we have
  * \f$r_1 = r_0 - A e_0\f$. This leads to the following residual
  * correction scheme
  *
    \f$
    r_i = (I - A P^{-1})^i r_0.
    \f$
  *
  * An additional application of the preconditioner will give the
  * $i^{th}$ correction
  *
    \f$
    e_i = P^{-1} \sum_{j=0}^i r_j.
    \f$
  *
  * This factory takes a preconditioner (or inverse) factory and
  * constructs \f$P^{-1}\f$, and then applies it, as outlined above,
  * a user specified number of times.
  *
  * Similar to the other preconditioner in Teko this class can be
  * constructed through an XML file. For example:
  \code
    <ParameterList name="Iterative Solve">
      <Parameter name="Type" type="string" value="Iterative Preconditioner"/>
      <Parameter name="Preconditioner Type" type="string" value="ML"/>
      <Parameter name="Iteration Count" type="int" value="3"/>
    </ParameterList>
  \endcode
  */
class IterativePreconditionerFactory : public virtual Teko::PreconditionerFactory {
 public:
  //! Default constructor, for use with the AutoClone class.
  IterativePreconditionerFactory();

  /** Construct a preconditioner factory that applies a specified
   * preconditioner, a fixed number of times.
   *
   * \param[in] correctionNum The correction number to be returned
   *                          by the preconditioner operator. If the
   *                          preconditioner is only applied once than
   *                          <code>correctionNum=0</code> and \f$e_0\f$
   *                          is returned.
   * \param[in] precFactory Factory used to construct the preconditioner
   */
  IterativePreconditionerFactory(unsigned int correctionNum,
                                 const Teuchos::RCP<Teko::InverseFactory>& precFactory);

  /** Construct a preconditioner factory that applies a specified
   * preconditioner, a fixed number of times.
   *
   * \param[in] correctionNum The correction number to be returned
   *                          by the preconditioner operator. If the
   *                          preconditioner is only applied once than
   *                          <code>correctionNum=0</code> and \f$e_0\f$
   *                          is returned.
   * \param[in] precFactory Factory used to construct the preconditioner
   */
  IterativePreconditionerFactory(unsigned int correctionNum,
                                 const Teuchos::RCP<Teko::PreconditionerFactory>& precFactory);

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
  virtual LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const;

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
   */
  virtual void initializeFromParameterList(const Teuchos::ParameterList& settings);

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
   */
  virtual bool updateRequestedParameters(const Teuchos::ParameterList& pl);

  //@}

 protected:
  unsigned int correctionNum_;
  Teuchos::RCP<Teko::InverseFactory> precFactory_;
};

}  // end namespace Teko

#endif
