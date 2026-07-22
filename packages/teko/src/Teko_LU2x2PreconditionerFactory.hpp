// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_LU2x2PreconditionerFactory_hpp__
#define __Teko_LU2x2PreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_LU2x2Strategy.hpp"

namespace Teko {

/** \brief Construct a preconditioner using a LDU dcomposition of a block
 *  2x2 matrix.
 *
 * This produces a preconditioner using the block-LDU decomposition of
 * the matrix. The general assumption made is that the matrix is 2x2
 * and the block factorization can be constructed (i.e. assumptions about
 * the invertability of some blocks). The pattern used, and the one you
 * should follow if you want to use this software is
 *
 * \f$
 * A = \left[
 * \begin{array}{cc}
 * A_{00} & A_{01} \\
 * A_{10} & A_{11}
 * \end{array}
 * \right]
 * = \left[
 * \begin{array}{cc}
 * I & 0  \\
 * A_{10} A_{00}^{-1} & I
 * \end{array}
 * \right]
 * \left[
 * \begin{array}{cc}
 * A_{00} & 0  \\
 * 0 & -S
 * \end{array}
 * \right]
 * \left[
 * \begin{array}{cc}
 * I &  A_{00}^{-1} A_{01} \\
 * 0 & I
 * \end{array}
 * \right]
 * \f$
 *
 * where the Schur complement is \f$ S=-A_{11}+A_{10} A_{00}^{-1} A_{01} \f$ .
 *
 * To use an LDU approximation 2 evaluations of \f$ A_{00}^{-1} \f$ and a single
 * evalution of \f$ S^{-1} \f$ are needed. For increased flexibility both
 * evaluations of \f$A_{00}^{-1}\f$ can be specified independently.
 * For righthand side vector \f$[f, g]^T\f$ and solution vector \f$[u,v]^T\f$
 * the two inverses (\f$A\f$-hat and \f$A\f$-tilde) are needed to evaluate
 *
 * \f$\hat{A}_{00} u^* = f\f$,
 *
 * \f$\tilde{A}_{00} v = A_{01} v\f$
 *
 * where \f$u^*\f$ is an intermediate step.
 *
 * In order to facilate using this class in a nonlinear solve (or for a
 * time-dependent problem) the additional abstraction of a ``Strategy''
 * has been added. This strategy, abstractly represented as the LU2x2Strategy,
 * provides the \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ operators. Typical usage for this class
 * is to build a LU2x2Strategy and pass it into the primary constructor.
 * Additional constructors are provided both for convenience and to ease
 * adoption. Underneath the hood all these constructors do is invoke the
 * corresponding strategy object.
 *
 * For example, assume that you have the particularly nice case that
 * your approximations of \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ are independent of the source
 * operator. Then, one way to instantiate a LU2x2PreconditionerFactory
 * is

   <code>
      RCP<LinearOpBase<double> > invA00 = buildInvA00(...);\n
      RCP<LinearOpBase<double> > invS   = buildInvS(...);\n
      RCP<LU2x2PreconditionerFactory> precFact = rcp(new LU2x2PreconditionerFactory(invA00,invS));
   </code>

 * Now using the strategy constructor, an entirely equivalent factory
 * object can be constructed by

   <code>
      RCP<LinearOpBase<double> > invA00 = buildInvA00(...);\n
      RCP<LinearOpBase<double> > invS   = buildInvS(...);\n
      RCP<LU2x2Strateghy> precStrat = rcp(new StaticLU2x2Strategy(invA00,invS));\n
      RCP<LU2x2PreconditionerFactory> precFact = rcp(new LU2x2PreconditionerFactory(precStrat));
   </code>

 * Notice that the StaticLU2x2Strategy takes the same objects
 * as the original constructor, it acts as an intermediary to tell the
 * LU2x2PreconditionerFactory what those operators are.
 **/
class LU2x2PreconditionerFactory : public BlockPreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  /** @brief Build a simple static LU2x2 preconditioner */
  LU2x2PreconditionerFactory(LinearOp& invA00, LinearOp& invS);

  /** @brief Build a simple static LU2x2 preconditioner */
  LU2x2PreconditionerFactory(LinearOp& hatInvA00, LinearOp& tildeInvA00, LinearOp& invS);

  /** @brief Constructor that permits the most generality in building \f$A_{00}^{-1}\f$ and
   *        \f$S^{-1}\f$.
   *
   * Constructor that permits the most generality in building \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$.
   *
   * @param[in] strategy  Strategy object that takes a 2x2 block matrix and
   *                      and constructs the \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ objects.
   */
  LU2x2PreconditionerFactory(const Teuchos::RCP<LU2x2Strategy>& strategy);

  /** \brief Default constructor for use with AutoClone
   *
   * Default constructor for use with AutoClone
   */
  LU2x2PreconditionerFactory();

  //@}

  /** \brief Create the LU 2x2 preconditioner operator.
   *
   * This method breaks apart the BlockLinearOp and builds a block
   * LU preconditioner. This will require two applications of the inverse
   * of the (0,0) block and one application of the inverse Schur complement.
   */
  LinearOp buildPreconditionerOperator(BlockedLinearOp& blo, BlockPreconditionerState& state) const;

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

  /** \brief Determine the type of inverse operator to build.
   *
   * Determine the type of inverse operator to build. If true
   * use the full LDU decomposition. If false only the upper
   * triangular solve should be used. Motivation for doing this
   * can be found in Murphy, Golub and Wathen, SISC 2000.
   *
   * \returns A boolean indicating the type of inverse operator
   *          to use.
   *
   * \note Default behavior for this class is to return true.
   */
  virtual bool useFullLDU() const { return useFullLDU_; }

  /** \brief Set the type of inverse operation to use.
   *
   * Set the type of inverse operator to use. If true
   * use the full LDU decomposition. If false only the upper
   * triangular solve should be used. Motivation for doing this
   * can be found in Murphy, Golub and Wathen, SISC 2000.
   *
   * \param[in] value Boolean indicating type of inverse operator
   *                  to build.
   */
  virtual void setFullLDU(bool value) { useFullLDU_ = value; }

 protected:
  using Teko::BlockPreconditionerFactory::buildPreconditionerOperator;

  //! some members
  Teuchos::RCP<LU2x2Strategy> invOpsStrategy_;

  /** If true, use full LDU decomposition, otherwise
   * use the Golub & Wathen style upper block. This is
   * true by default.
   */
  bool useFullLDU_;

 public:
  /** \brief Builder function for creating strategies.
   *
   * Builder function for creating strategies.
   *
   * \param[in] name     String name of strategy to build
   * \param[in] settings Parameter list describing the parameters for the
   *                     strategy to build
   * \param[in] invLib   Inverse library for the strategy to use.
   *
   * \returns If the name is associated with a strategy
   *          a pointer is returned, otherwise Teuchos::null is returned.
   */
  static RCP<LU2x2Strategy> buildStrategy(const std::string& name,
                                          const Teuchos::ParameterList& settings,
                                          const RCP<const InverseLibrary>& invLib,
                                          const RCP<RequestHandler>& rh);

  /** \brief Add a strategy to the builder. This is done using the
   *        clone pattern.
   *
   * Add a strategy to the builder. This is done using the
   * clone pattern. If your class does not support the Cloneable interface then
   * you can use the AutoClone class to construct your object.
   *
   * \note If this method is called twice with the same string, the latter clone pointer
   *       will be used.
   *
   * \param[in] name String to associate with this object
   * \param[in] clone Pointer to Cloneable object
   */
  static void addStrategy(const std::string& name, const RCP<Cloneable>& clone);

 private:
  //! for creating the strategy objects
  static CloneFactory<LU2x2Strategy> strategyBuilder_;

  //! This is where the default objects are put into the strategyBuilder_
  static void initializeStrategyBuilder();
};

}  // end namespace Teko

#endif
