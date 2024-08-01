// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_LSCStrategy_hpp__
#define __Teko_LSCStrategy_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"

namespace Teko {
namespace NS {

class LSCPrecondState;  // forward declaration

/** \brief Strategy for driving LSCPreconditionerFactory.
 *
 * Strategy for driving the LSCPreconditionerFactory. This
 * class provides all the pieces required by the LSC preconditioner.
 * The intent is that the user can overide them and build
 * there own implementation. Though a fairly substantial implementation
 * is provided in <code>InvLSCStrategy</code>.
 *
 * The basics of this method can be found in
 *
 * [1] Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
 *     for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
 *
 * [2] Elman, and Tuminaro, "Boundary Conditions in Approximate Commutator
 *     Preconditioners for the Navier-Stokes Equations," In press (8/2009)?
 *
 * The Least Squares Commuator preconditioner provides a (nearly) Algebraic approximation
 * of the Schur complement of the (Navier-)Stokes system
 *
 * \f$ A = \left[\begin{array}{cc}
 *        F & B^T \\
 *        B & C
 *     \end{array}\right] \f$
 *
 * The approximation to the Schur complement is
 *
 * \f$ C - B F^{-1} B^T \approx (B \hat{Q}_u^{-1} B^T - \gamma C)^{-1}
 *       (B \hat{Q}_u^{-1} F H B^T+C_I) (B H B^T - \gamma C)^{-1}
 *     + C_O \f$.
 *
 * Where \f$\hat{Q}_u\f$ is typically a diagonal approximation of the mass matrix,
 * and \f$H\f$ is an appropriate diagonal scaling matrix (see [2] for details).
 * The scalars \f$\alpha\f$ and \f$\gamma\f$ are chosen to stabilize an unstable
 * discretization (for the case of \f$C\neq 0\f$). If the system is stable then
 * they can be set to \f$0\f$ (see [1] for more details).
 *
 * In order to approximate \f$A\f$ two decompositions can be chosen, a full LU
 * decomposition and a purely upper triangular version. A full LU decomposition
 * requires that the velocity convection-diffusion operator (\f$F\f$) is inverted
 * twice, while an upper triangular approximation requires only a single inverse.
 *
 * The methods of this strategy provide the different pieces. For instance
 * <code>getInvF</code> provides \f$F^{-1}\f$. Similarly there are calls to get
 * the inverses of \f$B \hat{Q}_u^{-1} B^T - \gamma C\f$,
 * \f$B \hat{Q}_u^{-1} B^T - \gamma C\f$, and \f$\hat{Q}_u^{-1}\f$ as well as
 * the \f$H\f$ operator. All these methods are required by the
 * <code>LSCPreconditionerFactory</code>. Additionally there is a
 * <code>buildState</code> method that is called everytime a preconditiner is
 * (re)constructed. This is to allow for any preprocessing neccessary to be
 * handled.
 *
 * The final set of methods help construct a LSCStrategy object, they are
 * primarily used by the parameter list construction inteface. They are
 * more advanced and can be ignored by initial implementations of this
 * class.
 */
class LSCStrategy {
 public:
  virtual ~LSCStrategy() {}

  /** This informs the strategy object to build the state associated
   * with this operator.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   */
  virtual void buildState(BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

  /** Get the inverse of \f$B Q_u^{-1} B^T - \gamma C\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B Q_u^{-1} B^T - \gamma C\f$.
   */
  virtual LinearOp getInvBQBt(const BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

  /** Get the inverse of \f$B H B^T - \gamma C\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B H B^T - \gamma C\f$.
   */
  virtual LinearOp getInvBHBt(const BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

  /** Get the inverse of the \f$F\f$ block.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$F\f$.
   */
  virtual LinearOp getInvF(const BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

#if 0
   /** Get the inverse for stabilizing the whole Schur complement approximation.
     *
     * \param[in] A The linear operator to be preconditioned by LSC.
     * \param[in] state State object for storying reusable information about
     *                  the operator A.
     *
     * \returns The operator to stabilize the whole Schur complement (\f$\alpha D^{-1} \f$).
     */
   virtual LinearOp getInvAlphaD(const BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;
#endif

  /** Get the inverse to stablized stabilizing the Schur complement approximation using
   * a placement on the ``outside''.  That is what is the value for \f$C_O\f$. This quantity
   * may be null.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The operator to stabilize the whole Schur complement (originally \f$\alpha D^{-1}
   * \f$).
   */
  virtual LinearOp getOuterStabilization(const BlockedLinearOp& A,
                                         BlockPreconditionerState& state) const = 0;

  /** Get the inverse to stablized stabilizing the Schur complement approximation using
   * a placement on the ``inside''.  That is what is the value for \f$C_I\f$. This quantity
   * may be null.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The operator to stabilize the whole Schur complement.
   */
  virtual LinearOp getInnerStabilization(const BlockedLinearOp& A,
                                         BlockPreconditionerState& state) const = 0;

  /** Get the inverse mass matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The inverse of the mass matrix \f$Q_u\f$.
   */
  virtual LinearOp getInvMass(const BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

  /** Get the \f$H\f$ scaling matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The \f$H\f$ scaling matrix.
   */
  virtual LinearOp getHScaling(const BlockedLinearOp& A, BlockPreconditionerState& state) const = 0;

  /** Should the approximation of the inverse use a full LDU decomposition, or
   * is a upper triangular approximation sufficient.
   *
   * \returns True if the full LDU decomposition should be used, otherwise
   *          only an upper triangular version is used.
   */
  virtual bool useFullLDU() const = 0;

  /** Tell strategy that this operator is supposed to be symmetric.
   * Behavior of LSC is slightly different for non-symmetric case.
   *
   * \param[in] isSymmetric Is this operator symmetric?
   */
  virtual void setSymmetric(bool isSymmetric) = 0;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& /* pl */,
                                           const InverseLibrary& /* invLib */) {}

  //! For assiting in construction of the preconditioner
  virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const {
    return Teuchos::null;
  }

  //! For assiting in construction of the preconditioner
  virtual bool updateRequestedParameters(const Teuchos::ParameterList& /* pl */) { return true; }

  //! This method sets the request handler for this object
  void setRequestHandler(const Teuchos::RCP<RequestHandler>& rh) { requestHandler_ = rh; }

  //! This method gets the request handler uses by this object
  Teuchos::RCP<RequestHandler> getRequestHandler() const { return requestHandler_; }

 private:
  Teuchos::RCP<RequestHandler> requestHandler_;
};

}  // end namespace NS
}  // end namespace Teko

#endif
