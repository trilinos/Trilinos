// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_PresLaplaceLSCStrategy_hpp__
#define __Teko_PresLaplaceLSCStrategy_hpp__

#include "Teko_LSCStrategy.hpp"

namespace Teko {
namespace NS {

class LSCPrecondState;  // forward declration

/** \brief A strategy that takes a single inverse factory and
 *        uses that for all inverses. If no mass matrix is
 *        passed in the diagonal of the 1,1 block is used.
 *
 * A strategy that takes a single inverse factory and uses that
 * for all inverses. Optionally the mass matrix can be passed
 * in, if it is the diagonal is extracted and that is used to
 * form the inverse approximation.
 */
class PresLaplaceLSCStrategy : public LSCStrategy {
 public:
  //! \name Constructors
  //@{
  PresLaplaceLSCStrategy();
  PresLaplaceLSCStrategy(const Teuchos::RCP<InverseFactory> &factory);
  PresLaplaceLSCStrategy(const Teuchos::RCP<InverseFactory> &invFactF,
                         const Teuchos::RCP<InverseFactory> &invFactS);
  //@}

  virtual ~PresLaplaceLSCStrategy() {}

  //! Functions inherited from LSCStrategy
  //@{

  /** This informs the strategy object to build the state associated
   * with this operator.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   */
  virtual void buildState(BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of \f$B Q_u^{-1} B^T\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B Q_u^{-1} B^T\f$.
   */
  virtual LinearOp getInvBQBt(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of \f$B H B^T - \gamma C\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B H B^T - \gamma C\f$.
   */
  virtual LinearOp getInvBHBt(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of the \f$F\f$ block.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$F\f$.
   */
  virtual LinearOp getInvF(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse for stabilizing the whole schur complement approximation.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The operator to stabilize the whole Schur complement.
   */
  // virtual LinearOp getInvAlphaD(const BlockedLinearOp & A,BlockPreconditionerState & state)
  // const;
  virtual LinearOp getOuterStabilization(const BlockedLinearOp &A,
                                         BlockPreconditionerState &state) const;

  virtual LinearOp getInnerStabilization(const BlockedLinearOp & /* A */,
                                         BlockPreconditionerState & /* state */) const {
    return Teuchos::null;
  }

  /** Get the inverse mass matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The inverse of the mass matrix \f$Q_u\f$.
   */
  virtual LinearOp getInvMass(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the \f$H\f$ scaling matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The \f$H\f$ scaling matrix.
   */
  virtual LinearOp getHScaling(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Should the approximation of the inverse use a full LDU decomposition, or
   * is a upper triangular approximation sufficient.
   *
   * \returns True if the full LDU decomposition should be used, otherwise
   *          only an upper triangular version is used.
   */
  virtual bool useFullLDU() const { return useFullLDU_; }

  /** Tell strategy that this operator is supposed to be symmetric.
   * Behavior of LSC is slightly different for non-symmetric case.
   *
   * \param[in] isSymmetric Is this operator symmetric?
   */
  virtual void setSymmetric(bool isSymmetric) { isSymmetric_ = isSymmetric; }

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList &pl,
                                           const InverseLibrary &invLib);

  //! For assiting in construction of the preconditioner
  virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

  //! For assiting in construction of the preconditioner
  virtual bool updateRequestedParameters(const Teuchos::ParameterList &pl);
  //@}

  //! Initialize the state object using this blocked linear operator
  virtual void initializeState(const BlockedLinearOp &A, LSCPrecondState *state) const;

  /** Compute the inverses required for the LSC Schur complement
   *
   * \note This method assumes that the BQBt and BHBt operators have
   *       been constructed.
   */
  void computeInverses(const BlockedLinearOp &A, LSCPrecondState *state) const;

  //! Set the number of power series iterations to use when finding the spectral radius
  virtual void setEigSolveParam(int sz) { eigSolveParam_ = sz; }

  //! Return the number of power series iterations to use when finding the spectral radius
  virtual int getEigSolveParam() { return eigSolveParam_; }

  //! Set to true to use the Full LDU decomposition, false otherwise
  virtual void setUseFullLDU(bool val) { useFullLDU_ = val; }

 protected:
  // how to invert the matrices
  Teuchos::RCP<InverseFactory> invFactoryV_;
  Teuchos::RCP<InverseFactory> invFactoryP_;

  // flags for handling various options
  bool isSymmetric_;
  int eigSolveParam_;
  bool useFullLDU_;

  // scaling operator parameters
  bool useMass_;
  DiagonalType scaleType_;

 private:
  PresLaplaceLSCStrategy(const PresLaplaceLSCStrategy &);

 public:
  // some static functions for determining strings

  static std::string getPressureLaplaceString() { return "Pressure Laplace Operator"; }
  static std::string getVelocityMassString() { return "Velocity Mass Operator"; }
};

}  // end namespace NS
}  // end namespace Teko

#endif
