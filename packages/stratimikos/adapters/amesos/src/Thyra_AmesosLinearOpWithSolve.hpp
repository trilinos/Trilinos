// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Epetra_LinearProblem.h"
#include "Amesos_BaseSolver.h"


namespace Thyra {


/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass that adapts any
 * <tt>Amesos_BaseSolver</tt> object.
 *
 * See the <tt>LinearOpWithSolveBase</tt> interface for a description of how
 * to use objects of this type.
 *
 * <b>Note:</b> Clients should not generally directly create objects of this
 * type but instead should use <tt>AmesosLinearOpWithSolveFactory</tt>.  Only
 * very sophisticated users should ever directly interact with an object
 * through this subclass interface.
 *
 * \ingroup Amesos_Thyra_adapters_grp
 */
class AmesosLinearOpWithSolve : virtual public LinearOpWithSolveBase<double>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AmesosLinearOpWithSolve();

  /** \brief Calls <tt>this->initialize()</tt>. */
  AmesosLinearOpWithSolve(
    const Teuchos::RCP<const LinearOpBase<double> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
    const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
    const EOpTransp amesosSolverTransp,
    const double amesosSolverScalar
    );

  /** \brief First initialization.
   *
   * \param  fwdOp
   *           [in] The forward operator for which the factorization exists.
   * \param  epetraLP
   *           [in] The <tt>Epetra_LinearProblem</tt> object that was used to
   *           create the <tt>Amesos_BaseSolver</tt> object
   *           <tt>*amesosSolver</tt>.  Note that the RHS and the LHS
   *           multi-vector pointers in this object will be set and unset
   *           here.
   * \param  amesosSolver
   *           [in] Contains the factored, and ready to go,
   *           <tt>Amesos_BaseSolver</tt> object ready to solve linear system.
   * \param  amesosSolverTransp
   *           [in] Determines if the %Amesos solver should be used as its
   *           transpose or not.
   * \param  amesosSolverScalar
   *           [in] Determines the scaling factor associated with the %Amesos
   *           solver.  The solution to the linear solve is scaled by
   *           <tt>1/amesosSolverScalar</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>epetraLP.get()!=NULL</tt>
   * <li><tt>amesosSolver.get()!=NULL</tt>
   * <li><tt>*epetraLP->GetOperator()</tt> is compatible with <tt>*fwdOp</tt>
   * <li><tt>epetraLP->GetLHS()==NULL</tt>
   * <li><tt>epetraLP->GetRHS()==NULL</tt>
   * <li><tt>*amesosSolver</tt> contains the factorization of <tt>*fwdOp</tt> and is
   *     ready to solve linear systems!
   * </ul>
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_fwdOp().get() == fwdOp.get()</tt>
   * <li><tt>this->get_epetraLP().get() == epetraLP.get()</tt>
   * <li><tt>this->get_amesosSolver().get() == amesosSolver.get()</tt>
   * <li><tt>this->get_amesosSolverTransp() == amesosSolverTransp</tt>
   * <li><tt>this->get_amesosSolverScalar() == amesosSolverScalar</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const LinearOpBase<double> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
    const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
    const EOpTransp amesosSolverTransp,
    const double amesosSolverScalar
    );

  /** \brief Extract the <tt>LinearOpSourceBase<double></tt> object so that it can be modified.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_fwdOpSrc().get()</tt> before call.
   * <li><tt><tt>this->get_fwdOpSrc().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const LinearOpSourceBase<double> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<double> > get_fwdOp() const;

  /** \brief . */
  Teuchos::RCP<const LinearOpSourceBase<double> > get_fwdOpSrc() const;

  /** \brief . */
  Teuchos::RCP<Epetra_LinearProblem> get_epetraLP() const;

  /** \brief . */
  Teuchos::RCP<Amesos_BaseSolver> get_amesosSolver() const;

  /** \brief . */
  EOpTransp get_amesosSolverTransp() const;

  /** \brief . */
  double get_amesosSolverScalar() const;

  /** \brief Uninitialize.
   */
  void uninitialize(
    Teuchos::RCP<const LinearOpBase<double> > *fwdOp = NULL,
    Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOpSrc = NULL,
    Teuchos::RCP<Epetra_LinearProblem> *epetraLP = NULL,
    Teuchos::RCP<Amesos_BaseSolver> *amesosSolver = NULL,
    EOpTransp *amesosSolverTransp = NULL,
    double *amesosSolverScalar = NULL
    );
  
  //@}

  /** @name Overridden public functions from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<double> > range() const;
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<double> > domain() const;
  /** \brief. */
  Teuchos::RCP<const LinearOpBase<double> > clone() const;
  //@}

  /** @name Overridden public functions from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  //@}

protected:

  /** @name Overridden from LinearOpBase  */
  //@{
  /** \brief . */
  virtual bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X,
    const Ptr<MultiVectorBase<double> > &Y,
    const double alpha,
    const double beta
    ) const;
  //@}

  /** @name Overridden from LinearOpWithSolveBase. */
  //@{
  /** \brief . */
  virtual bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  SolveStatus<double> solveImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &B,
    const Ptr<MultiVectorBase<double> > &X,
    const Ptr<const SolveCriteria<double> > solveCriteria
    ) const;
  //@}

private:

  Teuchos::RCP<const LinearOpBase<double> > fwdOp_;
  Teuchos::RCP<const LinearOpSourceBase<double> > fwdOpSrc_;
  Teuchos::RCP<Epetra_LinearProblem> epetraLP_;
  Teuchos::RCP<Amesos_BaseSolver> amesosSolver_;
  EOpTransp amesosSolverTransp_;
  double amesosSolverScalar_;

  void assertInitialized() const;

};

// ///////////////////////////
// Inline members

inline
Teuchos::RCP<const LinearOpBase<double> >
AmesosLinearOpWithSolve::get_fwdOp() const
{
  return fwdOp_;
}

inline
Teuchos::RCP<const LinearOpSourceBase<double> >
AmesosLinearOpWithSolve::get_fwdOpSrc() const
{
  return fwdOpSrc_;
}

inline
Teuchos::RCP<Epetra_LinearProblem>
AmesosLinearOpWithSolve::get_epetraLP() const
{
  return epetraLP_;
}

inline
Teuchos::RCP<Amesos_BaseSolver>
AmesosLinearOpWithSolve::get_amesosSolver() const
{
  return amesosSolver_;
}

inline
EOpTransp AmesosLinearOpWithSolve::get_amesosSolverTransp() const
{
  return amesosSolverTransp_;
}

inline
double AmesosLinearOpWithSolve::get_amesosSolverScalar() const
{
  return amesosSolverScalar_;
}

} // namespace Thyra

#endif	// THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
