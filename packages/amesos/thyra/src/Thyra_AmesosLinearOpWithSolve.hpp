/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_SingleScalarLinearOpWithSolveBase.hpp"
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
class AmesosLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<double>                  // Public interface
  , virtual protected SingleScalarLinearOpWithSolveBase<double>   // Implementation detail
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AmesosLinearOpWithSolve();

  /** \brief Calls <tt>this->initialize()</tt>. */
  AmesosLinearOpWithSolve(
    const Teuchos::RefCountPtr<const LinearOpBase<double> > &fwdOp,
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RefCountPtr<Epetra_LinearProblem> &epetraLP,
    const Teuchos::RefCountPtr<Amesos_BaseSolver> &amesosSolver,
    const ETransp amesosSolverTransp,
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
    const Teuchos::RefCountPtr<const LinearOpBase<double> > &fwdOp,
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RefCountPtr<Epetra_LinearProblem> &epetraLP,
    const Teuchos::RefCountPtr<Amesos_BaseSolver> &amesosSolver,
    const ETransp amesosSolverTransp,
    const double amesosSolverScalar
    );

  /** \brief Extract the <tt>LinearOpSourceBase<double></tt> object so that it can be modified.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_fwdOpSrc().get()</tt> before call.
   * <li><tt><tt>this->get_fwdOpSrc().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr<const LinearOpSourceBase<double> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<double> > get_fwdOp() const;

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpSourceBase<double> > get_fwdOpSrc() const;

  /** \brief . */
  Teuchos::RefCountPtr<Epetra_LinearProblem> get_epetraLP() const;

  /** \brief . */
  Teuchos::RefCountPtr<Amesos_BaseSolver> get_amesosSolver() const;

  /** \brief . */
  ETransp get_amesosSolverTransp() const;

  /** \brief . */
  double get_amesosSolverScalar() const;

  /** \brief Uninitialize.
   */
  void uninitialize(
    Teuchos::RefCountPtr<const LinearOpBase<double> > *fwdOp = NULL,
    Teuchos::RefCountPtr<const LinearOpSourceBase<double> > *fwdOpSrc = NULL,
    Teuchos::RefCountPtr<Epetra_LinearProblem> *epetraLP = NULL,
    Teuchos::RefCountPtr<Amesos_BaseSolver> *amesosSolver = NULL,
    ETransp *amesosSolverTransp = NULL,
    double *amesosSolverScalar = NULL
    );
  
  //@}

  /** @name Overridden public functions from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > range() const;
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > domain() const;
  /** \brief. */
  Teuchos::RefCountPtr<const LinearOpBase<double> > clone() const;
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

  /** @name Overridden protected functions from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  /** \brief . */
  void apply(
    const ETransp M_trans,
    const MultiVectorBase<double> &X,
    MultiVectorBase<double> *Y,
    const double alpha,
    const double beta
    ) const;
  //@}

  /** @name Overridden protected functions from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(
    ETransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  void solve(
    const ETransp M_trans,
    const MultiVectorBase<double> &B,
    MultiVectorBase<double> *X,
    const int numBlocks,
    const BlockSolveCriteria<double> blockSolveCriteria[],
    SolveStatus<double> blockSolveStatus[]
    ) const;
  //@}

private:

  Teuchos::RefCountPtr<const LinearOpBase<double> > fwdOp_;
  Teuchos::RefCountPtr<const LinearOpSourceBase<double> > fwdOpSrc_;
  Teuchos::RefCountPtr<Epetra_LinearProblem> epetraLP_;
  Teuchos::RefCountPtr<Amesos_BaseSolver> amesosSolver_;
  ETransp amesosSolverTransp_;
  double amesosSolverScalar_;

  void assertInitialized() const;

};

// ///////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<const LinearOpBase<double> >
AmesosLinearOpWithSolve::get_fwdOp() const
{
  return fwdOp_;
}

inline
Teuchos::RefCountPtr<const LinearOpSourceBase<double> >
AmesosLinearOpWithSolve::get_fwdOpSrc() const
{
  return fwdOpSrc_;
}

inline
Teuchos::RefCountPtr<Epetra_LinearProblem>
AmesosLinearOpWithSolve::get_epetraLP() const
{
  return epetraLP_;
}

inline
Teuchos::RefCountPtr<Amesos_BaseSolver>
AmesosLinearOpWithSolve::get_amesosSolver() const
{
  return amesosSolver_;
}

inline
ETransp AmesosLinearOpWithSolve::get_amesosSolverTransp() const
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
