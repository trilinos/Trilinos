/*
// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#ifndef THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "Amesos2_Solver.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass in terms of
 * <tt>Amesos2</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Amesos2_Thyra_adapters_grp
 */
template<typename Scalar>
class Amesos2LinearOpWithSolve: virtual public LinearOpWithSolveBase<Scalar>
{
public:
  using MAT = Tpetra::CrsMatrix<Scalar>;
  using Op = Tpetra::Operator<Scalar>;
  using MV = Tpetra::MultiVector<Scalar>;
  using Solver = ::Amesos2::Solver<MAT, MV>;
  using ConverterT = TpetraOperatorVectorExtraction<Scalar>;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  Amesos2LinearOpWithSolve();

  /** \brief Calls <tt>this->initialize()</tt>. */
  Amesos2LinearOpWithSolve(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP< Solver > &amesos2Solver,
    const EOpTransp amesos2SolverTransp,
    const Scalar amesos2SolverScalar
    );

  /** \brief Initialize after construction. */
  void initialize(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP< Solver > &amesos2Solver
    );

  /** \brief Extract the forward <tt>LinearOpSourceBase<double></tt> object so
   * that it can be modified and remove it from this object.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_fwdOpSrc().get()</tt> before call.
   * <li><tt><tt>this->get_fwdOpSrc().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > get_fwdOp() const;

  /** \brief . */
  Teuchos::RCP<Solver> get_amesos2Solver();

  /** \brief . */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > get_fwdOpSrc() const;

  /** @name Overridden public functions from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief. */
  Teuchos::RCP<const LinearOpBase<Scalar> > clone() const;
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
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
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
  SolveStatus<Scalar> solveImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}

private:

  Teuchos::RCP<const LinearOpBase<Scalar> > fwdOp_;
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc_;
  Teuchos::RCP< Solver > amesos2Solver_;

  void assertInitialized() const;

};

// ///////////////////////////
// Inline members

template<typename Scalar>
inline
Teuchos::RCP<const LinearOpBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::get_fwdOp() const
{
  return fwdOp_;
}

template<typename Scalar>
inline
Teuchos::RCP<typename Amesos2LinearOpWithSolve<Scalar>::Solver>
Amesos2LinearOpWithSolve<Scalar>::get_amesos2Solver()
{
  return amesos2Solver_;
}

template<typename Scalar>
inline
Teuchos::RCP<const LinearOpSourceBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::get_fwdOpSrc() const
{
  return fwdOpSrc_;
}

} // namespace Thyra

#endif	// THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_HPP
