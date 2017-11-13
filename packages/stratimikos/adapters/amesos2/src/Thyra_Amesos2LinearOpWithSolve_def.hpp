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

#include "Thyra_Amesos2LinearOpWithSolve_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <Teuchos_RCP.hpp>


namespace Thyra {


// Constructors/initializers/accessors


template<typename Scalar>
Amesos2LinearOpWithSolve<Scalar>::Amesos2LinearOpWithSolve()
{}


template<typename Scalar>
Amesos2LinearOpWithSolve<Scalar>::Amesos2LinearOpWithSolve(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const Teuchos::RCP< Solver > &amesos2Solver,
  const EOpTransp amesos2SolverTransp,
  const Scalar amesos2SolverScalar
  )
{
  this->initialize(fwdOp,fwdOpSrc,amesos2Solver);
}


template<typename Scalar>
void Amesos2LinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const Teuchos::RCP< Solver > &amesos2Solver
  )
{
  fwdOp_ = fwdOp;
  fwdOpSrc_ = fwdOpSrc;
  amesos2Solver_ = amesos2Solver;
  const std::string fwdOpLabel = fwdOp_->getObjectLabel();
  if(fwdOpLabel.length())
    this->setObjectLabel( "lows("+fwdOpLabel+")" );
}


template<typename Scalar>
Teuchos::RCP<const LinearOpSourceBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::extract_fwdOpSrc()
{
  Teuchos::RCP<const LinearOpSourceBase<Scalar> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}


// Overridden from LinearOpBase


template<typename Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::range() const
{
  return ( fwdOp_.get() ? fwdOp_->range() : Teuchos::null );
}


template<typename Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::domain() const
{
  return  ( fwdOp_.get() ? fwdOp_->domain() : Teuchos::null );
}


template<typename Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}


// Overridden from Teuchos::Describable


template<typename Scalar>
std::string Amesos2LinearOpWithSolve<Scalar>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if(!is_null(amesos2Solver_)) {
    oss
      << "{fwdOp="<<fwdOp_->description()
      << ",amesos2Solver="<<typeName(*amesos2Solver_)<<"}";
  }
  return oss.str();
}


template<typename Scalar>
void Amesos2LinearOpWithSolve<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::typeName;
  using Teuchos::describe;
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab(out);
      if(!is_null(fwdOp_)) {
        out << "fwdOp = " << describe(*fwdOp_,verbLevel);
      }
      if(!is_null(amesos2Solver_)) {
        out << "amesos2Solver=" << typeName(*amesos2Solver_) << "\n";
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from LinearOpBase


template<typename Scalar>
bool Amesos2LinearOpWithSolve<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return ::Thyra::opSupported(*fwdOp_,M_trans);
}


template<typename Scalar>
void Amesos2LinearOpWithSolve<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply( *fwdOp_, M_trans, X, Y, alpha, beta );
}


// Overridden from LinearOpWithSolveBase

template<typename Scalar>
bool Amesos2LinearOpWithSolve<Scalar>::solveSupportsImpl(EOpTransp M_trans) const
{
  if (Thyra::real_trans(M_trans) == Thyra::NOTRANS) {
    // Assume every amesos2 solver supports a basic forward solve!
    return true;
  }
  // dai: not sure how to query amesos2 transpose capability
  return false;
}


template<typename Scalar>
bool Amesos2LinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  return true; // I am a direct solver so I should be able to do it all!
}


template<typename Scalar>
SolveStatus<Scalar>
Amesos2LinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{
  auto Btpetra = ConverterT::getConstTpetraMultiVector(Teuchos::rcpFromRef(B));

  auto Xtpetra = ConverterT::getTpetraMultiVector(Teuchos::rcpFromPtr(X));

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nSolving system using Amesos2 solver "
         << typeName(*amesos2Solver_) << " ...\n\n";

  amesos2Solver_->solve(Xtpetra.ptr(), Btpetra.ptr());

  SolveStatus<Scalar> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_CONVERGED; 
  return solveStatus;
}


}	// end namespace Thyra
