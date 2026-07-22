// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_LINEARSYSTEM_H
#define NLNML_LINEARSYSTEM_H

#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_MatrixFree.H"

#include "nlnml_coarselevelnoxinterface.H"

using namespace Teuchos;
namespace NLNML
{

class NLNML_LinearSystem : public virtual NOX::Epetra::LinearSystem
{

public:

  NLNML_LinearSystem(
              RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac,
              RefCountPtr<Epetra_Operator> J,
              RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec,
              RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
              RefCountPtr<Epetra_Operator> Prec,
              //const Epetra_Vector& soln,
              bool ismatrixfree,
              int level,
              int printlevel);

  virtual ~NLNML_LinearSystem() {}

  bool applyJacobian(const NOX::Epetra::Vector& input,
	             NOX::Epetra::Vector& result) const;

  bool applyJacobianTranspose(const NOX::Epetra::Vector& input,
		              NOX::Epetra::Vector& result) const;

  bool applyJacobianInverse(Teuchos::ParameterList &params,
		            const NOX::Epetra::Vector &input,
		            NOX::Epetra::Vector &result);

  bool applyRightPreconditioning(bool useTranspose,
				 Teuchos::ParameterList& params,
				 const NOX::Epetra::Vector& input,
				 NOX::Epetra::Vector& result) const;

  Teuchos::RefCountPtr<NOX::Epetra::Scaling> getScaling() { return null; }

  void resetScaling(const Teuchos::RefCountPtr<NOX::Epetra::Scaling>& s) { return; }

  bool computeJacobian(const NOX::Epetra::Vector& x);

  bool createPreconditioner(const NOX::Epetra::Vector& x,
			    Teuchos::ParameterList& p,
			    bool recomputeGraph) const;

  // Because of reuse features I prefer to do this myself
  bool destroyPreconditioner() const { return true; }

  // To be consistent with preventing preconditioner re-use, we flag this option
  // as invalid for now.  RWH

  bool recomputePreconditioner(const NOX::Epetra::Vector& x,
			Teuchos::ParameterList& linearSolverParams) const
  { return false; }

  // For now, we do not support reuse of the preconditioner.  This could change in th future. RWH
  PreconditionerReusePolicyType
  getPreconditionerPolicy(bool advanceReuseCounter)
  { return PRPT_REBUILD; }

  // Because of reuse features I prefer to do this myself
  bool checkPreconditionerReuse() { return true; }

  // Always...
  bool isPreconditionerConstructed() const { return true; }

  // Always...
  bool hasPreconditioner() const { return true; }

  Teuchos::RefCountPtr<const Epetra_Operator> getJacobianOperator() const
    { return J_; }

  Teuchos::RefCountPtr<Epetra_Operator> getJacobianOperator()
    { return J_; }

  Teuchos::RefCountPtr<const Epetra_Operator> getGeneratedPrecOperator() const
    { return Prec_; }

  Teuchos::RefCountPtr<Epetra_Operator> getGeneratedPrecOperator()
    { return Prec_; }

  void setJacobianOperatorForSolve(
                const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp);

  void setPrecOperatorForSolve(
                 const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp);


protected:

RefCountPtr<NOX::Epetra::Interface::Jacobian>       iJac_;
RefCountPtr<Epetra_Operator>                        J_;
RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec_;
RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface>   coarseinterface_;
RefCountPtr<Epetra_Operator>                        Prec_;
//const Epetra_Vector&                                soln_;
bool                                                ismatrixfree_;
int                                                 level_;
int                                                 printlevel_;


} ; // class NLNML_LinearSystem
} // namespace NLNML

#endif // ML_NOX_LINEARSYSTEM_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
