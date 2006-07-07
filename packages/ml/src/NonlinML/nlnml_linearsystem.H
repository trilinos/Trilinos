/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_LINEARSYSTEM_H
#define NLNML_LINEARSYSTEM_H

#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_MatrixFree.H"

#include "nlnml_coarselevelnoxinterface.H"

using namespace std;
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
