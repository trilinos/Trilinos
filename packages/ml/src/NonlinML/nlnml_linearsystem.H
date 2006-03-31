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

namespace NLNML
{

class NLNML_LinearSystem : public virtual NOX::Epetra::LinearSystem
{

public:

  NLNML_LinearSystem();

  virtual ~NLNML_LinearSystem() {}

  bool applyJacobian(const NOX::Epetra::Vector& input, 
	             NOX::Epetra::Vector& result) const;
  
  bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
		              NOX::Epetra::Vector& result) const;

  bool applyJacobianInverse(NOX::Parameter::List &params, 
		            const NOX::Epetra::Vector &input, 
		            NOX::Epetra::Vector &result);

  bool applyRightPreconditioning(bool useTranspose,
				 NOX::Parameter::List& params, 
				 const NOX::Epetra::Vector& input, 
				 NOX::Epetra::Vector& result) const;

  Teuchos::RefCountPtr<NOX::Epetra::Scaling> getScaling();

  void resetScaling(const Teuchos::RefCountPtr<NOX::Epetra::Scaling>& s);

  bool computeJacobian(const NOX::Epetra::Vector& x);

  bool createPreconditioner(const NOX::Epetra::Vector& x, 
			    NOX::Parameter::List& p,
			    bool recomputeGraph);

  bool destroyPreconditioner();
  
  bool checkPreconditionerReuse();

  bool isPreconditionerConstructed() const;
  
  bool hasPreconditioner() const;
  
  Teuchos::RefCountPtr<const Epetra_Operator> getJacobianOperator() const;
  
  Teuchos::RefCountPtr<const Epetra_Operator> getGeneratedPrecOperator() const;
  
  Teuchos::RefCountPtr<Epetra_Operator> getGeneratedPrecOperator();
  
  void setJacobianOperatorForSolve(
                const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp);

  void setPrecOperatorForSolve(
                 const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp);


protected:



} ; // class NLNML_LinearSystem
} // namespace NLNML

#endif // ML_NOX_LINEARSYSTEM_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
