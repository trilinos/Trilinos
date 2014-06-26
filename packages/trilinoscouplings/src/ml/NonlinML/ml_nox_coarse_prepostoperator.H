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
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
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
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NOX_COARSE_PREPOSTOPERATOR_H
#define NOX_COARSE_PREPOSTOPERATOR_H

#include <iostream>

// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX headers
#include "NOX_Utils.H"          
#include "NOX_Common.H"         
#include "NOX_Solver_Generic.H"	
#include "NOX_Abstract_Group.H"	
#include "NOX_Epetra_Vector.H"	
#include "NOX_EpetraNew_LinearSystem.H"
#include "NOX_Parameter_PrePostOperator.H"
#include "NOX_EpetraNew_LinearSystemAztecOO.H"
#include "NOX_EpetraNew_Group.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_Solver_Manager.H"

// ML-headers
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"

// User Defined
#include "ml_nox_fineinterface.H"
#include "ml_nox_coarseinterface.H"

// forward declarations

namespace ML_NOX
{


class Ml_Nox_CoarsePrePostOperator : public NOX::Parameter::PrePostOperator {

public:
  
  //! %Abstract %Vector constructor (does nothing)
  Ml_Nox_CoarsePrePostOperator(ML_NOX::Nox_CoarseProblem_Interface& coarseinterface,
                               ML_NOX::Ml_Nox_Fineinterface& fineinterface);

  //! Copy constructor (doesnothing)
  Ml_Nox_CoarsePrePostOperator(const ML_NOX::Ml_Nox_CoarsePrePostOperator& source);

  //! %Abstract %Vector destructor (does nothing)
  virtual ~Ml_Nox_CoarsePrePostOperator();

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
  void runPreIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
  void runPostIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  void runPreSolve(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  void runPostSolve(const NOX::Solver::Generic& solver);

  const string& getType() const { return type_; }
  NOX::Parameter::Arbitrary* clone() const;
  
  
private:

  string                               type_;
  ML_NOX::Nox_CoarseProblem_Interface& coarseinterface_;
  ML_NOX::Ml_Nox_Fineinterface&        fineinterface_;

}; // class Ml_Nox_CoarsePrePostOperator


} // namespace ML_NOX
#endif // NOX_COARSE_PREPOSTOPERATOR_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
