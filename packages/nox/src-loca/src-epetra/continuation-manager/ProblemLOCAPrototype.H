// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PROBLEM_LOCA_PROTOTYPE_H
#define PROBLEM_LOCA_PROTOTYPE_H

#include "ProblemNOXPrototype.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Abstract_MultiVector.H"

class ProblemLOCAPrototype :
  public ProblemNOXPrototype {

    public:

      // Constructor
      ProblemLOCAPrototype();

      // Destructor
      virtual ~ProblemLOCAPrototype();

      // Pointer to the continuation parameter list
      virtual LOCA::ParameterVector GetContinuableParams() const =0 ;

      // Setting one continuation parameter
      virtual bool SetContinuableParameter(std::string label,double value)=0;

      // Printing the solution at each step
      virtual bool PrintSolutionFile(const std::string & fileName, const Epetra_Vector & x,
      const Teuchos::ParameterList & xParams);

      // Setting the solution file parameters
      virtual bool SetSolutionFileParameters(const Epetra_Vector & x);

      // Getting the solution file parameters
      virtual Teuchos::RCP <Teuchos::ParameterList> GetSolutionFileParameters();

      // Updates the continuation file
      virtual bool UpdateContinuationFile( const std::string & fileName,
          const int & idStep,
      const Teuchos::ParameterList & continuationFileParams);

      // Setting the continuation file parameters
      virtual bool SetContinuationFileParameters(const Epetra_Vector & x);

      // Getting the continuation file parameters
      virtual Teuchos::RCP <Teuchos::ParameterList>
    GetContinuationFileParameters();

      // Compute derivative with respect to periodic direction
      virtual bool ComputePeriodicDirectionDerivative(const Epetra_Vector & x,
      Epetra_Vector & dx);

//      // Printing eigenvalues at each step
//      virtual bool PrintEigenValues(
//      const Teuchos::RCP < std::vector<double> > & evals_r,
//      const Teuchos::RCP < std::vector<double> > & evals_i,
//      const Teuchos::RCP < NOX::Abstract::MultiVector > & evecs_r,
//      const Teuchos::RCP < NOX::Abstract::MultiVector > & evecs_i
//      const int & stepId );
//
//      // Compute derivative with respect to time
//      virtual bool ComputeTimeDerivative(const Epetra_Vector & x,
//          Epetra_Vector & dxdt);
//
//      // Compute derivative with respect to space
//      virtual bool ComputeSpaceDerivative(const Epetra_Vector & x,
//          Epetra_Vector & dxdt);
//
//      // Compute a shifted version of the Jacobian, namely A = alpha*J + beta*M
//      // necessary for linear stability
//      virtual bool ComputeShiftedJacobian (double alpha, double beta) = 0;
//
//      // Return a pointer to the shifted jacobian
//      virtual Teuchos::RCP< Epetra_CrsMatrix > GetShiftedJacF() const = 0;
//
};
#endif
