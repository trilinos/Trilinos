/*
//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
*/

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
//	  const Teuchos::RCP < std::vector<double> > & evals_r,
//	  const Teuchos::RCP < std::vector<double> > & evals_i,
//	  const Teuchos::RCP < NOX::Abstract::MultiVector > & evecs_r,
//	  const Teuchos::RCP < NOX::Abstract::MultiVector > & evecs_i
//	  const int & stepId );
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
