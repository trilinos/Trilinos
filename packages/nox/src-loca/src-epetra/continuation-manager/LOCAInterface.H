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

#ifndef LOCA_INTERFACE_H
#define LOCA_INTERFACE_H

// Trilinos headers
#include <Epetra_Vector.h>
#include <LOCA_Epetra.H>

// ParaCont headers
#include "ProblemLOCAPrototype.H"
#include "ContinuationManager.H"

class  LOCAInterface :
//  public LOCA::Epetra::Interface::Required,
  public LOCA::Epetra::Interface::TimeDependent,
  public NOX::Epetra::Interface::Jacobian {

    public:

      // Constructor
      LOCAInterface( Teuchos::RCP <ProblemLOCAPrototype> & aProblem ,
      Teuchos::RCP <ContinuationManager> aContinuationManager);

      // Destructor
      ~LOCAInterface();

      // Compute and return F
      bool computeF(const Epetra_Vector& x, Epetra_Vector& f,
            const NOX::Epetra::Interface::Required::FillType F);

      // Compute an explicit Jacobian
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      // Set a parameter in the user's code.
      void setParameters(const LOCA::ParameterVector& params);

      // Print the solution
      void printSolution (const Epetra_Vector &x, const double conParam);

      // Compute the Shifted matrix needed for linear stability
      virtual bool computeShiftedMatrix (double alpha, double beta,
                                 const Epetra_Vector &x, Epetra_Operator &A);

      virtual void setXdot(const Epetra_Vector& xdot, const double time);

    private:

      // The continuation manager
      Teuchos::RCP <ContinuationManager> continuationManager;

      // The Problem
      Teuchos::RCP <ProblemLOCAPrototype> problem;

      double t; // time, for non-autonomous problems

};

#endif
