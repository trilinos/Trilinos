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

#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

// Trilinos headers
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

// ParaCont headers
#include "ProblemLOCAPrototype.H"

class LinearSystem:
  public ProblemLOCAPrototype{

    public:

      // Constructor
#ifdef HAVE_MPI
      LinearSystem( const Teuchos::RefCountPtr <Epetra_MpiComm> & aComm );
#else
      LinearSystem( const Teuchos::RefCountPtr <Epetra_SerialComm> & aComm );
#endif

      // Destructor
      ~LinearSystem();

      // Implementation of ProblemNOXPrototype
      // virtual methods  ---------------------------------------
      //
      // Compute F
      virtual bool ComputeF(const Epetra_Vector & x, 
          Epetra_Vector & f);
      
      // Compute the Jacobian of F
      virtual bool ComputeJacF(const Epetra_Vector & x);
      
      // Returns a pointer to the Jacobian of F
      virtual Teuchos::RefCountPtr <Epetra_CrsMatrix> 
	GetJacF() const;

      // Returns a pointer to the initial guess
      virtual Teuchos::RefCountPtr <Epetra_Vector> 
	GetInitialGuess() const;

      // Implementation of ProblemLOCAPrototype
      // virtual methods  --------------------------------------
      //                                                         
      // Pointer to the continuable parameter list
      virtual LOCA::ParameterVector 
	GetContinuableParams() const;

      // Setting one continuation parameter
      virtual bool SetContinuableParameter(std::string label,
	  double value);

      // Updates the continuation file
      virtual bool UpdateContinuationFile( 
	  const std::string & fileName,
          const int & idStep,
          const Teuchos::ParameterList & continuationFileParams);
  
      // Setting the continuation file parameters
      virtual bool SetContinuationFileParameters(const Epetra_Vector & x);
  
      // Getting the continuation file parameters
      virtual Teuchos::RefCountPtr <Teuchos::ParameterList> 
	GetContinuationFileParameters();

      // Printing the solution at each step
      virtual bool PrintSolutionFile(const std::string & fileName, 
	  const Epetra_Vector & x,
	  const Teuchos::ParameterList & xParams);

    private:

      // Internally stored variables ----------------------
      //
      // Communicator 
#ifdef HAVE_MPI
      Teuchos::RefCountPtr <Epetra_MpiComm> comm;
#else
      Teuchos::RefCountPtr <Epetra_SerialComm> comm;
#endif

      // Continuable parameter list
      LOCA::ParameterVector continuableParams;

      // The continuation file parameters
      Teuchos::RefCountPtr <Teuchos::ParameterList> 
	continuationFileParams;

      // Initial Guess
      Teuchos::RefCountPtr <Epetra_Vector> initialGuess;

      // Jacobian
      Teuchos::RefCountPtr <Epetra_CrsMatrix>  jacobian;

      // From Local to global indices (vector map)
      int * myGlobalElements;

      // Local number of elements (vector map)
      int numMyElements;

      // Epetra Map of the vector's components
      Teuchos::RefCountPtr <Epetra_Map>  vectorMap;
  };

#endif
