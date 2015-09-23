//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef AZTECOO_CONDITIONNUMBER_H
#define AZTECOO_CONDITIONNUMBER_H

class Epetra_Map;
class Epetra_Vector;
class Epetra_Operator;
class AztecOO;

/*!

 \brief Condition number estimator using AztecOO.

 This object will estimate the condition number of an Epetra_Operator.
 */
class AztecOOConditionNumber {

 public:
  
  //! Solver type to use.
  enum SolverType {
    //! CG for symmetric matrices
    CG_, 
    //! GMRES for nonsymmetric
    GMRES_
  };
  
  //! Constructor.  
  AztecOOConditionNumber();
  
  //! Destructor
  ~AztecOOConditionNumber();
  
  //! Initialization
  void initialize(const Epetra_Operator& op,
		  SolverType solverType=GMRES_,
		  int krylovSubspaceSize=100, 
		  bool printSolve=false);
  
  //! Estimates the condition number. 
  int computeConditionNumber(int maxIters, double tol);

  //! Return condition number computed by last call to computeConditionNumber.
  double getConditionNumber();

 protected:

  //! Frees all memory allocated with new by this object.
  void freeMemory();

 protected:

  //! Condition number calculated in computeConditionNumber.
  double conditionNumber_;

  //! Map to create left hand side vector.  
  Epetra_Map* domainMap_;

  //! Map to create right hand side vector.  
  Epetra_Map* rangeMap_;

  //! Operator supplied by user in initialization.
  Epetra_Operator* operator_;

  //! RHS vector.  This is initializaed to a random vector.
  Epetra_Vector* rhs_;

  //! Dummy vector.  Initializaed to zero.
  Epetra_Vector* dummy_;

  //! solver object.
  AztecOO* solver_;

  //! Conditional for printing solve to output.
  bool printSolve_;

};


#endif
