
/*
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
*/

#ifndef AZTECOO_STATUSTESTMAXITERS_H
#define AZTECOO_STATUSTESTMAXITERS_H

#include "AztecOO_StatusTest.h"
class Epetra_MultiVector;

//! AztecOO_StatusTestMaxIters: An AztecOO_StatusTest class specifying a maximum number of iterations.

/* This implementation of the AztecOO_StatusTest base class tests the number of iterations performed
   against a maximum number allowed.

  \warning Presently it is not valid to associate one status test instance with two different AztecOO objects.
*/

class AztecOO_StatusTestMaxIters: public AztecOO_StatusTest {

 public:

  //@{ \name Constructors/destructors.
  //! Constructor
  AztecOO_StatusTestMaxIters(int MaxIters);

  //! Destructor
  virtual ~AztecOO_StatusTestMaxIters() {};
  //@}

  //@{ \name Methods that implement the AztecOO_StatusTest interface.

  //! Indicates if residual vector is required by this convergence test: returns false for this class.
  bool ResidualVectorRequired() const {return(false);} ;

  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met..

    \param CurrentIter (In) Current iteration of iterative method.  Compared against MaxIters value
    passed in at construction.  If CurrentIter < MaxIters, we return with StatusType = Unconverged.
    Otherwise, StatusType will be set to Failed.

    \param CurrentResVector (In) Ignored by this class.

    \param CurrentResNormEst (In) Ignored by this class.

    \param SolutionUpdated (In) Ignored by this class.


    \return StatusType Unconverged if CurrentIter<MaxIters, Failed if CurrentIters>=MaxIters.
  */
  AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
					 double CurrentResNormEst,
				 bool SolutionUpdated);
  AztecOO_StatusType GetStatus() const {return(status_);};

  ostream& Print(ostream& stream, int indent = 0) const;
  //@}
  
  //@{ \name Methods to access data members.

  //! Returns the maximum number of iterations set in the constructor.
  int GetMaxIters() const {return(MaxIters_);};

  //! Returns the current number of iterations from the most recent StatusTest call.
  int GetNumIters() const {return(Niters_);};

  //@}

private:

  //@{ \name Private data members.
  //! Maximum number of iterations allowed
  int MaxIters_;

  //! Current number of iterations
  int Niters_;

  //! Status
  AztecOO_StatusType status_;
  //@}

};

#endif /* AZTECOO_STATUSTESTMAXITERS_H */
