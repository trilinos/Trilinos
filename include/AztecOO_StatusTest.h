
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

#ifndef AZTECOO_STATUSTEST_H
#define AZTECOO_STATUSTEST_H

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

class Epetra_MultiVector;
#include "AztecOO_StatusType.h"
#include "Epetra_ConfigDefs.h"
#include <iomanip>

//! AztecOO_StatusTest: A pure virtual class for extending the status testing capabilities of AztecOO.

/* AztecOO_StatusTest is an interface that can be implemented to extend the convergence testing
   capabilities of AztecOO.  Almost any kind of test can be expressed using this mechanism, 
   including composite tests. In this situation, two existing AztecOO_StatusTest objects 
   test1 and test2 can be used to create a new test.  See AztecOO_StatusTestCombo for details

  \warning Presently it is not valid to associate one status test instance with two different AztecOO objects.

*/

class AztecOO_StatusTest {

 public:
  //@{ \name Constructors/destructors.

  //! Constructor
  AztecOO_StatusTest() {};

  //! Destructor
  virtual ~AztecOO_StatusTest() {};
  //@}

  //@{ \name Methods that must be implemented by any AztecOO_StatusTest implementation.
  //! Indicates if residual vector is required by this convergence test.
  /*! If this method returns true, then the ResidualVector argument to the Converged() method will
    defined.  If this method returns false, then the ResidualVector may not be defined when Converged() is
    called.  Some iterative methods do not explicitly construct the residual vector at each iteration.  Thus,
    if this vector is not required, this vector will not need to be constructed if this method returns false.
  */
  virtual bool ResidualVectorRequired() const = 0;

  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  The calling routine may pass in the
    current native residual vector (the one naturally produced as part of the iterative method) or a
    pre-computed estimate of the two-norm of the current residual, or both or neither.  The calling routine
    should also indicate if the solution of the linear problem has been updated to be compatible with
    the residual.  Some methods, such as GMRES do update the solution at each iteration.

    \param CurrentIter (In) Current iteration of iterative method.

    \param CurrentResVector (In) The current residuals of the iterative process.  These values are 
    assumed to be the residuals that are a "natural" by-product of the iterative process.  
    Typically this means they are not explicitly generated and are therefore subject to round-off 
    error.  These values will also reflect the influence of any Any rigorous use of stopping criteria 
    should not
    rely solely on results using this vector.  Instead, tests can be performed using this vector 
    and, after convergence
    is reached with this vector.

    \param CurrentResNormEst (In) If the iterative method can cheaply provide this estimate, 
    as an alternative
    or in addition to providing the CurrentResVector, this value will contain that estimate.  
    The value will be
    set to -1.0 if no estimate is available.

    \param SolutionUpdated (In) If this argument is true, then the solution vector that is part 
    of the Epetra_LinearProblem
    object being solved is consistent with the residual.  Some iterative methods do not keep the 
    solution vector
    updated with the residual at each iteration.  For example GMRES does not generate the solution 
    until the least-
    square problem is solved at the end of the Arnoldi process.


    \return AztecOO_StatusType: Unconverged, Converged or Failed.
  */
  virtual  AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
				 double CurrentResNormEst,
				 bool SolutionUpdated) = 0;

  //! Return the result of the most recent checkStatus call
  virtual  AztecOO_StatusType GetStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual std::ostream& Print(std::ostream& stream, int indent = 0) const = 0;
  //@}

  virtual void PrintStatus(std::ostream& os, AztecOO_StatusType type) const {
    os << std::setiosflags(std::ios::left) << std::setw(13) << std::setfill('.');
    switch (type) {
    case  Failed:
      os << "Failed";
      break;
    case  Converged:
      os << "Converged";
      break;
    case  Unconverged:
    default:
      os << "**";
      break;
    }
    os << std::setiosflags(std::ios::right) << std::setfill(' ');
    return;
  }
};


#endif /* AZTECOO_STATUSTEST_H */
