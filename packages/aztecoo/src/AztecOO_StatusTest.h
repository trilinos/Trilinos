
/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef AZTECOO_STATUSTEST_H
#define AZTECOO_STATUSTEST_H

class Epetra_MultiVector;
#include "AztecOO_StatusType.h"
//#include <iostream>
#include <iomanip>

//! AztecOO_StatusTest: A pure virtual class for extending the status testing capabilities of AztecOO.

/* AztecOO_StatusTest is an interface that can be implemented to extend the convergence testing
   capabilities of AztecOO.  Almost any kind of test can be expressed using this mechanism, 
   including composite tests. In this situation, two existing AztecOO_StatusTest objects 
   test1 and test2 can be used to create a new test.  See AztecOO_StatusTestCombo for details

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
  virtual ostream& Print(ostream& stream, int indent = 0) const = 0;
  //@}

  virtual void PrintStatus(ostream& os, AztecOO_StatusType type) const {
    os << setiosflags(ios::left) << setw(13) << setfill('.');
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
    os << setiosflags(ios::right) << setfill(' ');
    return;
  };
};


#endif /* AZTECOO_STATUSTEST_H */
