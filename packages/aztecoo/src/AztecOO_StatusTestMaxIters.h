
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

#ifndef AZTECOO_STATUSTESTMAXITERS_H
#define AZTECOO_STATUSTESTMAXITERS_H

#include "AztecOO_StatusTest.h"
class Epetra_MultiVector;

//! AztecOO_StatusTestMaxIters: An AztecOO_StatusTest class specifying a maximum number of iterations.

/* This implementation of the AztecOO_StatusTest base class tests the number of iterations performed
   against a maximum number allowed.
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
  virtual bool ResidualVectorRequired() const {return(false);} ;

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
  virtual AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
					 double CurrentResNormEst,
				 bool SolutionUpdated);
  virtual AztecOO_StatusType GetStatus() const;

  virtual ostream& Print(ostream& stream, int indent = 0) const;
  //@}
  
  //@{ \name Methods to access data members.

  //! Returns the maximum number of iterations set in the constructor.
  virtual int GetMaxIters() const;

  //! Returns the current number of iterations from the most recent StatusTest call.
  virtual int GetNumIters() const;

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
