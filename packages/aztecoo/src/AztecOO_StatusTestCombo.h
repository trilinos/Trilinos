
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

#ifndef AZTECOO_STATUSTESTCOMBO_H
#define AZTECOO_STATUSTESTCOMBO_H

#include "AztecOO_StatusTest.h"
#include <vector>
class Epetra_MultiVector;

//! AztecOO_StatusTestCombo: A pure virtual class for extending the status testing capabilities of AztecOO.

/* AztecOO_StatusTestCombo is an interface that can be implemented to extend the convergence testing
   capabilities of AztecOO.  This class supports composite tests.  In this situation,
   two existing AztecOO_StatusTestCombo objects test1 and test2 can be used to create a new test.

*/

class AztecOO_StatusTestCombo: public AztecOO_StatusTest {

 public:

  //@{ \name Enums.
  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests.
  */
  enum ComboType {AND, /*!< Require both subtests to be satisfied. */
		  OR   /*!< Require one or the other subtests to be satisfied. */
  };
  //@}

  //@{ \name Constructors/destructors.
  //! Constructor
  AztecOO_StatusTestCombo(ComboType t);

  //! Constructor with a single test.
  AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a);

  //! Constructor with two tests.
  AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b);

  //! Add another test to this combination.
  virtual AztecOO_StatusTestCombo& AddStatusTest(AztecOO_StatusTest& a);

  //! Destructor
  virtual ~AztecOO_StatusTestCombo() {};
  //@}

  //@{ \name Methods that implement the AztecOO_StatusTest interface.
  //! Indicates if residual vector is required by this convergence test.
  /*! If this method returns true, then one or more of the AztecOO_StatusTest objects that make up this combined
    test requires the Residual Vector to perform its test.
  */
  virtual bool ResidualVectorRequired() const;

  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  Depending on how the combined test
    is constructed this method will return the appropriate status type using common logic principals.  
    However, if any subtest returns with a Failed status type, the combined test will return a status
    type of Failed.

    \param CurrentIter (In) Current iteration of iterative method.

    \param CurrentResVector (In) The current residuals of the iterative process.  

    \param CurrentResNormEst (In) Estimate of the two-norm of the residual.  The value will be
    set to -1.0 if no estimate is available.

    \param SolutionUpdated (In) If this argument is true, then the solution vector that is part 
    of the Epetra_LinearProblem
    object being solved is consistent with the residual. 

    \return AztecOO_StatusType: Unconverged, Converged or Failed.
  */
  virtual AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
				 double CurrentResNormEst,
				 bool SolutionUpdated);
  virtual AztecOO_StatusType GetStatus() const;

  virtual ostream& Print(ostream& stream, int indent = 0) const;
  //@}
protected:

  //@{ \name Internal methods.
  //! Use this for checkStatus when this is an OR type combo. Updates status.
  virtual void OrOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
		    bool SolutionUpdated);

  //! Use this for checkStatus when this is an AND type combo. Updates status.
  virtual void AndOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
		     bool SolutionUpdated);

  //! Check whether or not it is safe to add a to the list of
  //! tests. This is necessary to avoid any infinite recursions.
  bool IsSafe(AztecOO_StatusTest& a);
  //@}

 private:

  //@{ \name Private data members.
  //! Type of test
  ComboType type_;

  //! Vector of generic status tests
  vector<AztecOO_StatusTest*> tests_;

  //! Status
   AztecOO_StatusType status_;
  //@}

};

#endif /* AZTECOO_STATUSTESTCOMBO_H */
