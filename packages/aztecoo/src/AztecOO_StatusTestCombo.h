
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

//! AztecOO_StatusTestCombo: A  class for extending the status testing capabilities of AztecOO via logical combinations.

/*! AztecOO_StatusTestCombo is an interface that can be implemented to extend the convergence testing
   capabilities of AztecOO.  This class supports composite tests.  In this situation,
   two or more existing AztecOO_StatusTestCombo objects test1 and test2 can be used to create a new test.
   For all combinations, if any tests returns Failed or returns not-a-number (NaN) status, then the combination test 
   returns Failed.
   There are three possible combinations:
   <ol>
   <li> OR combination:
   If an OR combination is selected, the status returns Converged if any one of the subtest returns
   as Converged.  
   <li> AND combination:
   If an AND combination is selected, the status returns Converged only when all subtests return as Converged.
   <li> SEQ combination:
   SEQ is a form of AND that will perform subtests in sequence.  If the first test returns Unconverged, Failed or NaN,
   no other subtests are done, and the status is returned as Unconverged if the first test was Unconverged, or as
   Failed if the first test was Failed or NaN.  If the first test returns Converged, the second test is checked in 
   the same fashion as the first.  If the second test is Converged, the third one is tested, and so on.

   The purpose of the SEQ combination is to allow the addition of expensive but more rigorous convergence tests.  For
   example, we could define a test that used the implicit residual vector (the one produced by the iterative method)
   as the first subtest and define a second test using the explicitly computed residual vector.  Explicitly computing
   the residual requires a matrix multiplication with the original matrix operator, an expensive operation.  By using
   the SEQ combination, we can avoid the matrix multiplication associated with the explicit residual calculation
   until the implicit residual is small.
   </ol>
   

*/

class AztecOO_StatusTestCombo: public AztecOO_StatusTest {

 public:

  //@{ \name Enums.
  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests, or a sequential AND (SEQ).
  */
  enum ComboType {AND,  /*!< Require all subtests to be satisfied. */
		  OR,   /*!< Require one or the other subtests to be satisfied. */
		  SEQ   /*!< Requires all subtests to be satisfied, but stops check after the first failed 
			  or unconverged status. */
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
  AztecOO_StatusTestCombo& AddStatusTest(AztecOO_StatusTest& a);

  //! Destructor
  virtual ~AztecOO_StatusTestCombo() {};
  //@}

  //@{ \name Methods that implement the AztecOO_StatusTest interface.
  //! Indicates if residual vector is required by this convergence test.
  /*! If this method returns true, then one or more of the AztecOO_StatusTest objects that make up this combined
    test requires the Residual Vector to perform its test.
  */
  bool ResidualVectorRequired() const;

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
  AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
				 double CurrentResNormEst,
				 bool SolutionUpdated);
  AztecOO_StatusType GetStatus() const {return(status_);};

  ostream& Print(ostream& stream, int indent = 0) const;

  //@}

  //@{ \name Methods to access data members.

  //! Returns the maximum number of iterations set in the constructor.
  ComboType GetComboType() const {return(type_);};

  //@}
protected:

  //@{ \name Internal methods.
  //! Use this for checkStatus when this is an OR type combo. Updates status.
  void OrOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
		    bool SolutionUpdated);

  //! Use this for checkStatus when this is an AND type combo. Updates status.
  void AndOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
		     bool SolutionUpdated);

  //! Use this for checkStatus when this is a sequential AND type combo. Updates status.
  void SeqOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
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
  std::vector<AztecOO_StatusTest*> tests_;

  //! Status
   AztecOO_StatusType status_;
  //@}

};

#endif /* AZTECOO_STATUSTESTCOMBO_H */
