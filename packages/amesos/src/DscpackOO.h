// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef _DSCPACKOO_H_
#define _DSCPACKOO_H_

#include "Amesos_ConfigDefs.h"
#include "Epetra_MpiComm.h"

#ifndef DSC_LBLAS1
#define DBL_R_NUM
extern "C" {
#include "dscmain.h"
}
#endif

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
//  #include "Epetra_LinearProblemRedistor.h"
#include "Epetra_Object.h"
//! DscpackOO:  An object-oriented wrapper for Dscpack.
/*!  DscpackOO will solve a linear systems of equations: \f$ AX=B
  \f$, using Epetra objects and the Dscpack solver library, where
  \f$A\f$ is an Epetra_RowMatrix and \f$X\f$ and \f$B\f$ are 
  Epetra_MultiVector objects.

  Dscpack execution can be tuned through a variety of parameters.
  DscpackOO.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
*/
class DscpackOO {
    
  public:
  //@{ \name Constructor methods
  //! DscpackOO Constructor.
  /*! Creates a DscpackOO instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. The
      Epetra_LinearProblem class is the preferred method for passing
      in the linear problem to DscpackOO because this class
      provides scaling capabilities and self-consistency checks that
      are not available when using other constructors.

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  DscpackOO(const Epetra_LinearProblem& LinearProblem);

  //! DscpackOO Destructor.
  /*! Completely deletes a DscpackOO object.  
  */
  virtual ~DscpackOO(void);
  //@}

  //@{ \name Post-construction setup methods.

  //!  Setting the transpose flag to true causes Solve() to compute A^t x = b 
  void SetTrans( bool trans ) { } ; 


  //@}
  //@{ \name Check/Attribute Access Methods.
    
  //!  Return the transpose flag 
  bool GetTrans( ) const { return false; } ;

  //! Prints a summary of solver parameters, performs simple sanity checks.
  /*!
    Not supported in release 0.1;
   */
  int CheckInput() const ;

  //@}

  //@{ \name Setting and Clearing Compact Representations of the pre-factorzation transformtaions - Not yet defined
  /*! Member functions to set and clear the compact representations of the pre-factorizations are not included in the
    interface yet, in part because I am not sure what the interface
    to them should be.  Here is a description of these functions, whose interface we will define later.

    The next two member functions have no impact on the factorization unless 
    "EquilibrationReuse" is set to DSOL_USE_STORED_EQUILIBRATION
    SetLeftEquilibrationVector - Sets the left equilibration vector for use
    on the next matrix factorization as specified by "EquilibrationReuse"
    ClearLeftEquilibrationVector - Forces the next matrix factorization to 
    compute an equilibration as specified by "EquilibrationType"
    GetLeftEquilibrationVector - 

    The following member functions are analogous to the ones for Left Equilibration.  
    SetRightEquilibrationVector - 
    ClearRightEquilibrationVector -
    GetRightEquilibrationVector - 

    The next two member functions have no impact on the factorization unless 
    "RowPermutationReuse" is set to DSOL_USE_STORED_ROW_PERM
    SetRowPermutationVector - Sets the left row permutation vector for use
    on the next matrix factorization as specified by "Row PermutationReuse"
    ClearRowPermutationVector - Forces the next matrix factorization to 
    compute a row permutation as specified by "RowPermutationType"
    GetRowPermutationVector - 

    The following member functions are analogous to the ones for Row Permutation.  
    SetColumnPermutationVector - 
    ClearColumnPermutationVector -
    GetColumnPermutationVector - 

   */
  //@}

  //@{ \name Post-solve access functions

  //! Returns the condition number estimate for the current problem, if one exists, returns -1.0 if no estimate
    /*! Not supported in release 0.1
     */
  double Condest() const;

  //@}


  //@{ \name Solve method
  //!  All computation is performed during the call to Solve() 
  /*!  Factor controls whether or not the matrix should be factored prior to the solve.
       Default is true.
   */
  int Solve(bool Factor) ;

  //@}
 protected:

  //
  //  These are not used in release 0.1
  //
  const Epetra_LinearProblem * Problem_;
  
  DSC_Solver	MyDSCObject;
  MPI_Comm MPIC ; 

  bool Factored_;
  bool FirstCallToSolve_;
  bool A_and_LU_built ;            // Tells us whether to free them 
  int *GlobalStructNewColNum ; 
  int *GlobalStructNewNum ;  
  int *GlobalStructOwner ; 
  int *LocalStructOldNum ; 

  int MyDscRank ; 
  int DscNumProcs ; 
  int NumLocalCols; 
  //  I don't see why the next three belong here (if they do)
  int NumGlobalCols;
  int NumLocalPPPStructs ; 
#if 0
  int NumLocalNonz ;
#endif

};


#endif /* _DSCPACKOO_H_ */

