
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_ConfigDefs.h"

#include "superlu_ddefs.h"
#include "supermatrix.h"
//  SuperLU defines Reduce to be a macro in util.h, this conflicts with Reduce() in Epetra_MultiVector.h
#undef Reduce

#ifndef _SUPERLUDIST2_OO_H_
#define _SUPERLUDIST2_OO_H_

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_LinearProblemRedistor.h"
#include "Epetra_Object.h"
//! Superludist2_OO:  An object-oriented wrapper for Superludist.
/*!  Superludist2_OO will solve a linear systems of equations: \f$ AX=B
  \f$, using Epetra objects and the Superludist solver library, where
  \f$A\f$ is an Epetra_RowMatrix and \f$X\f$ and \f$B\f$ are 
  Epetra_MultiVector objects.

  SuperLUdist execution can be tuned through a variety of parameters.
  Superludist2_OO.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.

  Superludist2_OO consists of five steps.  The first three listed below are referred
  to as the pre-factorization transformations.
  1)  Equilibration - to reduce the condition number of the problem
  2)  Row Permutation - to make the diagonal values larger
  3)  Col permutation - to help maintain sparsity 
  4)  Factorization - Compute the L and U factors
  5)  Solve - Perform the forward and back solves

  In other solvers, the steps are: symbolic factorization (Step 3 above), 
  numeric factorization (Step 4 above) and Solve (Step 5 above).  Step 2, row 
  permutation can be considered static pivoting and equilibration is akin
  to left and right scaling.  

  NO MECHANISM EXISTS TODAY TO SET OPTIONS OR PARAMETERS.  

  The following parameters are generic, i.e. potentially applicable 
  to all direct solvers.  The ints are passed as options.  The doubles 
  are passed as parameters.

    "EquilibrationType"       - int    - enum DsolEquilibrationOption{ 
                                         DSOL_DO_NOT_EQUILIBRATE, 
					 DSOL_ROW_EQUILIBRATION,
					 DSOL_COL_EQUILIBRATION,
					 DSOL_ROW_COL_EQUILIBRATION }
					 
    "EquilibrationReuse"      - int    - enum DsolEquilibrationReuseOption{ 
                                         DSOL_EQUILIBRATE_NEW_MATRICES, 
					 DSOL_USE_STORED_EQUILIBRATION }
					 If DSOL_USE_STORED_EQUILIBRATION 
					 is set, new equilibration 
					 arrays are only computed if 
					 none have been stored.  There are
					 two ways that equilibration constants
					 can be stored.  1)  Each time that 
					 the matrices are equilibrated, 
					 the equilibrations are stored.  
					 2)  The user can store equilibration
					 constants by calling SetLeftEquilibration() 
					 and SetRightEquilibration().
					 
    "ColumnPermutationType"   - int    - enum DsolColumnPermutationTypeOption{
                                         DSOL_NO_COLUMN_PERMUTATION,
                                         DSOL_MMD_AT_times_A, 
                                         DSOL_MMD_AT_plus_A, 
					 DSOL_COLAMD }


     "ColumnPermutationReuse" -  int   - enum DsolColumnPermutationReuseOption{ 
                                         DSOL_COL_PERM_NEW_MATRICES, 
					 DSOL_USE_STORED_COL_PERM }
					 If DSOL_USE_STORED_COL_PERM 
					 is set, a new column permutation
					 is only computed if 
					 none has been stored.  There are
					 two ways that a column permutation
					 can be stored.  1)  Each time that 
					 a column permutation is computed
					 it is stored.  
					 2)  The user can store a column permutation
					 by calling SetColumnPermutation().
				

     "RowPermutationType"   - int      - enum DsolRowPermutationTypeOption{
                                         DSOL_NO_ROW_PERMUTATION,
                                         DSOL_DUFF_KOSTER }

     "RowPermutationReuse" -  int      - enum DsolRowPermutationReuseOption{ 
                                         DSOL_ROW_PERM_NEW_MATRICES, 
					 DSOL_USE_STORED_ROW_PERM }
					 If DSOL_USE_STORED_ROW_PERM 
					 is set, a new row permutation
					 is only computed if 
					 none has been stored.  There are
					 two ways that a row permutation
					 can be stored.  1)  Each time that 
					 a row permutation is computed
					 it is stored.  
					 2)  The user can store a row permutation
					 by calling SetRowPermutation().
				
    "FactorType"            - int -      enum DsolFactorTypeOption{
                                         DSOL_DO_FACTOR,
                                         DSOL_DO_NOT_FACTOR } 
					 "FactorType" is not fully supported in release 0.1

    "FactorTypeReuse"       - int -      enum DsolFactorTypeOption{
                                         DSOL_REUSE_FACTOR,
                                         DSOL_DO_NOT_REUSE_FACTOR } 
					 "FactorTypeReuse" is not fully supported in release 0.1

    "BLAS block size"    - int -         BLAS block size 

    
*/
class Superludist2_OO {
    
  public:
  //@{ \name Constructor methods
  //! Superludist2_OO Constructor.
  /*! Creates a Superludist2_OO instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. The
      Epetra_LinearProblem class is the preferred method for passing
      in the linear problem to Superludist2_OO because this class
      provides scaling capabilities and self-consistency checks that
      are not available when using other constructors.

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Superludist2_OO(const Epetra_LinearProblem& LinearProblem);

  //! Superludist2_OO Destructor.
  /*! Completely deletes a Superludist2_OO object.  
  */
  virtual ~Superludist2_OO(void);
  //@}

  //@{ \name Post-construction setup methods.

  //!  Setting the transpose flag to true causes Solve() to compute A^t x = b 
  /*!
    Transpose not supported in release 0.1;
   */
  void SetTrans( bool trans ) { 
    assert( trans == false) ; 
    Transpose_ = trans ;} ; 


  //@}
  //@{ \name Check/Attribute Access Methods.
    
  //!  Return the transpose flag 
  bool GetTrans( ) const { return Transpose_ ;} ;

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
  Epetra_LinearProblemRedistor *redistor;
  Epetra_LinearProblem *redistProblem;
  //
  //  Here are the values returned by ExtractHbData
  //
  int M,N,nz;
  int *ptr, *ind;
  double *val, *rhs, *lhs;
  int Nrhs, ldrhs, ldlhs;
  

  bool Transpose_ ;
  bool Factored_;
  bool FirstCallToSolve_;
  //
  //  Here are the SuperLU data structures for A, L and U:
  //
  SOLVEstruct_t SOLVEstruct; 
  int numprocs;
  int nprow;
  int npcol;
  gridinfo_t grid;                 // SuperLU's grid information
  superlu_options_t options;
  SuperMatrix A;
  ScalePermstruct_t ScalePermstruct;
  SuperLUStat_t stat;
  LUstruct_t LUstruct;
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;
  bool A_and_LU_built ;            // Tells us whether to free them 


  //  This is needed by the old Superludist2_OO.cpp
  int numrows ; 

};


#endif /* _SUPERLUDIST2_OO_H_ */
