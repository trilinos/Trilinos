
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

#include <map>
#include <string>
using namespace std;
typedef map< string, int >  OptionsMap ;
typedef map< string, double >  ParametersMap ;

#ifndef _SUPERLUDISTOO_H_
#define _SUPERLUDISTOO_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
//! SuperludistOO:  An object-oriented wrapper for Superludist.
/*!  SuperludistOO will solve a linear systems of equations: \f$ AX=B
  \f$, using Epetra objects and the Superludist solver library, where
  \f$A\f$ is an Epetra_CrsMatrix and   \f$X\f$ and \f$B\f$ are 
  Epetra_MultiVector objects.

  SuperLUdist execution can be tuned through a variety of parameters.  
  SuperludistOO.h allows control of these parameters through the
  following named parameters.  SuperludistOO() will ignore parameters 
  with names that it does not recognize, allowing the same list of 
  parameters to be used in calls to other direct solvers.  Where possible, the 
  parameters are common to all direct solvers (although some may ignore 
  them).  However, some parameters, in particular tuning parameters, are
  unique to each solver.

  SuperludistOO consists of five steps.  The first three listed below are referred
  to as the pre-factorization transformations.
  1)  Equilibration - to reduce the condition number of the problem
  2)  Row Permutation - to make the diagonal values larger
  3)  Col permutation - to help maintain sparsity 
  4)  Factorization - Compute the L and U factors
  5)  Solve - Perform the forward and back solves

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
				

     "RowPermutationType"   - int      - enum DsolColumnPermutationTypeOption{
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
class SuperludistOO {
    
  public:
  /*! Creates a SuperludistOO instance, passing in already-defined objects for the matrix
      (as an Epetra_RowMatrix),
      left-hand-side and right-hand-side. 

      Note:  The matrix A must be an Epetra_CrsMatrix.
  */

  SuperludistOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);

  //! SuperludistOO Constructor.
  /*! Creates a SuperludistOO instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. The
      Epetra_LinearProblem class is the preferred method for passing
      in the linear problem to SuperludistOO because this class
      provides scaling capabilities and self-consistency checks that
      are not available when using other constructors.

      Note: The operator in an Epetra_LinearProblem must be an
      Epetra_CrsMatrix.

  */
  SuperludistOO(const Epetra_LinearProblem& LinearProblem);

  //! SuperludistOO Default constructor.
  SuperludistOO();

  //! SuperludistOO Copy Constructor.
  /*! Makes copy of an existing SuperludistOO instance.
  */
  SuperludistOO(const SuperludistOO& Solver);

  //! SuperludistOO Destructor.
  /*! Completely deletes a SuperludistOO object.  
  */
  virtual ~SuperludistOO(void);
  //@}

  //@{ \name Post-construction setup methods.

  //! SuperludistOO Epetra_LinearProblem Set
  /*! Associates an already defined Epetra_LinearProblem as the
      problem that will be solved during iterations.  This method
      allows the user to change which problem is being solved by an
      existing SuperludistOO object.
   */
  int SetProblem(const Epetra_LinearProblem& prob);

  //! SuperludistOO User Matrix Set
  /*! Associates an already defined Epetra_RowMatrix as the matrix
      that will be used by SuperludistOO as the linear operator when
      solving the linear system.  Only Epetra_CrsMatrix objects can be
      passed in through this method.
   */
  int SetUserMatrix(Epetra_RowMatrix * UserMatrix);

  //! SuperludistOO LHS Set
  /*! Associates an already defined Epetra_MultiVector (or
      Epetra_Vector) as the location where the solution will be return
      and, optionally, the initial guess.
   */
  int SetLHS(Epetra_MultiVector * X);

  //! SuperludistOO RHS Set
  /*! Associates an already defined Epetra_MultiVector (or
      Epetra_Vector) as the right-hand-side of the linear system.  */

  int SetRHS(Epetra_MultiVector * B);

  //! SuperludistOO Options Set
  /*! Adds these options to the set of options used to control SuperLUdist
   */
  int SetOptions( OptionsMap Options ) ; 

  //! SuperludistOO Parameters Set
  /*! Adds these parameters to the set of parameters used to control SuperLUdist
   */
  int SetParameters( ParametersMap Parameters ) ; 

  //! SuperludistOO Label Matrix for Superludist
  /*! This is used to label individual matrices within Superludist. This might
      be useful if several Superludist invocations, corresponding
      to different matrices, are involved.  */
  int  SetMatrixName(int label);
  //@}

  //@{ \name Check/Attribute Access Methods.
    
  //! Prints a summary of solver parameters, performs simple sanity checks.
  int CheckInput() const { };

  //! Get a copy of the Parameters
  ParametersMap GetParameters( ) const { return(Parameters_); } ; 

  //! Get a copy of the Options
  OptionsMap GetOptions( ) const { return(Options_); } ; 

  //! Get a pointer to the user matrix A.
  Epetra_RowMatrix * GetUserMatrix() const {return(UserMatrix_);};
  //! Get a pointer to the left-hand-side X.
  Epetra_MultiVector * GetLHS() const {return(X_);};
  //! Get a pointer to the right-hand-side B.
  Epetra_MultiVector * GetRHS() const {return(B_);};
  //@}

  //@{ \name Unsupported functions
  /*! The following member functions are not included in the
    interface yet, in part because I am not sure what the interface
    to them should be.

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

    //! Returns the total number of iterations performed on this problem.
    int NumIters() const {};
    
    //! Returns the true unscaled residual for this problem.
    double TrueResidual() const {};
    
    //! Returns the true scaled residual for this problem.
    double ScaledResidual() const {};
    
    //! SuperludistOO status extraction function.
    /*! Extract Superludist status array into user-provided array.  The array must be of 
        length AZ_STATUS_SIZE as defined in the az_superludist.h header file.

	KENTODO - Figure out how this translates to SuperludistOO.

	NOT supported in release 0.1
     */
    int GetAllSuperludistStatus(double * status) {} ;

  //! Returns the condition number estimate for the current, if one exists, returns -1.0 if no estimate
    /*! Not supported in release 0.1
     */
  double Condest() const {};

  //@}


  //
  //  
  //

#if 0
  bool GetPermc( ) const { return Permc_ ;} ;

  void SetPermc( bool permc ) { Permc_ = permc ;} ; 

#endif
  int SetSuperludistDefaults();

  bool GetTrans( ) const { return Transpose_ ;} ;

  void SetTrans( bool trans ) { Transpose_ = trans ;} ; 

  int SymbolicFactor() ;

  int NumericFactor() ;

  int Solve() ;

 protected:

  Epetra_Operator * UserOperator_;
  Epetra_RowMatrix * UserMatrix_;
  //  Epetra_Operator * PrecOperator_;
  //  Epetra_RowMatrix * PrecMatrix_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  bool Transpose_ ;
  int Permc_ ; 

  int x_LDA_;
  double *x_;
  int b_LDA_;
  double *b_;

  OptionsMap Options_;
  ParametersMap Parameters_;

  bool inConstructor_;
};


#endif /* _SUPERLUDISTOO_H_ */

