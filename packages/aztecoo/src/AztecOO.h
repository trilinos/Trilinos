
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

#ifndef _AZTECOO_H_
#define _AZTECOO_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "az_aztec.h"


//! AztecOO:  An object-oriented wrapper for Aztec.
/*! Currently it accepts a Petra matrix, initial guess and RHS as
  separate arguments, or alternatively, accepts a Epetra_LinearProblem.
  If constructed using a Epetra_LinearProblem, AztecOO will infer some
  solver/preconditioner, etc., options and parameters. Users may override
  these choices and manually choose from among the full set of Aztec options
  using the SetAztecOption() and SetAztecParam() functions.
*/

class AztecOO: public Epetra_Object{
    
  public:
  //!  AztecOO Constructor.
  /*! Creates a AztecOO instance, using separate A,X,B. 
  */
  AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X,
			 Epetra_MultiVector * B);

  //! AztecOO Constructor.
  /*! Creates a AztecOO instance, using a Epetra_LinearProblem.
  */
  AztecOO(const Epetra_LinearProblem& tlp);

  //! AztecOO Default constructor.
  AztecOO();

  //! AztecOO Copy Constructor.
  /*! Makes copy of an existing AztecOO instance.
  */
  AztecOO(const AztecOO& Solver);

  //! AztecOO External Scaling Set
  /*! Associates an already defined Aztec scaling object with this solve.
   */
  int SetPreconditioner(struct AZ_SCALING * Scaling) {Scaling_ = Scaling; return(0);};


  //! AztecOO External Preconditioner Set (object)
  /*! Associates an already defined Aztec preconditioner with this solve.
   */
  int SetPreconditioner(AZ_PRECOND * Prec) {Prec_ = Prec; return(0);};


  //! AztecOO External Preconditioner Set (function and data)
  /*! Associates an external function and data pointer with preconditioner
   */
  int SetPreconditioner(void  (*prec_function)(double *, int *, int *, double *,
					   struct AZ_MATRIX_STRUCT  *, 
					   struct AZ_PREC_STRUCT *),
		    void *prec_data);

  //! Forces explicit construction and retention of preconditioner.
  /*! Aztec typically constructs the precondition on the first call to the solve function.
      However, there are situations where we would like to compute the preconditioner ahead
      of time.  One particular case is when we want to confirm that the preconditioner 
      well-conditioned.  This method allows us to precompute the preconditioner.  It also
      provides a estimate of the condition number of the preconditioner.  If \it condest is
      large, e.g., > 1.0e+14, it is likely the preconditioner will fail.  In this case, using 
      threshold values (available in the incomplete factorizations) can be used to reduce
      the condition number.
  */
  int ConstructPreconditioner(double & condest);

  //! Destroys a preconditioner computed using ConstructPreconditioner().
  /*! The ConstructPreconditioner() method creates a persistent preconditioner.
      In other words the preconditioner will be used by all calls to the Iterate() 
      method.  DestroyPreconditioner() deletes the current preconditioner and restores
      AztecOO to a state where the preconditioner will computed on first use of the
      preconditioner solve.
  */
  int DestroyPreconditioner();

  //! Returns the condition number estimate for the current, if one exists, returns -1.0 if no estimate
  double Condest() const {return(condest_);};


  //! AztecOO Label Matrix for Aztec
  /*! This is used to label individual matrices within Aztec. This might
      be useful if several Aztec invokations are involved corresponding
      to different matrices.
   */
  int  SetMatrixName(int label);


  //! AztecOO iteration functions.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached..
      This one should be suitable for recursive invokations of Aztec.
  */
  int recursiveIterate(int MaxIters, double Tolerance);

  //! AztecOO iteration function.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached..
  */
  int Iterate(int MaxIters, double Tolerance);

  //! AztecOO iteration function.
  /*! Iterates on the specified matrix and vectors until MaxIters or Tolerance
      is reached..
  */
  int Iterate(Epetra_RowMatrix * A,
                      Epetra_MultiVector * X,
                      Epetra_MultiVector * B,
                      int MaxIters, double Tolerance);

  //! AztecOO option setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecOption(AZ_precond, AZ_Jacobi)
   */
  int SetAztecOption(int option, int value)
    {options_[option] = value; return(0);};

  //! AztecOO param setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecParam(AZ_drop, 1.0E-6)
   */
    int SetAztecParam(int param, double value)
    {params_[param] = value; return(0);};

  //! AztecOO option setting function.
  /*! Set all Aztec option values using an existing Aztec options array.
   */
    int SetAllAztecOptions(int * options)
      {for (int i=0; i<AZ_OPTIONS_SIZE; i++) options_[i] = options[i]; return(0);};

  //! AztecOO param setting function.
  /*! Set all Aztec parameter values using an existing Aztec params array. 
   */
    int SetAllAztecParams(double * params)
      {for (int i=0; i<AZ_PARAMS_SIZE; i++) params_[i] = params[i]; return(0);};
 
  //! AztecOO function to restore default options/parameter settings.
  /*! This function is called automatically within AztecOO's constructor,
   but if constructed using a Epetra_LinearProblem object, some options are
   reset based on the ProblemDifficultyLevel associated with the 
   Epetra_LinearProblem. 
  */
    int SetAztecDefaults();

    // Adaptive Solve methods

    //! Force the AdaptiveIterate() method to use default adaptive strategy.
    int SetUseAdaptiveDefaultsTrue(){useAdaptiveDefaults_ = true;return(0);};

  //! Set the parameter that control the AdaptiveIterate() method.
  /*! The AdaptiveIterate() method attempts to solve a given problem using multiple preconditioner
      and iterative method tuning parameters.  There are defaults that are coded into AdaptiveIterate()
      method, but the defaults can be over-ridden by the use of the SetAdaptiveParams() method. Details of 
      condition number management follow:
\verbinclude Managing_conditioning_howto.txt

    \param NumTrials In
           The number of Athresh and Rthresh pairs that should be tried when attempting to stabilize 
	   the preconditioner.
    \param athresholds In
           The list of absolute threshold values that should be tried when attempting to stabilize 
	   the preconditioner.
    \param rthresholds In
           The list of relative threshold values that should be tried when attempting to stabilize 
	   the preconditioner.
    \param condestThreshold In
           If the condition number estimate of the preconditioner is above this number, no attempt will be
	   made to try iterations.  Instead a new preconditioner will be computed using the next threshold 
	   pair.
    \param maxFill In
           In addition to managing the condest, the AdaptiveIterate() method will also try to increase
	   the preconditioner fill if it is determined that this might help.  maxFill specifies the
	   maximum fill allowed.
    \param maxKspace In
           In addition to managing the condest, the AdaptiveIterate() method will also try to increase
	   the Krylov subspace size if GMRES is being used and it is determined that this might help.
	   maxKspace specifies the maximum Krylov subspace allowed.
	   
  */
    int SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds,
			  double condestThreshold, double maxFill, int maxKspace);

    //! Attempts to solve the given linear problem using an adaptive strategy.
    int AdaptiveIterate(int MaxIters, double Tolerance);

    // Post-solve access functions

    //! Returns the total number of iterations performed on this problem.
    int NumIters() const {return((int) status_[AZ_its]);};

  //! Returns the true unscaled residual for this problem.
  double TrueResidual() const {return(status_[AZ_r]);};

  //! Returns the true scaled residual for this problem.
  double ScaledResidual() const {return(status_[AZ_scaled_r]);};

  //! AztecOO status extraction function.
  /*! Extract Aztec status array into user-provided array. 
   */
  int GetAllAztecStatus(double * status)
    {for (int i=0; i<AZ_STATUS_SIZE; i++) status[i] = status_[i]; return(0);};

  //! Prints a summary of solver parameters, performs simple sanity checks.
  int CheckInput() const {
    return(AZ_check_input(Amat_->data_org, options_, params_, proc_config_));};

  //! AztecOO Destructor.
  /*! Completely deletes a AztecOO object.  
  */
  virtual ~AztecOO(void);

 protected:

  int AllocAzArrays();
  void DeleteAzArrays();
  int SetAztecVariables();
  int SetProblemOptions(ProblemDifficultyLevel PDL,
                            bool ProblemSymmetric);
  void DeleteMemory();

  Epetra_RowMatrix * A_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  int N_update_;
  int N_local_;
  int * update_;
  int x_LDA_;
  double *x_;
  int b_LDA_;
  double *b_;
  int * proc_config_;
  int * options_;
  double * params_;
  double * status_;
  AZ_MATRIX *Amat_;
  AZ_PRECOND *Prec_;
  struct AZ_SCALING * Scaling_;

  bool azVarsAllocated_;
  double condest_;
  bool useAdaptiveDefaults_;
  int NumTrials_;
  double maxFill_;
  int maxKspace_;
  double * athresholds_;
  double * rthresholds_;
  double condestThreshold_;
};

// External prototypes
extern "C" void Epetra_Aztec_matvec(double *x, double *y, AZ_MATRIX *Amat, int proc_config[]);
extern "C" int Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[], 
				  AZ_MATRIX *Amat, int N_requested_rows,
				  int requested_rows[], int allocated_space);
extern "C" int Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat);

#ifdef AZTEC_OO_WITH_ML

extern "C" int Epetra_ML_matvec(void *data, int in, double *p, int out, double *ap);

extern "C" int Epetra_ML_getrow(void *data, int N_requested_rows, int requested_rows[],
                    int allocated_space, int columns[], double values[],
                    int row_lengths[]);

extern "C" int Epetra_ML_comm_wrapper(double vec[], void *data);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //! Trilinos_LinearProblem Convert Petra Matrix To ML matrix.
  /*! Take a Petra matrix and convert it to an ML matrix and stick it
      in a multigrid hierarchy at 'level'. This function should probably
      not sit here but be inside some ML_Petra wrapper object.
   */
  int  PetraMatrix2MLMatrix(ML *ml_handle, int level,
			    Epetra_RowMatrix * Amat);

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // AZTEC_OO_WITH_ML

#endif /* _AZTECOO_H_ */

