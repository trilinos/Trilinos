#ifndef _Aztec_OO_h_
#define _Aztec_OO_h_

//! Aztec_OO:  An object-oriented wrapper for Aztec.
/*! Currently it accepts a Petra matrix, initial guess and RHS as
  separate arguments, or alternatively, accepts a Petra_RDP_LinearProblem.
  If constructed using a Petra_RDP_LinearProblem, Aztec_OO will infer some
  solver/preconditioner, etc., options and parameters. Users may override
  these choices and manually choose from among the full set of Aztec options
  using the SetAztecOption() and SetAztecParam() functions.
*/

#include "Petra_Object.h"
#include "Petra_Comm.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_RowMatrix.h"
#include "Petra_RDP_LinearProblem.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifdef PETRA_MPI
#include "mpi.h"
#define AZ_MPI
#define AZTEC_MPI
#endif
#include "az_aztec.h"
#endif

class Aztec_OO: public Petra_Object{
    
  public:
  //!  Aztec_OO Constructor.
  /*! Creates a Aztec_OO instance, using separate A,X,B. 
  */
  Aztec_OO(Petra_RDP_RowMatrix * A, Petra_RDP_MultiVector * X,
			 Petra_RDP_MultiVector * B);

  //! Aztec_OO Constructor.
  /*! Creates a Aztec_OO instance, using a Petra_RDP_LinearProblem.
  */
  Aztec_OO(const Petra_RDP_LinearProblem& tlp);

  //! Aztec_OO Default constructor.
  Aztec_OO();

  //! Aztec_OO Copy Constructor.
  /*! Makes copy of an existing Aztec_OO instance.
  */
  Aztec_OO(const Aztec_OO& Solver);

  //! Aztec_OO External Scaling Set
  /*! Associates an already defined Aztec scaling object with this solve.
   */
  int SetPreconditioner(struct AZ_SCALING * Scaling) {Scaling_ = Scaling; return(0);};


  //! Aztec_OO External Preconditioner Set (object)
  /*! Associates an already defined Aztec preconditioner with this solve.
   */
  int SetPreconditioner(AZ_PRECOND * Prec) {Prec_ = Prec; return(0);};


  //! Aztec_OO External Preconditioner Set (function and data)
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
      Aztec_OO to a state where the preconditioner will computed on first use of the
      preconditioner solve.
  */
  int DestroyPreconditioner();

  //! Returns the condition number estimate for the current, if one exists, returns -1.0 if no estimate
  double Condest() const {return(condest_);};


  //! Aztec_OO Label Matrix for Aztec
  /*! This is used to label individual matrices within Aztec. This might
      be useful if several Aztec invokations are involved corresponding
      to different matrices.
   */
  int  SetMatrixName(int label);


  //! Aztec_OO iteration functions.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached..
      This one should be suitable for recursive invokations of Aztec.
  */
  int recursiveIterate(int MaxIters, double Tolerance);

  //! Aztec_OO iteration function.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached..
  */
  int Iterate(int MaxIters, double Tolerance);

  //! Aztec_OO iteration function.
  /*! Iterates on the specified matrix and vectors until MaxIters or Tolerance
      is reached..
  */
  int Iterate(Petra_RDP_RowMatrix * A,
                      Petra_RDP_MultiVector * X,
                      Petra_RDP_MultiVector * B,
                      int MaxIters, double Tolerance);

  //! Aztec_OO option setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecOption(AZ_precond, AZ_Jacobi)
   */
  int SetAztecOption(int option, int value)
    {options_[option] = value; return(0);};

  //! Aztec_OO param setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecParam(AZ_drop, 1.0E-6)
   */
    int SetAztecParam(int param, double value)
    {params_[param] = value; return(0);};

  //! Aztec_OO option setting function.
  /*! Set all Aztec option values using an existing Aztec options array.
   */
    int SetAllAztecOptions(int * options)
      {for (int i=0; i<AZ_OPTIONS_SIZE; i++) options_[i] = options[i]; return(0);};

  //! Aztec_OO param setting function.
  /*! Set all Aztec parameter values using an existing Aztec params array. 
   */
    int SetAllAztecParams(double * params)
      {for (int i=0; i<AZ_PARAMS_SIZE; i++) params_[i] = params[i]; return(0);};
 
  //! Aztec_OO function to restore default options/parameter settings.
  /*! This function is called automatically within Aztec_OO's constructor,
   but if constructed using a Petra_RDP_LinearProblem object, some options are
   reset based on the ProblemDifficultyLevel associated with the 
   Petra_RDP_LinearProblem. 
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

  //! Aztec_OO status extraction function.
  /*! Extract Aztec status array into user-provided array. 
   */
  int GetAllAztecStatus(double * status)
    {for (int i=0; i<AZ_STATUS_SIZE; i++) status[i] = status_[i]; return(0);};

  //! Prints a summary of solver parameters, performs simple sanity checks.
  int CheckInput() const {
    return(AZ_check_input(Amat_->data_org, options_, params_, proc_config_));};

  //! Aztec_OO Destructor.
  /*! Completely deletes a Aztec_OO object.  
  */
  virtual ~Aztec_OO(void);

 protected:

  int AllocAzArrays();
  void DeleteAzArrays();
  int SetAztecVariables();
  int SetProblemOptions(ProblemDifficultyLevel PDL,
                            bool ProblemSymmetric);
  void DeleteMemory();

  Petra_RDP_RowMatrix * A_;
  Petra_RDP_MultiVector * X_;
  Petra_RDP_MultiVector * B_;

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
extern "C" void Petra_Aztec_matvec(double *x, double *y, AZ_MATRIX *Amat, int proc_config[]);
extern "C" int Petra_Aztec_getrow(int columns[], double values[], int row_lengths[], 
				  AZ_MATRIX *Amat, int N_requested_rows,
				  int requested_rows[], int allocated_space);
extern "C" int Petra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat);

#ifdef AZTEC_OO_WITH_ML

extern "C" int Petra_ML_matvec(void *data, int in, double *p, int out, double *ap);

extern "C" int Petra_ML_getrow(void *data, int N_requested_rows, int requested_rows[],
                    int allocated_space, int columns[], double values[],
                    int row_lengths[]);

extern "C" int Petra_ML_comm_wrapper(double vec[], void *data);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //! Trilinos_LinearProblem Convert Petra Matrix To ML matrix.
  /*! Take a Petra matrix and convert it to an ML matrix and stick it
      in a multigrid hierarchy at 'level'. This function should probably
      not sit here but be inside some ML_Petra wrapper object.
   */
  int  PetraMatrix2MLMatrix(ML *ml_handle, int level,
			    Petra_RDP_RowMatrix * Amat);

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // AZTEC_OO_WITH_ML

#endif /* _Aztec_OO_h_ */

