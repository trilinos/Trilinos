/*!
 * \file ml_epetra_preconditioner.h
 *
 * \class MultiLevelPreconditioner
 *
 * \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 * ML offers two preconditioners suitable for the solution of 
 * Epetra_LinearProblem objects. This file define one the two, called
 * MultiLevelOperator (in the ML_Epetra namespace). This preconditioner is
 * simple wrapper of the ML_Solve() function, so that ML can be applied to
 * Epetra_MultiVector's. 
 *
 * When you should use MultiLevelOperator:
 * - when your code already defines the required ML objects, with the optimal
 *   choice of parameters, and you want to use ML for Epetra_LinearProblem or
 *   AztecOO problems;
 *
 * When you should use MultiLevelPreconditioner:  
 * - when you have an Epetra_RowMatrix, and you don't want to code the
 *   conversion to ML_Operator, the creation of the hierarchy and the
 *   aggregates, and/or you want to experiment various combinations of the
 *   parameters, simply changing some parameters in a Teuchos::ParameterList.
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

#ifndef _ML_EPETRA_PRECONDITIONER_H_
#define _ML_EPETRA_PRECONDITIONER_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_Comm;
class Epetra_CrsMatrix;
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

#define ML_MEM_SIZE      20
#define ML_MEM_INITIAL    0
#define ML_MEM_FINAL      1
#define ML_MEM_SMOOTHER   2
#define ML_MEM_COARSE     3
#define ML_MEM_HIERARCHY  4
#define ML_MEM_PREC_FIRST 5
#define ML_MEM_PREC_OTHER 6
#define ML_MEM_TOT1       7
#define ML_MEM_TOT2       8
#define ML_MEM_INITIAL_USED    10
#define ML_MEM_FINAL_USED      11
#define ML_MEM_SMOOTHER_USED   12
#define ML_MEM_COARSE_USED     13
#define ML_MEM_HIERARCHY_USED  14
#define ML_MEM_PREC_FIRST_USED 15
#define ML_MEM_PREC_OTHER_USED 16
#define ML_MEM_TOT1_USED       17
#define ML_MEM_TOT2_USED       18

#ifdef HAVE_ML_TRIUTILS
#include "Trilinos_Util_CommandLineParser.h"
#endif

#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

namespace ML_Epetra
{

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  /*! This function, defined in the namespace ML_Epetra, can be used to set
    default values in a user's defined Teuchos::ParameterList.
    \param ProblemType (In) : a string, whose possible values are:
       - "DD" : defaults for 2-level domain decomposition preconditioners based
       on aggregation;
       - "DD-ML" : 3-level domain decomposition preconditioners, with coarser
       spaces defined by aggregation;
       - "SA" : classical smoothed aggregation preconditioners;
       - "maxwell" : default values for Maxwell.
    \param List (Out) : list which will populated by the default parameters
    \param options (In) : integer array, of size \c AZ_OPTIONS_SIZE, that will be
    populated with suitable values. A pointer to \c options will be stick into
    the parameters list. Note that this array is still required to apply the
    preconditioner! Do not delete options, nor let it go out of scope. The default value is 
    0, meaning that \c SetDefaults() will allocate the array. It is
    responsibility of the user to free this memory.
    \param params (Out) : double array, of size \c AZ_PARAMS_SIZE. See comments
    for \c options.    
    \param Prefix (In) : a string value, defaulted to "". All parameters will have
    \c Prefix as prefix. For example, the maximum number of level is defined as
    \c "aggregation: max levels". If \c Prefix == "vel prob: ", than the maximum
    number of levels will be inserted as \c "vel prob: aggregation: max levels".
    An ML_Epetra::MultiLevelPreconditioner can be created with a specified
    prefix. This is useful when more than one preconditioner must be created, and
    the user wants to put all the parameters in the same parameters list.
   */
  int SetDefaults(string ProblemType, Teuchos::ParameterList & List,
		  int * options = 0, double * params = 0, const string Prefix = "");
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsDD(Teuchos::ParameterList & List, const string Prefix = "",
		    int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners, using LU on each subdomain
  int SetDefaultsDD_LU(Teuchos::ParameterList & List, const string Prefix = "",
		       int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.  
  int SetDefaultsDD_3Levels(Teuchos::ParameterList & List, const string Prefix = "",
			    int * options = 0, double * params = 0);
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners with LU.
  int SetDefaultsDD_3Levels_LU(Teuchos::ParameterList & List, const string Prefix = "",
			       int * options = 0, double * params = 0);

  //! Sets default parameters for Maxwell's equations.
  int SetDefaultsMaxwell(Teuchos::ParameterList & List, const string Prefix = "",
			 int * options = 0, double * params = 0);
  
  //! Sets classical smoothed aggregation.
  int SetDefaultsSA(Teuchos::ParameterList & List, const string Prefix = "",
		    int * options = 0, double * params = 0);

#ifdef HAVE_ML_TRIUTILS
  //! Add values set in Command Line Parser to the internally stored parameter list object.
  int Set(Teuchos::ParameterList & List, Trilinos_Util::CommandLineParser & CLP);  
#endif

//! MultiLevelPreconditioner: An implementation of the Epetra_RowMatrix class.
/*! MultiLevelPreconditioner class implements Epetra_RowMatrix using a
    an Epetra_RowMatrix, and possibly a Teuchos parameters list, that
    specifies how to construct the preconditioner.
    The resulting preconditioner is completely black-box. The user needs
    to prive the linear system matrix, and specify in the parameters list 
    the required options.

    The code accepts any Epetra_RowMatrix-derived class. Some code can take
    advantage if Epetra_RowMatrix is an Epetra_CrsMatrix, or and 
    Epetra_VbrMatrix.

    This file requires ML to be configured with the following options:
    - \c --enable-epetra
    - \c --enable-teuchos
    
    The following option is suggested:
    - \c --enable-amesos

    Some part of this class needs the following options:
    - \c --enable-aztecoo
    - \c --enable-anasazi
    
    It is important to note that ML is more restrictive than Epetra for
    the definition of maps. It is required that RowMatrixRowMap() is equal 
    to OperatorRangeMap(). This is because ML needs to perform matrix-vector
    product, as well as getrow() functions, on the same data distribution.
    
    Also, for square matrices, OperatorDomainMap() 
    must be as OperatorRangeMap(). 

    Several examples are provided in the \c examples subdirectories:
    - ml_example_epetra_preconditioner.cpp is an introductory 
      example;
    - ml_example_epetra_preconditioner_2level.cpp shows how to
      define a 2-level domain decomposition preconditioner using 
      this class;
    - ml_example_epetra_preconditioner_viz.cpp details how to
      visualize the aggregates;
    - ml_example_epetra_preconditioner_vbr.cpp is an example for
      VBR matrices;
    - ml_example_epetra_preconditioner_Maxwell.cpp reports how to
      use this class for Maxwell problems.
    - ml_example_epetra_preconditioner_AztecMSR.cpp shows how to
      convert an Aztec matrix (in MSR format, but the example can be
      easily modified for VBR matrix formats) into Epetra matrices.
      
    \warning The Maxwell interface is still under development. 

    \author Marzio Sala, SNL 9214
*/  
class MultiLevelPreconditioner : public virtual Epetra_RowMatrix {
      
public:  

  //@{ \name Destructor.

  //! Constructs an MultiLevelPreconditioner with default values.

  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
                           const bool ComputePrec );

#ifdef HAVE_ML_TRIUTILS
  //! Constructs an MultiLevelPreconditioner, input parameters are specific in the command line

  MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
			   Trilinos_Util::CommandLineParser & CLP,
                           const bool ComputePrec );
#endif
  
  //! Constructs an MultiLevelPreconditioner. Retrives parameters (with prefix \c Prefix) from \c List.
  
  MultiLevelPreconditioner( const Epetra_RowMatrix & RowMatrix,
			    const Teuchos::ParameterList & List,
			    const bool ComputePrec=true, const char Prefix[]="");

  //! Constructs an MultiLevelPreconditioner from an ML_Operator. Retrives parameters from \c List.
  
  MultiLevelPreconditioner( ML_Operator * Operator,
			    const Teuchos::ParameterList & List,
			    const bool ComputePrec=true,
			    const char Prefix[]="" );
  
  //! Constructs an MultiLevelPreconditioner for Maxwell equations. Retrives parameters from \c List.
  /*! Constructs an MultiLevelPreconditioner for Maxwell equations. The constructor
    requires the edge matrix, the connectivity matrix T, the nodal matrix.
  */
  MultiLevelPreconditioner( const Epetra_RowMatrix & EdgeMatrix,
			    const Epetra_RowMatrix & TMatrix,
			    const Epetra_RowMatrix & NodeMatrix,
			    const Teuchos::ParameterList & List,
			    const bool ComputePrec=true,
			    const char Prefix[]="");

  //@}
  
  //@{ \name Destructor.

  ~MultiLevelPreconditioner() {
    if( IsComputePreconditionerOK_ ) Destroy_ML_Preconditioner(); 
  }

  //@}
  
  //@{ \name Query functions

  //! Prints label associated to this object.
  char * Label() const{return(Label_);};  
  
  //! Prints unused parameters in the input ParameterList on standard output. */
  void PrintUnused() const
  {
    List_.unused(std::cout);
  }

  //! Prints unused parameters in the input ParameterList on the specified stream.
  void PrintUnused(ostream & os) const
  {
    List_.unused(os);
  }

  //! Prints unused parameters in the input ParameterList to cout on proc \c MyPID. 
  /*! Mispelled parameters are simply ignored. Therefore, it is often the best
   * choice to print out the parameters that have not been used in the
   * construction phase. 
   * - \param MyPID (In) : ID of process that should print the unused parameters.
   */
  void PrintUnused(const int MyPID) const;

  //! Gets a reference to the internally stored parameters' list.
  Teuchos::ParameterList & GetList() 
  {
    return List_;
  }

  // Get a copy of the internally stored output list.
  Teuchos::ParameterList GetOutputList() 
  {
    return OutputList_;
  }

  //! Prints on \c cout the values of the internally stored parameter list for processor \c MyPID
  void PrintList(int MyPID);

  //! Copies \c List into the internally stored parameter list object.
  int SetParameterList(const Teuchos::ParameterList & List);

  //@}
  
  //@{ \name Mathematical functions.

  //! Apply the inverse of the preconditioner to an Epetra_MultiVector (NOT AVAILABLE)
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  //! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //@}
  
  //@{ \name Atribute access functions


  //! Computes the multilevel hierarchy.
  /*! Computes the multilevel hierarchy. This function retrives the user's defines parameters (as
    specified in the input ParameterList), or takes default values otherwise, and creates the ML
    objects for aggregation and hierarchy. Allocated data can be freed used DestroyPreconditioner(),
    or by the destructor,

    In a Newton-type procedure, several linear systems have to be solved, Often, these systems
    are not too different. In this case, it might be convenient to keep the already 
    computed preconditioner (with hierarchy, coarse solver, smoothers), and use it to
    precondition the next linear system. ML offers a way to determine whether the 
    already available preconditioner is "good enough" for the next linear system. 
    The user should proceed as follows:
    - define \c "adaptive: enable" == \c true
    - solve the first linear system. ML tries to estimate the rate of convergence, and record it;
    - change the values of the linear system matrix (but NOT its structure)
    - compute the new preconditioner as \c ComputePreconditioner(true)
    It is supposed that the pointer to the Epetra_RowMatrix remains constant. Currently,
    it is not possible to modify this pointer (other than creating a new preconditioner)
  */
  
  int ComputePreconditioner(const bool CheckFiltering = false);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  int ComputeFilteringPreconditioner();
#endif
  
  //! Queries whether multilevel hierarchy has been computed or not.
  int IsPreconditionerComputed()  const
  {
    return( IsComputePreconditionerOK_ );
  }

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose){return(-1);}

  //! Returns the infinity norm (not implemented).
  double NormInf() const {return(0.0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};

  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(*Comm_);};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(*DomainMap_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(*RangeMap_);};
  //@}

  //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
  void DestroyPreconditioner() 
  {
    if( IsComputePreconditionerOK_ ) Destroy_ML_Preconditioner();
  }

  //! Destroyes all structures allocated in \c ComputePreconditioner().   
  void Destroy_ML_Preconditioner();

  //!  Returns a reference to RowMatrix->Map().
  const Epetra_BlockMap & Map() const
  {
    return RowMatrix_->Map();
  }

  int NumMyRowEntries(int row, int &rows) const
  {
    return( RowMatrix_->NumMyRowEntries(row,rows) );
  }

  int MaxNumEntries() const
  {
    return( RowMatrix_->MaxNumEntries() );
  }

  int ExtractMyRowCopy(int i, int j, int & k, double * vals, int * indices) const
  {
    return( RowMatrix_->ExtractMyRowCopy(i,j,k,vals,indices) );
  }

  int ExtractDiagonalCopy(Epetra_Vector & v) const 
  {
    return( RowMatrix_->ExtractDiagonalCopy(v) );
  }

   int Multiply(bool flag, const Epetra_MultiVector & x,
		Epetra_MultiVector & y) const
  {
    return( RowMatrix_->Multiply(flag,x,y) );
  }
  
  int Solve(bool f1, bool f2, bool f3, const Epetra_MultiVector & x,
	    Epetra_MultiVector & y) const
  {
    return( RowMatrix_->Solve(f1,f2,f3,x,y) );
  }
  
  int InvRowSums(Epetra_Vector & x) const
  {
    return( RowMatrix_->InvRowSums(x) );;
  }

  int LeftScale(const Epetra_Vector & x)
  {
    //return( RowMatrix_->LeftScale(x) );
    return -1;
    
  }

  int InvColSums(Epetra_Vector & x) const
  {
    return( RowMatrix_->InvColSums(x) );
  }

  int  RightScale(const Epetra_Vector & x)
  {
    //    return( RowMatrix_->RightScale(x) );
    return -1;
  }
  
  bool Filled() const 
  {
    return( RowMatrix_->Filled() );
  }
  
  double NormOne() const 
  {
    return( RowMatrix_->NormOne() );
  }
  
  int NumGlobalNonzeros() const
  {
    return( RowMatrix_->NumGlobalNonzeros() );
  }
  
  int NumGlobalRows() const
  {
    return( RowMatrix_->NumGlobalRows() );
  }
  
  int NumGlobalCols() const
  {
    return( RowMatrix_->NumGlobalCols() );
  }
  
  int NumGlobalDiagonals() const
  {
    return( RowMatrix_->NumGlobalDiagonals() );
  }

  int NumMyNonzeros() const
  {
    return( RowMatrix_->NumMyNonzeros() );
  }
  
  int NumMyRows() const
  {
    return( RowMatrix_->NumMyRows() );
  }
  
  int NumMyCols() const
  {
    return( RowMatrix_->NumMyCols() );
  }
  
  int NumMyDiagonals() const
  {
    return( RowMatrix_->NumMyDiagonals() );
  }
  
  bool UpperTriangular() const
  {
    return( RowMatrix_->UpperTriangular() );
  }
  
  bool LowerTriangular() const
  {
    return( RowMatrix_->LowerTriangular());
  }

  const Epetra_Import * RowMatrixImporter() const
  {
    return RowMatrix_->RowMatrixImporter();
  }
  
  const Epetra_Map& RowMatrixRowMap() const
  {
    return RowMatrix_->RowMatrixRowMap();
  }
  
  const Epetra_Map& RowMatrixColMap() const
  {
    return RowMatrix_->RowMatrixColMap();
  }

  //@{ \name debugging and other utilities

  //! Stops the code, waiting for a debugger to attach
  /*! BreakForDebugger() is useful when the user wants to attach to the running
   * process(es). This is a very easy task for serial runs -- just run gdb.
   * Parallel runs may result more problematic. In this case, one can proceed as
   * follows:
   * - define the enviromental variable ML_BREAK_FOR_DEBUGGER (example, in BASH,
   *   \c export \c ML_BREAK_FOR_DEBUGGER=1 )
   * - run the parallel code on a terminal (example, \c mpirun \c -np \c 4 \c
   *   ml_example.exe )
   * - the code will stop in the first call to ComputePreconditioner(). This may
   *   occur in the construction phase. Few information about the ID number of
   *   each process will be showed.
   * - in another terminal, attach to the desired process.
   * - insert one character to let the code continue, and debug as required.
   */
  int BreakForDebugger();

  //! Prints the computational stencil for the specified row and equation (for 2D Cartesian grids only)
  /*! For problems defined on 2D Cartesian grids (with node numbering increasing
   * along the x-axis), this function prints out the stencil in an intelligible
   * form.
   * \param nx (In) : number of nodes along the X-axis
   * \param ny (In) : number of nodes along the Y-axis
   * \param NodeID (In) : (local) ID of node that will be used to print the
   *   stencil. If set to -1, the code will automatically chose an internal node.
   *   Default: -1.
   * \param EquationID (In) : ID of the equation that will be used to print the
   *   stencil (default = 0)  
   */
  int PrintStencil2D(const int nx, const int ny, 
		     int NodeID = -1,
		     const int EquationID = 0);
  //! Analyze the linear system matrix by testing several ML parameters
  int AnalyzeMatrix(char * Defaults, bool IsSymmetric = false);
  
  //@}

private:

  //! Copy constructor, should not be used
  MultiLevelPreconditioner(const MultiLevelPreconditioner & rhs) 
  {};

  //! operator =, should not be used.
  MultiLevelPreconditioner & operator = (const MultiLevelPreconditioner & rhs)
  {
    return *this;
  };

  //@{ \name Internal setting functions
  //! Initializes object with defauls values.
  void Initialize();

  //! Sets smoothers for non-Maxwell equations.
  void SetSmoothers();

  //! Sets smoothers for Maxwell equations.
  void SetSmoothersMaxwell();

  //! Sets coarse level solvers.
  void SetCoarse();

  //! Sets aggregation schemes.
  void SetAggregation();

  //! Sets preconditioner type (usually, V-cycle).
  void SetPreconditioner();

  //! Sets the null space for non-Maxwell problems.
  void SetNullSpace();

  //! Sets the null space for Maxwell equations.
  void SetNullSpaceMaxwell();

  //! Set parameters for eigen-computations.
  void SetEigenList();

  //! Sets prolongator smoother parameters.
  void SetSmoothingDamping();

  //! Sets damping parameter for classical smoothed aggregation.
  void SetSmoothingDampingClassic();
  
  //! Prints a line on cout.
  void PrintLine() const;

  //! Creates label for this object (printed out by AztecOO)
  int CreateLabel();

  void CreateAuxiliaryMatrix(Epetra_FECrsMatrix * & FakeMatrix);

  void PrintMem(char *fmt, int size, int, int);

  void PrintMemoryUsage();

  int SetFiltering();

  void VizMePleaze();
  
  void RandomAndZero(double *, double *, int);
  
  //! Checks whether the previously computed preconditioner is still valuable for the newly available linear system.
  /*! Used only when \c "adaptive: enable" is false and \c "filtering: enable" is true.
   * \warning: still under development
   */
  bool CheckPreconditionerFiltering();
  
  //! Checks whether the previously computed preconditioner is still valuable for the newly available linear system.
  /*! Used only when \c "adaptive: enable" is \c true, and
   * ComputePreconditioner(true) is called. */
  bool CheckPreconditionerKrylov();

  void AnalyzeMatrixProperties(char * Defaults, bool IsSymmetric = false);
  void AnalyzeMatrixEigenvalues(char * Defaults, bool IsSymmetric = false);
  void AnalyzeMatrixSmoothers(char * Defaults, bool IsSymmetric = false);
  void AnalyzeMatrixCoarsening(char * Defaults, bool IsSymmetric = false);

  //@}

  //@{ \name Internal data
  
  //! Pointer to ML_Struct
  ML * ml_;
  //! ML_Aggregate, contains aggregate information
  ML_Aggregate *agg_;
  //! Label for this object
  char * Label_;

  //! pointer to linear system matrix
  const Epetra_RowMatrix * RowMatrix_;
  //! specifies whether a hierarchy already exists or not.
  bool IsComputePreconditionerOK_;
  
  //! Number of levels
  int NumLevels_;
  //! Domain Map
  const Epetra_Map * DomainMap_;
  //! Range Map
  const Epetra_Map * RangeMap_;
  //! Epetra communicator object
  const Epetra_Comm * Comm_;
  bool  ownership_;
  //! proc_config for Aztec smoothers
  int   ProcConfig_[AZ_PROC_SIZE];
  //! options for Aztec smoothers
  int   SmootherOptions_[AZ_OPTIONS_SIZE];
  //! params for Aztec smoothers
  double SmootherParams_[AZ_PARAMS_SIZE];
  //! status for Aztec smoothers
  double SmootherStatus_[AZ_STATUS_SIZE];

  //! List containing all input parameters.
  Teuchos::ParameterList List_;
  //! List containing all output parameters
  Teuchos::ParameterList OutputList_;      
  //! List containing all the parameters for eigen-computations.
  Teuchos::ParameterList EigenList_;       

  //! Maximum number of levels
  int MaxLevels_;
  //! Integer array used to easily handle ML_INCREASING and ML_DECREASING
  /*! Integer array, of size MaxLevels_, that contain the ML level ID
    for the first logical level, and so on for all levels. The ML level ID
    of logical level L is LevelID_[L].
    In this interface, all levels move from 0 to MaxLevels-1.
    ML's level for interface's level i is LevelID_[i]
  */
  int * LevelID_;

  //! If not NULL, contains the allocated null space vector 
  double * NullSpaceToFree_;              

  //! All user's defined input data have this prefix
  char Prefix_[80];
  //! all cout's have this prefix (default'd in Initialize() )
  string PrintMsg_;
  //! all cerr's have this prefix (default'd in Initialize() )
  char ErrorMsg_[80];
  //! true of information has to be printed on this process
  bool verbose_;
  //! Number of PDE equations.
  int NumPDEEqns_;

  //@}

  //@{ \name Maxwell variables

  //! true if Maxwell equations are used
  bool SolvingMaxwell_;
  //! Main matrix for Maxwell
  const Epetra_RowMatrix * EdgeMatrix_;
  //! aux matrix for Maxwell
  const Epetra_RowMatrix * NodeMatrix_;
  //! T matrix for Maxwell
  const Epetra_RowMatrix * TMatrix_;
  ML_Operator * TMatrixML_;
  ML_Operator * TMatrixTransposeML_;
  ML_Operator ** Tmat_array, ** Tmat_trans_array;
  //! ML structures for Maxwell
  ML * ml_edges_, * ml_nodes_;

  void ** nodal_args_, ** edge_args_;

  //@}

  //@{ \name Variables for Timing
  //! Number of applications
  int NumApplications_;
  //! CPU time for all applications of the preconditioner
  double ApplicationTime_;
  bool FirstApplication_;
  //@ CPU time for first application
  double FirstApplicationTime_;
  //! Number of construction phases
  int NumConstructions_;
  //! CPU time for construction of the preconditioner.
  double ConstructionTime_;

  //@}
  
  // other stuff for old ML's compatibility
  Epetra_CrsMatrix * RowMatrixAllocated_;

  bool AnalyzeMemory_;
  
  int memory_[ML_MEM_SIZE];

  // filtering stuff

  Epetra_MultiVector * flt_R_;
  mutable Epetra_SerialDenseMatrix flt_A_;
  mutable Epetra_SerialDenseVector flt_rhs_, flt_lhs_;
  mutable Epetra_SerialDenseSolver flt_solver_;
  double * flt_NullSpace_;
  struct ML_CSR_MSRdata * flt_MatrixData_;
  ML * flt_ml_;
  ML_Aggregate * flt_agg_;
  
  // CheckPreconditioner related stuff
  Epetra_MultiVector       * SchurDecomposition_;
  Epetra_SerialDenseMatrix SchurMatrix_;
  Epetra_SerialDenseSolver SchurSolver_;
  double RateOfConvergence_;

}; // class MultiLevelPreconditioner
 
} // namespace ML_Epetra

#endif /* defined ML_EPETRA and ML_TEUCHOS */

#endif /* _ML_EPETRA_PRECONDITIONER_H_ */
