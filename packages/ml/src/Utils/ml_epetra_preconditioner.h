#ifndef _ML_EPETRA_PRECONDITIONER_H_
#define _ML_EPETRA_PRECONDITIONER_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_Comm;
class Epetra_CrsMatrix;

#ifdef HAVE_ML_TRIUTILS
#include "Trilinos_Util_CommandLineParser.h"
#endif

#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//! ML_Epetra: default namespace for all Epetra interfaces.

namespace ML_Epetra
{

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaults(string ProblemType, Teuchos::ParameterList & List, char * Prefix = "");
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsDD(Teuchos::ParameterList & List, char * Prefix = "");
  
  int SetDefaultsDD_LU(Teuchos::ParameterList & List, char * Prefix = "");
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.  
  int SetDefaultsDD_3Levels(Teuchos::ParameterList & List, char * Prefix = "");
  
  //! Sets default parameters for Maxwell's equations.
  int SetDefaultsMaxwell(Teuchos::ParameterList & List, char * Prefix = "");
  
  //! Sets classical smoothed aggregation.
  int SetDefaultsSA(Teuchos::ParameterList & List, char * Prefix = "");


//! MultiLevelPreconditioner: An implementation of the Epetra_RowMatrix class.
/*! MultiLevelPreconditioner class implements Epetra_RowMatrix using a
    an Epetra_RowMatrix, and possibly a Teuchos parameters list, that
    specifies how to construct the preconditioner.
    The resulting preconditioner is completely black-box.
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
  
  //! Prints unused parameters in the input ParameterList.
  void PrintUnused() const
  {
    List_.unused(std::cout);
  }

  //! Prints unused parameters in the input ParameterList.
  void PrintUnused(ostream & os) const
  {
    List_.unused(os);
  }

  //! Prints unused parameters in the input ParameterList to cout on proc \c MyPID. 
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

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //@}
  
  //@{ \name Atribute access functions


  //! Computes the multilevel hierarchy.
  /*! Computes the multilevel hierarchy. This function retrives the user's defines parameters (as
    specified in the input ParameterList), or takes default values otherwise, and creates the ML
    objects for aggregation and hierarchy. Allocated data can be freed used DestroyPreconditioner(). */
  int ComputePreconditioner();

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

  //! Sets level smoothers for non-Maxwell equations.
  void SetSmoothers();

  //! Sets level smoothers for Maxwell equations.
  void SetSmoothersMaxwell();

  //! Sets coarse level solvers.
  void SetCoarse();

  //! Sets aggregation schemes.
  void SetAggregation();

  //! Sets preconditioner type (usually, V-cycle).
  void SetPreconditioner();

  //! Sets the null space.
  void SetNullSpace();

  //! Sets the null space for Maxwell equations.
  void SetNullSpaceMaxwell();

  //! Set parameters for eigen-computations.
  void SetEigenList();

  //! Sets prolongator smoother parameters.
  void SetSmoothingDamping();

  //! Prints a line on cout.
  void PrintLine() const;

  //! Creates label for this object (printed out by AztecOO)
  int CreateLabel();

  //@}

  //@{ \name Internal data
  
  ML * ml_;                                 //! Pointer to ML_Struct
  ML_Aggregate *agg_;                       //! ML_Aggregate, contains aggregate information
  
  char * Label_;                            //! Label for this object

  const Epetra_RowMatrix * RowMatrix_;      //! pointer to linear system matrix
  //! specifies whether a hierarchy already exists or not.
  bool IsComputePreconditionerOK_;
  
  //! Number of levels
  int NumLevels_;
  const Epetra_Map * DomainMap_;            //! Domain Map
  const Epetra_Map * RangeMap_;             //! Range Map
  const Epetra_Comm * Comm_;                //! Epetra communicator object.
  bool  ownership_;
  int   ProcConfig_[AZ_PROC_SIZE];          //! proc_config for Aztec smoothers
  int   SmootherOptions_[AZ_OPTIONS_SIZE];  //! options for Aztec smoothers
  double SmootherParams_[AZ_PARAMS_SIZE];   //! params for Aztec smoothers
  double SmootherStatus_[AZ_STATUS_SIZE];   //! status for Aztec smoothers

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
  
}; // class MultiLevelPreconditioner
 
} // namespace ML_Epetra

#endif /* defined ML_EPETRA and ML_TEUCHOS */

#endif /* _ML_EPETRA_PRECONDITIONER_H_ */
