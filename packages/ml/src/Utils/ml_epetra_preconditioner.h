#ifndef _ML_EPETRA_PRECONDITIONER_H_
#define _ML_EPETRA_PRECONDITIONER_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Comm;

#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

namespace Epetra_ML
{

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaults(string ProblemType, Teuchos::ParameterList & List, char * Prefix = "",
		  int SmootherOptions[] = NULL, double SmootherParams[] = NULL);
  
  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  int SetDefaultsDD(Teuchos::ParameterList & List, char * Prefix = "",
		    int SmootherOptions[] = NULL, double SmootherParams[] = NULL);
  
  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.  
  int SetDefaultsDD_3Levels(Teuchos::ParameterList & List, char * Prefix = "",
			    int SmootherOptions[] = NULL, double SmootherParams[] = NULL);
  
  //! Sets default parameters for Maxwell's equations.
  int SetDefaultsMaxwell(Teuchos::ParameterList & List, char * Prefix = "" ,
			 int SmootherOptions[] = NULL, double SmootherParams[] = NULL);
  
  //! Sets classical smoothed aggregation.
  int SetDefaultsSA(Teuchos::ParameterList & List, char * Prefix = ""  ,
		    int SmootherOptions[] = NULL, double SmootherParams[] = NULL);

  
class MultiLevelPreconditioner : public virtual Epetra_RowMatrix {
      
public:  

  //! Constructs an MultiLevelPreconditioner with default values.

  MultiLevelPreconditioner::MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
						     const bool ComputePrec );

    //! Constructs an MultiLevelPreconditioner. Retrives parameters (with prefix \c Prefix) from \c List.
  
  MultiLevelPreconditioner( const Epetra_RowMatrix & RowMatrix,
			    const Teuchos::ParameterList & List,
			    const bool ComputePrec=true, const char Prefix[]="");

  //! Constructs an MultiLevelPreconditioner from an ML_Operator. Retrives parameters from \c List.
  
  MultiLevelPreconditioner::MultiLevelPreconditioner( ML_Operator * Operator,
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

  ~MultiLevelPreconditioner() {
    if( IsComputePreconditionerOK_ ) Destroy_ML_Preconditioner(); 
  }

  //! Prints unused parameters in the input ParameterList.
  void PrintUnused() 
  {
    List_.unused(std::cout);
  }

  //! Prints unused parameters in the input ParameterList.
  void PrintUnused(ostream & os) 
  {
    List_.unused(os);
  }

  void PrintUnused(int MyPID);

  Teuchos::ParameterList & GetList() 
  {
    return List_;
  }

  Teuchos::ParameterList & GetOutputList() 
  {
    return OutputList_;
  }

  //! Prints on \c cout the values of the internally stored parameter list for processor \c MyPID
  void PrintList(int MyPID);

  //! Copies \c List into the internally stored parameter list object.
  int SetParameterList(const Teuchos::ParameterList & List);

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
  
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};

  int SetUseTranspose(bool UseTranspose){return(-1);}

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  double NormInf() const {return(0.0);};

  char * Label() const{return(Label_);};
  
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

  //! Destroies all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
  void DestroyPreconditioner() 
  {
    if( IsComputePreconditionerOK_ ) Destroy_ML_Preconditioner();
  }

  //! Destroies all structures allocated in \c ComputePreconditioner().   
  void Destroy_ML_Preconditioner();

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
  
    
private:

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

  //! Set preconditioner type (usually, V-cycle).
  void SetPreconditioner();

  void SetNullSpace();

  void SetEigenList();
  
  void PrintLine();
  
  ML * ml_;                                 // ML_Struct
  ML_Aggregate *agg_;                       // ML_Aggregate, contains aggregate information
  
  char * Label_;

  const Epetra_RowMatrix * RowMatrix_;      // pointer to linear system matrix
  bool IsComputePreconditionerOK_;
  
  int CreateLabel();
  
  int NumLevels_;
  const Epetra_Map * DomainMap_;
  const Epetra_Map * RangeMap_;
  const Epetra_Comm * Comm_;
  bool  ownership_;
  int   ProcConfig_[AZ_PROC_SIZE];          // some Aztec's vectors
  int   SmootherOptions_[AZ_OPTIONS_SIZE];
  double SmootherParams_[AZ_PARAMS_SIZE];
  double SmootherStatus_[AZ_STATUS_SIZE];
  Teuchos::ParameterList List_;             // all input parameters are here
  Teuchos::ParameterList OutputList_;       // various informations

  Teuchos::ParameterList EigenList_;        // for advanced R constructions
  
  int MaxLevels_;
  int * LevelID_;                           // used to easily handle ML_INCREASING and ML_DECREASING.
                                            // In this interface, all levels move from 0 to MaxLevels-1.
                                            // ML's level for interface's level i is LevelID_[i]

  double * NullSpaceToFree_;                // not NULL if null space vectors have been allocated
  
  char Prefix_[80];                         // all user's defined input data have this prefix
  string PrintMsg_;                         // all cout's have this prefix (default'd in Initialize() )
  char ErrorMsg_[80];                       // all cerr's have this prefix (default'd in Initialize() )
  bool verbose_;
  int NumPDEEqns_;
  
  // Maxwell stuff
  bool SolvingMaxwell_;                     // true if Maxwell equations are used
  const Epetra_RowMatrix * EdgeMatrix_;     // Main matrix for Maxwell
  const Epetra_RowMatrix * NodeMatrix_;     // aux matrix for Maxwell
  const Epetra_RowMatrix * TMatrix_;
  const Epetra_RowMatrix * TMatrixTranspose_;
  ML_Operator ** Tmat_array, ** Tmat_trans_array;
  ML * ml_edges_, * ml_nodes_;

  // variables for Timing
  int NumApplications_;
  double ApplicationTime_; // all applications
  bool FirstApplication_;
  double FirstApplicationTime_; // only for first application
  int NumConstructions_;
  double ConstructionTime_;

  // other stuff for old ML's compatibility
  Epetra_CrsMatrix * RowMatrixAllocated_;
  
};
 
}

#endif /* defined ML_EPETRA and ML_TEUCHOS */

#endif /* _ML_EPETRA_PRECONDITIONER_H_ */
