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

using namespace Teuchos;

class Epetra_ML_Time 
{
public:
  Epetra_ML_Time()
  {
    Reset();
  }
  
  inline void Update(double time, bool flag)
  {
    ApplicationTime_ += time;
    if( flag ) FirstApplicationTime_ = time;
    ++NumApplications_;
  }

  inline void Reset() 
  {
    ApplicationTime_ = 0.0;
    FirstApplicationTime_ = 0.0;
    NumApplications_ = 0;
  }

  inline double Application() const
  {
    return ApplicationTime_;
  }

  inline double FirstApplication() const
  {
    return FirstApplicationTime_;
  }

  inline int NumApplications() const
  {
    return NumApplications_;
  }
  
private:
  double ApplicationTime_;
  double FirstApplicationTime_;
  int NumApplications_;
};

  
class Epetra_ML_Preconditioner: public virtual Epetra_RowMatrix {
      
public:

  Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
			    ParameterList & List, bool);

  ~Epetra_ML_Preconditioner() {
    Destroy_ML_Preconditioner(); 
 }

  int ComputePreconditioner();

  int IsPreconditionerComputed()  const
  {
    return( IsComputePreconditionerOK_ );
  }
  
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
  const Epetra_Comm & Comm() const{return(Comm_);};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(DomainMap_);};
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(RangeMap_);};
  //@}

  void DestroyPreconditioner() 
  {
    Destroy_ML_Preconditioner();
  }
  
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
  
    
protected:

  ML * ml_;
  ML_Aggregate *agg_;
  
  char * Label_;

 private:

  const Epetra_RowMatrix * RowMatrix_;
  bool IsComputePreconditionerOK_;
  
  int CreateLabel();
  
  int NumLevels_;
  const Epetra_Map & DomainMap_;
  const Epetra_Map & RangeMap_;
  const Epetra_Comm & Comm_;
  bool  ownership_;
  int   ProcConfig_[AZ_PROC_SIZE];
  int   SmootherOptions_[AZ_OPTIONS_SIZE];
  double SmootherParams_[AZ_PARAMS_SIZE];
  double SmootherStatus_[AZ_STATUS_SIZE];
  ParameterList & List_;
  int MaxLevels_;
  Epetra_ML_Time * ML_Time_;
  
};

#endif /* defined ML_EPETRA and ML_TEUCHOS */

#endif /* _ML_EPETRA_PRECONDITIONER_H_ */
