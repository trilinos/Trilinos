#include "Ifpack_SORa.h"
#include "Ifpack.h"
#include "Ifpack_Utils.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Import.h"


using Teuchos::RefCountPtr;
using Teuchos::rcp;


#ifdef HAVE_IFPACK_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"


#define ABS(x) ((x)>=0?(x):-(x))

Ifpack_SORa::Ifpack_SORa(Epetra_RowMatrix* A):
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  Alpha_(1.5),
  Gamma_(1.0),
  NumSweeps_(1),
  IsParallel_(false),
  Time_(A->Comm())
{
  Epetra_CrsMatrix *Acrs=dynamic_cast<Epetra_CrsMatrix*>(A);
  TEST_FOR_EXCEPT(!Acrs) 
  A_=rcp(Acrs,false);
}

void Ifpack_SORa::Destroy(){
}



int Ifpack_SORa::Initialize(){
  Alpha_            = List_.get("sora: alpha", Alpha_);
  Gamma_            = List_.get("sora: gamma",Gamma_);
  NumSweeps_        = List_.get("sora: sweeps",NumSweeps_);

  if (A_->Comm().NumProc() != 1) IsParallel_ = true;
  else IsParallel_ = false;

  // Counters
  IsInitialized_=true;
  NumInitialize_++;
  return 0;
}

int Ifpack_SORa::SetParameters(Teuchos::ParameterList& parameterlist){
  List_=parameterlist;
  return 0;
}


int Ifpack_SORa::Compute(){
  if(!IsInitialized_) Initialize();

  
  Epetra_CrsMatrix *Askew2=0, *Aherm2=0,*W=0;
  int *rowptr_s,*colind_s,*rowptr_h,*colind_h;
  double *vals_s,*vals_h;
  Epetra_Vector Adiag(A_->RowMap());

  // Label
  sprintf(Label_, "IFPACK SORa (alpha=%5.2f gamma=%5.2f)",Alpha_,Gamma_); 

  // Create Askew2
  // Note: Here I let EpetraExt do the thinking for me.  Since it gets all the maps correct for the E + F^T stencil.
  // There are probably more efficient ways to do this but this method has the bonus of allowing code reuse.
  IFPACK_CHK_ERR(EpetraExt::MatrixMatrix::Add(*A_,false,1,*A_,true,-1,Askew2));
  Askew2->FillComplete();
  int nnz2=Askew2->NumMyNonzeros();
  int N=Askew2->NumMyRows();

  // Create Aherm2
  IFPACK_CHK_ERR(EpetraExt::MatrixMatrix::Add(*A_,false,1,*A_,true,1,Aherm2));
  Aherm2->FillComplete();

  // Grab pointers
  IFPACK_CHK_ERR(Askew2->ExtractCrsDataPointers(rowptr_s,colind_s,vals_s));
  IFPACK_CHK_ERR(Aherm2->ExtractCrsDataPointers(rowptr_h,colind_h,vals_h));

  // Sanity checking: Make sure the sparsity patterns are the same
#define SANITY_CHECK
#ifdef SANITY_CHECK
  for(int i=0;i<N;i++)
    if(rowptr_s[i]!=rowptr_h[i]) IFPACK_CHK_ERR(-2);
  for(int i=0;i<nnz2;i++)
    if(colind_s[i]!=colind_h[i]) IFPACK_CHK_ERR(-3);
#endif

  // Grab diagonal of A
  A_->ExtractDiagonalCopy(Adiag);

  // Allocate the diagonal for W
  Epetra_Vector *Wdiag = new Epetra_Vector(A_->RowMap());

  // Build the W matrix (strict lower triangle only)
  // Note: Relies on EpetraExt giving me identical sparsity patterns for both Askew2 and Aherm2 (see sanity check above)
  int maxentries=Askew2->MaxNumEntries();
  int* gids=new int [maxentries];
  double* newvals=new double[maxentries];
  W=new Epetra_CrsMatrix(Copy,A_->RowMap(),0);
  for(int i=0;i<N;i++){
    // Build the - (1+alpha)/2 E - (1-alpha)/2 F part of the W matrix
    int rowgid=A_->GRID(i);
    double c_data=0.0;
    int idx=0;
    for(int j=rowptr_s[i];j<rowptr_s[i+1];j++){      
      int colgid=Askew2->GCID(colind_s[j]);
      c_data+=fabs(vals_s[j]);
      if(rowgid>colgid){
	gids[idx]=colgid;
	newvals[idx]=vals_h[j]/2 + Alpha_ * vals_s[j]/2;
	idx++;
      }
    }
    IFPACK_CHK_ERR(W->InsertGlobalValues(rowgid,idx,newvals,gids));
    // Do the diagonal
    double w_val= c_data*Alpha_*Gamma_/4 + Adiag[A_->LRID(rowgid)];

    // Note: I drop a zero on the diagonal just to make sure that my rowmap is a subset of my column map.
    double zero=0.0;
    W->InsertGlobalValues(rowgid,1,&zero,&rowgid);//HAQ
    IFPACK_CHK_ERR(Wdiag->ReplaceGlobalValues(1,&w_val,&rowgid));
  }
  W->FillComplete(A_->DomainMap(),A_->RangeMap());      
  W_=rcp(W);
  Wdiag_=rcp(Wdiag);

  // Cleanup
  delete [] gids;
  delete [] newvals;
  delete Aherm2;
  delete Askew2;

  // Counters
  IsComputed_=true;
  NumCompute_++;
  ComputeTime_ += Time_.ElapsedTime();
  return 0;
}


int Ifpack_SORa::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{  
  if(!IsComputed_) return -1;
  Time_.ResetStartTime();
  bool initial_guess_is_zero=false;  
  int NumMyRows=W_->NumMyRows();
  int NumVectors = X.NumVectors();

  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0]){
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
    // Since the user didn't give us anything better, our initial guess is zero.
    Y.Scale(0.0);
    initial_guess_is_zero=true;
  }
  else
    Xcopy = Teuchos::rcp( &X, false );

  Epetra_MultiVector Temp(Y);
  
  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  // Note: Assuming that the matrix has an importer.  Ifpack_PointRelaxation doesn't do this, but given that 
  // I have a CrsMatrix, I'm probably OK.
  if (IsParallel_)  Y2 = Teuchos::rcp( new Epetra_MultiVector(W_->Importer()->TargetMap(),NumVectors));
  else Y2 = Teuchos::rcp( &Y, false );

  // Pointer grabs
  int* rowptr,*colind;
  double *values;
  double **t_ptr,** y_ptr, ** y2_ptr, **x_ptr,*d_ptr;
  Y2->ExtractView(&y2_ptr);
  Y.ExtractView(&y_ptr);
  Temp.ExtractView(&t_ptr);
  Xcopy->ExtractView(&x_ptr);
  Wdiag_->ExtractView(&d_ptr);
  IFPACK_CHK_ERR(W_->ExtractCrsDataPointers(rowptr,colind,values));

  double norm;

  for(int i=0; i<NumSweeps_; i++){
    // Import Y2 if parallel
    if(IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*W_->Importer(),Insert));
    
    // Calculate Ax 
    if(!initial_guess_is_zero  || i > 0) A_->Apply(Y,Temp);
    else Temp.Scale(0.0);

    // Do backsolve & update
    // x = x  + W^{-1} (b - A x)
    // NTS: This works since W_ has only the strict lower triangle and Wdiag_ has the diagonal
    for(int j=0; j<NumMyRows; j++){
      double diag=d_ptr[j];
      for (int m=0 ; m<NumVectors; m++) {
	double dtmp=0.0;
	for(int k=rowptr[j];k<rowptr[j+1];k++){
	  dtmp+= values[k]*y2_ptr[m][colind[k]];
	}
	// Yes, we need to update both of these.
	y2_ptr[m][j] += (x_ptr[m][j] - t_ptr[m][j] - dtmp)/diag;     
	y_ptr[m][j] = y2_ptr[m][j];
      }
    }
  }
  
  // Counter update
  NumApplyInverse_++;
  ApplyInverseTime_ += Time_.ElapsedTime();
  return 0;
}


ostream& Ifpack_SORa::Print(ostream& os) const{
  os<<"Ifpack_SORa"<<endl;
  os<<" alpha = "<<Alpha_<<endl;
  os<<" gamma = "<<Gamma_<<endl;
  os<<endl;
  return os;
}


double Ifpack_SORa::Condest(const Ifpack_CondestType CT, 
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  return -1.0;
}

#endif
