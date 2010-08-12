#include "Ifpack_SORa.h"
#include "Ifpack.h"
#include "Ifpack_Utils.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_IntVector.h"

using Teuchos::RefCountPtr;
using Teuchos::rcp;


#ifdef HAVE_IFPACK_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"


#define ABS(x) ((x)>=0?(x):-(x))

// Useful functions horked from ML
int* FindLocalDiricheltRowsFromOnesAndZeros(const Epetra_CrsMatrix & Matrix, int &numBCRows);
void Apply_BCsToMatrixRowsAndColumns(const int *dirichletRows, int numBCRows,const Epetra_IntVector &dirichletColumns,const Epetra_CrsMatrix & Matrix);
Epetra_IntVector * FindLocalDirichletColumnsFromRows(const int *dirichletRows, int numBCRows,const Epetra_CrsMatrix & Matrix);

Ifpack_SORa::Ifpack_SORa(Epetra_RowMatrix* A):
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  Alpha_(1.5),
  Gamma_(1.0),
  NumSweeps_(1),
  IsParallel_(false),
  HaveOAZBoundaries_(false),
  UseInterprocDamping_(false),
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
  HaveOAZBoundaries_= List_.get("sora: oaz boundaries", HaveOAZBoundaries_);
  UseInterprocDamping_ = List_.get("sora: use interproc damping",UseInterprocDamping_);

  if (A_->Comm().NumProc() != 1) IsParallel_ = true;
  else {
    IsParallel_ = false;    
    // Don't use interproc damping, for obvious reasons
    UseInterprocDamping_=false;
  }

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

  // Dirichlet Detection & Nuking of Aherm2 and Askew2
  // Note: Relies on Aherm2/Askew2 having identical sparsity patterns (see sanity check above)
  if(HaveOAZBoundaries_){ 
    int numBCRows;
    int* dirRows=FindLocalDiricheltRowsFromOnesAndZeros(*A_,numBCRows);
    Epetra_IntVector* dirCols=FindLocalDirichletColumnsFromRows(dirRows,numBCRows,*Aherm2);
    Apply_BCsToMatrixRowsAndColumns(dirRows,numBCRows,*dirCols,*Aherm2);
    Apply_BCsToMatrixRowsAndColumns(dirRows,numBCRows,*dirCols,*Askew2);
    delete [] dirRows;
    delete dirCols;
  }

  // Grab diagonal of A
  A_->ExtractDiagonalCopy(Adiag);

  // Allocate the diagonal for W
  Epetra_Vector *Wdiag = new Epetra_Vector(A_->RowMap());

  // Build the W matrix (lower triangle only)
  // Note: Relies on EpetraExt giving me identical sparsity patterns for both Askew2 and Aherm2 (see sanity check above)
  int maxentries=Askew2->MaxNumEntries();
  int* gids=new int [maxentries];
  double* newvals=new double[maxentries];
  W=new Epetra_CrsMatrix(Copy,A_->RowMap(),0);
  for(int i=0;i<N;i++){
    // Build the - (1+alpha)/2 E - (1-alpha)/2 F part of the W matrix
    int rowgid=A_->GRID(i);
    double c_data=0.0;
    double ipdamp=0.0;
    int idx=0;

    for(int j=rowptr_s[i];j<rowptr_s[i+1];j++){      
      int colgid=Askew2->GCID(colind_s[j]);
      c_data+=fabs(vals_s[j]);
      if(rowgid>colgid){
	// Rely on the fact that off-diagonal entries are always numbered last, dropping the entry entirely.
	if(colind_s[j] < N) {       
	  gids[idx]=colgid;
	  newvals[idx]=vals_h[j]/2 + Alpha_ * vals_s[j]/2;
	  idx++;
	}
	else{
	  ipdamp+=fabs(vals_h[j]/2 + Alpha_ * vals_s[j]/2);
	}
      }   
    }
    if(idx>0)
      IFPACK_CHK_ERR(W->InsertGlobalValues(rowgid,idx,newvals,gids));

    // Do the diagonal
    double w_val= c_data*Alpha_*Gamma_/4 + Adiag[A_->LRID(rowgid)];
    if(UseInterprocDamping_) w_val+=ipdamp;

    W->InsertGlobalValues(rowgid,1,&w_val,&rowgid);
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
  Epetra_MultiVector Temp(A_->RowMap(),NumVectors);

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

  Teuchos::RefCountPtr< Epetra_MultiVector > T2;
  // Note: Assuming that the matrix has an importer.  Ifpack_PointRelaxation doesn't do this, but given that 
  // I have a CrsMatrix, I'm probably OK.
  // Note: This is the lazy man's version sacrificing a few extra flops for avoiding if statements to determine 
  // if things are on or off processor.
  // Note: T2 must be zero'd out
  if (IsParallel_ && W_->Importer())  T2 = Teuchos::rcp( new Epetra_MultiVector(W_->Importer()->TargetMap(),NumVectors,true));
  else T2 = Teuchos::rcp( new Epetra_MultiVector(A_->RowMap(),NumVectors,true));

  // Pointer grabs
  int* rowptr,*colind;
  double *values;
  double **t_ptr,** y_ptr, ** t2_ptr, **x_ptr,*d_ptr;
  T2->ExtractView(&t2_ptr);
  Y.ExtractView(&y_ptr);
  Temp.ExtractView(&t_ptr);
  Xcopy->ExtractView(&x_ptr);
  Wdiag_->ExtractView(&d_ptr);
  IFPACK_CHK_ERR(W_->ExtractCrsDataPointers(rowptr,colind,values));


  for(int i=0; i<NumSweeps_; i++){
    // Calculate b-Ax 
    if(!initial_guess_is_zero  || i > 0) {      
      A_->Apply(Y,Temp);
      Temp.Update(1.0,*Xcopy,-1.0);
    }
    else 
      Temp.Update(1.0,*Xcopy,0.0);

    // Note: The off-processor entries of T2 never get touched (they're always zero) and the other entries are updated 
    // in this sweep before they are used, so we don't need to reset T2 to zero here.

    // Do backsolve & update
    // x = x  + W^{-1} (b - A x)
    for(int j=0; j<NumMyRows; j++){
      double diag=d_ptr[j];
      for (int m=0 ; m<NumVectors; m++) {
	double dtmp=0.0;
	// Note: Since the diagonal is in the matrix, we need to zero that entry of T2 here to make sure it doesn't contribute.
	t2_ptr[m][j]=0.0;
	for(int k=rowptr[j];k<rowptr[j+1];k++){
	  dtmp+= values[k]*t2_ptr[m][colind[k]];
	}
	// Yes, we need to update both of these.
	t2_ptr[m][j] = (t_ptr[m][j]- dtmp)/diag;     
	y_ptr[m][j] += t2_ptr[m][j];
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







// ============================================================================
inline int* FindLocalDiricheltRowsFromOnesAndZeros(const Epetra_CrsMatrix & Matrix, int &numBCRows){
  int *dirichletRows = new int[Matrix.NumMyRows()];
  numBCRows = 0;
  for (int i=0; i<Matrix.NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j]!=0) nz++;      
      if (nz==1) dirichletRows[numBCRows++] = i;           
      // EXPERIMENTAL: Treat Special Inflow Boundaries as Dirichlet Boundaries
      if(nz==2) dirichletRows[numBCRows++] = i;        
    }/*end if*/
  }/*end for*/
  return dirichletRows;
}/*end FindLocalDiricheltRowsFromOnesAndZeros*/


// ====================================================================== 
 //! Finds Dirichlet the local Dirichlet columns, given the local Dirichlet rows
inline Epetra_IntVector * FindLocalDirichletColumnsFromRows(const int *dirichletRows, int numBCRows,const Epetra_CrsMatrix & Matrix){
  // Zero the columns corresponding to the Dirichlet rows.  Completely ignore the matrix maps.
  
  // Build rows
  Epetra_IntVector ZeroRows(Matrix.RowMap());
  Epetra_IntVector *ZeroCols=new Epetra_IntVector(Matrix.ColMap());
  ZeroRows.PutValue(0);  
  ZeroCols->PutValue(0);  

  // Easy Case: We're all on one processor
  if(Matrix.RowMap().SameAs(Matrix.ColMap())){
    for (int i=0; i < numBCRows; i++)
      (*ZeroCols)[dirichletRows[i]]=1;
    return ZeroCols;
  }

  // Flag the rows which are zero locally
  for (int i=0; i < numBCRows; i++)
    ZeroRows[dirichletRows[i]]=1;

  // Boundary exchange to move the data  
  if(Matrix.RowMap().SameAs(Matrix.DomainMap())){
    // Use A's Importer if we have one
    ZeroCols->Import(ZeroRows,*Matrix.Importer(),Insert);     
  }
  else{
    // Use custom importer if we don't
    Epetra_Import Importer(Matrix.ColMap(),Matrix.RowMap());
    ZeroCols->Import(ZeroRows,Importer,Insert);      
  }
  return ZeroCols;
}


// ====================================================================== 
inline void Apply_BCsToMatrixRowsAndColumns(const int *dirichletRows, int numBCRows,const Epetra_IntVector &dirichletColumns,const Epetra_CrsMatrix & Matrix){
  /* This function zeros out rows & columns of Matrix.
     Comments: The graph of Matrix is unchanged.
  */
  // Nuke the rows
  for(int i=0;i<numBCRows;i++){
    int numEntries, *cols;
    double *vals;
    int ierr = Matrix.ExtractMyRowView(dirichletRows[i],numEntries,vals,cols);
    for (int j=0; j<numEntries; j++) vals[j]=0.0; 
  }/*end for*/

  // Nuke the columns
  for (int i=0; i < Matrix.NumMyRows(); i++) {
    int numEntries;
    double *vals;
    int *cols;
    Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    for (int j=0; j < numEntries; j++) {
      if (dirichletColumns[ cols[j] ] > 0)  vals[j] = 0.0;
    }/*end for*/
  }/*end for*/
}/* end Apply_BCsToMatrixColumns */

#endif
