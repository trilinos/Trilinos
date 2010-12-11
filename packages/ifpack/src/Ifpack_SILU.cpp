//@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_SILU.h"
#ifdef HAVE_IFPACK_SUPERLU

#include "Ifpack_CondestType.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

// SuperLU includes
extern "C" {int dfill_diag(int n, NCformat *Astore);}

using Teuchos::RefCountPtr;
using Teuchos::rcp;

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
#  include "Teuchos_TimeMonitor.hpp"
#endif

//==============================================================================
Ifpack_SILU::Ifpack_SILU(Epetra_RowMatrix* Matrix_in) :
  A_(rcp(Matrix_in,false)),
  Comm_(Matrix_in->Comm()),
  UseTranspose_(false),
  NumMyDiagonals_(0),
  DropTol_(1e-4),
  FillTol_(1e-2),
  FillFactor_(10.0),
  DropRule_(9), 
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  Time_(Comm()),
  etree_(0),
  perm_r_(0),
  perm_c_(0)
{
  Teuchos::ParameterList List;
  SetParameters(List);
  SY_.Store=0;
}



//==============================================================================
void Ifpack_SILU::Destroy()
{
  if(IsInitialized_){
    // SuperLU cleanup
    StatFree(&stat_);

    Destroy_CompCol_Permuted(&SAc_);
    Destroy_SuperNode_Matrix(&SL_);
    Destroy_CompCol_Matrix(&SU_);

    // Make sure NOT to clean up Epetra's memory
    Destroy_SuperMatrix_Store(&SA_);
    if(SY_.Store) Destroy_SuperMatrix_Store(&SY_);
    SY_.Store=0;

    // Cleanup stuff I allocated
    delete [] etree_;etree_=0;
    delete [] perm_r_;perm_r_=0;
    delete [] perm_c_;perm_c_=0;  
  }
}

//==========================================================================
int Ifpack_SILU::SetParameters(Teuchos::ParameterList& List)
{
  DropTol_=List.get("fact: drop tolerance",DropTol_);
  FillTol_=List.get("fact: zero pivot threshold",FillTol_);
  FillFactor_=List.get("fact: maximum fill factor",FillFactor_);
  DropRule_=List.get("fact: silu drop rule",DropRule_);

  // set label
  sprintf(Label_, "IFPACK SILU (drop=%d, zpv=%f, ffact=%f, rthr=%f)",
	  DropRule(),FillTol(),FillFactor(),DropTol());
  return(0);
}

//==========================================================================
int Ifpack_SILU::Initialize() 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_SILU::Initialize");
#endif

  Time_.ResetStartTime();

  // reset this object
  Destroy();

  IsInitialized_ = false;

  Epetra_CrsMatrix* CrsMatrix;
  CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(&*A_);

  if(CrsMatrix && CrsMatrix->RowMap().SameAs(CrsMatrix->ColMap()) && CrsMatrix->IndicesAreContiguous()){
    // Case #1: Use original matrix
    Aover_=rcp(CrsMatrix,false);
  }
  else if(CrsMatrix && CrsMatrix->IndicesAreContiguous()){
    // Case #2: Extract using CrsDataPointers w/ contiguous indices
    int size = A_->MaxNumEntries();
    int N=A_->NumMyRows();
    Aover_ = rcp(new Epetra_CrsMatrix(Copy,A_->RowMatrixRowMap(), size));
    vector<int> Indices(size);
    vector<double> Values(size);

    int i,j,ct,*rowptr,*colind;
    double *values;
    IFPACK_CHK_ERR(CrsMatrix->ExtractCrsDataPointers(rowptr,colind,values));

    // Use the fact that EpetraCrsMatrices always number the off-processor columns *LAST*   
    for(i=0;i<N;i++){
      for(j=rowptr[i],ct=0;j<rowptr[i+1];j++){
	if(colind[j]<N){
	  Indices[ct]=CrsMatrix->GCID(colind[j]);
	  Values[ct]=values[j];
	  ct++;
	}
      }
      Aover_->InsertGlobalValues(CrsMatrix->GRID(i),ct,&Values[0],&Indices[0]);
    }
    IFPACK_CHK_ERR(Aover_->FillComplete(CrsMatrix->RowMap(),CrsMatrix->RowMap()));  
  }
  else{
    // Case #3: Extract using copys
    int size = A_->MaxNumEntries();
    Aover_ = rcp(new Epetra_CrsMatrix(Copy,A_->RowMatrixRowMap(), size));
    if (Aover_.get() == 0) IFPACK_CHK_ERR(-5); // memory allocation error

    vector<int> Indices1(size),Indices2(size);
    vector<double> Values1(size),Values2(size);

    // extract each row at-a-time, and insert it into
    // the graph, ignore all off-process entries
    int N=A_->NumMyRows();
    for (int i = 0 ; i < N ; ++i) {
      int NumEntries;
      int GlobalRow = A_->RowMatrixRowMap().GID(i);
      IFPACK_CHK_ERR(A_->ExtractMyRowCopy(i, size, NumEntries, 
					  &Values1[0], &Indices1[0]));

      // convert to global indices, keeping only on-proc entries
      int ct=0;
      for (int j=0; j < NumEntries ; ++j) {
	if(Indices1[j] < N){
	  Indices2[ct] = A_->RowMatrixColMap().GID(Indices1[j]);
	  Values2[ct]=Values1[j];
	  ct++;
	} 
      }
      IFPACK_CHK_ERR(Aover_->InsertGlobalValues(GlobalRow,ct,
						&Values2[0],&Indices2[0]));
    }    
    IFPACK_CHK_ERR(Aover_->FillComplete(A_->RowMatrixRowMap(),A_->RowMatrixRowMap()));
  }

  // Finishing touches
  Aover_->OptimizeStorage();
  Graph_=rcp(const_cast<Epetra_CrsGraph*>(&Aover_->Graph()),false); 

  IsInitialized_ = true;
  NumInitialize_++;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack_SILU::Compute() 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_SILU::Compute");
#endif

  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  // Initialize the SuperLU statistics & options variables.
  StatInit(&stat_);
  ilu_set_default_options(&options_);
  options_.ILU_DropTol=DropTol_;
  options_.ILU_FillTol=FillTol_;
  options_.ILU_DropRule=DropRule_;
  options_.ILU_FillFactor=FillFactor_;

  // Grab pointers from Aover_
  int *rowptr,*colind;
  double *values;
  IFPACK_CHK_ERR(Aover_->ExtractCrsDataPointers(rowptr,colind,values));
  int N=Aover_->NumMyRows();

  // Copy the data over to SuperLU land - mark as a transposed CompCol Matrix
  dCreate_CompCol_Matrix(&SA_,N,N,Aover_->NumMyNonzeros(),
			 values,colind,rowptr,SLU_NC,SLU_D,SLU_GE);

  // Fill any zeros on the diagonal
  // Commented out for now
  dfill_diag(N, (NCformat*)SA_.Store);

  // Allocate SLU memory
  etree_=new int [N];
  perm_r_=new int[N];
  perm_c_=new int[N];

  // Get column permutation
  int permc_spec=options_.ColPerm;
  if ( permc_spec != MY_PERMC && options_.Fact == DOFACT )
    get_perm_c(permc_spec,&SA_,perm_c_);

  // Preorder by column permutation
  sp_preorder(&options_, &SA_, perm_c_, etree_, &SAc_);

  // Call the factorization
  int panel_size = sp_ienv(1);
  int relax      = sp_ienv(2);
  int info=0;
  dgsitrf(&options_,&SAc_,relax,panel_size,etree_,NULL,0,perm_c_,perm_r_,&SL_,&SU_,&stat_,&info);
  if(info<0) IFPACK_CHK_ERR(info);

  IsComputed_ = true;
  NumCompute_++;
  ComputeTime_ += Time_.ElapsedTime();
  return 0;
}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_SILU::Solve(bool Trans, const Epetra_MultiVector& X, 
                      Epetra_MultiVector& Y) const 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_SILU::ApplyInverse - Solve");
#endif

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);
  int nrhs=X.NumVectors();
  int N=A_->NumMyRows();

  // Copy X over to Y
  Y=X;

  // Move to SuperLU land
  // NTS: Need to do epetra-style realloc-if-nrhs-changes thing
  if(SY_.Store) Destroy_SuperMatrix_Store(&SY_);
  SY_.Store=0;
  dCreate_Dense_Matrix(&SY_,N,nrhs,Y[0],N,SLU_DN,SLU_D,SLU_GE);

  int info;
  dgstrs(TRANS,&SL_,&SU_,perm_c_,perm_r_,&SY_,&stat_,&info);
  if(!info) IFPACK_CHK_ERR(info);

  return(info);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_SILU::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const 
{

  if (!IsComputed())
    IFPACK_CHK_ERR(-1);

  return(0);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_SILU::ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_SILU::ApplyInverse");
#endif

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_.ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr< const Epetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  IFPACK_CHK_ERR(Solve(Ifpack_SILU::UseTranspose(), *Xcopy, Y));

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
double Ifpack_SILU::Condest(const Ifpack_CondestType CT, 
                           const int MaxIters, const double Tol,
                              Epetra_RowMatrix* Matrix_in)
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_SILU::Condest");
#endif

  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_SILU::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_SILU: " << Label() << endl << endl;
    os << "Dropping rule      = "<< DropRule() << endl;
    os << "Zero pivot thresh  = "<< FillTol() << endl;
    os << "Max fill factor    = "<< FillFactor() << endl;
    os << "Drop tolerance     = "<< DropTol() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows     = " << A_->NumGlobalRows() << endl;
    if (IsComputed_) {
      // Internal SuperLU info
      int fnnz=0;
      if(SL_.Store) fnnz+=((SCformat*)SL_.Store)->nnz;
      if(SU_.Store) fnnz+=((NCformat*)SU_.Store)->nnz;
      int lufill=fnnz - A_->NumMyRows();
      os << "No. of nonzeros in L+U    = "<< lufill<<endl;
      os << "Fill ratio: nnz(F)/nnz(A) = "<<(double)lufill / (double)A_->NumMyNonzeros()<<endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
    if (ComputeTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}

#endif
