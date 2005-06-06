// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Amesos_Superlu.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"

namespace SLU
{
extern "C" {
#include "dsp_defs.h"
}
}
struct SLUData
{
  SLU::SuperMatrix A, B, X, L, U;
#ifdef USE_DGSTRF
  SLU::SuperMatrix AC;
#endif
  SLU::superlu_options_t SLU_options;
  SLU::mem_usage_t mem_usage;
  SLU::fact_t refactor_option ;         //  SamePattern or SamePattern_SameRowPerm 
};

using namespace SLU;
using namespace Teuchos;

//=============================================================================
Amesos_Superlu::Amesos_Superlu(const Epetra_LinearProblem &prob ):
  NumGlobalRows_(-1),
  NumGlobalNonzeros_(-1),
  UseTranspose_(false),
  FactorizationOK_(false),
  FactorizationDone_(false),
  iam_(0),
  SerialMap_(Teuchos::null), 
  SerialCrsMatrixA_(Teuchos::null), 
  ImportToSerial_(Teuchos::null),
  SerialMatrix_(0),
  RowMatrixA_(0)
{
  data_ = new SLUData();
  ferr_.resize(1);
  berr_.resize(1);

  Problem_ = &prob ; 

  dCreate_Dense_Matrix( &(data_->X), 
			0, 
			0, 
			&DummyArray[0],
			0, 
			SLU_DN, SLU_D, SLU_GE);
    
  dCreate_Dense_Matrix( &(data_->B), 
			0, 
			0, 
			&DummyArray[0],
			0, 
			SLU_DN, SLU_D, SLU_GE);
}

//=============================================================================
Amesos_Superlu::~Amesos_Superlu(void) 
{
  if (PrintTiming_) PrintTiming();
  if (PrintStatus_) PrintStatus();

  Destroy_SuperMatrix_Store(&data_->B);
  Destroy_SuperMatrix_Store(&data_->X);

  if (iam_ == 0) { 
    if (FactorizationDone_) { 
      Destroy_SuperMatrix_Store(&data_->A);
      Destroy_SuperNode_Matrix(&data_->L);
      Destroy_CompCol_Matrix(&data_->U);
    }
  }

  delete data_; 
}

// ====================================================================== 
int Amesos_Superlu::SetParameters(Teuchos::ParameterList &ParameterList) 
{
  // print some statistics (on process 0). Do not include timing
  if (ParameterList.isParameter("PrintStatus"))
    PrintStatus_ = ParameterList.get("PrintStatus", PrintStatus_);

  // print some statistics (on process 0). Do not include timing
  if (ParameterList.isParameter("PrintTiming"))
    PrintTiming_ = ParameterList.get("PrintTiming", PrintTiming_);

  if (ParameterList.isParameter("ComputeTrueResidual"))
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",ComputeTrueResidual_);

  if (ParameterList.isParameter("ComputeVectorNorms"))
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",ComputeVectorNorms_);

  /* the list is empty now
  if (ParameterList.isSublist("Superlu") ) {
    Teuchos::ParameterList SuperluParams = ParameterList.sublist("Superlu") ;
  }  
  */
  return 0;
}

// ====================================================================== 
bool Amesos_Superlu::MatrixShapeOK() const 
{ 
  if (GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
      GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints()) 
    return(false);

  return(true);
}

// ======================================================================
int Amesos_Superlu::ConvertToSerial() 
{ 
  RowMatrixA_ = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA_ == 0)
    AMESOS_CHK_ERR(-1); // cannot cast to RowMatrix

  iam_ = Comm().MyPID() ;

  const Epetra_Map &OriginalMap = RowMatrixA_->RowMatrixRowMap() ; 

  NumGlobalRows_ = RowMatrixA_->NumGlobalRows();
  NumGlobalNonzeros_ = RowMatrixA_->NumGlobalNonzeros();
  if (NumGlobalRows_ != RowMatrixA_->NumGlobalCols())
    AMESOS_CHK_ERR(-1); // only square matrices
  if (NumGlobalNonzeros_ <= 0)
    AMESOS_CHK_ERR(-2); // empty matrix

  int NumMyElements_ = 0;
  if (iam_ == 0) NumMyElements_ = NumGlobalRows_;

  // If the matrix is distributed, then brings it to processor zero.
  // This also requires the construction of a serial map, and
  // the exporter from distributed to serial.
  // Otherwise, simply take the pointer of RowMatrixA_ and
  // set it to SerialMatrix_.

  if (Comm().NumProc() == 1) // Bug #1411 - should recognize serial matrices even when NumProc() > 1
  { 
     SerialMatrix_ = RowMatrixA_;
  } 
  else 
  {
    SerialMap_ = rcp(new Epetra_Map(NumGlobalRows_, NumMyElements_, 0, Comm()));

    ImportToSerial_ = rcp(new Epetra_Import(SerialMap(), OriginalMap));

    SerialCrsMatrixA_ = rcp(new Epetra_CrsMatrix(Copy, SerialMap(), 0));
    SerialCrsMatrixA_->Import(*RowMatrixA_, ImportToSerial(), Add); 
    
    // MS // Set zero element if not present, possibly add
    // MS // something to the diagonal
    double AddToDiag = 0.0;
#if 0
    if (iam_ == 0)
    {
      int this_res ;
      for (int i = 0 ; i < NumGlobalRows_; i++ ) {
	//  I am not sure what the intent is here, 
	//  but InsertGlobalValues returns 1 meaning "out of room" 
	//  and as a result, we sum AddToDiag_ into this diagonal element
	//  a second time.
        if (this_res=SerialCrsMatrixA_->InsertGlobalValues(i, 1, &AddToDiag, &i)) {
	  cout << __FILE__ << "::" << __LINE__ 
	       << " this_res = " <<  this_res 
	       << endl ; 
          SerialCrsMatrixA_->SumIntoGlobalValues(i, 1, &AddToDiag, &i);
	}
      }
    }
#endif

    SerialCrsMatrixA_->FillComplete();
    SerialMatrix_ = SerialCrsMatrixA_.get();
  }
  
  return(0);
}

//
//  See also pre and post conditions in Amesos_Superlu.h
//  Preconditions:  
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalRows_ has been set to the dimension of the matrix
//    NumGlobalNonzeros_ has been set to the number of non-zeros in the matrix
//      (i.e. ConvertToSerial() has been called) 
//    FactorizationDone_ and FactorizationOK_  must be accurate
//
//  Postconditions:
//    data->A 
//    FactorizationDone_ = true
//
//  Issues:
//
//
int Amesos_Superlu::Factor()
{
  // Convert matrix to the form that Superlu expects (Ap, Ai, Aval) 
  // I suppose that the matrix has already been formed on processor 0
  // as SerialMatrix_.

  if (iam_ == 0) 
  {
    if (NumGlobalRows_ != SerialMatrix_->NumGlobalRows() ||
        NumGlobalRows_ != SerialMatrix_->NumGlobalCols() ||
        NumGlobalRows_ != SerialMatrix_->NumMyRows() ||
        NumGlobalRows_ != SerialMatrix_->NumMyCols() ||
        NumGlobalNonzeros_ != SerialMatrix_->NumGlobalNonzeros())
    {
      AMESOS_CHK_ERR(-1); // something fishy here
    }

    NumGlobalNonzeros_ = SerialMatrix_->NumGlobalNonzeros();
    Ap_.resize(NumGlobalRows_ + 1);
    Ai_.resize(EPETRA_MAX( NumGlobalRows_, NumGlobalNonzeros_)); 
    Aval_.resize(EPETRA_MAX( NumGlobalRows_, NumGlobalNonzeros_)); 

    int NzThisRow ;
    int Ai_index = 0 ; 
    int MyRow;
    double *RowValues;
    int *ColIndices;
    vector<int> ColIndicesV_;
    vector<double> RowValuesV_;
    int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();

    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    if ( SuperluCrs == 0 ) {
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
    }
    for ( MyRow = 0; MyRow < NumGlobalRows_ ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	EPETRA_CHK_ERR( SerialMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					  NzThisRow, &RowValuesV_[0], 
					  &ColIndicesV_[0] ) != 0 );
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }
      Ap_[MyRow] = Ai_index ; 
      for ( int j = 0; j < NzThisRow; j++ ) { 
	Ai_[Ai_index] = ColIndices[j] ; 
	Aval_[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    assert( NumGlobalRows_ == MyRow );
    Ap_[ NumGlobalRows_ ] = Ai_index ; 

    if ( FactorizationDone_ ) { 
      Destroy_SuperMatrix_Store(&data_->A);
      Destroy_SuperNode_Matrix(&data_->L);
      Destroy_CompCol_Matrix(&data_->U);
    }
      
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix( &(data_->A), NumGlobalRows_, NumGlobalRows_,
			    NumGlobalNonzeros_, &Aval_[0],
			    &Ai_[0], &Ap_[0], SLU_NR, SLU_D, SLU_GE );
  }
  return 0;
}   

// ====================================================================== 
int Amesos_Superlu::ReFactor()
{
  // Convert matrix to the form that Superlu expects (Ap, Ai, Aval) 
  // I suppose that ConvertToSerial() has already been called, 
  // there I have SerialMatrix_ stored at it should. 
  // Only processor 0 does something useful here.

  if (iam_ == 0) 
  {
    if (NumGlobalRows_ != SerialMatrix_->NumGlobalRows() ||
        NumGlobalRows_ != SerialMatrix_->NumGlobalCols() ||
        NumGlobalRows_ != SerialMatrix_->NumMyRows() ||
        NumGlobalRows_ != SerialMatrix_->NumMyCols() ||
        NumGlobalNonzeros_ != SerialMatrix_->NumGlobalNonzeros())
    {
      AMESOS_CHK_ERR(-1);
    }

    Ap_.resize(NumGlobalRows_+ 1);
    Ai_.resize(EPETRA_MAX( NumGlobalRows_, NumGlobalNonzeros_)); 
    Aval_.resize(EPETRA_MAX(NumGlobalRows_, NumGlobalNonzeros_)); 

    int NzThisRow ;
    int Ai_index = 0 ; 
    int MyRow;
    double *RowValues;
    int *ColIndices;
    vector<int> ColIndicesV_;
    vector<double> RowValuesV_;
    int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();

    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    if ( SuperluCrs == 0 ) {
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
    }
    for ( MyRow = 0; MyRow < NumGlobalRows_ ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	EPETRA_CHK_ERR( SerialMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					  NzThisRow, &RowValuesV_[0], 
					  &ColIndicesV_[0] ) != 0 );
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }

      if (Ap_[MyRow] != Ai_index) 
        AMESOS_CHK_ERR(-4);

      for (int j = 0; j < NzThisRow; j++) { 
	assert(Ai_[Ai_index] == ColIndices[j]) ;   // FIXME this may not work.   
	Aval_[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    assert( NumGlobalRows_ == MyRow );
    Ap_[ NumGlobalRows_ ] = Ai_index ; 
    
    assert ( FactorizationDone_ ) ; 
    Destroy_SuperMatrix_Store(&data_->A);
    Destroy_SuperNode_Matrix(&data_->L);
    Destroy_CompCol_Matrix(&data_->U);

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix( &(data_->A), NumGlobalRows_, NumGlobalRows_,
			    NumGlobalNonzeros_, &Aval_[0],
			    &Ai_[0], &Ap_[0], SLU_NR, SLU_D, SLU_GE );
  }
  return 0;
}
//
// ====================================================================== 
int Amesos_Superlu::SymbolicFactorization() 
{
  // nothing happens here, all it is done in NumericFactorization()
  FactorizationOK_ = false; 
  return 0;
}

// ======================================================================
int Amesos_Superlu::NumericFactorization() 
{
  InitTime(Comm());
  ResetTime();

  ConvertToSerial(); 

  perm_r_.resize(NumGlobalRows_);
  perm_c_.resize(NumGlobalRows_);
  etree_.resize(NumGlobalRows_);
  R_.resize(NumGlobalRows_);
  C_.resize(NumGlobalRows_);

  SLU::superlu_options_t& SLUopt =  data_->SLU_options ; 

  set_default_options( &SLUopt ) ; 
  if (FactorizationOK_) {
    ReFactor() ; 
    SLUopt.Fact = data_->refactor_option ;
  }  else { 
    Factor() ; 
    FactorizationOK_ = true;
    SLUopt.Fact = DOFACT;
  }

  int Ierr[1];
  if ( iam_ == 0 ) { 
    double rpg, rcond;
    equed_ = 'N';

#if 0
    if ( ! UseTranspose() )        // FIXME - I doubt we need this here.
      assert( SLUopt.Trans == NOTRANS ) ; 
    else
      SLUopt.Trans = TRANS ; 
    

    //    SLUopt.ColPerm  = COLAMD ;

    cout << " SLUopt.ColPerm  = " << SLUopt.ColPerm  << endl ; 
    cout << " SLUopt.Equil  = " << SLUopt.Equil  << endl ; 
    cout << " SLUopt.Fact  = " << SLUopt.Fact  << endl ; 
    cout << " SLUopt.IterRefine  = " << SLUopt.IterRefine  << endl ; 
    cout << " data_->A.Stype  = " << data_->A.Stype  
	 << " SLU_NC = " << SLU_NC 
	 << " SLU_NR = " << SLU_NR 
	 << endl ; 
    cout << " SLUopt.ColPerm  = " << SLUopt.ColPerm  << endl ; 
#endif

    data_->B.nrow = NumGlobalRows_; 
    data_->B.ncol = 0;
    SLU::DNformat* Bstore = (SLU::DNformat *) (data_->B.Store) ; 
    Bstore->lda = NumGlobalRows_; 
    Bstore->nzval = &DummyArray[0];
    data_->X.nrow = NumGlobalRows_; 
    data_->X.ncol = 0;
    SLU::DNformat* Xstore = (SLU::DNformat *) (data_->X.Store) ; 
    Xstore->lda = NumGlobalRows_; 
    Xstore->nzval = &DummyArray[0];

    SuperLUStat_t SLU_stat ;
    StatInit( &SLU_stat ) ; 
    assert( SLUopt.Fact == DOFACT); 
    dgssvx( &(SLUopt), &(data_->A), 
	    &perm_c_[0], &perm_r_[0], &etree_[0], &equed_, &R_[0], 
	    &C_[0], &(data_->L), &(data_->U), NULL, 0, 
	    &(data_->B), &(data_->X), &rpg, &rcond, &ferr_[0], 
	    &berr_[0], &(data_->mem_usage), &SLU_stat,
	    &Ierr[0] );
    StatFree( &SLU_stat ) ; 
  }

  FactorizationDone_ = true; 

  AddTime("numeric");

  ++NumNumericFact_;

  return(0);
}

// ====================================================================== 
int Amesos_Superlu::Solve() 
{ 
  if (!FactorizationDone_) {
    FactorizationOK_ = false;
    AMESOS_CHK_ERR(NumericFactorization());
  }

  ResetTime();

  Epetra_MultiVector* vecX = Problem_->GetLHS(); 
  Epetra_MultiVector* vecB = Problem_->GetRHS(); 
  int Ierr;

  if (vecX == 0 || vecB == 0)  
    AMESOS_CHK_ERR(-1); // vectors not set

  int nrhs = vecX->NumVectors(); 
  if (nrhs != vecB->NumVectors())
    AMESOS_CHK_ERR(-2); // vectors not compatible

  ferr_.resize(nrhs); 
  berr_.resize(nrhs); 

  Epetra_MultiVector* SerialB = 0;
  Epetra_MultiVector* SerialX = 0;

  double *SerialXvalues ;
  double *SerialBvalues ;

  if (Comm().NumProc() == 1) 
  { 
    SerialB = vecB; 
    SerialX = vecX; 
  } 
  else 
  { 
    SerialX = new Epetra_MultiVector(SerialMap(), nrhs); 
    SerialB = new Epetra_MultiVector(SerialMap(), nrhs); 

    SerialB->Import(*vecB, ImportToSerial(), Insert);
    SerialB = SerialB; 
    SerialX = SerialX; 
  } 

  // Call SUPERLU's dgssvx to perform the solve
  // At this point I have, on processor 0, the solution vector
  // in SerialX and the right-hand side in SerialB
  
  if (iam_ == 0) 
  {
    int SerialXlda;
    int SerialBlda;

    int ierr;
    ierr = SerialX->ExtractView(&SerialXvalues, &SerialXlda);
    assert (ierr == 0);
    assert (SerialXlda == NumGlobalRows_ ) ; 

    ierr = SerialB->ExtractView(&SerialBvalues, &SerialBlda);
    assert (ierr == 0);
    assert (SerialBlda == NumGlobalRows_ ) ; 
    
    SLU::SuperMatrix& dataX = (data_->X) ;
    dataX.nrow =   NumGlobalRows_; 
    dataX.ncol =   nrhs ; 
    data_->X.nrow = NumGlobalRows_; 
    SLU::DNformat* Xstore = (SLU::DNformat *) (data_->X.Store) ; 
    Xstore->lda = SerialXlda; 
    Xstore->nzval = &SerialXvalues[0]; 

    SLU::SuperMatrix& dataB = (data_->B) ;
    dataB.nrow =   NumGlobalRows_; 
    dataB.ncol =   nrhs ; 
    data_->B.nrow = NumGlobalRows_; 
    SLU::DNformat* Bstore = (SLU::DNformat *) (data_->B.Store) ; 
    Bstore->lda = SerialBlda; 
    Bstore->nzval = &SerialBvalues[0]; 

    double rpg, rcond;
#if 0
    char fact, trans, refact;
    fact = 'F';
    refact = 'N';
    trans = 'N' ;    //  This will allow us to try trans and not trans - see if either works.
    //    equed = 'N';
#endif
#if 0
    dgssvx( &fact, &trans, &refact, &(data_->A), &(data_->iparam), &perm_c_[0],
	    &perm_r_[0], &etree_[0], &equed_, &R_[0], &C_[0], &(data_->L), &(data_->U),
	    NULL, 0, &(data_->B), &(data_->X), &rpg, &rcond,
	    &ferr_[0], &berr_[0], &(data_->mem_usage), &Ierr[0] );
#endif

    SuperLUStat_t SLU_stat ;
    StatInit( &SLU_stat ) ;//    Copy the scheme used in dgssvx1.c 
    data_->SLU_options.Fact = FACTORED ; 
    SLU::superlu_options_t& SLUopt =  data_->SLU_options ; 
    if (UseTranspose()) 
      SLUopt.Trans = TRANS; 
    else
      SLUopt.Trans = NOTRANS; 

    //#ifdef USE_DGSTRF
    //    assert( equed_ == 'N' ) ; 
    dgssvx( &(SLUopt), &(data_->A), 
	    &perm_c_[0], &perm_r_[0], &etree_[0], &equed_, &R_[0], 
	    &C_[0], &(data_->L), &(data_->U), NULL, 0, 
	    &(data_->B), &(data_->X), &rpg, &rcond, &ferr_[0], 
	    &berr_[0], &(data_->mem_usage), &SLU_stat,
	    &Ierr);
    //    assert( equed_ == 'N' ) ; 
    StatFree( &SLU_stat ) ; 
  }

  //
  //  Copy X back to the original vector
  // 
  if (Comm().NumProc() != 1) 
  { 
    vecX->Export(*SerialX, ImportToSerial(), Insert);
    delete SerialB;
    delete SerialX;
  } 

  //  All processes should return the same error code
  if (Comm().NumProc() != 1)
    Comm().Broadcast(&Ierr, 1, 0); 

#if 0 
  // MS // compute vector norms, as done in Amesos_Mumps
  if (ComputeVectorNorms_ == true) 
  {
    double NormLHS, NormRHS;
    for (int i = 0 ; i < nrhs ; ++i) 
    {
      AMESOS_CHK_ERR((*vecX)(i)->Norm2(&NormLHS));
      AMESOS_CHK_ERR((*vecB)(i)->Norm2(&NormRHS));
      if (Comm().MyPID() == 0) 
      {
	cout << "Amesos_Superlu : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }
  
  // MS // add compute true residual, as done in Amesos_Mumps
  if (ComputeTrueResidual_) 
  {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for (int i = 0 ; i < nrhs ; ++i) 
    {
      (RowMatrixA_->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));
      
      //
      //  In Amesos_Klu, verbose_ turns this print statement on 
      //
      if (false && Comm().MyPID() == 0) 
      {
	cout << "Amesos_Superlu : Vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }
  
  SolTime_+= Time_->ElapsedTime();
#else
  AddTime("solve");

  if (ComputeTrueResidual_)
    ComputeTrueResidual(*(GetProblem()->GetMatrix()), *vecX, *vecB, 
                        UseTranspose(), "Amesos_Superlu");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Superlu");

#endif
  ++NumSolve_;

  return(Ierr);
}

// ================================================ ====== ==== ==== == =

void Amesos_Superlu::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  string p = "Amesos_Superlu : ";
  PrintLine();

  int n = GetProblem()->GetMatrix()->NumGlobalRows();
  int nnz = GetProblem()->GetMatrix()->NumGlobalNonzeros();

  cout << p << "Matrix has " << n << " rows"
       << " and " << nnz << " nonzeros" << endl;
  cout << p << "Nonzero elements per row = "
       << 1.0 *  nnz / n << endl;
  cout << p << "Percentage of nonzero elements = "
       << 100.0 * nnz /(pow(n,2.0)) << endl;
  cout << p << "Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;
}

// ================================================ ====== ==== ==== == =
void Amesos_Superlu::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Superlu : ";
  PrintLine();

  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime << " (s), avg = " << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime << " (s), avg = " << SolTime << " (s)" << endl;

  PrintLine();

  return;
}
