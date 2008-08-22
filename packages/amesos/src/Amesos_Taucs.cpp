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

#include "Amesos_Taucs.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
extern "C" {
#include "taucs.h"
}

using namespace Teuchos;


//
//  deaalocFunctorHandleDelete requires a fucntion that takes a ** argument, these
//  functions meet that requirement
// 
void taucs_dccs_free_ptr( taucs_ccs_matrix** taucs_ccs_matrix_in ) { 
  taucs_dccs_free( *taucs_ccs_matrix_in ); 
}

void taucs_supernodal_factor_free_ptr( taucs_ccs_matrix** taucs_ccs_matrix_in ) { 
  taucs_supernodal_factor_free( *taucs_ccs_matrix_in ); 
}

class Amesos_Taucs_Pimpl {
public:

#define USE_REF_COUNT_PTR_FOR_A_AND_L
#ifdef USE_REF_COUNT_PTR_FOR_A_AND_L
   Teuchos::RCP<taucs_ccs_matrix> A_ ;
   Teuchos::RCP<taucs_ccs_matrix> L_ ;
#else
   taucs_ccs_matrix* A_ ;
   taucs_ccs_matrix* L_ ;

  Amesos_Taucs_Pimpl():
    A_(0),
    L_(0)
  {}

  ~Amesos_Taucs_Pimpl(void){
  if (A_ != 0)
  {
    taucs_dccs_free(A_);
    A_ = 0;
  }

  if (L_ != 0)
  {
    taucs_supernodal_factor_free(L_);
    L_ = 0;
  }
  }

#endif
} ;

//=============================================================================
Amesos_Taucs::Amesos_Taucs(const Epetra_LinearProblem &prob) :
  PrivateTaucsData_( rcp( new Amesos_Taucs_Pimpl() ) ),
  Matrix_(0),
  Problem_(&prob),
  MtxConvTime_(-1),
  MtxRedistTime_(-1),
  VecRedistTime_(-1),
  SymFactTime_(-1),
  NumFactTime_(-1),
  SolveTime_(-1)
{ }

//=============================================================================
Amesos_Taucs::~Amesos_Taucs(void) 
{

  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();

#if 0
  if (A_ != 0)
  {
    taucs_dccs_free(A_);
    A_ = 0;
  }

  if (L_ != 0)
  {
    taucs_supernodal_factor_free(L_);
    L_ = 0;
  }
#endif
}

//=============================================================================
int Amesos_Taucs::ConvertToSerial() 
{
  ResetTimer();

  int NumGlobalRows = Matrix_->NumGlobalRows();

  // create a serial map
  int NumMyRows = 0;
  if (Comm().MyPID() == 0) 
    NumMyRows = NumGlobalRows;

  SerialMap_ = rcp(new Epetra_Map(-1, NumMyRows, 0, Comm()));
  if (SerialMap_.get() == 0)
    AMESOS_CHK_ERR(-1);

  Importer_ = rcp(new Epetra_Import(SerialMap(),Map()));
  if (Importer_.get() == 0)
    AMESOS_CHK_ERR(-1);

  SerialCrsMatrix_ = rcp(new Epetra_CrsMatrix(Copy, SerialMap(), 0));
  if (SerialCrsMatrix_.get() == 0)
    AMESOS_CHK_ERR(-1);

  AMESOS_CHK_ERR(SerialCrsMatrix().Import(Matrix(), Importer(), Add));

  AMESOS_CHK_ERR(SerialCrsMatrix().FillComplete());

  SerialMatrix_ = rcp(SerialCrsMatrix_.get(), false);

  MtxRedistTime_ = AddTime("Total matrix redistribution time", MtxRedistTime_);

  return 0;
}

//=============================================================================
int Amesos_Taucs::ConvertToTaucs()
{
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    int n = SerialMatrix().NumMyRows();
    int nnz = SerialMatrix().NumMyNonzeros();
    int nnz_sym = (nnz) / 2 + n;

#ifndef USE_REF_COUNT_PTR_FOR_A_AND_L
    if (PrivateTaucsData_->A_ != 0)
    {
      taucs_dccs_free(PrivateTaucsData_->A_);
      PrivateTaucsData_->A_ = 0;
    }
#endif

#ifdef USE_REF_COUNT_PTR_FOR_A_AND_L
    PrivateTaucsData_->A_ = 
      rcp( taucs_ccs_create(n, n, nnz_sym, TAUCS_DOUBLE),
		 deallocFunctorHandleDelete<taucs_ccs_matrix>(taucs_dccs_free_ptr), true );	   
#else
    PrivateTaucsData_->A_ = taucs_ccs_create(n, n, nnz_sym, TAUCS_DOUBLE) ;
#endif
    PrivateTaucsData_->A_->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;

    int count = 0;
    int MaxNumEntries = SerialMatrix().MaxNumEntries();
    std::vector<int>    Indices(MaxNumEntries);
    std::vector<double> Values(MaxNumEntries);

    PrivateTaucsData_->A_->colptr[0] = 0;

    for (int i = 0 ; i < n ; ++i)
    {
      int count2 = 0;
      int ierr, NumEntriesThisRow;
      ierr = SerialMatrix().ExtractMyRowCopy(i, MaxNumEntries,
                                             NumEntriesThisRow,
                                             &Values[0], &Indices[0]);
      if (ierr < 0)
        AMESOS_CHK_ERR(ierr);

      for (int j = 0 ; j < NumEntriesThisRow ; ++j)
      {
        if (Indices[j] < i)
          continue;

        if (Indices[j] == i)
          Values[j] += AddToDiag_;

        PrivateTaucsData_->A_->rowind[count] = Indices[j];
        PrivateTaucsData_->A_->values.d[count] = Values[j];
        ++count;
        ++count2;
      }

      PrivateTaucsData_->A_->colptr[i + 1] = PrivateTaucsData_->A_->colptr[i] + count2;
    }

    if (count > SerialMatrix().NumMyNonzeros())
      AMESOS_CHK_ERR(-1); // something wrong here
  }

  MtxConvTime_ = AddTime("Total matrix conversion time", MtxConvTime_);

  return 0;
}

//=============================================================================
int Amesos_Taucs::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  return 0;
}

//=============================================================================
int Amesos_Taucs::PerformSymbolicFactorization() 
{
  ResetTimer();

  if (Comm().MyPID() == 0)
  {

#ifdef USE_REF_COUNT_PTR_FOR_A_AND_L
    PrivateTaucsData_->L_ = 
      rcp( ((taucs_ccs_matrix*)taucs_ccs_factor_llt_symbolic(&*PrivateTaucsData_->A_)),
	   deallocFunctorHandleDelete<taucs_ccs_matrix>(taucs_supernodal_factor_free_ptr), true );	  
#else
    if (PrivateTaucsData_->L_ != 0)
    {
      taucs_supernodal_factor_free(PrivateTaucsData_->L_);
      PrivateTaucsData_->L_ = 0;
    }
    PrivateTaucsData_->L_ = 
      (taucs_ccs_matrix*)taucs_ccs_factor_llt_symbolic(&*PrivateTaucsData_->A_) ;
#endif 
  }

  SymFactTime_ = AddTime("Total symbolic factorization time", SymFactTime_);

  return 0;
}

//=============================================================================
int Amesos_Taucs::PerformNumericFactorization( ) 
{
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    taucs_supernodal_factor_free_numeric( &*PrivateTaucsData_->L_);

    int ierr = taucs_ccs_factor_llt_numeric(&*PrivateTaucsData_->A_, &*PrivateTaucsData_->L_);

    if (ierr != 0) 
    {
      std::cerr << "Amesos_Taucs: error during numeric factorization ("
        << ierr << ")" << std::endl;
      AMESOS_CHK_ERR(-1);
    }
  }

  NumFactTime_ = AddTime("Total numeric factorization time", NumFactTime_);

  return 0;
}

//=============================================================================
bool Amesos_Taucs::MatrixShapeOK() const 
{
  bool OK = true;

  if (GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
      GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) 
  {
    OK = false;
  }
  return OK;
}

//=============================================================================
int Amesos_Taucs::SymbolicFactorization() 
{
  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__ << " Entering SymbolicFactorization()" << std::endl ; 
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  CreateTimer(Comm());

  ++NumSymbolicFact_;

  Matrix_ = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  Map_ = &(Matrix_->RowMatrixRowMap());

  // =========================================================== //
  // redistribute and create all import/export objects only      //
  // if more than one processor is used. Otherwise simply set    //
  // dummy pointers to Matrix() and Map(), without giving the    //
  // ownership to the smart pointer.                             //
  // =========================================================== //

  if (Comm().NumProc() != 1) 
    ConvertToSerial();
  else
  {
    SerialMap_ = rcp(const_cast<Epetra_Map*>(&Map()), false);
    SerialMatrix_ = rcp(const_cast<Epetra_RowMatrix*>(&Matrix()), false);
  }

  // =========================================================== //
  // Only on processor zero, convert the matrix into CSR format, //
  // as required by TAUCS.                                       //
  // =========================================================== //

  ConvertToTaucs();

  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;

  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__  << " Leaving SymbolicFactorization()" << std::endl ; 
  return(0);
}

//=============================================================================
int Amesos_Taucs::NumericFactorization() 
{
  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__ << " Entering NumericFactorization()" << std::endl ; 
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  ++NumNumericFact_;

  // FIXME: this must be checked, now all the matrix is shipped twice here
  ConvertToSerial();
  ConvertToTaucs();

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;

  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__
			   << " Leaving NumericFactorization()" << std::endl ; 
  return(0);
}

//=============================================================================
int Amesos_Taucs::Solve() 
{
  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__
			   << " Entering Solve()" << std::endl ; 
  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  ++NumSolve_;

  Epetra_MultiVector* X = Problem_->GetLHS();
  Epetra_MultiVector* B = Problem_->GetRHS();

  if ((X == 0) || (B == 0))
    AMESOS_CHK_ERR(-1); 

  int NumVectors = X->NumVectors();
  if (NumVectors != B->NumVectors())
    AMESOS_CHK_ERR(-1); 

  // vectors with SerialMap_
  RCP<Epetra_MultiVector> SerialB;
  RCP<Epetra_MultiVector> SerialX;

  ResetTimer();

  if (Comm().NumProc() == 1) 
  {
    SerialB = rcp(B,false);
    SerialX = rcp(X,false);
  } 
  else 
  {
    SerialX = rcp(new Epetra_MultiVector(SerialMap(),NumVectors));
    SerialB = rcp(new Epetra_MultiVector(SerialMap(),NumVectors));

    SerialB->Import(*B,Importer(),Insert);
  }

  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    double* SerialXValues;
    double* SerialBValues;
    int LDA;

    AMESOS_CHK_ERR(SerialX->ExtractView(&SerialXValues,&LDA));

    // FIXME: check LDA
    AMESOS_CHK_ERR(SerialB->ExtractView(&SerialBValues,&LDA));

    for (int i = 0 ; i < NumVectors ; ++i)
    {
      int ierr = taucs_supernodal_solve_llt(&*PrivateTaucsData_->L_, 
                                            SerialXValues + i * LDA,
                                            SerialBValues + i * LDA);

      if (ierr != TAUCS_SUCCESS)
      {
        std::cerr << "Error occurred in taucs_ccs_solve()" << std::endl;
        AMESOS_CHK_ERR(-1);
      }
    }
  }

  SolveTime_ = AddTime("Total solve time", SolveTime_);

  ResetTimer();

  if (Comm().NumProc() != 1) 
    X->Export(*SerialX, Importer(), Insert);

  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *X, *B, false, "Amesos_Taucs");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_Taucs");

  if ( debug_ > 0 ) std::cout << __FILE__ << "::" << __LINE__
			   << " Leaving Solve()" << std::endl ; 
  return(0) ;
}

// ====================================================================== 
void Amesos_Taucs::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  std::string p = "Amesos_Taucs : ";

  if (Matrix_ == 0) return;

  
  PrintLine();
  int n = Matrix().NumGlobalRows();
  int nnz = Matrix().NumGlobalNonzeros();
    
  std::cout << p << "Matrix has " << n << " rows"
       << " and " << nnz << " nonzeros" << std::endl;
  if (n > 0) 
  { 
    std::cout << p << "Nonzero elements per row = "
         << 1.0 *  nnz / n << std::endl;
    std::cout << p << "Percentage of nonzero elements = "
         << 100.0 * nnz /(pow(n,2.0)) << std::endl;
  }
  PrintLine();

  return;
}

// ====================================================================== 
void Amesos_Taucs::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime(MtxConvTime_);
  double MatTime = GetTime(MtxRedistTime_);
  double VecTime = GetTime(VecRedistTime_);
  double SymTime = GetTime(SymFactTime_);
  double NumTime = GetTime(NumFactTime_);
  double SolTime = GetTime(SolveTime_);

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  std::string p = "Amesos_Taucs : ";
  PrintLine();

  std::cout << p << "Time to convert matrix to Taucs format = "
       << ConTime << " (s)" << std::endl;
  std::cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << std::endl;
  std::cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << std::endl;
  std::cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << std::endl;
  std::cout << p << "Time for sym fact = "
       << SymTime << " (s), avg = " << SymTime << " (s)" << std::endl;
  std::cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << std::endl;
  std::cout << p << "Time for num fact = "
       << NumTime << " (s), avg = " << NumTime << " (s)" << std::endl;
  std::cout << p << "Number of solve phases = "
       << NumSolve_ << std::endl;
  std::cout << p << "Time for solve = "
       << SolTime << " (s), avg = " << SolTime << " (s)" << std::endl;

  PrintLine();

  return;
}
