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

using namespace Teuchos;

//=============================================================================
Amesos_Taucs::Amesos_Taucs(const Epetra_LinearProblem &prob) :
  Problem_(&prob),
  Matrix_(0),
  A_(0),
  L_(0)
{ }

//=============================================================================
Amesos_Taucs::~Amesos_Taucs(void) 
{

  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();

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

//=============================================================================
int Amesos_Taucs::ConvertToSerial() 
{
  ResetTime();

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

  AddTime("matrix redistribution");

  return 0;
}

//=============================================================================
int Amesos_Taucs::ConvertToTaucs()
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    int n = SerialMatrix().NumMyRows();
    int nnz = SerialMatrix().NumMyNonzeros();
    int nnz_sym = (nnz) / 2 + n;

    if (A_ != 0)
    {
      taucs_dccs_free(A_);
      A_ = 0;
    }

    A_ = taucs_ccs_create(n, n, nnz_sym, TAUCS_DOUBLE);
    A_->flags = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;

    int count = 0;
    int MaxNumEntries = SerialMatrix().MaxNumEntries();
    vector<int>    Indices(MaxNumEntries);
    vector<double> Values(MaxNumEntries);

    A_->colptr[0] = 0;

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

        A_->rowind[count] = Indices[j];
        A_->values.d[count] = Values[j];
        ++count;
        ++count2;
      }

      A_->colptr[i + 1] = A_->colptr[i] + count2;
    }

    if (count > SerialMatrix().NumMyNonzeros())
      AMESOS_CHK_ERR(-1); // something wrong here
  }

  AddTime("conversion");

  return 0;
}

//=============================================================================
int Amesos_Taucs::SetParameters(Teuchos::ParameterList &ParameterList) 
{
  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);

  // compute norms of some vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",false);

  // compute the true residual Ax-b after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",false);

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",1);

  return 0;
}

//=============================================================================
int Amesos_Taucs::PerformSymbolicFactorization() 
{
  ResetTime();

  if (Comm().MyPID() == 0)
  {
    if (L_ != 0)
    {
      taucs_supernodal_factor_free(L_);
      L_ = 0;
    }

    L_ = (taucs_ccs_matrix*)taucs_ccs_factor_llt_symbolic(A_);
  }

  AddTime("symbolic");

  return 0;
}

//=============================================================================
int Amesos_Taucs::PerformNumericFactorization( ) 
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    taucs_supernodal_factor_free_numeric(L_);

    int ierr = taucs_ccs_factor_llt_numeric(A_, L_);

    if (ierr != 0) 
    {
      cerr << "Amesos_Taucs: error during symbolic factorization ("
        << ierr << ")" << endl;
      AMESOS_CHK_ERR(-1);
    }
  }

  AddTime("numeric");

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
  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Entering SymbolicFactorization()" << endl ; 
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  InitTime(Comm());

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

  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Leaving SymbolicFactorization()" << endl ; 
  return(0);
}

//=============================================================================
int Amesos_Taucs::NumericFactorization() 
{
  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Entering NumericFactorization()" << endl ; 
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  ++NumNumericFact_;

  // FIXME: this must be checked, now all the matrix is shipped twice here
  ConvertToSerial();
  ConvertToTaucs();

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;

  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Leaving NumericFactorization()" << endl ; 
  return(0);
}

//=============================================================================
int Amesos_Taucs::Solve() 
{
  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Entering Solve()" << endl ; 
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
  RefCountPtr<Epetra_MultiVector> SerialB;
  RefCountPtr<Epetra_MultiVector> SerialX;

  ResetTime();

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

  AddTime("vector redistribution");

  ResetTime();

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
      int ierr = taucs_supernodal_solve_llt(L_, 
                                            SerialXValues + i * LDA,
                                            SerialBValues + i * LDA);

      if (ierr != TAUCS_SUCCESS)
      {
        cerr << "Error occurred in taucs_ccs_solve()" << endl;
        AMESOS_CHK_ERR(-1);
      }
    }
  }

  AddTime("solve");

  ResetTime();

  if (Comm().NumProc() != 1) 
    X->Export(*SerialX, Importer(), Insert);

  AddTime("vector redistribution");

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *X, *B, false, "Amesos_Taucs");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_Taucs");

  if ( verbose_ > 1 ) cout << __FILE__ << "::" << __LINE__ << " Leaving Solve()" << endl ; 
  return(0) ;
}

// ====================================================================== 
void Amesos_Taucs::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  string p = "Amesos_Taucs : ";
  PrintLine();

  const Epetra_RowMatrix& A = Matrix() ;
  if ( Matrix_ != 0 ) { 
    
    int n = Matrix().NumGlobalRows();
    int nnz = Matrix().NumGlobalNonzeros();
    
    cout << p << "Matrix has " << n << " rows"
	 << " and " << nnz << " nonzeros" << endl;
    if ( n > 0 ) { 
      cout << p << "Nonzero elements per row = "
	   << 1.0 *  nnz / n << endl;
      cout << p << "Percentage of nonzero elements = "
	   << 100.0 * nnz /(pow(n,2.0)) << endl;
    }
    
    PrintLine();
  }

  return;
}

// ====================================================================== 
void Amesos_Taucs::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime("conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Taucs : ";
  PrintLine();

  cout << p << "Time to convert matrix to Taucs format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = "
       << SymTime << " (s), avg = " << SymTime << " (s)" << endl;
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
