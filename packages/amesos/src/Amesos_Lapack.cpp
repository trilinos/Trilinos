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

#include "Amesos_Lapack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"

//=============================================================================
Amesos_Lapack::Amesos_Lapack(const Epetra_LinearProblem &Problem) :
  SerialMatrix_(0),
  SerialMap_(0),
  UseTranspose_(false),
  Problem_(&Problem),
  PrintTiming_(false),
  PrintStatus_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  verbose_(1),
  Threshold_(0.0),
  AddToDiag_(0.0),
  RowImporter_(0),
  IsSymbolicFactorizationOK_(false),
  IsNumericFactorizationOK_(false),
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0),
  Time_(0)
{
  Teuchos::ParameterList ParamList;
  SetParameters(ParamList);
}

//=============================================================================
Amesos_Lapack::~Amesos_Lapack(void) {

  if (SerialMatrix_) 
    delete SerialMatrix_ ;

  if (SerialMap_) 
    delete SerialMap_ ;

  if (RowImporter_) 
    delete RowImporter_;

  if (Time_)
    delete Time_;

  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || (verbose_ == 2)) 
    PrintTiming();
  if ((verbose_ && PrintStatus_) || (verbose_ == 2)) 
    PrintStatus();
}


//=============================================================================
int Amesos_Lapack::SetParameters( Teuchos::ParameterList &ParameterList ) {

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  // ========================================= //
  // retrive parameters from list.             //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

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

  if( ParameterList.isParameter("Threshold") )
    Threshold_ = ParameterList.get("Threshold", 0.0);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);
 
  // MS // now comment it out, if we have parameters for LAPACK sublist
  // MS // uncomment it
  /*
  if (ParameterList.isSublist("Lapack") ) {
    Teuchos::ParameterList KluParams = ParameterList.sublist("Lapack") ;
  }
  */

  return 0;
}

//=============================================================================
bool Amesos_Lapack::MatrixShapeOK() const 
{
  if (GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints()) {
    return(false);
  }
  return(true);
}


//=============================================================================
int Amesos_Lapack::SymbolicFactorization() 
{
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  if (Time_ == 0) 
    Time_ = new Epetra_Time( Comm() );

  Time_->ResetStartTime();

  IsSymbolicFactorizationOK_ = true;
  SymTime_ += Time_->ElapsedTime();
  ++NumSymbolicFact_;

  return(0);
}

//=============================================================================
int Amesos_Lapack::NumericFactorization() 
{

  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK() == false)
    AMESOS_CHK_ERR(SymbolicFactorization()); // do nothing anyway

  if (MyPID() == 0) {
    AMESOS_CHK_ERR(DenseMatrix_.Reshape(1,1));
    AMESOS_CHK_ERR(DenseMatrix_.Reshape(NumGlobalRows(),NumGlobalRows()));
  }

  // convert the matrix from sparse to dense
  if (NumProc() == 1) {
    AMESOS_CHK_ERR(SerialToDense());
  }
  else {
    AMESOS_CHK_ERR(DistributedToDense());
  }

  Time_->ResetStartTime();

  if (MyPID() == 0) {
    AMESOS_CHK_ERR(DenseSolver_.SetMatrix(DenseMatrix_));
    AMESOS_CHK_ERR(DenseSolver_.Factor());
  }

  IsNumericFactorizationOK_ = true;
  NumTime_ += Time_->ElapsedTime();
  ++NumNumericFact_;
  
  return(0);
}


//=============================================================================
int Amesos_Lapack::Solve() 
{

  if (IsNumericFactorizationOK() == false)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector* X = Problem_->GetLHS() ;
  const Epetra_MultiVector* B = Problem_->GetRHS() ;

  if (X == 0)
    AMESOS_CHK_ERR(-1); // something wrong with user's input
  if (B == 0)
    AMESOS_CHK_ERR(-1); // something wrong with user's input

  if (X->NumVectors() != B->NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong with user's input

  // Timing are set inside each Solve().
  if (NumProc() == 1) {
    AMESOS_CHK_ERR(SolveSerial(*X,*B));
  }
  else {
    AMESOS_CHK_ERR(SolveDistributed(*X,*B));
  }

  return(0);
}

//=============================================================================
int Amesos_Lapack::SolveSerial(Epetra_MultiVector& X,
			       const Epetra_MultiVector& B) 
{
  assert(NumProc() == 1);
  Time_->ResetStartTime();
  
  int NumVectors = X.NumVectors();

  Epetra_SerialDenseMatrix DenseX(NumGlobalRows(),NumVectors);
  Epetra_SerialDenseMatrix DenseB(NumGlobalRows(),NumVectors);

  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
      DenseB(i,j) = B[j][i];

  DenseSolver_.SetVectors(DenseX,DenseB);
  (DenseSolver_.Solve());

  for (int i = 0 ; i < NumMyRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
       X[j][i] = DenseX(i,j);

  SolTime_ += Time_->ElapsedTime();
  ++NumSolve_;

  return(0) ;
}

//=============================================================================
int Amesos_Lapack::SolveDistributed(Epetra_MultiVector& X,
				    const Epetra_MultiVector& B) 
{
  assert(NumProc() != 1);
  Time_->ResetStartTime();
  
  int NumVectors = X.NumVectors();

  // vector that contains RHS, overwritten by solution,
  // with all elements on on process 0.
  Epetra_MultiVector SerialVector(SerialMap(),NumVectors);
  // import off-process data
  AMESOS_CHK_ERR(SerialVector.Import(B,RowImporter(),Insert));

  VecTime_ += Time_->ElapsedTime();

  if (MyPID() == 0) {
    Epetra_SerialDenseMatrix DenseX(NumGlobalRows(),NumVectors);
    Epetra_SerialDenseMatrix DenseB(NumGlobalRows(),NumVectors);

    for (int i = 0 ; i < NumGlobalRows() ; ++i)
      for (int j = 0 ; j < NumVectors ; ++j)
	DenseB(i,j) = SerialVector[j][i];

    DenseSolver_.SetVectors(DenseX,DenseB);
    (DenseSolver_.Solve());

    for (int i = 0 ; i < NumGlobalRows() ; ++i)
      for (int j = 0 ; j < NumVectors ; ++j)
	SerialVector[j][i] = DenseX(i,j);
  }

  AMESOS_CHK_ERR(X.Export(SerialVector,RowImporter(),Insert));

  SolTime_ += Time_->ElapsedTime();
  ++NumSolve_;

  return(0) ;
}

//=============================================================================
int Amesos_Lapack::SerialToDense()
{

  assert (NumProc() == 1);
  Time_->ResetStartTime();

  if (Matrix() == 0)
    AMESOS_CHK_ERR(-1); // something is wrong in the user's input

  // zero out matrix. This may be needed if the user
  // has changed the values of the matrix.
  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumGlobalRows() ; ++j)
      DenseMatrix_(i,j) = 0.0;

  // allocate storage to extract matrix rows.
  int Length = Matrix()->MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  int count = 0;
  for (int j = 0 ; j < Matrix()->NumMyRows() ; ++j) {

    int NumEntries;
    int ierr = Matrix()->ExtractMyRowCopy(j, Length, NumEntries,
					&Values[0], &Indices[0]);
    AMESOS_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {
      if ( abs(Values[k]) >= Threshold_)
	if (UseTranspose())
	  DenseMatrix_(Indices[k],j) = Values[k];
	else
	  DenseMatrix_(j,Indices[k]) = Values[k];
      if (Indices[k] == j)
	DenseMatrix_(j,j) = Values[k] + AddToDiag_;
    }
  }

  ConTime_ += Time_->ElapsedTime();
  return 0;
}

//=============================================================================
int Amesos_Lapack::DistributedToDense()
{

  assert(NumProc() != 1);
  Time_->ResetStartTime();

  if (SerialMatrix_) {
    delete SerialMatrix_;
    SerialMatrix_ = 0;
  }

  // check whether allocated 
  if (!RowImporter().SourceMap().SameAs(Matrix()->RowMatrixRowMap())) {
    delete SerialMap_;
    SerialMap_ = 0;
    delete RowImporter_;
    RowImporter_ = 0;
  }

  AMESOS_CHK_ERR(SerialMatrix().Import(*Matrix(),RowImporter(),Insert));
  AMESOS_CHK_ERR(SerialMatrix().FillComplete());

  MatTime_ += Time_->ElapsedTime();

  if (MyPID())
    return(0);

  // zero out matrix. This may be needed if the user
  // has changed the values of the matrix.
  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumGlobalRows() ; ++j)
      DenseMatrix_(i,j) = 0.0;

  // allocate storage to extract matrix rows.
  int Length = SerialMatrix().MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  int count = 0;
  for (int j = 0 ; j < SerialMatrix().NumMyRows() ; ++j) {

    int NumEntries;
    int ierr = SerialMatrix().ExtractMyRowCopy(j, Length, NumEntries,
					       &Values[0], &Indices[0]);
    AMESOS_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {
      int col = Indices[k];
      if ( abs(Values[k]) >= Threshold_) {
	if (UseTranspose())
	  DenseMatrix_(col,j) = Values[k];
	else
	  DenseMatrix_(j,col) = Values[k];
      }
      if (col == j)
	DenseMatrix_(j,j) = Values[k] + AddToDiag_;
    }
  }

  ConTime_ += Time_->ElapsedTime();
  return 0;
} 

// ================================================ ====== ==== ==== == =
void Amesos_Lapack::PrintStatus()
{
  if (MyPID() != 0) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Lapack : Matrix has " << NumGlobalRows() << " rows"
       << " and " << Matrix()->NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Lapack : Nonzero elements per row = "
       << 1.0 * Matrix()->NumGlobalNonzeros() / NumGlobalRows() << endl;
  cout << "Amesos_Lapack : Percentage of nonzero elements = "
       << 100.0 * Matrix()->NumGlobalNonzeros() /
          (NumGlobalRows() * NumGlobalRows()) << endl;
  cout << "Amesos_Lapack : Use transpose = " << UseTranspose_ << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;

}

// ================================================ ====== ==== ==== == =
void Amesos_Lapack::PrintTiming()
{
  if (MyPID() == 0) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Lapack : Time to convert matrix to KLU format = "
       << ConTime_ << " (s)" << endl;
  cout << "Amesos_Lapack : Time to redistribute matrix = "
       << MatTime_ << " (s)" << endl;
  cout << "Amesos_Lapack : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Lapack : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Lapack : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Lapack : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Lapack : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Lapack : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Lapack : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
}
