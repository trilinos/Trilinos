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
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_LAPACK.h"
using namespace Teuchos;

//=============================================================================
Amesos_Lapack::Amesos_Lapack(const Epetra_LinearProblem &Problem) :
  UseTranspose_(false),
  Problem_(&Problem),
//  verbose_(0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0)
{
  Teuchos::ParameterList ParamList;
  SetParameters(ParamList);
  verbose_ = 0;
}

//=============================================================================
Amesos_Lapack::~Amesos_Lapack(void) 
{
  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || (verbose_ == 2)) 
    PrintTiming();
  if ((verbose_ && PrintStatus_) || (verbose_ == 2)) 
    PrintStatus();
}

//=============================================================================
// Default values are defined in the constructor.
int Amesos_Lapack::SetParameters( Teuchos::ParameterList &ParameterList ) 
{
  // retrive general parameters

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", PrintTiming_);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", PrintStatus_);

  // compute norms of solution and right-hand side vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",ComputeVectorNorms_);

  // compute the true residual ||b - Ax|| after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",ComputeTrueResidual_);

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",verbose_);

  // Drop all elements whose abolute value is below Threshold
  if( ParameterList.isParameter("Threshold") )
    Threshold_ = ParameterList.get("Threshold", Threshold_);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", AddToDiag_);
 
  bool Equilibrate = true;
  if (ParameterList.isSublist("Lapack") ) {
    Teuchos::ParameterList& LAPACKParams = ParameterList.sublist("Lapack") ;
    Equilibrate = LAPACKParams.get("Equilibrate", Equilibrate);
  }
  DenseSolver_.FactorWithEquilibration(Equilibrate);

  return(0);
}

//=============================================================================
bool Amesos_Lapack::MatrixShapeOK() const 
{
  if (GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
      GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints()) {
    return(false);
  }
  else
    return(true);
}

//=============================================================================
// Build SerialMap_ and Importer_. SerialMap_ is equal to RowMatrixRowMap()
// for the serial case, otherwise we actually need to build a new map.
// Importer_ is constructed only for the parallel case. All these objects are
// destroyed and rebuilt every time SymbolicFactorization() is called, to
// allow the possibility of completely different matrices given in input.
//
// NOTE: A possible bug is nasted here. If the user destroys the
// RowMatrixRowMap(), Importer() may not work properly any more. Anyway,
// this bug is due to misuse of the class, since it is supposed that the
// structure of A (and therefore all its components) remain untouched after a
// call to SymbolicFactorization().
// ============================================================================
int Amesos_Lapack::SymbolicFactorization() 
{
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  if (GetProblem()->GetMatrix() == 0)
    AMESOS_CHK_ERR(-1); // problem still not properly set

  InitTime(Comm()); // Initialize timer
  ResetTime();

  if (Comm().NumProc() == 1)
    SerialMap_ = rcp(const_cast<Epetra_Map*>(&(Matrix()->RowMatrixRowMap())), 
                     false);
  else
  {
    int NumElements = 0;
    if (Comm().MyPID() == 0)
      NumElements = NumGlobalRows();

    SerialMap_ = rcp(new Epetra_Map(-1, NumElements, 0, Comm()));
    if (SerialMap_.get() == 0)
      AMESOS_CHK_ERR(-1);

    Importer_ = rcp(new Epetra_Import(SerialMap(), Matrix()->RowMatrixRowMap()));
    if (Importer_.get() == 0)
      AMESOS_CHK_ERR(-1);
  }

  AddTime("symbolic");
  IsSymbolicFactorizationOK_ = true;
  ++NumSymbolicFact_;

  return(0);
}

//=============================================================================
int Amesos_Lapack::NumericFactorization() 
{
  IsNumericFactorizationOK_ = false;

  // perform the symbolic (that is, build map and importer) if not done yet.
  if (IsSymbolicFactorizationOK() == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  // Only on processor 0 define the dense matrix.
  if (Comm().MyPID() == 0)
    AMESOS_CHK_ERR(DenseMatrix_.Shape(NumGlobalRows(),NumGlobalRows()));

  AMESOS_CHK_ERR(DistributedToSerial());
  AMESOS_CHK_ERR(SerialToDense());
  AMESOS_CHK_ERR(DenseToFactored());

  IsNumericFactorizationOK_ = true;
  ++NumNumericFact_;
  
  return(0);
}

//=============================================================================
int Amesos_Lapack::Solve() 
{
  if (IsNumericFactorizationOK() == false)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector*       X = Problem_->GetLHS();
  const Epetra_MultiVector* B = Problem_->GetRHS();

  if (X == 0) AMESOS_CHK_ERR(-1); // something wrong with user's input
  if (B == 0) AMESOS_CHK_ERR(-1); // something wrong with user's input
  if (X->NumVectors() != B->NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong with user's input

  // Timing are set inside each Solve().
  int ierr;
  if (Comm().NumProc() == 1)
    ierr = SolveSerial(*X,*B);
  else
    ierr = SolveDistributed(*X,*B);

  AMESOS_RETURN(ierr);
}

//=============================================================================
int Amesos_Lapack::SolveSerial(Epetra_MultiVector& X,
			       const Epetra_MultiVector& B) 
{
  ResetTime();
  
  int NumVectors = X.NumVectors();

  Epetra_SerialDenseMatrix DenseX(NumGlobalRows(),NumVectors);
  Epetra_SerialDenseMatrix DenseB(NumGlobalRows(),NumVectors);

  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
      DenseB(i,j) = B[j][i];

  DenseSolver_.SetVectors(DenseX,DenseB);
  AMESOS_CHK_ERR(DenseSolver_.Solve());

  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
       X[j][i] = DenseX(i,j);

  AddTime("solve");
  ++NumSolve_;

  return(0) ;
}

//=============================================================================
int Amesos_Lapack::SolveDistributed(Epetra_MultiVector& X,
				    const Epetra_MultiVector& B) 
{
  ResetTime();
  
  int NumVectors = X.NumVectors();

  // vector that contains RHS, overwritten by solution,
  // with all elements on on process 0.
  Epetra_MultiVector SerialVector(SerialMap(),NumVectors);
  // import off-process data
  AMESOS_CHK_ERR(SerialVector.Import(B,Importer(),Insert));

  AddTime("vector redistribution");
  ResetTime();

  if (Comm().MyPID() == 0) {
    Epetra_SerialDenseMatrix DenseX(NumGlobalRows(),NumVectors);
    Epetra_SerialDenseMatrix DenseB(NumGlobalRows(),NumVectors);

    for (int i = 0 ; i < NumGlobalRows() ; ++i)
      for (int j = 0 ; j < NumVectors ; ++j)
	DenseB(i,j) = SerialVector[j][i];

    DenseSolver_.SetVectors(DenseX,DenseB);
    AMESOS_CHK_ERR(DenseSolver_.Solve());

    for (int i = 0 ; i < NumGlobalRows() ; ++i)
      for (int j = 0 ; j < NumVectors ; ++j)
	SerialVector[j][i] = DenseX(i,j);
  }

  AddTime("solve");
  ResetTime();

  AMESOS_CHK_ERR(X.Export(SerialVector,Importer(),Insert));

  AddTime("vector redistribution");
  ++NumSolve_;

  return(0) ;
}

//=============================================================================
int Amesos_Lapack::SerialToDense()
{
  if (Comm().MyPID())
    return(0);

  ResetTime();

  for (int i = 0 ; i < NumGlobalRows() ; ++i)
    for (int j = 0 ; j < NumGlobalRows() ; ++j)
      DenseMatrix_(i,j) = 0.0;

  // allocate storage to extract matrix rows.
  int Length = SerialMatrix().MaxNumEntries();
  vector<double> Values(Length);
  vector<int>    Indices(Length);

  for (int j = 0 ; j < SerialMatrix().NumMyRows() ; ++j) 
  {
    int NumEntries;
    int ierr = SerialMatrix().ExtractMyRowCopy(j, Length, NumEntries,
                                               &Values[0], &Indices[0]);
    AMESOS_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {
      if (fabs(Values[k]) >= Threshold_)
	if (UseTranspose())
	  DenseMatrix_(Indices[k],j) = Values[k];
	else
	  DenseMatrix_(j,Indices[k]) = Values[k];
      if (Indices[k] == j)
	DenseMatrix_(j,j) = Values[k] + AddToDiag_;
    }
  }

  AddTime("matrix conversion");
  return(0);
}

//=============================================================================
int Amesos_Lapack::DistributedToSerial()
{
  ResetTime();

  if (Comm().NumProc() == 1)
    SerialMatrix_ = rcp(const_cast<Epetra_RowMatrix*>(Matrix()), false);
  else
  {
    SerialCrsMatrix_ = rcp(new Epetra_CrsMatrix(Copy,SerialMap(), 0));
    SerialMatrix_    = rcp(&SerialCrsMatrix(), false);
    AMESOS_CHK_ERR(SerialCrsMatrix().Import(*Matrix(),Importer(),Insert));
    AMESOS_CHK_ERR(SerialCrsMatrix().FillComplete());
  }

  AddTime("matrix redistribution");

  return(0);
} 

// ====================================================================== 
int Amesos_Lapack::GEEV(Epetra_Vector& Er, Epetra_Vector& Ei)
{
  if (IsSymbolicFactorizationOK() == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  if (Comm().MyPID() == 0)
    AMESOS_CHK_ERR(DenseMatrix_.Shape(NumGlobalRows(),NumGlobalRows()));

  AMESOS_CHK_ERR(DistributedToSerial());
  AMESOS_CHK_ERR(SerialToDense());

  Teuchos::RefCountPtr<Epetra_Vector> LocalEr;
  Teuchos::RefCountPtr<Epetra_Vector> LocalEi;

  if (Comm().NumProc() == 1)
  {
    LocalEr = Teuchos::rcp(&Er, false);
    LocalEi = Teuchos::rcp(&Ei, false);
  }
  else
  {
    LocalEr = Teuchos::rcp(new Epetra_Vector(*SerialMap_));
    LocalEi = Teuchos::rcp(new Epetra_Vector(*SerialMap_));
  }

  if (Comm().MyPID() == 0) 
  {
    int n = NumGlobalRows();
    char jobvl = 'N'; /* V/N to calculate/not calculate left eigenvectors
                         of matrix H.*/
    char jobvr = 'N'; /* As above, but for right eigenvectors. */
    int info = 0;
    int ldvl = n;
    int ldvr = n;

    Er.PutScalar(0.0);
    Ei.PutScalar(0.0);

    Epetra_LAPACK LAPACK;

    vector<double> work(1);
    int lwork = -1;

    LAPACK.GEEV(jobvl, jobvr, n, DenseMatrix_.A(), n, 
                LocalEr->Values(), LocalEi->Values(), NULL,
                ldvl, NULL, 
                ldvr, &work[0], 
                lwork, &info);

    lwork = (int)work[0];
    work.resize(lwork);
    LAPACK.GEEV(jobvl, jobvr, n, DenseMatrix_.A(), n, 
                LocalEr->Values(), LocalEi->Values(), NULL,
                ldvl, NULL, 
                ldvr, &work[0], 
                lwork, &info);

    if (info)
      AMESOS_CHK_ERR(info);
  }

  if (Comm().NumProc() != 1)
  {
    // I am not really sure that exporting the results make sense... 
    // It is just to be coherent with the other parts of the code.
    Er.Export(*LocalEr, Importer(), Insert);
    Ei.Export(*LocalEi, Importer(), Insert);
  }

  return(0);
}

//=============================================================================
int Amesos_Lapack::DenseToFactored() 
{
  ResetTime();

  if (Comm().MyPID() == 0) {
    AMESOS_CHK_ERR(DenseSolver_.SetMatrix(DenseMatrix_));
    AMESOS_CHK_ERR(DenseSolver_.Factor());
  }

  AddTime("numeric");
  return(0);
}

// ================================================ ====== ==== ==== == =
void Amesos_Lapack::PrintStatus()
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  PrintLine();
  string p = "Amesos_Lapack : ";

  int NumGlobalElements = Matrix()->NumGlobalRows();
  int numentries        = Matrix()->NumGlobalNonzeros();

  cout << p << "Matrix has " << NumGlobalElements << " rows"
       << " and " << numentries << " nonzeros" << endl;
  cout << p << "Nonzero elements per row = "
       << 1.0 * numentries / NumGlobalElements << endl;
  cout << p << "Percentage of nonzero elements = "
       << 100.0 * numentries / (pow(NumGlobalElements,2.0)) << endl;
  cout << p << "Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;
}

// ================================================ ====== ==== ==== == =
void Amesos_Lapack::PrintTiming()
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime("matrix conversion");
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

  string p = "Amesos_Lapack : ";
  PrintLine();

  cout << p << "Time to convert matrix to Klu format = "
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
