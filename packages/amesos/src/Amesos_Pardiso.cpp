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

#include "Amesos_Pardiso.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"

#define F77_PARDISOINIT F77_FUNC(pardisoinit, PARDISOINIT)
#define F77_PARDISO F77_FUNC(pardiso, PARDISO)

/* PARDISO prototype. */
extern "C" int F77_PARDISOINIT
    (void *, int *, int *);

extern "C" int F77_PARDISO
    (void *, int *, int *, int *, int *, int *, 
     double *, int *, int *, int *, int *, int *, 
     int *, double *, double *, int *);


//=============================================================================
Amesos_Pardiso::Amesos_Pardiso(const Epetra_LinearProblem &prob) :
  UseTranspose_(false),
  Problem_(&prob),
  nrhs_(1)
{
}

//=============================================================================
Amesos_Pardiso::~Amesos_Pardiso(void) 
{
  int n = SerialMatrix().NumMyRows();
  int phase = -1;                 /* Release internal memory. */
  int error = 0;
  int idum;
  double ddum;

  if (Problem_->GetOperator() != 0 && Comm().MyPID() == 0)
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                &n, &ddum, &ia_[0], &ja_[0], &idum, &nrhs_,
                iparm_, &msglvl_, &ddum, &ddum, &error);

  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
}

//=============================================================================
int Amesos_Pardiso::ConvertToSerial() 
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
int Amesos_Pardiso::ConvertToPardiso()
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    ia_.resize(SerialMatrix().NumMyRows());
    ja_.resize(SerialMatrix().NumMyNonzeros());
    aa_.resize(SerialMatrix().NumMyNonzeros());

    int MaxNumEntries = SerialMatrix().MaxNumEntries();
    vector<int>    Indices(MaxNumEntries);
    vector<double> Values(MaxNumEntries);

    // Requires FORTRAN numbering (from 1)
    ia_[0] = 1;
    int count = 0;

    for (int i = 0 ; i < SerialMatrix().NumMyRows() ; ++i)
    {
      int ierr, NumEntriesThisRow;
      ierr = SerialMatrix().ExtractMyRowCopy(i, MaxNumEntries, 
                                             NumEntriesThisRow, 
                                             &Values[0], &Indices[0]);
      if (ierr < 0)
        AMESOS_CHK_ERR(ierr);

      ia_[i + 1] = ia_[i] + NumEntriesThisRow;

      for (int j = 0 ; j < NumEntriesThisRow ; ++j)
      {
        if (Indices[j] == i) 
          Values[j] += AddToDiag_;

        ja_[count] = Indices[j] + 1;
        aa_[count] = Values[j];
        ++count;
      }
    }
    
    if (count != SerialMatrix().NumMyNonzeros())
      AMESOS_CHK_ERR(-1); // something wrong here
  }

  AddTime("matrix conversion");

  return 0;
}

//=============================================================================
int Amesos_Pardiso::SetParameters(Teuchos::ParameterList &ParameterList) 
{
  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));

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
int Amesos_Pardiso::PerformSymbolicFactorization() 
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    int num_procs;

    // FIXME: at this point only read unsym matrix
    mtype_ = 11; 

    // ============================================================== //
    // Setup Pardiso control parameters und initialize the solvers    //
    // internal adress pointers. This is only necessary for the FIRST //
    // call of the PARDISO solver.                                    // 
    // Also get the number of processors (in terms of OMP threads),   //
    // and quit if the user doesn't specify them. This is as done in  //
    // one of the PARDISO examples.                                   //
    // ============================================================== //

    F77_PARDISOINIT(pt_,  &mtype_, iparm_);

    char* var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
      sscanf( var, "%d", &num_procs );
    else {
      cerr << "Please set the environment OMP_NUM_THREADS to either" << endl;
      cerr << "1 or the number of OMP processes you want to use" << endl;
      AMESOS_CHK_ERR(-1);
    }

    iparm_[2]  = num_procs;

    maxfct_ = 1;         /* Maximum number of numerical factorizations.  */
    mnum_   = 1;         /* Which factorization to use. */

    int phase = 11; 
    int error = 0;
    int n = SerialMatrix().NumMyRows();
    int idum;
    double ddum;

    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error);

    if (error != 0) 
    {
      cerr << "Amesos_Pardiso: error during symbolic factorization ("
        << error << ")" << endl;
      AMESOS_CHK_ERR(-1);
    }
  }

  AddTime("symbolic");

  return 0;
}

//=============================================================================
int Amesos_Pardiso::PerformNumericFactorization( ) 
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    int phase = 22;
    int error;
    int n = SerialMatrix().NumMyRows();
    int idum;
    double ddum;

    F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error);

    if (error != 0) 
    {
      cerr << "Amesos_Pardiso: error during symbolic factorization ("
        << error << ")" << endl;
      AMESOS_CHK_ERR(-1);
    }
  }

  AddTime("numeric");

  return 0;
}

//=============================================================================
bool Amesos_Pardiso::MatrixShapeOK() const 
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
int Amesos_Pardiso::SymbolicFactorization() 
{
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
  // as required by PARDISO.                                     //
  // =========================================================== //

  ConvertToPardiso();

  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_Pardiso::NumericFactorization() 
{
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  ++NumNumericFact_;

  // FIXME: this must be checked, now all the matrix is shipped twice here
  ConvertToSerial();
  ConvertToPardiso();

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_Pardiso::Solve() 
{
  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector* X = Problem_->GetLHS();
  Epetra_MultiVector* B = Problem_->GetRHS();

  if ((X == 0) || (B == 0))
    AMESOS_CHK_ERR(-1); 

  int NumVectors = X->NumVectors();
  if (NumVectors != B->NumVectors())
    AMESOS_CHK_ERR(-1); 

  // vectors with SerialMap_
  Epetra_MultiVector* SerialB;
  Epetra_MultiVector* SerialX;

  ResetTime();

  if (Comm().NumProc() == 1) 
  {
    SerialB = B;
    SerialX = X;
  } 
  else 
  {
    SerialX = new Epetra_MultiVector(SerialMap(),NumVectors);
    SerialB = new Epetra_MultiVector(SerialMap(),NumVectors);

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

    int error;
    int idum = 0;
    int n = SerialMatrix().NumMyRows();
    int phase = 33;

    for (int i = 0 ; i < NumVectors ; ++i)
      F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                         &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                         iparm_, &msglvl_, 
                         SerialBValues + i * n,
                         SerialXValues + i * n,
                         &error);

    if (error != 0) 
    {
      cerr << "Amesos_Pardiso: error during solution ("
        << error << ")" << endl;
      AMESOS_CHK_ERR(-1);
    }
  }

  AddTime("solve");

  //  Copy X back to the original vector

  ResetTime();

  if (Comm().NumProc() != 1) 
  {
    X->Export(*SerialX, Importer(), Insert);
    delete SerialB;
    delete SerialX;
  } // otherwise we are already in place.

  AddTime("vector redistribution");

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *X, *B, UseTranspose(), "Amesos_Taucs");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_Taucs");

  ++NumSolve_;

  return(0) ;
}

// ====================================================================== 
void Amesos_Pardiso::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  string p = "Amesos_Pardiso : ";
  PrintLine();

  int n = Matrix().NumGlobalRows();
  int nnz = Matrix().NumGlobalNonzeros();

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

// ====================================================================== 
void Amesos_Pardiso::PrintTiming() const
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

  string p = "Amesos_Pardiso : ";
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
