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

#define IPARM(I) iparm_[(I) - 1]

using namespace Teuchos;

//=============================================================================
Amesos_Pardiso::Amesos_Pardiso(const Epetra_LinearProblem &prob) :
  UseTranspose_(false),
  Problem_(&prob),
  maxfct_(1),
  mnum_(1),
  msglvl_(0),
  nrhs_(1),
  pardiso_initialized_(false)
{
  // Sets all parameters (many unused) to zero
  for (int i = 1 ; i < 64 ; ++i)
    IPARM(i) = 0; 

  // setting parameters from manual's default
  IPARM(1) = 0; // use default values
  IPARM(2) = 2; // Fill-in reduction reordering
  IPARM(3) = 1; // Number of processors
  IPARM(4) = 0; // Preconditioned CGS
  IPARM(6) = 0; // write solution on X
  IPARM(8) = 0; // number of iterative refinement steps
  IPARM(10) = 8; // pivot perturbation
  IPARM(11) = 1; // MPS scaling of the unsymmetric reordering
  IPARM(18) = -1; // number of nonzeros in factor
  IPARM(19) = 0; // MFlops of factorization
  IPARM(21) = 1; // pivoting for undefinite symmetric matrices
}

//=============================================================================
Amesos_Pardiso::~Amesos_Pardiso() 
{
  int phase = -1;                 /* Release internal memory. */
  int error = 0;
  int idum;
  double ddum;

  if (pardiso_initialized_ ) {
    int n = SerialMatrix().NumMyRows();
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                &n, &ddum, &ia_[0], &ja_[0], &idum, &nrhs_,
                iparm_, &msglvl_, &ddum, &ddum, &error);
  }

  AMESOS_CHK_ERRV(CheckError(error));
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
    ia_.resize(SerialMatrix().NumMyRows()+1);
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
int Amesos_Pardiso::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  // retrive PARDISO's specific parameters

  if (ParameterList.isSublist("pardiso")) 
  {
    const Teuchos::ParameterList& PardisoList = ParameterList.sublist("Pardiso");

    if (PardisoList.isParameter("MSGLVL"))
      msglvl_ = PardisoList.get<int>("MSGLVL");
    else
      if ( debug_ ) msglvl_ = 1 ; //  msglvl prints statistical information, but is the closest 
    //  thing I found to debug print statements - KSS

    if (PardisoList.isParameter("IPARM(1)"))
      IPARM(1) = PardisoList.get<int>("IPARM(1)");

    if (PardisoList.isParameter("IPARM(2)"))
      IPARM(2) = PardisoList.get<int>("IPARM(2)");

    if (PardisoList.isParameter("IPARM(3)"))
      IPARM(3) = PardisoList.get<int>("IPARM(3)");

    if (PardisoList.isParameter("IPARM(4)"))
      IPARM(4) = PardisoList.get<int>("IPARM(4)");

    if (PardisoList.isParameter("IPARM(8)"))
      IPARM(8) = PardisoList.get<int>("IPARM(8)");

    if (PardisoList.isParameter("IPARM(10)"))
      IPARM(10) = PardisoList.get<int>("IPARM(10)");

    if (PardisoList.isParameter("IPARM(11)"))
      IPARM(11) = PardisoList.get<int>("IPARM(11)");

    if (PardisoList.isParameter("IPARM(18)"))
      IPARM(18) = PardisoList.get<int>("IPARM(18)");

    if (PardisoList.isParameter("IPARM(19)"))
      IPARM(19) = PardisoList.get<int>("IPARM(19)");

    if (PardisoList.isParameter("IPARM(21)"))
      IPARM(21) = PardisoList.get<int>("IPARM(21)");
  }
  
  return 0;
}

//=============================================================================
int Amesos_Pardiso::PerformSymbolicFactorization() 
{
  ResetTime();

  if (Comm().MyPID() == 0) 
  {
    // at this point only read unsym matrix
    mtype_ = 11; 

    // ============================================================== //
    // Setup Pardiso control parameters und initialize the solvers    //
    // internal adress pointers. This is only necessary for the FIRST //
    // call of the PARDISO solver.                                    // 
    // The number of processors is specified by IPARM(2), in the      //
    // Pardiso sublist.                                               //
    // ============================================================== //

    F77_PARDISOINIT(pt_,  &mtype_, iparm_);
    pardiso_initialized_ = true; 
    /*
    char* var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
      sscanf( var, "%d", &num_procs );
    else {
      cerr << "Please set the environment OMP_NUM_THREADS to either" << endl;
      cerr << "1 or the number of OMP processes you want to use" << endl;
      AMESOS_CHK_ERR(-1);
    }

    iparm_[2]  = num_procs;
    */

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

    AMESOS_CHK_ERR(CheckError(error));
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

    AMESOS_CHK_ERR(CheckError(error));
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

    AMESOS_CHK_ERR(CheckError(error));
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
    ComputeTrueResidual(Matrix(), *X, *B, UseTranspose(), "Amesos_Pardiso");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_Pardiso");

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
  cout << p << "Nonzero elements per row       = "
       << 1.0 *  nnz / n << endl;
  cout << p << "Percentage of nonzero elements = "
       << 100.0 * nnz /(pow(n,2.0)) << endl;
  cout << p << "Use transpose                  = " << UseTranspose_ << endl;
  cout << p << "Number of performed iterative ref. steps = " << IPARM(9) << endl;
  cout << p << "Peak memory symbolic factorization       = " << IPARM(15) << endl;
  cout << p << "Permanent memory symbolic factorization  = " << IPARM(16) << endl;
  cout << p << "Memory numerical fact. and solution      = " << IPARM(17) << endl;
  cout << p << "Number of nonzeros in factors            = " << IPARM(18) << endl;
  cout << p << "MFlops of factorization                  = " << IPARM(19) << endl;
  cout << p << "CG/CGS diagnostic                        = " << IPARM(20) << endl;
  cout << p << "Inertia: Number of positive eigenvalues  = " << IPARM(22) << endl;
  cout << p << "Inertia: Number of negative eigenvalues  = " << IPARM(23) << endl;

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

  cout << p << "Time to convert matrix to Pardiso format = "
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

// ====================================================================== 
int Amesos_Pardiso::CheckError(const int error) const
{
  if (!error)
    return 0;
  
  cerr << "Amesos: PARDISO returned error code " << error << endl;
  cerr << "Amesos: Related message from manual is:" << endl;

  switch(error)
  {
  case -1:
    cerr << "Input inconsistent" << endl;
    break;
  case -2:
    cerr << "Not enough memory" << endl;
    break;
  case -3:
    cerr << "Reordering problems" << endl;
    break;
  case -4:
    cerr << "Zero pivot, numerical fact. or iterative refinement problem. " << endl;
    break;
  case -5:
    cerr << "Unclassified (internal) error" << endl;
    break;
  case -6:
    cerr << "Preordering failed (matrix types 11, 13 only)" << endl;
    break;
  case -7:
    cerr << "Diagonal matrix problem." << endl;
    break;
  }

  AMESOS_RETURN(error);
}
