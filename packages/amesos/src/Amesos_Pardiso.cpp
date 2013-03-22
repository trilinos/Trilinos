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

#ifdef HAVE_AMESOS_PARDISO_MKL

#include "mkl_pardiso.h"
#define F77_PARDISO PARDISO

#else

#define F77_PARDISOINIT F77_FUNC(pardisoinit, PARDISOINIT)
#define F77_PARDISO F77_FUNC(pardiso, PARDISO)

/* PARDISO prototype. */
extern "C" int F77_PARDISOINIT
    (void *, int *, int *, int *, double *, int *);

extern "C" int F77_PARDISO
    (void *, int *, int *, int *, int *, int *, 
     double *, int *, int *, int *, int *, int *, 
     int *, double *, double *, int *, double *);
#endif

#define IPARM(i) iparm_[i-1]

using namespace Teuchos;

//=============================================================================
Amesos_Pardiso::Amesos_Pardiso(const Epetra_LinearProblem &prob) :
  UseTranspose_(false),
  Problem_(&prob),
  MtxConvTime_(-1),
  MtxRedistTime_(-1),
  VecRedistTime_(-1),
  SymFactTime_(-1),
  NumFactTime_(-1),
  SolveTime_(-1),
  maxfct_(1),
  mnum_(1),
  msglvl_(0),
  nrhs_(1),
  pardiso_initialized_(false)
{
  for (int i = 0; i < 64; i++) {
    iparm_[i] = 0;
  }
  iparm_[0] = 1; /* No solver default */
  iparm_[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm_[2] = 1;
  iparm_[3] = 0; /* No iterative-direct algorithm */
  iparm_[4] = 0; /* No user fill-in reducing permutation */
  iparm_[5] = 0; /* Write solution into x */
  iparm_[6] = 0; /* Not in use */
  iparm_[7] = 0; /* Max numbers of iterative refinement steps */
  iparm_[8] = 0; /* Not in use */
  iparm_[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm_[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm_[11] = 0; /* Normal solve (0), or a transpose solve (1) */
  iparm_[12] = 0; /* Do not use (non-)symmetric matchings */
  iparm_[13] = 0; /* Output: Number of perturbed pivots */
  iparm_[14] = 0; /* Peak memory in KB during analysis */
  iparm_[15] = 0; /* Permanent mem in KB from anal. that is used in ph 2&3 */
  iparm_[16] = 0; /* Peak double prec  mem in KB incl one LU factor */
  iparm_[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm_[18] = -1; /* Output: Mflops for LU factorization */
  iparm_[19] = 0; /* Output: Numbers of CG Iterations */

  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++) {
    pt_[i] = 0;
  }
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
#ifdef HAVE_AMESOS_PARDISO_MKL
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                &n, &ddum, &ia_[0], &ja_[0], &idum, &nrhs_,
                iparm_, &msglvl_, &ddum, &ddum, &error);
#else
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                &n, &ddum, &ia_[0], &ja_[0], &idum, &nrhs_,
                iparm_, &msglvl_, &ddum, &ddum, &error, dparm_);
#endif
  }

  AMESOS_CHK_ERRV(CheckError(error));
  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
}

//=============================================================================
int Amesos_Pardiso::ConvertToSerial() 
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
int Amesos_Pardiso::ConvertToPardiso()
{
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    ia_.resize(SerialMatrix().NumMyRows()+1);
    ja_.resize(SerialMatrix().NumMyNonzeros());
    aa_.resize(SerialMatrix().NumMyNonzeros());

    int MaxNumEntries = SerialMatrix().MaxNumEntries();
    std::vector<int>    Indices(MaxNumEntries);
    std::vector<double> Values(MaxNumEntries);

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

  MtxConvTime_ = AddTime("Total matrix conversion time", MtxConvTime_);

  return 0;
}

//=============================================================================
int Amesos_Pardiso::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  // We fill iparm_ using named parameters
  if (ParameterList.isSublist("Pardiso"))
  {
	param_ = ParameterList.sublist("Pardiso");
  }

  msglvl_ = param_.get<int>("Message level", static_cast<int>(debug_));
  iparm_[0] = param_.get<int>("No default parameters", 1);
  iparm_[1] = param_.get<int>("Use METIS reordering" , 2);
  iparm_[2] = param_.get<int>("Number of processors", 1);
  iparm_[3] = param_.get<int>("Do preconditioned CGS iterations", 0);
  iparm_[4] = param_.get<int>("Use user permutation", 0);
  iparm_[5] = param_.get<int>("Solution on X/B", 0);
  iparm_[7] = param_.get<int>("Max num of iterative refinement steps", 0);
  iparm_[9] = param_.get<int>("Perturbation for pivot elements 10^-k", 13);
  iparm_[10] = param_.get<int>("Use (non-)symmetric scaling vectors", 1);
  iparm_[11] = param_.get<int>("Solve transposed", 0);
  iparm_[12] = param_.get<int>("Use (non-)symmetric matchings", 0);
  iparm_[17] = param_.get<int>("Number of non-zeros in LU; -1 to compute", -1);
  iparm_[18] = param_.get<int>("Mflops for LU fact; -1 to compute", -1);

  return 0;
}

//=============================================================================
int Amesos_Pardiso::PerformSymbolicFactorization() 
{
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    // at this point only read unsym matrix
    mtype_ = 11; 
    int error = 0;
    int solver = 0;

    // ============================================================== //
    // Setup Pardiso control parameters und initialize the solvers    //
    // internal adress pointers. This is only necessary for the FIRST //
    // call of the PARDISO solver.                                    // 
    // The number of processors is specified by IPARM(2), in the      //
    // Pardiso sublist.                                               //
    // ============================================================== //
#ifndef HAVE_AMESOS_PARDISO_MKL
    F77_PARDISOINIT(pt_,  &mtype_, &solver, iparm_, dparm_, &error);
#endif
    pardiso_initialized_ = true;

    int num_procs = 1;
    char* var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
      sscanf( var, "%d", &num_procs );
    IPARM(3) = num_procs;

    maxfct_ = 1;         /* Maximum number of numerical factorizations.  */
    mnum_   = 1;         /* Which factorization to use. */

    int phase = 11; 
    error = 0;
    int n = SerialMatrix().NumMyRows();
    int idum;
    double ddum;

#ifdef HAVE_AMESOS_PARDISO_MKL
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error);
#else
    F77_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error, dparm_);
#endif
    AMESOS_CHK_ERR(CheckError(error));
  }

  SymFactTime_ = AddTime("Total symbolic factorization time", SymFactTime_);

  return 0;
}

//=============================================================================
int Amesos_Pardiso::PerformNumericFactorization( ) 
{
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    int phase = 22;
    int error;
    int n = SerialMatrix().NumMyRows();
    int idum;
    double ddum;

#ifdef HAVE_AMESOS_PARDISO_MKL
    F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error);
#else
    F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_, &ddum, &ddum, &error, dparm_);
#endif

    AMESOS_CHK_ERR(CheckError(error));
  }

  NumFactTime_ = AddTime("Total numeric factorization time", NumFactTime_);

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

  ResetTimer();

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

    int error;
    int idum = 0;
    int n = SerialMatrix().NumMyRows();
    int phase = 33;

    for (int i = 0 ; i < NumVectors ; ++i)
#ifdef HAVE_AMESOS_PARDISO_MKL
      F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                         &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                         iparm_, &msglvl_, 
                         SerialBValues + i * n,
                         SerialXValues + i * n,
                         &error);
#else
    F77_PARDISO (pt_, &maxfct_, &mnum_, &mtype_, &phase,
                       &n, &aa_[0], &ia_[0], &ja_[0], &idum, &nrhs_,
                       iparm_, &msglvl_,
                       SerialBValues + i * n,
                       SerialXValues + i * n,
                       &error, dparm_);
#endif
    AMESOS_CHK_ERR(CheckError(error));
  }

  SolveTime_ = AddTime("Total solve time", SolveTime_);

  //  Copy X back to the original vector

  ResetTimer();

  if (Comm().NumProc() != 1) 
  {
    X->Export(*SerialX, Importer(), Insert);
    delete SerialB;
    delete SerialX;
  } // otherwise we are already in place.

  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

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

  std::string p = "Amesos_Pardiso : ";
  PrintLine();

  int n = Matrix().NumGlobalRows();
  int nnz = Matrix().NumGlobalNonzeros();

  std::cout << p << "Matrix has " << n << " rows"
       << " and " << nnz << " nonzeros" << std::endl;
  std::cout << p << "Nonzero elements per row       = "
       << 1.0 *  nnz / n << std::endl;
  std::cout << p << "Percentage of nonzero elements = "
       << 100.0 * nnz /(pow(n,2.0)) << std::endl;
  std::cout << p << "Use transpose                  = " << UseTranspose_ << std::endl;
  std::cout << p << "Number of performed iterative ref. steps = " << IPARM(9) << std::endl;
  std::cout << p << "Peak memory symbolic factorization       = " << IPARM(15) << std::endl;
  std::cout << p << "Permanent memory symbolic factorization  = " << IPARM(16) << std::endl;
  std::cout << p << "Memory numerical fact. and solution      = " << IPARM(17) << std::endl;
  std::cout << p << "Number of nonzeros in factors            = " << IPARM(18) << std::endl;
  std::cout << p << "MFlops of factorization                  = " << IPARM(19) << std::endl;
  std::cout << p << "CG/CGS diagnostic                        = " << IPARM(20) << std::endl;
  std::cout << p << "Inertia: Number of positive eigenvalues  = " << IPARM(22) << std::endl;
  std::cout << p << "Inertia: Number of negative eigenvalues  = " << IPARM(23) << std::endl;

  PrintLine();

  return;
}

// ====================================================================== 
void Amesos_Pardiso::PrintTiming() const
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

  std::string p = "Amesos_Pardiso : ";
  PrintLine();

  std::cout << p << "Time to convert matrix to Pardiso format = "
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

// ====================================================================== 
int Amesos_Pardiso::CheckError(const int error) const
{
  if (!error)
    return 0;
  
  std::cerr << "Amesos: PARDISO returned error code " << error << std::endl;
  std::cerr << "Amesos: Related message from manual is:" << std::endl;

  switch(error)
  {
  case -1:
    std::cerr << "Input inconsistent" << std::endl;
    break;
  case -2:
    std::cerr << "Not enough memory" << std::endl;
    break;
  case -3:
    std::cerr << "Reordering problems" << std::endl;
    break;
  case -4:
    std::cerr << "Zero pivot, numerical fact. or iterative refinement problem. " << std::endl;
    break;
  case -5:
    std::cerr << "Unclassified (internal) error" << std::endl;
    break;
  case -6:
    std::cerr << "Preordering failed (matrix types 11, 13 only)" << std::endl;
    break;
  case -7:
    std::cerr << "Diagonal matrix problem." << std::endl;
    break;
  }

  AMESOS_RETURN(error);
}
