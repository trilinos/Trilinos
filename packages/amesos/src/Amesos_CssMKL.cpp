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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Amesos_CssMKL.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"

#include "mkl_pardiso.h"
#include "mkl_cluster_sparse_solver.h"
#define IPARM(i) iparm_[i-1]

using namespace Teuchos;

//=============================================================================
Amesos_CssMKL::Amesos_CssMKL(const Epetra_LinearProblem &prob) :
  UseTranspose_(false),
  Problem_(&prob),
  MtxConvTime_(-1),
  MtxRedistTime_(-1),
  VecRedistTime_(-1),
  SymFactTime_(-1),
  NumFactTime_(-1),
  SolveTime_(-1),
  mtype_(11), 
  maxfct_(1),
  mnum_(1),
  msglvl_(0),
  nrhs_(1)
{
  for( int i = 0; i < 64; ++i ){
    pt_[i] = nullptr;
    iparm_[i] = 0;
  }
  iparm_[0] = 1; /* No solver default */
  // Reset some of the default parameters
  iparm_[1] = 10;  /* 2: Fill-in reordering from METIS, 3: thread dissection, 10: MPI version of the nested dissection and symbolic factorization*/
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm_[7] = 0;   /* Max numbers of iterative refinement steps */
  iparm_[9] = 13;  /* Perturb the pivot elements with 1E-13 */
  iparm_[10] = 0;  /* Disable nonsymmetric permutation and scaling MPS */
  iparm_[11] = 0;  /* Normal solve (0), or a transpose solve (1) */
  iparm_[12] = 0;  /* Do not use (non-)symmetric matchings */
  iparm_[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm_[20] = -1; /* Pivoting for symmetric indefinite matrices */
  iparm_[26] = 1;  /* Check input matrix is sorted */
  iparm_[27] = 0;  /* double-precision */
  iparm_[34] = 1;  /* Use zero-based indexing */

  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only       */
  /*    necessary for the FIRST call of the CSS solver.                   */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++) {
    pt_[i] = 0;
  }
}

//=============================================================================
Amesos_CssMKL::~Amesos_CssMKL() 
{
  /*
   * Free any memory allocated by the CssMKL library functions
   */
  int error = 0;
  if (IsSymbolicFactorizationOK_)
  {
    int phase = -1;         // release all internal solver memory
    int n = Matrix_->NumGlobalRows();
    void *bdummy, *xdummy;
    const MPI_Fint CssComm = CssComm_;
    cluster_sparse_solver( pt_, const_cast<int*>(&maxfct_),
                           const_cast<int*>(&mnum_), &mtype_, &phase, &n,
                           &aa_[0], &ia_[0], &ja_[0],
                           &perm_[0], &nrhs_, iparm_,
                           const_cast<int*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
    IsSymbolicFactorizationOK_ = false;
    IsNumericFactorizationOK_ = false;
  }

  AMESOS_CHK_ERRV(CheckError(error));
  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
}

//=============================================================================
int Amesos_CssMKL::ConvertToCssMKL()
{
  ResetTimer();

  // =========================================================== //
  // Load Matrix from Problem                                    //
  // =========================================================== //
  OriginalMatrix_ = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  if  ( Reindex_ ) {
#ifdef HAVE_AMESOS_EPETRAEXT
    auto CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(Problem_->GetOperator());
    if(!CrsMatrix) {
      std::cerr << "Amesos_CssMKL requires CsrMatrix to reindex matrices." << std::endl;
      AMESOS_CHK_ERR(-8);
    }
    const Epetra_Map& OriginalMap = CrsMatrix->RowMap();
    StdIndex_ = rcp( new Amesos_StandardIndex( OriginalMap  ) );
    Matrix_ = StdIndex_->StandardizeIndex( CrsMatrix );
    if(!Matrix_) {
      std::cerr << "Amesos_CssMKL reindexing failed" << std::endl;
      AMESOS_CHK_ERR(-8);
    }
#else
    std::cerr << "Amesos_CssMKL requires EpetraExt to reindex matrices." << std::endl;
    AMESOS_CHK_ERR(-8);
#endif
  } else {
    Matrix_ = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  }

  // =========================================================== //
  // Read Matrix into Crs (with Linear/Uniform index)            //
  // =========================================================== //
  int ierr;
  int MaxNumEntries_ = Matrix_->MaxNumEntries();
  int NumMyElements = Matrix_->NumMyRows(); 
  int nnz_loc = Matrix_->NumMyNonzeros();
  ia_.resize( NumMyElements+1 );
  ja_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ); 
  aa_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ); 

  std::vector<int> ColIndicesV_(MaxNumEntries_);
  std::vector<double> RowValuesV_(MaxNumEntries_);

  int *Global_Columns_ = Matrix_->RowMatrixColMap().MyGlobalElements();

  // for sorting column indexes
  typedef std::pair<int, double> Data;
  std::vector<Data> sort_array(MaxNumEntries_);

  // NOTE : Assumes input is distributed in contiguous 1D Block
  int NzThisRow;
  int Ai_index = 0;
  for (int MyRow = 0; MyRow < NumMyElements ; MyRow++) 
  {
    ierr = Matrix_->ExtractMyRowCopy(MyRow, MaxNumEntries_, NzThisRow,
                                            &RowValuesV_[0], &ColIndicesV_[0]);
    AMESOS_CHK_ERR(ierr);

    double *RowValues =  &RowValuesV_[0];
    int    *ColIndices = &ColIndicesV_[0];

    // sort column indexes (default std::sort by first elements of pair)
    for ( int j = 0; j < NzThisRow; j++ ) { 
      sort_array[j].first = Global_Columns_[ColIndices[j]];
      sort_array[j].second = RowValues[j];
    }
    std::sort(&sort_array[0], &sort_array[NzThisRow]);

    ia_[MyRow] = Ai_index ; 
    for ( int j = 0; j < NzThisRow; j++ ) { 
      ja_[Ai_index] = sort_array[j].first;  //Global_Columns_[ColIndices[j]]; 
      aa_[Ai_index] = sort_array[j].second; //RowValues[j];
      Ai_index++;
    }
  }
  ia_[ NumMyElements ] = Ai_index; 
  assert( NumMyElements == MyRow );

  MtxConvTime_ = AddTime("Total matrix conversion time", MtxConvTime_);

  return 0;
}

//=============================================================================
int Amesos_CssMKL::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  // We fill iparm_ using named parameters
  if (ParameterList.isSublist("CssMKL"))
  {
    param_ = ParameterList.sublist("CssMKL");
  }

  msglvl_ = param_.get<int>("Message level", 0);
  iparm_[0] = param_.get<int>("No default parameters", 1);
  iparm_[1] = param_.get<int>("Use METIS reordering" , 10);
  iparm_[7] = param_.get<int>("Max num of iterative refinement steps", 0);
  iparm_[9] = param_.get<int>("Perturbation for pivot elements 10^-k", 13);
  iparm_[10] = param_.get<int>("Use (non-)symmetric scaling vectors", 0);
  iparm_[11] = param_.get<int>("Solve transposed", 0);
  iparm_[12] = param_.get<int>("Use (non-)symmetric matchings", 0);
  iparm_[17] = param_.get<int>("Number of non-zeros in LU; -1 to compute", -1);
  iparm_[20] = param_.get<int>("Pivot for symmetric indefinite matrix", -1);

  return 0;
}

//=============================================================================
int Amesos_CssMKL::PerformSymbolicFactorization() 
{
  ResetTimer();

  int error = 0;
  {
    // Get communicator
    const Epetra_MpiComm& EMpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm());
    const MPI_Comm CssEComm = EMpiComm.Comm();
    CssComm_ = MPI_Comm_c2f(CssEComm);

    // Set matrix 1D block distribution
    auto rangeMap = Matrix_->RowMatrixRowMap();
    int n = Matrix_->NumGlobalRows();
    int first_row = rangeMap.MinMyGID();
    int last_row  = rangeMap.MaxMyGID();
    iparm_[39] = 2;  /* Matrix input format. */
    iparm_[40] = first_row;  /* > Beginning of input domain. */
    iparm_[41] = last_row;   /* > End of input domain. */

    // Allocate perm
    perm_.resize(n);

    // Perform Symbolic Analysis
    int phase = 11; // Analysis
    void *bdummy, *xdummy;
    const MPI_Fint CssComm = CssComm_;

    cluster_sparse_solver( pt_, const_cast<int*>(&maxfct_),
                           const_cast<int*>(&mnum_), &mtype_, &phase, &n,
                           &aa_[0], &ia_[0], &ja_[0],
                           &perm_[0], &nrhs_, iparm_,
                           const_cast<int*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
  }
  SymFactTime_ = AddTime("Total symbolic factorization time", SymFactTime_);

  AMESOS_CHK_ERR(CheckError(error));
  return error;
}

//=============================================================================
int Amesos_CssMKL::PerformNumericFactorization( ) 
{
  ResetTimer();

  int error = 0;
  {
    //int phase = 12; // Analysis, numerical factorization
    int phase = 22; // Analysis, numerical factorization
    int n = Matrix_->NumGlobalRows();
    void *bdummy, *xdummy;
    const MPI_Fint CssComm = CssComm_;
    cluster_sparse_solver( pt_, const_cast<int*>(&maxfct_),
                           const_cast<int*>(&mnum_), &mtype_, &phase, &n,
                           &aa_[0], &ia_[0], &ja_[0],
                           &perm_[0], &nrhs_, iparm_,
                           const_cast<int*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
  }
  NumFactTime_ = AddTime("Total numeric factorization time", NumFactTime_);
 
  // Any failure should be sent to all other processors.
  AMESOS_CHK_ERR(CheckError(error));
  return error;
}

//=============================================================================
bool Amesos_CssMKL::MatrixShapeOK() const 
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
int Amesos_CssMKL::SymbolicFactorization() 
{
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  CreateTimer(Comm());

  ++NumSymbolicFact_;

  // =========================================================== //
  // Convert the matrix into CSR format,                         //
  // =========================================================== //
  ConvertToCssMKL();

  // =========================================================== //
  // Perform Symbolic                                            //
  // =========================================================== //
  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_CssMKL::NumericFactorization() 
{
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  ++NumNumericFact_;

  // =========================================================== //
  // (Re)Convert the matrix into CSR format,                     //
  // =========================================================== //
  // FIXME: this must be checked, now all the matrix is shipped twice here
  ConvertToCssMKL();

  // =========================================================== //
  // Perform Numeric                                             //
  // =========================================================== //
  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_CssMKL::Solve() 
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

  ResetTimer();
  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

  ResetTimer();
  int error = 0;
  {
    double* XValues;
    double* BValues;
    int LDA;

    // FIXME: check LDA
    AMESOS_CHK_ERR(X->ExtractView(&XValues,&LDA));
    AMESOS_CHK_ERR(B->ExtractView(&BValues,&LDA));

    int idum = 0;
    int n = Matrix().NumGlobalRows();
    int phase = 33;

    const MPI_Fint CssComm = CssComm_;

    nrhs_ = NumVectors;
    cluster_sparse_solver( pt_, const_cast<int*>(&maxfct_),
                           const_cast<int*>(&mnum_), &mtype_, &phase, &n,
                           &aa_[0], &ia_[0], &ja_[0],
                           &perm_[0], &nrhs_, iparm_,
                           const_cast<int*>(&msglvl_), 
                           BValues,
                           XValues,
                           &CssComm, &error );
  }
  SolveTime_ = AddTime("Total solve time", SolveTime_);

  // Any failure should be sent to all other processors.
  Comm().Broadcast( &error, 1, 0 );
  if ( error ) {
    if (Comm().MyPID() == 0) {
      AMESOS_CHK_ERR(CheckError(error));
    } else
      return error;
  }

  //  Copy X back to the original vector
  ResetTimer();
  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *X, *B, UseTranspose(), "Amesos_CssMKL");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_CssMKL");

  ++NumSolve_;

  return(0) ;
}

// ====================================================================== 
void Amesos_CssMKL::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  std::string p = "Amesos_CssMKL : ";
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
void Amesos_CssMKL::PrintTiming() const
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

  std::string p = "Amesos_CssMKL : ";
  PrintLine();

  std::cout << p << "Time to convert matrix to CssMKL format = "
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
int Amesos_CssMKL::CheckError(const int error) const
{
  if (!error)
    return 0;
  
  std::cerr << "Amesos: CSS returned error code " << error << std::endl;
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
