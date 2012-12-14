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

#ifdef HAVE_AMESOS_CSPARSE
#include "Amesos_CSparse.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"


using namespace Teuchos;

//=============================================================================
Amesos_CSparse::Amesos_CSparse(const Epetra_LinearProblem &prob) :
  UseTranspose_(false),
  Problem_(&prob),
  MtxConvTime_(-1),
  MtxRedistTime_(-1),
  VecRedistTime_(-1),
  SymFactTime_(-1),
  NumFactTime_(-1),
  SolveTime_(-1)
{
    // Initialize data structures for CSparse
}

//=============================================================================
Amesos_CSparse::~Amesos_CSparse() 
{
    // Clean up

    //AMESOS_CHK_ERRV(CheckError(error));
    if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
    if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
}

//=============================================================================
int Amesos_CSparse::ConvertToSerial() 
{
    ResetTimer();
  const Epetra_RowMatrix* StdIndexMatrix_ ; 
  Epetra_CrsMatrix* CrsMatrixA_;

    if  ( Reindex_ ) {
        if ( Matrix_ == 0 ) AMESOS_CHK_ERR(-7);
#ifdef HAVE_AMESOS_EPETRAEXT
        const Epetra_Map& OriginalMap = Matrix_->RowMatrixRowMap();

        const Epetra_Map &OriginalDomainMap = 
        UseTranspose()?GetProblem()->GetOperator()->OperatorRangeMap():
        GetProblem()->GetOperator()->OperatorDomainMap();
        const Epetra_Map &OriginalRangeMap = 
        UseTranspose()?GetProblem()->GetOperator()->OperatorDomainMap():
        GetProblem()->GetOperator()->OperatorRangeMap();

        StdIndex_ = rcp( new Amesos_StandardIndex( OriginalMap  ) );
        StdIndexDomain_ = rcp( new Amesos_StandardIndex( OriginalDomainMap  ) );
        StdIndexRange_ = rcp( new Amesos_StandardIndex( OriginalRangeMap  ) );

        CrsMatrixA_ = dynamic_cast<Epetra_CrsMatrix *>(Problem_->GetOperator());
        StdIndexMatrix_ = StdIndex_->StandardizeIndex( CrsMatrixA_ );
#else
        std::cerr << "Amesos_CSparse requires EpetraExt to reindex matrices." << std::endl 
        AMESOS_CHK_ERR(-8);
#endif
    } else { 
        StdIndexMatrix_ = Matrix_ ;
    }

    int NumGlobalRows = Matrix_->NumGlobalRows();

    // create a serial map
    int NumMyRows = 0;
    if (Comm().MyPID() == 0) 
        NumMyRows = NumGlobalRows;

    SerialMap_ = rcp(new Epetra_Map(-1, NumMyRows, 0, Comm()));
    if (SerialMap_.get() == 0)
    AMESOS_CHK_ERR(-1);

    Importer_ = rcp(new Epetra_Import(SerialMap(),StdIndexMatrix_->RowMatrixRowMap()));
    if (Importer_.get() == 0)
    AMESOS_CHK_ERR(-1);

    SerialCrsMatrix_ = rcp(new Epetra_CrsMatrix(Copy, SerialMap(), 0));
    if (SerialCrsMatrix_.get() == 0)
    AMESOS_CHK_ERR(-1);

    AMESOS_CHK_ERR(SerialCrsMatrix().Import(*StdIndexMatrix_, Importer(), Add));

    AMESOS_CHK_ERR(SerialCrsMatrix().FillComplete());

    SerialMatrix_ = rcp(SerialCrsMatrix_.get(), false);

    MtxRedistTime_ = AddTime("Total matrix redistribution time", MtxRedistTime_);

    return 0;
}

//=============================================================================
int Amesos_CSparse::ConvertToCSparse()
{

#ifdef HAVE_AMESOS_CSPARSE

  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    csMatrix.p = (ptrdiff_t *) malloc((SerialMatrix().NumMyRows()+1)*sizeof(ptrdiff_t));
    csMatrix.i = (ptrdiff_t *) malloc(SerialMatrix().NumMyNonzeros()*sizeof(ptrdiff_t));
    csMatrix.x = (double *) malloc(SerialMatrix().NumMyNonzeros()*sizeof(double));
    csMatrix.nzmax = SerialMatrix().NumMyNonzeros();
    csMatrix.m = SerialMatrix().NumMyRows();
    csMatrix.n = SerialMatrix().NumMyRows();
    csMatrix.nz = -1;

    int MaxNumEntries = SerialMatrix().MaxNumEntries();
    std::vector<int>    Indices(MaxNumEntries);
    std::vector<double> Values(MaxNumEntries);

    csMatrix.p[0] = 0;
    int count = 0;

    for (int i = 0 ; i < SerialMatrix().NumMyRows() ; ++i)
    {
      int ierr, NumEntriesThisRow;
      ierr = SerialMatrix().ExtractMyRowCopy(i, MaxNumEntries, 
                                             NumEntriesThisRow, 
                                             &Values[0], &Indices[0]);
      if (ierr < 0)
        AMESOS_CHK_ERR(ierr);

      csMatrix.p[i + 1] = csMatrix.p[i] + NumEntriesThisRow;

      for (int j = 0 ; j < NumEntriesThisRow ; ++j)
      {
        if (Indices[j] == i) 
          Values[j] += AddToDiag_;

        csMatrix.i[count] = Indices[j];
        csMatrix.x[count] = Values[j];
        ++count;
      }
    }
    
    if (count != SerialMatrix().NumMyNonzeros())
      AMESOS_CHK_ERR(-1); // something wrong here*/

    // TODO: Avoid this transpose once we have cs_sptsolve()
    csTranMatrix = cs_transpose(&csMatrix, 1);

    free(csMatrix.p);
    free(csMatrix.i);
    free(csMatrix.x);
  }

  MtxConvTime_ = AddTime("Total matrix conversion time", MtxConvTime_);

  return 0;
#else
  AMESOS_CHK_ERR(-1); // Don't have CSPARSE
  return 1;
#endif
}

//=============================================================================
int Amesos_CSparse::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  // retrive CSparse's specific parameters

  if (ParameterList.isSublist("CSparse")) 
  {
      // set solver specific parameters here.

  }
  
  return 0;
}

//=============================================================================
int Amesos_CSparse::PerformSymbolicFactorization() 
{
#ifdef HAVE_AMESOS_CSPARSE

  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    // ============================================================== //
    // Setup CSparse parameteres and call symbolic factorization.
    // ============================================================== //
    int error = 0;

    // Call Csparse here.
    csSymbolic = cs_sqr(2, csTranMatrix, 0);

    AMESOS_CHK_ERR(CheckError(error));
  }

  SymFactTime_ = AddTime("Total symbolic factorization time", SymFactTime_);

  return 0;
#else
  AMESOS_CHK_ERR(-1); // Don't have CSPARSE
  return 1;
#endif
}

//=============================================================================
int Amesos_CSparse::PerformNumericFactorization( ) 
{
#ifdef HAVE_AMESOS_CSPARSE
  ResetTimer();

  if (Comm().MyPID() == 0) 
  {
    int error = 0;
    // Call Csparse here.
    csNumeric = cs_lu(csTranMatrix, csSymbolic, 1e-15);

    AMESOS_CHK_ERR(CheckError(error));
  }

  NumFactTime_ = AddTime("Total numeric factorization time", NumFactTime_);

  return 0;
#else
  AMESOS_CHK_ERR(-1); // Don't have CSPARSE
  return 1;
#endif
}

//=============================================================================
bool Amesos_CSparse::MatrixShapeOK() const 
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
int Amesos_CSparse::SymbolicFactorization() 
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
  // as required by CSparse.                                     //
  // =========================================================== //

  ConvertToCSparse();

  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_CSparse::NumericFactorization() 
{
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  ++NumNumericFact_;

  // FIXME: this must be checked, now all the matrix is shipped twice here
  ConvertToSerial();
  ConvertToCSparse();

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;

  return(0);
}

//=============================================================================
int Amesos_CSparse::Solve() 
{
#ifdef HAVE_AMESOS_CSPARSE
  Epetra_MultiVector* vecX = 0 ;
  Epetra_MultiVector* vecB = 0 ;

#ifdef HAVE_AMESOS_EPETRAEXT
  Teuchos::RCP<Epetra_MultiVector> vecX_rcp;
  Teuchos::RCP<Epetra_MultiVector> vecB_rcp;
#endif

  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector* X = Problem_->GetLHS();
  Epetra_MultiVector* B = Problem_->GetRHS();

  if ((X == 0) || (B == 0))
    AMESOS_CHK_ERR(-1); 

  int NumVectors = X->NumVectors();
  if (NumVectors != B->NumVectors())
    AMESOS_CHK_ERR(-1); 

    if ( Reindex_ ) { 
#ifdef HAVE_AMESOS_EPETRAEXT
      vecX_rcp = StdIndexDomain_->StandardizeIndex( *X ) ;
      vecB_rcp = StdIndexRange_->StandardizeIndex( *B ) ;

      vecX = &*vecX_rcp;
      vecB = &*vecB_rcp;
#else
      AMESOS_CHK_ERR( -13 ) ; // Amesos_CSparse can't handle non-standard indexing without EpetraExt 
#endif
    } else {
      vecX = X ;
      vecB = B ;
    } 
    
  // vectors with SerialMap_
  Epetra_MultiVector* SerialB;
  Epetra_MultiVector* SerialX;

  ResetTimer();

  SerialX = new Epetra_MultiVector(SerialMap(),NumVectors);
  SerialB = new Epetra_MultiVector(SerialMap(),NumVectors);

  SerialB->Import(*vecB,Importer(),Insert);

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

    int error = 0;
    int n = SerialMatrix().NumMyRows();

    double *x = (double *) malloc(n * sizeof(double));
    // TODO: OMP here ??
    for (int i = 0 ; i < NumVectors ; ++i)
    {
        // Call Csparse here
        cs_ipvec(csNumeric->pinv, SerialBValues+i*n, x, n);
        cs_lsolve(csNumeric->L, x);
        cs_usolve(csNumeric->U, x);
        cs_ipvec(csSymbolic->q, x, SerialXValues+i*n, n);
    }
    free(x);

    AMESOS_CHK_ERR(CheckError(error));
  }

  SolveTime_ = AddTime("Total solve time", SolveTime_);

  //  Copy X back to the original vector

  ResetTimer();

  vecX->Export(*SerialX, Importer(), Insert);
  delete SerialB;
  delete SerialX;

  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_);

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *X, *B, UseTranspose(), "Amesos_CSparse");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*X, *B, "Amesos_CSparse");

  ++NumSolve_;

  return(0) ;
#else
  AMESOS_CHK_ERR(-1); // Don't have CSPARSE
  return 1;
#endif
}

// ====================================================================== 
void Amesos_CSparse::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  std::string p = "Amesos_CSparse : ";
  PrintLine();

  PrintLine();

  return;
}

// ====================================================================== 
void Amesos_CSparse::PrintTiming() const
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

  std::string p = "Amesos_CSparse : ";
  PrintLine();

  std::cout << p << "Time to convert matrix to Csparse format = "
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
int Amesos_CSparse::CheckError(const int error) const
{
  if (!error)
    return 0;
  
  std::cerr << "Amesos: CSparse returned error code " << error << std::endl;

  AMESOS_RETURN(error);
}

#endif
