/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#include "Tifpack_ConfigDefs.hpp"
#include <iomanip>
#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_PointRelaxation.hpp"
#include "Tifpack_Utils.hpp"
#include "Tifpack_Condest.hpp"

static const int TIFPACK_JACOBI = 0;
static const int TIFPACK_GS = 1;
static const int TIFPACK_SGS = 2;

//==============================================================================
Tifpack_PointRelaxation::
Tifpack_PointRelaxation(const Tpetra_RowMatrix* Matrix_in) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  NumSweeps_(1),
  DampingFactor_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  PrecType_(TIFPACK_JACOBI),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Matrix_(Teuchos::rcp(Matrix_in,false)),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false)
{
}

//==============================================================================
int Tifpack_PointRelaxation::SetParameters(Teuchos::ParameterList& List)
{

  string PT;
  if (PrecType_ == TIFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == TIFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == TIFPACK_SGS)
    PT = "symmetric Gauss-Seidel";

  PT = List.get("relaxation: type", PT);

  if (PT == "Jacobi")
    PrecType_ = TIFPACK_JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = TIFPACK_GS;
  else if (PT == "symmetric Gauss-Seidel")
    PrecType_ = TIFPACK_SGS;
  else {
    TIFPACK_CHK_ERR(-2);
  }
  
  NumSweeps_            = List.get("relaxation: sweeps",NumSweeps_);
  DampingFactor_        = List.get("relaxation: damping factor", 
                                   DampingFactor_);
  MinDiagonalValue_     = List.get("relaxation: min diagonal value", 
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("relaxation: zero starting solution", 
                                   ZeroStartingSolution_);

  DoBackwardGS_         = List.get("relaxation: backward mode",DoBackwardGS_);
  
  SetLabel();

  return(0);
}

//==============================================================================
const Tpetra_Comm& Tifpack_PointRelaxation::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Tpetra_Map& Tifpack_PointRelaxation::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Tpetra_Map& Tifpack_PointRelaxation::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Tifpack_PointRelaxation::
Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    TIFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    TIFPACK_CHK_ERR(-2);

  TIFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  return(0);
}

//==============================================================================
int Tifpack_PointRelaxation::Initialize()
{
  IsInitialized_ = false;

  if (Matrix_ == Teuchos::null)
    TIFPACK_CHK_ERR(-2);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Tpetra_Time(Comm()) );

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    TIFPACK_CHK_ERR(-2); // only square matrices

  NumMyRows_ = Matrix_->NumMyRows();
  NumMyNonzeros_ = Matrix_->NumMyNonzeros();
  NumGlobalRows_ = Matrix_->NumGlobalRows();
  NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros();

  if (Comm().NumProc() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Tifpack_PointRelaxation::Compute()
{
  int ierr = 0;
  if (!IsInitialized())
    TIFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (NumSweeps_ == 0) ierr = 1; // Warning: no sweeps performed.
  if (NumSweeps_ < 0)
    TIFPACK_CHK_ERR(-2); // at least one application
  
  Diagonal_ = Teuchos::rcp( new Tpetra_Vector(Matrix().RowMatrixRowMap()) );

  if (Diagonal_ == Teuchos::null)
    TIFPACK_CHK_ERR(-5);

  TIFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*Diagonal_));

  // check diagonal elements, store the inverses, and verify that
  // no zeros are around. If an element is zero, then by default
  // its inverse is zero as well (that is, the row is ignored).
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double& diag = (*Diagonal_)[i];
    if (std::abs(diag) < MinDiagonalValue_)
      diag = MinDiagonalValue_;
    if (diag != 0.0)
      diag = 1.0 / diag;
  }
  ComputeFlops_ += NumMyRows_;

#if 0
  // some methods require the inverse of the diagonal, compute it
  // now and store it in `Diagonal_'
  if ((PrecType_ == TIFPACK_JACOBI) || (PrecType_ == TIFPACK_GS)) {
    Diagonal_->Reciprocal(*Diagonal_);
    // update flops
    ComputeFlops_ += NumMyRows_;
  }
#endif

  // We need to import data from external processors. Here I create an
  // Tpetra_Import object because I cannot assume that Matrix_ has one.
  // This is a bit of waste of resources (but the code is more robust).
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == TIFPACK_GS) || (PrecType_ == TIFPACK_SGS))) {
    Importer_ = Teuchos::rcp( new Tpetra_Import(Matrix().RowMatrixColMap(),
                                  Matrix().RowMatrixRowMap()) );
    if (Importer_ == Teuchos::null) TIFPACK_CHK_ERR(-5);
  }

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  TIFPACK_CHK_ERR(ierr);
  return(0);
}

//==============================================================================
ostream& Tifpack_PointRelaxation::Print(ostream & os) const
{

  double MyMinVal, MyMaxVal;
  double MinVal, MaxVal;

  if (IsComputed_) {
    Diagonal_->MinValue(&MyMinVal);
    Diagonal_->MaxValue(&MyMaxVal);
    Comm().MinAll(&MyMinVal,&MinVal,1);
    Comm().MinAll(&MyMaxVal,&MaxVal,1);
  }

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Tifpack_PointRelaxation" << endl;
    os << "Sweeps         = " << NumSweeps_ << endl;
    os << "damping factor = " << DampingFactor_ << endl;
    if (PrecType_ == TIFPACK_JACOBI)
      os << "Type           = Jacobi" << endl;
    else if (PrecType_ == TIFPACK_GS)
      os << "Type           = Gauss-Seidel" << endl;
    else if (PrecType_ == TIFPACK_SGS)
      os << "Type           = symmetric Gauss-Seidel" << endl;
    if (DoBackwardGS_) 
      os << "Using backward mode (GS only)" << endl;
    if (ZeroStartingSolution_) 
      os << "Using zero starting solution" << endl;
    else
      os << "Using input starting solution" << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Minimum value on stored diagonal = " << MinVal << endl;
      os << "Maximum value on stored diagonal = " << MaxVal << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_ 
       << "  " << std::setw(15) << InitializeTime_ 
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_ 
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_ 
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}

//==============================================================================
double Tifpack_PointRelaxation::
Condest(const Tifpack_CondestType CT, 
        const int MaxIters, const double Tol,
	Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always computes it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Tifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
void Tifpack_PointRelaxation::SetLabel()
{
  string PT;
  if (PrecType_ == TIFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == TIFPACK_GS){
    PT = "GS";
    if(DoBackwardGS_)
      PT = "Backward " + PT;
  }    
  else if (PrecType_ == TIFPACK_SGS)
    PT = "SGS";

  Label_ = "TIFPACK (" + PT + ", sweeps=" + Tifpack_toString(NumSweeps_)
    + ", damping=" + Tifpack_toString(DampingFactor_) + ")";
}

//==============================================================================
// Note that Tifpack_PointRelaxation and Jacobi is much faster than
// Tifpack_AdditiveSchwarz<Tifpack_PointRelaxation> (because of the
// way the matrix-vector produce is performed).
//
// Another ML-related observation is that the starting solution (in Y)
// is NOT supposed to be zero. This may slow down the application of just
// one sweep of Jacobi.
//
int Tifpack_PointRelaxation::
ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  if (!IsComputed())
    TIFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    TIFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr< const Tpetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );
    
  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  // Flops are updated in each of the following. 
  switch (PrecType_) {
  case TIFPACK_JACOBI:
    TIFPACK_CHK_ERR(ApplyInverseJacobi(*Xcopy,Y));
    break;
  case TIFPACK_GS:
    TIFPACK_CHK_ERR(ApplyInverseGS(*Xcopy,Y));
    break;
  case TIFPACK_SGS:
    TIFPACK_CHK_ERR(ApplyInverseSGS(*Xcopy,Y));
    break;
  default:
    TIFPACK_CHK_ERR(-1); // something wrong
  }

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
// This preconditioner can be much slower than AztecOO and ML versions
// if the matrix-vector product requires all ExtractMyRowCopy() 
// (as done through Tifpack_AdditiveSchwarz).
int Tifpack_PointRelaxation::
ApplyInverseJacobi(const Tpetra_MultiVector& RHS, Tpetra_MultiVector& LHS) const
{

  int NumVectors = LHS.NumVectors();
  Tpetra_MultiVector A_times_LHS( LHS.Map(),NumVectors );

  for (int j = 0; j < NumSweeps_ ; j++) {

    TIFPACK_CHK_ERR(Apply(LHS,A_times_LHS));
    TIFPACK_CHK_ERR(A_times_LHS.Update(1.0,RHS,-1.0));
    for (int v = 0 ; v < NumVectors ; ++v)
      TIFPACK_CHK_ERR(LHS(v)->Multiply(DampingFactor_, *(A_times_LHS(v)), 
                                     *Diagonal_, 1.0));

  }

  // Flops:
  // - matrix vector              (2 * NumGlobalNonzeros_)
  // - update                     (2 * NumGlobalRows_)
  // - Multiply:
  //   - DampingFactor            (NumGlobalRows_)
  //   - Diagonal                 (NumGlobalRows_)
  //   - A + B                    (NumGlobalRows_)
  //   - 1.0                      (NumGlobalRows_)
  ApplyInverseFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);

  return(0);
}

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseGS(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  const Tpetra_CrsMatrix* CrsMatrix = dynamic_cast<const Tpetra_CrsMatrix*>(&*Matrix_);
  // try to pick the best option; performances may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    if (CrsMatrix->StorageOptimized())
      return(ApplyInverseGS_FastCrsMatrix(CrsMatrix, X, Y));
    else
      return(ApplyInverseGS_CrsMatrix(CrsMatrix, X, Y));
  }
  else
    return(ApplyInverseGS_RowMatrix(X, Y));
} //ApplyInverseGS()

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseGS_RowMatrix(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  // extract views (for nicer and faster code)
  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int j = 0; j < NumSweeps_ ; j++) {

    // data exchange is here, once per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    // FIXME: do I really need this code below?
    if (NumVectors == 1) {

      double* y0_ptr = y_ptr[0];
      double* y20_ptr = y2_ptr[0];
      double* x0_ptr = x_ptr[0];

      if(!DoBackwardGS_){      
        /* Forward Mode */
        for (int i = 0 ; i < NumMyRows_ ; ++i) {

          int NumEntries;
          int col;
          TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));
          
          double dtemp = 0.0;
          for (int k = 0 ; k < NumEntries ; ++k) {
            
            col = Indices[k];
            dtemp += Values[k] * y20_ptr[col];
          }
          
          y20_ptr[i] += DampingFactor_ * d_ptr[i] * (x0_ptr[i] - dtemp);
        }
      }
      else {
        /* Backward Mode */
        for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

          int NumEntries;
          int col;
          TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));
          double dtemp = 0.0;
          for (int k = 0 ; k < NumEntries ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y20_ptr[i];
          }
          
          y20_ptr[i] += DampingFactor_ * d_ptr[i] * (x0_ptr[i] - dtemp);
        }
      }
      
      // using Export() sounded quite expensive
      if (IsParallel_)
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y0_ptr[i] = y20_ptr[i];

    }
    else {
      if(!DoBackwardGS_){      
        /* Forward Mode */
        for (int i = 0 ; i < NumMyRows_ ; ++i) {
          
          int NumEntries;
          int col;
          TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));
          
          for (int m = 0 ; m < NumVectors ; ++m) {
            
            double dtemp = 0.0;
            for (int k = 0 ; k < NumEntries ; ++k) {
              
              col = Indices[k];
              dtemp += Values[k] * y2_ptr[m][col];
            }
            
            y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
          }
        }
      }
      else {
        /* Backward Mode */
        for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {
          int NumEntries;
          int col;
          TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));

          for (int m = 0 ; m < NumVectors ; ++m) {
            
            double dtemp = 0.0;
            for (int k = 0 ; k < NumEntries ; ++k) {

              col = Indices[k];
              dtemp += Values[k] * y2_ptr[m][col];
            }

            y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);

          }
        }
      }

      // using Export() sounded quite expensive   
      if (IsParallel_)
        for (int m = 0 ; m < NumVectors ; ++m) 
          for (int i = 0 ; i < NumMyRows_ ; ++i)
            y_ptr[m][i] = y2_ptr[m][i];
      
    }
  }

  ApplyInverseFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);

  return(0);
} //ApplyInverseGS_RowMatrix()

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseGS_CrsMatrix(const Tpetra_CrsMatrix* A, const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int* Indices;
  double* Values;

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    if(!DoBackwardGS_){  
      /* Forward Mode */
      for (int i = 0 ; i < NumMyRows_ ; ++i) {

        int NumEntries;
        int col;
        double diag = d_ptr[i];
        
        TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));
        
        for (int m = 0 ; m < NumVectors ; ++m) {
          
          double dtemp = 0.0;
          
          for (int k = 0; k < NumEntries; ++k) {
            
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

        int NumEntries;
        int col;
        double diag = d_ptr[i];
        
        TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));
        
        for (int m = 0 ; m < NumVectors ; ++m) {
          
          double dtemp = 0.0;
          for (int k = 0; k < NumEntries; ++k) {
            
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
          
        }
      }
    }
    
    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
} //ApplyInverseGS_CrsMatrix()

//
//==============================================================================
// ApplyInverseGS_FastCrsMatrix() requires Tpetra_CrsMatrix + OptimizeStorage()

int Tifpack_PointRelaxation::
ApplyInverseGS_FastCrsMatrix(const Tpetra_CrsMatrix* A, const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  TIFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    if(!DoBackwardGS_){  
      /* Forward Mode */
      for (int i = 0 ; i < NumMyRows_ ; ++i) {
        
        int col;
        double diag = d_ptr[i];
        
        for (int m = 0 ; m < NumVectors ; ++m) {
          
          double dtemp = 0.0;
          
          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {
            
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }      
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

        int col;
        double diag = d_ptr[i];
        
        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

        }
      }
    }
    

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
} //ApplyInverseGS_FastCrsMatrix()

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseSGS(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  const Tpetra_CrsMatrix* CrsMatrix = dynamic_cast<const Tpetra_CrsMatrix*>(&*Matrix_);
  // try to pick the best option; performances may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    if (CrsMatrix->StorageOptimized())
      return(ApplyInverseSGS_FastCrsMatrix(CrsMatrix, X, Y));
    else
      return(ApplyInverseSGS_CrsMatrix(CrsMatrix, X, Y));
  }
  else
    return(ApplyInverseSGS_RowMatrix(X, Y));
}

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseSGS_RowMatrix(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      TIFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
}

//==============================================================================
int Tifpack_PointRelaxation::
ApplyInverseSGS_CrsMatrix(const Tpetra_CrsMatrix* A, const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int* Indices;
  double* Values;

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = 0; k < NumEntries; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = 0; k < NumEntries; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
}

//==============================================================================
// this requires Tpetra_CrsMatrix + OptimizeStorage()
int Tifpack_PointRelaxation::
ApplyInverseSGS_FastCrsMatrix(const Tpetra_CrsMatrix* A, const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  TIFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Tpetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      TIFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //TIFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
}
