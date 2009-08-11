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
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_IKLU.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Utils.hpp"
#include "Tifpack_HPPashTable.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Util.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <functional>

using namespace Teuchos;

//==============================================================================
// FIXME: allocate Comm_ and Time_ the first Initialize() call
Tifpack_IKLU::Tifpack_IKLU(const Tpetra_RowMatrix* A) :
  A_(*A),
  Comm_(A->Comm()),
  Condest_(-1.0),
  Relax_(0.),
  Athresh_(0.0),
  Rthresh_(1.0),
  LevelOfFill_(1.0),
  DropTolerance_(1e-12),
  IsInitialized_(false),
  IsComputed_(false),
  UseTranspose_(false),
  NumMyRows_(-1),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(Comm()),
  GlobalNonzeros_(0),
  csrA_(0),
  cssS_(0),
  csrnN_(0)
{
  // do nothing here..
}

//==============================================================================
Tifpack_IKLU::~Tifpack_IKLU()
{
  Destroy();
}

//==============================================================================
void Tifpack_IKLU::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
  if (csrA_)
    csr_spfree( csrA_ );
  if (cssS_)
    csr_sfree( cssS_ );
  if (csrnN_)
    csr_nfree( csrnN_ );
}

//==========================================================================
int Tifpack_IKLU::SetParameters(Teuchos::ParameterList& List)
{
  try 
  {
    LevelOfFill_ = List.get<double>("fact: ilut level-of-fill", LevelOfFill());
    if (LevelOfFill_ <= 0.0)
      TIFPACK_CHK_ERR(-2); // must be greater than 0.0

    Athresh_ = List.get<double>("fact: absolute threshold", Athresh_);
    Rthresh_ = List.get<double>("fact: relative threshold", Rthresh_);
    Relax_ = List.get<double>("fact: relax value", Relax_);
    DropTolerance_ = List.get<double>("fact: drop tolerance", DropTolerance_);

    Label_ = "TIFPACK IKLU (fill=" + Tifpack_toString(LevelOfFill())
      + ", relax=" + Tifpack_toString(RelaxValue())
      + ", athr=" + Tifpack_toString(AbsoluteThreshold())
      + ", rthr=" + Tifpack_toString(RelativeThreshold())
      + ", droptol=" + Tifpack_toString(DropTolerance())
      + ")";
    return(0);
  }
  catch (...)
  {
    cerr << "Caught an exception while parsing the parameter list" << endl;
    cerr << "This typically means that a parameter was set with the" << endl;
    cerr << "wrong type (for example, int instead of double). " << endl;
    cerr << "please check the documentation for the type required by each parameer." << endl;
    TIFPACK_CHK_ERR(-1);
  }

  return(0);
}

//==========================================================================
int Tifpack_IKLU::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.ResetStartTime();

  if (A_.Comm().NumProc() != 1) {
    cout << " There are too many processors !!! " << endl;
    cerr << "Tifpack_IKLU can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a subdomain solver for Tifpack_AdditiveSchwarz," << endl;
    cerr << "and it is currently not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }
  
  // check dimensions of input matrix only in serial
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    TIFPACK_CHK_ERR(-2);
    
  NumMyRows_ = Matrix().NumMyRows();
  NumMyNonzeros_ = Matrix().NumMyNonzeros();

  int RowNnz, Length = Matrix().MaxNumEntries();
  vector<int>    RowIndices(Length);
  vector<double> RowValues(Length);

  //cout << "Processor " << Comm().MyPID() << " owns " << NumMyRows_ << " rows and has " << NumMyNonzeros_ << " nonzeros " << endl;
  // get general symbolic structure of the matrix
  csrA_ = csr_spalloc( NumMyRows_, NumMyRows_, NumMyNonzeros_, 1, 0 );

  // copy the symbolic structure into csrA_
  int count = 0;
  csrA_->p[0] = 0;
  for (int i = 0; i < NumMyRows_; ++i ) {

    TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(i,Length,RowNnz,
				       &RowValues[0],&RowIndices[0]));
    for (int j = 0 ; j < RowNnz ; ++j) {
      csrA_->j[count++] = RowIndices[j];
      //cout << "Row = " << i << ", Column = " << RowIndices[j] << ", Value = " << RowValues[j] << endl;
    }
    csrA_->p[i+1] = csrA_->p[i] + RowNnz;
  }

  // Perform symbolic analysis on the current matrix structure
  int order = 1;
  cssS_ = csr_sqr( order, csrA_ );
  for (int i = 0; i < NumMyRows_; ++i ) {
     cout << "AMD Perm (from inside KLU) [" << i << "] = " << cssS_->q[i] << endl;
  }

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
class Tifpack_AbsComp 
{
 public:
  inline bool operator()(const double& x, const double& y) 
  {
    return(TIFPACK_ABS(x) > TIFPACK_ABS(y));
  }
};

//==========================================================================

int Tifpack_IKLU::Compute() 
{
  if (!IsInitialized()) 
    TIFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();

  bool distributed = (Comm().NumProc() > 1)?true:false;
  if (distributed)
  {
    SerialComm_ = rcp(new Tpetra_SerialComm);
    SerialMap_ = rcp(new Tpetra_Map(NumMyRows_, 0, *SerialComm_));
    assert (SerialComm_.get() != 0);
    assert (SerialMap_.get() != 0);
  }
  else
    SerialMap_ = rcp(const_cast<Tpetra_Map*>(&A_.RowMatrixRowMap()), false);

  int RowNnz;
  vector<int>    RowIndices(Length);
  vector<double> RowValues(Length);

  // copy the values from A_ into csrA_
  int count = 0;
  for (int i = 0; i < NumMyRows_; ++i ) {

    TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(i,Length,RowNnz,
				       &RowValues[0],&RowIndices[0]));
    // make sure each row has the same number of nonzeros
    if (RowNnz != (csrA_->p[i+1]-csrA_->p[i])) {
      cout << "The number of nonzeros for this row does not math the expected number of nonzeros!!!" << endl;
    }
    for (int j = 0 ; j < RowNnz ; ++j) {
      
      csrA_->x[count++] = RowValues[j];
      //cout << "Row = " << i << ", Column = " << RowIndices[j] << ", Value = " << RowValues[j] << endl;
    }
  }
  
  // compute the lu factors
  double tol = 0.1;
  csrnN_ = csr_lu( &*csrA_, &*cssS_, tol );

  // Create L and U as a view of the information stored in csrnN_->L and csrnN_->U
  csr* L_tmp = csrnN_->L;
  csr* U_tmp = csrnN_->U;
  vector<int> numEntriesL( NumMyRows_ ), numEntriesU( NumMyRows_ );
  for (int i=0; i < NumMyRows_; ++i) {
    numEntriesL[i] = ( L_tmp->p[i+1] - L_tmp->p[i] );
    numEntriesU[i] = ( U_tmp->p[i+1] - U_tmp->p[i] );
  }
  L_ = rcp(new Tpetra_CrsMatrix(View, *SerialMap_, &numEntriesL[0]));
  U_ = rcp(new Tpetra_CrsMatrix(View, *SerialMap_, &numEntriesU[0]));

  // Insert the values into L and U
  for (int i=0; i < NumMyRows_; ++i) {
    L_->InsertGlobalValues( i, numEntriesL[i], &(L_tmp->x[L_tmp->p[i]]), &(L_tmp->j[L_tmp->p[i]]) );
    U_->InsertGlobalValues( i, numEntriesU[i], &(U_tmp->x[U_tmp->p[i]]), &(U_tmp->j[U_tmp->p[i]]) );
  }

  TIFPACK_CHK_ERR(L_->FillComplete());
  TIFPACK_CHK_ERR(U_->FillComplete());

  int MyNonzeros = L_->NumGlobalNonzeros() + U_->NumGlobalNonzeros();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;

  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}
  
//=============================================================================
int Tifpack_IKLU::ApplyInverse(const Tpetra_MultiVector& X, 
			     Tpetra_MultiVector& Y) const
{
  if (!IsComputed())
    TIFPACK_CHK_ERR(-2); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    TIFPACK_CHK_ERR(-3); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  // NOTE: L_ and U_ are based on SerialMap_, while Xcopy is based
  // on A.Map()... which are in general different. However, Solve()
  // does not seem to care... which is fine with me.
  //
  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy and apply permutation.
  vector<int> invq( NumMyRows_ );

  for (int i=0; i<NumMyRows_; ++i ) {
    csrnN_->perm[ csrnN_->pinv[i] ] = i;
    invq[ cssS_->q[i] ] = i;
  }

  Teuchos::RefCountPtr<Tpetra_MultiVector> Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X.Map(),X.NumVectors()), false );
  Teuchos::RefCountPtr<Tpetra_MultiVector> Ytemp = Teuchos::rcp( new Tpetra_MultiVector(Y.Map(),Y.NumVectors()) );

  for (int i=0; i<NumMyRows_; ++i) {
    for (int j=0; j<X.NumVectors(); ++j) {
      Xcopy->ReplaceMyValue( invq[i], j, (*X(j))[i] );
    }
  }

  if (!UseTranspose_)
  {
    // solves LU Y = X 
    TIFPACK_CHK_ERR(L_->Solve(false,false,false,*Xcopy,*Ytemp));
    TIFPACK_CHK_ERR(U_->Solve(true,false,false,*Ytemp,*Ytemp));
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    TIFPACK_CHK_ERR(U_->Solve(true,true,false,*Xcopy,*Ytemp));
    TIFPACK_CHK_ERR(L_->Solve(false,true,false,*Ytemp,*Ytemp));
  }

  // Reverse the permutation.
  for (int i=0; i<NumMyRows_; ++i) {
    for (int j=0; j<Y.NumVectors(); ++j) {
      Y.ReplaceMyValue( csrnN_->perm[i], j, (*(*Ytemp)(j))[i] );
    }
  }

  ++NumApplyInverse_;
  ApplyInverseFlops_ += X.NumVectors() * 2 * GlobalNonzeros_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Tifpack_IKLU::Apply(const Tpetra_MultiVector& X, 
		      Tpetra_MultiVector& Y) const 
{

  return(-98);
}

//=============================================================================
double Tifpack_IKLU::Condest(const Tifpack_CondestType CT, 
                            const int MaxIters, const double Tol,
			    Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Tifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Tifpack_IKLU::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Tifpack_IKLU: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Relax value        = " << RelaxValue() << endl;
    os << "Condition number estimate       = " << Condest() << endl;
    os << "Global number of rows           = " << A_.NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros in A         = " << A_.NumGlobalNonzeros() << endl;
      os << "Number of nonzeros in L + U     = " << NumGlobalNonzeros() 
         << " ( = " << 100.0 * NumGlobalNonzeros() / A_.NumGlobalNonzeros() 
         << " % of A)" << endl;
      os << "nonzeros / rows                 = " 
        << 1.0 * NumGlobalNonzeros() / U_->NumGlobalRows() << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "               0.0            0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops(); 
    if (ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}
