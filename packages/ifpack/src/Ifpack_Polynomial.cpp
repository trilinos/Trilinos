/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include <iomanip>
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Ifpack_Polynomial.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include <complex>
#ifdef HAVE_IFPACK_AZTECOO
#include "Ifpack_DiagPreconditioner.h"
#include "AztecOO.h"
#endif

#ifdef HAVE_IFPACK_EPETRAEXT
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#endif


#define ABS(x) ((x)>0?(x):-(x))

// Helper function for normal equations
inline void Apply_Transpose(Teuchos::RCP<const Epetra_Operator> Operator_,const Epetra_MultiVector &X,Epetra_MultiVector &Y){
  Epetra_Operator * Operator=const_cast<Epetra_Operator*>(&*Operator_);
  Operator->SetUseTranspose(true);
  Operator->Apply(X,Y);
  Operator->SetUseTranspose(false);
}




//==============================================================================
// NOTE: any change to the default values should be committed to the other
//       constructor as well.
Ifpack_Polynomial::
Ifpack_Polynomial(const Epetra_Operator* Operator) :
  IsInitialized_(false),
  IsComputed_(false),
  IsIndefinite_(false),
  IsComplex_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  PolyDegree_(3),
  LSPointsReal_(10),
  LSPointsImag_(10),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  RealEigRatio_(10.0),
  ImagEigRatio_(10.0),
  Label_(),
  LambdaRealMin_(0.0),
  LambdaRealMax_(-1.0),
  LambdaImagMin_(0.0),
  LambdaImagMax_(0.0),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  UseBlockMode_(false),
  SolveNormalEquations_(false),
  IsRowMatrix_(false),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
// NOTE: This constructor has been introduced because SWIG does not appear
//       to appreciate dynamic_cast. An instruction of type
//       Matrix_ = dynamic_cast<const Epetra_RowMatrix*> in the
//       other construction does not work in PyTrilinos -- of course
//       it does in any C++ code (for an Epetra_Operator that is also
//       an Epetra_RowMatrix).
//
// FIXME: move declarations into a separate method?
Ifpack_Polynomial::
Ifpack_Polynomial(const Epetra_RowMatrix* Operator) :
  IsInitialized_(false),
  IsComputed_(false),
  IsIndefinite_(false),
  IsComplex_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  PolyDegree_(3),
  LSPointsReal_(10),
  LSPointsImag_(10),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  RealEigRatio_(10.0),
  ImagEigRatio_(10.0),
  EigMaxIters_(10),
  Label_(),
  LambdaRealMin_(0.0),
  LambdaRealMax_(-1.0),
  LambdaImagMin_(0.0),
  LambdaImagMax_(0.0),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  Matrix_(Teuchos::rcp(Operator,false)),
  UseBlockMode_(false),
  SolveNormalEquations_(false),
  IsRowMatrix_(true),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
int Ifpack_Polynomial::SetParameters(Teuchos::ParameterList& List)
{

  RealEigRatio_         = List.get("polynomial: real eigenvalue ratio", RealEigRatio_);
  ImagEigRatio_         = List.get("polynomial: imag eigenvalue ratio", ImagEigRatio_);
  LambdaRealMin_        = List.get("polynomial: min real part", LambdaRealMin_);
  LambdaRealMax_        = List.get("polynomial: max real part", LambdaRealMax_);
  LambdaImagMin_        = List.get("polynomial: min imag part", LambdaImagMin_);
  LambdaImagMax_        = List.get("polynomial: max imag part", LambdaImagMax_);
  PolyDegree_           = List.get("polynomial: degree",PolyDegree_);
  LSPointsReal_         = List.get("polynomial: real interp points",LSPointsReal_);
  LSPointsImag_         = List.get("polynomial: imag interp points",LSPointsImag_);
  IsIndefinite_         = List.get("polynomial: indefinite",IsIndefinite_);
  IsComplex_            = List.get("polynomial: complex",IsComplex_);
  MinDiagonalValue_     = List.get("polynomial: min diagonal value",
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("polynomial: zero starting solution",
                                   ZeroStartingSolution_);

  Epetra_Vector* ID     = List.get("polynomial: operator inv diagonal",
                                   (Epetra_Vector*)0);
  EigMaxIters_          = List.get("polynomial: eigenvalue max iterations",EigMaxIters_);

#ifdef HAVE_IFPACK_EPETRAEXT
  // This is *always* false if EpetraExt isn't enabled
  UseBlockMode_         = List.get("polynomial: use block mode",UseBlockMode_);
  if(!List.isParameter("polynomial: block list")) UseBlockMode_=false;
  else{
    BlockList_          = List.get("polynomial: block list",BlockList_);

    // Since we know we're doing a matrix inverse, clobber the block list
    // w/"invert" if it's set to multiply
    Teuchos::ParameterList Blist;
    Blist=BlockList_.get("blockdiagmatrix: list",Blist);
    string dummy("invert");
    string ApplyMode=Blist.get("apply mode",dummy);
    if(ApplyMode==string("multiply")){
      Blist.set("apply mode","invert");
      BlockList_.set("blockdiagmatrix: list",Blist);
    }
  }
#endif

  SolveNormalEquations_ = List.get("polynomial: solve normal equations",SolveNormalEquations_);

  if (ID != 0)
  {
    InvDiagonal_ = Teuchos::rcp( new Epetra_Vector(*ID) );
  }

  SetLabel();

  return(0);
}

//==============================================================================
const Epetra_Comm& Ifpack_Polynomial::Comm() const
{
  return(Operator_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_Polynomial::OperatorDomainMap() const
{
  return(Operator_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_Polynomial::OperatorRangeMap() const
{
  return(Operator_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_Polynomial::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  if (IsRowMatrix_)
  {
    IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  }
  else
  {
    IFPACK_CHK_ERR(Operator_->Apply(X,Y));
  }

  return(0);
}

//==============================================================================
int Ifpack_Polynomial::Initialize()
{
  IsInitialized_ = false;

  if (Operator_ == Teuchos::null)
    IFPACK_CHK_ERR(-2);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );

  if (IsRowMatrix_)
  {
    if (Matrix().NumGlobalRows64() != Matrix().NumGlobalCols64())
      IFPACK_CHK_ERR(-2); // only square matrices

    NumMyRows_ = Matrix_->NumMyRows();
    NumMyNonzeros_ = Matrix_->NumMyNonzeros();
    NumGlobalRows_ = Matrix_->NumGlobalRows64();
    NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros64();
  }
  else
  {
    if (Operator_->OperatorDomainMap().NumGlobalElements64() !=
        Operator_->OperatorRangeMap().NumGlobalElements64())
      IFPACK_CHK_ERR(-2); // only square operators
  }

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_Polynomial::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (PolyDegree_ <= 0)
    IFPACK_CHK_ERR(-2); // at least one application

#ifdef HAVE_IFPACK_EPETRAEXT
  // Check to see if we can run in block mode
  if(IsRowMatrix_ && InvDiagonal_ == Teuchos::null && UseBlockMode_){
    const Epetra_CrsMatrix *CrsMatrix=dynamic_cast<const Epetra_CrsMatrix*>(&*Matrix_);

    // If we don't have CrsMatrix, we can't use the block preconditioner
    if(!CrsMatrix) UseBlockMode_=false;
    else{
      int ierr;
      InvBlockDiagonal_=Teuchos::rcp(new EpetraExt_PointToBlockDiagPermute(*CrsMatrix));
      if(InvBlockDiagonal_==Teuchos::null) IFPACK_CHK_ERR(-6);

      ierr=InvBlockDiagonal_->SetParameters(BlockList_);
      if(ierr) IFPACK_CHK_ERR(-7);

      ierr=InvBlockDiagonal_->Compute();
      if(ierr) IFPACK_CHK_ERR(-8);
    }
  }
#endif
  if (IsRowMatrix_ && InvDiagonal_ == Teuchos::null && !UseBlockMode_)
  {
    InvDiagonal_ = Teuchos::rcp( new Epetra_Vector(Matrix().Map()) );

    if (InvDiagonal_ == Teuchos::null)
      IFPACK_CHK_ERR(-5);

    IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*InvDiagonal_));

    // Inverse diagonal elements
    // Replace zeros with 1.0
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      double diag = (*InvDiagonal_)[i];
      if (IFPACK_ABS(diag) < MinDiagonalValue_)
        (*InvDiagonal_)[i] = MinDiagonalValue_;
      else
        (*InvDiagonal_)[i] = 1.0 / diag;
    }
  }

  // Automatically compute maximum eigenvalue estimate of D^{-1}A if user hasn't provided one
  double lambda_real_min, lambda_real_max, lambda_imag_min, lambda_imag_max;
  if (LambdaRealMax_ == -1) {
    //PowerMethod(Matrix(), *InvDiagonal_, EigMaxIters_, lambda_max);
    GMRES(Matrix(), *InvDiagonal_, EigMaxIters_, lambda_real_min, lambda_real_max, lambda_imag_min, lambda_imag_max);
    LambdaRealMin_=lambda_real_min;  LambdaImagMin_=lambda_imag_min;
    LambdaRealMax_=lambda_real_max;  LambdaImagMax_=lambda_imag_max;
    //std::cout<<"LambdaRealMin: "<<LambdaRealMin_<<std::endl;
    //std::cout<<"LambdaRealMax: "<<LambdaRealMax_<<std::endl;
    //std::cout<<"LambdaImagMin: "<<LambdaImagMin_<<std::endl;
    //std::cout<<"LambdaImagMax: "<<LambdaImagMax_<<std::endl;
  }

  // find least squares polynomial for (LSPointsReal_*LSPointsImag_) zeros
  // on a rectangle in the complex plane defined as
  // [LambdaRealMin_,LambdaRealMax_] x [LambdaImagMin_,LambdaImagMax_]

  std::complex<double> zero(0.0,0.0);
  std::complex<double> one(1.0,0.0);

  // Compute points in complex plane
  double lenx = LambdaRealMax_-LambdaRealMin_;
  int      nx = ceil(lenx*((double) LSPointsReal_));
  if (nx<2) { nx = 2; }
  double   hx = lenx/((double) nx);
  std::vector<double> xs;
  if(abs(lenx)>1.0e-8) {
    for( int pt=0; pt<=nx; pt++ ) {
      xs.push_back(hx*pt+LambdaRealMin_);
    }
  }
  else {
    xs.push_back(LambdaRealMax_);
    nx=1;
  }
  double leny = LambdaImagMax_-LambdaImagMin_;
  int      ny = ceil(leny*((double) LSPointsImag_));
  if (ny<2) { ny = 2; }
  double   hy = leny/((double) ny);
  std::vector<double> ys;
  if(abs(leny)>1.0e-8) {
    for( int pt=0; pt<=ny; pt++ ) {
      ys.push_back(hy*pt+LambdaImagMin_);
    }
  }
  else {
    ys.push_back(LambdaImagMax_);
    ny=1;
  }
  std::vector< std::complex<double> > cpts;
  for( int jj=0; jj<ny; jj++ ) {
    for( int ii=0; ii<nx; ii++ ) {
      std::complex<double> cpt(xs[ii],ys[jj]);
      cpts.push_back(cpt);
    }
  }
  cpts.push_back(zero);

#ifdef HAVE_TEUCHOS_COMPLEX

  // Construct overdetermined Vandermonde matrix
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Vmatrix(cpts.size(),PolyDegree_+1);
  Vmatrix.putScalar(zero);
  for (int jj = 0; jj <= PolyDegree_; ++jj) {
    for (int ii = 0; ii < static_cast<int> (cpts.size ()) - 1; ++ii) {
      if (jj > 0) {
        Vmatrix(ii,jj) = pow(cpts[ii],jj);
      }
      else {
        Vmatrix(ii,jj) = one;
      }
    }
  }
  Vmatrix(cpts.size()-1,0)=one;

  // Right hand side: all zero except last entry
  Teuchos::SerialDenseMatrix< int,std::complex<double> > RHS(cpts.size(),1);
  RHS.putScalar(zero);
  RHS(cpts.size()-1,0)=one;

  // Solve least squares problem using LAPACK
  Teuchos::LAPACK< int, std::complex<double> > lapack;
  const int N = Vmatrix.numCols();
  Teuchos::Array<double> singularValues(N);
  Teuchos::Array<double> rwork(1);
  rwork.resize (std::max (1, 5 * N));
  std::complex<double> lworkScalar(1.0,0.0);
  int info = 0;
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(), RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(), RHS.values(),    RHS.numRows(),
                   &lworkScalar, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                             "_GELSS workspace query returned INFO = "
                             << info << " != 0.");
  const int lwork = static_cast<int> (real(lworkScalar));
  TEUCHOS_TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
                             "_GELSS workspace query returned LWORK = "
                             << lwork << " < 0.");
  // Allocate workspace.  Size > 0 means &work[0] makes sense.
  Teuchos::Array< std::complex<double> > work (std::max (1, lwork));
  // Solve the least-squares problem.
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(),  RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(),  RHS.values(),   RHS.numRows(),
                   &work[0], lwork, &info);

  coeff_.resize(PolyDegree_+1);
  std::complex<double> c0=RHS(0,0);
  for(int ii=0; ii<=PolyDegree_; ii++) {
    // test that the imaginary part is nonzero
    //TEUCHOS_TEST_FOR_EXCEPTION(abs(imag(RHS(ii,0))) > 1e-8, std::logic_error,
    //                         "imaginary part of polynomial coefficients is nonzero! coeff = "
    //                         << RHS(ii,0));
    coeff_[ii]=real(RHS(ii,0)/c0);
    //std::cout<<"coeff["<<ii<<"]="<<coeff_[ii]<<std::endl;
  }

#else

  // Construct overdetermined Vandermonde matrix
  Teuchos::SerialDenseMatrix< int, double > Vmatrix(xs.size()+1,PolyDegree_+1);
  Vmatrix.putScalar(0.0);
  for( int jj=0; jj<=PolyDegree_; jj++) {
    for( int ii=0; ii<xs.size(); ii++) {
      if(jj>0) {
        Vmatrix(ii,jj)=pow(xs[ii],jj);
      }
      else {
        Vmatrix(ii,jj)=1.0;
      }
    }
  }
  Vmatrix(xs.size(),0)=1.0;

  // Right hand side: all zero except last entry
  Teuchos::SerialDenseMatrix< int, double > RHS(xs.size()+1,1);
  RHS.putScalar(0.0);
  RHS(xs.size(),0)=1.0;

  // Solve least squares problem using LAPACK
  Teuchos::LAPACK< int, double > lapack;
  const int N = Vmatrix.numCols();
  Teuchos::Array<double> singularValues(N);
  Teuchos::Array<double> rwork(1);
  rwork.resize (std::max (1, 5 * N));
  double lworkScalar(1.0);
  int info = 0;
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(), RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(), RHS.values(),    RHS.numRows(),
                   &lworkScalar, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                             "_GELSS workspace query returned INFO = "
                             << info << " != 0.");
  const int lwork = static_cast<int> (lworkScalar);
  TEUCHOS_TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
                             "_GELSS workspace query returned LWORK = "
                             << lwork << " < 0.");
  // Allocate workspace.  Size > 0 means &work[0] makes sense.
  Teuchos::Array< double > work (std::max (1, lwork));
  // Solve the least-squares problem.
  lapack.GELS('N', Vmatrix.numRows(), Vmatrix.numCols(),  RHS.numCols(),
                   Vmatrix.values(),  Vmatrix.numRows(),  RHS.values(),   RHS.numRows(),
                   &work[0], lwork, &info);

  coeff_.resize(PolyDegree_+1);
  double c0=RHS(0,0);
  for(int ii=0; ii<=PolyDegree_; ii++) {
    // test that the imaginary part is nonzero
    //TEUCHOS_TEST_FOR_EXCEPTION(abs(imag(RHS(ii,0))) > 1e-8, std::logic_error,
    //                         "imaginary part of polynomial coefficients is nonzero! coeff = "
    //                         << RHS(ii,0));
    coeff_[ii]=RHS(ii,0)/c0;
  }

#endif

#ifdef IFPACK_FLOPCOUNTERS
  ComputeFlops_ += NumMyRows_;
#endif

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Ifpack_Polynomial::Print(ostream & os) const
{

  double MyMinVal, MyMaxVal;
  double MinVal, MaxVal;

  if (IsComputed_) {
    InvDiagonal_->MinValue(&MyMinVal);
    InvDiagonal_->MaxValue(&MyMaxVal);
    Comm().MinAll(&MyMinVal,&MinVal,1);
    Comm().MinAll(&MyMaxVal,&MaxVal,1);
  }

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Polynomial" << endl;
    os << "Degree of polynomial      = " << PolyDegree_ << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows     = " << Operator_->OperatorRangeMap().NumGlobalElements64() << endl;
    if (IsComputed_) {
      os << "Minimum value on stored inverse diagonal = " << MinVal << endl;
      os << "Maximum value on stored inverse diagonal = " << MaxVal << endl;
    }
    if (ZeroStartingSolution_)
      os << "Using zero starting solution" << endl;
    else
      os << "Using input starting solution" << endl;
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
double Ifpack_Polynomial::
Condest(const Ifpack_CondestType CT,
        const int MaxIters, const double Tol,
        Epetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always computes it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
void Ifpack_Polynomial::SetLabel()
{
  Label_ = "IFPACK (Least squares polynomial), degree=" + Ifpack_toString(PolyDegree_);
}

//==============================================================================
int Ifpack_Polynomial::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (PolyDegree_ == 0)
    return 0;

  int nVec = X.NumVectors();
  if (nVec != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  Epetra_MultiVector Xcopy(X);
  if(ZeroStartingSolution_==true) {
    Y.PutScalar(0.0);
  }

  // mfh 20 Mar 2014: IBD never gets used, so I'm commenting out the
  // following lines of code in order to forestall build warnings.
// #ifdef HAVE_IFPACK_EPETRAEXT
//   EpetraExt_PointToBlockDiagPermute* IBD=0;
//   if (UseBlockMode_) IBD=&*InvBlockDiagonal_;
// #endif

  Y.Update(-coeff_[1], Xcopy, 1.0);
  for (int ii = 2; ii < static_cast<int> (coeff_.size ()); ++ii) {
    const Epetra_MultiVector V(Xcopy);
    Operator_->Apply(V,Xcopy);
    Xcopy.Multiply(1.0, *InvDiagonal_, Xcopy, 0.0);
    // Update Y
    Y.Update(-coeff_[ii], Xcopy, 1.0);
  }

  // Flops are updated in each of the following.
  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_Polynomial::
PowerMethod(const Epetra_Operator& Operator,
            const Epetra_Vector& InvPointDiagonal,
            const int MaximumIterations,
            double& lambda_max)
{
  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Epetra_Vector x(Operator.OperatorDomainMap());
  Epetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) IFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator.Apply(x, y);
    IFPACK_CHK_ERR(y.Multiply(1.0, InvPointDiagonal, y, 0.0));
    IFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    IFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    IFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) IFPACK_CHK_ERR(-1);
    IFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}

//==============================================================================
int Ifpack_Polynomial::
CG(const Epetra_Operator& Operator,
   const Epetra_Vector& InvPointDiagonal,
   const int MaximumIterations,
   double& lambda_min, double& lambda_max)
{
#ifdef HAVE_IFPACK_AZTECOO
  Epetra_Vector x(Operator.OperatorDomainMap());
  Epetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);

  Epetra_LinearProblem LP(const_cast<Epetra_Operator*>(&Operator), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  Ifpack_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use the CG estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif
}

//==============================================================================
#ifdef HAVE_IFPACK_EPETRAEXT
int Ifpack_Polynomial::
PowerMethod(const int MaximumIterations,  double& lambda_max)
{

  if(!UseBlockMode_) IFPACK_CHK_ERR(-1);
  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Epetra_Vector x(Operator_->OperatorDomainMap());
  Epetra_Vector y(Operator_->OperatorRangeMap());
  Epetra_Vector z(Operator_->OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) IFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator_->Apply(x, z);
    InvBlockDiagonal_->ApplyInverse(z,y);
    if(SolveNormalEquations_){
      InvBlockDiagonal_->ApplyInverse(y,z);
      Apply_Transpose(Operator_,z, y);
    }

    IFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    IFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    IFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) IFPACK_CHK_ERR(-1);
    IFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}
#endif

//==============================================================================
#ifdef HAVE_IFPACK_EPETRAEXT
int Ifpack_Polynomial::
CG(const int MaximumIterations,
   double& lambda_min, double& lambda_max)
{
  IFPACK_CHK_ERR(-1);// NTS: This always seems to yield errors in AztecOO, ergo,
                     // I turned it off.

  if(!UseBlockMode_) IFPACK_CHK_ERR(-1);

#ifdef HAVE_IFPACK_AZTECOO
  Epetra_Vector x(Operator_->OperatorDomainMap());
  Epetra_Vector y(Operator_->OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);
  Epetra_LinearProblem LP(const_cast<Epetra_RowMatrix*>(&*Matrix_), &x, &y);

  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  solver.SetPrecOperator(&*InvBlockDiagonal_);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use the CG estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif
}
#endif

//==============================================================================
int Ifpack_Polynomial::
GMRES(const Epetra_Operator& Operator,
      const Epetra_Vector& InvPointDiagonal,
      const int MaximumIterations,
      double& lambda_real_min, double& lambda_real_max,
      double& lambda_imag_min, double& lambda_imag_max)
{
#ifdef HAVE_IFPACK_AZTECOO
  Epetra_Vector x(Operator_->OperatorDomainMap());
  Epetra_Vector y(Operator_->OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);
  Epetra_LinearProblem LP(const_cast<Epetra_RowMatrix*>(&*Matrix_), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);
  Ifpack_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);
  const double* status = solver.GetAztecStatus();
  lambda_real_min = status[AZ_lambda_real_min];
  lambda_real_max = status[AZ_lambda_real_max];
  lambda_imag_min = status[AZ_lambda_imag_min];
  lambda_imag_max = status[AZ_lambda_imag_max];
  return(0);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use the GMRES estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif
}
