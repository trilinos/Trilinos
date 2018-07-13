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
#ifdef HAVE_IFPACK_SPARSKIT
#include "Ifpack_Preconditioner.h"
#include "Ifpack_SPARSKIT.h"
#include "Ifpack_Condest.h"
#include "Ifpack_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"
#include "Teuchos_ParameterList.hpp"

#define F77_ILUT  F77_FUNC(ilut, ILUT)
#define F77_ILUD  F77_FUNC(ilud, ILUD)
#define F77_ILUTP F77_FUNC(ilutp, ILUTP)
#define F77_ILUDP F77_FUNC(iludp, ILUDP)
#define F77_ILUK  F77_FUNC(iluk, ILUK)
#define F77_ILU0  F77_FUNC(ilu0, ILU0)
#define F77_LUSOL F77_FUNC(lusol, LUSOL)

extern "C" {
  void F77_ILUT(int *, double*, int*, int*, int*, double*,
                double*, int*, int*, int*, double*, int*, int*);
  void F77_ILUD(int *, double*, int*, int*, double*, double*,
                double*, int*, int*, int*, double*, int*, int*);
  void F77_ILUTP(int *, double*, int*, int*, int*, double*, double*, int*,
                 double*, int*, int*, int*, double*, int*, int*, int*);
  void F77_ILUDP(int *, double*, int*, int*, double*, double*, double*, int*,
                 double*, int*, int*, int*, double*, int*, int*, int*);
  void F77_ILUK(int *, double*, int*, int*, int*, 
                double*, int*, int*, int*, int*, double*, int*, int*);
  void F77_ILU0(int*, double*, int*, int*, double*, int*, int*, int*, int*);
  void F77_LUSOL(int *, double*, double*, double*, int*, int*);
}


//==============================================================================
Ifpack_SPARSKIT::Ifpack_SPARSKIT(Epetra_RowMatrix* A) :
  A_(*A),
  Comm_(A->Comm()),
  UseTranspose_(false),
  lfil_(0),
  droptol_(0.0),
  tol_(0.0),
  permtol_(0.1),
  alph_(0.0),
  mbloc_(-1),
  Type_("ILUT"),
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0),
  ApplyInverseTime_(0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0)
{
  Teuchos::ParameterList List;
  SetParameters(List);
}

//==============================================================================
Ifpack_SPARSKIT::~Ifpack_SPARSKIT()
{
}

//==========================================================================
int Ifpack_SPARSKIT::SetParameters(Teuchos::ParameterList& List)
{
  lfil_    = List.get("fact: sparskit: lfil", lfil_);
  tol_     = List.get("fact: sparskit: tol", tol_);
  droptol_ = List.get("fact: sparskit: droptol", droptol_);
  permtol_ = List.get("fact: sparskit: permtol", permtol_);
  alph_    = List.get("fact: sparskit: alph", alph_);
  mbloc_   = List.get("fact: sparskit: mbloc", mbloc_);
  Type_    = List.get("fact: sparskit: type", Type_);

  // set label
  Label_ = "IFPACK SPARSKIT (Type=" + Type_ + ", fill=" +
    Ifpack_toString(lfil_) + ")";

  return(0);
}

//==========================================================================
int Ifpack_SPARSKIT::Initialize()
{
  IsInitialized_ = true;
  IsComputed_    = false;
  return(0);
}

//==========================================================================
int Ifpack_SPARSKIT::Compute() 
{
  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;

  // convert the matrix into SPARSKIT format. The matrix is then
  // free'd after method Compute() returns.

  // convert the matrix into CSR format. Note that nnz is an over-estimate,
  // since it contains the nonzeros corresponding to external nodes as well.
  int n   = Matrix().NumMyRows();
  int nnz = Matrix().NumMyNonzeros();

  std::vector<double> a(nnz);
  std::vector<int>    ja(nnz);
  std::vector<int>    ia(n + 1);

  const int MaxNumEntries = Matrix().MaxNumEntries();

  std::vector<double> Values(MaxNumEntries);
  std::vector<int>    Indices(MaxNumEntries);

  int count = 0;

  ia[0] = 1;
  for (int i = 0 ; i < n ; ++i)
  {
    int NumEntries;
    int NumMyEntries = 0;
    Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumEntries, &Values[0], 
                              &Indices[0]);

    // NOTE: There might be some issues here with the ILU(0) if the column indices aren't sorted.
    // The other factorizations are probably OK.

    for (int j = 0 ; j < NumEntries ; ++j)
    {
      if (Indices[j] < n) // skip non-local columns
      {
        a[count]  = Values[j];
        ja[count] = Indices[j] + 1; // SPARSKIT is FORTRAN
        ++count;
        ++NumMyEntries;
      }
    }
    ia[i + 1] = ia[i] + NumMyEntries;
  }

  if (mbloc_ == -1) mbloc_ = n;

  int iwk;

  if (Type_ == "ILUT" || Type_ == "ILUTP" || Type_ == "ILUD" ||
      Type_ == "ILUDP")
    iwk = nnz + 2 * lfil_ * n;
  else if (Type_ == "ILUK")
    iwk = (2 * lfil_ + 1) * nnz + 1;
  else if (Type_ == "ILU0")
    iwk = nnz+2;

  int ierr = 0;

  alu_.resize(iwk);
  jlu_.resize(iwk);
  ju_.resize(n + 1);

  std::vector<int>    jw(n + 1);
  std::vector<double> w(n + 1);

  if (Type_ == "ILUT")
  {
    jw.resize(2 * n);
    F77_ILUT(&n, &a[0], &ja[0], &ia[0], &lfil_, &droptol_,
             &alu_[0], &jlu_[0], &ju_[0], &iwk, &w[0], &jw[0], &ierr);
  }
  else if (Type_ == "ILUD")
  {
    jw.resize(2 * n);
    F77_ILUD(&n, &a[0], &ja[0], &ia[0], &alph_, &tol_,
             &alu_[0], &jlu_[0], &ju_[0], &iwk, &w[0], &jw[0], &ierr);
  }
  else if (Type_ == "ILUTP")
  {
    jw.resize(2 * n);
    iperm_.resize(2 * n);
    F77_ILUTP(&n, &a[0], &ja[0], &ia[0], &lfil_, &droptol_, &permtol_, 
              &mbloc_, &alu_[0], &jlu_[0], &ju_[0], &iwk, &w[0], &jw[0],
              &iperm_[0], &ierr);
    for (int i = 0 ; i < n ; ++i)
      iperm_[i]--;
  }
  else if (Type_ == "ILUDP")
  {
    jw.resize(2 * n);
    iperm_.resize(2 * n);
    F77_ILUDP(&n, &a[0], &ja[0], &ia[0], &alph_, &droptol_, &permtol_, 
              &mbloc_, &alu_[0], &jlu_[0], &ju_[0], &n, &w[0], &jw[0],
              &iperm_[0], &ierr);
    for (int i = 0 ; i < n ; ++i)
      iperm_[i]--;
  }
  else if (Type_ == "ILUK")
  {
    std::vector<int> levs(iwk);
    jw.resize(3 * n);
    F77_ILUK(&n, &a[0], &ja[0], &ia[0], &lfil_, 
             &alu_[0], &jlu_[0], &ju_[0], &levs[0], &iwk, &w[0], &jw[0], &ierr);
  }
  else if (Type_ == "ILU0")
  {
    // here w is only of size n
    jw.resize(2 * n);
    F77_ILU0(&n, &a[0], &ja[0], &ia[0], 
             &alu_[0], &jlu_[0], &ju_[0], &jw[0], &ierr);
  }
  IFPACK_CHK_ERR(ierr);

  IsComputed_ = true;
  return(0);
}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_SPARSKIT::ApplyInverse(const Epetra_MultiVector& X, 
                                  Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-3); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK_CHK_ERR(-2); // Return error: X and Y not the same size

  int n = Matrix().NumMyRows();

  for (int i = 0 ; i < X.NumVectors() ; ++i)
    F77_LUSOL(&n, (double*)X(i)->Values(), Y(i)->Values(), (double*)&alu_[0], 
                (int*)&jlu_[0], (int*)&ju_[0]);

  // still need to fix support for permutation
  if (Type_ == "ILUTP" || Type_ == "ILUDP")
  {
    std::vector<double> tmp(n);
    for (int j = 0 ; j < n ; ++j)
      tmp[iperm_[j]] = Y[0][j];
    for (int j = 0 ; j < n ; ++j)
      Y[0][j] = tmp[j];
  }

  ++NumApplyInverse_;
  return(0);

}

//=============================================================================
double Ifpack_SPARSKIT::Condest(const Ifpack_CondestType CT, 
                                const int MaxIters, const double Tol,
                                Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_SPARSKIT::Print(std::ostream& os) const
{
  using std::endl;
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_SPARSKIT: " << Label() << endl << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << A_.NumGlobalRows() << endl;
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

#endif // HAVE_IFPACK_SPARSKIT
