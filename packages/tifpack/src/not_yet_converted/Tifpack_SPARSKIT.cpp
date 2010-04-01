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
#ifdef HAVE_TIFPACK_SPARSKIT
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_SPARSKIT.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Utils.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Util.hpp"
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
Tifpack_SPARSKIT::Tifpack_SPARSKIT(Tpetra_RowMatrix* A) :
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
Tifpack_SPARSKIT::~Tifpack_SPARSKIT()
{
}

//==========================================================================
int Tifpack_SPARSKIT::SetParameters(Teuchos::ParameterList& List)
{
  lfil_    = List.get("fact: sparskit: lfil", lfil_);
  tol_     = List.get("fact: sparskit: tol", tol_);
  droptol_ = List.get("fact: sparskit: droptol", droptol_);
  permtol_ = List.get("fact: sparskit: permtol", permtol_);
  alph_    = List.get("fact: sparskit: alph", alph_);
  mbloc_   = List.get("fact: sparskit: mbloc", mbloc_);
  Type_    = List.get("fact: sparskit: type", Type_);

  // set label
  Label_ = "TIFPACK SPARSKIT (Type=" + Type_ + ", fill=" +
    Tifpack_toString(lfil_) + ")";

  return(0);
}

//==========================================================================
int Tifpack_SPARSKIT::Initialize()
{
  IsInitialized_ = true;
  IsComputed_    = false;
  return(0);
}

//==========================================================================
int Tifpack_SPARSKIT::Compute() 
{
  if (!IsInitialized()) 
    TIFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;

  // convert the matrix into SPARSKIT format. The matrix is then
  // free'd after method Compute() returns.

  // convert the matrix into CSR format. Note that nnz is an over-estimate,
  // since it contains the nonzeros corresponding to external nodes as well.
  int n   = Matrix().NumMyRows();
  int nnz = Matrix().NumMyNonzeros();

  vector<double> a(nnz);
  vector<int>    ja(nnz);
  vector<int>    ia(n + 1);

  const int MaxNumEntries = Matrix().MaxNumEntries();

  vector<double> Values(MaxNumEntries);
  vector<int>    Indices(MaxNumEntries);

  int count = 0;

  ia[0] = 1;
  for (int i = 0 ; i < n ; ++i)
  {
    int NumEntries;
    int NumMyEntries = 0;
    Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumEntries, &Values[0], 
                              &Indices[0]);

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
    iwk = nnz;

  int ierr = 0;

  alu_.resize(iwk);
  jlu_.resize(iwk);
  ju_.resize(n + 1);

  vector<int>    jw(n + 1);
  vector<double> w(n + 1);

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
    vector<int> levs(iwk);
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

  TIFPACK_CHK_ERR(ierr);

  IsComputed_ = true;
  return(0);
}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Tifpack_SPARSKIT::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                                  Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  if (!IsComputed())
    TIFPACK_CHK_ERR(-3); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors()) 
    TIFPACK_CHK_ERR(-2); // Return error: X and Y not the same size

  int n = Matrix().NumMyRows();

  for (int i = 0 ; i < X.NumVectors() ; ++i)
      F77_LUSOL(&n, (double*)Y(i)->Values(), X(i)->Values(), (double*)&alu_[0], 
                (int*)&jlu_[0], (int*)&ju_[0]);

#if 0 
  // still need to fix support for permutation
  if (Type_ == "ILUTP" || Type_ == "ILUDP")
  {
    vector<double> tmp(n);
    for (int j = 0 ; j < n ; ++j)
      tmp[j] = Y[0][iperm_[j]];
    for (int j = 0 ; j < n ; ++j)
      Y[0][j] = tmp[j];
  }
#endif

  ++NumApplyInverse_;
  return(0);

}

//=============================================================================
double Tifpack_SPARSKIT::Condest(const Tifpack_CondestType CT, 
                                const int MaxIters, const double Tol,
                                Tpetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Tifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//=============================================================================
std::ostream&
Tifpack_SPARSKIT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Tifpack_SPARSKIT: " << Label() << endl << endl;
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

#endif // HAVE_TIFPACK_SPARSKIT
