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

#ifndef TIFPACK_ILUT_HPP
#define TIFPACK_ILUT_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_CondestType.hpp"
#include "Tifpack_ScalingType.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Tifpack {

//! Tifpack::ILUT: A class for constructing and using an ILUT factorization
// of a given Tpetra::RowMatrix.

/*! The Tifpack::ILUT class computes a "Relaxed" ILUT factorization with level k fill 
    of a given Tpetra::RowMatrix. 

    <P> Please refer to \ref ifp_ilu for a general description of the ILU algorithm.

    <P>The complete list of supported parameters is reported in page \ref ifp_params. 

    \author Michael Heroux, SNL 9214.

    \date Last modified on 22-Jan-05.
*/    
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
class ILUT: public Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
      
public:
  // @{ Constructors and Destructors
  //! ILUT constuctor with variable number of indices per row.
  ILUT(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);
  
  //! ILUT Destructor
  virtual ~ILUT();

  // @}
  // @{ Construction methods
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes five parameter names: level_fill, drop_tolerance,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive. For level_fill the ParameterEntry must have type int, the 
     threshold entries must have type double and overlap_mode must have type
     Tpetra::CombineMode.
  */
  int SetParameters(Teuchos::ParameterList& parameterlis);

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
    \param In 
           A - User matrix to be factored.
    \warning The graph of A must be identical to the graph passed in to IlukGraph constructor.
             
   */
  void Initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Compute IC factor U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the RILU(k) factors L and U using the current:
    <ol>
    <li> IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  void Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const {return(IsComputed_);};

  // Mathematical functions.
  
  //! Returns the result of a ILUT forward/back solve on a Tpetra::MultiVector X in Y.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  void applyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Computed the estimated condition number and returns the value.
  double Condest(const CondestType CT = Cheap, 
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
		 Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  double Condest() const
  {
    return(Condest_);
  }

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
      does not support transpose use, this method should return a value of -1.
      
     \param
     UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

     \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & getDomainMap() const {return(A_.getDomainMap());};

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Tpetra::Map<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getRangeMap() const{return(A_.getRangeMap());};

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  const Teuchos::Comm & Comm() const{return(Comm_);};

  //! Returns a reference to the matrix to be preconditioned.
  const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Matrix() const
  {
    return(A_);
  }

  //! Returns a reference to the L factor.
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & L() const {return(*L_);};
  
  //! Returns a reference to the U factor.
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & U() const {return(*U_);};
    
  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  inline double LevelOfFill() const {
    return(LevelOfFill_);
  }

  //! Set relative threshold value
  inline double RelaxValue() const {
    return(Relax_);
  }

  //! Get absolute threshold value
  inline double AbsoluteThreshold() const
  {
    return(Athresh_);
  }

  //! Get relative threshold value
  inline double RelativeThreshold() const
  {
    return(Rthresh_);
  }
    
  //! Gets the dropping tolerance
  inline double DropTolerance() const
  {
    return(DropTolerance_);
  }
    
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {
    // FIXME: diagonal of L_ should not be stored
    return(L().NumGlobalNonzeros() + U().NumGlobalNonzeros() - L().NumGlobalRows());
  }
 
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {
    return(L().NumMyNonzeros() + U().NumMyNonzeros());
  }

private:
  
  // @}
  // @{ Internal methods

  //! Copy constructor (should never be used)
  ILUT(const ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) :
    A_(RHS.Matrix()),
    Comm_(RHS.Comm()),
    Time_(Comm());

  //! operator= (should never be used)
  ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS);

  //! Releases all allocated memory.
  void Destroy();

  // @}
  // @{ Internal data

  //! reference to the matrix to be preconditioned.
  const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A_;
  //! Reference to the communicator object.
  const Teuchos::Comm& Comm_;
  //! L factor
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > L_;
  //! U factor
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > U_;
  //! Condition number estimate.
  double Condest_;
  //! relaxation value
  double Relax_;
  //! Absolute threshold
  double Athresh_;
  //! Relative threshold
  double Rthresh_;
  //! Level-of-fill
  double LevelOfFill_;
  //! Discards all elements below this tolerance
  double DropTolerance_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;
  //! \c true if transpose has to be used.
  bool UseTranspose_;
  //! Number of local rows.
  int NumMyRows_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  //! Used for timing purposes
  mutable Teuchos::Time Time_;
  //! Global number of nonzeros in L and U factors
  int GlobalNonzeros_;
  Teuchos::RCP<Tpetra::SerialComm> SerialComm_;
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > SerialMap_;
}; // ILUT

ILUT::ILUT(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* A) :
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
  GlobalNonzeros_(0)
{
  // do nothing here..
}

//==============================================================================
ILUT::~ILUT()
{
  Destroy();
}

//==============================================================================
void ILUT::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int ILUT::SetParameters(Teuchos::ParameterList& List)
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

    Label_ = "TIFPACK ILUT (fill=" + Tifpack_toString(LevelOfFill())
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
int ILUT::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.ResetStartTime();

  // check only in serial
  if (Comm().NumProc() == 1 && Matrix().NumMyRows() != Matrix().NumMyCols())
    TIFPACK_CHK_ERR(-2);

  NumMyRows_ = Matrix().NumMyRows();

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
// MS // completely rewritten the algorithm on 25-Jan-05, using more STL
// MS // and in a nicer and cleaner way. Also, it is more efficient.
// MS // WARNING: Still not fully tested!
int ILUT::Compute()
{
  if (!IsInitialized())
    TIFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  NumMyRows_ = A_.NumMyRows();
  int Length = A_.MaxNumEntries();
  vector<int>    RowIndicesL(Length);
  vector<double> RowValuesL(Length);
  vector<int>    RowIndicesU(Length);
  vector<double> RowValuesU(Length);
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

  int RowNnzU;

  L_ = rcp(new Tpetra_CrsMatrix(Copy, *SerialMap_, 0));
  U_ = rcp(new Tpetra_CrsMatrix(Copy, *SerialMap_, 0));

  if ((L_.get() == 0) || (U_.get() == 0))
    TIFPACK_CHK_ERR(-5); // memory allocation error

  // insert first row in U_ and L_
  TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(0,Length,RowNnzU,
                                     &RowValuesU[0],&RowIndicesU[0]));

  if (distributed)
  {
    int count = 0;
    for (int i = 0 ;i < RowNnzU ; ++i)
    {
      if (RowIndicesU[i] < NumMyRows_){
        RowIndicesU[count] = RowIndicesU[i];
        RowValuesU[count] = RowValuesU[i];
        ++count;
      }
      else
        continue;
    }
    RowNnzU = count;
  }

  // modify diagonal
  for (int i = 0 ;i < RowNnzU ; ++i) {
    if (RowIndicesU[i] == 0) {
      double& v = RowValuesU[i];
      v = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      break;
    }
  }

  TIFPACK_CHK_ERR(U_->InsertGlobalValues(0,RowNnzU,&(RowValuesU[0]),
                                        &(RowIndicesU[0])));
   // FIXME: DOES IT WORK IN PARALLEL ??
  RowValuesU[0] = 1.0;
  RowIndicesU[0] = 0;
  TIFPACK_CHK_ERR(L_->InsertGlobalValues(0,1,&(RowValuesU[0]),
                                        &(RowIndicesU[0])));

  int hash_size = 128;
  while (hash_size < (int) 1.5 * A_.MaxNumEntries() * LevelOfFill())
    hash_size *= 2;

  Tifpack_HashTable SingleRowU(hash_size - 1, 1);
  Tifpack_HashTable SingleRowL(hash_size - 1, 1);

  vector<int> keys;      keys.reserve(hash_size * 10);
  vector<double> values; values.reserve(hash_size * 10);
  vector<double> AbsRow; AbsRow.reserve(hash_size * 10);

  // =================== //
  // start factorization //
  // =================== //

  double this_proc_flops = 0.0;

  for (int row_i = 1 ; row_i < NumMyRows_ ; ++row_i)
  {
    // get row `row_i' of the matrix, store in U pointers
    TIFPACK_CHK_ERR(A_.ExtractMyRowCopy(row_i,Length,RowNnzU,
                                       &RowValuesU[0],&RowIndicesU[0]));

    if (distributed)
    {
      int count = 0;
      for (int i = 0 ;i < RowNnzU ; ++i)
      {
        if (RowIndicesU[i] < NumMyRows_){
          RowIndicesU[count] = RowIndicesU[i];
          RowValuesU[count] = RowValuesU[i];
          ++count;
        }
        else
          continue;
      }
      RowNnzU = count;
    }

    int NnzLower = 0;
    int NnzUpper = 0;

    for (int i = 0 ;i < RowNnzU ; ++i) {
      if (RowIndicesU[i] < row_i)
        NnzLower++;
      else if (R wIndicesU[i] == row_i) {
        // add threshold
        NnzUpper++;
        double& v = RowValuesU[i];
        v = AbsoluteThreshold() * EPETRA_SGN(v) + RelativeThreshold() * v;
      }
      else
        NnzUpper++;
    }

    int FillL = (int)(LevelOfFill() * NnzLower);
    int FillU = (int)(LevelOfFill() * NnzUpper);
    if (FillL == 0) FillL = 1;
    if (FillU == 0) FillU = 1;

    // convert line `row_i' into STL map for fast access
    SingleRowU.reset();

    for (int i = 0 ; i < RowNnzU ; ++i) {
        SingleRowU.set(RowIndicesU[i], RowValuesU[i]);
    }

    // for the multipliers
    SingleRowL.reset();

    int start_col = NumMyRows_;
    for (int i = 0 ; i < RowNnzU ; ++i)
      start_col = EPETRA_MIN(start_col, RowIndicesU[i]);

    int flops = 0;

    for (int col_k = start_col ; col_k < row_i ; ++col_k) {

      int*    ColIndicesK;
      double* ColValuesK;
      int     ColNnzK;

      TIFPACK_CHK_ERR(U_->ExtractGlobalRowView(col_k, ColNnzK, ColValuesK,
                                              ColIndicesK));

      // FIXME: can keep trace of diagonals
      double DiagonalValueK = 0.0;
      for (int i = 0 ; i < ColNnzK ; ++i) {
        if (ColIndicesK[i] == col_k) {
          DiagonalValueK = ColValuesK[i];
          break;
        }
      }

      // FIXME: this should never happen!
      if (DiagonalValueK == 0.0)
        DiagonalValueK = AbsoluteThreshold();

      double xxx = SingleRowU.get(col_k);
      if (std::abs(xxx) > DropTolerance()) {
        SingleRowL.set(col_k, xxx / DiagonalValueK);
        ++flops;

        for (int j = 0 ; j < ColNnzK ; ++j) {
          int col_j = ColIndicesK[j];

          if (col_j < col_k) continue;

          double yyy = SingleRowL.get(col_k);
          if (yyy !=  0.0)
            SingleRowU.set(col_j, -yyy * ColValuesK[j], true);
          flops += 2;
        }
      }
    }

    this_proc_flops += 1.0 * flops;

    double cutoff = DropTolerance();
    double DiscardedElements = 0.0;
    int count;

    // drop elements to satisfy LevelOfFill(), start with L
    count = 0;
    int sizeL = SingleRowL.getNumEntries();
    keys.resize(sizeL);
    values.resize(sizeL);

    AbsRow.resize(sizeL);

    SingleRowL.arrayify(
      keys.size() ? &keys[0] : 0,
      values.size() ? &values[0] : 0
      );
    for (int i = 0; i < sizeL; ++i)
      if (std::abs(values[i]) > DropTolerance()) {
        AbsRow[count++] = std::abs(values[i]);
      }

    if (count > FillL) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillL, AbsRow.begin() + count,
                  greater<double>());
      cutoff = AbsRow[FillL];
    }

    for (int i = 0; i < sizeL; ++i) {
      if (std::abs(values[i]) >= cutoff) {
        TIFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &values[i], (int*)&keys[i]));
      }
      else
        DiscardedElements += values[i];
    }

    // FIXME: DOES IT WORK IN PARALLEL ???
    // add 1 to the diagonal
    double dtmp = 1.0;
    TIFPACK_CHK_ERR(L_->InsertGlobalValues(row_i,1, &dtmp, &row_i));

    // same business with U_
    count = 0;
    int sizeU = SingleRowU.getNumEntries();
    AbsRow.resize(sizeU + 1);
    keys.resize(sizeU + 1);
    values.resize(sizeU + 1);

    SingleRowU.arrayify(&keys[0], &values[0]);

    for (int i = 0; i < sizeU; ++i)
      if (keys[i] >= row_i && std::abs(values[i]) > DropTolerance())
      {
        AbsRow[count++] = std::abs(values[i]);
      }

    if (count > FillU) {
      nth_element(AbsRow.begin(), AbsRow.begin() + FillU, AbsRow.begin() + count,
                  greater<double>());
      cutoff = AbsRow[FillU];
    }

    // sets the factors in U_
    for (int i = 0; i < sizeU; ++i)
    {
      int col = keys[i];
      double val = values[i];

      if (col >= row_i) {
        if (std::abs(val) >= cutoff || row_i == col) {
          TIFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &val, &col));
        }
        else
          DiscardedElements += val;
      }
    }

    // FIXME: not so sure of that!
    if (RelaxValue() != 0.0) {
      DiscardedElements *= RelaxValue();
      TIFPACK_CHK_ERR(U_->InsertGlobalValues(row_i,1, &DiscardedElements,
                                            &row_i));
    }
  }

  double tf;
  Comm().SumAll(&this_proc_flops, &tf, 1);
  ComputeFlops_ += tf;

  TIFPACK_CHK_ERR(L_->FillComplete());
  TIFPACK_CHK_ERR(U_->FillComplete());

#if 0
  // to check the complete factorization
  Tpetra_Vector LHS(A_.RowMatrixRowMap());
  Tpetra_Vector RHS1(A_.RowMatrixRowMap());
  Tpetra_Vector RHS2(A_.RowMatrixRowMap());
  Tpetra_Vector RHS3(A_.RowMatrixRowMap());
  LHS.Random();

  cout << "A = " << A_.NumGlobalNonzeros() << endl;
  cout << "L = " << L_->NumGlobalNonzeros() << endl;
  cout << "U = " << U_->NumGlobalNonzeros() << endl;

  A_.Multiply(false,LHS,RHS1);
  U_->Multiply(false,LHS,RHS2);
  L_->Multiply(false,RHS2,RHS3);

  RHS1.Update(-1.0, RHS3, 1.0);
  double Norm;
  RHS1.Norm2(&Norm);
#endif

  int MyNonzeros = L_->NumGlobalNonzeros() + U_->NumGlobalNonzeros();
  Comm().SumAll(&MyNonzeros, &GlobalNonzeros_, 1);

  IsComputed_ = true;

  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
int ILUT::ApplyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
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
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (!UseTranspose_)
  {
    // solves LU Y = X
    TIFPACK_CHK_ERR(L_->Solve(false,false,false,*Xcopy,Y));
    TIFPACK_CHK_ERR(U_->Solve(true,false,false,Y,Y));
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    TIFPACK_CHK_ERR(U_->Solve(true,true,false,*Xcopy,Y));
    TIFPACK_CHK_ERR(L_->Solve(false,true,false,Y,Y));
  }

  ++NumApplyInverse_;
  ApplyInverseFlops_ += X.NumVectors() * 2 * GlobalNonzeros_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int ILUT::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
          Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  return(-98);
}

//=============================================================================
Scalar ILUT::Condest(const Tifpack::CondestType CT,
                     const LocalOrdinal MaxIters, const Scalar Tol,
          Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0)
    Condest_ = Tifpack::Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
ILUT::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "ILUT: " << Label() << endl << endl;
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

}//namespace Tifpack

#endif /* TIFPACK_ILUT_HPP */
