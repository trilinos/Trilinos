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
#include "Tifpack_Condest.hpp"
#include "Tifpack_DenseRow.hpp"
#include "Tifpack_ScalingType.hpp"
#include "Tifpack_Parameters.hpp"
#include "Tpetra_CrsMatrix.hpp"
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
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  // @{ Constructors and Destructors
  //! ILUT constuctor with RowMatrix input
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
  void SetParameters(Teuchos::ParameterList& params);

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

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  // Mathematical functions.
  
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Computed the estimated condition number and returns the value.
  //! Returns the result of a ILUT forward/back solve on a Tpetra::MultiVector X in Y.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  void applyInverse(
      const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
            Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  magnitudeType Condest(const CondestType CT = Cheap, 
                 const LocalOrdinal MaxIters = 1550,
                 const magnitudeType Tol = 1e-9,
        Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  magnitudeType Condest() const
  {
    return(Condest_);
  }

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the apply() and applyInverse() methods.  If the implementation of this interface 
      does not support transpose use, this method should return a value of -1.
      
     \param
     UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

     \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);}

  //! Returns 0.0 because this class cannot compute Inf-norm.
  magnitudeType NormInf() const {return(0.0);}

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);}

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);}

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  const Teuchos::Comm<int> & Comm() const{return(*Comm_);}

  //! Returns a reference to the matrix to be preconditioned.
  const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Matrix() const
  {
    return(*A_);
  }

  //! Returns a reference to the L factor.
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & L() const {return(*L_);}
  
  //! Returns a reference to the U factor.
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & U() const {return(*U_);}
    
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

  //! Returns the number of calls to applyInverse().
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

  //! Returns the time spent in applyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  inline double LevelOfFill() const {
    return(LevelOfFill_);
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

  //! Get the relax value
  inline magnitudeType RelaxValue() const
  {
    return(RelaxValue_);
  }

  //! Gets the dropping tolerance
  inline magnitudeType DropTolerance() const
  {
    return(DropTolerance_);
  }
    
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {
    // FIXME: diagonal of L_ should not be stored
    return(L().getGlobalNumEntries() + U().getGlobalNumEntries() - L().getGlobalNumRows());
  }
 
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {
    return(L().getNodeNumEntries() + U().getNodeNumEntries());
  }

private:
  
  // @}
  // @{ Internal methods

  //! Copy constructor (should never be used)
  ILUT(const ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS);

  //! operator= (should never be used)
  ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS);

  //! Releases all allocated memory.
  void Destroy();

  // @}
  // @{ Internal data

  //! reference to the matrix to be preconditioned.
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Reference to the communicator object.
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! L factor
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > L_;
  //! U factor
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > U_;
  //! Condition number estimate.
  magnitudeType Condest_;
  //! Absolute threshold
  double Athresh_;
  //! Relative threshold
  double Rthresh_;
  magnitudeType RelaxValue_;
  //! Level-of-fill
  double LevelOfFill_;
  //! Discards all elements below this tolerance
  magnitudeType DropTolerance_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;
  //! \c true if transpose has to be used.
  bool UseTranspose_;
  //! Number of local rows.
  GlobalOrdinal NumMyRows_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to applyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to applyInverse().
  mutable double ApplyInverseTime_;
  //! Used for timing purposes
  mutable Teuchos::Time Time_;
  //! Global number of nonzeros in L and U factors
  size_t GlobalNonzeros_;
}; // ILUT

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ILUT(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A) :
  A_(A),
  Comm_(A->getRowMap()->getComm()),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  RelaxValue_(0.0),
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
  Time_("Tifpack::ILUT"),
  GlobalNonzeros_(0)
{
  // do nothing here..
}

//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~ILUT()
{
  Destroy();
}

//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameters(Teuchos::ParameterList& params)
{
  Tifpack::GetParameter(params, "fact: ilut level-of-fill", LevelOfFill_);
  TEST_FOR_EXCEPTION(LevelOfFill_ <= 0.0, std::runtime_error,
    "Tifpack::ILUT::SetParameters ERROR, level-of-fill must be >= 0.");

  Tifpack::GetParameter(params, "fact: absolute threshold", Athresh_);
  Tifpack::GetParameter(params, "fact: relative threshold", Rthresh_);
  Tifpack::GetParameter(params, "fact: relax value", RelaxValue_);
  Tifpack::GetParameter(params, "fact: drop tolerance", DropTolerance_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Initialize()
{
  // delete previously allocated factorization
  Destroy();

  Time_.start(true);

  // check only in serial
  TEST_FOR_EXCEPTION(Comm().getSize() == 1 && Matrix().getNodeNumRows() != Matrix().getNodeNumCols(), std::runtime_error, "Tifpack::ILUT::Initialize ERROR, matrix must be square");

  NumMyRows_ = Matrix().getNodeNumRows();

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  Time_.stop();
  InitializeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Compute()
{
  if (!IsInitialized()) {
    Initialize();
  }

  Time_.start(true);
  IsComputed_ = false;

  NumMyRows_ = A_->getNodeNumRows();
  size_t Length = A_->getNodeMaxNumRowEntries();
  Teuchos::Array<LocalOrdinal>    RowIndicesL(Length);
  Teuchos::Array<Scalar> RowValuesL(Length);
  Teuchos::Array<LocalOrdinal>    RowIndicesU(Length);
  Teuchos::Array<Scalar> RowValuesU(Length);

  size_t RowNnzU;

  L_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap(), A_->getColMap(), 0));
  U_ = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap(), A_->getColMap(), 0));

  TEST_FOR_EXCEPTION(L_ == Teuchos::null || U_ == Teuchos::null, std::runtime_error,
     "Tifpack::ILUT::Compute ERROR, failed to allocate L_ or U_");

  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar one  = Teuchos::ScalarTraits<Scalar>::one();

  Teuchos::Array<Scalar> DiagU(NumMyRows_,zero);

  Teuchos::Array<Teuchos::Array<LocalOrdinal> > Uindices(NumMyRows_);
  Teuchos::Array<Teuchos::Array<Scalar> > Ucoefs(NumMyRows_);

  // insert first row in U_ and L_
  LocalOrdinal lrow0 = 0;
  A_->getLocalRowCopy(lrow0, RowIndicesU(), RowValuesU(), RowNnzU);

  int count = 0;
  for (size_t i = 0 ;i < RowNnzU ; ++i) {
    if (RowIndicesU[i] < NumMyRows_){
      if (RowIndicesU[i] == 0) {
        Scalar& v = RowValuesU[i];
        v = AbsoluteThreshold() * TIFPACK_SGN(v) + RelativeThreshold() * v;
        DiagU[0] = v;
      }
      Uindices[lrow0].push_back(RowIndicesU[i]);
      Ucoefs[lrow0].push_back(RowValuesU[i]);
      ++count;
    }
  }

  RowValuesU[0] = one;
  RowIndicesU[0] = 0;
  L_->insertLocalValues(0,RowIndicesU(0,1), RowValuesU(0,1));

  LocalOrdinal max_col = NumMyRows_;

  Teuchos::Array<magnitudeType> AbsRow;

  Teuchos::Array<LocalOrdinal> tmp_idx;
  Teuchos::Array<Scalar> tmpv;

  // =================== //
  // start factorization //
  // =================== //

  for (LocalOrdinal row_i = 1 ; row_i < NumMyRows_ ; ++row_i)
  {
    Teuchos::ArrayRCP<const LocalOrdinal> ColIndices;
    Teuchos::ArrayRCP<const Scalar> ColValues;

    A_->getLocalRowView(row_i, ColIndices, ColValues);
    RowNnzU = ColIndices.size();

    Tifpack::DenseRow<Scalar,LocalOrdinal> SingleRowU(max_col);

    size_t NnzLower = 0;
    size_t NnzUpper = 0;

    LocalOrdinal start_col = ColIndices[0];

    size_t count = 0;
    for(size_t i=0; i<RowNnzU; ++i) {
      if (ColIndices[i] < NumMyRows_) {
        Scalar v = ColValues[i];
        if (ColIndices[i] < row_i)
          NnzLower++;
        else if (ColIndices[i] == row_i) {
          // add threshold
          NnzUpper++;
          v = AbsoluteThreshold() * TIFPACK_SGN(v) + RelativeThreshold() * v;
        }
        else
          NnzUpper++;

        if (ColIndices[i] < start_col) start_col = ColIndices[i];

        SingleRowU[ColIndices[i]] = v;
        ++count;
      }
    }
    RowNnzU = count;

    //Can magnitudeType always be assign-converted to double?
    double lof = LevelOfFill();

    int FillL = (int)(lof * NnzLower);
    int FillU = (int)(lof * NnzUpper);
    if (FillL == 0) FillL = 1;
    if (FillU == 0) FillU = 1;

    // for the multipliers
    Tifpack::DenseRow<Scalar,LocalOrdinal> SingleRowL(max_col);

    for (LocalOrdinal col_k = start_col ; col_k < row_i ; ++col_k) {

      Scalar xxx = SingleRowU[col_k];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(xxx) <= DropTolerance()) continue;

      Teuchos::Array<LocalOrdinal>& ColIndices = Uindices[col_k];
      Teuchos::Array<Scalar>& ColValues = Ucoefs[col_k];
      size_t ColNnzK = ColIndices.size();

      Scalar DiagonalValueK = DiagU[col_k];

      // FIXME: this should never happen!
      if (DiagonalValueK == zero)
        DiagonalValueK = AbsoluteThreshold();

      Scalar yyy = xxx / DiagonalValueK;

      if (yyy !=  zero) {
        SingleRowL[col_k] = yyy;
        for (size_t j = 0 ; j < ColNnzK ; ++j) {
          SingleRowU[ColIndices[j]] += -yyy * ColValues[j];
        }
      }
    }


    // drop elements to satisfy LevelOfFill(), start with L

    SingleRowL.reduceToArraysL( AbsRow, FillL, DropTolerance(), row_i, tmp_idx, tmpv );

    tmp_idx.push_back(row_i);
    tmpv.push_back(one);

    L_->insertLocalValues(row_i, tmp_idx(), tmpv());

    // same business with U_

    DiagU[row_i] =
     SingleRowU.reduceToArraysU(AbsRow, FillU, DropTolerance(), row_i, Uindices[row_i], Ucoefs[row_i]);
  }

  for(LocalOrdinal i=0; i<NumMyRows_; ++i) {
    U_->insertLocalValues(i, Uindices[i](), Ucoefs[i]() );
  }

  L_->fillComplete();
  U_->fillComplete();

  size_t MyNonzeros = L_->getGlobalNumEntries() + U_->getGlobalNumEntries();
  Teuchos::reduceAll(Comm(), Teuchos::REDUCE_SUM, 1, &MyNonzeros, &GlobalNonzeros_);

  IsComputed_ = true;

  ++NumCompute_;
  Time_.stop();
  ComputeTime_ += Time_.totalElapsedTime();
}

//=============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse(
           const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                 Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                 Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(!IsComputed(), std::runtime_error,
    "Tifpack::ILUT::applyInverse ERROR, Computed() hasn't been called yet.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Tifpack::ILUT::applyInverse ERROR, X.getNumVectors() != Y.getNumVectors().");

  Time_.start(true);

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xcopy;
  if (&(X.get2dView()[0][0]) == &(Y.get2dView()[0][0]))
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (mode == Teuchos::NO_TRANS)
  {
    // solves LU Y = X
    L_->solve(*Xcopy,Y,mode);
    U_->solve(Y,Y,mode);
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    U_->solve(*Xcopy,Y,mode);
    L_->solve(Y,Y,mode);
  }

  ++NumApplyInverse_;
  Time_.stop();
  ApplyInverseTime_ += Time_.totalElapsedTime();
}
//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(
     const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
           Teuchos::ETransp mode) const
{
  throw std::runtime_error("Tifpack::ILUT::apply is not implemented. Do you mean 'applyInverse' instead?");
}

//=============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
typename ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::magnitudeType
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Condest(const Tifpack::CondestType CT,
                     const LocalOrdinal MaxIters, const magnitudeType Tol,
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
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
std::ostream&
ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Print(std::ostream& os) const
{
  if (!Comm().getRank()) {
    os << std::endl;
    os << "================================================================================" << std::endl;
    os << "ILUT: " << std::endl;
    os << "Level-of-fill      = " << LevelOfFill() << std::endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << std::endl;
    os << "Relative threshold = " << RelativeThreshold() << std::endl;
    os << "Relax value        = " << RelaxValue() << std::endl;
    os << "Condition number estimate       = " << Condest() << std::endl;
    os << "Global number of rows           = " << A_->getGlobalNumRows() << std::endl;
    if (IsComputed_) {
      os << "Number of nonzeros in A         = " << A_->getGlobalNumEntries() << std::endl;
      os << "Number of nonzeros in L + U     = " << NumGlobalNonzeros()
         << " ( = " << 100.0 * NumGlobalNonzeros() / A_->getGlobalNumEntries()
         << " % of A)" << std::endl;
      os << "nonzeros / rows                 = "
        << 1.0 * NumGlobalNonzeros() / U_->getGlobalNumRows() << std::endl;
    }
    os << std::endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << std::endl;
    os << "-----           -------   --------------       ------------     --------" << std::endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize()
       << "  " << std::setw(15) << InitializeTime()
       << "               0.0            0.0" << std::endl;
    os << "Compute()       "   << std::setw(5) << NumCompute()
       << "  " << std::setw(15) << ComputeTime();
    os << "================================================================================" << std::endl;
    os << std::endl;
  }

  return(os);
}

}//namespace Tifpack

#endif /* TIFPACK_ILUT_HPP */
