/*@HEADER
// ***********************************************************************
// 
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

//-----------------------------------------------------
// Ifpack2::ILUT is a translation of the Aztec ILUT
// implementation. The Aztec ILUT implementation was
// written by Ray Tuminaro.
// See notes below, in the Ifpack2::ILUT::Compute method.
// ABW.
//------------------------------------------------------

#ifndef IFPACK2_ILUT_DEF_HPP
#define IFPACK2_ILUT_DEF_HPP

namespace Ifpack2 {

//Definitions for the ILUT methods:

//==============================================================================
template <class MatrixType>
ILUT<MatrixType>::ILUT(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A) :
  A_(A),
  Comm_(A->getRowMap()->getComm()),
  Athresh_(0.0),
  Rthresh_(1.0),
  RelaxValue_(0.0),
  LevelOfFill_(1.0),
  DropTolerance_(1e-12),
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  Time_("Ifpack2::ILUT"),
  NumMyRows_(-1),
  NumGlobalNonzeros_(0)
{ 
  TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      Teuchos::typeName(*this) << "::ILUT(): input matrix reference was null.");
}

//==========================================================================
template <class MatrixType>
ILUT<MatrixType>::~ILUT() {
}

//==========================================================================
template <class MatrixType>
void ILUT<MatrixType>::setParameters(const Teuchos::ParameterList& params) {
  Ifpack2::getParameter(params, "fact: ilut level-of-fill", LevelOfFill_);
  TEST_FOR_EXCEPTION(LevelOfFill_ <= 0.0, std::runtime_error,
    "Ifpack2::ILUT::SetParameters ERROR, level-of-fill must be >= 0.");

  double tmp = -1;
  Ifpack2::getParameter(params, "fact: absolute threshold", tmp);
  if (tmp != -1) Athresh_ = tmp;
  tmp = -1;
  Ifpack2::getParameter(params, "fact: relative threshold", tmp);
  if (tmp != -1) Rthresh_ = tmp;
  tmp = -1;
  Ifpack2::getParameter(params, "fact: relax value", tmp);
  if (tmp != -1) RelaxValue_ = tmp;
  tmp = -1;
  Ifpack2::getParameter(params, "fact: drop tolerance", tmp);
  if (tmp != -1) DropTolerance_ = tmp;
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
ILUT<MatrixType>::getComm() const{
  return(Comm_);
}

//==========================================================================
template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
ILUT<MatrixType>::getMatrix() const {
  return(A_);
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
ILUT<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
ILUT<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==============================================================================
template <class MatrixType>
bool ILUT<MatrixType>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template <class MatrixType>
int ILUT<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template <class MatrixType>
int ILUT<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template <class MatrixType>
int ILUT<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template <class MatrixType>
double ILUT<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double ILUT<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double ILUT<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType>
global_size_t ILUT<MatrixType>::getGlobalNumEntries() const { 
  return(L_->getGlobalNumEntries() + U_->getGlobalNumEntries());
}

//==========================================================================
template<class MatrixType>
size_t ILUT<MatrixType>::getNodeNumEntries() const {
  return(L_->getNodeNumEntries() + U_->getNodeNumEntries());
}

//=============================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
ILUT<MatrixType>::computeCondEst(
                     CondestType CT,
                     LocalOrdinal MaxIters, 
                     magnitudeType Tol,
                     const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix) {
  if (!isComputed()) { // cannot compute right now
    return(-1.0);
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return(Condest_);
}

//==========================================================================
template<class MatrixType>
void ILUT<MatrixType>::initialize() {
  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;
  L_ = Teuchos::null;
  U_ = Teuchos::null;

  Time_.start(true);

  // check only in serial
  TEST_FOR_EXCEPTION(Comm_->getSize() == 1 && A_->getNodeNumRows() != A_->getNodeNumCols(), std::runtime_error, "Ifpack2::ILUT::Initialize ERROR, matrix must be square");

  NumMyRows_ = A_->getNodeNumRows();

  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  Time_.stop();
  InitializeTime_ += Time_.totalElapsedTime();
}

template<typename Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType scalar_mag(const Scalar& s)
{
  return Teuchos::ScalarTraits<Scalar>::magnitude(s);
}

//==========================================================================
template<class MatrixType>
void ILUT<MatrixType>::compute() {
  //--------------------------------------------------------------------------
  // Ifpack2::ILUT is a translation of the Aztec ILUT implementation. The Aztec
  // ILUT implementation was written by Ray Tuminaro.
  //
  // This isn't an exact translation of the Aztec ILUT algorithm, for the
  // following reasons:
  // 1. Minor differences result from the fact that Aztec factors a MSR format
  // matrix in place, while the code below factors an input CrsMatrix which
  // remains untouched and stores the resulting factors in separate L and U
  // CrsMatrix objects.
  // Also, the Aztec code begins by shifting the matrix pointers back
  // by one, and the pointer contents back by one, and then using 1-based
  // Fortran-style indexing in the algorithm. This Ifpack2 code uses C-style
  // 0-based indexing throughout.
  // 2. Aztec stores the inverse of the diagonal of U. This Ifpack2 code
  // stores the non-inverted diagonal in U.
  // The triangular solves (in Ifpack2::ILUT::apply()) are performed by
  // calling the Tpetra::CrsMatrix::solve method on the L and U objects, and
  // this requires U to contain the non-inverted diagonal.
  //
  // ABW.
  //--------------------------------------------------------------------------

  if (!isInitialized()) {
    initialize();
  }

  Time_.start(true);

  NumMyRows_ = A_->getNodeNumRows();

  L_ = Teuchos::rcp(new MatrixType(A_->getRowMap(), A_->getColMap(), 0));
  U_ = Teuchos::rcp(new MatrixType(A_->getRowMap(), A_->getColMap(), 0));

  TEST_FOR_EXCEPTION(L_ == Teuchos::null || U_ == Teuchos::null, std::runtime_error,
     "Ifpack2::ILUT::Compute ERROR, failed to allocate L_ or U_");

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();

  // CGB: note, this caching approach may not be necessary anymore
  // We will store ArrayView objects that are views of the rows of U, so that
  // we don't have to repeatedly retrieve the view for each row. These will
  // be populated row by row as the factorization proceeds.
  Teuchos::Array<Teuchos::ArrayView<const LocalOrdinal> > Uindices(NumMyRows_);
  Teuchos::Array<Teuchos::ArrayView<const Scalar> >       Ucoefs(NumMyRows_);

  // If this macro is defined, files containing the L and U factors
  // will be written. DON'T CHECK IN THE CODE WITH THIS MACRO ENABLED!!!
  // #define IFPACK2_WRITE_FACTORS
#ifdef IFPACK2_WRITE_FACTORS
  std::ofstream ofsL("L.tif.mtx", std::ios::out);
  std::ofstream ofsU("U.tif.mtx", std::ios::out);
#endif

  // Calculate how much fill will be allowed in addition to the space that
  // corresponds to the input matrix entries.
  double local_nnz = static_cast<double>(A_->getNodeNumEntries());
  double fill = ((getLevelOfFill()-1)*local_nnz)/(2*NumMyRows_);

  // std::ceil gives the smallest integer larger than the argument.
  // this may give a slightly different result than Aztec's fill value in
  // some cases.
  double fill_ceil=std::ceil(fill);

  typedef typename Teuchos::Array<LocalOrdinal>::size_type      Tsize_t;

  // Similarly to Aztec, we will allow the same amount of fill for each
  // row, half in L and half in U.
  Tsize_t fillL = static_cast<Tsize_t>(fill_ceil);
  Tsize_t fillU = static_cast<Tsize_t>(fill_ceil);

  Teuchos::Array<Scalar> InvDiagU(NumMyRows_,zero);

  Teuchos::Array<LocalOrdinal> tmp_idx;
  Teuchos::Array<Scalar> tmpv;

  enum { UNUSED, ORIG, FILL };
  LocalOrdinal max_col = NumMyRows_;

  Teuchos::Array<int> pattern(max_col, UNUSED);
  Teuchos::Array<Scalar> cur_row(max_col, zero);
  Teuchos::Array<magnitudeType> unorm(max_col);
  magnitudeType rownorm;
  Teuchos::Array<LocalOrdinal> L_cols_heap;
  Teuchos::Array<LocalOrdinal> U_cols;
  Teuchos::Array<LocalOrdinal> L_vals_heap;
  Teuchos::Array<LocalOrdinal> U_vals_heap;

  // A comparison object which will be used to create 'heaps' of indices
  // that are ordered according to the corresponding values in the
  // 'cur_row' array.
  greater_indirect<Scalar,LocalOrdinal> vals_comp(cur_row);

  // =================== //
  // start factorization //
  // =================== //

  for (LocalOrdinal row_i = 0 ; row_i < NumMyRows_ ; ++row_i) {
    Teuchos::ArrayView<const LocalOrdinal> ColIndicesA;
    Teuchos::ArrayView<const Scalar> ColValuesA;

    A_->getLocalRowView(row_i, ColIndicesA, ColValuesA);
    size_t RowNnz = ColIndicesA.size();

    // Always include the diagonal in the U factor. The value should get
    // set in the next loop below.
    U_cols.push_back(row_i);
    cur_row[row_i] = zero;
    pattern[row_i] = ORIG;

    Tsize_t L_cols_heaplen = 0;
    rownorm = (magnitudeType)0;
    for(size_t i=0; i<RowNnz; ++i) {
      if (ColIndicesA[i] < NumMyRows_) {
        if (ColIndicesA[i] < row_i) {
          add_to_heap(ColIndicesA[i], L_cols_heap, L_cols_heaplen);
        }
        else if (ColIndicesA[i] > row_i) {
          U_cols.push_back(ColIndicesA[i]);
        }

        cur_row[ColIndicesA[i]] = ColValuesA[i];
        pattern[ColIndicesA[i]] = ORIG;
        rownorm += scalar_mag(ColValuesA[i]);
      }
    }

    // Alter the diagonal according to the absolute-threshold and
    // relative-threshold values. If not set, those values default
    // to zero and one respectively.
    const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rthresh = getRelativeThreshold();
    const Scalar& v = cur_row[row_i];
    cur_row[row_i] = (Scalar)(getAbsoluteThreshold()*IFPACK2_SGN(v)) + rthresh*v;

    Tsize_t orig_U_len = U_cols.size();
    RowNnz = L_cols_heap.size() + orig_U_len;
    rownorm = getDropTolerance() * rownorm/RowNnz;

    // The following while loop corresponds to the 'L30' goto's in Aztec.
    Tsize_t L_vals_heaplen = 0;
    while(L_cols_heaplen > 0) {
      LocalOrdinal row_k = L_cols_heap.front();

      Scalar multiplier = cur_row[row_k] * InvDiagU[row_k];
      cur_row[row_k] = multiplier;
      magnitudeType mag_mult = scalar_mag(multiplier);
      if (mag_mult*unorm[row_k] < rownorm) {
        pattern[row_k] = UNUSED;
        rm_heap_root(L_cols_heap, L_cols_heaplen);
        continue;
      }
      if (pattern[row_k] != ORIG) {
        if (L_vals_heaplen < fillL) {
          add_to_heap(row_k, L_vals_heap, L_vals_heaplen, vals_comp);
        }
        else if (L_vals_heaplen==0 ||
                 mag_mult < scalar_mag(cur_row[L_vals_heap.front()])) {
          pattern[row_k] = UNUSED;
          rm_heap_root(L_cols_heap, L_cols_heaplen);
          continue;
        }
        else {
          pattern[L_vals_heap.front()] = UNUSED;
          rm_heap_root(L_vals_heap, L_vals_heaplen, vals_comp);
          add_to_heap(row_k, L_vals_heap, L_vals_heaplen, vals_comp);
        }
      }

      /* Reduce current row */

      Teuchos::ArrayView<const LocalOrdinal>& ColIndicesU = Uindices[row_k];
      Teuchos::ArrayView<const Scalar>& ColValuesU = Ucoefs[row_k];
      Tsize_t ColNnzU = ColIndicesU.size();

      for(Tsize_t j=0; j<ColNnzU; ++j) {
        if (ColIndicesU[j] > row_k) {
          Scalar tmp = multiplier * ColValuesU[j];
          LocalOrdinal col_j = ColIndicesU[j];
          if (pattern[col_j] != UNUSED) {
            cur_row[col_j] -= tmp;
          }
          else if (scalar_mag(tmp) > rownorm) {
            cur_row[col_j] = -tmp;
            pattern[col_j] = FILL;
            if (col_j > row_i) {
              U_cols.push_back(col_j);
            }
            else {
              add_to_heap(col_j, L_cols_heap, L_cols_heaplen);
            }
          }
        }
      }

      rm_heap_root(L_cols_heap, L_cols_heaplen);
    }//end of while(L_cols_heaplen) loop


    // Put indices and values for L into arrays and then into the L_ matrix.

    //   first, the original entries from the L section of A:
    for(Tsize_t i=0; i<ColIndicesA.size(); ++i) {
      if (ColIndicesA[i] < row_i) {
        tmp_idx.push_back(ColIndicesA[i]);
        tmpv.push_back(cur_row[ColIndicesA[i]]);
        pattern[ColIndicesA[i]] = UNUSED;
      }
    }

    //   next, the L entries resulting from fill:
    for(Tsize_t j=0; j<L_vals_heaplen; ++j) {
      tmp_idx.push_back(L_vals_heap[j]);
      tmpv.push_back(cur_row[L_vals_heap[j]]);
      pattern[L_vals_heap[j]] = UNUSED;
    }

    // L has a one on the diagonal, but we don't explicitly store it.
    // If we don't store it, then the Tpetra/Kokkos kernel which performs
    // the triangular solve can assume a unit diagonal, take a short-cut
    // and perform faster.

    L_->insertLocalValues(row_i, tmp_idx(), tmpv());
#ifdef IFPACK2_WRITE_FACTORS
    for(Tsize_t ii=0; ii<tmp_idx.size(); ++ii) {
      ofsL << row_i << " " << tmp_idx[ii] << " " << tmpv[ii] << std::endl;
    }
#endif

    tmp_idx.clear();
    tmpv.clear();

    // Pick out the diagonal element, store its reciprocal.
    if (cur_row[row_i] == zero) {
      std::cerr << "Ifpack2::ILUT::Compute: zero pivot encountered! Replacing with rownorm and continuing...(You may need to set the parameter 'fact: absolute threshold'.)" << std::endl;
      cur_row[row_i] = rownorm;
    }
    InvDiagU[row_i] = one / cur_row[row_i];

    // Non-inverted diagonal is stored for U:
    tmp_idx.push_back(row_i);
    tmpv.push_back(cur_row[row_i]);
    unorm[row_i] = scalar_mag(cur_row[row_i]);
    pattern[row_i] = UNUSED;

    // Now put indices and values for U into arrays and then into the U_ matrix.
    // The first entry in U_cols is the diagonal, which we just handled, so we'll
    // start our loop at j=1.

    Tsize_t U_vals_heaplen = 0;
    for(Tsize_t j=1; j<U_cols.size(); ++j) {
      LocalOrdinal col = U_cols[j];
      if (pattern[col] != ORIG) {
        if (U_vals_heaplen < fillU) {
          add_to_heap(col, U_vals_heap, U_vals_heaplen, vals_comp);
        }
        else if (U_vals_heaplen!=0 && scalar_mag(cur_row[col]) >
                 scalar_mag(cur_row[U_vals_heap.front()])) {
          rm_heap_root(U_vals_heap, U_vals_heaplen, vals_comp);
          add_to_heap(col, U_vals_heap, U_vals_heaplen, vals_comp);
        }
      }
      else {
        tmp_idx.push_back(col);
        tmpv.push_back(cur_row[col]);
        unorm[row_i] += scalar_mag(cur_row[col]);
      }
      pattern[col] = UNUSED;
    }

    for(Tsize_t j=0; j<U_vals_heaplen; ++j) {
      tmp_idx.push_back(U_vals_heap[j]);
      tmpv.push_back(cur_row[U_vals_heap[j]]);
      unorm[row_i] += scalar_mag(cur_row[U_vals_heap[j]]);
    }

    unorm[row_i] /= (orig_U_len + U_vals_heaplen);

    U_->insertLocalValues(row_i, tmp_idx(), tmpv() );
#ifdef IFPACK2_WRITE_FACTORS
    for(int ii=0; ii<tmp_idx.size(); ++ii) {
      ofsU <<row_i<< " " <<tmp_idx[ii]<< " " <<tmpv[ii]<< std::endl;
    }
#endif
    tmp_idx.clear();
    tmpv.clear();

    U_->getLocalRowView(row_i, Uindices[row_i], Ucoefs[row_i] );

    L_cols_heap.clear();
    U_cols.clear();
    L_vals_heap.clear();
    U_vals_heap.clear();
  } // end of for(row_i) loop

  L_->fillComplete();
  U_->fillComplete();

  global_size_t MyNonzeros = L_->getGlobalNumEntries() + U_->getGlobalNumEntries();
  Teuchos::reduceAll(*Comm_,Teuchos::REDUCE_SUM,1,&MyNonzeros,&NumGlobalNonzeros_);

  IsComputed_ = true;

  ++NumCompute_;
  Time_.stop();
  ComputeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template <class MatrixType>
void ILUT<MatrixType>::apply(
           const Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& X,
                 Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Y,
                 Teuchos::ETransp mode,
               typename MatrixType::scalar_type alpha,
               typename MatrixType::scalar_type beta) const
{
  TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::ILUT::apply() ERROR, compute() hasn't been called yet.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::ILUT::apply() ERROR, X.getNumVectors() != Y.getNumVectors().");

  Time_.start(true);

  //The following comment is from Ifpack's code. Leave it in Ifpack2 for now...
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
    L_->localSolve(*Xcopy,Y,Teuchos::NO_TRANS);
    U_->localSolve(Y,Y,Teuchos::NO_TRANS);
  }
  else
  {
    // solves U(trans) L(trans) Y = X
    U_->localSolve(*Xcopy,Y,mode);
    L_->localSolve(Y,Y,mode);
  }

  ++NumApply_;
  Time_.stop();
  ApplyTime_ += Time_.totalElapsedTime();
}

//=============================================================================
template <class MatrixType>
std::string ILUT<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//=============================================================================
template <class MatrixType>
void ILUT<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: 
  //    high: 
  // extreme: 
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << endl;
    out << "Level-of-fill      = " << getLevelOfFill()       << endl;
    out << "Absolute threshold = " << getAbsoluteThreshold() << endl;
    out << "Relative threshold = " << getRelativeThreshold() << endl;
    out << "Relax value        = " << getRelaxValue()        << endl;
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    if (isComputed()) {
      out << "Number of nonzeros in A         = " << A_->getGlobalNumEntries() << endl;
      out << "Number of nonzeros in L + U     = " << getGlobalNumEntries()
          << " ( = " << 100.0 * (double)getGlobalNumEntries() / (double)A_->getGlobalNumEntries() << " % of A)" << endl;
      out << "nonzeros / rows                 = " << 1.0 * getGlobalNumEntries() / U_->getGlobalNumRows() << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s) " << endl;
    out << "------------    -------    ---------------" << endl;
    out << "initialize()    " << setw(7) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << "compute()       " << setw(7) << getNumCompute()    << "    " << setw(15) << getComputeTime()    << endl;
    out << "apply()         " << setw(7) << getNumApply()      << "    " << setw(15) << getApplyTime()      << endl;
    out << "==============================================================================="                << endl;
    out << endl;
  }
}


}//namespace Ifpack2

#endif /* IFPACK2_ILUT_DEF_HPP */

