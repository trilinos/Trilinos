/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_ILUT_DEF_HPP
#define IFPACK2_ILUT_DEF_HPP

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#include <Ifpack2_Heap.hpp>
#include <Ifpack2_Condest.hpp>
#include <Ifpack2_LocalFilter.hpp>
#include <Ifpack2_Parameters.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#include <Ifpack2_ILUT.hpp>

#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>


namespace Ifpack2 {

  namespace {

    /// \brief Default drop tolerance for ILUT.
    ///
    /// \tparam ScalarType The "scalar type"; the type of entries in
    ///   the input sparse matrix to ILUT.  This is the same as the
    ///   scalar_type typedef of ILUT.
    ///
    /// \warning This is an implementation detail of Ifpack2.  Do NOT
    ///   depend on this function or use it in your code.  It may go
    ///   away entirely or change interface or behavior without
    ///   warning.
    ///
    /// This function preserves the previous default drop tolerance
    /// (1e-12, independent of scalar type), thus ensuring backwards
    /// compatibility for the common case of ScalarType=double.
    /// However, it provides a more reasonable default for other
    /// scalar types of possibly lower or higher precision than
    /// double.
    ///
    /// This function is templated on ScalarType, rather than its
    /// magnitude type, so that we can handle complex numbers
    /// specially if desired.
    ///
    /// In order to override the default, just specialize this
    /// function for your particular ScalarType.
    template<class ScalarType>
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    ilutDefaultDropTolerance () {
      using Teuchos::as;
      typedef Teuchos::ScalarTraits<ScalarType> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;

      // 1/2.  Hopefully this can be represented in magnitude_type.
      const magnitude_type oneHalf = STM::one() / (STM::one() + STM::one());

      // The min ensures that in case magnitude_type has very low
      // precision, we'll at least get some value strictly less than
      // one.
      return std::min (as<magnitude_type> (1000) * STS::magnitude (STS::eps ()), oneHalf);
    }

    // Full specialization for ScalarType = double.
    // This specialization preserves ILUT's previous default behavior.
    template<>
    Teuchos::ScalarTraits<double>::magnitudeType
    ilutDefaultDropTolerance<double> () {
      return 1e-12;
    }

  } // namespace (anonymous)


template <class MatrixType>
ILUT<MatrixType>::ILUT (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
  Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ()),
  RelaxValue_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
  LevelOfFill_ (1),
  DropTolerance_ (ilutDefaultDropTolerance<scalar_type> ()),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{}

template <class MatrixType>
ILUT<MatrixType>::~ILUT()
{}

template <class MatrixType>
void ILUT<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  int fillLevel = 1;
  magnitude_type absThresh = STM::zero ();
  magnitude_type relThresh = STM::one ();
  magnitude_type relaxValue = STM::zero ();
  magnitude_type dropTol = ilutDefaultDropTolerance<scalar_type> ();

  bool gotFillLevel = false;
  try {
    // Try getting the fill level as an int.
    fillLevel = params.get<int> ("fact: ilut level-of-fill");
    gotFillLevel = true;
  }
  catch (InvalidParameterName&) {
    // We didn't really get it, but this tells us to stop looking.
    gotFillLevel = true;
  }
  catch (InvalidParameterType&) {
    // The name is right, but the type is wrong; try different types.
    // We don't have to check InvalidParameterName again, since we
    // checked it above, and it has nothing to do with the type.
  }

  if (! gotFillLevel) {
    // Try magnitude_type, for backwards compatibility.
    try {
      fillLevel = as<int> (params.get<magnitude_type> ("fact: ilut level-of-fill"));
    }
    catch (InvalidParameterType&) {}
  }
  if (! gotFillLevel) {
    // Try double, for backwards compatibility.
    try {
      fillLevel = as<int> (params.get<double> ("fact: ilut level-of-fill"));
    }
    catch (InvalidParameterType&) {}
  }
  // If none of the above attempts succeed, accept the default value.

  TEUCHOS_TEST_FOR_EXCEPTION(
    fillLevel <= 0, std::runtime_error,
    "Ifpack2::ILUT: The \"fact: ilut level-of-fill\" parameter must be "
    "strictly greater than zero, but you specified a value of " << fillLevel
    << ".  Remember that for ILUT, the fill level p means something different "
    "than it does for ILU(k).  ILU(0) produces factors with the same sparsity "
    "structure as the input matrix A; ILUT with p = 0 always produces a "
    "diagonal matrix, and is thus probably not what you want.");

  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double> ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double> ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relaxValue = params.get<magnitude_type> ("fact: relax value");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relaxValue = as<magnitude_type> (params.get<double> ("fact: relax value"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    dropTol = params.get<magnitude_type> ("fact: drop tolerance");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    dropTol = as<magnitude_type> (params.get<double> ("fact: drop tolerance"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  // "Commit" the values only after validating all of them.  This
  // ensures that there are no side effects if this routine throws an
  // exception.

  // mfh 28 Nov 2012: The previous code would not assign Athresh_,
  // Rthresh_, RelaxValue_, or DropTolerance_ if the read-in value was
  // -1.  I don't know if keeping this behavior is correct, but I'll
  // keep it just so as not to change previous behavior.

  LevelOfFill_ = fillLevel;
  if (absThresh != -STM::one ()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one ()) {
    Rthresh_ = relThresh;
  }
  if (relaxValue != -STM::one ()) {
    RelaxValue_ = relaxValue;
  }
  if (dropTol != -STM::one ()) {
    DropTolerance_ = dropTol;
  }
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
ILUT<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const typename ILUT<MatrixType>::row_matrix_type>
ILUT<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const typename ILUT<MatrixType>::map_type>
ILUT<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const typename ILUT<MatrixType>::map_type>
ILUT<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool ILUT<MatrixType>::hasTransposeApply () const {
  return true;
}


template <class MatrixType>
int ILUT<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int ILUT<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int ILUT<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double ILUT<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template<class MatrixType>
double ILUT<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template<class MatrixType>
double ILUT<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}


template<class MatrixType>
global_size_t ILUT<MatrixType>::getGlobalNumEntries () const {
  return L_->getGlobalNumEntries () + U_->getGlobalNumEntries ();
}


template<class MatrixType>
size_t ILUT<MatrixType>::getNodeNumEntries () const {
  return L_->getNodeNumEntries () + U_->getNodeNumEntries ();
}


template<class MatrixType>
typename ILUT<MatrixType>::magnitude_type
ILUT<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) {
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template<class MatrixType>
void ILUT<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != A_.getRawPtr ()) {
    // Check in serial or one-process mode if the matrix is square.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! A.is_null () && A->getComm ()->getSize () == 1 &&
      A->getNodeNumRows () != A->getNodeNumCols (),
      std::runtime_error, "Ifpack2::ILUT::setMatrix: If A's communicator only "
      "contains one process, then A must be square.  Instead, you provided a "
      "matrix A with " << A->getNodeNumRows () << " rows and "
      << A->getNodeNumCols () << " columns.");

    // It's legal for A to be null; in that case, you may not call
    // initialize() until calling setMatrix() with a nonnull input.
    // Regardless, setting the matrix invalidates any previous
    // factorization.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;
    L_ = Teuchos::null;
    U_ = Teuchos::null;
    Condest_ = -Teuchos::ScalarTraits<magnitude_type>::one();
    A_ = A;
  }
}


template<class MatrixType>
void ILUT<MatrixType>::initialize ()
{
  Teuchos::Time timer ("ILUT::initialize");
  {
    Teuchos::TimeMonitor timeMon (timer);

    // Check that the matrix is nonnull.
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::ILUT::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;
    L_ = Teuchos::null;
    U_ = Teuchos::null;

    A_local_ = makeLocalFilter (A_); // Compute the local filter.

    IsInitialized_ = true;
    ++NumInitialize_;
  }
  InitializeTime_ += timer.totalElapsedTime ();
}


template<typename ScalarType>
typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
scalar_mag (const ScalarType& s)
{
  return Teuchos::ScalarTraits<ScalarType>::magnitude(s);
}


template<class MatrixType>
void ILUT<MatrixType>::compute ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::reduceAll;

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

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("ILUT::compute");
  { // Timer scope for timing compute()
    Teuchos::TimeMonitor timeMon (timer, true);
    const scalar_type zero = STS::zero ();
    const scalar_type one  = STS::one ();

    const local_ordinal_type myNumRows = A_local_->getNodeNumRows ();
    L_ = rcp (new crs_matrix_type (A_local_->getRowMap (), A_local_->getColMap (), 0));
    U_ = rcp (new crs_matrix_type (A_local_->getRowMap (), A_local_->getColMap (), 0));

    // CGB: note, this caching approach may not be necessary anymore
    // We will store ArrayView objects that are views of the rows of U, so that
    // we don't have to repeatedly retrieve the view for each row. These will
    // be populated row by row as the factorization proceeds.
    Array<ArrayView<const local_ordinal_type> > Uindices (myNumRows);
    Array<ArrayView<const scalar_type> >       Ucoefs (myNumRows);

    // If this macro is defined, files containing the L and U factors
    // will be written. DON'T CHECK IN THE CODE WITH THIS MACRO ENABLED!!!
    // #define IFPACK2_WRITE_FACTORS
#ifdef IFPACK2_WRITE_FACTORS
    std::ofstream ofsL("L.tif.mtx", std::ios::out);
    std::ofstream ofsU("U.tif.mtx", std::ios::out);
#endif

    // The code here uses double for fill calculations, even though
    // the fill level is actually an integer.  The point is to avoid
    // rounding and overflow for integer calculations.  If int is <=
    // 32 bits, it can never overflow double, so this cast is always
    // OK as long as int has <= 32 bits.

    // Calculate how much fill will be allowed in addition to the
    // space that corresponds to the input matrix entries.
    double local_nnz = static_cast<double> (A_local_->getNodeNumEntries ());
    double fill;
    {
      const double fillLevel = as<double> (getLevelOfFill ());
      fill = ((fillLevel - 1) * local_nnz) / (2 * myNumRows);
    }

    // std::ceil gives the smallest integer larger than the argument.
    // this may give a slightly different result than Aztec's fill value in
    // some cases.
    double fill_ceil=std::ceil(fill);

    // Similarly to Aztec, we will allow the same amount of fill for each
    // row, half in L and half in U.
    size_type fillL = static_cast<size_type>(fill_ceil);
    size_type fillU = static_cast<size_type>(fill_ceil);

    Array<scalar_type> InvDiagU (myNumRows, zero);

    Array<local_ordinal_type> tmp_idx;
    Array<scalar_type> tmpv;

    enum { UNUSED, ORIG, FILL };
    local_ordinal_type max_col = myNumRows;

    Array<int> pattern(max_col, UNUSED);
    Array<scalar_type> cur_row(max_col, zero);
    Array<magnitude_type> unorm(max_col);
    magnitude_type rownorm;
    Array<local_ordinal_type> L_cols_heap;
    Array<local_ordinal_type> U_cols;
    Array<local_ordinal_type> L_vals_heap;
    Array<local_ordinal_type> U_vals_heap;

    // A comparison object which will be used to create 'heaps' of indices
    // that are ordered according to the corresponding values in the
    // 'cur_row' array.
    greater_indirect<scalar_type,local_ordinal_type> vals_comp(cur_row);

    // =================== //
    // start factorization //
    // =================== //

    ArrayRCP<local_ordinal_type> ColIndicesARCP;
    ArrayRCP<scalar_type>       ColValuesARCP;
    if (! A_local_->supportsRowViews ()) {
      const size_t maxnz = A_local_->getNodeMaxNumRowEntries ();
      ColIndicesARCP.resize (maxnz);
      ColValuesARCP.resize (maxnz);
    }

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      ArrayView<const local_ordinal_type> ColIndicesA;
      ArrayView<const scalar_type> ColValuesA;
      size_t RowNnz;

      if (A_local_->supportsRowViews ()) {
        A_local_->getLocalRowView (row_i, ColIndicesA, ColValuesA);
        RowNnz = ColIndicesA.size ();
      }
      else {
        A_local_->getLocalRowCopy (row_i, ColIndicesARCP (), ColValuesARCP (), RowNnz);
        ColIndicesA = ColIndicesARCP (0, RowNnz);
        ColValuesA = ColValuesARCP (0, RowNnz);
      }

      // Always include the diagonal in the U factor. The value should get
      // set in the next loop below.
      U_cols.push_back(row_i);
      cur_row[row_i] = zero;
      pattern[row_i] = ORIG;

      size_type L_cols_heaplen = 0;
      rownorm = STM::zero ();
      for (size_t i = 0; i < RowNnz; ++i) {
        if (ColIndicesA[i] < myNumRows) {
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
      const magnitude_type rthresh = getRelativeThreshold();
      const scalar_type& v = cur_row[row_i];
      cur_row[row_i] = as<scalar_type> (getAbsoluteThreshold() * IFPACK2_SGN(v)) + rthresh*v;

      size_type orig_U_len = U_cols.size();
      RowNnz = L_cols_heap.size() + orig_U_len;
      rownorm = getDropTolerance() * rownorm/RowNnz;

      // The following while loop corresponds to the 'L30' goto's in Aztec.
      size_type L_vals_heaplen = 0;
      while (L_cols_heaplen > 0) {
        local_ordinal_type row_k = L_cols_heap.front();

        scalar_type multiplier = cur_row[row_k] * InvDiagU[row_k];
        cur_row[row_k] = multiplier;
        magnitude_type mag_mult = scalar_mag(multiplier);
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

        ArrayView<const local_ordinal_type>& ColIndicesU = Uindices[row_k];
        ArrayView<const scalar_type>& ColValuesU = Ucoefs[row_k];
        size_type ColNnzU = ColIndicesU.size();

        for(size_type j=0; j<ColNnzU; ++j) {
          if (ColIndicesU[j] > row_k) {
            scalar_type tmp = multiplier * ColValuesU[j];
            local_ordinal_type col_j = ColIndicesU[j];
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
      for (size_type i = 0; i < ColIndicesA.size (); ++i) {
        if (ColIndicesA[i] < row_i) {
          tmp_idx.push_back(ColIndicesA[i]);
          tmpv.push_back(cur_row[ColIndicesA[i]]);
          pattern[ColIndicesA[i]] = UNUSED;
        }
      }

      //   next, the L entries resulting from fill:
      for (size_type j = 0; j < L_vals_heaplen; ++j) {
        tmp_idx.push_back(L_vals_heap[j]);
        tmpv.push_back(cur_row[L_vals_heap[j]]);
        pattern[L_vals_heap[j]] = UNUSED;
      }

      // L has a one on the diagonal, but we don't explicitly store
      // it.  If we don't store it, then the kernel which performs the
      // triangular solve can assume a unit diagonal, take a short-cut
      // and perform faster.

      L_->insertLocalValues (row_i, tmp_idx (), tmpv ());
#ifdef IFPACK2_WRITE_FACTORS
      for (size_type ii = 0; ii < tmp_idx.size (); ++ii) {
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

      size_type U_vals_heaplen = 0;
      for(size_type j=1; j<U_cols.size(); ++j) {
        local_ordinal_type col = U_cols[j];
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

      for(size_type j=0; j<U_vals_heaplen; ++j) {
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

    // FIXME (mfh 03 Apr 2013) Do we need to supply a domain and range Map?
    L_->fillComplete();
    U_->fillComplete();
  }
  ComputeTime_ += timer.totalElapsedTime ();
  IsComputed_ = true;
  ++NumCompute_;
}


template <class MatrixType>
void ILUT<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  Teuchos::Time timer ("ILUT::apply");
  { // Timer scope for timing apply()
    Teuchos::TimeMonitor timeMon (timer, true);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::ILUT::apply: You must call compute() to compute the incomplete "
      "factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::ILUT::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    if (alpha == Teuchos::ScalarTraits<scalar_type>::zero ()) {
      // alpha == 0, so we don't need to apply the operator.
      //
      // The special case for beta == 0 ensures that if Y contains Inf
      // or NaN values, we replace them with 0 (following BLAS
      // convention), rather than multiplying them by 0 to get NaN.
      if (beta == Teuchos::ScalarTraits<scalar_type>::zero ()) {
        Y.putScalar (beta);
      } else {
        Y.scale (beta);
      }
      return;
    }

    // If beta != 0, create a temporary multivector Y_temp to hold the
    // contents of alpha*M^{-1}*X.  Otherwise, alias Y_temp to Y.
    RCP<MV> Y_temp;
    if (beta == Teuchos::ScalarTraits<scalar_type>::zero ()) {
      Y_temp = rcpFromRef (Y);
    } else {
      Y_temp = rcp (new MV (Y.getMap (), Y.getNumVectors ()));
    }

    // If X and Y are pointing to the same memory location, create an
    // auxiliary vector, X_temp, so that we don't clobber the input
    // when computing the output.  Otherwise, alias X_temp to X.
    RCP<const MV> X_temp;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      X_temp = rcp (new MV (X, Teuchos::Copy));
    } else {
      X_temp = rcpFromRef (X);
    }

    // Create a temporary multivector Y_mid to hold the intermediate
    // between the L and U (or U and L, for the transpose or conjugate
    // transpose case) solves.
    RCP<MV> Y_mid = rcp (new MV (Y.getMap (), Y.getNumVectors ()));

    if (mode == Teuchos::NO_TRANS) { // Solve L U Y = X
      L_->template localSolve<scalar_type, scalar_type> (*X_temp, *Y_mid, mode);
      // FIXME (mfh 20 Aug 2013) Is it OK to use Y_temp for both the
      // input and the output?
      U_->template localSolve<scalar_type, scalar_type> (*Y_mid, *Y_temp, mode);
    }
    else { // Solve U^* L^* Y = X
      U_->template localSolve<scalar_type, scalar_type> (*X_temp, *Y_mid, mode);
      // FIXME (mfh 20 Aug 2013) Is it OK to use Y_temp for both the
      // input and the output?
      L_->template localSolve<scalar_type, scalar_type> (*Y_mid, *Y_temp, mode);
    }

    if (beta == Teuchos::ScalarTraits<scalar_type>::zero ()) {
      Y.scale (alpha);
    } else { // beta != 0
      Y.update (alpha, *Y_temp, beta);
    }
  }
  ++NumApply_;
  ApplyTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
std::string ILUT<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::ILUT\": {";
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  os << "Level-of-fill: " << getLevelOfFill() << ", "
     << "absolute threshold: " << getAbsoluteThreshold() << ", "
     << "relative threshold: " << getRelativeThreshold() << ", "
     << "relaxation value: " << getRelaxValue() << ", ";

  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]"
       << ", Global nnz: " << A_->getGlobalNumEntries();
  }

  os << "}";
  return os.str ();
}


template <class MatrixType>
void
ILUT<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using std::endl;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
  OSTab tab0 (out);

  if (vl > VERB_NONE) {
    out << "\"Ifpack2::ILUT\":" << endl;
    OSTab tab1 (out);
    out << "MatrixType: " << TypeNameTraits<MatrixType>::name () << endl;
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
    }
    out << "Initialized: " << (isInitialized () ? "true" : "false")
        << endl
        << "Computed: " << (isComputed () ? "true" : "false")
        << endl
        << "Level of fill: " << getLevelOfFill () << endl
        << "Absolute threshold: " << getAbsoluteThreshold () << endl
        << "Relative threshold: " << getRelativeThreshold () << endl
        << "Relax value: " << getRelaxValue () << endl;

    if (isComputed () && vl >= VERB_HIGH) {
      const double fillFraction =
        (double) getGlobalNumEntries () / (double) A_->getGlobalNumEntries ();
      const double nnzToRows =
        (double) getGlobalNumEntries () / (double) U_->getGlobalNumRows ();

      out << "Dimensions of L: [" << L_->getGlobalNumRows () << ", "
          << L_->getGlobalNumRows () << "]" << endl
          << "Dimensions of U: [" << U_->getGlobalNumRows () << ", "
          << U_->getGlobalNumRows () << "]" << endl
          << "Number of nonzeros in factors: " << getGlobalNumEntries () << endl
          << "Fill fraction of factors over A: " << fillFraction << endl
          << "Ratio of nonzeros to rows: " << nnzToRows << endl;
    }

    out << "Number of initialize calls: " << getNumInitialize () << endl
        << "Number of compute calls: " << getNumCompute () << endl
        << "Number of apply calls: " << getNumApply () << endl
        << "Total time in seconds for initialize: " << getInitializeTime () << endl
        << "Total time in seconds for compute: " << getComputeTime () << endl
        << "Total time in seconds for apply: " << getApplyTime () << endl;

    out << "Local matrix:" << endl;
    A_local_->describe (out, vl);
  }
}

template <class MatrixType>
Teuchos::RCP<const typename ILUT<MatrixType>::row_matrix_type>
ILUT<MatrixType>::makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A->getComm ()->getSize () > 1) {
    return Teuchos::rcp (new LocalFilter<MatrixType> (A));
  } else {
    return A;
  }
}

}//namespace Ifpack2


// FIXME (mfh 16 Sep 2014) We should really only use RowMatrix here!
// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_ILUT_INSTANT(S,LO,GO,N) \
  template class Ifpack2::ILUT< Tpetra::RowMatrix<S, LO, GO, N> >; \
  template class Ifpack2::ILUT< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif /* IFPACK2_ILUT_DEF_HPP */
