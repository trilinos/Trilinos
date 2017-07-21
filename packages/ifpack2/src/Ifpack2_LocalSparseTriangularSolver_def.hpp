/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifdef HAVE_IFPACK2_SHYLUHTS
# include "shylu_hts.hpp"
#endif

namespace Ifpack2 {

namespace Details {
struct TrisolverType {
  enum Enum {
    Internal, //!< Tpetra::CrsMatrix::localSolve
    HTS       //!< Multicore ShyLU/HTS
  };

  static void loadPLTypeOption (Teuchos::Array<std::string>& type_strs, Teuchos::Array<Enum>& type_enums) {
    type_strs.resize(2);
    type_strs[0] = "Internal";
    type_strs[1] = "HTS";
    type_enums.resize(2);
    type_enums[0] = Internal;
    type_enums[1] = HTS;
  }
};
}

template<class MatrixType>
class LocalSparseTriangularSolver<MatrixType>::HtsImpl {
public:
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> crs_matrix_type;

  void reset () {
#ifdef HAVE_IFPACK2_SHYLUHTS
    Timpl_ = Teuchos::null;
    levelset_block_size_ = 1;
#endif
  }

  void setParameters (const Teuchos::ParameterList& pl) {
#ifdef HAVE_IFPACK2_SHYLUHTS
    const char* block_size_s = "trisolver: block size";
    if (pl.isParameter(block_size_s)) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG( ! pl.isType<int>(block_size_s),
                                   "The parameter \"" << block_size_s << "\" must be of type int.");
      levelset_block_size_ = pl.get<int>(block_size_s);
    }
    if (levelset_block_size_ < 1)
      levelset_block_size_ = 1;
#endif
  }

  // HTS has the phases symbolic+numeric, numeric, and apply. Hence the first
  // call to compute() will trigger the symbolic+numeric phase, and subsequent
  // calls (with the same Timpl_) will trigger the numeric phase. In the call to
  // initialize(), essentially nothing happens.
  void initialize (const crs_matrix_type& /* unused */) {
#ifdef HAVE_IFPACK2_SHYLUHTS
    reset();
    transpose_ = conjugate_ = false;
#endif
  }

  void compute (const crs_matrix_type& T_in, const Teuchos::RCP<Teuchos::FancyOStream>& out) {
#ifdef HAVE_IFPACK2_SHYLUHTS
    using Teuchos::ArrayRCP;

    Teuchos::ArrayRCP<const size_t> rowptr;
    Teuchos::ArrayRCP<const local_ordinal_type> colidx;
    Teuchos::ArrayRCP<const scalar_type> val;
    T_in.getAllValues(rowptr, colidx, val);

    Teuchos::RCP<HtsCrsMatrix> T_hts = Teuchos::rcpWithDealloc(
      HTST::make_CrsMatrix(rowptr.size() - 1,
                           rowptr.getRawPtr(), colidx.getRawPtr(), val.getRawPtr(),
                           transpose_, conjugate_),
      HtsCrsMatrixDeleter());

    if (Teuchos::nonnull(Timpl_)) {
      // Reuse the nonzero pattern.
      HTST::reprocess_numeric(Timpl_.get(), T_hts.get());
    } else {
      // Build from scratch.
      if (T_in.getCrsGraph().is_null()) {
        if (Teuchos::nonnull(out))
          *out << "HTS compute failed because T_in.getCrsGraph().is_null().\n";
        return;
      }
      if ( ! T_in.getCrsGraph()->isSorted()) {
        if (Teuchos::nonnull(out))
          *out << "HTS compute failed because ! T_in.getCrsGraph().isSorted().\n";
        return;
      }
      if ( ! T_in.isStorageOptimized()) {
        if (Teuchos::nonnull(out))
          *out << "HTS compute failed because ! T_in.isStorageOptimized().\n";
        return;
      }

      typename HTST::PreprocessArgs args;
      args.T = T_hts.get();
      args.max_nrhs = 1;
      args.nthreads = omp_get_max_threads();
      args.save_for_reprocess = true;
      typename HTST::Options opts;
      opts.levelset_block_size = levelset_block_size_;
      args.options = &opts;

      try {
        Timpl_ = Teuchos::rcpWithDealloc(HTST::preprocess(args), TImplDeleter());
      } catch (const std::exception& e) {
        if (Teuchos::nonnull(out))
          *out << "HTS preprocess threw: " << e.what() << "\n";
      }
    }
#endif
  }

  // HTS may not be able to handle a matrix, so query whether compute()
  // succeeded.
  bool isComputed () {
#ifdef HAVE_IFPACK2_SHYLUHTS
    return Teuchos::nonnull(Timpl_);
#else
    return false;
#endif
  }

  // Y := beta * Y + alpha * (M * X)
  void localApply (const MV& X, MV& Y,
                   const Teuchos::ETransp mode,
                   const scalar_type& alpha, const scalar_type& beta) const {
#ifdef HAVE_IFPACK2_SHYLUHTS
    const auto& X_view = X.template getLocalView<Kokkos::HostSpace>();
    const auto& Y_view = Y.template getLocalView<Kokkos::HostSpace>();
    // Only does something if #rhs > current capacity.
    HTST::reset_max_nrhs(Timpl_.get(), X_view.dimension_1());
    // Switch alpha and beta because of HTS's opposite convention.
    HTST::solve_omp(Timpl_.get(), X_view.data(), X_view.dimension_1(), Y_view.data(), beta, alpha);
#endif
  }

private:
#ifdef HAVE_IFPACK2_SHYLUHTS
  typedef ::Experimental::HTS<local_ordinal_type, size_t, scalar_type> HTST;
  typedef typename HTST::Impl TImpl;
  typedef typename HTST::CrsMatrix HtsCrsMatrix;

  struct TImplDeleter {
    void free (TImpl* impl) {
      HTST::delete_Impl(impl);
    }
  };

  struct HtsCrsMatrixDeleter {
    void free (HtsCrsMatrix* T) {
      HTST::delete_CrsMatrix(T);
    }
  };

  Teuchos::RCP<TImpl> Timpl_;
  bool transpose_, conjugate_;
  int levelset_block_size_;
#endif
};

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
LocalSparseTriangularSolver (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A)
{
  initializeState();
  typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  if (! A.is_null ()) {
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.is_null (), std::invalid_argument,
       "Ifpack2::LocalSparseTriangularSolver constructor: "
       "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
LocalSparseTriangularSolver (const Teuchos::RCP<const row_matrix_type>& A,
                             const Teuchos::RCP<Teuchos::FancyOStream>& out) :
  A_ (A),
  out_ (out)
{
  initializeState();
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG Ifpack2::LocalSparseTriangularSolver constructor"
          << std::endl;
  }
  typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  if (! A.is_null ()) {
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.is_null (), std::invalid_argument,
       "Ifpack2::LocalSparseTriangularSolver constructor: "
       "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
LocalSparseTriangularSolver ()
{
  initializeState();
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
LocalSparseTriangularSolver (const bool /* unused */, const Teuchos::RCP<Teuchos::FancyOStream>& out) :
  out_ (out)
{
  initializeState();
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG Ifpack2::LocalSparseTriangularSolver constructor"
          << std::endl;
  }
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::initializeState ()
{
  isInitialized_ = false;
  isComputed_ = false;
  numInitialize_ = 0;
  numCompute_ = 0;
  numApply_ = 0;
  initializeTime_ = 0.0;
  computeTime_ = 0.0;
  applyTime_ = 0.0;
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
~LocalSparseTriangularSolver ()
{}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
setParameters (const Teuchos::ParameterList& pl)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::Array;

  Details::TrisolverType::Enum trisolverType = Details::TrisolverType::Internal;
  do {
    static const char typeName[] = "trisolver: type";

    if ( ! pl.isType<std::string>(typeName)) break;

    // Map std::string <-> TrisolverType::Enum.
    Array<std::string> trisolverTypeStrs;
    Array<Details::TrisolverType::Enum> trisolverTypeEnums;
    Details::TrisolverType::loadPLTypeOption (trisolverTypeStrs, trisolverTypeEnums);
    Teuchos::StringToIntegralParameterEntryValidator<Details::TrisolverType::Enum>
      s2i(trisolverTypeStrs (), trisolverTypeEnums (), typeName, false);

    trisolverType = s2i.getIntegralValue(pl.get<std::string>(typeName));
  } while (0);

  if (trisolverType == Details::TrisolverType::HTS) {
    htsImpl_ = Teuchos::rcp (new HtsImpl ());
    htsImpl_->setParameters (pl);
  }
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
initialize ()
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::initialize: ";
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG " << prefix << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "You must call "
     "setMatrix() with a nonnull input matrix before you may call "
     "initialize() or compute().");
  if (A_crs_.is_null ()) {
    typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> crs_matrix_type;
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.is_null (), std::invalid_argument, prefix <<
       "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
  auto G = A_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_'s RowGraph G is null.  "
     "Please report this bug to the Ifpack2 developers.");
  // At this point, the graph MUST be fillComplete.  The "initialize"
  // (symbolic) part of setup only depends on the graph structure, so
  // the matrix itself need not be fill complete.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! G->isFillComplete (), std::runtime_error, "If you call this method, "
     "the matrix's graph must be fill complete.  It is not.");

  if (Teuchos::nonnull (htsImpl_))
    htsImpl_->initialize (*A_crs_);

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
compute ()
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::compute: ";
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG " << prefix << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "You must call "
     "setMatrix() with a nonnull input matrix before you may call "
     "initialize() or compute().");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_crs_.is_null (), std::logic_error, prefix << "A_ is nonnull, but "
     "A_crs_ is null.  Please report this bug to the Ifpack2 developers.");
  // At this point, the matrix MUST be fillComplete.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_crs_->isFillComplete (), std::runtime_error, "If you call this "
     "method, the matrix must be fill complete.  It is not.");

  if (! isInitialized_) {
    initialize ();
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (! isInitialized_, std::logic_error, prefix << "initialize() should have "
     "been called by this point, but isInitialized_ is false.  "
     "Please report this bug to the Ifpack2 developers.");

  if (Teuchos::nonnull (htsImpl_))
    htsImpl_->compute (*A_crs_, out_);

  auto A   = A_crs_->getLocalMatrix();
  auto ptr = A.graph.row_map;
  auto ind = A.graph.entries;
  auto val = A.values;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type annz    = A.nnz ();

  typename local_matrix_type::row_map_type::non_const_type  newptr ("ptr", ptr.dimension_0 ());
  typename local_matrix_type::index_type::non_const_type    newind ("ind", ind.dimension_0 ());
  typename local_matrix_type::values_type                   newval ("val", val.dimension_0 ());

  local_ordinal_type lclNumRows = A.numRows ();
  local_ordinal_type rowStart = 0;
  for (local_ordinal_type lclRow = 0; lclRow < lclNumRows; ++lclRow) {
     auto A_r = A.row (lclNumRows-lclRow-1);
     const local_ordinal_type numEnt = A_r.length;
     for (local_ordinal_type k = 0; k < numEnt; ++k) {
        newval(rowStart + k) = A_r.value(numEnt - k - 1);
        newind(rowStart + k) = A_r.colidx(numEnt - k - 1);
     }
     newptr(lclRow) = rowStart;
     rowStart = rowStart + numEnt;
  }

  Ulocal_ = Teuchos::rcp(new local_matrix_type("Ulocal", numRows, numCols, annz, newval, newptr, newind));

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type,
         global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type,
         global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::apply: ";
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG " << prefix;
    if (A_crs_.is_null ()) {
      *out_ << "A_crs_ is null!" << std::endl;
    }
    else {
      const std::string uplo = A_crs_->isUpperTriangular () ? "U" :
        (A_crs_->isLowerTriangular () ? "L" : "N");
      const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
        (mode == Teuchos::TRANS ? "T" : "N");
      const std::string diag =
        (A_crs_->getNodeNumDiags () < A_crs_->getNodeNumRows ()) ? "U" : "N";
      *out_ << "uplo=\"" << uplo
            << "\", trans=\"" << trans
            << "\", diag=\"" << diag << "\"" << std::endl;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (! isComputed (), std::runtime_error, prefix << "If compute() has not yet "
     "been called, or if you have changed the matrix via setMatrix(), you must "
     "call compute() before you may call this method.");
  // If isComputed() is true, it's impossible for the matrix to be
  // null, or for it not to be a Tpetra::CrsMatrix.
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::logic_error, prefix << "A_ is null.  "
     "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_crs_.is_null (), std::logic_error, prefix << "A_crs_ is null.  "
     "Please report this bug to the Ifpack2 developers.");
  // However, it _is_ possible that the user called resumeFill() on
  // the matrix, after calling compute().  This is NOT allowed.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_crs_->isFillComplete (), std::runtime_error, "If you call this "
     "method, the matrix must be fill complete.  It is not.  This means that "
     " you must have called resumeFill() on the matrix before calling apply(). "
     "This is NOT allowed.  Note that this class may use the matrix's data in "
     "place without copying it.  Thus, you cannot change the matrix and expect "
     "the solver to stay the same.  If you have changed the matrix, first call "
     "fillComplete() on it, then call compute() on this object, before you call"
     " apply().  You do NOT need to call setMatrix, as long as the matrix "
     "itself (that is, its address in memory) is the same.");

  auto G = A_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_'s RowGraph G is null.  "
     "Please report this bug to the Ifpack2 developers.");
  auto importer = G->getImporter ();
  auto exporter = G->getExporter ();

  if (! importer.is_null ()) {
    if (X_colMap_.is_null () || X_colMap_->getNumVectors () != X.getNumVectors ()) {
      X_colMap_ = rcp (new MV (importer->getTargetMap (), X.getNumVectors ()));
    }
    else {
      X_colMap_->putScalar (STS::zero ());
    }
    // See discussion of Github Issue #672 for why the Import needs to
    // use the ZERO CombineMode.  The case where the Export is
    // nontrivial is likely never exercised.
    X_colMap_->doImport (X, *importer, Tpetra::ZERO);
  }
  RCP<const MV> X_cur = importer.is_null () ? rcpFromRef (X) :
    Teuchos::rcp_const_cast<const MV> (X_colMap_);

  if (! exporter.is_null ()) {
    if (Y_rowMap_.is_null () || Y_rowMap_->getNumVectors () != Y.getNumVectors ()) {
      Y_rowMap_ = rcp (new MV (exporter->getSourceMap (), Y.getNumVectors ()));
    }
    else {
      Y_rowMap_->putScalar (STS::zero ());
    }
    Y_rowMap_->doExport (Y, *importer, Tpetra::ADD);
  }
  RCP<MV> Y_cur = exporter.is_null () ? rcpFromRef (Y) : Y_rowMap_;

  localApply (*X_cur, *Y_cur, mode, alpha, beta);

  if (! exporter.is_null ()) {
    Y.putScalar (STS::zero ());
    Y.doExport (*Y_cur, *exporter, Tpetra::ADD);
  }

  ++numApply_;
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
localApply (const MV& X,
            MV& Y,
            const Teuchos::ETransp mode,
            const scalar_type& alpha,
            const scalar_type& beta) const
{
  if (mode == Teuchos::NO_TRANS && Teuchos::nonnull (htsImpl_) &&
      htsImpl_->isComputed ()) {
    htsImpl_->localApply (X, Y, mode, alpha, beta);
    return;
  }

  using Teuchos::RCP;
  typedef scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;

  if (beta == STS::zero ()) {
    if (alpha == STS::zero ()) {
      Y.putScalar (STS::zero ()); // Y := 0 * Y (ignore contents of Y)
    }
    else { // alpha != 0
      A_crs_->template localSolve<ST, ST> (X, Y, mode);
      if (alpha != STS::one ()) {
        Y.scale (alpha);
      }
    }
  }
  else { // beta != 0
    if (alpha == STS::zero ()) {
      Y.scale (beta); // Y := beta * Y
    }
    else { // alpha != 0
      MV Y_tmp (Y, Teuchos::Copy);
      A_crs_->template localSolve<ST, ST> (X, Y_tmp, mode); // Y_tmp := M * X
      Y.update (alpha, Y_tmp, beta); // Y := beta * Y + alpha * Y_tmp
    }
  }
}


template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumInitialize () const {
  return numInitialize_;
}

template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumCompute () const {
  return numCompute_;
}

template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumApply () const {
  return numApply_;
}

template <class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getInitializeTime () const {
  return initializeTime_;
}

template<class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getComputeTime () const {
  return computeTime_;
}

template<class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getApplyTime() const {
  return applyTime_;
}

template <class MatrixType>
std::string
LocalSparseTriangularSolver<MatrixType>::
description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::LocalSparseTriangularSolver\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Matrix: not null"
       << ", Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", "
       << A_->getGlobalNumCols () << "]";
  }

  if (Teuchos::nonnull (htsImpl_))
    os << ", HTS computed: " << (htsImpl_->isComputed () ? "true" : "false");

  os << "}";
  return os.str ();
}

template <class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  // Default verbosity level is VERB_LOW, which prints only on Process
  // 0 of the matrix's communicator.
  const Teuchos::EVerbosityLevel vl
    = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl != Teuchos::VERB_NONE) {
    // Print only on Process 0 in the matrix's communicator.  If the
    // matrix is null, though, we have to get the communicator from
    // somewhere, so we ask Tpetra for its default communicator.  If
    // MPI is enabled, this wraps MPI_COMM_WORLD or a clone thereof.
    auto comm = A_.is_null () ?
      Tpetra::DefaultPlatform::getDefaultPlatform ().getComm () :
      A_->getComm ();

    // Users aren't supposed to do anything with the matrix on
    // processes where its communicator is null.
    if (! comm.is_null () && comm->getRank () == 0) {
      // By convention, describe() should always begin with a tab.
      Teuchos::OSTab tab0 (out);
      // Output is in YAML format.  We have to escape the class name,
      // because it has a colon.
      out << "\"Ifpack2::LocalSparseTriangularSolver\":" << endl;
      Teuchos::OSTab tab1 (out);
      out << "Scalar: " << Teuchos::TypeNameTraits<scalar_type>::name () << endl
          << "LocalOrdinal: " << Teuchos::TypeNameTraits<local_ordinal_type>::name () << endl
          << "GlobalOrdinal: " << Teuchos::TypeNameTraits<global_ordinal_type>::name () << endl
          << "Node: " << Teuchos::TypeNameTraits<node_type>::name () << endl;
    }
  }
}

template <class MatrixType>
Teuchos::RCP<const typename LocalSparseTriangularSolver<MatrixType>::map_type>
LocalSparseTriangularSolver<MatrixType>::
getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error,
     "Ifpack2::LocalSparseTriangularSolver::getDomainMap: "
     "The matrix is null.  Please call setMatrix() with a nonnull input "
     "before calling this method.");
  return A_->getDomainMap ();
}

template <class MatrixType>
Teuchos::RCP<const typename LocalSparseTriangularSolver<MatrixType>::map_type>
LocalSparseTriangularSolver<MatrixType>::
getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error,
     "Ifpack2::LocalSparseTriangularSolver::getRangeMap: "
     "The matrix is null.  Please call setMatrix() with a nonnull input "
     "before calling this method.");
  return A_->getRangeMap ();
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::setMatrix: ";

  // If the pointer didn't change, do nothing.  This is reasonable
  // because users are supposed to call this method with the same
  // object over all participating processes, and pointer identity
  // implies object identity.
  if (A.getRawPtr () != A_.getRawPtr ()) {
    // Check in serial or one-process mode if the matrix is square.
    TEUCHOS_TEST_FOR_EXCEPTION
      (! A.is_null () && A->getComm ()->getSize () == 1 &&
       A->getNodeNumRows () != A->getNodeNumCols (),
       std::runtime_error, prefix << "If A's communicator only contains one "
       "process, then A must be square.  Instead, you provided a matrix A with "
       << A->getNodeNumRows () << " rows and " << A->getNodeNumCols ()
       << " columns.");

    // It's legal for A to be null; in that case, you may not call
    // initialize() until calling setMatrix() with a nonnull input.
    // Regardless, setting the matrix invalidates the preconditioner.
    isInitialized_ = false;
    isComputed_ = false;

    typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> crs_matrix_type;
    if (A.is_null ()) {
      A_crs_ = Teuchos::null;
      A_ = Teuchos::null;
    }
    else { // A is not null
      Teuchos::RCP<const crs_matrix_type> A_crs =
        Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A);
      TEUCHOS_TEST_FOR_EXCEPTION
        (A_crs.is_null (), std::invalid_argument, prefix <<
         "The input matrix A is not a Tpetra::CrsMatrix.");
      A_crs_ = A_crs;
      A_ = A;
    }

    if (Teuchos::nonnull (htsImpl_))
      htsImpl_->reset ();
  } // pointers are not the same
}

} // namespace Ifpack2

#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::LocalSparseTriangularSolver< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
