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
#include "Tpetra_Core.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Tpetra_Details_determineLocalTriangularStructure.hpp"

#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
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
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    Timpl_ = Teuchos::null;
    levelset_block_size_ = 1;
#endif
  }

  void setParameters (const Teuchos::ParameterList& pl) {
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
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
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    reset();
    transpose_ = conjugate_ = false;
#endif
  }

  void compute (const crs_matrix_type& T_in, const Teuchos::RCP<Teuchos::FancyOStream>& out) {
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
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
#ifdef _OPENMP
      args.nthreads = omp_get_max_threads();
#else
      args.nthreads = 1;
#endif
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
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    return Teuchos::nonnull(Timpl_);
#else
    return false;
#endif
  }

  // Y := beta * Y + alpha * (M * X)
  void localApply (const MV& X, MV& Y,
                   const Teuchos::ETransp mode,
                   const scalar_type& alpha, const scalar_type& beta) const {
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    const auto& X_view = X.template getLocalView<Kokkos::HostSpace>();
    const auto& Y_view = Y.template getLocalView<Kokkos::HostSpace>();
    // Only does something if #rhs > current capacity.
    HTST::reset_max_nrhs(Timpl_.get(), X_view.extent(1));
    // Switch alpha and beta because of HTS's opposite convention.
    HTST::solve_omp(Timpl_.get(),
                    // For std/Kokkos::complex.
                    reinterpret_cast<const scalar_type*>(X_view.data()),
                    X_view.extent(1),
                    // For std/Kokkos::complex.
                    reinterpret_cast<scalar_type*>(Y_view.data()),
                    beta, alpha);
#endif
  }

private:
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
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
  reverseStorage_ = false;
  isInternallyChanged_ = false;
  numInitialize_ = 0;
  numCompute_ = 0;
  numApply_ = 0;
  initializeTime_ = 0.0;
  computeTime_ = 0.0;
  applyTime_ = 0.0;
  uplo_ = "N";
  diag_ = "N";
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

  if (pl.isParameter("trisolver: reverse U"))
    reverseStorage_ = pl.get<bool>("trisolver: reverse U");

  TEUCHOS_TEST_FOR_EXCEPTION
    (reverseStorage_ && trisolverType == Details::TrisolverType::HTS,
     std::logic_error, "Ifpack2::LocalSparseTriangularSolver::setParameters: "
     "You are not allowed to enable both HTS and the \"trisolver: reverse U\" "
     "options.  See GitHub issue #2647.");
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
initialize ()
{
  using crs_matrix_type = Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type>;
  using local_matrix_type = typename crs_matrix_type::local_matrix_type;

  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::initialize: ";
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG " << prefix << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "You must call "
     "setMatrix() with a nonnull input matrix before you may call "
     "initialize() or compute().");
  if (A_crs_.is_null ()) {
    auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.get () == nullptr, std::invalid_argument,
       prefix << "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
  auto G = A_crs_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_crs_'s RowGraph G is null.  "
     "Please report this bug to the Ifpack2 developers.");
  // At this point, the graph MUST be fillComplete.  The "initialize"
  // (symbolic) part of setup only depends on the graph structure, so
  // the matrix itself need not be fill complete.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! G->isFillComplete (), std::runtime_error, "If you call this method, "
     "the matrix's graph must be fill complete.  It is not.");

  // FIXME (mfh 01,02 Jun 2018) isUpperTriangular has been DEPRECATED.
  // See GitHub Issue #2630.  I'm using isUpperTriangularImpl ONLY to
  // avoid deprecated warnings.  Users may NOT call this method.
  //
  // FIXME (mfh 02 Jun 2018) Move the
  // determineLocalTriangularStructure call above this test, so we can
  // use that result, rather than the deprecated method.
  if (reverseStorage_ && A_crs_->isUpperTriangularImpl() && htsImpl_.is_null()) {
    // Reverse the storage for an upper triangular matrix
    auto Alocal = A_crs_->getLocalMatrix();
    auto ptr    = Alocal.graph.row_map;
    auto ind    = Alocal.graph.entries;
    auto val    = Alocal.values;

    auto numRows = Alocal.numRows();
    auto numCols = Alocal.numCols();
    auto numNnz = Alocal.nnz();

    typename decltype(ptr)::non_const_type  newptr ("ptr", ptr.extent (0));
    typename decltype(ind)::non_const_type  newind ("ind", ind.extent (0));
    decltype(val)                           newval ("val", val.extent (0));

    // FIXME: The code below assumes UVM
    crs_matrix_type::execution_space::fence();
    newptr(0) = 0;
    for (local_ordinal_type row = 0, rowStart = 0; row < numRows; ++row) {
      auto A_r = Alocal.row(numRows-1 - row);

      auto numEnt = A_r.length;
      for (local_ordinal_type k = 0; k < numEnt; ++k) {
        newind(rowStart + k) = numCols-1 - A_r.colidx(numEnt-1 - k);
        newval(rowStart + k) = A_r.value (numEnt-1 - k);
      }
      rowStart += numEnt;
      newptr(row+1) = rowStart;
    }
    crs_matrix_type::execution_space::fence();

    // Reverse maps
    using map_type = typename crs_matrix_type::map_type;
    Teuchos::RCP<map_type> newRowMap, newColMap;
    {
      // Reverse row map
      auto rowMap = A_->getRowMap();
      auto numElems = rowMap->getNodeNumElements();
      auto rowElems = rowMap->getNodeElementList();

      Teuchos::Array<global_ordinal_type> newRowElems(rowElems.size());
      for (size_t i = 0; i < numElems; i++)
        newRowElems[i] = rowElems[numElems-1 - i];

      newRowMap = Teuchos::rcp(new map_type(rowMap->getGlobalNumElements(), newRowElems, rowMap->getIndexBase(), rowMap->getComm()));
    }
    {
      // Reverse column map
      auto colMap = A_->getColMap();
      auto numElems = colMap->getNodeNumElements();
      auto colElems = colMap->getNodeElementList();

      Teuchos::Array<global_ordinal_type> newColElems(colElems.size());
      for (size_t i = 0; i < numElems; i++)
        newColElems[i] = colElems[numElems-1 - i];

      newColMap = Teuchos::rcp(new map_type(colMap->getGlobalNumElements(), newColElems, colMap->getIndexBase(), colMap->getComm()));
    }

    // Construct new matrix
    local_matrix_type newLocalMatrix("Upermuted", numRows, numCols, numNnz, newval, newptr, newind);

    A_crs_ = Teuchos::rcp(new crs_matrix_type(newLocalMatrix, newRowMap, newColMap, A_crs_->getDomainMap(), A_crs_->getRangeMap()));

    isInternallyChanged_ = true;
  }

  if (Teuchos::nonnull (htsImpl_))
  {
    htsImpl_->initialize (*A_crs_);
    isInternallyChanged_ = true;
  }

  auto lclMatrix = A_crs_->getLocalMatrix ();
  auto lclRowMap = A_crs_->getRowMap ()->getLocalMap ();
  auto lclColMap = A_crs_->getColMap ()->getLocalMap ();
  using Tpetra::Details::determineLocalTriangularStructure;
  // mfh 30 Apr 2018: See GitHub Issue #2658 for why this is false.
  constexpr bool ignoreMapsForTriangularStructure = true;
  auto result =
    determineLocalTriangularStructure (lclMatrix.graph, lclRowMap, lclColMap,
                                       ignoreMapsForTriangularStructure);
  using LO = local_ordinal_type;
  const LO lclNumRows = lclRowMap.getNodeNumElements ();
  const LO lclNumCols = lclColMap.getNodeNumElements ();
  // NOTE (mfh 30 Apr 2018) Original test for implicit unit diagonal was
  //
  // (A_crs_->getNodeNumDiags () < A_crs_->getNodeNumRows ()) ? "U" : "N";
  //
  // I don't agree with this test -- it's not an implicitly stored
  // unit diagonal if there are SOME entries -- but I'm leaving it for
  // backwards compatibility.
  this->diag_ = (result.diagCount < lclNumRows) ? "U" : "N";
  this->uplo_ = result.couldBeLowerTriangular ? "L" :
    (result.couldBeUpperTriangular ? "U" : "N");

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
      Teuchos::RCP<const crs_matrix_type> A_crs =
          Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_);
      const std::string uplo = this->uplo_;
      const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
        (mode == Teuchos::TRANS ? "T" : "N");
      const std::string diag = this->diag_;
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

  auto G = A_crs_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_crs_'s RowGraph G is null.  "
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
localTriangularSolve (const MV& Y,
                      MV& X,
                      const Teuchos::ETransp mode) const
{
  using Teuchos::CONJ_TRANS;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  typedef Kokkos::HostSpace host_memory_space;
  using device_type = typename MV::device_type;
  using dev_memory_space = typename device_type::memory_space;
  const char tfecfFuncName[] = "localTriangularSolve: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (! A_crs_->isFillComplete (), std::runtime_error,
     "The matrix is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (! X.isConstantStride () || ! Y.isConstantStride (), std::invalid_argument,
     "X and Y must be constant stride.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    ( A_crs_->getNodeNumRows() > 0 && this->uplo_ == "N", std::runtime_error,
      "The matrix is neither upper triangular or lower triangular.  "
      "You may only call this method if the matrix is triangular.  "
      "Remember that this is a local (per MPI process) property, and that "
      "Tpetra only knows how to do a local (per process) triangular solve.");
  using STS = Teuchos::ScalarTraits<scalar_type>;
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (STS::isComplex && mode == TRANS, std::logic_error, "This method does "
     "not currently support non-conjugated transposed solve (mode == "
     "Teuchos::TRANS) for complex scalar types.");

  // FIXME (mfh 19 May 2016) This makes some Ifpack2 tests fail.
  //
  // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
  //   (Y.template need_sync<device_type> () && !
  //    Y.template need_sync<Kokkos::HostSpace> (), std::runtime_error,
  //    "Y must be sync'd to device memory before you may call this method.");

  const std::string uplo = this->uplo_;
  const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
    (mode == Teuchos::TRANS ? "T" : "N");
  const std::string diag = this->diag_;
  auto A_lcl = this->A_crs_->getLocalMatrix ();

  // NOTE (mfh 20 Aug 2017): KokkosSparse::trsv currently is a
  // sequential, host-only code.  See
  // https://github.com/kokkos/kokkos-kernels/issues/48.  This
  // means that we need to sync to host, then sync back to device
  // when done.
  X.template sync<host_memory_space> ();
  const_cast<MV&> (Y).template sync<host_memory_space> ();
  X.template modify<host_memory_space> (); // we will write to X

  if (X.isConstantStride () && Y.isConstantStride ()) {
    auto X_lcl = X.template getLocalView<host_memory_space> ();
    auto Y_lcl = Y.template getLocalView<host_memory_space> ();
    KokkosSparse::trsv (uplo.c_str (), trans.c_str (), diag.c_str (),
                        A_lcl, Y_lcl, X_lcl);
  }
  else {
    const size_t numVecs =
      std::min (X.getNumVectors (), Y.getNumVectors ());
    for (size_t j = 0; j < numVecs; ++j) {
      auto X_j = X.getVector (j);
      auto Y_j = X.getVector (j);
      auto X_lcl = X_j->template getLocalView<host_memory_space> ();
      auto Y_lcl = Y_j->template getLocalView<host_memory_space> ();
      KokkosSparse::trsv (uplo.c_str (), trans.c_str (),
                          diag.c_str (), A_lcl, Y_lcl, X_lcl);
    }
  }

  X.template sync<dev_memory_space> ();
  const_cast<MV&> (Y).template sync<dev_memory_space> ();
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
      this->localTriangularSolve (X, Y, mode);
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
      this->localTriangularSolve (X, Y_tmp, mode); // Y_tmp := M * X
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
      Tpetra::getDefaultComm () :
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
  if (A.getRawPtr () != A_.getRawPtr () || isInternallyChanged_) {
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
