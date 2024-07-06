// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP

#include "Ifpack2_LocalSparseTriangularSolver_decl.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Tpetra_Details_determineLocalTriangularStructure.hpp"
#include "KokkosSparse_trsv.hpp"

#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
# include "shylu_hts.hpp"
#endif

namespace Ifpack2 {

namespace Details {
struct TrisolverType {
  enum Enum {
    Internal, //!< Tpetra::CrsMatrix::localSolve
    HTS,      //!< Multicore ShyLU/HTS
    KSPTRSV   //!< Multicore/GPU KokkosKernels sptrsv
  };

  static void loadPLTypeOption (Teuchos::Array<std::string>& type_strs, Teuchos::Array<Enum>& type_enums) {
    type_strs.resize(3);
    type_strs[0] = "Internal";
    type_strs[1] = "HTS";
    type_strs[2] = "KSPTRSV";
    type_enums.resize(3);
    type_enums[0] = Internal;
    type_enums[1] = HTS;
    type_enums[2] = KSPTRSV;
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
    (void)pl;
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
    (void)T_in;
    (void)out;
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    using Teuchos::ArrayRCP;

    auto rowptr = T_in.getLocalRowPtrsHost();
    auto colidx = T_in.getLocalIndicesHost();
    auto val = T_in.getLocalValuesHost(Tpetra::Access::ReadOnly);
    Kokkos::fence();

    Teuchos::RCP<HtsCrsMatrix> T_hts = Teuchos::rcpWithDealloc(
      HTST::make_CrsMatrix(rowptr.size() - 1,
                           rowptr.data(), colidx.data(), 
                           // For std/Kokkos::complex.
                           reinterpret_cast<const scalar_type*>(val.data()),
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
                   const Teuchos::ETransp /* mode */,
                   const scalar_type& alpha, const scalar_type& beta) const {
    (void)X;
    (void)Y;
    (void)alpha;
    (void)beta;
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
    const auto& X_view = X.getLocalViewHost (Tpetra::Access::ReadOnly);
    const auto& Y_view = Y.getLocalViewHost (Tpetra::Access::ReadWrite);

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
  isKokkosKernelsSptrsv_ = false;
  isKokkosKernelsStream_ = false;
  num_streams_ = 0;
  uplo_ = "N";
  diag_ = "N";
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
~LocalSparseTriangularSolver ()
{
  if (!isKokkosKernelsStream_) {
    if (Teuchos::nonnull (kh_)) {
      kh_->destroy_sptrsv_handle();
    }
  }
  else {
    for (int i = 0; i < num_streams_; i++) {
      if (Teuchos::nonnull (kh_v_[i])) {
        kh_v_[i]->destroy_sptrsv_handle();
      }
    }
  }
}

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

  if (trisolverType == Details::TrisolverType::KSPTRSV) {
    this->isKokkosKernelsSptrsv_ = true;
  }
  else {
    this->isKokkosKernelsSptrsv_ = false;
  }

  if (pl.isParameter("trisolver: reverse U"))
    reverseStorage_ = pl.get<bool>("trisolver: reverse U");

  TEUCHOS_TEST_FOR_EXCEPTION
    (reverseStorage_ && (trisolverType == Details::TrisolverType::HTS || trisolverType == Details::TrisolverType::KSPTRSV),
     std::logic_error, "Ifpack2::LocalSparseTriangularSolver::setParameters: "
     "You are not allowed to enable both HTS or KSPTRSV and the \"trisolver: reverse U\" "
     "options.  See GitHub issue #2647.");
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
initialize ()
{
  using Tpetra::Details::determineLocalTriangularStructure;

  using local_matrix_type = typename crs_matrix_type::local_matrix_device_type;
  using LO = local_ordinal_type;

  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::initialize: ";
  if (! out_.is_null ()) {
    *out_ << ">>> DEBUG " << prefix << std::endl;
  }

  if (!isKokkosKernelsStream_) {
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

    // mfh 30 Apr 2018: See GitHub Issue #2658.
    constexpr bool ignoreMapsForTriStructure = true;
    auto lclTriStructure = [&] {
      auto lclMatrix = A_crs_->getLocalMatrixDevice ();
      auto lclRowMap = A_crs_->getRowMap ()->getLocalMap ();
      auto lclColMap = A_crs_->getColMap ()->getLocalMap ();
      auto lclTriStruct =
        determineLocalTriangularStructure (lclMatrix.graph,
                                           lclRowMap,
                                           lclColMap,
                                           ignoreMapsForTriStructure);
      const LO lclNumRows = lclRowMap.getLocalNumElements ();
      this->diag_ = (lclTriStruct.diagCount < lclNumRows) ? "U" : "N";
      this->uplo_ = lclTriStruct.couldBeLowerTriangular ? "L" :
        (lclTriStruct.couldBeUpperTriangular ? "U" : "N");
      return lclTriStruct;
    } ();
    
    if (reverseStorage_ && lclTriStructure.couldBeUpperTriangular &&
        htsImpl_.is_null ()) {
      // Reverse the storage for an upper triangular matrix
      auto Alocal = A_crs_->getLocalMatrixDevice();
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
      typename crs_matrix_type::execution_space().fence();
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
      typename crs_matrix_type::execution_space().fence();
    
      // Reverse maps
      Teuchos::RCP<map_type> newRowMap, newColMap;
      {
        // Reverse row map
        auto rowMap = A_->getRowMap();
        auto numElems = rowMap->getLocalNumElements();
        auto rowElems = rowMap->getLocalElementList();
    
        Teuchos::Array<global_ordinal_type> newRowElems(rowElems.size());
        for (size_t i = 0; i < numElems; i++)
          newRowElems[i] = rowElems[numElems-1 - i];
    
        newRowMap = Teuchos::rcp(new map_type(rowMap->getGlobalNumElements(), newRowElems, rowMap->getIndexBase(), rowMap->getComm()));
      }
      {
        // Reverse column map
        auto colMap = A_->getColMap();
        auto numElems = colMap->getLocalNumElements();
        auto colElems = colMap->getLocalElementList();
    
        Teuchos::Array<global_ordinal_type> newColElems(colElems.size());
        for (size_t i = 0; i < numElems; i++)
          newColElems[i] = colElems[numElems-1 - i];
    
        newColMap = Teuchos::rcp(new map_type(colMap->getGlobalNumElements(), newColElems, colMap->getIndexBase(), colMap->getComm()));
      }
    
      // Construct new matrix
      local_matrix_type newLocalMatrix("Upermuted", numRows, numCols, numNnz, newval, newptr, newind);
    
      A_crs_ = Teuchos::rcp(new crs_matrix_type(newLocalMatrix, newRowMap, newColMap, A_crs_->getDomainMap(), A_crs_->getRangeMap()));
    
      isInternallyChanged_ = true;
    
      // FIXME (mfh 18 Apr 2019) Recomputing this is unnecessary, but I
      // didn't want to break any invariants, especially considering
      // that this branch is likely poorly tested.
      auto newLclTriStructure =
        determineLocalTriangularStructure (newLocalMatrix.graph,
                                           newRowMap->getLocalMap (),
                                           newColMap->getLocalMap (),
                                           ignoreMapsForTriStructure);
      const LO newLclNumRows = newRowMap->getLocalNumElements ();
      this->diag_ = (newLclTriStructure.diagCount < newLclNumRows) ? "U" : "N";
      this->uplo_ = newLclTriStructure.couldBeLowerTriangular ? "L" :
        (newLclTriStructure.couldBeUpperTriangular ? "U" : "N");
    }
  }
  else {
    for (int i = 0; i < num_streams_; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (A_crs_v_[i].is_null (), std::runtime_error, prefix << "You must call "
         "setMatrix() with a nonnull input matrix before you may call "
         "initialize() or compute().");
      auto G = A_crs_v_[i]->getGraph ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (G.is_null (), std::logic_error, prefix << "A_crs_ are nonnull, "
         "but A_crs_'s RowGraph G is null.  "
         "Please report this bug to the Ifpack2 developers.");
      // At this point, the graph MUST be fillComplete.  The "initialize"
      // (symbolic) part of setup only depends on the graph structure, so
      // the matrix itself need not be fill complete.
      TEUCHOS_TEST_FOR_EXCEPTION
        (! G->isFillComplete (), std::runtime_error, "If you call this method, "
         "the matrix's graph must be fill complete.  It is not.");
	  
      // mfh 30 Apr 2018: See GitHub Issue #2658.
      constexpr bool ignoreMapsForTriStructure = true;
      std::string prev_uplo_ = this->uplo_;
      std::string prev_diag_ = this->diag_;
      auto lclMatrix = A_crs_v_[i]->getLocalMatrixDevice ();
      auto lclRowMap = A_crs_v_[i]->getRowMap ()->getLocalMap ();
      auto lclColMap = A_crs_v_[i]->getColMap ()->getLocalMap ();
      auto lclTriStruct =
        determineLocalTriangularStructure (lclMatrix.graph,
                                           lclRowMap,
                                           lclColMap,
                                           ignoreMapsForTriStructure);
      const LO lclNumRows = lclRowMap.getLocalNumElements ();
      this->diag_ = (lclTriStruct.diagCount < lclNumRows) ? "U" : "N";
      this->uplo_ = lclTriStruct.couldBeLowerTriangular ? "L" :
        (lclTriStruct.couldBeUpperTriangular ? "U" : "N");
      if (i > 0) {
        TEUCHOS_TEST_FOR_EXCEPTION
          ((this->diag_ != prev_diag_) || (this->uplo_ != prev_uplo_),
           std::logic_error, prefix << "A_crs_'s structures in streams "
           "are different. Please report this bug to the Ifpack2 developers.");
      }
    }
  }

  if (Teuchos::nonnull (htsImpl_))
  {
    htsImpl_->initialize (*A_crs_);
    isInternallyChanged_ = true;
  }

  const bool ksptrsv_valid_uplo = (this->uplo_ != "N");
  kh_v_nonnull_ = false;
  if (this->isKokkosKernelsSptrsv_ && ksptrsv_valid_uplo && this->diag_ != "U")
  {
    if (!isKokkosKernelsStream_) {
      kh_ = Teuchos::rcp (new k_handle());
      const bool is_lower_tri = (this->uplo_ == "L") ? true : false;

      auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_, true);
      auto Alocal = A_crs->getLocalMatrixDevice();
      auto ptr    = Alocal.graph.row_map;
      auto ind    = Alocal.graph.entries;
      auto val    = Alocal.values;

      auto numRows = Alocal.numRows();    
      kh_->create_sptrsv_handle(kokkosKernelsAlgorithm(), numRows, is_lower_tri);
      KokkosSparse::Experimental::sptrsv_symbolic(kh_.getRawPtr(), ptr, ind, val);
    } else {
      kh_v_ = std::vector< Teuchos::RCP<k_handle> >(num_streams_);
      for (int i = 0; i < num_streams_; i++) {
        kh_v_[i] = Teuchos::rcp (new k_handle ());
        auto A_crs_i = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_crs_v_[i], true);
        auto Alocal_i = A_crs_i->getLocalMatrixDevice();
        auto ptr_i    = Alocal_i.graph.row_map;
        auto ind_i    = Alocal_i.graph.entries;
        auto val_i    = Alocal_i.values;
    
        auto numRows_i = Alocal_i.numRows();

        const bool is_lower_tri = (this->uplo_ == "L") ? true : false;
        kh_v_[i]->create_sptrsv_handle(kokkosKernelsAlgorithm(), numRows_i, is_lower_tri);
        KokkosSparse::Experimental::sptrsv_symbolic(kh_v_[i].getRawPtr(), ptr_i, ind_i, val_i);
      }
      kh_v_nonnull_ = true;
    }
  }

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
KokkosSparse::Experimental::SPTRSVAlgorithm 
LocalSparseTriangularSolver<MatrixType>::kokkosKernelsAlgorithm() const
{
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(KOKKOS_ENABLE_CUDA)
  // CuSparse only supports int type ordinals 
  // and scalar types of float, double, float complex and double complex
  if constexpr (std::is_same<Kokkos::Cuda, HandleExecSpace>::value &&
      std::is_same<int, local_ordinal_type>::value &&
    (std::is_same<scalar_type, float>::value ||
      std::is_same<scalar_type, double>::value ||
      std::is_same<scalar_type, Kokkos::complex<float>>::value ||
      std::is_same<scalar_type, Kokkos::complex<double>>::value))
  {
    return KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE;
  }
#endif
  return KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1;
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

  if (!isKokkosKernelsStream_) {
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
  }
  else {
    for(int i = 0; i < num_streams_; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (A_crs_v_[i].is_null (), std::runtime_error, prefix << "You must call "
         "setMatrices() with a nonnull input matrix before you may call "
         "initialize() or compute().");
      // At this point, the matrix MUST be fillComplete.
      TEUCHOS_TEST_FOR_EXCEPTION
        (! A_crs_v_[i]->isFillComplete (), std::runtime_error, "If you call this "
         "method, the matrix must be fill complete.  It is not.");
    }
  }

  if (! isInitialized_) {
    initialize ();
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (! isInitialized_, std::logic_error, prefix << "initialize() should have "
     "been called by this point, but isInitialized_ is false.  "
     "Please report this bug to the Ifpack2 developers.");

// NOTE (Nov-09-2022):
// For Cuda >= 11.3 (using cusparseSpSV), always call symbolic during compute
// even when matrix values are changed with the same sparsity pattern.
// For Cuda >= 12.1 has a new cusparseSpSV_updateMatrix function just for updating the
// values that is substantially faster.
// This would all be much better handled via a KokkosSparse::Experimental::sptrsv_numeric(...)
// that could hide the Cuda implementation details.
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && (CUDA_VERSION >= 11030)
  if constexpr ( std::is_same_v<typename crs_matrix_type::execution_space, Kokkos::Cuda> )
  {
    if (this->isKokkosKernelsSptrsv_) {
      if (Teuchos::nonnull(kh_) && !isKokkosKernelsStream_) {
        auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_crs_, true);
        auto Alocal = A_crs->getLocalMatrixDevice();
        auto val    = Alocal.values;
  #if (CUSPARSE_VERSION >= 12100)
        auto *sptrsv_handle = kh_->get_sptrsv_handle();
        auto cusparse_handle = sptrsv_handle->get_cuSparseHandle();
        cusparseSpSV_updateMatrix(cusparse_handle->handle,
                          cusparse_handle->spsvDescr,
                          val.data(),
                          CUSPARSE_SPSV_UPDATE_GENERAL);
  #else
        auto ptr    = Alocal.graph.row_map;
        auto ind    = Alocal.graph.entries;
        KokkosSparse::Experimental::sptrsv_symbolic(kh_.getRawPtr(), ptr, ind, val);
  #endif
      } else if (kh_v_nonnull_) {
        for (int i = 0; i < num_streams_; i++) {
          auto A_crs_i = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_crs_v_[i], true);
          auto Alocal_i = A_crs_i->getLocalMatrixDevice();
          auto val_i    = Alocal_i.values;
  #if (CUSPARSE_VERSION >= 12100)
          auto *sptrsv_handle = kh_v_[i]->get_sptrsv_handle();
          auto cusparse_handle = sptrsv_handle->get_cuSparseHandle();
          KOKKOS_CUSPARSE_SAFE_CALL(
              cusparseSetStream(cusparse_handle->handle, exec_space_instances_[i].cuda_stream()));
          cusparseSpSV_updateMatrix(cusparse_handle->handle,
                            cusparse_handle->spsvDescr,
                            val_i.data(),
                            CUSPARSE_SPSV_UPDATE_GENERAL);
  #else
          auto ptr_i    = Alocal_i.graph.row_map;
          auto ind_i    = Alocal_i.graph.entries;
          KokkosSparse::Experimental::sptrsv_symbolic(exec_space_instances_[i], kh_v_[i].getRawPtr(), ptr_i, ind_i, val_i);
  #endif
        }
      }
    }
  }
#endif

  if (! isComputed_) {//Only compute if not computed before
    if (Teuchos::nonnull (htsImpl_))
      htsImpl_->compute (*A_crs_, out_);

    isComputed_ = true;
    ++numCompute_;
  }
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
    if (!isKokkosKernelsStream_) {
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
    else {
      for (int i = 0; i < num_streams_; i++) {
        if (A_crs_v_[i].is_null ()) {
          *out_ << "A_crs_v_[" << i << "]" << " is null!" << std::endl;
        }
        else {
          const std::string uplo = this->uplo_;
          const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
            (mode == Teuchos::TRANS ? "T" : "N");
          const std::string diag = this->diag_;
          *out_ << "A_crs_v_[" << i << "]: "
                << "uplo=\"" << uplo
                << "\", trans=\"" << trans
                << "\", diag=\"" << diag << "\"" << std::endl;
        }
      }
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (! isComputed (), std::runtime_error, prefix << "If compute() has not yet "
     "been called, or if you have changed the matrix via setMatrix(), you must "
     "call compute() before you may call this method.");
  // If isComputed() is true, it's impossible for the matrix to be
  // null, or for it not to be a Tpetra::CrsMatrix.
  if (!isKokkosKernelsStream_) {
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
  }
  else {
    for (int i = 0; i < num_streams_; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (A_crs_v_[i].is_null (), std::logic_error, prefix << "A_crs_ is null.  "
         "Please report this bug to the Ifpack2 developers.");
      // However, it _is_ possible that the user called resumeFill() on
      // the matrix, after calling compute().  This is NOT allowed.
      TEUCHOS_TEST_FOR_EXCEPTION
        (! A_crs_v_[i]->isFillComplete (), std::runtime_error, "If you call this "
         "method, the matrix must be fill complete.  It is not.  This means that "
         " you must have called resumeFill() on the matrix before calling apply(). "
         "This is NOT allowed.  Note that this class may use the matrix's data in "
         "place without copying it.  Thus, you cannot change the matrix and expect "
         "the solver to stay the same.  If you have changed the matrix, first call "
         "fillComplete() on it, then call compute() on this object, before you call"
         " apply().  You do NOT need to call setMatrix, as long as the matrix "
         "itself (that is, its address in memory) is the same.");
    }
  }

  RCP<const MV> X_cur;
  RCP<MV> Y_cur;

  if (!isKokkosKernelsStream_) {
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
    X_cur = importer.is_null () ? rcpFromRef (X) :
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
    Y_cur = exporter.is_null () ? rcpFromRef (Y) : Y_rowMap_;
  }
  else {
    // Currently assume X and Y are local vectors (same sizes as A_crs_).
    // Should do a better job here!!!
    X_cur = rcpFromRef (X);
    Y_cur = rcpFromRef (Y);
  }

  localApply (*X_cur, *Y_cur, mode, alpha, beta);

  if (!isKokkosKernelsStream_) {
    auto G = A_crs_->getGraph ();
    auto exporter = G->getExporter ();
    if (! exporter.is_null ()) {
      Y.putScalar (STS::zero ());
      Y.doExport (*Y_cur, *exporter, Tpetra::ADD);
    }
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
  const char tfecfFuncName[] = "localTriangularSolve: ";

  if (!isKokkosKernelsStream_)
  {
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! A_crs_->isFillComplete (), std::runtime_error,
       "The matrix is not fill complete.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      ( A_crs_->getLocalNumRows() > 0 && this->uplo_ == "N", std::runtime_error,
        "The matrix is neither upper triangular or lower triangular.  "
        "You may only call this method if the matrix is triangular.  "
        "Remember that this is a local (per MPI process) property, and that "
        "Tpetra only knows how to do a local (per process) triangular solve.");
  }
  else
  {
    for (int i = 0; i < num_streams_; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! A_crs_v_[i]->isFillComplete (), std::runtime_error,
         "The matrix is not fill complete.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        ( A_crs_v_[i]->getLocalNumRows() > 0 && this->uplo_ == "N", std::runtime_error,
          "The matrix is neither upper triangular or lower triangular.  "
          "You may only call this method if the matrix is triangular.  "
          "Remember that this is a local (per MPI process) property, and that "
          "Tpetra only knows how to do a local (per process) triangular solve.");
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (! X.isConstantStride () || ! Y.isConstantStride (), std::invalid_argument,
     "X and Y must be constant stride.");
  using STS = Teuchos::ScalarTraits<scalar_type>;
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (STS::isComplex && mode == TRANS, std::logic_error, "This method does "
     "not currently support non-conjugated transposed solve (mode == "
     "Teuchos::TRANS) for complex scalar types.");

  const std::string uplo = this->uplo_;
  const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
    (mode == Teuchos::TRANS ? "T" : "N");
  const size_t numVecs = std::min (X.getNumVectors (), Y.getNumVectors ());

  if (Teuchos::nonnull(kh_) && this->isKokkosKernelsSptrsv_ && trans == "N")
  {
    auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (this->A_);
    auto A_lclk = A_crs->getLocalMatrixDevice ();
    auto ptr    = A_lclk.graph.row_map;
    auto ind    = A_lclk.graph.entries;
    auto val    = A_lclk.values;
	
    for (size_t j = 0; j < numVecs; ++j) {
      auto X_j = X.getVectorNonConst (j);
      auto Y_j = Y.getVector (j);
      auto X_lcl = X_j->getLocalViewDevice (Tpetra::Access::ReadWrite);
      auto Y_lcl = Y_j->getLocalViewDevice (Tpetra::Access::ReadOnly);
      auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::ALL (), 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::ALL (), 0);
      KokkosSparse::Experimental::sptrsv_solve(kh_.getRawPtr(), ptr, ind, val, Y_lcl_1d, X_lcl_1d);
      // TODO is this fence needed...
      typename k_handle::HandleExecSpace().fence();
    }
  } // End using regular interface of Kokkos Kernels Sptrsv
  else if (kh_v_nonnull_ && this->isKokkosKernelsSptrsv_ && trans == "N")
  {
    std::vector<lno_row_view_t>        ptr_v(num_streams_);
    std::vector<lno_nonzero_view_t>    ind_v(num_streams_);
    std::vector<scalar_nonzero_view_t> val_v(num_streams_);
    std::vector<k_handle *>            KernelHandle_rawptr_v_(num_streams_);
    for (size_t j = 0; j < numVecs; ++j) {
      auto X_j = X.getVectorNonConst (j);
      auto Y_j = Y.getVector (j);
      auto X_lcl = X_j->getLocalViewDevice (Tpetra::Access::ReadWrite);
      auto Y_lcl = Y_j->getLocalViewDevice (Tpetra::Access::ReadOnly);
      auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::ALL (), 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::ALL (), 0);
      std::vector<decltype(X_lcl_1d)> x_v(num_streams_);
      std::vector<decltype(Y_lcl_1d)> y_v(num_streams_);
      local_ordinal_type stream_begin = 0;
      local_ordinal_type stream_end;
      for (int i = 0; i < num_streams_; i++) {
        auto A_crs_i = Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_crs_v_[i]);
        auto Alocal_i = A_crs_i->getLocalMatrixDevice();
        ptr_v[i] = Alocal_i.graph.row_map;
        ind_v[i] = Alocal_i.graph.entries;
        val_v[i] = Alocal_i.values;
        stream_end = stream_begin + Alocal_i.numRows();
        x_v[i] = Kokkos::subview (X_lcl, Kokkos::make_pair(stream_begin, stream_end), 0);
        y_v[i] = Kokkos::subview (Y_lcl, Kokkos::make_pair(stream_begin, stream_end), 0);;
        KernelHandle_rawptr_v_[i] = kh_v_[i].getRawPtr();
        stream_begin = stream_end;
      }
      Kokkos::fence();
      KokkosSparse::Experimental::sptrsv_solve_streams( exec_space_instances_, KernelHandle_rawptr_v_,
                                                        ptr_v, ind_v, val_v, y_v, x_v );
      Kokkos::fence();
    }
  } // End using stream interface of Kokkos Kernels Sptrsv
  else
  {
    const std::string diag = this->diag_;
    // NOTE (mfh 20 Aug 2017): KokkosSparse::trsv currently is a
    // sequential, host-only code.  See
    // https://github.com/kokkos/kokkos-kernels/issues/48. 

    auto A_lcl = this->A_crs_->getLocalMatrixHost ();

    if (X.isConstantStride () && Y.isConstantStride ()) {
      auto X_lcl = X.getLocalViewHost (Tpetra::Access::ReadWrite);
      auto Y_lcl = Y.getLocalViewHost (Tpetra::Access::ReadOnly);
      KokkosSparse::trsv (uplo.c_str (), trans.c_str (), diag.c_str (),
          A_lcl, Y_lcl, X_lcl);
    }
    else {
      for (size_t j = 0; j < numVecs; ++j) {
        auto X_j = X.getVectorNonConst (j);
        auto Y_j = Y.getVector (j);
        auto X_lcl = X_j->getLocalViewHost (Tpetra::Access::ReadWrite);
        auto Y_lcl = Y_j->getLocalViewHost (Tpetra::Access::ReadOnly);
        KokkosSparse::trsv (uplo.c_str (), trans.c_str (),
            diag.c_str (), A_lcl, Y_lcl, X_lcl);
      }
    }
  }
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

  if(isKokkosKernelsSptrsv_) os << "KK-SPTRSV, ";
  if(isKokkosKernelsStream_) os << "KK-SolveStream, ";
	
  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Matrix dimensions: ["
       << A_->getGlobalNumRows () << ", "
       << A_->getGlobalNumCols () << "]"
       << ", Number of nonzeros: " << A_->getGlobalNumEntries();
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
       A->getLocalNumRows () != A->getLocalNumCols (),
       std::runtime_error, prefix << "If A's communicator only contains one "
       "process, then A must be square.  Instead, you provided a matrix A with "
       << A->getLocalNumRows () << " rows and " << A->getLocalNumCols ()
       << " columns.");

    // It's legal for A to be null; in that case, you may not call
    // initialize() until calling setMatrix() with a nonnull input.
    // Regardless, setting the matrix invalidates the preconditioner.
    isInitialized_ = false;
    isComputed_ = false;

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

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
setStreamInfo (const bool& isKokkosKernelsStream, const int& num_streams,
               const std::vector<HandleExecSpace>& exec_space_instances)
{
  isKokkosKernelsStream_ = isKokkosKernelsStream;
  num_streams_ = num_streams;
  exec_space_instances_ = exec_space_instances;
  A_crs_v_ = std::vector< Teuchos::RCP<crs_matrix_type> >(num_streams_);
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
setMatrices (const std::vector< Teuchos::RCP<crs_matrix_type> >& A_crs_v)
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::setMatrixWithStreams: ";

  for(int i = 0; i < num_streams_; i++) {
    // If the pointer didn't change, do nothing.  This is reasonable
    // because users are supposed to call this method with the same
    // object over all participating processes, and pointer identity
    // implies object identity.
    if (A_crs_v[i].getRawPtr () != A_crs_v_[i].getRawPtr () || isInternallyChanged_) {
      // Check in serial or one-process mode if the matrix is square.
      TEUCHOS_TEST_FOR_EXCEPTION
        (! A_crs_v[i].is_null () && A_crs_v[i]->getComm ()->getSize () == 1 &&
         A_crs_v[i]->getLocalNumRows () != A_crs_v[i]->getLocalNumCols (),
         std::runtime_error, prefix << "If A's communicator only contains one "
         "process, then A must be square.  Instead, you provided a matrix A with "
         << A_crs_v[i]->getLocalNumRows () << " rows and " << A_crs_v[i]->getLocalNumCols ()
         << " columns.");
    
      // It's legal for A to be null; in that case, you may not call
      // initialize() until calling setMatrix() with a nonnull input.
      // Regardless, setting the matrix invalidates the preconditioner.
      isInitialized_ = false;
      isComputed_ = false;
    
      if (A_crs_v[i].is_null ()) {
        A_crs_v_[i] = Teuchos::null;
      }
      else { // A is not null
        Teuchos::RCP<crs_matrix_type> A_crs =
          Teuchos::rcp_dynamic_cast<crs_matrix_type> (A_crs_v[i]);
        TEUCHOS_TEST_FOR_EXCEPTION
          (A_crs.is_null (), std::invalid_argument, prefix <<
           "The input matrix A is not a Tpetra::CrsMatrix.");
        A_crs_v_[i] = A_crs;
      }
    } // pointers are not the same
  }
}

} // namespace Ifpack2

#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::LocalSparseTriangularSolver< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
