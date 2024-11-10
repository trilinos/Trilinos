// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_ILUT_DEF_HPP
#define IFPACK2_ILUT_DEF_HPP

#include <type_traits>
#include "Kokkos_StaticCrsGraph.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Time.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "KokkosSparse_par_ilut.hpp"

#include "Ifpack2_Heap.hpp"
#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_Details_getParamTryingTypes.hpp"

namespace Ifpack2 {

  namespace {

    struct IlutImplType {
      enum Enum {
        Serial,
        PAR_ILUT   //!< KokkosKernels par_ilut (parallel ILUT)
      };

      static void loadPLTypeOption (Teuchos::Array<std::string>& type_strs, Teuchos::Array<Enum>& type_enums) {
        type_strs.resize(2);
        type_strs[0] = "serial";
        type_strs[1] = "par_ilut";
        type_enums.resize(2);
        type_enums[0] = Serial;
        type_enums[1] = PAR_ILUT;
      }
    };


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
    inline typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    ilutDefaultDropTolerance () {
      typedef Teuchos::ScalarTraits<ScalarType> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;

      // 1/2.  Hopefully this can be represented in magnitude_type.
      const magnitude_type oneHalf = STM::one() / (STM::one() + STM::one());

      // The min ensures that in case magnitude_type has very low
      // precision, we'll at least get some value strictly less than
      // one.
      return std::min (static_cast<magnitude_type> (1000) * STS::magnitude (STS::eps ()), oneHalf);
    }

    // Full specialization for ScalarType = double.
    // This specialization preserves ILUT's previous default behavior.
    template<>
    inline Teuchos::ScalarTraits<double>::magnitudeType
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
  LevelOfFill_ (1.0),
  DropTolerance_ (ilutDefaultDropTolerance<scalar_type> ()),
  par_ilut_options_{1, 0., -1, -1, 0.75, false},
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false),
  useKokkosKernelsParILUT_(false)
  
{
  allocateSolvers();
}

template<class MatrixType>
void ILUT<MatrixType>::allocateSolvers ()
{
  L_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> ());
  L_solver_->setObjectLabel("lower");
  U_solver_ = Teuchos::rcp (new LocalSparseTriangularSolver<row_matrix_type> ());
  U_solver_->setObjectLabel("upper");
}

template <class MatrixType>
void ILUT<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  using Ifpack2::Details::getParamTryingTypes;
  const char prefix[] = "Ifpack2::ILUT: ";

  // Don't actually change the instance variables until we've checked
  // all parameters.  This ensures that setParameters satisfies the
  // strong exception guarantee (i.e., is transactional).

  // Parsing implementation type
  IlutImplType::Enum ilutimplType = IlutImplType::Serial;
  do {
    static const char typeName[] = "fact: type";

    if ( ! params.isType<std::string>(typeName)) break;

    // Map std::string <-> IlutImplType::Enum.
    Teuchos::Array<std::string> ilutimplTypeStrs;
    Teuchos::Array<IlutImplType::Enum> ilutimplTypeEnums;
    IlutImplType::loadPLTypeOption (ilutimplTypeStrs, ilutimplTypeEnums);
    Teuchos::StringToIntegralParameterEntryValidator<IlutImplType::Enum>
      s2i(ilutimplTypeStrs (), ilutimplTypeEnums (), typeName, false);

    ilutimplType = s2i.getIntegralValue(params.get<std::string>(typeName));
  } while (0);

  if (ilutimplType == IlutImplType::PAR_ILUT) {
    this->useKokkosKernelsParILUT_ = true;
  }
  else {
    this->useKokkosKernelsParILUT_ = false;
  }

  // Fill level in ILUT is a double, not a magnitude_type, because it
  // depends on LO and GO, not on Scalar.  Also, you can't cast
  // arbitrary magnitude_type (e.g., Sacado::MP::Vector) to double.
  double fillLevel = LevelOfFill_;
  {
    const std::string paramName ("fact: ilut level-of-fill");
    TEUCHOS_TEST_FOR_EXCEPTION(
      (params.isParameter(paramName) && this->useKokkosKernelsParILUT_), std::runtime_error,
      "Ifpack2::ILUT: Parameter " << paramName << " is meaningless for algorithm par_ilut.");
    getParamTryingTypes<double, double, float>
      (fillLevel, params, paramName, prefix);
    TEUCHOS_TEST_FOR_EXCEPTION
      (fillLevel < 1.0, std::runtime_error,
       "Ifpack2::ILUT: The \"" << paramName << "\" parameter must be >= "
       "1.0, but you set it to " << fillLevel << ".  For ILUT, the fill level "
       "means something different than it does for ILU(k).  ILU(0) produces "
       "factors with the same sparsity structure as the input matrix A. For "
       "ILUT, level-of-fill = 1.0 will produce factors with nonzeros matching "
       "the sparsity structure of A. level-of-fill > 1.0 allows for additional "
       "fill-in.");
  }

  magnitude_type absThresh = Athresh_;
  {
    const std::string paramName ("fact: absolute threshold");
    getParamTryingTypes<magnitude_type, magnitude_type, double>
      (absThresh, params, paramName, prefix);
  }

  magnitude_type relThresh = Rthresh_;
  {
    const std::string paramName ("fact: relative threshold");
    getParamTryingTypes<magnitude_type, magnitude_type, double>
      (relThresh, params, paramName, prefix);
  }

  magnitude_type relaxValue = RelaxValue_;
  {
    const std::string paramName ("fact: relax value");
    getParamTryingTypes<magnitude_type, magnitude_type, double>
      (relaxValue, params, paramName, prefix);
  }

  magnitude_type dropTol = DropTolerance_;
  {
    const std::string paramName ("fact: drop tolerance");
    getParamTryingTypes<magnitude_type, magnitude_type, double>
      (dropTol, params, paramName, prefix);
  }

  int par_ilut_max_iter=20;
  magnitude_type par_ilut_residual_norm_delta_stop=1e-2;
  int par_ilut_team_size=0;
  int par_ilut_vector_size=0;
  float par_ilut_fill_in_limit=0.75;
  bool par_ilut_verbose=false;
  if (this->useKokkosKernelsParILUT_) {
    par_ilut_max_iter = par_ilut_options_.max_iter;
    par_ilut_residual_norm_delta_stop = par_ilut_options_.residual_norm_delta_stop;
    par_ilut_team_size = par_ilut_options_.team_size;
    par_ilut_vector_size = par_ilut_options_.vector_size;
    par_ilut_fill_in_limit = par_ilut_options_.fill_in_limit;
    par_ilut_verbose = par_ilut_options_.verbose;

    std::string par_ilut_plist_name("parallel ILUT options");
    if (params.isSublist(par_ilut_plist_name)) {
      Teuchos::ParameterList const &par_ilut_plist = params.sublist(par_ilut_plist_name);

      std::string paramName("maximum iterations");
      getParamTryingTypes<int, int>(par_ilut_max_iter, par_ilut_plist, paramName, prefix);

      paramName = "residual norm delta stop";
      getParamTryingTypes<magnitude_type, magnitude_type, double>(par_ilut_residual_norm_delta_stop, par_ilut_plist, paramName, prefix);

      paramName = "team size";
      getParamTryingTypes<int, int>(par_ilut_team_size, par_ilut_plist, paramName, prefix);

      paramName = "vector size";
      getParamTryingTypes<int, int>(par_ilut_vector_size, par_ilut_plist, paramName, prefix);

      paramName = "fill in limit";
      getParamTryingTypes<float, float, double>(par_ilut_fill_in_limit, par_ilut_plist, paramName, prefix);

      paramName = "verbose";
      getParamTryingTypes<bool, bool>(par_ilut_verbose, par_ilut_plist, paramName, prefix);

    } // if (params.isSublist(par_ilut_plist_name))

    par_ilut_options_.max_iter = par_ilut_max_iter;
    par_ilut_options_.residual_norm_delta_stop = par_ilut_residual_norm_delta_stop;
    par_ilut_options_.team_size = par_ilut_team_size;
    par_ilut_options_.vector_size = par_ilut_vector_size;
    par_ilut_options_.fill_in_limit = par_ilut_fill_in_limit;
    par_ilut_options_.verbose = par_ilut_verbose;

  } //if (this->useKokkosKernelsParILUT_)

  // Forward to trisolvers.
  L_solver_->setParameters(params);
  U_solver_->setParameters(params);

  LevelOfFill_ = fillLevel;
  Athresh_ = absThresh;
  Rthresh_ = relThresh;
  RelaxValue_ = relaxValue;
  DropTolerance_ = dropTol;
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
size_t ILUT<MatrixType>::getNodeSmootherComplexity() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getNodeSmootherComplexity: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  // ILUT methods cost roughly one apply + the nnz in the upper+lower triangles
  return A_->getLocalNumEntries() + getLocalNumEntries();
}


template<class MatrixType>
global_size_t ILUT<MatrixType>::getGlobalNumEntries () const {
  return L_->getGlobalNumEntries () + U_->getGlobalNumEntries ();
}


template<class MatrixType>
size_t ILUT<MatrixType>::getLocalNumEntries () const {
  return L_->getLocalNumEntries () + U_->getLocalNumEntries ();
}


template<class MatrixType>
void ILUT<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != A_.getRawPtr ()) {
    // Check in serial or one-process mode if the matrix is square.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! A.is_null () && A->getComm ()->getSize () == 1 &&
      A->getLocalNumRows () != A->getLocalNumCols (),
      std::runtime_error, "Ifpack2::ILUT::setMatrix: If A's communicator only "
      "contains one process, then A must be square.  Instead, you provided a "
      "matrix A with " << A->getLocalNumRows () << " rows and "
      << A->getLocalNumCols () << " columns.");

    // It's legal for A to be null; in that case, you may not call
    // initialize() until calling setMatrix() with a nonnull input.
    // Regardless, setting the matrix invalidates any previous
    // factorization.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;

    // The sparse triangular solvers get a triangular factor as their
    // input matrix.  The triangular factors L_ and U_ are getting
    // reset, so we reset the solvers' matrices to null.  Do that
    // before setting L_ and U_ to null, so that latter step actually
    // frees the factors.
    if (! L_solver_.is_null ()) {
      L_solver_->setMatrix (Teuchos::null);
    }
    if (! U_solver_.is_null ()) {
      U_solver_->setMatrix (Teuchos::null);
    }

    L_ = Teuchos::null;
    U_ = Teuchos::null;
    A_ = A;
  }
}

template <class MatrixType>
Teuchos::RCP<const typename ILUT<MatrixType>::row_matrix_type>
ILUT<MatrixType>::makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;

  // If A_'s communicator only has one process, or if its column and
  // row Maps are the same, then it is already local, so use it
  // directly.
  if (A->getRowMap ()->getComm ()->getSize () == 1 ||
      A->getRowMap ()->isSameAs (* (A->getColMap ()))) {
    return A;
  }

  // If A_ is already a LocalFilter, then use it directly.  This
  // should be the case if RILUT is being used through
  // AdditiveSchwarz, for example.
  RCP<const LocalFilter<row_matrix_type> > A_lf_r =
    rcp_dynamic_cast<const LocalFilter<row_matrix_type> > (A);
  if (! A_lf_r.is_null ()) {
    return rcp_implicit_cast<const row_matrix_type> (A_lf_r);
  }
  else {
    // A_'s communicator has more than one process, its row Map and
    // its column Map differ, and A_ is not a LocalFilter.  Thus, we
    // have to wrap it in a LocalFilter.
    return rcp (new LocalFilter<row_matrix_type> (A));
  }
}


template<class MatrixType>
void ILUT<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::rcp_const_cast;
  Teuchos::Time timer ("ILUT::initialize");
  double startTime = timer.wallTime();
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

    A_local_ = makeLocalFilter(A_); // Compute the local filter.
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local_.is_null(), std::logic_error, "Ifpack2::RILUT::initialize: "
      "makeLocalFilter returned null; it failed to compute A_local.  "
      "Please report this bug to the Ifpack2 developers.");

    if (this->useKokkosKernelsParILUT_) {
      this->KernelHandle_ = Teuchos::rcp(new kk_handle_type());
      KernelHandle_->create_par_ilut_handle();
      auto par_ilut_handle = KernelHandle_->get_par_ilut_handle();
      par_ilut_handle->set_residual_norm_delta_stop(par_ilut_options_.residual_norm_delta_stop);
      par_ilut_handle->set_team_size(par_ilut_options_.team_size);
      par_ilut_handle->set_vector_size(par_ilut_options_.vector_size);
      par_ilut_handle->set_fill_in_limit(par_ilut_options_.fill_in_limit);
      par_ilut_handle->set_verbose(par_ilut_options_.verbose);
      par_ilut_handle->set_async_update(false);

      RCP<const crs_matrix_type> A_local_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A_local_);
      if (A_local_crs.is_null()) {
         // the result is a host-based matrix, which is the same as what happens in RILUK
        local_ordinal_type numRows = A_local_->getLocalNumRows();
        Array<size_t> entriesPerRow(numRows);
        for(local_ordinal_type i = 0; i < numRows; i++) {
          entriesPerRow[i] = A_local_->getNumEntriesInLocalRow(i);
        }
        RCP<crs_matrix_type> A_local_crs_nc =
          rcp (new crs_matrix_type (A_local_->getRowMap (),
                                    A_local_->getColMap (),
                                    entriesPerRow()));
        // copy entries into A_local_crs
        nonconst_local_inds_host_view_type indices("indices",A_local_->getLocalMaxNumRowEntries());
        nonconst_values_host_view_type values("values",A_local_->getLocalMaxNumRowEntries());
        for(local_ordinal_type i = 0; i < numRows; i++) {
          size_t numEntries = 0;
          A_local_->getLocalRowCopy(i, indices, values, numEntries);
          A_local_crs_nc->insertLocalValues(i, numEntries, reinterpret_cast<scalar_type*>(values.data()), indices.data());
        }
        A_local_crs_nc->fillComplete (A_local_->getDomainMap (), A_local_->getRangeMap ());
        A_local_crs = rcp_const_cast<const crs_matrix_type> (A_local_crs_nc);
      }
      auto A_local_crs_device = A_local_crs->getLocalMatrixDevice();

      //KokkosKernels requires unsigned
      typedef typename Kokkos::View<usize_type*, array_layout, device_type> ulno_row_view_t;
      const int NumMyRows = A_local_crs->getRowMap()->getLocalNumElements();
      L_rowmap_ = ulno_row_view_t("L_row_map", NumMyRows + 1);
      U_rowmap_ = ulno_row_view_t("U_row_map", NumMyRows + 1);
      L_rowmap_orig_ = ulno_row_view_t("L_row_map_orig", NumMyRows + 1);
      U_rowmap_orig_ = ulno_row_view_t("U_row_map_orig", NumMyRows + 1);

      KokkosSparse::Experimental::par_ilut_symbolic(KernelHandle_.getRawPtr(),
                                                    A_local_crs_device.graph.row_map, A_local_crs_device.graph.entries,
                                                    L_rowmap_,
                                                    U_rowmap_);

      Kokkos::deep_copy(L_rowmap_orig_, L_rowmap_);
      Kokkos::deep_copy(U_rowmap_orig_, U_rowmap_);
    }

    IsInitialized_ = true;
    ++NumInitialize_;
  } //timer scope
  InitializeTime_ += (timer.wallTime() - startTime);
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
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("ILUT::compute");
  double startTime = timer.wallTime();
  { // Timer scope for timing compute()
  Teuchos::TimeMonitor timeMon (timer, true);

  if (!this->useKokkosKernelsParILUT_)
  {
    //--------------------------------------------------------------------------
    // Ifpack2::ILUT's serial version is a translation of the Aztec ILUT
    // implementation. The Aztec ILUT implementation was written by Ray Tuminaro.
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

    const scalar_type zero = STS::zero ();
    const scalar_type one  = STS::one ();

    const local_ordinal_type myNumRows = A_local_->getLocalNumRows ();

    // If this macro is defined, files containing the L and U factors
    // will be written. DON'T CHECK IN THE CODE WITH THIS MACRO ENABLED!!!
    // #define IFPACK2_WRITE_ILUT_FACTORS
#ifdef IFPACK2_WRITE_ILUT_FACTORS
    std::ofstream ofsL("L.ifpack2_ilut.mtx", std::ios::out);
    std::ofstream ofsU("U.ifpack2_ilut.mtx", std::ios::out);
#endif

    // Calculate how much fill will be allowed in addition to the
    // space that corresponds to the input matrix entries.
    double local_nnz = static_cast<double> (A_local_->getLocalNumEntries ());
    double fill = ((getLevelOfFill () - 1.0) * local_nnz) / (2 * myNumRows);

    // std::ceil gives the smallest integer larger than the argument.
    // this may give a slightly different result than Aztec's fill value in
    // some cases.
    double fill_ceil=std::ceil(fill);

    // Similarly to Aztec, we will allow the same amount of fill for each
    // row, half in L and half in U.
    size_type fillL = static_cast<size_type>(fill_ceil);
    size_type fillU = static_cast<size_type>(fill_ceil);

    Array<scalar_type> InvDiagU (myNumRows, zero);

    Array<Array<local_ordinal_type> > L_tmp_idx(myNumRows);
    Array<Array<scalar_type> > L_tmpv(myNumRows);
    Array<Array<local_ordinal_type> > U_tmp_idx(myNumRows);
    Array<Array<scalar_type> > U_tmpv(myNumRows);

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
    nonconst_local_inds_host_view_type ColIndicesARCP;
    nonconst_values_host_view_type ColValuesARCP;
    if (! A_local_->supportsRowViews ()) {
      const size_t maxnz = A_local_->getLocalMaxNumRowEntries ();
      Kokkos::resize(ColIndicesARCP,maxnz);
      Kokkos::resize(ColValuesARCP,maxnz);
    }

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      local_inds_host_view_type  ColIndicesA;
      values_host_view_type ColValuesA;
      size_t RowNnz;

      if (A_local_->supportsRowViews ()) {
        A_local_->getLocalRowView (row_i, ColIndicesA, ColValuesA);
        RowNnz = ColIndicesA.size ();
      }
      else {
        A_local_->getLocalRowCopy (row_i, ColIndicesARCP, ColValuesARCP, RowNnz);
        ColIndicesA = Kokkos::subview(ColIndicesARCP,std::make_pair((size_t)0, RowNnz));
        ColValuesA  = Kokkos::subview(ColValuesARCP,std::make_pair((size_t)0, RowNnz));
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

        ArrayView<local_ordinal_type> ColIndicesU = U_tmp_idx[row_k]();
        ArrayView<scalar_type> ColValuesU = U_tmpv[row_k]();
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
      for (size_type i = 0; i < (size_type)ColIndicesA.size (); ++i) {
        if (ColIndicesA[i] < row_i) {
          L_tmp_idx[row_i].push_back(ColIndicesA[i]);
          L_tmpv[row_i].push_back(cur_row[ColIndicesA[i]]);
          pattern[ColIndicesA[i]] = UNUSED;
        }
      }

      //   next, the L entries resulting from fill:
      for (size_type j = 0; j < L_vals_heaplen; ++j) {
        L_tmp_idx[row_i].push_back(L_vals_heap[j]);
        L_tmpv[row_i].push_back(cur_row[L_vals_heap[j]]);
        pattern[L_vals_heap[j]] = UNUSED;
      }

      // L has a one on the diagonal, but we don't explicitly store
      // it.  If we don't store it, then the kernel which performs the
      // triangular solve can assume a unit diagonal, take a short-cut
      // and perform faster.

#ifdef IFPACK2_WRITE_ILUT_FACTORS
      for (size_type ii = 0; ii < L_tmp_idx[row_i].size (); ++ii) {
        ofsL << row_i << " " << L_tmp_idx[row_i][ii] << " " 
                             << L_tmpv[row_i][ii] << std::endl;
      }
#endif


      // Pick out the diagonal element, store its reciprocal.
      if (cur_row[row_i] == zero) {
        std::cerr << "Ifpack2::ILUT::Compute: zero pivot encountered! "
                  << "Replacing with rownorm and continuing..."
                  << "(You may need to set the parameter "
                  << "'fact: absolute threshold'.)" << std::endl;
        cur_row[row_i] = rownorm;
      }
      InvDiagU[row_i] = one / cur_row[row_i];

      // Non-inverted diagonal is stored for U:
      U_tmp_idx[row_i].push_back(row_i);
      U_tmpv[row_i].push_back(cur_row[row_i]);
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
          U_tmp_idx[row_i].push_back(col);
          U_tmpv[row_i].push_back(cur_row[col]);
          unorm[row_i] += scalar_mag(cur_row[col]);
        }
        pattern[col] = UNUSED;
      }

      for(size_type j=0; j<U_vals_heaplen; ++j) {
        U_tmp_idx[row_i].push_back(U_vals_heap[j]);
        U_tmpv[row_i].push_back(cur_row[U_vals_heap[j]]);
        unorm[row_i] += scalar_mag(cur_row[U_vals_heap[j]]);
      }

      unorm[row_i] /= (orig_U_len + U_vals_heaplen);

#ifdef IFPACK2_WRITE_ILUT_FACTORS
      for(int ii=0; ii<U_tmp_idx[row_i].size(); ++ii) {
        ofsU <<row_i<< " " <<U_tmp_idx[row_i][ii]<< " " 
                           <<U_tmpv[row_i][ii]<< std::endl;
      }
#endif

      L_cols_heap.clear();
      U_cols.clear();
      L_vals_heap.clear();
      U_vals_heap.clear();
    } // end of for(row_i) loop

    // Now allocate and fill the matrices
    Array<size_t> nnzPerRow(myNumRows);

    // Make sure to release the old memory for L & U prior to recomputing to
    // avoid bloating the high-water mark.
    L_ = Teuchos::null;
    U_ = Teuchos::null;
    L_solver_->setMatrix(Teuchos::null);
    U_solver_->setMatrix(Teuchos::null);

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      nnzPerRow[row_i] = L_tmp_idx[row_i].size();
    }

    L_ = rcp (new crs_matrix_type (A_local_->getRowMap(), A_local_->getColMap(),
                                   nnzPerRow()));

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      L_->insertLocalValues (row_i, L_tmp_idx[row_i](), L_tmpv[row_i]());
    }

    L_->fillComplete(); 

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      nnzPerRow[row_i] = U_tmp_idx[row_i].size();
    }

    U_ = rcp (new crs_matrix_type (A_local_->getRowMap(), A_local_->getColMap(),
                                   nnzPerRow()));

    for (local_ordinal_type row_i = 0 ; row_i < myNumRows ; ++row_i) {
      U_->insertLocalValues (row_i, U_tmp_idx[row_i](), U_tmpv[row_i]());
    }

    U_->fillComplete();

    L_solver_->setMatrix(L_);
    L_solver_->initialize ();
    L_solver_->compute ();

    U_solver_->setMatrix(U_);
    U_solver_->initialize ();
    U_solver_->compute ();

  } //if (!this->useKokkosKernelsParILUT_)
  else {
    // Set L, U rowmaps back to original state. Par_ilut can change them, which invalidates them
    // if compute is called again.
    if (this->isComputed()) {
      Kokkos::resize(L_rowmap_, L_rowmap_orig_.size());
      Kokkos::resize(U_rowmap_, U_rowmap_orig_.size());
      Kokkos::deep_copy(L_rowmap_, L_rowmap_orig_);
      Kokkos::deep_copy(U_rowmap_, U_rowmap_orig_);
    }

    RCP<const crs_matrix_type> A_local_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A_local_);
    {//Make sure values in A is picked up even in case of pattern reuse
      if(A_local_crs.is_null()) {
        local_ordinal_type numRows = A_local_->getLocalNumRows();
        Array<size_t> entriesPerRow(numRows);
        for(local_ordinal_type i = 0; i < numRows; i++) {
          entriesPerRow[i] = A_local_->getNumEntriesInLocalRow(i);
        }
        RCP<crs_matrix_type> A_local_crs_nc =
          rcp (new crs_matrix_type (A_local_->getRowMap (),
                                    A_local_->getColMap (),
                                    entriesPerRow()));
        // copy entries into A_local_crs
        nonconst_local_inds_host_view_type indices("indices",A_local_->getLocalMaxNumRowEntries());
        nonconst_values_host_view_type values("values",A_local_->getLocalMaxNumRowEntries());
        for(local_ordinal_type i = 0; i < numRows; i++) {
          size_t numEntries = 0;
          A_local_->getLocalRowCopy(i, indices, values, numEntries);
          A_local_crs_nc->insertLocalValues(i, numEntries, reinterpret_cast<scalar_type*>(values.data()),indices.data());
        }
        A_local_crs_nc->fillComplete (A_local_->getDomainMap (), A_local_->getRangeMap ());
        A_local_crs = rcp_const_cast<const crs_matrix_type> (A_local_crs_nc);
      }
      auto lclMtx = A_local_crs->getLocalMatrixDevice();
      A_local_rowmap_  = lclMtx.graph.row_map;
      A_local_entries_ = lclMtx.graph.entries;
      A_local_values_  = lclMtx.values;
    }

    //JHU TODO Should allocation of L & U's column (aka entry) and value arrays occur here or in init()?
    auto par_ilut_handle = KernelHandle_->get_par_ilut_handle();
    auto nnzL = par_ilut_handle->get_nnzL();
    static_graph_entries_t L_entries_ = static_graph_entries_t("L_entries", nnzL);
    local_matrix_values_t L_values_ = local_matrix_values_t("L_values", nnzL);

    auto nnzU = par_ilut_handle->get_nnzU();
    static_graph_entries_t U_entries_ = static_graph_entries_t("U_entries", nnzU);
    local_matrix_values_t U_values_ = local_matrix_values_t("U_values", nnzU);

    KokkosSparse::Experimental::par_ilut_numeric(KernelHandle_.getRawPtr(),
                                                   A_local_rowmap_, A_local_entries_, A_local_values_,
                                                   L_rowmap_, L_entries_, L_values_, U_rowmap_, U_entries_, U_values_);

    auto L_kokkosCrsGraph = local_graph_device_type(L_entries_, L_rowmap_);
    auto U_kokkosCrsGraph = local_graph_device_type(U_entries_, U_rowmap_);

    local_matrix_device_type L_localCrsMatrix_device;
    L_localCrsMatrix_device = local_matrix_device_type("L_Factor_localmatrix",
                                                       A_local_->getLocalNumRows(),
                                                       L_values_,
                                                       L_kokkosCrsGraph);

    L_ = rcp (new crs_matrix_type (L_localCrsMatrix_device,
                                   A_local_crs->getRowMap(),
                                   A_local_crs->getColMap(),
                                   A_local_crs->getDomainMap(),
                                   A_local_crs->getRangeMap(),
                                   A_local_crs->getGraph()->getImporter(),
                                   A_local_crs->getGraph()->getExporter()));

    local_matrix_device_type U_localCrsMatrix_device;
    U_localCrsMatrix_device = local_matrix_device_type("U_Factor_localmatrix",
                                                       A_local_->getLocalNumRows(),
                                                       U_values_,
                                                       U_kokkosCrsGraph);

    U_ = rcp (new crs_matrix_type (U_localCrsMatrix_device,
                                   A_local_crs->getRowMap(),
                                   A_local_crs->getColMap(),
                                   A_local_crs->getDomainMap(),
                                   A_local_crs->getRangeMap(),
                                   A_local_crs->getGraph()->getImporter(),
                                   A_local_crs->getGraph()->getExporter()));

    L_solver_->setMatrix (L_);
    L_solver_->compute ();//NOTE: Only do compute if the pointer changed. Otherwise, do nothing
    U_solver_->setMatrix (U_);
    U_solver_->compute ();//NOTE: Only do compute if the pointer changed. Otherwise, do nothing
  } //if (!this->useKokkosKernelsParILUT_) ... else ...

  } // Timer scope for timing compute()
  ComputeTime_ += (timer.wallTime() - startTime);
  IsComputed_ = true;
  ++NumCompute_;
} //compute()


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

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::ILUT::apply: You must call compute() to compute the incomplete "
    "factorization, before calling apply().");

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::ILUT::apply: X and Y must have the same number of columns.  "
    "X has " << X.getNumVectors () << " columns, but Y has "
    << Y.getNumVectors () << " columns.");

  const scalar_type one = STS::one ();
  const scalar_type zero = STS::zero ();

  Teuchos::Time timer ("ILUT::apply");
  double startTime = timer.wallTime();
  { // Start timing
    Teuchos::TimeMonitor timeMon (timer, true);

    if (alpha == one && beta == zero) {
      if (mode == Teuchos::NO_TRANS) { // Solve L (U Y) = X for Y.
        // Start by solving L Y = X for Y.
        L_solver_->apply (X, Y, mode);

        // Solve U Y = Y.
        U_solver_->apply (Y, Y, mode);
      }
      else { // Solve U^P (L^P Y)) = X for Y (where P is * or T).

        // Start by solving U^P Y = X for Y.
        U_solver_->apply (X, Y, mode);

        // Solve L^P Y = Y.
        L_solver_->apply (Y, Y, mode);
      }
    }
    else { // alpha != 1 or beta != 0
      if (alpha == zero) {
        // The special case for beta == 0 ensures that if Y contains Inf
        // or NaN values, we replace them with 0 (following BLAS
        // convention), rather than multiplying them by 0 to get NaN.
        if (beta == zero) {
          Y.putScalar (zero);
        } else {
          Y.scale (beta);
        }
      } else { // alpha != zero
        MV Y_tmp (Y.getMap (), Y.getNumVectors ());
        apply (X, Y_tmp, mode);
        Y.update (alpha, Y_tmp, beta);
      }
    }
  }//end timing

  ++NumApply_;
  ApplyTime_ += (timer.wallTime() - startTime);
} //apply()


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

}//namespace Ifpack2


// FIXME (mfh 16 Sep 2014) We should really only use RowMatrix here!
// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_ILUT_INSTANT(S,LO,GO,N) \
  template class Ifpack2::ILUT< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif /* IFPACK2_ILUT_DEF_HPP */
