// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_SolverCore_def.hpp
  \class  Amesos2::SolverCore
  \author Eric T Bavier <etbavie@sandia.gov>
  \date   Thu May 27 14:02:35 CDT 2010

  \brief  Templated core-functionality class for Amesos2 solvers.  Definitions.
*/

#ifndef AMESOS2_SOLVERCORE_DEF_HPP
#define AMESOS2_SOLVERCORE_DEF_HPP

#include "Kokkos_ArithTraits.hpp"

#include "Amesos2_MatrixAdapter_def.hpp"
#include "Amesos2_MultiVecAdapter_def.hpp"

#include "Amesos2_Util.hpp"

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas.hpp"

namespace Amesos2 {


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverCore<ConcreteSolver,Matrix,Vector>::SolverCore(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : matrixA_(createConstMatrixAdapter<Matrix>(A))
  , multiVecX_(X) // may be null
  , multiVecB_(B) // may be null
  , globalNumRows_(matrixA_->getGlobalNumRows())
  , globalNumCols_(matrixA_->getGlobalNumCols())
  , globalNumNonZeros_(matrixA_->getGlobalNNZ())
  , rowIndexBase_(matrixA_->getRowIndexBase())
  , columnIndexBase_(matrixA_->getColumnIndexBase())
  , rank_(Teuchos::rank(*this->getComm()))
  , root_(rank_ == 0)
  , nprocs_(Teuchos::size(*this->getComm()))
{
    TEUCHOS_TEST_FOR_EXCEPTION(
    !matrixShapeOK(),
    std::invalid_argument,
    "Matrix shape inappropriate for this solver");
}


/// Destructor
template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverCore<ConcreteSolver,Matrix,Vector>::~SolverCore( )
{
  // Nothing
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::preOrdering()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  loadA(PREORDERING);

  int error_code = static_cast<solver_type*>(this)->preOrdering_impl();
  if (error_code == EXIT_SUCCESS){
    ++status_.numPreOrder_;
    status_.last_phase_ = PREORDERING;
  }

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::symbolicFactorization()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( !status_.preOrderingDone() ){
    preOrdering();
    if( !matrix_loaded_ ) loadA(SYMBFACT);
  } else {
    loadA(SYMBFACT);
  }

  int error_code = static_cast<solver_type*>(this)->symbolicFactorization_impl();
  if (error_code == EXIT_SUCCESS){
    ++status_.numSymbolicFact_;
    status_.last_phase_ = SYMBFACT;
  }

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::numericFactorization()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( !status_.symbolicFactorizationDone() ){
    symbolicFactorization();
    if( !matrix_loaded_ ) loadA(NUMFACT);
  } else {
    loadA(NUMFACT);
  }

  int error_code = static_cast<solver_type*>(this)->numericFactorization_impl();
  if (error_code == EXIT_SUCCESS){
    ++status_.numNumericFact_;
    status_.last_phase_ = NUMFACT;
  }

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve()
{
  solve(multiVecX_.ptr(), multiVecB_.ptr());
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve(const Teuchos::Ptr<Vector> X,
                                                const Teuchos::Ptr<const Vector> B) const
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  X.assert_not_null();
  B.assert_not_null();

  if (control_.useIterRefine_) {
    solve_ir(X, B, control_.maxNumIterRefines_, control_.verboseIterRefine_);
    return;
  }

  const Teuchos::RCP<MultiVecAdapter<Vector> > x =
    createMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(X));
  const Teuchos::RCP<const MultiVecAdapter<Vector> > b =
    createConstMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(B));

#ifdef HAVE_AMESOS2_DEBUG
  // Check some required properties of X and B
  TEUCHOS_TEST_FOR_EXCEPTION
    (x->getGlobalLength() != matrixA_->getGlobalNumCols(),
     std::invalid_argument,
     "MultiVector X must have length equal to the number of "
     "global columns in A.  X->getGlobalLength() = "
     << x->getGlobalLength() << " != A->getGlobalNumCols() = "
     << matrixA_->getGlobalNumCols() << ".");

  TEUCHOS_TEST_FOR_EXCEPTION(b->getGlobalLength() != matrixA_->getGlobalNumRows(),
                     std::invalid_argument,
                     "MultiVector B must have length equal to the number of "
                     "global rows in A");

  TEUCHOS_TEST_FOR_EXCEPTION(x->getGlobalNumVectors() != b->getGlobalNumVectors(),
                     std::invalid_argument,
                     "X and B MultiVectors must have the same number of vectors");
#endif  // HAVE_AMESOS2_DEBUG

  if( !status_.numericFactorizationDone() ){
    // This casting-away of constness is probably OK because this
    // function is meant to be "logically const"
    const_cast<type&>(*this).numericFactorization();
  }

  int error_code = static_cast<const solver_type*>(this)->solve_impl(Teuchos::outArg(*x), Teuchos::ptrInArg(*b));
  if (error_code == EXIT_SUCCESS){
    ++status_.numSolve_;
    status_.last_phase_ = SOLVE;
  }
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve(Vector* X, const Vector* B) const
{
  solve(Teuchos::ptr(X), Teuchos::ptr(B));
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
int
SolverCore<ConcreteSolver,Matrix,Vector>::solve_ir(const int maxNumIters, const bool verbose)
{
  return solve_ir(multiVecX_.ptr(), multiVecB_.ptr(), maxNumIters, verbose);
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
int
SolverCore<ConcreteSolver,Matrix,Vector>::solve_ir(Vector* X, const Vector* B, const int maxNumIters, const bool verbose) const
{
  return solve_ir(Teuchos::ptr(X), Teuchos::ptr(B), maxNumIters, verbose);
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
int
SolverCore<ConcreteSolver,Matrix,Vector>::solve_ir(const Teuchos::Ptr<      Vector> x,
                                                   const Teuchos::Ptr<const Vector> b,
                                                   const int maxNumIters,
                                                   const bool verbose) const
{
  using KAT              = Kokkos::ArithTraits<scalar_type>;
  using impl_scalar_type = typename KAT::val_type;
  using magni_type       = typename KAT::mag_type;
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  using host_crsmat_t    = KokkosSparse::CrsMatrix<impl_scalar_type, int, host_execution_space, void, int>;
  using host_graph_t     = typename host_crsmat_t::StaticCrsGraphType;
  using host_values_t    = typename host_crsmat_t::values_type::non_const_type;
  using host_row_map_t   = typename host_graph_t::row_map_type::non_const_type;
  using host_colinds_t   = typename host_graph_t::entries_type::non_const_type;
  using host_mvector_t   = Kokkos::View<impl_scalar_type **, Kokkos::LayoutLeft, host_execution_space>;
  using host_vector_t    = Kokkos::View<impl_scalar_type *,  Kokkos::LayoutLeft, host_execution_space>;
  using host_magni_view  = Kokkos::View<magni_type  *,  Kokkos::LayoutLeft, host_execution_space>;

  const impl_scalar_type one(1.0);
  const impl_scalar_type mone = impl_scalar_type(-one);
  const magni_type eps = KAT::eps ();

  // get data needed for IR
  using MVAdapter = MultiVecAdapter<Vector>;
  Teuchos::RCP<      MVAdapter> X = createMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(x));
  Teuchos::RCP<const MVAdapter> B = createConstMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(b));

  auto r_ = B->clone();
  auto e_ = X->clone();
  auto r = r_.ptr();
  auto e = e_.ptr();
  Teuchos::RCP<      MVAdapter> R = createMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(r));
  Teuchos::RCP<      MVAdapter> E = createMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(e));

  const size_t nrhs = X->getGlobalNumVectors();
  const int nnz   = this->matrixA_->getGlobalNNZ();
  const int nrows = this->matrixA_->getGlobalNumRows();

  // get local matriix
  host_crsmat_t crsmat;
  host_graph_t static_graph;
  host_row_map_t rowmap_view;
  host_colinds_t colind_view;
  host_values_t  values_view;
  if (this->root_) {
    Kokkos::resize(rowmap_view, 1+nrows);
    Kokkos::resize(colind_view, nnz);
    Kokkos::resize(values_view, nnz);
  } else {
    Kokkos::resize(rowmap_view, 1);
    Kokkos::resize(colind_view, 0);
    Kokkos::resize(values_view, 0);
  }

  int nnz_ret = 0;
  Util::get_crs_helper_kokkos_view<
    MatrixAdapter<Matrix>, host_values_t, host_colinds_t, host_row_map_t>::do_get(
      this->matrixA_.ptr(),
      values_view, colind_view, rowmap_view,
      nnz_ret, ROOTED, ARBITRARY, this->rowIndexBase_);

  if (this->root_) {
    static_graph = host_graph_t(colind_view, rowmap_view);
    crsmat = host_crsmat_t("CrsMatrix", nrows, values_view, static_graph);
  }

  //
  // ** First Solve **
  static_cast<const solver_type*>(this)->solve_impl(Teuchos::outArg(*X), Teuchos::ptrInArg(*B));


  // auxiliary scalar Kokkos views
  const int ldx = (this->root_ ? X->getGlobalLength() : 0);
  const int ldb = (this->root_ ? B->getGlobalLength() : 0);
  const int ldr = (this->root_ ? R->getGlobalLength() : 0);
  const int lde = (this->root_ ? E->getGlobalLength() : 0);
  const bool     initialize_data = true;
  const bool not_initialize_data = true;
  host_mvector_t X_view;
  host_mvector_t B_view;
  host_mvector_t R_view;
  host_mvector_t E_view;

  global_size_type rowIndexBase = this->rowIndexBase_;
  auto Xptr = Teuchos::Ptr<      MVAdapter>(X.ptr());
  auto Bptr = Teuchos::Ptr<const MVAdapter>(B.ptr());
  auto Rptr = Teuchos::Ptr<      MVAdapter>(R.ptr());
  auto Eptr = Teuchos::Ptr<      MVAdapter>(E.ptr());
  Util::get_1d_copy_helper_kokkos_view<MVAdapter, host_mvector_t>::
    do_get(    initialize_data, Xptr, X_view, ldx, CONTIGUOUS_AND_ROOTED, rowIndexBase);
  Util::get_1d_copy_helper_kokkos_view<MVAdapter, host_mvector_t>::
    do_get(    initialize_data, Bptr, B_view, ldb, CONTIGUOUS_AND_ROOTED, rowIndexBase);
  Util::get_1d_copy_helper_kokkos_view<MVAdapter, host_mvector_t>::
    do_get(not_initialize_data, Rptr, R_view, ldr, CONTIGUOUS_AND_ROOTED, rowIndexBase);
  Util::get_1d_copy_helper_kokkos_view<MVAdapter, host_mvector_t>::
    do_get(not_initialize_data, Eptr, E_view, lde, CONTIGUOUS_AND_ROOTED, rowIndexBase);


  host_magni_view x0norms("x0norms", nrhs);
  host_magni_view bnorms("bnorms", nrhs);
  host_magni_view enorms("enorms", nrhs);
  if (this->root_) {
    // compute initial solution norms (used for stopping criteria)
    for (size_t j = 0; j < nrhs; j++) { 
      auto x_subview = Kokkos::subview(X_view, Kokkos::ALL(), j);
      host_vector_t x_1d (const_cast<impl_scalar_type*>(x_subview.data()), x_subview.extent(0));
      x0norms(j) = KokkosBlas::nrm2(x_1d);
    }
    if (verbose) {
      std::cout << std::endl
                << " SolverCore :: solve_ir (maxNumIters = " << maxNumIters
                << ", tol = " << x0norms(0) << " * " << eps << " = " << x0norms(0)*eps
                << ")" << std::endl;
    }

    // compute residual norm
    if (verbose) {
      std::cout << " bnorm = ";
      for (size_t j = 0; j < nrhs; j++) { 
        auto b_subview = Kokkos::subview(B_view, Kokkos::ALL(), j);
        host_vector_t b_1d (const_cast<impl_scalar_type*>(b_subview.data()), b_subview.extent(0));
        bnorms(j) = KokkosBlas::nrm2(b_1d);
        std::cout << bnorms(j) << ", ";
      }
      std::cout << std::endl;
    }
  }


  //
  // ** Iterative Refinement **
  int numIters = 0;
  int converged = 0; // 0 = has not converged, 1 = converged
  for (numIters = 0; numIters < maxNumIters && converged == 0; ++numIters) {
    // r = b - Ax (on rank-0)
    if (this->root_) {
      Kokkos::deep_copy(R_view, B_view);
      KokkosSparse::spmv("N", mone, crsmat, X_view, one, R_view);
      Kokkos::fence();
    
      if (verbose) {
        // compute residual norm
        std::cout << " > " << numIters << " : norm(r,x,e) = ";
        for (size_t j = 0; j < nrhs; j++) { 
          auto r_subview = Kokkos::subview(R_view, Kokkos::ALL(), j);
          auto x_subview = Kokkos::subview(X_view, Kokkos::ALL(), j);
          host_vector_t r_1d (const_cast<impl_scalar_type*>(r_subview.data()), r_subview.extent(0));
          host_vector_t x_1d (const_cast<impl_scalar_type*>(x_subview.data()), x_subview.extent(0));
          impl_scalar_type rnorm = KokkosBlas::nrm2(r_1d);
          impl_scalar_type xnorm = KokkosBlas::nrm2(x_1d);
          std::cout << rnorm << " -> " << rnorm/bnorms(j) << " " << xnorm << " " << enorms(j) << ", ";
        }
        std::cout << std::endl;
      }
    }

    // e = A^{-1} r 
    Util::put_1d_data_helper_kokkos_view<MVAdapter, host_mvector_t>::
      do_put(Rptr, R_view, ldr, CONTIGUOUS_AND_ROOTED, rowIndexBase);
    static_cast<const solver_type*>(this)->solve_impl(Teuchos::outArg(*E), Teuchos::ptrInArg(*R));
    Util::get_1d_copy_helper_kokkos_view<MVAdapter, host_mvector_t>::
      do_get(initialize_data, Eptr, E_view, lde, CONTIGUOUS_AND_ROOTED, rowIndexBase);
    
    // x = x + e (on rank-0)
    if (this->root_) {
      KokkosBlas::axpy(one, E_view, X_view);

      if (numIters < maxNumIters-1) {
        // compute norm of corrections for "convergence" check
        converged = 1;
        for (size_t j = 0; j < nrhs; j++) { 
          auto e_subview = Kokkos::subview(E_view, Kokkos::ALL(), j);
          host_vector_t e_1d (const_cast<impl_scalar_type*>(e_subview.data()), e_subview.extent(0));
          enorms(j) = KokkosBlas::nrm2(e_1d);
          if (enorms(j) > eps * x0norms(j)) {
            converged = 0;
          }
        }
        if (verbose && converged) {
          std::cout << " converged " << std::endl;
        }
      }
    }

    // broadcast "converged"
    Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &converged);
  } // end of for-loop for IR iteration

  if (verbose && this->root_) {
    // r = b - Ax
    Kokkos::deep_copy(R_view, B_view);
    KokkosSparse::spmv("N", mone, crsmat, X_view, one, R_view);
    Kokkos::fence();
    std::cout << " > final residual norm = ";
    for (size_t j = 0; j < nrhs; j++) { 
      auto r_subview = Kokkos::subview(R_view, Kokkos::ALL(), j);
      host_vector_t r_1d (const_cast<impl_scalar_type*>(r_subview.data()), r_subview.extent(0));
      scalar_type rnorm = KokkosBlas::nrm2(r_1d);
      std::cout << rnorm << " -> " << rnorm/bnorms(j) << ", ";
    }
    std::cout << std::endl << std::endl;
  }

  // put X for output
  Util::put_1d_data_helper_kokkos_view<MVAdapter, host_mvector_t>::
    do_put(Xptr, X_view, ldx, CONTIGUOUS_AND_ROOTED, rowIndexBase);

  return numIters;
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
bool
SolverCore<ConcreteSolver,Matrix,Vector>::matrixShapeOK()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  return( static_cast<solver_type*>(this)->matrixShapeOK_impl() );
}

  // The RCP should probably be to a const Matrix, but that throws a
  // wrench in some of the traits mechanisms that aren't expecting
  // const types
template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::setA( const Teuchos::RCP<const Matrix> a,
                                                EPhase keep_phase )
{
  matrixA_ = createConstMatrixAdapter(a);

#ifdef HAVE_AMESOS2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (keep_phase != CLEAN) &&
                      (globalNumRows_ != matrixA_->getGlobalNumRows() ||
                       globalNumCols_ != matrixA_->getGlobalNumCols()),
                      std::invalid_argument,
                      "Dimensions of new matrix be the same as the old matrix if "
                      "keeping any solver phase" );
#endif

  status_.last_phase_ = keep_phase;

  // Reset phase counters
  switch( status_.last_phase_ ){
  case CLEAN:
    status_.numPreOrder_ = 0;
    // Intentional fallthrough.
  case PREORDERING:
    status_.numSymbolicFact_ = 0;
    // Intentional fallthrough.
  case SYMBFACT:
    status_.numNumericFact_ = 0;
    // Intentional fallthrough.
  case NUMFACT:                 // probably won't ever happen by itself
    status_.numSolve_ = 0;
    // Intentional fallthrough.
  case SOLVE:                   // probably won't ever happen
    break;
  }

  // Re-get the matrix dimensions in case they have changed
  globalNumNonZeros_ = matrixA_->getGlobalNNZ();
  globalNumCols_     = matrixA_->getGlobalNumCols();
  globalNumRows_     = matrixA_->getGlobalNumRows();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::setParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( parameterList->name() == "Amesos2" ){
    // First, validate the parameterList
    Teuchos::RCP<const Teuchos::ParameterList> valid_params = getValidParameters();
    parameterList->validateParameters(*valid_params);

    // Do everything here that is for generic status and control parameters
    control_.setControlParameters(parameterList);

    // Finally, hook to the implementation's parameter list parser
    // First check if there is a dedicated sublist for this solver and use that if there is
    if( parameterList->isSublist(name()) ){
      // Have control look through the solver's parameter list to see if
      // there is anything it recognizes (mostly the "Transpose" parameter)
      control_.setControlParameters(Teuchos::sublist(parameterList, name()));

      static_cast<solver_type*>(this)->setParameters_impl(Teuchos::sublist(parameterList, name()));
    }
  }

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Teuchos::RCP<const Teuchos::ParameterList>
SolverCore<ConcreteSolver,Matrix,Vector>::getValidParameters() const
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1( timers_.totalTime_ );
#endif

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  //RCP<ParameterList> control_params = rcp(new ParameterList("Amesos2 Control"));
  RCP<ParameterList> control_params = rcp(new ParameterList("Amesos2"));
  control_params->set("Transpose", false, "Whether to solve with the matrix transpose");
  control_params->set("Iterative refinement", false, "Whether to solve with iterative refinement");
  control_params->set("Number of iterative refinements", 2, "Number of iterative refinements");
  control_params->set("Verboes for iterative refinement", false, "Verbosity for iterative refinements");
  //  control_params->set("AddToDiag", "");
  //  control_params->set("AddZeroToDiag", false);
  //  control_params->set("MatrixProperty", "general");
  //  control_params->set("Reindex", false);

  RCP<const ParameterList>
    solver_params = static_cast<const solver_type*>(this)->getValidParameters_impl();
  // inject the "Transpose" parameter into the solver's valid parameters
  Teuchos::rcp_const_cast<ParameterList>(solver_params)->set("Transpose", false,
                                                             "Whether to solve with the "
                                                             "matrix transpose");

  RCP<ParameterList> amesos2_params = rcp(new ParameterList("Amesos2"));
  amesos2_params->setParameters(*control_params);
  amesos2_params->set(name(), *solver_params);

  return amesos2_params;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
SolverCore<ConcreteSolver,Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << name() << " solver interface";
  return oss.str();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if( matrixA_.is_null() || (rank_ != 0) ){ return; }
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
  Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm();
  size_t width = 1;
  for( size_t dec = 10; dec < globalNumRows_; dec *= 10 ) {
    ++width;
  }
  width = std::max<size_t>(width,size_t(11)) + 2;
  Teuchos::OSTab tab(out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: print O(P) info, num entries per node
  //    high: print O(N) info, num entries per row
  // extreme: print O(NNZ) info: print indices and values
  //
  // for medium and higher, print constituent objects at specified verbLevel
  if( vl != VERB_NONE ) {
    std::string p = name();
    Util::printLine(out);
    out << this->description() << std::endl << std::endl;

    out << p << "Matrix has " << globalNumRows_ << " rows"
        << " and " << globalNumNonZeros_ << " nonzeros"
        << std::endl;
    if( vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Nonzero elements per row = "
          << globalNumNonZeros_ / globalNumRows_
          << std::endl;
      out << p << "Percentage of nonzero elements = "
          << 100.0 * globalNumNonZeros_ / (globalNumRows_ * globalNumCols_)
          << std::endl;
    }
    if( vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Use transpose = " << control_.useTranspose_
          << std::endl;
      out << p << "Use iterative refinement = " << control_.useIterRefine_
          << std::endl;
    }
    if ( vl == VERB_EXTREME ){
      printTiming(out,vl);
    }
    Util::printLine(out);
  }
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::printTiming(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if( matrixA_.is_null() || (rank_ != 0) ){ return; }

  double preTime  = timers_.preOrderTime_.totalElapsedTime();
  double symTime  = timers_.symFactTime_.totalElapsedTime();
  double numTime  = timers_.numFactTime_.totalElapsedTime();
  double solTime  = timers_.solveTime_.totalElapsedTime();
  double totTime  = timers_.totalTime_.totalElapsedTime();
  double overhead = totTime - (preTime + symTime + numTime + solTime);

  std::string p = name() + " : ";
  Util::printLine(out);

  if(verbLevel != Teuchos::VERB_NONE)
    {
      out << p << "Time to convert matrix to implementation format = "
          << timers_.mtxConvTime_.totalElapsedTime() << " (s)"
          << std::endl;
      out << p << "Time to redistribute matrix = "
          << timers_.mtxRedistTime_.totalElapsedTime() << " (s)"
          << std::endl;

      out << p << "Time to convert vectors to implementation format = "
          << timers_.vecConvTime_.totalElapsedTime() << " (s)"
          << std::endl;
      out << p << "Time to redistribute vectors = "
          << timers_.vecRedistTime_.totalElapsedTime() << " (s)"
          << std::endl;

      out << p << "Number of pre-orderings = "
          << status_.getNumPreOrder()
          << std::endl;
      out << p << "Time for pre-ordering = "
          << preTime << " (s), avg = "
          << preTime / status_.getNumPreOrder() << " (s)"
          << std::endl;

      out << p << "Number of symbolic factorizations = "
          << status_.getNumSymbolicFact()
          << std::endl;
      out << p << "Time for sym fact = "
          << symTime << " (s), avg = "
          << symTime / status_.getNumSymbolicFact() << " (s)"
          << std::endl;

      out << p << "Number of numeric factorizations = "
          << status_.getNumNumericFact()
          << std::endl;
      out << p << "Time for num fact = "
          << numTime << " (s), avg = "
          << numTime / status_.getNumNumericFact() << " (s)"
          << std::endl;

      out << p << "Number of solve phases = "
          << status_.getNumSolve()
          << std::endl;
      out << p << "Time for solve = "
          << solTime << " (s), avg = "
          << solTime / status_.getNumSolve() << " (s)"
          << std::endl;

      out << p << "Total time spent in Amesos2 = "
          << totTime << " (s)"
          << std::endl;
      out << p << "Total time spent in the Amesos2 interface = "
          << overhead << " (s)"
          << std::endl;
      out << p << "  (the above time does not include solver time)"
          << std::endl;
      out << p << "Amesos2 interface time / total time = "
          << overhead / totTime
          << std::endl;
      Util::printLine(out);
    }
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::getTiming(
  Teuchos::ParameterList& timingParameterList) const
{
  Teuchos::ParameterList temp;
  timingParameterList = temp.setName("NULL");
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
SolverCore<ConcreteSolver,Matrix,Vector>::name() const
{
  std::string solverName = solver_type::name;
  return solverName;
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::loadA(EPhase current_phase)
{
  matrix_loaded_ = static_cast<solver_type*>(this)->loadA_impl(current_phase);
}


} // end namespace Amesos2

#endif  // AMESOS2_SOLVERCORE_DEF_HPP
