// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_TACHO_DEF_HPP
#define AMESOS2_TACHO_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Tacho_decl.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

template <class Matrix, class Vector>
TachoSolver<Matrix,Vector>::TachoSolver(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::TachoSolver,Matrix,Vector>(A, X, B)
{
  data_.method    = 1;      // Cholesky
  data_.variant   = 2;      // solver variant
  data_.streams   = 1;      // # of streams
  data_.dofs_per_node = 1;  // DoFs / node
  data_.pivot_pert = false; // Diagonal pertubation
  data_.verbose    = false; // verbose
}


template <class Matrix, class Vector>
TachoSolver<Matrix,Vector>::~TachoSolver( )
{
  if ( this->root_ ) {
    data_.solver.release();
  }
}

template <class Matrix, class Vector>
std::string
TachoSolver<Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << "Tacho solver interface";
  return oss.str();
}

template<class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::preOrdering_impl()
{
  return(0);
}

template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::symbolicFactorization_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
#endif

  int status = 0;
  if ( this->root_ ) {
    if(do_optimization()) {
      this->matrixA_->returnRowPtr_kokkos_view(host_row_ptr_view_);
      this->matrixA_->returnColInd_kokkos_view(host_cols_view_);
    }

    data_.solver.setSolutionMethod(data_.method);
    data_.solver.setLevelSetOptionAlgorithmVariant(data_.variant);
    data_.solver.setSmallProblemThresholdsize(data_.small_problem_threshold_size);
    data_.solver.setVerbose(data_.verbose);
    data_.solver.setLevelSetOptionNumStreams(data_.streams);
    // TODO: Confirm param options
    // data_.solver.setMaxNumberOfSuperblocks(data_.max_num_superblocks);

    // Symbolic factorization currently must be done on host
    if (data_.dofs_per_node > 1) {
      data_.solver.analyze(this->globalNumCols_, data_.dofs_per_node, host_row_ptr_view_, host_cols_view_);
    } else {
      data_.solver.analyze(this->globalNumCols_, host_row_ptr_view_, host_cols_view_);
    }
    data_.solver.initialize();
  }
  return status;
}


template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::numericFactorization_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

  int status = 0;
  if ( this->root_ ) {
    if(do_optimization()) {
     // instead of holding onto the device poinster
     //  this->matrixA_->returnValues_kokkos_view(device_nzvals_view_);
     // make an explicit copy
     device_value_type_array device_nzvals_temp;
     this->matrixA_->returnValues_kokkos_view(device_nzvals_temp);
     Kokkos::deep_copy(device_nzvals_view_, device_nzvals_temp);
    }
    if (data_.pivot_pert) {
      data_.solver.useDefaultPivotTolerance();
    } else {
      data_.solver.useNoPivotTolerance();
    }
    data_.solver.factorize(device_nzvals_view_);
  }
  return status;
}

template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  // don't allocate b since it's handled by the copy manager and might just be
  // be assigned, not copied anyways.
  // also don't allocate x since we will also use do_get to allocate this if
  // necessary. When a copy is not necessary we'll solve directly to the x
  // values in the MV.
  bool bDidAssignX;
  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
#endif
    const bool initialize_data = true;
    const bool do_not_initialize_data = false;
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
                             device_solve_array_t>::do_get(initialize_data, B, this->bValues_,
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
    bDidAssignX = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
                             device_solve_array_t>::do_get(do_not_initialize_data, X, this->xValues_,
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
  }

  int ierr = 0; // returned error code

  if ( this->root_ ) {  // Do solve!
    // Bump up the workspace size if needed
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
    if (workspace_.extent(0) < this->globalNumRows_ || workspace_.extent(1) < nrhs) {
      workspace_ = device_solve_array_t(
        Kokkos::ViewAllocateWithoutInitializing("t"), this->globalNumRows_, nrhs);
    }

    data_.solver.solve(xValues_, bValues_, workspace_);

    int status = 0; // TODO: determine what error handling will be
    if(status != 0) {
      ierr = status;
    }
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  TEUCHOS_TEST_FOR_EXCEPTION( ierr != 0, std::runtime_error,
    "tacho_solve has error code: " << ierr );

  /* Update X's global values */

  // if bDidAssignX, then we solved straight to the adapter's X memory space without
  // requiring additional memory allocation, so the x data is already in place.
  if(!bDidAssignX) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    // This will do nothing is if the target view matches the src view, which
    // can be the case if the memory spaces match. See comments above for do_get.
    Util::template put_1d_data_helper_kokkos_view<
      MultiVecAdapter<Vector>,device_solve_array_t>::do_put(X, xValues_,
                                        as<size_t>(ld_rhs),
                                        ROOTED, this->rowIndexBase_);
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::matrixShapeOK_impl() const
{
  // Tacho can only apply the solve routines to square matrices
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
TachoSolver<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  // TODO: Confirm param options

  // factorization type
  auto method_name = parameterList->get<std::string> ("method", "chol");
  if (method_name == "chol")
    data_.method = 1;
  else if (method_name == "ldl")
    data_.method = 2;
  else if (method_name == "lu")
    data_.method = 3;
  else {
    std::cout << "Error: not supported solution method\n";
  }
  // solver type
  data_.variant = parameterList->get<int> ("variant", 2);
  // small problem threshold
  data_.small_problem_threshold_size = parameterList->get<int> ("small problem threshold size", 1024);
  // verbosity
  data_.verbose = parameterList->get<bool> ("verbose", false);
  // # of streams
  data_.streams = parameterList->get<int> ("num-streams", 1);
  // DoFs / node
  data_.dofs_per_node = parameterList->get<int> ("dofs-per-node", 1);
  // Perturb tiny pivots
  data_.pivot_pert = parameterList->get<bool> ("perturb-pivot", false);
  // TODO: Confirm param options
  // data_.num_kokkos_threads = parameterList->get<int>("kokkos-threads", 1);
  // data_.max_num_superblocks = parameterList->get<int>("max-num-superblocks", 4);
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
TachoSolver<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set("method", "chol", "Type of factorization, chol, ldl, or lu");
    pl->set("variant", 2, "Type of solver variant, 0, 1, or 2");
    pl->set("small problem threshold size", 1024, "Problem size threshold below with Tacho uses LAPACK.");
    pl->set("verbose", false, "Verbosity");
    pl->set("num-streams", 1, "Number of GPU streams");
    pl->set("dofs-per-node", 1, "DoFs per node");
    pl->set("perturb-pivot", false, "Perturb tiny pivots");

    // TODO: Confirm param options
    // pl->set("kokkos-threads", 1, "Number of threads");
    // pl->set("max-num-superblocks", 4, "Max number of superblocks");

    valid_params = pl;
  }

  return valid_params;
}

template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::do_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1));
}

template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::loadA_impl(EPhase current_phase)
{

  if(current_phase == SOLVE) {
    return(false);
  }

  if(!do_optimization()) {
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Note views are allocated but eventually we should remove this.
    // The internal copy manager will decide if we can assign or deep_copy
    // and then allocate if necessary. However the GPU solvers are serial right
    // now so I didn't complete refactoring the matrix code for the parallel
    // case. If we added that later, we should have it hooked up to the copy
    // manager and then these allocations can go away.
    {
      if( this->root_ ) {
        if (device_nzvals_view_.extent(0) != this->globalNumNonZeros_)
          Kokkos::resize(device_nzvals_view_, this->globalNumNonZeros_);
        if (host_cols_view_.extent(0) != this->globalNumNonZeros_)
          Kokkos::resize(host_cols_view_, this->globalNumNonZeros_);
        if (host_row_ptr_view_.extent(0) != this->globalNumRows_ + 1)
          Kokkos::resize(host_row_ptr_view_, this->globalNumRows_ + 1);
      } else {
        Kokkos::resize(device_nzvals_view_, 0);
        Kokkos::resize(host_cols_view_, 0);
        Kokkos::resize(host_row_ptr_view_, 1);
      }
    }

    typename host_size_type_array::value_type nnz_ret = 0;
    {
  #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
  #endif

      TEUCHOS_TEST_FOR_EXCEPTION( this->rowIndexBase_ != this->columnIndexBase_,
                          std::runtime_error,
                          "Row and column maps have different indexbase ");

      Util::get_crs_helper_kokkos_view<MatrixAdapter<Matrix>,
        device_value_type_array, host_ordinal_type_array, host_size_type_array>::do_get(
                                                      this->matrixA_.ptr(),
                                                      device_nzvals_view_,
                                                      host_cols_view_,
                                                      host_row_ptr_view_,
                                                      nnz_ret,
                                                      ROOTED, ARBITRARY,
                                                      this->columnIndexBase_);
    }
  }
  else {
    if( this->root_ ) {
      // instead of holding onto the device poinster (which could cause issue)
      // make an explicit copy
      device_nzvals_view_ = device_value_type_array(
        Kokkos::ViewAllocateWithoutInitializing("nzvals"), this->globalNumNonZeros_);
    }
  }

  return true;
}


template<class Matrix, class Vector>
const char* TachoSolver<Matrix,Vector>::name = "Tacho";


} // end namespace Amesos2

#endif  // AMESOS2_TACHO_DEF_HPP
