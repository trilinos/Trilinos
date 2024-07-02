// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_STRUMPACK_DEF_HPP
#define AMESOS2_STRUMPACK_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Util.hpp"

#include <memory>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include "StrumpackSparseSolverMPIDist.hpp"
#else
#include "StrumpackSparseSolver.hpp"
#endif

namespace Amesos2 {


  template <class Matrix, class Vector>
  STRUMPACK<Matrix,Vector>::STRUMPACK(Teuchos::RCP<const Matrix> A,
                                      Teuchos::RCP<Vector> X,
                                      Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::STRUMPACK,Matrix,Vector>(A, X, B)

  {
    using Teuchos::Comm;
#ifdef HAVE_MPI
    using Teuchos::MpiComm;
#endif
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    typedef global_ordinal_type GO;
#ifdef HAVE_MPI
    typedef Tpetra::Map<local_ordinal_type, GO, node_type> map_type;
#endif
    RCP<const Comm<int> > comm = this->getComm ();

#ifdef HAVE_MPI
    RCP<const MpiComm<int> > mpiComm =
      rcp_dynamic_cast<const MpiComm<int> > (comm);
    TEUCHOS_TEST_FOR_EXCEPTION
      (mpiComm.is_null (), std::logic_error, "Amesos2::STRUMPACK "
       "constructor: The matrix's communicator is not an MpiComm!");
    MPI_Comm rawMpiComm = (* (mpiComm->getRawMpiComm ())) ();

    sp_ = Teuchos::RCP<strumpack::StrumpackSparseSolverMPIDist<scalar_type,GO>>
      // (new strumpack::StrumpackSparseSolverMPIDist<scalar_type,GO>(rawMpiComm, this->control_.verbose_));
      (new strumpack::StrumpackSparseSolverMPIDist<scalar_type,GO>(rawMpiComm, true));
#else
    sp_ = Teuchos::RCP<strumpack::StrumpackSparseSolver<scalar_type,GO>>
      (new strumpack::StrumpackSparseSolver<scalar_type,GO>(this->control_.verbose_, this->root_));

#endif

/*
    Do we need this?
    (What parameters do we set here that are not already provided?)
*/
    RCP<ParameterList> default_params =
      parameterList (* (this->getValidParameters ()));
    this->setParameters (default_params);

#ifdef HAVE_MPI
    const size_t myNumRows = this->matrixA_->getLocalNumRows();
    const GO indexBase = this->matrixA_->getRowMap ()->getIndexBase ();
    strumpack_rowmap_ =
      rcp (new map_type (this->globalNumRows_, myNumRows, indexBase, comm));
#endif
  }


////////////////////////////////////////////////////////////////////////////////
//                    DELETE                                                  //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  STRUMPACK<Matrix,Vector>::~STRUMPACK( )
  {
  }


////////////////////////////////////////////////////////////////////////////////
//                     PRE-ORDERING                                           //
////////////////////////////////////////////////////////////////////////////////
  template<class Matrix, class Vector>
  int
  STRUMPACK<Matrix,Vector>::preOrdering_impl()
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTime( this->timers_.preOrderTime_ );
#endif

    // nothing to do: reordering and symbolic factorization are done
    // together in call to ->reorder

    return EXIT_SUCCESS;
  }



////////////////////////////////////////////////////////////////////////////////
//                     SYMBOLIC-FACTORIZATION                                 //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  int
  STRUMPACK<Matrix,Vector>::symbolicFactorization_impl()
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
#endif

    strumpack::ReturnCode ret = sp_->reorder();
    TEUCHOS_TEST_FOR_EXCEPTION( ret == strumpack::ReturnCode::MATRIX_NOT_SET,
                                std::runtime_error,
                                "STRUMPACK matrix reordering failed: "
                                "matrix was not set." );
    TEUCHOS_TEST_FOR_EXCEPTION( ret == strumpack::ReturnCode::REORDERING_ERROR,
                                std::runtime_error,
                                "STRUMPACK matrix reordering failed." );

    return EXIT_SUCCESS;
  }


////////////////////////////////////////////////////////////////////////////////
//                               NUMERIC-FACTORIZATION                        //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  int
  STRUMPACK<Matrix,Vector>::numericFactorization_impl()
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

    strumpack::ReturnCode ret = sp_->factor();
    // Check output
    TEUCHOS_TEST_FOR_EXCEPTION( ret != strumpack::ReturnCode::SUCCESS,
                                std::runtime_error,
                                "Error in STRUMPACK factorization." );

    return EXIT_SUCCESS;
  }


////////////////////////////////////////////////////////////////////////////////
//                                SOLVE                                       //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  int
  STRUMPACK<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                       const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {

#ifdef HAVE_MPI
    // local_len_rhs is how many of the multivector rows belong to
    // this processor
    const size_t local_len_rhs = strumpack_rowmap_->getLocalNumElements();
    const global_size_type nrhs = X->getGlobalNumVectors();

    // make sure our multivector storage is sized appropriately
    bvals_.resize(nrhs * local_len_rhs);
    xvals_.resize(nrhs * local_len_rhs);

    {

#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer(this->timers_.vecConvTime_);
#endif

      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif
        //get dristributed data from Trilinos
       typedef Util::get_1d_copy_helper<MultiVecAdapter<Vector>,scalar_type> copy_helper;
        copy_helper::do_get(B,
                            bvals_(),
                            local_len_rhs,
                            Teuchos::ptrInArg(*strumpack_rowmap_));
      }
    }         // end block for conversion time

    {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
       strumpack::DenseMatrixWrapper<scalar_type>
       Bsp(local_len_rhs, nrhs, bvals_().getRawPtr(), local_len_rhs),
       Xsp(local_len_rhs, nrhs, xvals_().getRawPtr(), local_len_rhs);
       strumpack::ReturnCode ret =sp_->solve(Bsp, Xsp);

       TEUCHOS_TEST_FOR_EXCEPTION( ret != strumpack::ReturnCode::SUCCESS,
                                    std::runtime_error,
                                    "Error in STRUMPACK solve" );
    } // end block for solve time


    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      //get dristributed data from STRUMPACK after solving
      typedef Util::put_1d_data_helper<MultiVecAdapter<Vector>,scalar_type> put_helper;
      put_helper::do_put(X,
                         xvals_(),
                         local_len_rhs,
                         Teuchos::ptrInArg(*strumpack_rowmap_));
    }
#else //NO MPI
    using Teuchos::as;
    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    const size_t nrhs = X->getGlobalNumVectors();
    bvals_.resize(nrhs * ld_rhs);
    xvals_.resize(nrhs * ld_rhs);

    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

      strumpack::DenseMatrixWrapper<scalar_type>
         Bsp(ld_rhs, nrhs, bvals_().getRawPtr(), ld_rhs),
         Xsp(ld_rhs, nrhs, xvals_().getRawPtr(), ld_rhs);
      strumpack::ReturnCode ret =sp_->solve(Bsp, Xsp);

      TEUCHOS_TEST_FOR_EXCEPTION( ret != strumpack::ReturnCode::SUCCESS,
                                    std::runtime_error,
                                    "Error in STRUMPACK solve" );
    } // end block for solve time

    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    Util::put_1d_data_helper<
        MultiVecAdapter<Vector>,scalar_type>::do_put(X, xvals_(),
            as<size_t>(ld_rhs),
            ROOTED, this->rowIndexBase_);
  }
#endif
    return EXIT_SUCCESS;
  }


  template <class Matrix, class Vector>
  bool
  STRUMPACK<Matrix,Vector>::matrixShapeOK_impl() const
  {
#ifdef HAVE_MPI
    // STRUMPACK requires square matrices
    return( this->globalNumRows_ == this->globalNumCols_ );
#else
    return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
#endif
  }



////////////////////////////////////////////////////////////////////////////////
//                              SET_PARAMETERS                                //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  void
  STRUMPACK<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterEntryValidator;

    RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

    if( parameterList->isParameter("Matching") ){
       RCP<const ParameterEntryValidator> matching_validator = valid_params->getEntry("Matching").validator();
       parameterList->getEntry("Matching").setValidator(matching_validator);

       sp_->options().set_matching(getIntegralValue<strumpack::MatchingJob>(*parameterList, "Matching"));
     }

    if( parameterList->isParameter("Ordering") ){
       RCP<const ParameterEntryValidator> reordering_validator = valid_params->getEntry("Ordering").validator();
       parameterList->getEntry("Ordering").setValidator(reordering_validator);

       sp_->options().set_reordering_method(getIntegralValue<strumpack::ReorderingStrategy>(*parameterList, "Ordering"));
     }

    if( parameterList->isParameter("ReplaceTinyPivot") ){
       RCP<const ParameterEntryValidator> replacepivot_validator = valid_params->getEntry("ReplaceTinyPivot").validator();
       parameterList->getEntry("ReplaceTinyPivot").setValidator(replacepivot_validator);

       if( replacepivot_validator) {
         sp_->options().enable_replace_tiny_pivots();
       }
       else{
         sp_->options().disable_replace_tiny_pivots();
       }
    }

    if( parameterList->isParameter("IterRefine") ){
       RCP<const ParameterEntryValidator> iter_refine_validator = valid_params->getEntry("IterRefine").validator();
       parameterList->getEntry("IterRefine").setValidator(iter_refine_validator);

       sp_->options().set_Krylov_solver(getIntegralValue<strumpack::KrylovSolver>(*parameterList, "IterRefine"));
    }

    if( parameterList->isParameter("Compression") ){
       RCP<const ParameterEntryValidator> compression_validator = valid_params->getEntry("Compression").validator();
       parameterList->getEntry("Compression").setValidator(compression_validator);

       sp_->options().set_compression(getIntegralValue<strumpack::CompressionType>(*parameterList, "Compression"));
    }

    TEUCHOS_TEST_FOR_EXCEPTION( this->control_.useTranspose_,
                        std::invalid_argument,
                        "STRUMPACK does not support solving the tranpose system" );

  }

  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  STRUMPACK<Matrix,Vector>::getValidParameters_impl() const
  {
    using std::string;
    using Teuchos::tuple;
    using Teuchos::ParameterList;
    using Teuchos::EnhancedNumberValidator;
    using Teuchos::setStringToIntegralParameter;
    using Teuchos::stringToIntegralParameterEntryValidator;

    static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

    if( is_null(valid_params) ){
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

        setStringToIntegralParameter<strumpack::MatchingJob>("Matching", "NONE",
                                                    "Specifies how to permute the "
                                                    "matrix for numerical stability",
                                                    tuple<string>("NONE", "MAX_CARDINALITY", "MAX_SMALLEST_DIAGONAL", "MAX_SMALLEST_DIAGONAL_2", "MAX_DIAGONAL_SUM", "MAX_DIAGONAL_PRODUCT_SCALING", "COMBBLAS"),
                                                    tuple<string>("NONE", "MAX_CARDINALITY", "MAX_SMALLEST_DIAGONAL", "MAX_SMALLEST_DIAGONAL_2", "MAX_DIAGONAL_SUM", "MAX_DIAGONAL_PRODUCT_SCALING", "COMBBLAS"),
                                                    tuple<strumpack::MatchingJob>(strumpack::MatchingJob::NONE,
                                                                           strumpack::MatchingJob::MAX_CARDINALITY,
                                                                           strumpack::MatchingJob::MAX_SMALLEST_DIAGONAL,
                                                                           strumpack::MatchingJob::MAX_SMALLEST_DIAGONAL_2,
                                                                           strumpack::MatchingJob::MAX_DIAGONAL_SUM,
                                                                           strumpack::MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING, 
                                                                           strumpack::MatchingJob::COMBBLAS),
                                                    pl.getRawPtr());

#if defined(STRUMPACK_USE_PARMETIS)
        std::string default_ordering("PARMETIS");
#else
        std::string default_ordering("METIS");
#endif
        setStringToIntegralParameter<strumpack::ReorderingStrategy>("Ordering", default_ordering,
                                                    "Specifies how to permute the "
                                                    "matrix for sparsity preservation",
                                                    tuple<string>("NATURAL", "PARMETIS", "METIS", "SCOTCH", "GEOMETRIC", "PTSCOTCH", "RCM"),
                                                    tuple<string>("Natural ordering",
                                                                  "ParMETIS ordering",
                                                                  "Metis ordering",
                                                                  "Scotch ordering",
                                                                  "Geometric ordering",
                                                                  "PTScotch ordering",
                                                                  "RCM"),
                                                    tuple<strumpack::ReorderingStrategy>(strumpack::ReorderingStrategy::NATURAL,
                                                                           strumpack::ReorderingStrategy::PARMETIS,
                                                                           strumpack::ReorderingStrategy::METIS,
                                                                           strumpack::ReorderingStrategy::SCOTCH,
                                                                           strumpack::ReorderingStrategy::GEOMETRIC,
                                                                           strumpack::ReorderingStrategy::PTSCOTCH,
                                                                           strumpack::ReorderingStrategy::RCM),
                                                    pl.getRawPtr());

        pl->set("ReplaceTinyPivot", true, "Specifies whether to replace tiny diagonals during LU factorization");


//   There are multiple options available for an iterative refinement,
//   however we recommend the use of "DIRECT" within the Amesos2 interface
        setStringToIntegralParameter<strumpack::KrylovSolver>("IterRefine", "DIRECT",
                                                    "Type of iterative refinement to use",
                                                    tuple<string>("AUTO", "DIRECT", "REFINE", "PREC_GMRES", "GMRES", "PREC_BICGSTAB", "BICGSTAB"),
                                                    tuple<string>("Use iterative refinement if no compression is used, otherwise use GMRes.",
                                                                  "Single application of the multifrontal solver.",
                                                                  "Iterative refinement.",
                                                                  "Preconditioned GMRes.",
                                                                  "UN-preconditioned GMRes.",
                                                                  "Preconditioned BiCGStab.",
                                                                  "UN-preconditioned BiCGStab."),
                                                    tuple<strumpack::KrylovSolver>(strumpack::KrylovSolver::AUTO,
                                                                           strumpack::KrylovSolver::DIRECT,
                                                                           strumpack::KrylovSolver::REFINE,
                                                                           strumpack::KrylovSolver::PREC_GMRES,
                                                                           strumpack::KrylovSolver::GMRES,
                                                                           strumpack::KrylovSolver::PREC_BICGSTAB,
                                                                           strumpack::KrylovSolver::BICGSTAB),
                                                    pl.getRawPtr());

//   There are multiple options available for the compression of the matrix,
//   we recommend the use of "NONE" within the Amesos2 interface
        setStringToIntegralParameter<strumpack::CompressionType>("Compression", "NONE",
                                                    "Type of compression to use",
                                                    tuple<string>("NONE", "HSS", "BLR", "HODLR", "LOSSLESS", "LOSSY"),
                                                    tuple<string>("No compression, purely direct solver.",
                                                                  "HSS compression of frontal matrices.",
                                                                  "Block low-rank compression of fronts.",
                                                                  "Hierarchically Off-diagonal Low-Rank compression of frontal matrices.",
                                                                  "Lossless compresssion.",
                                                                  "Lossy compresssion."),
                                                    tuple<strumpack::CompressionType>(strumpack::CompressionType::NONE,
                                                                           strumpack::CompressionType::HSS,
                                                                           strumpack::CompressionType::BLR,
                                                                           strumpack::CompressionType::HODLR,
                                                                           strumpack::CompressionType::LOSSLESS,
                                                                           strumpack::CompressionType::LOSSY),
                                                    pl.getRawPtr());




      valid_params = pl;
    }

    return valid_params;
  }



////////////////////////////////////////////////////////////////////////////////
//                           LOAD_DATA                                        //
////////////////////////////////////////////////////////////////////////////////

  template <class Matrix, class Vector>
  bool
  STRUMPACK<Matrix,Vector>::loadA_impl(EPhase current_phase){
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ptrInArg;
    using Teuchos::as;
    using Teuchos::rcp_dynamic_cast; // Do I need this?

    using Teuchos::Comm;
#ifdef HAVE_MPI
    using Teuchos::MpiComm;

    #ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
    #endif

    Teuchos::RCP<const MatrixAdapter<Matrix> > redist_mat
      = this->matrixA_->get(ptrInArg(*strumpack_rowmap_));

    typedef global_ordinal_type GO;
    GO l_nnz, l_rows;
    l_nnz  = as<GO>(redist_mat->getLocalNNZ());
    l_rows = as<GO>(redist_mat->getLocalNumRows());

    RCP<const Comm<int> > comm = this->getComm ();
    RCP<const MpiComm<int> > mpiComm =
      rcp_dynamic_cast<const MpiComm<int> > (comm);

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    Array<GO> dist(numProcs+1);
    dist[0] = 0;
    dist[myRank+1] = as<GO>(strumpack_rowmap_->getMaxGlobalIndex()) + 1;

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dist.data()+1, sizeof(GO), MPI_BYTE,
                   mpiComm->getRawMpiComm()->operator()());
    Kokkos::resize(nzvals_view_, l_nnz);
    Kokkos::resize(colind_view_, l_nnz);
    Kokkos::resize(rowptr_view_, l_rows + 1);


    GO nnz_ret = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

      Util::get_crs_helper_kokkos_view<MatrixAdapter<Matrix>,
        host_value_type_array, host_ordinal_type_array, host_ordinal_type_array >::do_get(
                                         redist_mat.ptr(),
                                         nzvals_view_, colind_view_, rowptr_view_,
                                         nnz_ret,
                                         ptrInArg(*strumpack_rowmap_),
                                         ROOTED,
                                         ARBITRARY);
    }


    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != l_nnz,
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");

    // Get the csr data type for this type of matrix
    sp_->set_distributed_csr_matrix
      (l_rows, rowptr_view_.data(), colind_view_.data(),
       nzvals_view_.data(), dist.getRawPtr(), false);

#else
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    typedef global_ordinal_type GO;
    GO nnz_ret = 0;

    if( this->root_ ){
      Kokkos::resize(nzvals_view_, this->globalNumNonZeros_);
      Kokkos::resize(colind_view_, this->globalNumNonZeros_);
      Kokkos::resize(rowptr_view_, this->globalNumRows_ + 1);
    }
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    Util::get_crs_helper_kokkos_view<MatrixAdapter<Matrix>,
      host_value_type_array, host_ordinal_type_array, host_ordinal_type_array >::do_get(
                                       this->matrixA_.ptr(),
                                       nzvals_view_, colind_view_, rowptr_view_,
                                       nnz_ret,
                                       ROOTED,
                                       ARBITRARY, this->rowIndexBase_);
   }

    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != this->globalNumNonZeros_,
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");

    // Get the csr data type for this type of matrix
    sp_->set_csr_matrix(this->globalNumRows_, rowptr_view_.data(), colind_view_.data(),
       nzvals_view_.data(), false);

#endif
    return true;
  }


  template<class Matrix, class Vector>
  const char* STRUMPACK<Matrix,Vector>::name = "STRUMPACK";


} // end namespace Amesos2
#endif  // AMESOS2_STRUMPACK_DEF_HPP
