// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
// @HEADER

/**
   \file   Amesos2_Superludist_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Tue May 24 08:37:17 MDT 2011

   \brief  Definitions for the Amesos2 SuperLU_DIST solver interface
*/

#ifndef AMESOS2_SUPERLUDIST_DEF_HPP
#define AMESOS2_SUPERLUDIST_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Superludist_TypeMap.hpp"
#include "Amesos2_Util.hpp"


namespace Amesos2 {


  template <class Matrix, class Vector>
  Superludist<Matrix,Vector>::Superludist(Teuchos::RCP<const Matrix> A,
                                          Teuchos::RCP<Vector> X,
                                          Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::Superludist,Matrix,Vector>(A, X, B)
    , nzvals_()                 // initialization to empty arrays
    , colind_()
    , rowptr_()
    , bvals_()
    , xvals_()
    , in_grid_(false)
  {
    ////////////////////////////////////////////
    // Set up the SuperLU_DIST processor grid //
    ////////////////////////////////////////////

    int nprocs = this->getComm()->getSize();
    SLUD::int_t nprow, npcol;
    get_default_grid_size(nprocs, nprow, npcol);
    data_.mat_comm = dynamic_cast<const Teuchos::MpiComm<int>* >(this->matrixA_->getComm().getRawPtr())->getRawMpiComm()->operator()();
    SLUD::superlu_gridinit(data_.mat_comm, nprow, npcol, &(data_.grid));

    ////////////////////////////////////////////////////////
    // Set Some default parameters.                       //
    //                                                    //
    // Must do this after grid has been created in        //
    // case user specifies the nprow and npcol parameters //
    ////////////////////////////////////////////////////////
    Teuchos::RCP<Teuchos::ParameterList> default_params
      = Teuchos::parameterList( *(this->getValidParameters()) );
    this->setParameters(default_params);

    // Set some internal options
    data_.options.Fact = SLUD::DOFACT;
    data_.equed = SLUD::NOEQUIL; // No equilibration has yet been performed
    data_.options.SolveInitialized  = SLUD::NO;
    data_.options.RefineInitialized = SLUD::NO;
    data_.rowequ = false;
    data_.colequ = false;
    data_.perm_r.resize(this->globalNumRows_);
    data_.perm_c.resize(this->globalNumCols_);

    ////////////////////////////////////////////////////////////////
    // Set up a communicator for the parallel column ordering and //
    // parallel symbolic factorization.                           //
    ////////////////////////////////////////////////////////////////
    data_.symb_comm = MPI_COMM_NULL;
    int color = MPI_UNDEFINED;
    int my_rank = this->rank_;

    /* domains is the next power of 2 less than nprow*npcol.  This
     * value will be used for creating an MPI communicator for the
     * pre-ordering and symbolic factorization methods.
     */
    data_.domains = (int) ( pow(2.0, floor(log10((double)nprow*npcol)/log10(2.0))) );

    if( this->rank_ < data_.domains ) color = 0;
    MPI_Comm_split (data_.mat_comm, color, my_rank, &(data_.symb_comm));

    //////////////////////////////////////////////////////////////////////
    // Set up a row map that maps to only processors that are in the    //
    // SuperLU processor grid.  This will be used for redistributing A. //
    //////////////////////////////////////////////////////////////////////

    int my_weight = 0;
    if( this->rank_ < nprow * npcol ){
      in_grid_ = true; my_weight = 1; // I am in the grid, and I get some of the matrix rows
    }
    // TODO: might only need to initialize if parallel symbolic factorization is requested.
    // TODO: Need to fix this map for indexbase ?
    superlu_rowmap_
      = Tpetra::createWeightedContigMapWithNode<local_ordinal_type,
      global_ordinal_type,
      node_type>(my_weight,
                 this->globalNumRows_,
                 this->getComm(),
                 Kokkos::DefaultNode::getDefaultNode());
    // TODO: the node above should technically come from the matrix
    // itself.  Might need to add a getNode method to the matrix
    // adapter.

    //////////////////////////////////
    // Do some other initialization //
    //////////////////////////////////

    data_.A.Store = NULL;
    function_map::LUstructInit(this->globalNumRows_, this->globalNumCols_, &(data_.lu));
    SLUD::PStatInit(&(data_.stat));
    // We do not use ScalePermstructInit because we will use our own
    // arrays for storing perm_r and perm_c
    data_.scale_perm.perm_r = data_.perm_r.getRawPtr();
    data_.scale_perm.perm_c = data_.perm_c.getRawPtr();
  }


  template <class Matrix, class Vector>
  Superludist<Matrix,Vector>::~Superludist( )
  {
    /* Free SuperLU_DIST data_types
     * - Matrices
     * - Vectors
     * - Stat object
     * - ScalePerm, LUstruct, grid, and solve objects
     *
     * Note: the function definitions are the same regardless whether
     * complex or real, so we arbitrarily use the D namespace
     */
    if ( this->status_.getNumPreOrder() > 0 ){
      free( data_.sizes );
      free( data_.fstVtxSep );
    }

    // Cleanup old matrix store memory if it's non-NULL.  Our
    // Teuchos::Array's will destroy rowind, colptr, and nzval for us
    if( data_.A.Store != NULL ){
      SLUD::Destroy_SuperMatrix_Store_dist( &(data_.A) );
    }

    // LU data is initialized in numericFactorization_impl()
    if ( this->status_.getNumNumericFact() > 0 ){
      function_map::Destroy_LU(this->globalNumRows_, &(data_.grid), &(data_.lu));
    }
    function_map::LUstructFree(&(data_.lu));

    // If a symbolic factorization is ever performed without a
    // follow-up numericfactorization, there are some arrays in the
    // Pslu_freeable struct which will never be free'd by
    // SuperLU_DIST.
    if ( this->status_.symbolicFactorizationDone() &&
         !this->status_.numericFactorizationDone() ){
      if ( data_.pslu_freeable.xlsub != NULL ){
        free( data_.pslu_freeable.xlsub );
        free( data_.pslu_freeable.lsub );
      }
      if ( data_.pslu_freeable.xusub != NULL ){
        free( data_.pslu_freeable.xusub );
        free( data_.pslu_freeable.usub );
      }
      if ( data_.pslu_freeable.supno_loc != NULL ){
        free( data_.pslu_freeable.supno_loc );
        free( data_.pslu_freeable.xsup_beg_loc );
        free( data_.pslu_freeable.xsup_end_loc );
      }
      free( data_.pslu_freeable.globToLoc );
    }

    SLUD::PStatFree( &(data_.stat) ) ;

    // Teuchos::Arrays will free R, C, perm_r, and perm_c
    // SLUD::D::ScalePermstructFree(&(data_.scale_perm));

    if ( data_.options.SolveInitialized == SLUD::YES )
      function_map::SolveFinalize(&(data_.options), &(data_.solve_struct));

    SLUD::superlu_gridexit(&(data_.grid)); // TODO: are there any
                                           // cases where grid
                                           // wouldn't be initialized?

    if ( data_.symb_comm != MPI_COMM_NULL ) MPI_Comm_free(&(data_.symb_comm));
  }

  template<class Matrix, class Vector>
  int
  Superludist<Matrix,Vector>::preOrdering_impl()
  {
    // We will always use the NATURAL row ordering to avoid the
    // sequential bottleneck present when doing any other row
    // ordering scheme from SuperLU_DIST
    //
    // Set perm_r to be the natural ordering
    SLUD::int_t slu_rows_ub = Teuchos::as<SLUD::int_t>(this->globalNumRows_);
    for( SLUD::int_t i = 0; i < slu_rows_ub; ++i ) data_.perm_r[i] = i;

    // loadA_impl();                    // Refresh matrix values

    if( in_grid_ ){
      // If this function has been called at least once, then the
      // sizes, and fstVtxSep arrays were allocated in
      // get_perm_c_parmetis.  Delete them before calling that
      // function again.  These arrays will also be dealloc'd in the
      // deconstructor.
      if( this->status_.getNumPreOrder() > 0 ){
        free( data_.sizes );
        free( data_.fstVtxSep );
      }
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor preOrderTime( this->timers_.preOrderTime_ );
#endif

      float info = 0.0;
      info = SLUD::get_perm_c_parmetis( &(data_.A),
                                        data_.perm_r.getRawPtr(), data_.perm_c.getRawPtr(),
                                        data_.grid.nprow * data_.grid.npcol, data_.domains,
                                        &(data_.sizes), &(data_.fstVtxSep),
                                        &(data_.grid), &(data_.symb_comm) );

      TEUCHOS_TEST_FOR_EXCEPTION( info > 0.0,
                          std::runtime_error,
                          "SuperLU_DIST pre-ordering ran out of memory after allocating "
                          << info << " bytes of memory" );
    }

    // Ordering will be applied directly before numeric factorization,
    // after we have a chance to get updated coefficients from the
    // matrix

    return EXIT_SUCCESS;
  }



  template <class Matrix, class Vector>
  int
  Superludist<Matrix,Vector>::symbolicFactorization_impl()
  {
    // loadA_impl();                    // Refresh matrix values

    if( in_grid_ ){

#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
#endif

      float info = 0.0;
      info = SLUD::symbfact_dist((data_.grid.nprow) * (data_.grid.npcol),
                                 data_.domains, &(data_.A), data_.perm_c.getRawPtr(),
                                 data_.perm_r.getRawPtr(), data_.sizes,
                                 data_.fstVtxSep, &(data_.pslu_freeable),
                                 &(data_.grid.comm), &(data_.symb_comm),
                                 &(data_.mem_usage));

      TEUCHOS_TEST_FOR_EXCEPTION( info > 0.0,
                          std::runtime_error,
                          "SuperLU_DIST symbolic factorization ran out of memory after"
                          " allocating " << info << " bytes of memory" );
    }
    same_symbolic_ = false;
    same_solve_struct_ = false;

    return EXIT_SUCCESS;
  }


  template <class Matrix, class Vector>
  int
  Superludist<Matrix,Vector>::numericFactorization_impl(){
    using Teuchos::as;

    // loadA_impl();                    // Refresh the matrix values

    // if( data_.options.Equil == SLUD::YES ){
    //   // Apply the scalings computed in preOrdering
    //   function_map::laqgs(&(data_.A), data_.R.getRawPtr(),
    //                    data_.C.getRawPtr(), data_.rowcnd, data_.colcnd,
    //                    data_.amax, &(data_.equed));

    //   data_.rowequ = (data_.equed == SLUD::ROW) || (data_.equed == SLUD::BOTH);
    //   data_.colequ = (data_.equed == SLUD::COL) || (data_.equed == SLUD::BOTH);
    // }

    if( in_grid_ ){
      // Apply the column ordering, so that AC is the column-permuted A, and compute etree
      size_t nnz_loc = ((SLUD::NRformat_loc*)data_.A.Store)->nnz_loc;
      for( size_t i = 0; i < nnz_loc; ++i ) colind_[i] = data_.perm_c[colind_[i]];

      // Distribute data from the symbolic factorization
      if( same_symbolic_ ){
        // Note: with the SamePattern_SameRowPerm options, it does not
        // matter that the glu_freeable member has never been
        // initialized, because it is never accessed.  It is a
        // placeholder arg.  The real work is done in data_.lu
        function_map::pdistribute(SLUD::SamePattern_SameRowPerm,
                                  as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                  &(data_.A), &(data_.scale_perm),
                                  &(data_.glu_freeable), &(data_.lu),
                                  &(data_.grid));
      } else {
        function_map::dist_psymbtonum(SLUD::DOFACT,
                                      as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                      &(data_.A), &(data_.scale_perm),
                                      &(data_.pslu_freeable), &(data_.lu),
                                      &(data_.grid));
      }

      // Retrieve the normI of A (required by gstrf).
      double anorm = function_map::plangs((char *)"I", &(data_.A), &(data_.grid));

      int info = 0;
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

        function_map::gstrf(&(data_.options), this->globalNumRows_,
                            this->globalNumCols_, anorm, &(data_.lu),
                            &(data_.grid), &(data_.stat), &info);
      }

      // Check output
      TEUCHOS_TEST_FOR_EXCEPTION( info > 0,
                          std::runtime_error,
                          "L and U factors have been computed but U("
                          << info << "," << info << ") is exactly zero "
                          "(i.e. U is singular)");
    }

    // The other option, that info_st < 0, denotes invalid parameters
    // to the function, but we'll assume for now that that won't
    // happen.

    data_.options.Fact = SLUD::FACTORED;
    same_symbolic_ = true;

    return EXIT_SUCCESS;
  }


  template <class Matrix, class Vector>
  int
  Superludist<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                         const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    using Teuchos::as;

    // local_len_rhs is how many of the multivector rows belong to
    // this processor in the SuperLU_DIST processor grid.
    const size_t local_len_rhs = superlu_rowmap_->getNodeNumElements();
    const global_size_type nrhs = X->getGlobalNumVectors();
    const global_ordinal_type first_global_row_b = superlu_rowmap_->getMinGlobalIndex();

    // make sure our multivector storage is sized appropriately
    bvals_.resize(nrhs * local_len_rhs);
    xvals_.resize(nrhs * local_len_rhs);

    // We assume the global length of the two vectors have already been
    // checked for compatibility

    {                           // get the values from B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer(this->timers_.vecConvTime_);
#endif

      {
        // The input dense matrix for B should be distributed in the
        // same manner as the superlu_dist matrix.  That is, if a
        // processor has m_loc rows of A, then it should also have
        // m_loc rows of B (and the same rows).  We accomplish this by
        // distributing the multivector rows with the same Map that
        // the matrix A's rows are distributed.
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

        // get grid-distributed mv data.  The multivector data will be
        // distributed across the processes in the SuperLU_DIST grid.
        typedef Util::get_1d_copy_helper<MultiVecAdapter<Vector>,slu_type> copy_helper;
        copy_helper::do_get(B,
                            bvals_(),
                            local_len_rhs,
                            Teuchos::ptrInArg(*superlu_rowmap_));
      }
    }         // end block for conversion time

    if( in_grid_ ){
      // if( data_.options.trans == SLUD::NOTRANS ){
      //   if( data_.rowequ ){            // row equilibration has been done on AC
      //  // scale bxvals_ by diag(R)
      //  Util::scale(bxvals_(), as<size_t>(len_rhs), ldbx_, data_.R(),
      //              SLUD::slu_mt_mult<slu_type,magnitude_type>());
      //   }
      // } else if( data_.colequ ){       // column equilibration has been done on AC
      //   // scale bxvals_ by diag(C)
      //   Util::scale(bxvals_(), as<size_t>(len_rhs), ldbx_, data_.C(),
      //            SLUD::slu_mt_mult<slu_type,magnitude_type>());
      // }

      // Initialize the SOLVEstruct_t.
      //
      // We are able to reuse the solve struct if we have not changed
      // the sparsity pattern of L and U since the last solve
      if( !same_solve_struct_ ){
        if( data_.options.SolveInitialized == SLUD::YES ){
          function_map::SolveFinalize(&(data_.options), &(data_.solve_struct));
        }
        function_map::SolveInit(&(data_.options), &(data_.A), data_.perm_r.getRawPtr(),
                                data_.perm_c.getRawPtr(), as<SLUD::int_t>(nrhs), &(data_.lu),
                                &(data_.grid), &(data_.solve_struct));
        // Flag that we can reuse this solve_struct unless another
        // symbolicFactorization is called between here and the next
        // solve.
        same_solve_struct_ = true;
      }

      int ierr = 0; // returned error code
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

        function_map::gstrs(as<SLUD::int_t>(this->globalNumRows_), &(data_.lu),
                            &(data_.scale_perm), &(data_.grid), bvals_.getRawPtr(),
                            as<SLUD::int_t>(local_len_rhs), as<SLUD::int_t>(first_global_row_b),
                            as<SLUD::int_t>(local_len_rhs), as<int>(nrhs),
                            &(data_.solve_struct), &(data_.stat), &ierr);
      } // end block for solve time

      TEUCHOS_TEST_FOR_EXCEPTION( ierr < 0,
                          std::runtime_error,
                          "Argument " << -ierr << " to gstrs had an illegal value" );

      // "Un-scale" the solution so that it is a solution of the original system
      // if( data_.options.trans == SLUD::NOTRANS ){
      //   if( data_.colequ ){    // column equilibration has been done on AC
      //  // scale bxvals_ by diag(C)
      //  Util::scale(bxvals_(), as<size_t>(len_rhs), ldbx_, data_.C(),
      //              SLUD::slu_mt_mult<slu_type,magnitude_type>());
      //   }
      // } else if( data_.rowequ ){               // row equilibration has been done on AC
      //   // scale bxvals_ by diag(R)
      //   Util::scale(bxvals_(), as<size_t>(len_rhs), ldbx_, data_.R(),
      //            SLUD::slu_mt_mult<slu_type,magnitude_type>());
      // }
      {                         // permute B to a solution of the original system
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif
        SLUD::int_t ld = as<SLUD::int_t>(local_len_rhs);
        function_map::permute_Dense_Matrix(as<SLUD::int_t>(first_global_row_b),
                                           as<SLUD::int_t>(local_len_rhs),
                                           data_.solve_struct.row_to_proc,
                                           data_.solve_struct.inv_perm_c,
                                           bvals_.getRawPtr(), ld,
                                           xvals_.getRawPtr(), ld,
                                           as<int>(nrhs),
                                           &(data_.grid));
      }
    }

    /* Update X's global values */
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      typedef Util::put_1d_data_helper<MultiVecAdapter<Vector>,slu_type> put_helper;
      put_helper::do_put(X,
                         xvals_(),
                         local_len_rhs,
                         Teuchos::ptrInArg(*superlu_rowmap_));
    }

    return EXIT_SUCCESS;
  }


  template <class Matrix, class Vector>
  bool
  Superludist<Matrix,Vector>::matrixShapeOK_impl() const
  {
    // SuperLU_DIST requires square matrices
    return( this->globalNumRows_ == this->globalNumCols_ );
  }


  template <class Matrix, class Vector>
  void
  Superludist<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterEntryValidator;

    RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

    if( parameterList->isParameter("npcol") || parameterList->isParameter("nprow") ){
      TEUCHOS_TEST_FOR_EXCEPTION( !(parameterList->isParameter("nprow") &&
                            parameterList->isParameter("npcol")),
                          std::invalid_argument,
                          "nprow and npcol must be set together" );

      SLUD::int_t nprow = parameterList->template get<SLUD::int_t>("nprow");
      SLUD::int_t npcol = parameterList->template get<SLUD::int_t>("npcol");

      TEUCHOS_TEST_FOR_EXCEPTION( nprow * npcol > this->getComm()->getSize(),
                          std::invalid_argument,
                          "nprow and npcol combination invalid" );

      if( (npcol != data_.grid.npcol) || (nprow != data_.grid.nprow) ){
        // De-allocate the default grid that was initialized in the constructor
        SLUD::superlu_gridexit(&(data_.grid));
        // Create a new grid
        SLUD::superlu_gridinit(data_.mat_comm, nprow, npcol, &(data_.grid));
      } // else our grid has not changed size since the last initialization
    }

    TEUCHOS_TEST_FOR_EXCEPTION( this->control_.useTranspose_,
                        std::invalid_argument,
                        "SuperLU_DIST does not support solving the tranpose system" );

    data_.options.Trans = SLUD::NOTRANS; // should always be set this way;

    // TODO: Uncomment when supported
    // bool equil = parameterList->get<bool>("Equil", true);
    // data_.options.Equil = equil ? SLUD::YES : SLUD::NO;
    data_.options.Equil = SLUD::NO;

    if( parameterList->isParameter("ColPerm") ){
      RCP<const ParameterEntryValidator> colperm_validator = valid_params->getEntry("ColPerm").validator();
      parameterList->getEntry("ColPerm").setValidator(colperm_validator);

      data_.options.ColPerm = getIntegralValue<SLUD::colperm_t>(*parameterList, "ColPerm");
    }

    // Always use the "NOROWPERM" option to avoid a serial bottleneck
    // with the weighted bipartite matching algorithm used for the
    // "LargeDiag" RowPerm.  Note the inconsistency with the SuperLU
    // User guide (which states that the value should be "NATURAL").
    data_.options.RowPerm = SLUD::NOROWPERM;

    // TODO: Uncomment when supported
    // if( parameterList->isParameter("IterRefine") ){
    //   RCP<const ParameterEntryValidator> iter_refine_validator = valid_params->getEntry("IterRefine").validator();
    //   parameterList->getEntry("IterRefine").setValidator(iter_refine_validator);

    //   data_.options.IterRefine = getIntegralValue<SLUD::IterRefine_t>(*parameterList, "IterRefine");
    // }
    data_.options.IterRefine = SLUD::NOREFINE;

    bool replace_tiny = parameterList->get<bool>("ReplaceTinyPivot", true);
    data_.options.ReplaceTinyPivot = replace_tiny ? SLUD::YES : SLUD::NO;
  }


  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  Superludist<Matrix,Vector>::getValidParameters_impl() const
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

      Teuchos::RCP<EnhancedNumberValidator<SLUD::int_t> > col_row_validator
        = Teuchos::rcp( new EnhancedNumberValidator<SLUD::int_t>() );
      col_row_validator->setMin(1);

      pl->set("npcol", data_.grid.npcol,
              "Number of columns in the processor grid. "
              "Must be set with nprow", col_row_validator);
      pl->set("nprow", data_.grid.nprow,
              "Number of rows in the SuperLU_DIST processor grid. "
              "Must be set together with npcol", col_row_validator);

      // validator will catch any value besides NOTRANS
      setStringToIntegralParameter<SLUD::trans_t>("Trans", "NOTRANS",
                                                  "Solve for the transpose system or not",
                                                  tuple<string>("NOTRANS"),
                                                  tuple<string>("Do not solve with transpose"),
                                                  tuple<SLUD::trans_t>(SLUD::NOTRANS),
                                                  pl.getRawPtr());

      // TODO: uncomment when supported
      // pl->set("Equil", false, "Whether to equilibrate the system before solve");

      // TODO: uncomment when supported
      // setStringToIntegralParameter<SLUD::IterRefine_t>("IterRefine", "NOREFINE",
      //                                                     "Type of iterative refinement to use",
      //                                                     tuple<string>("NOREFINE", "DOUBLE"),
      //                                                     tuple<string>("Do not use iterative refinement",
      //                                                                   "Do double iterative refinement"),
      //                                                     tuple<SLUD::IterRefine_t>(SLUD::NOREFINE,
      //                                                                               SLUD::DOUBLE),
      //                                                     pl.getRawPtr());

      pl->set("ReplaceTinyPivot", true,
              "Specifies whether to replace tiny diagonals during LU factorization");

      setStringToIntegralParameter<SLUD::colperm_t>("ColPerm", "PARMETIS",
                                                    "Specifies how to permute the columns of the "
                                                    "matrix for sparsity preservation",
                                                    tuple<string>("NATURAL", "PARMETIS"),
                                                    tuple<string>("Natural ordering",
                                                                  "ParMETIS ordering on A^T + A"),
                                                    tuple<SLUD::colperm_t>(SLUD::NATURAL,
                                                                           SLUD::PARMETIS),
                                                    pl.getRawPtr());

      valid_params = pl;
    }

    return valid_params;
  }


  template <class Matrix, class Vector>
  void
  Superludist<Matrix,Vector>::get_default_grid_size(int nprocs,
                                                    SLUD::int_t& nprow,
                                                    SLUD::int_t& npcol) const {
    TEUCHOS_TEST_FOR_EXCEPTION( nprocs < 1,
                        std::invalid_argument,
                        "Number of MPI processes must be at least 1" );
    SLUD::int_t c, r = 1;
    while( r*r <= nprocs ) r++;
    nprow = npcol = --r;                // fall back to square grid
    c = nprocs / r;
    while( (r--)*c != nprocs ){
      c = nprocs / r;           // note integer division
    }
    ++r;
    // prefer the square grid over a single row (which will only happen
    // in the case of a prime nprocs
    if( r > 1 || nprocs < 9){   // nprocs < 9 is a heuristic for the small cases
      nprow = r;
      npcol = c;
    }
  }


  template <class Matrix, class Vector>
  bool
  Superludist<Matrix,Vector>::loadA_impl(EPhase current_phase){
    // Extract the necessary information from mat and call SLU function
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ptrInArg;
    using Teuchos::as;

    using SLUD::int_t;

#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Cleanup old store memory if it's non-NULL
    if( data_.A.Store != NULL ){
      SLUD::Destroy_SuperMatrix_Store_dist( &(data_.A) );
      data_.A.Store = NULL;
    }

    Teuchos::RCP<const MatrixAdapter<Matrix> > redist_mat
      = this->matrixA_->get(ptrInArg(*superlu_rowmap_));

    int_t l_nnz, l_rows, g_rows, g_cols, fst_global_row;
    l_nnz  = as<int_t>(redist_mat->getLocalNNZ());
    l_rows = as<int_t>(redist_mat->getLocalNumRows());
    g_rows = as<int_t>(redist_mat->getGlobalNumRows());
    g_cols = g_rows;            // we deal with square matrices
    fst_global_row = as<int_t>(superlu_rowmap_->getMinGlobalIndex());

    nzvals_.resize(l_nnz);
    colind_.resize(l_nnz);
    rowptr_.resize(l_rows + 1);

    int_t nnz_ret = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

      Util::get_crs_helper<
      MatrixAdapter<Matrix>,
        slu_type, int_t, int_t >::do_get(redist_mat.ptr(),
                                         nzvals_(), colind_(), rowptr_(),
                                         nnz_ret,
                                         ptrInArg(*superlu_rowmap_),
                                         ARBITRARY);
  }

    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != l_nnz,
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");

  // Get the SLU data type for this type of matrix
  SLUD::Dtype_t dtype = type_map::dtype;

  if( in_grid_ ){
    function_map::create_CompRowLoc_Matrix(&(data_.A),
                                           g_rows, g_cols,
                                           l_nnz, l_rows, fst_global_row,
                                           nzvals_.getRawPtr(),
                                           colind_.getRawPtr(),
                                           rowptr_.getRawPtr(),
                                           SLUD::SLU_NR_loc,
                                           dtype, SLUD::SLU_GE);
  }

  return true;
}


  template<class Matrix, class Vector>
  const char* Superludist<Matrix,Vector>::name = "SuperLU_DIST";


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUDIST_DEF_HPP
