// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include <Teuchos_Details_MpiTypeTraits.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Superludist_TypeMap.hpp"
#include "Amesos2_Util.hpp"


namespace Amesos2 {


  template <class Matrix, class Vector>
  Superludist<Matrix,Vector>::Superludist(Teuchos::RCP<const Matrix> A,
                                          Teuchos::RCP<Vector> X,
                                          Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::Superludist,Matrix,Vector>(A, X, B)
    , bvals_()
    , xvals_()
    , in_grid_(false)
    , force_symbfact_(false)
    , is_contiguous_(true)
  {
    using Teuchos::Comm;
    // It's OK to depend on MpiComm explicitly here, because
    // SuperLU_DIST requires MPI anyway.
    using Teuchos::MpiComm;
    using Teuchos::outArg;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    typedef global_ordinal_type GO;
    typedef Tpetra::Map<local_ordinal_type, GO, node_type> map_type;

    ////////////////////////////////////////////
    // Set up the SuperLU_DIST processor grid //
    ////////////////////////////////////////////

    RCP<const Comm<int> > comm = this->getComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    SLUD::int_t nprow, npcol;
    get_default_grid_size (numProcs, nprow, npcol);

    {
      // FIXME (mfh 16 Dec 2014) getComm() just returns
      // matrixA_->getComm(), so it's not clear why we need to ask for
      // the matrix's communicator separately here.
      RCP<const Comm<int> > matComm = this->matrixA_->getComm ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        matComm.is_null (), std::logic_error, "Amesos2::Superlustdist "
        "constructor: The matrix's communicator is null!");
      RCP<const MpiComm<int> > matMpiComm =
        rcp_dynamic_cast<const MpiComm<int> > (matComm);
      // FIXME (mfh 16 Dec 2014) If the matrix's communicator is a
      // SerialComm, we probably could just use MPI_COMM_SELF here.
      // I'm not sure if SuperLU_DIST is smart enough to handle that
      // case, though.
      TEUCHOS_TEST_FOR_EXCEPTION(
        matMpiComm.is_null (), std::logic_error, "Amesos2::Superlustdist "
        "constructor: The matrix's communicator is not an MpiComm!");
      TEUCHOS_TEST_FOR_EXCEPTION(
        matMpiComm->getRawMpiComm ().is_null (), std::logic_error, "Amesos2::"
        "Superlustdist constructor: The matrix's communicator claims to be a "
        "Teuchos::MpiComm<int>, but its getRawPtrComm() method returns "
        "Teuchos::null!  This means that the underlying MPI_Comm doesn't even "
        "exist, which likely implies that the Teuchos::MpiComm was constructed "
        "incorrectly.  It means something different than if the MPI_Comm were "
        "MPI_COMM_NULL.");
      MPI_Comm rawMpiComm = (* (matMpiComm->getRawMpiComm ())) ();
      data_.mat_comm = rawMpiComm;
      // This looks a bit like ScaLAPACK's grid initialization (which
      // technically takes place in the BLACS, not in ScaLAPACK
      // proper). See http://netlib.org/scalapack/slug/node34.html.
      // The main difference is that SuperLU_DIST depends explicitly
      // on MPI, while the BLACS hides its communication protocol.
      SLUD::superlu_gridinit(data_.mat_comm, nprow, npcol, &(data_.grid));
    }

    ////////////////////////////////////////////////////////
    // Set some default parameters.                       //
    //                                                    //
    // Must do this after grid has been created in        //
    // case user specifies the nprow and npcol parameters //
    ////////////////////////////////////////////////////////
    SLUD::set_default_options_dist(&data_.options);

    RCP<ParameterList> default_params =
      parameterList (* (this->getValidParameters ()));
    this->setParameters (default_params);

    // Set some internal options
    data_.options.Fact = SLUD::DOFACT;
    data_.equed = SLUD::NOEQUIL; // No equilibration has yet been performed
    data_.options.SolveInitialized  = SLUD::NO;
    data_.options.RefineInitialized = SLUD::NO;
    data_.rowequ = false;
    data_.colequ = false;
    data_.perm_r.resize(this->globalNumRows_);
    data_.perm_c.resize(this->globalNumCols_);
    data_.largediag_mc64_job = 4;
    for (global_size_type i = 0; i < this->globalNumRows_; i++)
      data_.perm_r[i] = i;
    for (global_size_type i = 0; i < this->globalNumCols_; i++)
      data_.perm_c[i] = i;

    ////////////////////////////////////////////////////////////////
    // Set up a communicator for the parallel column ordering and //
    // parallel symbolic factorization.                           //
    ////////////////////////////////////////////////////////////////
    data_.symb_comm = MPI_COMM_NULL;

    // domains is the next power of 2 less than nprow*npcol.  This
    // value will be used for creating an MPI communicator for the
    // pre-ordering and symbolic factorization methods.
    data_.domains = (int) ( pow(2.0, floor(log10((double)nprow*npcol)/log10(2.0))) );

    const int color = (myRank < data_.domains) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split (data_.mat_comm, color, myRank, &(data_.symb_comm));

    //////////////////////////////////////////////////////////////////////
    // Set up a row Map that only includes processes that are in the
    // SuperLU process grid.  This will be used for redistributing A.
    //////////////////////////////////////////////////////////////////////

    // mfh 16 Dec 2014: We could use createWeightedContigMapWithNode
    // with myProcParticipates as the weight, but that costs an extra
    // all-reduce.

    // Set to 1 if I am in the grid, and I get some of the matrix rows.
    int myProcParticipates = 0;
    if (myRank < nprow * npcol) {
      in_grid_ = true;
      myProcParticipates = 1;
    }

    // Compute how many processes in the communicator belong to the
    // process grid.
    int numParticipatingProcs = 0;
    reduceAll<int, int> (*comm, REDUCE_SUM, myProcParticipates,
                         outArg (numParticipatingProcs));
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->globalNumRows_ != 0 && numParticipatingProcs == 0,
      std::logic_error, "Amesos2::Superludist constructor: The matrix has "
      << this->globalNumRows_ << " > 0 global row(s), but no processes in the "
      "communicator want to participate in its factorization!  nprow = "
      << nprow << " and npcol = " << npcol << ".");

    // Divide up the rows among the participating processes.
    size_t myNumRows = 0;
    {
      const GO GNR = static_cast<GO> (this->globalNumRows_);
      const GO quotient = (numParticipatingProcs == 0) ? static_cast<GO> (0) :
        GNR / static_cast<GO> (numParticipatingProcs);
      const GO remainder =
        GNR - quotient * static_cast<GO> (numParticipatingProcs);
      const GO lclNumRows = (static_cast<GO> (myRank) < remainder) ?
        (quotient + static_cast<GO> (1)) : quotient;
      myNumRows = static_cast<size_t> (lclNumRows);
    }

    // TODO: might only need to initialize if parallel symbolic factorization is requested.
    const GO indexBase = this->matrixA_->getRowMap ()->getIndexBase ();
    superlu_rowmap_ =
      rcp (new map_type (this->globalNumRows_, myNumRows, indexBase, comm));

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
      /// These are created by superlu malloc and should be deallocated by superlu free
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
      SUPERLU_FREE( data_.sizes );
      SUPERLU_FREE( data_.fstVtxSep );
#else
      free( data_.sizes );
      free( data_.fstVtxSep );
#endif
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
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
        SUPERLU_FREE( data_.pslu_freeable.xlsub );
        SUPERLU_FREE( data_.pslu_freeable.lsub );
#else
        free( data_.pslu_freeable.xlsub );
        free( data_.pslu_freeable.lsub );
#endif
      }
      if ( data_.pslu_freeable.xusub != NULL ){
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
        SUPERLU_FREE( data_.pslu_freeable.xusub );
        SUPERLU_FREE( data_.pslu_freeable.usub );
#else
        free( data_.pslu_freeable.xusub );
        free( data_.pslu_freeable.usub );
#endif
      }
      if ( data_.pslu_freeable.supno_loc != NULL ){
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
        SUPERLU_FREE( data_.pslu_freeable.supno_loc );
        SUPERLU_FREE( data_.pslu_freeable.xsup_beg_loc );
        SUPERLU_FREE( data_.pslu_freeable.xsup_end_loc );
#else
        free( data_.pslu_freeable.supno_loc );
        free( data_.pslu_freeable.xsup_beg_loc );
        free( data_.pslu_freeable.xsup_end_loc );
#endif
      }
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
      SUPERLU_FREE( data_.pslu_freeable.globToLoc );
#else
      free( data_.pslu_freeable.globToLoc );
#endif
    }

    SLUD::PStatFree( &(data_.stat) ) ;

    // Teuchos::Arrays will free R, C, perm_r, and perm_c
    // SLUD::D::ScalePermstructFree(&(data_.scale_perm));

    if ( data_.options.SolveInitialized == SLUD::YES )
      function_map::SolveFinalize(&(data_.options), &(data_.solve_struct));

    // gridexit of an older version frees SuperLU_MPI_DOUBLE_COMPLE,
    // which could cause an issue if there are still active instances of superludist?
    SLUD::superlu_gridexit(&(data_.grid)); // TODO: are there any
                                           // cases where grid
                                           // wouldn't be initialized?

    if ( data_.symb_comm != MPI_COMM_NULL ) MPI_Comm_free(&(data_.symb_comm));
  }

  template<class Matrix, class Vector>
  void
  Superludist<Matrix,Vector>::computeRowPermutationLargeDiagMC64(SLUD::SuperMatrix& GA)
  {
    int job = data_.largediag_mc64_job;
    if (job == 5)
    {
      data_.R1.resize(data_.A.nrow);
      data_.C1.resize(data_.A.ncol);
    }

    SLUD::NCformat *GAstore = (SLUD::NCformat*) GA.Store;
    SLUD::int_t* colptr = GAstore->colptr;
    SLUD::int_t* rowind = GAstore->rowind;
    SLUD::int_t nnz = GAstore->nnz;
    slu_type *a_GA = (slu_type *) GAstore->nzval;
    MPI_Datatype mpi_dtype = Teuchos::Details::MpiTypeTraits<magnitude_type>::getType();
    MPI_Datatype mpi_itype = Teuchos::Details::MpiTypeTraits<SLUD::int_t>::getType();

    int iinfo = 0;
    if ( !data_.grid.iam ) { /* Process 0 finds a row permutation */
      iinfo = function_map::ldperm_dist(job, data_.A.nrow, nnz, colptr, rowind, a_GA,
              data_.perm_r.getRawPtr(), data_.R1.getRawPtr(), data_.C1.getRawPtr());

      MPI_Bcast( &iinfo, 1, MPI_INT, 0, data_.grid.comm );
      if ( iinfo == 0 ) {
          MPI_Bcast( data_.perm_r.getRawPtr(), data_.A.nrow, mpi_itype, 0, data_.grid.comm );
          if ( job == 5 && data_.options.Equil ) {
              MPI_Bcast( data_.R1.getRawPtr(), data_.A.nrow, mpi_dtype, 0, data_.grid.comm );
              MPI_Bcast( data_.C1.getRawPtr(), data_.A.ncol, mpi_dtype, 0, data_.grid.comm );
          }
      }
    } else {
      MPI_Bcast( &iinfo, 1, mpi_int_t, 0, data_.grid.comm );
      if ( iinfo == 0 ) {
        MPI_Bcast( data_.perm_r.getRawPtr(), data_.A.nrow, mpi_itype, 0, data_.grid.comm );
        if ( job == 5 && data_.options.Equil ) {
            MPI_Bcast( data_.R1.getRawPtr(), data_.A.nrow, mpi_dtype, 0, data_.grid.comm );
            MPI_Bcast( data_.C1.getRawPtr(), data_.A.ncol, mpi_dtype, 0, data_.grid.comm );
        }
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION( iinfo != 0,
                        std::runtime_error,
                        "SuperLU_DIST pre-ordering failed to compute row perm with "
                        << iinfo << std::endl);

    if (job == 5)
    {
      for (SLUD::int_t i = 0; i < data_.A.nrow; ++i) data_.R1[i] = exp(data_.R1[i]);
      for (SLUD::int_t i = 0; i < data_.A.ncol; ++i) data_.C1[i] = exp(data_.C1[i]);
    }
  }


  template<class Matrix, class Vector>
  int
  Superludist<Matrix,Vector>::preOrdering_impl()
  {
    if (data_.options.RowPerm == SLUD::NOROWPERM) {
      SLUD::int_t slu_rows_ub = Teuchos::as<SLUD::int_t>(this->globalNumRows_);
      for( SLUD::int_t i = 0; i < slu_rows_ub; ++i ) data_.perm_r[i] = i;
    }
    else if (data_.options.RowPerm == SLUD::LargeDiag_MC64) {
      if (!force_symbfact_)
        // defer to numerical factorization because row permutation requires the matrix values
        return (EXIT_SUCCESS + 1);
    }
    // loadA_impl();                    // Refresh matrix values

    if( in_grid_ ){
      // If this function has been called at least once, then the
      // sizes, and fstVtxSep arrays were allocated in
      // get_perm_c_parmetis.  Delete them before calling that
      // function again.  These arrays will also be dealloc'd in the
      // deconstructor.
      if( this->status_.getNumPreOrder() > 0 ){
#if defined(AMESOS2_ENABLES_SUPERLUDIST_VERSION5_AND_HIGHER)
        SUPERLU_FREE( data_.sizes );
        SUPERLU_FREE( data_.fstVtxSep );
#else
        free( data_.sizes );
        free( data_.fstVtxSep );
#endif
      }
      float info = 0.0;
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor preOrderTime( this->timers_.preOrderTime_ );
#endif
        info = SLUD::get_perm_c_parmetis( &(data_.A),
                                          data_.perm_r.getRawPtr(), data_.perm_c.getRawPtr(),
                                          data_.grid.nprow * data_.grid.npcol, data_.domains,
                                          &(data_.sizes), &(data_.fstVtxSep),
                                          &(data_.grid), &(data_.symb_comm) );
      }
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
    if (!force_symbfact_) {
       if (data_.options.RowPerm == SLUD::LargeDiag_MC64) {
          // defer to numerical factorization because row permutation requires the matrix values
          return (EXIT_SUCCESS + 1);
       }
    }

    if( in_grid_ ){

      float info = 0.0;
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
#endif

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
        info = SLUD::symbfact_dist(&(data_.options), (data_.grid.nprow) * (data_.grid.npcol),
                                   data_.domains, &(data_.A), data_.perm_c.getRawPtr(),
                                   data_.perm_r.getRawPtr(), data_.sizes,
                                   data_.fstVtxSep, &(data_.pslu_freeable),
                                   &(data_.grid.comm), &(data_.symb_comm),
                                   &(data_.mem_usage));

#else
        info = SLUD::symbfact_dist((data_.grid.nprow) * (data_.grid.npcol),
                                   data_.domains, &(data_.A), data_.perm_c.getRawPtr(),
                                   data_.perm_r.getRawPtr(), data_.sizes,
                                   data_.fstVtxSep, &(data_.pslu_freeable),
                                   &(data_.grid.comm), &(data_.symb_comm),
                                   &(data_.mem_usage));
#endif
      }
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
    SLUD::SuperMatrix GA;      /* Global A in NC format */
    bool need_value = false;

    if( in_grid_ ) {
      if( data_.options.Equil == SLUD::YES ) {
        SLUD::int_t info = 0;

        // Compute scaling
        data_.R.resize(this->globalNumRows_);
        data_.C.resize(this->globalNumCols_);
        function_map::gsequ_loc(&(data_.A), data_.R.getRawPtr(), data_.C.getRawPtr(),
                                &(data_.rowcnd), &(data_.colcnd), &(data_.amax), &info, &(data_.grid));

        // Apply the scalings
        function_map::laqgs_loc(&(data_.A), data_.R.getRawPtr(), data_.C.getRawPtr(),
                                data_.rowcnd, data_.colcnd, data_.amax,
                                &(data_.equed));

        data_.rowequ = (data_.equed == SLUD::ROW) || (data_.equed == SLUD::BOTH);
        data_.colequ = (data_.equed == SLUD::COL) || (data_.equed == SLUD::BOTH);

        // Compute and apply the row permutation
        if (data_.options.RowPerm == SLUD::LargeDiag_MC64) {
          // Create a column-order copy of A
          need_value = true;
          SLUD::D::pdCompRow_loc_to_CompCol_global(true, &data_.A, &data_.grid, &GA);

          // Compute row permutation
          computeRowPermutationLargeDiagMC64(GA);

          // Here we do symbolic factorization
          force_symbfact_ = true;
          preOrdering_impl();
          symbolicFactorization_impl();
          force_symbfact_ = false;

          // Apply row-permutation scaling for job=5
          // Here we do it manually to bypass the threshold check in laqgs_loc
          if (data_.largediag_mc64_job == 5)
          {
            SLUD::NRformat_loc *Astore  = (SLUD::NRformat_loc*) data_.A.Store;
            slu_type *a = (slu_type*) Astore->nzval;
            SLUD::int_t m_loc   = Astore->m_loc;
            SLUD::int_t fst_row = Astore->fst_row;
            SLUD::int_t i, j, irow = fst_row, icol;

            /* Scale the distributed matrix further.
             A <-- diag(R1)*A*diag(C1)            */
            SLUD::slu_dist_mult<slu_type, magnitude_type> mult_op;
            for (j = 0; j < m_loc; ++j) {
              for (i = rowptr_view_.data()[j]; i < rowptr_view_.data()[j+1]; ++i) {
                  icol = colind_view_.data()[i];
                  a[i] = mult_op(a[i], data_.R1[irow] * data_.C1[icol]);
              }
              ++irow;
            }

            /* Multiply together the scaling factors */
            if ( data_.rowequ ) for (i = 0; i < data_.A.nrow; ++i) data_.R[i] *= data_.R1[i];
            else for (i = 0; i < data_.A.nrow; ++i) data_.R[i] = data_.R1[i];
            if ( data_.colequ ) for (i = 0; i < data_.A.ncol; ++i) data_.C[i] *= data_.C1[i];
            else for (i = 0; i < data_.A.ncol; ++i) data_.C[i] = data_.C1[i];

            data_.rowequ = data_.colequ = 1;
          }
        }
      }

      // Apply the column ordering, so that AC is the column-permuted A, and compute etree
      size_t nnz_loc = ((SLUD::NRformat_loc*)data_.A.Store)->nnz_loc;
      for( size_t i = 0; i < nnz_loc; ++i ) colind_view_(i) = data_.perm_c[colind_view_(i)];

      // Distribute data from the symbolic factorization
      if( same_symbolic_ ){
        // Note: with the SamePattern_SameRowPerm options, it does not
        // matter that the glu_freeable member has never been
        // initialized, because it is never accessed.  It is a
        // placeholder arg.  The real work is done in data_.lu
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
        data_.options.Fact = SLUD::SamePattern_SameRowPerm;
        function_map::pdistribute(&(data_.options),
                                  as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                  &(data_.A), &(data_.scale_perm),
                                  &(data_.glu_freeable), &(data_.lu),
                                  &(data_.grid));
#else
        function_map::pdistribute(SLUD::SamePattern_SameRowPerm,
                                  as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                  &(data_.A), &(data_.scale_perm),
                                  &(data_.glu_freeable), &(data_.lu),
                                  &(data_.grid));
#endif
      } else {
#if (SUPERLU_DIST_MAJOR_VERSION > 7)
        data_.options.Fact = SLUD::DOFACT;
        function_map::dist_psymbtonum(&(data_.options),
                                      as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                      &(data_.A), &(data_.scale_perm),
                                      &(data_.pslu_freeable), &(data_.lu),
                                      &(data_.grid));
#else
        function_map::dist_psymbtonum(SLUD::DOFACT,
                                      as<SLUD::int_t>(this->globalNumRows_), // aka "n"
                                      &(data_.A), &(data_.scale_perm),
                                      &(data_.pslu_freeable), &(data_.lu),
                                      &(data_.grid));
#endif
      }

      // Retrieve the normI of A (required by gstrf).
      bool notran = (data_.options.Trans == SLUD::NOTRANS);
      magnitude_type anorm = function_map::plangs((notran ? (char *)"1" : (char *)"I"), &(data_.A), &(data_.grid));

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

    if (need_value)
      SLUD::Destroy_CompCol_Matrix_dist(&GA);

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
    const size_t local_len_rhs = superlu_rowmap_->getLocalNumElements();
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

      // Apply row-scaling if requested
      if (data_.options.Equil == SLUD::YES && data_.rowequ) {
        SLUD::int_t ld = as<SLUD::int_t>(local_len_rhs);
        SLUD::slu_dist_mult<slu_type, magnitude_type> mult_op;
        for(global_size_type j = 0; j < nrhs; ++j) {
          for(size_t i = 0; i < local_len_rhs; ++i) {
            bvals_[i + j*ld] = mult_op(bvals_[i + j*ld], data_.R[first_global_row_b + i]);
          }
        }
      }

      // Solve
      int ierr = 0; // returned error code
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

#if (SUPERLU_DIST_MAJOR_VERSION > 7)
        function_map::gstrs(&(data_.options), as<SLUD::int_t>(this->globalNumRows_), &(data_.lu),
                            &(data_.scale_perm), &(data_.grid), bvals_.getRawPtr(),
                            as<SLUD::int_t>(local_len_rhs), as<SLUD::int_t>(first_global_row_b),
                            as<SLUD::int_t>(local_len_rhs), as<int>(nrhs),
                            &(data_.solve_struct), &(data_.stat), &ierr);
#else
        function_map::gstrs(as<SLUD::int_t>(this->globalNumRows_), &(data_.lu),
                            &(data_.scale_perm), &(data_.grid), bvals_.getRawPtr(),
                            as<SLUD::int_t>(local_len_rhs), as<SLUD::int_t>(first_global_row_b),
                            as<SLUD::int_t>(local_len_rhs), as<int>(nrhs),
                            &(data_.solve_struct), &(data_.stat), &ierr);
#endif
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

      // Apply col-scaling if requested
      if (data_.options.Equil == SLUD::YES && data_.colequ) {
        SLUD::int_t ld = as<SLUD::int_t>(local_len_rhs);
        SLUD::slu_dist_mult<slu_type, magnitude_type> mult_op;
        for(global_size_type j = 0; j < nrhs; ++j) {
          for(size_t i = 0; i < local_len_rhs; ++i) {
            xvals_[i + j*ld] = mult_op(xvals_[i + j*ld], data_.C[first_global_row_b + i]);
          }
        }
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

    // Equilbration option
    bool equil = parameterList->get<bool>("Equil", false);
    data_.options.Equil = equil ? SLUD::YES : SLUD::NO;

    if( parameterList->isParameter("RowPerm") ){
      RCP<const ParameterEntryValidator> rowperm_validator = valid_params->getEntry("RowPerm").validator();
      parameterList->getEntry("RowPerm").setValidator(rowperm_validator);

      data_.options.RowPerm = getIntegralValue<SLUD::rowperm_t>(*parameterList, "RowPerm");
    }

    if( parameterList->isParameter("LargeDiag_MC64-Options") ){
      data_.largediag_mc64_job = parameterList->template get<int>("LargeDiag_MC64-Options");
    }

    if( parameterList->isParameter("ColPerm") ){
      RCP<const ParameterEntryValidator> colperm_validator = valid_params->getEntry("ColPerm").validator();
      parameterList->getEntry("ColPerm").setValidator(colperm_validator);

      data_.options.ColPerm = getIntegralValue<SLUD::colperm_t>(*parameterList, "ColPerm");
    }

    // TODO: Uncomment when supported
    // if( parameterList->isParameter("IterRefine") ){
    //   RCP<const ParameterEntryValidator> iter_refine_validator = valid_params->getEntry("IterRefine").validator();
    //   parameterList->getEntry("IterRefine").setValidator(iter_refine_validator);
    //   data_.options.IterRefine = getIntegralValue<SLUD::IterRefine_t>(*parameterList, "IterRefine");
    // }
    data_.options.IterRefine = SLUD::NOREFINE;

    bool replace_tiny = parameterList->get<bool>("ReplaceTinyPivot", true);
    data_.options.ReplaceTinyPivot = replace_tiny ? SLUD::YES : SLUD::NO;

    if( parameterList->isParameter("IsContiguous") ){
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }
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
    using Teuchos::setIntParameter;
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

      // Equillbration
      pl->set("Equil", false, "Whether to equilibrate the system before solve");

      // TODO: uncomment when supported
      // setStringToIntegralParameter<SLUD::IterRefine_t>("IterRefine", "NOREFINE",
      //                                                     "Type of iterative refinement to use",
      //                                                     tuple<string>("NOREFINE", "DOUBLE"),
      //                                                     tuple<string>("Do not use iterative refinement",
      //                                                                   "Do double iterative refinement"),
      //                                                     tuple<SLUD::IterRefine_t>(SLUD::NOREFINE,
      //                                                                               SLUD::DOUBLE),
      //                                                     pl.getRawPtr());

      // Tiny pivot handling
      pl->set("ReplaceTinyPivot", true,
              "Specifies whether to replace tiny diagonals during LU factorization");

      // Row permutation
      setStringToIntegralParameter<SLUD::rowperm_t>("RowPerm", "NOROWPERM",
                                                    "Specifies how to permute the rows of the "
                                                    "matrix for sparsity preservation",
                                                    tuple<string>("NOROWPERM", "LargeDiag_MC64"),
                                                    tuple<string>("Natural ordering",
                                                                  "Duff/Koster algorithm"),
                                                    tuple<SLUD::rowperm_t>(SLUD::NOROWPERM,
                                                                           SLUD::LargeDiag_MC64),
                                                    pl.getRawPtr());

      setIntParameter("LargeDiag_MC64-Options", 4, "Options for RowPerm-LargeDiag_MC64", pl.getRawPtr());

      // Column permutation
      setStringToIntegralParameter<SLUD::colperm_t>("ColPerm", "PARMETIS",
                                                    "Specifies how to permute the columns of the "
                                                    "matrix for sparsity preservation",
                                                    tuple<string>("NATURAL", "PARMETIS"),
                                                    tuple<string>("Natural ordering",
                                                                  "ParMETIS ordering on A^T + A"),
                                                    tuple<SLUD::colperm_t>(SLUD::NATURAL,
                                                                           SLUD::PARMETIS),
                                                    pl.getRawPtr());

      pl->set("IsContiguous", true, "Whether GIDs contiguous");

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

    Kokkos::resize(nzvals_view_, l_nnz);
    Kokkos::resize(colind_view_, l_nnz);
    Kokkos::resize(rowptr_view_, l_rows + 1);
    int_t nnz_ret = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

      Util::get_crs_helper_kokkos_view<MatrixAdapter<Matrix>,
        host_value_type_array,host_ordinal_type_array, host_size_type_array >::do_get(
                                         redist_mat.ptr(),
                                         nzvals_view_, colind_view_, rowptr_view_,
                                         nnz_ret,
                                         ptrInArg(*superlu_rowmap_),
                                         ROOTED,
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
                                           nzvals_view_.data(),
                                           colind_view_.data(),
                                           rowptr_view_.data(),
                                           SLUD::SLU_NR_loc,
                                           dtype, SLUD::SLU_GE);
  }

  return true;
}


  template<class Matrix, class Vector>
  const char* Superludist<Matrix,Vector>::name = "SuperLU_DIST";


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUDIST_DEF_HPP
