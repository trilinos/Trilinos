// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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

#include "Amesos2_Util.hpp"


namespace Amesos {


  template <class Matrix, class Vector>
  Superludist<Matrix,Vector>::Superludist(Teuchos::RCP<Matrix> A,
					  Teuchos::RCP<Vector> X,
					  Teuchos::RCP<Vector> B)
    : SolverCore<Amesos::Superludist,Matrix,Vector>(A, X, B)
    , nzvals_()                 // initialization to empty arrays
    , colind_()
    , rowptr_()
    , bvals_()
    , xvals_()
    , in_grid_(false)
  {
    // Set some default parameters
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

    ////////////////////////////////////////////
    // Set up the SuperLU_DIST processor grid //
    ////////////////////////////////////////////

    int nprocs = this->getComm()->getSize();
    SLUD::int_t nprow, npcol;
    get_default_grid_size(nprocs, nprow, npcol);
    data_.mat_comm = dynamic_cast<const Teuchos::MpiComm<int>* >(this->matrixA_->getComm().getRawPtr())->getRawMpiComm()->operator()();
    SLUD::superlu_gridinit(data_.mat_comm, nprow, npcol, &(data_.grid));

    // Set up a communicator for the parallel column ordering and
    // parallel symbolic factorization.
    data_.symb_comm = MPI_COMM_NULL;
    int color = MPI_UNDEFINED;
    int my_rank = this->status_.myPID_;

    /* domains is the next power of 2 less than nprow*npcol.  This
     * value will be used for creating an MPI communicator for the
     * pre-ordering and symbolic factorization methods.
     */
    data_.domains = (int) ( pow(2, ((int) (log10((double)nprow*npcol)/log10(2.0)))) );

    if( this->status_.myPID_ < data_.domains ) color = 0;
    MPI_Comm_split (data_.mat_comm, color, my_rank, &(data_.symb_comm));

    //////////////////////////////////////////////////////////////////////
    // Set up a row map that maps to only processors that are in the    //
    // SuperLU processor grid.  This will be used for redistributing A. //
    //////////////////////////////////////////////////////////////////////
  
    int my_weight = 0;
    if( this->status_.myPID_ < nprow * npcol ){
      in_grid_ = true; my_weight = 1; // I am in the grid, and I get some of the matrix rows
    }
    // TODO: might only need to initialize if parallel symbolic factorization is requested.
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
    /* Free SuperLU_MT data_types
     * - Matrices
     * - Vectors
     * - Stat object
     * - ScalePerm, LUstruct, grid, and solve objects
     *
     * Note: the function definitions are the same regardless whether
     * complex or real, so we arbitrarily use the D namespace
     */
    if ( this->getNumPreOrder() > 0 ){
      free( data_.sizes );
      free( data_.fstVtxSep );
    }

    // Storage is initialized in numericFactorization_impl()
    if ( this->getNumNumericFact() > 0 ){
      // Cleanup old matrix store memory if it's non-NULL
      if( data_.A.Store != NULL ){
	SLUD::Destroy_SuperMatrix_Store_dist( &(data_.A) );
      }

      // Our Teuchos::Array's will destroy rowind, colptr, and nzval for us
      function_map::Destroy_LU(this->globalNumRows_, &(data_.grid), &(data_.lu));
      function_map::LUstructFree(&(data_.lu));
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

    // The request for symbolic factorization will dictate what type
    // of storage is allocated for the matrix.
    if( data_.options.ParSymbFact == SLUD::NO ){
      // // First, do equilibration if requested
      // if( data_.options.Equil == SLUD::YES ){
      //        data_.R.resize(this->globalNumRows_);
      //        data_.C.resize(this->globalNumCols_);

      //        int info;

      //        function_map::gsequ(&(data_.A), data_.R.getRawPtr(),
      //                            data_.C.getRawPtr(), &(data_.rowcnd),
      //                            &(data_.colcnd), &(data_.amax), &info);
      //        TEST_FOR_EXCEPTION( info != 0,
      //                            std::runtime_error,
      //                            "SuperLU_MT gsequ returned with status " << info );
      //        // Scalings will be applied in numericFactorization
      // }

      // // Next find a row permutation.
      // if( data_.options.RowPerm == SLUD::NATURAL ){
      //        for( int i = 0; i < this->globalNumRows_; ++i ) data_.perm_r[i] = i;
      // } else if( data_.options.RowPerm == SLUD::LargeDiag ){
      //        // Currently this must be done serially.
      //        matrix_helper::createCcsMatrix(this->matrixA_.ptr(),
      //                                       nzvals_(), colind_(), rowptr_(),
      //                                       Teuchos::ptrFromRef(data_.A),
      //                                       this->timers_.mtxRedistTime_);
      //        // Everyone has a global matrix copy, so we just all compute
      //        // the row permutation
      //        if( data_.options.Equil == SLUD::YES ){
      //          data_.R1.resize(this->globalNumRows_);
      //          data_.C1.resize(this->globalNumCols_);
      //        }
      //        // TODO: Do we need to apply R and C to A's nzvals before
      //        // finding a row permutation?
      //        function_map::ldperm(5, this->globalNumRows_, this->globalNumNonZeros_,
      //                             colind_, rowptr_, nzvals_, data_.perm_r.getRawPtr(),
      //                             data_.R1.getRawPtr(), data_.C1.getRawPtr());

      //        // Note: if equilibration was not requested, then R1 and C1
      //        // scalings will never be applied to the matrix values
      //        // (indeed, the arrays will never be accessed)

      // }      // else it must be MY_PERMR, and we expect perm_r to have been initialized

      // // Finally, find a column permutation if requested
      // // Use either the column-ordering found in the users perm_c or
      // // the requested computed ordering
      // int perm_spec = data_.options.ColPerm;
      // if( perm_spec != SLUD::MY_PERMC ){
      //        {                           // start matrix conversion block
      //          Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);

      //          if( data_.options.RowPerm != SLUD::LargeDiag ){
      //            matrix_helper::createCcsMatrix(this->matrixA_.ptr(),
      //                                           nzvals_(), colind_(), rowptr_(),
      //                                           Teuchos::ptrFromRef(data_.A),
      //                                           this->timers_.mtxRedistTime_);
      //          } // else we've already created the supermatrix
      //        } // end matrix conversion block

      //        TEST_FOR_EXCEPTION( perm_spec == SLUD::PARMETIS,
      //                            std::invalid_argument,
      //                            "Parmetis orderings not yet supported in Amesos2" );

      //        // TODO: Do we need to apply R1 and C1 to A before finding a
      //        // column permutation?
      //        SLUD::get_perm_c_dist(this->status_.myPID,
      //                              perm_spec,
      //                              &(data_.A),
      //                              data_.perm_c.getRawPtr());
      // }
      // // Else the user's perm_c will be applied later
    } else {                    // ParSymbFact == YES (default, required)

      loadA();			// Refresh matrix values

      if( in_grid_ ){
	// If this function has been called at least once, then the
	// sizes, and fstVtxSep arrays were allocated in
	// get_perm_c_parmetis.  Delete them before calling that
	// function again.  These arrays will also be dealloc'd in the
	// deconstructor.
	if( this->getNumPreOrder() > 0 ){
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
	
	TEST_FOR_EXCEPTION( info > 0.0,
			    std::runtime_error,
			    "SuperLU_DIST pre-ordering ran out of memory after allocating "
			    << info << " bytes of memory" );
      }
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
    loadA();			// Refresh matrix values
    
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
      
      TEST_FOR_EXCEPTION( info > 0.0,
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

    loadA();			// Refresh the matrix values

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
      TEST_FOR_EXCEPTION( info > 0,
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
  Superludist<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
					 const Teuchos::Ptr<MultiVecAdapter<Vector> > B) const
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
      
      TEST_FOR_EXCEPTION( ierr < 0,
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
      {				// permute B to a solution of the original system
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

    if( parameterList->isParameter("npcol") || parameterList->isParameter("nprow") ){
      TEST_FOR_EXCEPTION( parameterList->isParameter("nprow") &&
			  parameterList->isParameter("npcol"),
			  std::invalid_argument,
			  "nprow and npcol must be set together" );
      SLUD::int_t nprow = parameterList->template get<SLUD::int_t>("nprow");
      SLUD::int_t npcol = parameterList->template get<SLUD::int_t>("npcol");

      TEST_FOR_EXCEPTION( nprow * npcol > this->getComm()->getSize(),
			  std::invalid_argument,
			  "nprow and npcol combination invalid" );

      // De-allocate the default grid that was initialized in the constructor
      SLUD::superlu_gridexit(&(data_.grid));
      // Create a new grid
      SLUD::superlu_gridinit(data_.mat_comm, nprow, npcol, &(data_.grid));
    }

    TEST_FOR_EXCEPTION( this->control_.useTranspose_,
			std::invalid_argument,
			"SuperLU_DIST does not support solving the tranpose system" );

    if ( parameterList->isParameter("Trans") ){
      std::string fact = parameterList->template get<std::string>("Trans");
      TEST_FOR_EXCEPTION( fact == "TRANS" || fact == "CONJ",
			  std::invalid_argument,
			  "SuperLU_DIST does not support solving the transpose system" );
    }
    // The 'fact' parameter is in the options struct, but it seems to
    // be more for adherance to the same interface as the other
    // SuperLU's.  The documentation and the code suggest that the
    // option is not supported.
    data_.options.Trans = SLUD::NOTRANS; // Should always be set this way


    if( parameterList->isParameter("Equil") ){
      if ( parameterList->template isType<bool>("Equil") ){
	bool equil = parameterList->template get<bool>("Equil");
	if( equil ){
	  data_.options.Equil = SLUD::YES;
	} else {
	  data_.options.Equil = SLUD::NO;
	}
      } else if ( parameterList->template isType<std::string>("Equil") ) {
	std::string equil = parameterList->template get<std::string>("Equil");
	if ( equil == "YES" || equil == "yes" ){
	  data_.options.Equil = SLUD::YES;
	} else if ( equil == "NO" || equil == "no" ) {
	  data_.options.Equil = SLUD::NO;
	}
      }
    }

    if( parameterList->isParameter("ParSymbFact") ){
      if ( parameterList->template isType<bool>("ParSymbFact") ){
	bool parsymbfact = parameterList->template get<bool>("ParSymbFact");
	if( parsymbfact ){
	  data_.options.ParSymbFact = SLUD::YES;
	} else {
	  data_.options.ParSymbFact = SLUD::NO;
	}
      } else if ( parameterList->template isType<std::string>("ParSymbFact") ) {
	std::string parsymbfact = parameterList->template get<std::string>("ParSymbFact");
	if ( parsymbfact == "YES" || parsymbfact == "yes" ){
	  data_.options.ParSymbFact = SLUD::YES;
	} else if ( parsymbfact == "NO" || parsymbfact == "no" ) {
	  data_.options.ParSymbFact = SLUD::NO;
	}
      }

      TEST_FOR_EXCEPTION( data_.options.ParSymbFact == SLUD::NO,
			  std::invalid_argument,
			  "Amesos::Superludist requires parallel symbolic factorization" );
    }


    if( parameterList->isParameter("SymmetricMode") ){
      if ( parameterList->template isType<bool>("SymmetricMode") ){
	bool sym = parameterList->template get<bool>("SymmetricMode");
	if( sym ){
	  data_.options.SymmetricMode = SLUD::YES;
	} else {
	  data_.options.SymmetricMode = SLUD::NO;
	}
      } else if ( parameterList->template isType<std::string>("SymmetricMode") ) {
	std::string sym = parameterList->template get<std::string>("SymmetricMode");
	if ( sym == "YES" || sym == "yes" ){
	  data_.options.SymmetricMode = SLUD::YES;
	} else if ( sym == "NO" || sym == "no" ) {
	  data_.options.SymmetricMode = SLUD::NO;
	}
      }
    }

    if( parameterList->isParameter("PrintStat") ){
      if ( parameterList->template isType<bool>("PrintStat") ){
	bool sym = parameterList->template get<bool>("PrintStat");
	if( sym ){
	  data_.options.PrintStat = SLUD::YES;
	} else {
	  data_.options.PrintStat = SLUD::NO;
	}
      } else if ( parameterList->template isType<std::string>("PrintStat") ) {
	std::string sym = parameterList->template get<std::string>("PrintStat");
	if ( sym == "YES" || sym == "yes" ){
	  data_.options.PrintStat = SLUD::YES;
	} else if ( sym == "NO" || sym == "no" ) {
	  data_.options.PrintStat = SLUD::NO;
	}
      }
    }

    if( parameterList->isParameter("ColPerm") ){
      std::string col_perm_method = parameterList->template get<std::string>("ColPerm");

      TEST_FOR_EXCEPTION( col_perm_method != "PARMETIS" && col_perm_method != "NATURAL",
			  std::invalid_argument,
			  "Amesos::Superludist only supports natural column orderings "
			  "or parallel ordering with ParMETIS" );

      if( col_perm_method == "NATURAL" ){
	data_.options.ColPerm = SLUD::NATURAL;
      } else if ( col_perm_method == "MMD_AT_PLUS_A" ) {
	data_.options.ColPerm = SLUD::MMD_AT_PLUS_A;
      } else if ( col_perm_method == "MMD_ATA" ) {
	data_.options.ColPerm = SLUD::MMD_ATA;
      } else if ( col_perm_method == "COLAMD" ) {
	data_.options.ColPerm = SLUD::COLAMD;
      } else if ( col_perm_method == "METIS_AT_PLUS_A" ) {
	data_.options.ColPerm = SLUD::METIS_AT_PLUS_A;
      } else if ( col_perm_method == "PARMETIS" ) {
	data_.options.ColPerm = SLUD::PARMETIS;
      } else if ( col_perm_method == "MY_PERMC" ) {
	data_.options.ColPerm = SLUD::MY_PERMC;

	// Now we also expect to find a parameter in parameterList called
	// "perm_c"
	TEST_FOR_EXCEPTION(
			   !parameterList->isParameter("perm_c"),
			   std::invalid_argument,
			   "MY_PERMC option specified without accompanying 'perm_c' parameter.");

	data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

	TEST_FOR_EXCEPTION(
			   as<global_size_type>(data_.perm_c.size()) != this->globalNumCols_,
			   std::length_error,
			   "'perm_c' parameter not of correct length.");
      } else {
	TEST_FOR_EXCEPTION( true, std::invalid_argument,
			    "Unrecognized value for 'ColPerm' parameter.");
      }
    }

    // We also recognize a lone 'perm_c' parameter, assuming that ColPerm = MY_PERMC
    if( parameterList->isParameter("perm_c") ){
      data_.options.ColPerm = SLUD::MY_PERMC;
      data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

      TEST_FOR_EXCEPTION(
			 as<global_size_type>(data_.perm_c.size()) == this->globalNumCols_,
			 std::length_error,
			 "'perm_c' parameter not of correct length.");
    }

    if( parameterList->isParameter("IterRefine") ){
      std::string refine = parameterList->template get<std::string>("IterRefine");
      if( refine == "NO" ){
	data_.options.IterRefine = SLUD::NOREFINE;
      } else if ( refine == "SINGLE" ) {
	data_.options.IterRefine = SLUD::SINGLE;
      } else if ( refine == "DOUBLE" ) {
	data_.options.IterRefine = SLUD::DOUBLE;
      } else if ( refine == "EXTRA" ) {
	data_.options.IterRefine = SLUD::EXTRA;
      } else {
	TEST_FOR_EXCEPTION(
			   true,
			   std::invalid_argument,
			   "Unrecognized value for 'IterRefine' parameter.");
      }
    }

    if( parameterList->isParameter("RowPerm") ){
      std::string method = parameterList->template get<std::string>("RowPerm");
      if( method == "NO" || method == "NATURAL" ){
	data_.options.RowPerm = SLUD::NOROWPERM;
      } else if( method == "LargeDiag" ){
	data_.options.RowPerm = SLUD::LargeDiag;
      } else if( method == "MY_PERMR" ){
	data_.options.RowPerm = SLUD::MY_PERMR;

	// Expect to find a parameter called "perm_r"
	TEST_FOR_EXCEPTION( !parameterList->isParameter("perm_r"),
			    std::invalid_argument,
			    "MY_PERMR specified without accompanying 'perm_r' parameter." );

	data_.perm_r = parameterList->template get<Teuchos::Array<int> >("perm_r");
	TEST_FOR_EXCEPTION( as<global_size_type>(data_.perm_r.size()) != this->globalNumRows_,
			    std::invalid_argument,
			    "'perm_r' parameter not of correct length." );
      } else {
	TEST_FOR_EXCEPTION( true, std::invalid_argument,
			    "Unrecognized value for 'RowPerm' parameter." );
      }
      TEST_FOR_EXCEPTION( data_.options.RowPerm != SLUD::NOROWPERM,
			  std::invalid_argument,
			  "Amesos::Superludist does not support row permutations" );
    }

    if( parameterList->isParameter("ReplaceTinyPivot") ){
      if( parameterList->template isType<bool>("ReplaceTinyPivot") ){
	bool replace = parameterList->template get<bool>("ReplaceTinyPivot");
	if( replace ){
	  data_.options.ReplaceTinyPivot = SLUD::YES;
	} else {
	  data_.options.ReplaceTinyPivot = SLUD::NO;
	}
      } else if ( parameterList->template isType<std::string>("ReplaceTinyPivot") ){
	std::string replace = parameterList->template get<std::string>("ReplaceTinyPivot");
	if( replace == "YES" || replace == "yes" ){
	  data_.options.ReplaceTinyPivot = SLUD::YES;
	} else if( replace == "NO" || replace == "no" ){
	  data_.options.ReplaceTinyPivot = SLUD::NO;
	}
      }
    }
  }


  /** \internal
   *
   * Using the "NO" (don't do any row permutation) option for RowPerm
   * eliminates a sequential bottleneck in the SuperLU_DIST code.  If
   * the other option ("LargeDiag") is used, then SuperLU_DIST must
   * gather the matrix to processor 0 in order to perform a sequential
   * bipartite matching algorithm.  Disregard the SuperLU_DIST
   * documentation that says that "NATURAL" is a valid option.
   */
  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  Superludist<Matrix,Vector>::getValidParameters_impl() const
  {
    using Teuchos::ParameterList;

    Teuchos::RCP<ParameterList> valid_params = Teuchos::rcp(new ParameterList());

    valid_params->set("Equil", false); // TODO: change to true when supported
    valid_params->set("ParSymbFact", true);
    valid_params->set("ReplaceTinyPivot", true);
    valid_params->set("ColPerm", "PARMETIS");
    valid_params->set("RowPerm", "NO");
    valid_params->set("IterRefine", "DOUBLE");
    valid_params->set("Trans", "NOTRANS");
    valid_params->set("PrintStat", false);

    return valid_params;
  }


  template <class Matrix, class Vector>
  void
  Superludist<Matrix,Vector>::get_default_grid_size(int nprocs,
						    SLUD::int_t& nprow,
						    SLUD::int_t& npcol) const {
    TEST_FOR_EXCEPTION( nprocs < 1,
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
  void
  Superludist<Matrix,Vector>::loadA(){
    // Extract the necessary information from mat and call SLU function
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ptrInArg;

#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Cleanup old store memory if it's non-NULL
    if( data_.A.Store != NULL ){
      SLUD::Destroy_SuperMatrix_Store_dist( &(data_.A) );
    }
    
    Teuchos::RCP<const MatrixAdapter<Matrix> > redist_mat
      = this->matrixA_->get(ptrInArg(*superlu_rowmap_));
    
    SLUD::int_t l_nnz, l_rows, g_rows, g_cols, fst_global_row;
    l_nnz  = Teuchos::as<SLUD::int_t>(redist_mat->getLocalNNZ());
    l_rows = Teuchos::as<SLUD::int_t>(redist_mat->getLocalNumRows());
    g_rows = Teuchos::as<SLUD::int_t>(redist_mat->getGlobalNumRows());
    // g_cols = Teuchos::as<SLUD::int_t>(redist_mat->getGlobalNumCols());
    g_cols = g_rows;		// should be a square matrix anyhow
    fst_global_row = Teuchos::as<SLUD::int_t>(superlu_rowmap_->getMinGlobalIndex());

    nzvals_.resize(l_nnz);
    colind_.resize(l_nnz);
    rowptr_.resize(l_rows + 1);
    
    SLUD::int_t nnz_ret = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif
      
      Util::get_crs_helper<
        MatrixAdapter<Matrix>,
	slu_type,
	SLUD::int_t,
	SLUD::int_t >::do_get(redist_mat.ptr(), nzvals_(), colind_(),
			      rowptr_(), nnz_ret, ptrInArg(*superlu_rowmap_),
			      Util::Arbitrary);
    }
    
    TEST_FOR_EXCEPTION( nnz_ret != l_nnz,
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
  }


  template<class Matrix, class Vector>
  const char* Superludist<Matrix,Vector>::name = "SuperLU_DIST";


} // end namespace Amesos

#endif  // AMESOS2_SUPERLUDIST_DEF_HPP
