// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_MUMPS_def.hpp
   \author Joshua Dennis Booth <jdbooth21@sandia.gov>

   \brief  Definitions for the Amesos2 MUMPS solver interface
*/


#ifndef AMESOS2_MUMPS_DEF_HPP
#define AMESOS2_MUMPS_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#endif

#include <limits>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_MUMPS_decl.hpp"

namespace Amesos2 
{


  template <class Matrix, class Vector>
  MUMPS<Matrix,Vector>::MUMPS(
                              Teuchos::RCP<const Matrix> A,
                              Teuchos::RCP<Vector>       X,
                              Teuchos::RCP<const Vector> B )
    : SolverCore<Amesos2::MUMPS,Matrix,Vector>(A, X, B)
    , is_contiguous_(true)
  {

    typedef FunctionMap<MUMPS,scalar_type>  function_map;
    
    MUMPS_MATRIX_LOAD = false;
    MUMPS_STRUCT = false;
    MUMPS_MATRIX_LOAD_PREORDERING = false;
      
    #ifdef HAVE_MPI
    using Teuchos::Comm;
    using Teuchos::MpiComm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
  
    //Comm
    mumps_par.comm_fortran = -987654;
    RCP<const Comm<int> > matComm = this->matrixA_->getComm();
    //Add Exception Checking
    TEUCHOS_TEST_FOR_EXCEPTION(
                               matComm.is_null(), std::logic_error, "Amesos2::Comm");
    RCP<const MpiComm<int> > matMpiComm = 
      rcp_dynamic_cast<const MpiComm<int> >(matComm);
    //Add Exception Checking
    TEUCHOS_TEST_FOR_EXCEPTION(
                               matMpiComm->getRawMpiComm().is_null(), 
                               std::logic_error, "Amesos2::MPI");
    MPI_Comm rawMpiComm = (* (matMpiComm->getRawMpiComm()) )();
    mumps_par.comm_fortran = (int) MPI_Comm_c2f(rawMpiComm);
    #endif  
  
    mumps_par.job = -1;
    mumps_par.par = 1;
    mumps_par.sym = 0;
    
    function_map::mumps_c(&(mumps_par)); //init
    MUMPS_ERROR();
    
    mumps_par.n = this->globalNumCols_;

    /*Default parameters*/
    mumps_par.icntl[0] = -1; // Turn off error messages
    mumps_par.icntl[1] = -1; // Turn off diagnositic printing
    mumps_par.icntl[2] = -1; // Turn off global information messages
    mumps_par.icntl[3] = 1; // No messages
    mumps_par.icntl[4] = 0;  // Matrix given in assembled (Triplet) form
    mumps_par.icntl[5] = 7;  // Choose column permuation automatically
    mumps_par.icntl[6] = 7;  // Choose ordering methods automatically
    mumps_par.icntl[7] = 7; // Choose scaling automatically
    mumps_par.icntl[8] = 1;  // Compuate Ax = b;
    mumps_par.icntl[9] = 0;  // iterative refinement
    mumps_par.icntl[10] = 0; //Do not collect statistics
    mumps_par.icntl[11] = 0; // Automatic choice of ordering strategy
    mumps_par.icntl[12] = 0; //Use ScaLAPACK for root node
    mumps_par.icntl[13] = 20; // Increase memory allocation 20% at a time
    mumps_par.icntl[17] = 0; 
    mumps_par.icntl[18] = 0; // do not provide back the Schur complement
    mumps_par.icntl[19] = 0; //  RHS is given in dense form
    mumps_par.icntl[20] = 0; //RHS is in dense form
    mumps_par.icntl[21] = 0; // Do all compuations in-core
    mumps_par.icntl[22] = 0; // No max MB for work
    mumps_par.icntl[23] = 0; // Do not perform null pivot detection
    mumps_par.icntl[24] = 0; // No null space basis compuation
    mumps_par.icntl[25] = 0; // Do not condense/reduce Schur RHS
    mumps_par.icntl[27] = 1; // sequential analysis
    mumps_par.icntl[28] = 0; // 
    mumps_par.icntl[29] = 0; //
    mumps_par.icntl[30] = 0; //
    mumps_par.icntl[31] = 0; 
    mumps_par.icntl[32] = 0; 
  }
  
  template <class Matrix, class Vector>
  MUMPS<Matrix,Vector>::~MUMPS( )
  {
    /* Clean up the struc*/
    typedef FunctionMap<MUMPS,scalar_type> function_map;

    if(MUMPS_STRUCT == true)
      {
        free(mumps_par.a);
        free(mumps_par.jcn);
        free(mumps_par.irn);
      }
      mumps_par.job = -2;
      if (this->rank_ < this->nprocs_) {
          function_map::mumps_c(&(mumps_par));
      }

  }
  
  template<class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::preOrdering_impl()
  {
    /* TODO: Define what it means for MUMPS
     */
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
    #endif
    
    return(0);
  }//end preOrdering_impl()

  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::symbolicFactorization_impl()
  {
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
    #endif  

    typedef FunctionMap<MUMPS,scalar_type> function_map;
    if ( this->globalNumRows_ > 0 ) {
      mumps_par.par = 1;
      mumps_par.job = 1; // sym factor
      function_map::mumps_c(&(mumps_par));
      MUMPS_ERROR();
    }
    return(0);
  }//end symblicFactortion_impl()
  
  
  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::numericFactorization_impl()
  {
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
    #endif

    typedef FunctionMap<MUMPS,scalar_type> function_map;
    if ( this->globalNumRows_ > 0 ) {
      mumps_par.job = 2;
      function_map::mumps_c(&(mumps_par));
      MUMPS_ERROR();
    }
    return(0);
  }//end numericFactorization_impl()


  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::solve_impl(
                         const Teuchos::Ptr<MultiVecAdapter<Vector> >  X,
                         const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    typedef FunctionMap<MUMPS,scalar_type> function_map;

    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    const size_t nrhs = X->getGlobalNumVectors();
    const size_t val_store_size = Teuchos::as<size_t>(ld_rhs * nrhs);

    xvals_.resize(val_store_size);
    bvals_.resize(val_store_size);

    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
    #endif

    Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
      mumps_type>::do_get(B, bvals_(), Teuchos::as<size_t>(ld_rhs),
        (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
        this->rowIndexBase_);
 
    int ierr = 0; // returned error code
    if ( this->globalNumRows_ > 0 ) {
      mumps_par.nrhs = nrhs;
      mumps_par.lrhs = mumps_par.n;
      mumps_par.job = 3;

      if ( this->root_ ) 
        {
          mumps_par.rhs = bvals_.getRawPtr();
        }

      #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
      #endif

      function_map::mumps_c(&(mumps_par));
      MUMPS_ERROR();
    }

    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer2(this->timers_.vecRedistTime_);
    #endif

    Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,mumps_type>::do_put(X, bvals_(), Teuchos::as<size_t>(ld_rhs),
        (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED);

    // ch: see function loadA_impl()
    MUMPS_MATRIX_LOAD_PREORDERING = false;
    return(ierr);
  }//end solve()
  
  
  template <class Matrix, class Vector>
  bool
  MUMPS<Matrix,Vector>::matrixShapeOK_impl() const
  {
    // The MUMPS can only handle square for right now
    return( this->globalNumRows_ == this->globalNumCols_ );
  }
  
  
  template <class Matrix, class Vector>
  void
  MUMPS<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    using Teuchos::RCP;
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterEntryValidator;
    
    RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();
    /*To Do --- add support for parameters */
    if(parameterList->isParameter("ICNTL(1)")){
        mumps_par.icntl[0] = parameterList->get<int>("ICNTL(1)", -1);
    }
    if(parameterList->isParameter("ICNTL(2)")){
        mumps_par.icntl[1] = parameterList->get<int>("ICNTL(2)", -1);
    }
    if(parameterList->isParameter("ICNTL(3)")){
        mumps_par.icntl[2] = parameterList->get<int>("ICNTL(3)", -1);
    }
    if(parameterList->isParameter("ICNTL(4)")){
        mumps_par.icntl[3] = parameterList->get<int>("ICNTL(4)", 1);
    }
    if(parameterList->isParameter("ICNTL(6)")){
        mumps_par.icntl[5] = parameterList->get<int>("ICNTL(6)", 0);
    }
    if(parameterList->isParameter("ICNTL(7)")){
        mumps_par.icntl[6] = parameterList->get<int>("ICNTL(7)", 7);
    }
    if(parameterList->isParameter("ICNTL(9)")){
        mumps_par.icntl[8] = parameterList->get<int>("ICNTL(9)", 1);
    }
    if(parameterList->isParameter("ICNTL(11)")){
        mumps_par.icntl[10] = parameterList->get<int>("ICNTL(11)", 0);
    }
    if(parameterList->isParameter("ICNTL(14)")){
          mumps_par.icntl[13] = parameterList->get<int>("ICNTL(14)", 20);
    }
    if( parameterList->isParameter("IsContiguous") ){
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }
  }//end set parameters()

  
  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  MUMPS<Matrix,Vector>::getValidParameters_impl() const
  {
    using Teuchos::ParameterList;
    
    static Teuchos::RCP<const Teuchos::ParameterList> valid_params;
    
    if( is_null(valid_params) ){
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    
       pl->set("ICNTL(1)", -1, "See Manual" );
       pl->set("ICNTL(2)", -1, "See Manual" );
       pl->set("ICNTL(3)", -1, "See Manual" );
       pl->set("ICNTL(4)", 1, "See Manual" );
       pl->set("ICNTL(6)", 0, "See Manual" );
       pl->set("ICNTL(9)", 1, "See Manual" );
       pl->set("ICNTL(11)", 0, "See Manual" );
       pl->set("ICNTL(14)", 20, "See Manual" );
       pl->set("IsContiguous", true, "Whether GIDs contiguous");
      
       valid_params = pl;
    }
    
    return valid_params;
  }//end getValidParmaeters_impl()
  
  
  template <class Matrix, class Vector>
  bool
  MUMPS<Matrix,Vector>::loadA_impl(EPhase current_phase)
  {
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
    #endif
    if(MUMPS_MATRIX_LOAD == false || current_phase==NUMFACT)
      {
        // Only the root image needs storage allocated
        if( !MUMPS_MATRIX_LOAD && this->root_ ){
          Kokkos::resize(host_nzvals_view_, this->globalNumNonZeros_);
          Kokkos::resize(host_rows_view_, this->globalNumNonZeros_);
          Kokkos::resize(host_col_ptr_view_, this->globalNumRows_ + 1);
        }
  
        local_ordinal_type nnz_ret = 0;
        
        #ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
        #endif

        Util::get_ccs_helper_kokkos_view<
          MatrixAdapter<Matrix>,host_value_type_view,host_ordinal_type_view,host_ordinal_type_view>
          ::do_get(this->matrixA_.ptr(), host_nzvals_view_, host_rows_view_, host_col_ptr_view_, nnz_ret,
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            ARBITRARY,
            this->rowIndexBase_);
  
        if( this->root_ ){
                  TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != Teuchos::as<local_ordinal_type>(this->globalNumNonZeros_),
                                std::runtime_error,
                        "Did not get the expected number of non-zero vals");
        }
  
        
        if( this->root_ ){
          ConvertToTriplet();
        }
        /* ch: In general, the matrix is loaded during the preordering phase.
           However, if the matrix pattern has not changed during consecutive calls of numeric factorizations
           we can reuse the previous symbolic factorization. In this case, the matrix is not loaded in the preordering phase,
           because it is not called. Therefore, we need to load the matrix in the numeric factorization phase.
         */
          if (current_phase==PREORDERING){
              MUMPS_MATRIX_LOAD_PREORDERING = true;
          }
      }
      
      
    
    MUMPS_MATRIX_LOAD = true;
    return (true);
  }//end loadA_impl()
  
  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::ConvertToTriplet()
  {
    if ( !MUMPS_STRUCT ) {
      MUMPS_STRUCT = true;
      mumps_par.n =  this->globalNumCols_;
      mumps_par.nz = this->globalNumNonZeros_;
      mumps_par.a = (mumps_type*)malloc(mumps_par.nz * sizeof(mumps_type));
      mumps_par.irn = (MUMPS_INT*)malloc(mumps_par.nz *sizeof(MUMPS_INT));
      mumps_par.jcn = (MUMPS_INT*)malloc(mumps_par.nz * sizeof(MUMPS_INT));
    }
    if((mumps_par.a == NULL) || (mumps_par.irn == NULL) 
       || (mumps_par.jcn == NULL))
      {
        return -1;
      }
    /* Going from full CSC to full Triplet */
    /* Will have to add support for symmetric case*/
    local_ordinal_type tri_count = 0;
    local_ordinal_type i,j;
    local_ordinal_type max_local_ordinal = 0;
    
    for(i = 0; i < (local_ordinal_type)this->globalNumCols_; i++)
      {
        for( j = host_col_ptr_view_(i); j < host_col_ptr_view_(i+1)-1; j++)
          {
            mumps_par.jcn[tri_count] = (MUMPS_INT)i+1; //Fortran index
            mumps_par.irn[tri_count] = (MUMPS_INT)host_rows_view_(j)+1; //Fortran index
            mumps_par.a[tri_count] = host_nzvals_view_(j);
            
            tri_count++;
          }
        
        j = host_col_ptr_view_(i+1)-1;
        mumps_par.jcn[tri_count] = (MUMPS_INT)i+1; //Fortran index
        mumps_par.irn[tri_count] = (MUMPS_INT)host_rows_view_(j)+1; //Fortran index
        mumps_par.a[tri_count] = host_nzvals_view_(j);

        tri_count++;
        
        if(host_rows_view_(j) > max_local_ordinal)
          {
            max_local_ordinal = host_rows_view_(j);
          }
      }
    TEUCHOS_TEST_FOR_EXCEPTION(std::numeric_limits<MUMPS_INT>::max() <= max_local_ordinal,
                               std::runtime_error,
                               "Matrix index larger than MUMPS_INT");

    return 0;
  }//end Convert to Trip()

  template<class Matrix, class Vector>
  void
  MUMPS<Matrix,Vector>::MUMPS_ERROR()const
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    bool Wrong = ((mumps_par.info[0] != 0) || (mumps_par.infog[0] != 0)) && (this->rank_ < this->nprocs_);
    if(Wrong){
      if (this->rank_==0) {
        std::cerr << "Amesos2_Mumps : ERROR" << std::endl;
        if ( this->status_.getNumSolve() > 0) {
            std::cerr << " Last Phase : SOLVE" << std::endl;
        } else if( this->status_.numericFactorizationDone() ){
          std::cerr << " Last Phase : NUMFACT" << std::endl;
        } else if( this->status_.symbolicFactorizationDone() ){
          std::cerr << " Last Phase : SYMBFACT" << std::endl;
        } else if( this->status_.preOrderingDone() ){
          std::cerr << " Last Phase : PREORDERING" << std::endl;
        } else {
          std::cerr << " Last Phase : CLEAN" << std::endl;
        }
        std::cerr << "Amesos2_Mumps : INFOG(1) = " << mumps_par.infog[0] << std::endl;
        std::cerr << "Amesos2_Mumps : INFOG(2) = " << mumps_par.infog[1] << std::endl;
      }
      if (mumps_par.info[0] != 0  && Wrong) {
        std::cerr << "Amesos2_Mumps : On process " << this->matrixA_->getComm()->getRank()
        << ", INFO(1) = " << mumps_par.info[0] << std::endl;
        std::cerr << "Amesos2_Mumps : On process " << this->matrixA_->getComm()->getRank()
        << ", INFO(2) = " << mumps_par.info[1] << std::endl;
      }
      
      
    }
    // Throw on all ranks
    int WrongInt = Wrong;
    RCP<const Comm<int> > matComm = this->matrixA_->getComm();
    Teuchos::broadcast<int,int>(*matComm,0,1,&WrongInt);
    TEUCHOS_TEST_FOR_EXCEPTION(WrongInt>0,
                               std::runtime_error,
                               "MUMPS error");
    
  }//end MUMPS_ERROR()


  template<class Matrix, class Vector>
  const char* MUMPS<Matrix,Vector>::name = "MUMPS";
  
} // end namespace Amesos2

#endif  // AMESOS2_MUMPS_DEF_HPP
