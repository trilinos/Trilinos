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
    , nzvals_()                   // initialize to empty arrays
    , rowind_()
    , colptr_()
    , is_contiguous_(true)
  {

    typedef FunctionMap<MUMPS,scalar_type>  function_map;
    
    MUMPS_MATRIX_LOAD = false;
    MUMPS_STRUCT = false;
    
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
    if(MUMPS_STRUCT == true)
      {
        free(mumps_par.a);
        free(mumps_par.jcn);
        free(mumps_par.irn);
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
    
    typedef FunctionMap<MUMPS,scalar_type> function_map;
    
    mumps_par.par = 1;  
    mumps_par.job = 1; // sym factor
    function_map::mumps_c(&(mumps_par));
    MUMPS_ERROR();
    
    return(0);
  }//end symblicFactortion_impl()
  
  
  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::numericFactorization_impl()
  {
    using Teuchos::as;
    typedef FunctionMap<MUMPS,scalar_type> function_map;
    
    if ( this->root_ )
      {
        { // Do factorization
          #ifdef HAVE_AMESOS2_TIMERS
          Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
          #endif
          
          #ifdef HAVE_AMESOS2_VERBOSE_DEBUG
          std::cout << "MUMPS:: Before numeric factorization" << std::endl;
          std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
          std::cout << "rowind_ : " << rowind_.toString() << std::endl;
          std::cout << "colptr_ : " << colptr_.toString() << std::endl;
          #endif
        }
      }
    mumps_par.job = 2;
    function_map::mumps_c(&(mumps_par));
    MUMPS_ERROR();
    
    return(0);
  }//end numericFactorization_impl()


  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::solve_impl(
                         const Teuchos::Ptr<MultiVecAdapter<Vector> >  X,
                         const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    
    typedef FunctionMap<MUMPS,scalar_type> function_map;
    
    using Teuchos::as;
    
    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    const size_t nrhs = X->getGlobalNumVectors();
    
    const size_t val_store_size = as<size_t>(ld_rhs * nrhs);
    
    xvals_.resize(val_store_size);
    bvals_.resize(val_store_size);
    
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
    #endif

    if ( is_contiguous_ == true ) {
      Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
        slu_type>::do_get(B, bvals_(),as<size_t>(ld_rhs), ROOTED, this->rowIndexBase_);
    }
    else {
      Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
        slu_type>::do_get(B, bvals_(),as<size_t>(ld_rhs), CONTIGUOUS_AND_ROOTED, this->rowIndexBase_);
    }
 
    
    int ierr = 0; // returned error code
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
    
    
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
    #endif

    if ( is_contiguous_ == true ) {
      Util::put_1d_data_helper<
        MultiVecAdapter<Vector>,slu_type>::do_put(X, bvals_(),
            as<size_t>(ld_rhs),
            ROOTED);
    }
    else {
      Util::put_1d_data_helper<
        MultiVecAdapter<Vector>,slu_type>::do_put(X, bvals_(),
            as<size_t>(ld_rhs),
            CONTIGUOUS_AND_ROOTED);
    }
    
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
    
    if(parameterList->isParameter("ICNTL(1)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(1)");
      }
    if(parameterList->isParameter("ICNTL(2)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(2)");
      }
    if(parameterList->isParameter("ICNTL(3)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(3)");
      }
    if(parameterList->isParameter("ICNTL(4)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(4)");
      }
    if(parameterList->isParameter("ICNTL(6)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(6)");
      }
    if(parameterList->isParameter("ICNTL(9)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(9)");
      }
    if(parameterList->isParameter("ICNTL(11)"))
      {
        mumps_par.icntl[0] = getIntegralValue<local_ordinal_type>(*parameterList, 
                                                                  "ICNTL(11)");
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
    
      pl->set("ICNTL(1)", "no", "See Manual" );
      pl->set("ICNTL(2)", "no", "See Manual" );
      pl->set("ICNTL(3)", "no", "See Manual" );
      pl->set("ICNTL(4)", "no", "See Manual" );
      pl->set("ICNTL(6)", "no", "See Manual" );
      pl->set("ICNTL(9)", "no", "See Manual" );
      pl->set("ICNTL(11)", "no", "See Manual" );

      pl->set("IsContiguous", true, "Whether GIDs contiguous");
      
      valid_params = pl;
    }
    
    return valid_params;
  }//end getValidParmaeters_impl()
  
  
  template <class Matrix, class Vector>
  bool
  MUMPS<Matrix,Vector>::loadA_impl(EPhase current_phase)
  {
    using Teuchos::as;
    
    #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
    #endif
    
    if(MUMPS_MATRIX_LOAD == false)
      {
        // Only the root image needs storage allocated
        if( this->root_ ){
          nzvals_.resize(this->globalNumNonZeros_);
          rowind_.resize(this->globalNumNonZeros_);
          colptr_.resize(this->globalNumCols_ + 1);
        }
  
        local_ordinal_type nnz_ret = 0;
        
        #ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
        #endif
    
        if ( is_contiguous_ == true ) {
          Util::get_ccs_helper<
            MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
            ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
                nnz_ret, ROOTED, ARBITRARY, this->rowIndexBase_);
        }
        else {
          Util::get_ccs_helper<
            MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
            ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
                nnz_ret, CONTIGUOUS_AND_ROOTED, ARBITRARY, this->rowIndexBase_);
        }
  
        if( this->root_ ){
                  TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
                                std::runtime_error,
                        "Did not get the expected number of non-zero vals");
        }
  
        
        if( this->root_ ){
          ConvertToTriplet();
        }
      }
    
    MUMPS_MATRIX_LOAD = true;
    return (true);
  }//end loadA_impl()
  
  template <class Matrix, class Vector>
  int
  MUMPS<Matrix,Vector>::ConvertToTriplet()
  {
    MUMPS_STRUCT = true;
    mumps_par.n =  this->globalNumCols_;
    mumps_par.nz = this->globalNumNonZeros_;
    mumps_par.a = (magnitude_type*)malloc(mumps_par.nz * sizeof(magnitude_type));
    mumps_par.irn = (MUMPS_INT*)malloc(mumps_par.nz *sizeof(MUMPS_INT));
    mumps_par.jcn = (MUMPS_INT*)malloc(mumps_par.nz * sizeof(MUMPS_INT));

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
        for( j = colptr_[i]; j < colptr_[i+1]-1; j++)
          {
            mumps_par.jcn[tri_count] = (MUMPS_INT)i+1; //Fortran index
            mumps_par.irn[tri_count] = (MUMPS_INT)rowind_[j]+1; //Fortran index
            mumps_par.a[tri_count] = nzvals_[j];
            
            tri_count++;
          }
        
        j = colptr_[i+1]-1;
        mumps_par.jcn[tri_count] = (MUMPS_INT)i+1; //Fortran index
        mumps_par.irn[tri_count] = (MUMPS_INT)rowind_[j]+1; //Fortran index
        mumps_par.a[tri_count] = nzvals_[j];

        tri_count++;
        
        if(rowind_[j] > max_local_ordinal)
          {
            max_local_ordinal = rowind_[j];
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
    if(mumps_par.info[0] < 0)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(false,
                                   std::runtime_error,
                                   "MUMPS error");
      }
  }//end MUMPS_ERROR()


  template<class Matrix, class Vector>
  const char* MUMPS<Matrix,Vector>::name = "MUMPS";
  
} // end namespace Amesos2

#endif  // AMESOS2_MUMPS_DEF_HPP
