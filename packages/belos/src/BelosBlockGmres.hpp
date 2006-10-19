// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This file contains an implementation of the Block GMRES algorithm
// for solving real nonsymmetric linear systems of equations AX = B,
// where B is a matrix containing one or more right-hand sides, and X is
// the matrix of corresponding solutions. This implementation allows the 
// user to solve systems involving any number of right-hand sides. The
// block size used in the solver is user specified, and is independent
// of the number of right-hand sides. Thus, a system involving many 
// right-hand sides can be processed by solving for only some number  
// (the block size) of the right-hand sides simultaneously. Several passes
// through the block solver are used to solve for all of them. A single
// right-hand side system can be solved in the traditional way by choosing
// the block size equal to one, or it can be solved using a block 
// implementation (choosing a block size greater than one).
//   
//
#ifndef BELOS_BLOCK_GMRES_HPP
#define BELOS_BLOCK_GMRES_HPP

/*!
  \file BelosBlockGmres.hpp

  \brief Belos concrete class for solving nonsymmetric linear systems with the Generalized Miminum Residual (GMRES) method.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosDGKSOrthoManager.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::ParameterList;
using Teuchos::RefCountPtr;

/*!	
  \class Belos::BlockGmres
  
  \brief This class implements the Restarted Block GMRES algorithm
  for solving real nonsymmetric linear systems of equations AX = B,
  where B is a matrix containing one or more right-hand sides, and 
  X is the matrix of corresponding solutions.
  
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType, class MV, class OP>
  class BlockGmres : public IterativeSolver<ScalarType,MV,OP> { 
  public:
    
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    
    /** @name Constructor/Destructor */
    //@{
    //! %Belos::BlockGmres constructor.
    BlockGmres(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp, 
	       const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest,
               const RefCountPtr<OutputManager<ScalarType> > &om,
	       const RefCountPtr<ParameterList> &pl
	       );
    
    //! %Belos::BlockGmres destructor.
    virtual ~BlockGmres() {};
    //@}
    
    /** @name Accessor methods */
    //@{
    
    //! Get the iteration count for the current block of linear systems.
    int GetNumIters() const { return( _totaliter ); }
    
    //! Get the restart count of the iteration method for the current block of linear systems.
    int GetNumRestarts() const { return( _restartiter ); }
    
    //! Get the solvers native residuals for the current block of linear systems.
    /*! For GMRES this is not the same as the actual residual of the linear system and the
      residual is not in MultiVec form, so the normvec will be populated with the residual norm.
    */
    RefCountPtr<const MV> GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const;
    
    //! Get the true residuals for the current block of linear systems.
    /*! For GMRES this will force the solver to compute a current residual for its linear 
      systems, the current solution is not stored. <b> This is an expensive computation, 
      so a convergence test using these residuals should be secondary to using the native 
      residuals. </b>
    */
    RefCountPtr<MV> GetCurrentSoln();
    
    //! Get a pointer to the current linear problem.  
    /*! This may include a current solution, if the solver has recently restarted or completed.
     */
    RefCountPtr<LinearProblem<ScalarType,MV,OP> > GetLinearProblem() const { return( _lp ); }
    
    //! Get a pointer to the current status test.
    RefCountPtr<StatusTest<ScalarType,MV,OP> > GetStatusTest() const { return( _stest ); }

    //! Reset the solver, can pass in a new parameter list to change solver parameters.
    int Reset( const RefCountPtr<ParameterList>& pl = null );
    
    //@} 

    /** \name Solver application method. */
    //@{
    
    /*! \brief This method uses the iterative method to compute approximate
      solutions to the original problem.  This method can return unconverged if the
      maximum number of iterations is reached, or numerical breakdown is observed.
    */
    void Solve();

    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the block GMRES solver */
    std::string description() const;

    //@}
    
  private:

    //! Method for performing the block Krylov decomposition.
    bool BlockReduction(bool&);

    //! Method for updating QR factorization of upper Hessenberg matrix 
    void UpdateLSQR(Teuchos::SerialDenseMatrix<int,ScalarType>&,Teuchos::SerialDenseMatrix<int,ScalarType>&);

    //! Reference to the linear problem being solver for with the solver. [passed in by user]
    RefCountPtr<LinearProblem<ScalarType,MV,OP> > _lp; 

    //! Reference to the status test, which provides the stopping criteria for the solver. [passed in by user]
    RefCountPtr<StatusTest<ScalarType,MV,OP> > _stest; 

    //! Reference to the output manager for this linear solver. [passed in by user]
    RefCountPtr<OutputManager<ScalarType> > _om;

    //! Reference to the orthogonalization manager for this linear solver.
    RefCountPtr<OrthoManager<ScalarType, MV> > _ortho;

    //! Parameter list containing information for configuring the linear solver. [passed in by user]
    RefCountPtr<ParameterList> _pl;     

    //! Pointers to the Krylov basis constructed by the solver.
    RefCountPtr<MV> _basisvecs;

    //! Pointers to the preconditioned Krylov basis constructed by the solver. [ used in Flexible GMRES ]
    RefCountPtr<MV> _z_basisvecs; 

    //! Pointers to the current right-hand side and solution multivecs being solved for.
    RefCountPtr<MV> _cur_block_rhs, _cur_block_sol;

    //! Dense matrices for holding the upper Hessenberg matrix (H) of the Arnoldi factorization 
    Teuchos::SerialDenseMatrix<int,ScalarType> _hessmatrix;

    //! Dense vector for holding the right-hand side of the least squares problem.
    Teuchos::SerialDenseMatrix<int,ScalarType> _z;

    //! The output stream for sending solver information.
    RefCountPtr<ostream> _os;
    
    //! Length of the Krylov factorization.
    int _length;

    //! Current blocksize, number of restarts, total iteration number, current iteration number.
    int _blocksize, _restartiter, _totaliter, _iter;

    //! Whether this solver is using the flexible variant or not.
    bool _flexible;

    //! Storage for QR factorization of the least-squares system.
    Teuchos::SerialDenseVector<int,ScalarType> beta, sn;
    Teuchos::SerialDenseVector<int,MagnitudeType> cs;

    //! Restart the timers each time Solve() is called.
    bool _restartTimers;
    
    //! Internal timers
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerPrec, 
                                        _timerOrtho, _timerTotal;
  };

  //
  // Implementation
  //

  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class ScalarType, class MV, class OP>
  BlockGmres<ScalarType,MV,OP>::BlockGmres(const RefCountPtr< LinearProblem<ScalarType,MV,OP> >& lp, 
					   const RefCountPtr< StatusTest<ScalarType,MV,OP> >& stest,
					   const RefCountPtr< OutputManager<ScalarType> >& om,
					   const RefCountPtr< ParameterList > &pl ):
    _lp(lp),
    _stest(stest),
    _om(om),
    _pl(pl),
    _os(om->GetOStream()),
    _length(_pl->get("Length",25)),
    _blocksize(0), 
    _restartiter(0), 
    _totaliter(0),
    _iter(0),
    _flexible( (_pl->isParameter("Variant"))&&(Teuchos::getParameter<std::string>(*_pl, "Variant")=="Flexible") ),
    _restartTimers(true),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerPrec(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time"))    
  {
  }
    
  template <class ScalarType, class MV, class OP>
  RefCountPtr<const MV> 
  BlockGmres<ScalarType,MV,OP>::GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
  {
    //
    // If this is the first iteration for a new right-hand side return the
    // residual for the current block rhs and solution.
    // NOTE: Make sure the incoming vector is the correct size!
    //
    if ( normvec && (int)normvec->size() < _blocksize )                         
      normvec->resize( _blocksize );                                          

    if (_totaliter == 0) {
      RefCountPtr<MV> temp_res = MVT::Clone( *_cur_block_rhs, _blocksize );
      _lp->ComputeResVec( &*temp_res, &*_cur_block_sol, &*_cur_block_rhs);
      MVT::MvNorm( *temp_res, normvec );
      return temp_res;
    } else {
      if (normvec) {
        Teuchos::BLAS<int,ScalarType> blas;
        for (int j=0; j<_blocksize; j++)
          (*normvec)[j] = blas.NRM2( _blocksize, &_z(_iter*_blocksize, j ), 1);
      }
    }
    return null;
  }

  template <class ScalarType, class MV, class OP>
  RefCountPtr<MV> 
  BlockGmres<ScalarType,MV,OP>::GetCurrentSoln()
  {    
    //
    // If this is the first iteration of the Arnoldi factorization, 
    // return the current solution.  It has either been updated recently, 
    // if there was a restart, or we haven't computed anything yet.
    //
    RefCountPtr<MV> cur_sol_copy = MVT::CloneCopy(*_cur_block_sol);
    if (_iter==0) { 
        return cur_sol_copy; 
    } else {
      const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
      int i, m = _iter*_blocksize;
      Teuchos::BLAS<int,ScalarType> blas;
      //
      //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> y( Teuchos::Copy, _z, m, _blocksize );
      //
      //  Solve the least squares problem and compute current solutions.
      //
      blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	       Teuchos::NON_UNIT_DIAG, m, _blocksize, one,  
	       _hessmatrix.values(), _hessmatrix.stride(), y.values(), y.stride() );
    
      std::vector<int> index(m);
      for ( i=0; i<m; i++ ) {   
        index[i] = i;
      }
      if (_flexible) {
	RefCountPtr<const MV> Zjp1 = MVT::CloneView( *_z_basisvecs, index );
	MVT::MvTimesMatAddMv( one, *Zjp1, y, one, *cur_sol_copy );
      } 
      else {
	RefCountPtr<const MV> Vjp1 = MVT::CloneView( *_basisvecs, index );
	MVT::MvTimesMatAddMv( one, *Vjp1, y, one, *cur_sol_copy );
      }
    }
    return cur_sol_copy;
  }
    
  template <class ScalarType, class MV, class OP>
  void 
  BlockGmres<ScalarType,MV,OP>::Solve() 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);

    if ( _restartTimers ) {
      _timerOp->reset();
      _timerPrec->reset();
      _timerOrtho->reset();
    }

    int i=0;
    std::vector<int> index;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::RefCountPtr<MV> U_vec;
    bool dep_flg = false, ortho_flg = false;
    bool restart_flg = true;
    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Create the orthogonalization manager.
    //
    _ortho = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>() );
    //
    // Obtain the output stream from the OutputManager.
    //
    _os = _om->GetOStream();
    //
    // Obtain the first block linear system form the linear problem manager.
    //
    _cur_block_sol = _lp->GetCurrLHSVec();
    _cur_block_rhs = _lp->GetCurrRHSVec();
    //
    //  Start executable statements. 
    //
    while ( _cur_block_rhs.get() && _cur_block_sol.get() ) {
      //
      // Get the blocksize from the linear problem
      //
      _blocksize = _lp->GetCurrBlockSize();
      //
      if (_om->isVerbosityAndPrint( IterationDetails )) {
        *_os << endl;
        *_os << "===================================================" << endl;
        *_os << "Solving linear system(s):  "
             << _lp->GetRHSIndex() << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve()-1
             << "  [block size = " << _blocksize << "]\n";
        *_os << endl;
      }
      //
      // Reset the iteration counter for this block of right-hand sides.
      //
      _totaliter = 0;
      //
      // Make room for the Arnoldi vectors and F.
      //
      _basisvecs = MVT::Clone(*_cur_block_rhs,(_length+1)*_blocksize);
      if (_flexible)
        _z_basisvecs = MVT::Clone(*_cur_block_rhs, (_length+1)*_blocksize);
      //
      // Create the rectangular Hessenberg matrix and right-hand side of least squares problem.
      //
      _hessmatrix.shape((_length+1)*_blocksize, _length*_blocksize);
      _z.shape((_length+1)*_blocksize, _blocksize); 
      //
      //
      for ( _restartiter=0; _stest->CheckStatus(this) == Unconverged && restart_flg; ++_restartiter ) {
	//
	// Associate the initial block of _basisvecs with U_vec
	// Reset the index vector (this might have been changed if there was a restart)
	//
	index.resize(_blocksize);
	for (i=0; i < _blocksize; i++) { index[i] = i; }
	U_vec = MVT::CloneView( *_basisvecs, index );
	//
	// Compute current residual and place into 1st block
	//
	_lp->ComputeResVec( &*U_vec, &*_cur_block_sol, &*_cur_block_rhs );
	//
	// Reset orthogonalization failure flags
	//
	dep_flg = false; ortho_flg = false;
	//
	// Re-initialize RHS of the least squares system and create a view.
	//
	_z.putScalar();
	Teuchos::RefCountPtr< Teuchos::SerialDenseMatrix<int,ScalarType> > G10
	  = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, _z, _blocksize, _blocksize) );
	int rank = _ortho->normalize( *U_vec, G10 );
	//
	if (rank != _blocksize){
	  if (_om->isVerbosityAndPrint( Errors )){
	    *_os << "Exiting Block GMRES" << endl;
	    *_os << "  Restart iteration# " << _restartiter
		 << "  Iteration# " << _iter << endl;
	    *_os << "  ERROR: Failed to compute initial block of orthonormal basis vectors"
		 << endl << endl;
	  }
	}
	//
	// NOTE:  U_vec is a view into the set of basis vectors (_basisvecs), which is not guaranteed to
	// be updated when U_vec is changed.  Thus, to ensure consistency between U_vec and _basisvecs,
	// U_vec must be destroyed at this point!
	//
	U_vec = Teuchos::null;
	//
	//	
	for (_iter=0; _iter<_length && _stest->CheckStatus(this) == Unconverged && !ortho_flg; ++_iter, ++_totaliter) {
	  //
	  // Compute a length _length block Arnoldi Reduction (one step at a time),
	  // the exit_flg indicates if we cannot extend the Arnoldi Reduction.
          // If exit_flg is true, then we need to leave this loop and compute the latest solution.
          //
	  // NOTE:  There's an exception here if the blocksize is equal to one, then a lucky
	  // breakdown has occurred, so we should update the least-squares problem.
	  //
	  //dep_flg = true;
	  ortho_flg = BlockReduction(dep_flg);
	  if (ortho_flg && _blocksize > 1){ 
	    break;
	  }
	  //
	  // Update the new Least-Squares solution through the QR factorization of the new
	  // block in the upper Hessenberg matrix.
	  //
	  UpdateLSQR( _hessmatrix, _z );
	  //
	} // end for (_iter=0;...
	//
	// Update the solutions by solving the triangular system to get the Krylov weights.
	//
        if (_iter) {
	  //
	  // Make a copy of _z since it may be used in the convergence test to compute native residuals.
	  Teuchos::SerialDenseMatrix<int,ScalarType> _z_copy( Teuchos::Copy,_z, _iter*_blocksize, _blocksize );	
	  //
	  blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		     Teuchos::NON_UNIT_DIAG, _iter*_blocksize, _blocksize, one,
		     _hessmatrix.values(), _hessmatrix.stride(), _z_copy.values(), _z_copy.stride() ); 
	  //                                                                    
          // Update Solution.                                                   
	  //                                                                      
          // 1)  For the flexible variant the solution will be updated directly,
	  // Check status if _iter==_length                                       
          // otherwise the updated residual will be passed back to the linear problem.
	  //                                                                      
          // 2)  Inform the linear problem that the solution was updated, pass updated residual if necessary.
	  //
	  index.resize( _iter*_blocksize );
	  for (i=0; i < _iter*_blocksize; i++) { index[i] = i; }
	  if (_flexible) {
	    RefCountPtr<const MV> Zjp1 = MVT::CloneView( *_z_basisvecs, index );
	    MVT::MvTimesMatAddMv( one, *Zjp1, _z_copy, one, *_cur_block_sol );
	    _lp->SolutionUpdated();
	  } else {
	    RefCountPtr<const MV> Vjp1 = MVT::CloneView( *_basisvecs, index );
	    RefCountPtr<MV> solnUpdate = MVT::Clone( *_cur_block_sol, _blocksize ); 
	    MVT::MvTimesMatAddMv( one, *Vjp1, _z_copy, zero, *solnUpdate );
	    //
	    // Update the solution held by the linear problem.
	    //	  
	    _lp->SolutionUpdated( solnUpdate.get() );
	    //
	  }
	} // if (_iter) 
	//
	// Check status if _iter==_length
	//
	if (_iter == _length) {
	  _stest->CheckStatus(this);
	}
	//
	// Determine if we are going to allow a restart
	// NOTE:  We will try to restart if we actually computed a subspace of nontrivial dimension (iter!=0)
	//
	restart_flg = (_iter!=0 && restart_flg);
	if (ortho_flg && restart_flg) {
	  if (_om->isVerbosityAndPrint( Warnings )) {
	    *_os << "WARNING: Orthogonalization failure detected at local iteration "
		 << _iter<<", total iteration "<<_totaliter<<", restart will be performed!\n";
	  }
	} 
	//
	// Break out of this loop before the _restartiter is incremented if we are finished.
	//
        if ( _stest->GetStatus() != Unconverged || !restart_flg ) { break; }
        //
      } // end for (_restartiter=0;...
      //
      // Print out solver status
      //
      if (_om->isVerbosityAndPrint( FinalSummary )) {
        *_os << endl;
        _stest->Print(*_os);
        if (ortho_flg && _stest->GetStatus()!=Converged) {
          *_os << " Exiting Block GMRES --- " << endl;
          *_os << "  ERROR: Failed to compute new block of orthonormal basis vectors" << endl;
          *_os << "  ***Solution from previous step will be returned***"<< endl<< endl;
        }
      } 
      //
      // Inform the linear problem that we are finished with this block linear system.
      //	  
      _lp->SetCurrLSVec();
      //
      // Obtain the next block linear system from the linear problem manager.
      //
      _cur_block_sol = _lp->GetCurrLHSVec();
      _cur_block_rhs = _lp->GetCurrRHSVec();
      //
    } // end while( _cur_block_sol && _cur_block_rhs )
    //
    // Print timing details 

    // Stop timer.
    _timerTotal->stop();

    // Reset format that will be used to print the summary
    Teuchos::TimeMonitor::format().setPageWidth(54);

    if (_om->isVerbosity( Belos::TimingDetails )) {
      if (_om->doPrint())
        *_os <<"********************TIMING DETAILS********************"<<endl;
      Teuchos::TimeMonitor::summarize( *_os );
      if (_om->doPrint())
        *_os <<"******************************************************"<<endl;
    }
  } // end Solve()
  
  
  template<class ScalarType, class MV, class OP>
  int 
  BlockGmres<ScalarType,MV,OP>::Reset( const RefCountPtr<ParameterList>& pl )
  {
    // Set new parameter list if one is passed in.
    if (pl.get() != 0 )  
      _pl = pl;
    _length = _pl->get("Length",25);
    _blocksize = 0;
    _restartiter = 0; 
    _totaliter = 0;
    _iter = 0;
    //
    // If there is a "Variant" parameter in the list, check to see if it's "Flexible" (i.e. flexible GMRES)
    //
    if (_pl->isParameter("Variant"))
      _flexible = (Teuchos::getParameter<std::string>(*_pl, "Variant")=="Flexible");
    return 0;
  }

  template<class ScalarType, class MV, class OP>
  bool 
  BlockGmres<ScalarType,MV,OP>::BlockReduction ( bool& dep_flg ) 
  {
    //
    int i;	
    bool flg = false;
    std::vector<int> index( _blocksize );
    //
    // Associate the j-th block of _basisvecs with U_vec.
    //
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = _iter*_blocksize+i;
    }
    RefCountPtr<MV> U_vec = MVT::CloneView( *_basisvecs, index );
    //
    RefCountPtr<MV> Z_vec;
    if (_flexible) {                                                            
      Z_vec = MVT::CloneView( *_z_basisvecs, index );             
    }
    //
    // Associate the j+1-th block of _basisvecs with AU_vec.
    //
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = (_iter+1)*_blocksize+i;
    }
    RefCountPtr<MV> AU_vec = MVT::CloneView( *_basisvecs, index );
    //                                                                          
    // If this is the flexible variant apply operator separately, else apply composite operator.
    //                                                                          
    if (_flexible) {                                                            
      //
      //  Apply right preconditioning and store it in _z_basisvecs.             
      //                                                                        
      {
	Teuchos::TimeMonitor PrecTimer(*_timerPrec);
	_lp->ApplyRightPrec( *U_vec, *Z_vec );                                    
      }
      //                                                                        
      //  Apply operator and store it in AU_vec.                                
      //                                                                        
      {
	Teuchos::TimeMonitor OpTimer(*_timerOp);
	_lp->ApplyOp( *Z_vec, *AU_vec );                                          
      }
    } 
    else                                                                      
      { 
	Teuchos::TimeMonitor OpTimer(*_timerOp);
	_lp->Apply( *U_vec, *AU_vec ); 
      }

    // Clean up pointer to the j-th block of _basisvecs
    U_vec = Teuchos::null;

    //
    // Resize index and get previous vectors.
    //
    int num_prev = (_iter+1)*_blocksize;
    index.resize( num_prev );
    for (i=0; i<num_prev; i++) { 
      index[i] = i; 
    }
    RefCountPtr<MV> V_prev = MVT::CloneView( *_basisvecs, index );
    Teuchos::Array< RefCountPtr<const MV> > V_array( 1, V_prev );
    //
    // Get a view of the current part of the upper-hessenberg matrix.
    //
    RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > h_new 
      = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, _hessmatrix, num_prev, _blocksize, 0, _iter*_blocksize ) );
    Teuchos::Array< RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > h_array( 1, h_new );

    RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > r_new
      = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, _hessmatrix, _blocksize, _blocksize, num_prev, _iter*_blocksize ) );
    //
    // Orthonormalize the new block of the Krylov expansion
    // 
    int rank = _ortho->projectAndNormalize( *AU_vec, h_array, r_new, V_array );

    if (rank != _blocksize) {
      flg = true;
    } 

    return flg;

    //
  } // end BlockReduction()
    
  template<class ScalarType, class MV, class OP>
  void 
  BlockGmres<ScalarType,MV,OP>::UpdateLSQR( Teuchos::SerialDenseMatrix<int,ScalarType>& R, 
					    Teuchos::SerialDenseMatrix<int,ScalarType>& z )
  {
    int i, j, maxidx;
    ScalarType sigma, mu, vscale, maxelem;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Apply previous transformations and compute new transformation to reduce upper-Hessenberg
    // system to upper-triangular form.  
    // NOTE:  When _iter==0, storage will be created for the transformations.
    //
    if (_blocksize == 1) 
      {
	if (_iter==0) {
	  cs.resize( _length+1 );
	  sn.resize( _length+1 );
	}
	//
	// QR factorization of Least-Squares system with Givens rotations
	//
	for (i=0; i<_iter; i++) {
	  //
	  // Apply previous Givens rotations to new column of Hessenberg matrix
	  //
	  blas.ROT( 1, &R(i,_iter), 1, &R(i+1, _iter), 1, &cs[i], &sn[i] );
	}
	//
	// Calculate new Givens rotation
	//
	blas.ROTG( &R(_iter,_iter), &R(_iter+1,_iter), &cs[_iter], &sn[_iter] );
	R(_iter+1,_iter) = zero;
	//
	// Update RHS w/ new transformation
	//
	blas.ROT( 1, &z(_iter,0), 1, &z(_iter+1,0), 1, &cs[_iter], &sn[_iter] );
      } 
    else
      {
	if (_iter==0) {
	  beta.size((_length+1)*_blocksize);
	}
	//
	// QR factorization of Least-Squares system with Householder reflectors
	//
	for (j=0; j<_blocksize; j++) {
	  //
	  // Apply previous Householder reflectors to new block of Hessenberg matrix
	  //
	  for (i=0; i<_iter*_blocksize+j; i++) {
	    sigma = blas.DOT( _blocksize, &R(i+1,i), 1, &R(i+1,_iter*_blocksize+j), 1);
	    sigma += R(i,_iter*_blocksize+j);
	    sigma *= beta[i];
	    blas.AXPY(_blocksize, -sigma, &R(i+1,i), 1, &R(i+1,_iter*_blocksize+j), 1);
	    R(i,_iter*_blocksize+j) -= sigma;
	  }
	  //
	  // Compute new Householder reflector
	  //
	  maxidx = blas.IAMAX( _blocksize+1, &R(_iter*_blocksize+j,_iter*_blocksize+j), 1 );
	  maxelem = R(_iter*_blocksize+j+maxidx-1,_iter*_blocksize+j);
	  for (i=0; i<_blocksize+1; i++) 
	    R(_iter*_blocksize+j+i,_iter*_blocksize+j) /= maxelem;
	  sigma = blas.DOT( _blocksize, &R(_iter*_blocksize+j+1,_iter*_blocksize+j), 1, 
			    &R(_iter*_blocksize+j+1,_iter*_blocksize+j), 1 );
	  if (sigma == zero) {
	    beta[_iter*_blocksize + j] = zero;
	  } else {
	    mu = sqrt(R(_iter*_blocksize+j,_iter*_blocksize+j)*R(_iter*_blocksize+j,_iter*_blocksize+j)+sigma);
	    if ( Teuchos::ScalarTraits<ScalarType>::real(R(_iter*_blocksize+j,_iter*_blocksize+j)) < Teuchos::ScalarTraits<MagnitudeType>::zero() ) {
	      vscale = R(_iter*_blocksize+j,_iter*_blocksize+j) - mu;
	    } else {
	      vscale = -sigma / (R(_iter*_blocksize+j,_iter*_blocksize+j) + mu);
	    }
	    beta[_iter*_blocksize+j] = 2.0*vscale*vscale/(sigma + vscale*vscale);
	    R(_iter*_blocksize+j,_iter*_blocksize+j) = maxelem*mu;
	    for (i=0; i<_blocksize; i++)
	      R(_iter*_blocksize+j+1+i,_iter*_blocksize+j) /= vscale;
	  }
	  //
	  // Apply new Householder reflector to rhs
	  //
	  for (i=0; i<_blocksize; i++) {
	    sigma = blas.DOT( _blocksize, &R(_iter*_blocksize+j+1,_iter*_blocksize+j), 
			      1, &z(_iter*_blocksize+j+1,i), 1);
	    sigma += z(_iter*_blocksize+j,i);
	    sigma *= beta[_iter*_blocksize+j];
	    blas.AXPY(_blocksize, -sigma, &R(_iter*_blocksize+j+1,_iter*_blocksize+j), 
		      1, &z(_iter*_blocksize+j+1,i), 1);
	    z(_iter*_blocksize+j,i) -= sigma;
	  }
	}
      } // end if (_blocksize == 1)
  } // end UpdateLSQR
  
  // Overridden from Teuchos::Describable

  template<class ScalarType, class MV, class OP>
  std::string BlockGmres<ScalarType,MV,OP>::description() const
  {
    std::ostringstream oss;
    oss << "Belos::BlockGmres<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
    oss << "{";
    oss << "Variant=\'"<<(_flexible?"Flexible":"Standard")<<"\'";
    oss << "}";
    return oss.str();
  }

} // end namespace Belos

#endif
// End of file BelosBlockGmres.hpp


