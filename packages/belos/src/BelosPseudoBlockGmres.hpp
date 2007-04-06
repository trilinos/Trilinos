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
#ifndef BELOS_PSEUDO_BLOCK_GMRES_HPP
#define BELOS_PSEUDO_BLOCK_GMRES_HPP

/*!
  \file BelosPseudoBlockGmres.hpp

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
#include "BelosICGSOrthoManager.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Array.hpp"

using Teuchos::ParameterList;
using Teuchos::RefCountPtr;

/*!	
  \class Belos::PseudoBlockGmres
  
  \brief This class implements the Restarted Block GMRES algorithm
  for solving real nonsymmetric linear systems of equations AX = B,
  where B is a matrix containing one or more right-hand sides, and 
  X is the matrix of corresponding solutions.  Each linear system
  is solved independently, where the matrix-vector products are 
  aggregated for all the right-hand sides. 
  
  \author Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType, class MV, class OP>
  class PseudoBlockGmres : public IterativeSolver<ScalarType,MV,OP> { 
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
    //! %Belos::PseudoBlockGmres constructor.
    PseudoBlockGmres(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp, 
		     const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest,
		     const RefCountPtr<OutputManager<ScalarType> > &om,
		     const RefCountPtr<ParameterList> &pl
		     );
    
    //! %Belos::PseudoBlockGmres destructor.
    virtual ~PseudoBlockGmres() {};
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
    void UpdateLSQR(Teuchos::Array<Teuchos::SerialDenseMatrix<int,ScalarType> >&,
		    Teuchos::Array<Teuchos::SerialDenseVector<int,ScalarType> >&);

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
    Teuchos::Array<RefCountPtr<MV> > _basisvecs;

    //! Pointers to a work vector used to improve aggregate performance.
    RefCountPtr<MV> _U_vec, _AU_vec;

    //! Pointers to the current right-hand side and solution multivecs being solved for.
    RefCountPtr<MV> _cur_block_rhs, _cur_block_sol;

    //! Dense matrices for holding the upper Hessenberg matrix (H) of the Arnoldi factorization 
    Teuchos::Array<Teuchos::SerialDenseMatrix<int,ScalarType> > _hessmatrix;

    //! Dense vector for holding the right-hand side of the least squares problem.
    Teuchos::Array<Teuchos::SerialDenseVector<int,ScalarType> > _z;

    //! The output stream for sending solver information.
    RefCountPtr<ostream> _os;
    
    //! Length of the Krylov factorization.
    int _length;

    //! Current blocksize, number of restarts, total iteration number, current iteration number.
    int _blocksize, _restartiter, _totaliter, _iter;

    //! What type of orthogonalization being used in this solver.
    string _orthoType;

    //! Storage for QR factorization of the least-squares system.
    Teuchos::Array<Teuchos::SerialDenseVector<int,ScalarType> > beta, sn;
    Teuchos::Array<Teuchos::SerialDenseVector<int,MagnitudeType> > cs;

    //! Restart the timers each time Solve() is called.
    bool _restartTimers;
    
    //! Internal timers
    Teuchos::RefCountPtr<Teuchos::Time> _timerOrtho, _timerTotal;
  };

  //
  // Implementation
  //

  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class ScalarType, class MV, class OP>
  PseudoBlockGmres<ScalarType,MV,OP>::PseudoBlockGmres(const RefCountPtr< LinearProblem<ScalarType,MV,OP> >& lp, 
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
    _orthoType( "DGKS" ),
    _restartTimers(true),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Belos: Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Belos: Total time"))    
  {
    // Check the parameter list to see if another type of orthogonalization is specified.
    if (_pl->isParameter("Ortho Type")) {
      _orthoType = Teuchos::getParameter<std::string>(*_pl, "Ortho Type" );
    }
    //
    // Create the orthogonalization manager.
    //
    if (_orthoType=="ICGS")
      _ortho = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>() );
    else 
      _ortho = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>() );

    // Check the parameter list to see if the timers should be restarted each solve.
    if (_pl->isParameter("Restart Timers")) {
      _restartTimers = Teuchos::getParameter<bool>(*_pl, "Restart Timers" );
    }
  }
    
  template <class ScalarType, class MV, class OP>
  RefCountPtr<const MV> 
  PseudoBlockGmres<ScalarType,MV,OP>::GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
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
          (*normvec)[j] = blas.NRM2( 1, &_z[j](_iter), 1 );
      }
    }
    return null;
  }
  
  template <class ScalarType, class MV, class OP>
  RefCountPtr<MV> 
  PseudoBlockGmres<ScalarType,MV,OP>::GetCurrentSoln()
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
      std::vector<int> index(1), index2(_iter);
      for (int i=0; i<_iter; ++i) {
        index2[i] = i;
      }
      const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
      Teuchos::BLAS<int,ScalarType> blas;
      
      for (int i=0; i<_blocksize; ++i) {
        index[0] = i;
        RefCountPtr<MV> cur_block_copy_vec = MVT::CloneView( *cur_sol_copy, index );
        //
        //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
        //
        Teuchos::SerialDenseVector<int,ScalarType> y( Teuchos::Copy, _z[i].values(), _iter );
        //
        //  Solve the least squares problem and compute current solutions.
        //
        blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	           Teuchos::NON_UNIT_DIAG, _iter, 1, one,  
		   _hessmatrix[i].values(), _hessmatrix[i].stride(), y.values(), y.stride() );
	
	RefCountPtr<const MV> Vjp1 = MVT::CloneView( *_basisvecs[i], index2 );
	MVT::MvTimesMatAddMv( one, *Vjp1, y, one, *cur_block_copy_vec );
      }
    }
    return cur_sol_copy;
  }
    
  template <class ScalarType, class MV, class OP>
  void 
  PseudoBlockGmres<ScalarType,MV,OP>::Solve() 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);
    
    if ( _restartTimers ) {
      _timerOrtho->reset();
    }
    
    int i=0;
    std::vector<int> index, index2;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::RefCountPtr<MV> tmp_vec;
    bool dep_flg = false, ortho_flg = false;
    bool restart_flg = true;
    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
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
      _basisvecs.resize(_blocksize);
      for (i=0; i<_blocksize; ++i) {
        _basisvecs[i] = MVT::Clone(*_cur_block_rhs,_length+1);
      }
      //
      // Create the rectangular Hessenberg matrix and right-hand side of least squares problem.
      //
      _hessmatrix.resize(_blocksize);
      _z.resize(_blocksize);
      for (i=0; i<_blocksize; ++i) {
	_hessmatrix[i].shapeUninitialized(_length+1, _length);
        _z[i].shapeUninitialized(_length+1, 1); 
      }
      //
      // Outer restart loop begins here.
      //
      for ( _restartiter=0; _stest->CheckStatus(this) == Unconverged && restart_flg; ++_restartiter ) {
        //
	// Associate each initial block of _basisvecs[i] with U_vec[i]
	// Reset the index vector (this might have been changed if there was a restart)
	//
	index.resize(1);
	index2.resize(1);
	index[0] = 0;
	_U_vec = MVT::Clone( *_basisvecs[0], _blocksize );
        for (i=0; i<_blocksize; ++i) {
	  index2[0] = i;
	  tmp_vec = MVT::CloneView( *_basisvecs[i], index );
	  MVT::MvAddMv( one, *tmp_vec, zero, *tmp_vec, *MVT::CloneView( *_U_vec, index2 ) );
        }
	//
	// Compute current residual and place into 1st block
	//
	_lp->ComputeResVec( &*_U_vec, &*_cur_block_sol, &*_cur_block_rhs );
	//
	// Reset orthogonalization failure flags
	//
	dep_flg = false; ortho_flg = false;
	//
	// Re-initialize RHS of the least squares system and create a view.
	//
        for (i=0; i<_blocksize; ++i) {
	  _z[i].putScalar();
        }
        for (i=0; i<_blocksize; ++i) {
	  index2[0] = i;
	  tmp_vec = MVT::CloneView( *_U_vec, index2 );
          _ortho->normalize( *tmp_vec, Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>
						     (Teuchos::View, _z[i].values(), 1) ) );
	  MVT::MvAddMv( one, *tmp_vec, zero, *tmp_vec, *MVT::CloneView( *_basisvecs[i], index ) );
	}
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
	  index.resize( _iter );
	  index2.resize( 1 );
          // 
          //  Create a multivector to store the solution updates
          //
	  RefCountPtr<MV> solnUpdate = MVT::Clone( *_cur_block_sol, _blocksize ); 
          //
	  for (i=0; i<_blocksize; ++i) {
	    //
	    // Make a copy of _z since it may be used in the convergence test to compute native residuals.
	    Teuchos::SerialDenseVector<int,ScalarType> _z_copy( Teuchos::Copy, _z[i].values(), _iter );	
	    //
	    blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		       Teuchos::NON_UNIT_DIAG, _iter, 1, one, _hessmatrix[i].values(), 
		       _hessmatrix[i].stride(), _z_copy.values(), _z_copy.stride() ); 
	    //                                                                    
	    // Update Solution.                                                   
	    //                                                                      
	    // 1)  The updated residual will be passed back to the linear problem.
	    //                                                                      
	    // 2)  Inform the linear problem that the solution was updated, pass updated residual if necessary.
	    //
	    for (int j=0; j < _iter; j++) { index[j] = j; }
	    RefCountPtr<const MV> Vjp1 = MVT::CloneView( *_basisvecs[i], index );
            index2[0] = i;
	    MVT::MvTimesMatAddMv( one, *Vjp1, _z_copy, zero, *MVT::CloneView( *solnUpdate, index2 ) );
	  }
	  //
	  // Update the solution held by the linear problem.
	  //	  
	  _lp->SolutionUpdated( solnUpdate.get() );
	  //
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
  PseudoBlockGmres<ScalarType,MV,OP>::Reset( const RefCountPtr<ParameterList>& pl )
  {
    // Set new parameter list if one is passed in.
    if (pl.get() != 0 )  
      _pl = pl;
    _blocksize = 0;
    _restartiter = 0; 
    _totaliter = 0;
    _iter = 0;

    _length = _pl->get("Length", 25);
    if (_pl->isParameter("Ortho Type")) {
      std::string newOrthoType = Teuchos::getParameter<std::string>(*_pl,"Ortho Type");
      if (newOrthoType != _orthoType) {
        if (newOrthoType=="ICGS")
          _ortho = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>() );
        else 
          _ortho = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>() );
        _orthoType = newOrthoType;
      }
    }
    if (_pl->isParameter("Restart Timers")) {
      _restartTimers = Teuchos::getParameter<bool>(*_pl,"Restart Timers"); 
    }
    return 0;
  }

  template<class ScalarType, class MV, class OP>
  bool 
  PseudoBlockGmres<ScalarType,MV,OP>::BlockReduction ( bool& dep_flg ) 
  {
    //
    int i;	
    bool flg = false;
    std::vector<int> index, index2;
    RefCountPtr<MV> tmp_vec;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    //
    // Associate each of the _iter-th block of _basisvecs[i] with _U_vec[i]
    //
    index.resize(1);
    index2.resize(1);
    index[0] = _iter;
    for (i=0; i<_blocksize; ++i) {
      index2[0] = i;
      tmp_vec = MVT::CloneView( *_basisvecs[i], index );
      MVT::MvAddMv( one, *tmp_vec, zero, *tmp_vec, *MVT::CloneView( *_U_vec, index2 ) );
    }
    //
    // Create _AU_vec to hold A*_U_vec.
    //
    if (_AU_vec == Teuchos::null) {
      _AU_vec = MVT::Clone( *_basisvecs[0], _blocksize );
    }

    // Apply the operator to _work_vector
    _lp->Apply( *_U_vec, *_AU_vec );
    //
    // Resize index.
    //
    int num_prev = _iter+1;
    index.resize( num_prev );
    for (i=0; i<num_prev; ++i) { 
      index[i] = i; 
    }
    for (i=0; i<_blocksize; ++i) {
      //
      // Get previous Krylov vectors.
      //
      RefCountPtr<MV> V_prev = MVT::CloneView( *_basisvecs[i], index );
      Teuchos::Array< RefCountPtr<const MV> > V_array( 1, V_prev );
      //
      // Get a view of the new candidate vector.
      //
      index2[0] = i;
      RefCountPtr<MV> V_new = MVT::CloneView( *_AU_vec, index2 );
      //
      // Get a view of the current part of the upper-hessenberg matrix.
      //
      RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > h_new 
	= Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			( Teuchos::View, _hessmatrix[i], num_prev, 1, 0, _iter ) );
      Teuchos::Array< RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > h_array( 1, h_new );

      RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > r_new
	= Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			( Teuchos::View, _hessmatrix[i], 1, 1, num_prev, _iter ) );
      //
      // Orthonormalize the new block of the Krylov expansion
      // 
      int rank = 0;
      {
	Teuchos::TimeMonitor OrthoTimer(*_timerOrtho);
	rank = _ortho->projectAndNormalize( *V_new, h_array, r_new, V_array );
      }
      //
      // NOTE:  V_new is a copy of the iter+1 vector in _basisvecs[i], so the normalized vector has to be
      // be copied back in when V_new is changed.  
      //
      index2[0] = _iter+1;
      tmp_vec = MVT::CloneView( *_basisvecs[i], index2 );
      MVT::MvAddMv( one, *V_new, zero, *V_new, *tmp_vec );
    }
    // 
    // Now _AU*vec is the new _U*vec 
    // Note: Right now the solver shouldn't do this because it already copies in the _U_vec above.
    //	
    //_U_vec = _AU*vec;
    //
    return flg;
    //
  } // end BlockReduction()
    
  
  template<class ScalarType, class MV, class OP>
  void 
  PseudoBlockGmres<ScalarType,MV,OP>::UpdateLSQR( Teuchos::Array<Teuchos::SerialDenseMatrix<int,ScalarType> >& R, 
						  Teuchos::Array<Teuchos::SerialDenseVector<int,ScalarType> >& z )
  {
    int i, j;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Apply previous transformations and compute new transformation to reduce upper-Hessenberg
    // system to upper-triangular form.  
    // NOTE:  When _iter==0, storage will be created for the transformations.
    //
    if (_iter==0 && (int)cs.size()!=_blocksize) {
      cs.resize(_blocksize);
      sn.resize(_blocksize);
    }

    for (i=0; i<_blocksize; ++i) {
      //
      //  Update the least-squares QR for each linear system.
      //
      if ( _iter==0 ) {
	cs[i].shapeUninitialized(_length+1, 1);
	sn[i].shapeUninitialized(_length+1, 1);
      }
      //
      // QR factorization of Least-Squares system with Givens rotations
      //
      for (j=0; j<_iter; j++) {
	//
	// Apply previous Givens rotations to new column of Hessenberg matrix
	//
	blas.ROT( 1, &R[i](j,_iter), 1, &R[i](j+1, _iter), 1, &(cs[i])[j], &(sn[i])[j] );
      }
      //
      // Calculate new Givens rotation
      //
      blas.ROTG( &R[i](_iter,_iter), &R[i](_iter+1,_iter), &(cs[i])[_iter], &(sn[i])[_iter] );
      R[i](_iter+1,_iter) = zero;
      //
      // Update RHS w/ new transformation
      //
      blas.ROT( 1, &z[i](_iter), 1, &z[i](_iter+1), 1, &cs[i][_iter], &sn[i][_iter] );
    }

  } // end UpdateLSQR
  
  // Overridden from Teuchos::Describable

  template<class ScalarType, class MV, class OP>
  std::string PseudoBlockGmres<ScalarType,MV,OP>::description() const
  {
    std::ostringstream oss;
    oss << "Belos::PseudoBlockGmres<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
    oss << "{";
    oss << "Ortho Type=\'"<<_orthoType<<"\'";
    oss << "}";
    return oss.str();
  }

} // end namespace Belos

#endif
// End of file BelosPseudoBlockGmres.hpp


