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
// This file contains an implementation of the Block CG algorithm
// for solving real symmetric positive definite linear systems of 
// equations AX = B, where B is a matrix containing one or more right-hand
// sides, and X is the matrix of corresponding solutions. This implementation 
// allows the user to solve systems involving any number of right-hand sides.
// The block size used in the solver is user specified, and is independent
// of the number of right-hand sides. Thus, a system involving many 
// right-hand sides can be processed by solving for only some number  
// (the block size) of the right-hand sides simultaneously. Several passes
// through the block solver are used to solve for all of them. A single
// right-hand side system can be solved in the traditional way by choosing
// the block size equal to one, or it can be solved using a block 
// implementation (choosing a block size greater than one).
//
//
//
#ifndef BELOS_BLOCK_CG_HPP
#define BELOS_BLOCK_CG_HPP

/*!
  \file BelosBlockCG.hpp
  
  \brief Belos concrete class for solving symmetric positive definite linear systems with the
  preconditioned Conjugate Gradient (CG) method.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOperator.hpp"
#include "BelosStatusTest.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosCG.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	\class Belos::BlockCG

	\brief This class implements the preconditioned Conjugate Gradient algorithm for
	solving real symmetric positive definite linear systems of equations
	AX = B, where B is a matrix containing one or more right-hand sides.

	\author Teri Barth and Heidi Thornquist
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class BlockCG : public IterativeSolver<ScalarType,MV,OP> { 
public:
  //
  // Convenience typedefs
  //
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

  //@{ \name Constructor/Destructor.

  //! %Belos::BlockCG constructor.
  BlockCG(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp,
	  const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest,
	  const RefCountPtr<OutputManager<ScalarType> > &om,	  
	  const RefCountPtr<ParameterList> &pl = Teuchos::null
	  );

  //! %BlockCG destructor.
  virtual ~BlockCG() {};
  //@}
  
  //@{ \name Accessor methods
  
  //! Get the iteration count for the current block of linear systems.
  int GetNumIters() const { return( _iter ); }
  
  //! Get the restart count of the iteration method for the current block of linear systems [not valid for CG].
  int GetNumRestarts() const { return(0); }
   
  //! Get the solvers native residuals for the current block of linear systems.
  /*! 
      \note The memory for the residual MultiVec must be handled by the calling routine.
   */
  RefCountPtr<const MV> GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const;
  
  //! Get the actual residual vectors for the current block of linear systems.
  /*! This may force the solver to compute a current residual for its linear
  	systems.  For CG, this method is not useful since the linear problem
	manager always has the current solution (even when the blocksize is larger
	than the current number of linear systems being solved for).  
  */
  RefCountPtr<MV> GetCurrentSoln() { return MVT::CloneCopy( *_cur_block_sol ); };
  
  //! Get a constant reference to the current linear problem.  
  /*! This may include a current solution, if the solver has recently restarted or completed.
   */
  RefCountPtr<LinearProblem<ScalarType,MV,OP> > GetLinearProblem() const { return( _lp ); }

  RefCountPtr<StatusTest<ScalarType,MV,OP> > GetStatusTest() const { return( _stest ); }

  //@} 
  
  //@{ \name Solver application method.
  
  /*! \brief This method uses the iterative method to compute approximate solutions
    to the original problem.  This method can return unconverged if the maximum number
    of iterations is reached, or numerical breakdown is observed.
  */
  void Solve();
  //@}
    
private:

  void SetCGBlkTols();
  bool QRFactorDef(MV&, Teuchos::SerialDenseMatrix<int,ScalarType>&, std::vector<int>* );
  void CheckOrthogonality(MV&, MV&);
  void CheckAOrthogonality(MV&, MV&);
  void PrintCGIterInfo(const std::vector<int> &ind);
  void BlockIteration();

  //! Linear problem manager [ must be passed in by the user ]
  RefCountPtr<LinearProblem<ScalarType,MV,OP> > _lp; 

  //! Status test [ must be passed in by the user ]
  RefCountPtr<StatusTest<ScalarType,MV,OP> > _stest; 

  //! Output manager [ must be passed in by the user ]
  RefCountPtr<OutputManager<ScalarType> > _om;

  //! Parameter list containing information for configuring the linear solver. [ must be passed in by the user ]
  RefCountPtr<ParameterList> _pl;     
  
  //! Pointer to current linear systems block of solution vectors [obtained from linear problem manager]
  RefCountPtr<MV> _cur_block_sol;

  //! Pointer to current linear systems block of right-hand sides [obtained from linear problem manager]
  RefCountPtr<MV> _cur_block_rhs; 

  //! Pointer to block of the current residual vectors.
  RefCountPtr<MV> _residvecs;

  //! Output stream.
  RefCountPtr<ostream> _os;

  //! Current blocksize, iteration number, and basis pointer.
  int _blocksize, _iter, _new_blk;

  //! Numerical breakdown tolerances.
  MagnitudeType _prec, _dep_tol;
  
  //! Restart the timers each time Solve() is called.
  bool _restartTimers;

  //! Internal timers
  Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerPrec, 
                                      _timerOrtho, _timerTotal;
};
  
  //
  // Implementation
  //
  
  template <class ScalarType, class MV, class OP>
  BlockCG<ScalarType,MV,OP>::BlockCG(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp,
				     const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest,
				     const RefCountPtr<OutputManager<ScalarType> > &om,
				     const RefCountPtr<ParameterList> &pl
				     ) : 
    _lp(lp), 
    _stest(stest),
    _om(om),
    _pl(pl),
    _os(om->GetOStream()),
    _blocksize(0), 
    _iter(0),
    _new_blk(1),
    _prec(1.0), 
    _dep_tol(1.0),
    _restartTimers(true),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerPrec(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time"))
  { 
    //
    // Set the block orthogonality tolerances
    //
    SetCGBlkTols();
  }
  
  template <class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::SetCGBlkTols() 
  {
    typedef typename Teuchos::ScalarTraits<MagnitudeType> MGT;
    const MagnitudeType two = 2.0;
    const MagnitudeType eps = SCT::eps();
    _prec = eps;
    _dep_tol = MGT::one()/MGT::squareroot(two);
  }
  
  template <class ScalarType, class MV, class OP>
  RefCountPtr<const MV> 
  BlockCG<ScalarType,MV,OP>::GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
  {
    int i;
    std::vector<int> index( _blocksize );
    if (_new_blk)
      for (i=0; i<_blocksize; i++) { index[i] = i; }
    else
      for (i=0; i<_blocksize; i++) { index[i] = _blocksize + i; }
    RefCountPtr<MV> ResidMV = MVT::CloneView( *_residvecs, index );
    return ResidMV;
  }
  

  template <class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::Solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);

    if ( _restartTimers ) {
      _timerOp->reset();
      _timerPrec->reset();
      _timerOrtho->reset();
    }
    //
    // Get the current ostream from the OutputManager
    //
    _os = _om->GetOStream();
    //
    // Retrieve the first linear system to be solved.
    //
    _cur_block_sol = _lp->GetCurrLHSVec();
    _cur_block_rhs = _lp->GetCurrRHSVec();
    //
    //  Start executable statements. 
    //
    while (_cur_block_sol.get() && _cur_block_rhs.get() ) {
      //
      // Get the blocksize for this set of linear systems.
      //
      _blocksize = _lp->GetCurrBlockSize();
      if (_blocksize == 1 ) {
	//
	// Create single vector CG solver for this linear system.
	//
	Belos::CG<ScalarType,MV,OP> CGSolver(_lp, _stest, _om);
	CGSolver.Solve();
	//
      } else {
	BlockIteration();
      }
      //
      // Get the next block of linear systems, if it returns the null pointer we are done.
      //
      _cur_block_sol = _lp->GetCurrLHSVec();
      _cur_block_rhs = _lp->GetCurrRHSVec();
      //
    } // end while ( _cur_block_sol && _cur_block_rhs )
    // **********************************************************************************
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
    //
  } // end Solve()
  //  
  
  template <class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::BlockIteration ( ) 
  {
    //
    int i, j, k, info, num_ind;
    int ind_blksz, prev_ind_blksz;
    bool exit_flg = false;
    bool isPrec = (_lp->GetLeftPrec().get()!=NULL);
    char UPLO = 'U';
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::LAPACK<int,ScalarType> lapack;
    RefCountPtr<MV> P_prev, R_prev, AP_prev, P_new, R_new;
    //
    // Make additional space needed during iteration
    //
    std::vector<int> index2( _blocksize );
    std::vector<int> index( _blocksize );
    std::vector<int> ind_idx( _blocksize );
    std::vector<int> cols( _blocksize );
    //
    // Make room for the direction, residual, and operator applied to direction vectors.
    // We save 3 blocks of these vectors.  Also create 2 vectors of _blocksize to store
    // The residual norms used to determine valid direction/residual vectors.
    //
    RefCountPtr<MV> _basisvecs = MVT::Clone( *_cur_block_sol, 2*_blocksize ); 
    RefCountPtr<MV> temp_blk = MVT::Clone( *_cur_block_sol, _blocksize );
    std::vector<MagnitudeType> _cur_resid_norms( _blocksize );
    std::vector<MagnitudeType> _init_resid_norms( _blocksize );
    _residvecs = MVT::Clone( *_cur_block_sol, 2*_blocksize ); 
    //
    ind_blksz = _blocksize;
    prev_ind_blksz = _blocksize;
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha( _blocksize, _blocksize );
    Teuchos::SerialDenseMatrix<int,ScalarType> beta( _blocksize, _blocksize );
    Teuchos::SerialDenseMatrix<int,ScalarType> T2( _blocksize, _blocksize );
    //
    for (i=0;i<_blocksize;i++){
      index[i] = i;
      ind_idx[i] = i; 
      cols[i] = i;
    }
    //
    if (_om->isVerbosityAndPrint( IterationDetails )) {
      *_os << endl;
      *_os << "===================================================" << endl;
      *_os << "Solving linear system(s):  " << _lp->GetRHSIndex() 
	   << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve() << endl;
      *_os << endl;
    }	
    //
    //
    // ************ Compute the initial residuals ********************************
    //
    // Associate the first block of _basisvecs with P_prev and the
    // first block of _residvecs with R_prev
    //
    P_prev = MVT::CloneView( *_basisvecs, ind_idx );
    R_prev = MVT::CloneView( *_residvecs, ind_idx );
    AP_prev = MVT::CloneView( *temp_blk, ind_idx );
    //
    // Store initial guesses to AX = B in 1st block of _basisvecs
    //         P_prev = one*cur_block_sol + zero*P_prev
    //
    MVT::MvAddMv(one, *_cur_block_sol, zero, *_cur_block_sol, *P_prev);
    //
    // Multiply by A and store in AP_prev
    //       AP_prev = A*P_prev
    //
    {
      Teuchos::TimeMonitor OpTimer(*_timerOp);
      _lp->ApplyOp( *P_prev, *AP_prev );
    }
    //
    // Compute initial residual block and store in 1st block of _residvecs
    //     R_prev = cur_block_rhs - A*P_prev
    //
    MVT::MvAddMv(one, *_cur_block_rhs, -one, *AP_prev, *R_prev );
    //
    //-------Compute and save the initial residual norms----------
    //
    MVT::MvNorm(*R_prev, &_init_resid_norms);
    //
    // Update indices of current (independent) blocks.
    // If a residual is too small, it will be dropped from
    // the current block, thus, from future computations
    //
    j = 0;
    ind_idx.resize( 0 );
    for (i=0; i<_blocksize; i++){
      _cur_resid_norms[i] = _init_resid_norms[i];
      if (_init_resid_norms[i] > _prec)
	ind_idx.push_back( i );
    }
    ind_blksz = ind_idx.size(); 
    //
    if (ind_blksz > 0) { 
      //
      // All initial residuals have not converged -- continue Block CG	
      // Compute the initial block of direciton vectors
      //
      // Associate current blocks of residuals, directions, and solution block
      // with R_prev, P_prev, and cur_sol
      //
      R_prev = MVT::CloneView( *_residvecs, ind_idx );
      P_prev = MVT::CloneView( *_basisvecs, ind_idx );
      //
      //----------------Compute initial direction vectors--------------------------
      // Initially, they are set to the preconditioned residuals
      //
      if ( isPrec ) {
	Teuchos::TimeMonitor PrecTimer(*_timerPrec);
	_lp->ApplyLeftPrec( *R_prev, *P_prev );
      } else {
	MVT::MvAddMv( one , *R_prev, zero, *R_prev, *P_prev); 
      }
      //
      // Compute an orthonormal block of initial direction vectors,
      // and check for dependencies, adjusting indices of independent
      // vectors if needed
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> G(ind_blksz, ind_blksz);
      exit_flg = false;
      {
	Teuchos::TimeMonitor OrthoTimer(*_timerOrtho);
	exit_flg = QRFactorDef(*P_prev, G, &cols );
      }
      num_ind = cols.size();
      //
      if ( exit_flg ) {
	if (_om->isVerbosityAndPrint( Errors )) {
	  *_os << " Exiting Block CG iteration " << endl; 
	  *_os << " ERROR: No independent initial direction vectors" << endl;
	}
      }		
      if (num_ind < ind_blksz) {
	// The initial block of direction vectors are linearly dependent
	if (_om->isVerbosityAndPrint( Warnings )) {
	  *_os << " Initial direction vectors are dependent" << endl;
	  *_os << " Adjusting blocks and indices for iteration" << endl;
	}
	//
	// Adjust block and indices for iteration.
	//
	ind_blksz = num_ind;
	std::vector<int> temp_idx( ind_blksz );
	for (i=0; i< ind_blksz; i++)
	  temp_idx[i] = ind_idx[cols[i]];
	ind_idx.resize( ind_blksz );
	for (i=0; i< ind_blksz; i++)      
	  ind_idx[i] = temp_idx[i];
	
      }  // end if (num < ind_blksz)
    }  // end if (ind_blksz > 0)
    //		
    else {  // all initial residuals have converged
      if (_om->isVerbosityAndPrint( Warnings )) {
	*_os << " Exiting Block CG iteration " 
	     << " -- Iteration# " << _iter << endl;
	*_os << "Reason: All initial residuals have converged" << endl;
      }
      exit_flg = true;
    }
    
    // ***************************************************************************
    // ************************Main CG Loop***************************************
    // ***************************************************************************
    // 
    _new_blk = 1;
    for (_iter=0; _stest->CheckStatus(this) == Unconverged && !exit_flg; _iter++) {
      //
      //----------------Compute the new blocks of iterates and residuals------------------
      //
      // Get views of the previous blocks of residuals, direction vectors, etc.
      //
      if (_new_blk){
	P_prev = MVT::CloneView( *_basisvecs, ind_idx );
      }
      else {
	index2.resize( ind_blksz );
	for (i=0; i< ind_blksz; i++) {
	  index2[i] = _blocksize + ind_idx[i];
	} 
	P_prev = MVT::CloneView( *_basisvecs, index2 );
      }
      //
      index2.resize( _blocksize );
      for (i=0; i < _blocksize; i++){
	index2[i] = _blocksize + i;
      }
      if (_new_blk){
	R_prev = MVT::CloneView( *_residvecs, index );
	R_new = MVT::CloneView( *_residvecs, index2 );
      }
      else {
	R_prev = MVT::CloneView( *_residvecs, index2 );
	R_new = MVT::CloneView( *_residvecs, index );
      }
      //
      // Compute the coefficient matrix alpha
      //
      // P_prev^T * A * P_prev * alpha = P_prev^T * R_prev
      // 1) Compute P_prev^T * A * P_prev = T2 and P_prev^T * R_prev = T1
      // 2) Compute the Cholesky Factorization of T2
      // 3) Back and Forward Solves for alpha
      //
      // 1)
      if ( ind_blksz < prev_ind_blksz ) {  
	//
	// The number of independent direction vectors has changed,
	// so the dimension of the application multivectors needs to be resized.
	//
	AP_prev = MVT::CloneView( *temp_blk, ind_idx ); 
	alpha.reshape( ind_blksz, _blocksize );
	T2.reshape( ind_blksz, ind_blksz );
      }
      {
	Teuchos::TimeMonitor OpTimer(*_timerOp);
	_lp->ApplyOp( *P_prev, *AP_prev );
      }
      MVT::MvTransMv( one, *P_prev, *R_prev, alpha);   
      MVT::MvTransMv( one, *P_prev, *AP_prev, T2);
      //
      // 2)
      lapack.POTRF(UPLO, ind_blksz, T2.values(), ind_blksz, &info);
      if (info != 0) {
	if(_om->isVerbosityAndPrint( Errors )){
	  *_os << " Exiting Block CG iteration "
	       << " -- Iteration# " << _iter << endl;
	  *_os << " ERROR: Cannot compute coefficient matrix alpha" << endl;
	  *_os << " P_prev'* A*P_prev is singular" << endl;
	  *_os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      // 3)
      lapack.POTRS(UPLO, ind_blksz, _blocksize, T2.values(), ind_blksz, 
		   alpha.values(), ind_blksz, &info);
      // Note: solution returned in alpha
      if (info != 0) {
	if(_om->isVerbosityAndPrint( Errors  )){
	  *_os << " Exiting Block CG iteration "
	       << " -- Iteration# " << _iter << endl;
	  *_os << " ERROR: Cannot compute coefficient matrix alpha" << endl;
	  *_os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Update the solution: cur_block_sol = cur_block_sol + P_prev*alpha
      // 
      MVT::MvTimesMatAddMv( one, *P_prev, alpha, one, *_cur_block_sol );
      _lp->SolutionUpdated();
      //
      // Update the residual vectors: R_new = R_prev - A*P_prev*alpha
      //
      MVT::MvAddMv(one, *R_prev, zero, *R_prev, *R_new);
      MVT::MvTimesMatAddMv(-one, *AP_prev, alpha, one, *R_new);
      //
      // ****Compute the Current Relative Residual Norms and the Block Error****
      //
      MVT::MvNorm( *R_new, &_cur_resid_norms );
      //
      prev_ind_blksz = ind_blksz; // Save old ind_blksz of P_prev
      //
      // Update the number of current residuals that correspond
      // to linearly independent direction vectors. Note that
      // ind_idx are a subset of cur_idx.
      //
      k = 0;
      std::vector<int> temp_idx;
      for (i=0; i< ind_blksz; i++){
	if (_cur_resid_norms[ ind_idx[i] ] / _init_resid_norms[ ind_idx[i] ] > _prec)
	  temp_idx.push_back( ind_idx[i] );
      }
      ind_blksz = temp_idx.size();
      if (ind_blksz < prev_ind_blksz ) {
	ind_idx.resize( ind_blksz );
	for (i=0; i< ind_blksz; i++)
	  ind_idx[i] = temp_idx[i];
      }
      //
      // ****************Print iteration information*****************************
      //
      if (_om->isVerbosity( IterationDetails )) {
	PrintCGIterInfo( ind_idx );
      }
      //
      // ****************Test for breakdown*************************************
      //
      if (ind_blksz <= 0){
	if (_om->isVerbosityAndPrint( Errors )) {
	  *_os << " Exiting Block CG iteration " 
	       << " -- Iteration# " << _iter << endl;
	  *_os << " ERROR: No more independent direction vectors" << endl;
	  *_os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // **************Compute the new block of direction vectors****************
      //
      // Get views of the new blocks of independent direction vectors and
      // the corresponding residuals. 
      //
      index2.resize( ind_blksz );
      for (i=0; i<ind_blksz; i++){
	index2[i] = _blocksize + ind_idx[i];
      }
      
      if (_new_blk) {
	R_new = MVT::CloneView( *_residvecs, index2 );
	P_new = MVT::CloneView( *_basisvecs, index2 );
      }
      else {
	R_new = MVT::CloneView( *_residvecs, ind_idx );
	P_new = MVT::CloneView( *_basisvecs, ind_idx );
      }
      //
      // Put the current preconditioned initial residual into P_new 
      // since P_new = precond_resid + P_prev * beta
      //
      if ( isPrec ) {
	Teuchos::TimeMonitor PrecTimer(*_timerPrec);
	_lp->ApplyLeftPrec( *R_new, *P_new );
      } else { 
	MVT::MvAddMv( one, *R_new, zero, *R_new, *P_new ); 
      }
      // 
      // Compute coefficient matrix beta
      //
      // P_prev^T A * P_prev * beta = P_prev^T A * precond_resid
      // 1) Compute P_prev^T A * P_prev = T2 and P_prev^T * A * precond_resid = T3
      //                                     or (A*P_prev)^T * precond_resid (A SPD)
      // 2) Compute the Cholesky Factorization of T2
      // 3) Back and Forward Solves for beta
      //
      beta.reshape(prev_ind_blksz,ind_blksz);
      //
      // 1 & 2)  Note: we already have computed T2 and its Cholesky
      //         factorization during computation of alpha
      MVT::MvTransMv(-one, *AP_prev, *P_new, beta);
      // 3)
      lapack.POTRS(UPLO, prev_ind_blksz, ind_blksz, T2.values(), prev_ind_blksz, 
		   beta.values(), prev_ind_blksz, &info);
      // Note: Solution returned in beta
      if (info != 0) {
	if (_om->isVerbosityAndPrint( Errors )) {
	  *_os << " Exiting Block CG iteration " 
	       << " -- Iteration# " << _iter << endl;
	  *_os << "ERROR: Cannot compute coefficient matrix beta" << endl;
	  *_os << "Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Compute: P_new = precond_resid + P_prev * beta
      // ( remember P_new already = precond_resid )
      MVT::MvTimesMatAddMv(one, *P_prev, beta, one, *P_new);
      //
      // Check A-orthogonality of new and previous blocks of direction vectors
      //
      if (_om->isVerbosity( OrthoDetails )) {
	if (_om->doPrint()) *_os << "Orthogonality check" << endl;
	CheckAOrthogonality(*P_prev, *P_new);   
      }
      //
      // Compute orthonormal block of direction vectors,
      // and check for dependencies, adjusting indices of
      // independent vectors if needed
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> G(ind_blksz,ind_blksz);
      {
	Teuchos::TimeMonitor OrthoTimer(*_timerOrtho);
	exit_flg = QRFactorDef(*P_new, G, &cols);
      }
      //
      // Check if the orthogonalization has failed.
      //
      num_ind = cols.size();
      if ( exit_flg ) {
	ind_blksz = num_ind;
	if (_om->isVerbosityAndPrint( Errors )) {
	  *_os  << " Exiting Block CG iteration "  
		<< " -- Iteration# " << _iter << endl;
	  *_os << " ERROR: No more linearly independent direction vectors" << endl;
	  *_os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Check if the new block of direction vectors are linearly dependent
      //
      if (num_ind < ind_blksz) {
	if (_om->isVerbosityAndPrint( OrthoDetails )) {
	  *_os << "The new block of direction vectors are dependent " << endl;
	  *_os << "# independent direction vectors: " << num_ind << endl;
	  *_os << " Independent indices: " << endl;
	  for (i=0; i<num_ind ; i++)
	    *_os << cols[i] << " ";
	  *_os << endl << endl;
	}
	ind_blksz = num_ind;
	std::vector<int> temp_idx( ind_blksz );
	for (i=0; i<ind_blksz; i++)
	  temp_idx[i] = ind_idx[cols[i]];
	ind_idx.resize( ind_blksz );
	for (i=0; i<ind_blksz; i++)
	  ind_idx[i] = temp_idx[i];
      }
      //
      // Check A-orthogonality after orthonormalization
      //
      if (_om->isVerbosity( OrthoDetails )) {
	if (_om->doPrint())  *_os << "Orthogonality check after orthonormalization " << endl; 
	CheckAOrthogonality(*P_prev, *P_new);
      }
      // 
      // *****Update index of new blocks*****************************************
      //
      if (_new_blk)
	_new_blk--;
      else
	_new_blk++;
      //
    } // end of the main CG loop -- for(_iter = 0;...)
    //
    // Inform the linear problem manager that we are done with the current block of linear systems.
    //
    _lp->SetCurrLSVec();
    //
    // Print out solver status.
    //
    if (_om->isVerbosityAndPrint( FinalSummary )) {
      _stest->Print(*_os);
  }  
    //
  } // end BlockIteration()
  //
  
  template<class ScalarType, class MV, class OP>
  bool 
  BlockCG<ScalarType,MV,OP>::QRFactorDef (MV& VecIn, 
					  Teuchos::SerialDenseMatrix<int,ScalarType>& R, 
					  std::vector<int>* cols) 
  {
    int i, j, k;
    int num_orth, num_dep = 0;
    int nb = MVT::GetNumberVecs( VecIn );
    std::vector<int> index, dep_idx;
    std::vector<int> index2( 1 );
    RefCountPtr<MV> qj, Qj;
    Teuchos::SerialDenseVector<int,ScalarType> rj;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    std::vector<MagnitudeType> NormVecIn( nb );
    std::vector<MagnitudeType> qNorm( 1 );
    //
    // Zero out the array that will contain the Fourier coefficients.
    // 
    R.putScalar();
    //
    // Reset the size of the independent columns back to zero.
    //
    cols->resize( 0 );
    //
    // Compute the Frobenius norm of VecIn -- this will be used to
    // determine rank deficiency of VecIn
    //
    MVT::MvNorm( VecIn, &NormVecIn );
    ScalarType FroNorm = zero;
    for (j=0; j<nb; j++) {
      FroNorm += NormVecIn[j] * NormVecIn[j];
    }
    FroNorm = Teuchos::ScalarTraits<ScalarType>::squareroot(FroNorm);
    //
    // Start the loop to orthogonalize the nb columns of VecIn.
    //
    //
    for ( j=0; j<nb; j++ ) {
      //
      // Grab the j-th column of VecIn (the first column is indexed to 
      // be the zero-th one).
      //
      index2[0] = j;
      qj = MVT::CloneView( VecIn, index2 );
      if ( j ) {
	//
	// Grab the first j columns of VecIn (that are now an orthogonal
	// basis for first j columns of the entering VecIn).
	//
	index.resize( j );
	for (i=0; i<j; i++)
	  index[i] = i;
	Qj = MVT::CloneView( VecIn, index );
	rj.size(j);
	//
	// Enter a for loop that does two (num_orth) steps of classical 
	// Gram-Schmidt orthogonalization.
	//
	for (num_orth=0; num_orth<2; num_orth++) {
	  //
	  // Determine the Fourier coefficients for orthogonalizing column
	  // j of VecIn against columns 0:j-1 of VecIn. In other words, 
	  // result = trans(Qj)*qj.
	  //
	  MVT::MvTransMv( one, *Qj, *qj, rj );
	  //
	  // Sum result[0:j-1] into column j of R.
	  //
	  for (k=0; k<num_dep; k++) {
	    rj[dep_idx[k]] = zero;
	  }
	  //
	  for ( k=0; k<j; k++ ) {
	    R(k,j) += rj[k];
	  }
	  //
	  //   Compute qj <- qj - Qj * rj.
	  //
	  MVT::MvTimesMatAddMv(-one, *Qj, rj, one, *qj);
	}
	
      } // if (j)
      //
      // Compute the norm of column j of VecIn (=qj).
      //
      MVT::MvNorm( *qj, &qNorm );
      R(j,j) = qNorm[0];
      //
      if ( NormVecIn[j] > _prec && SCT::magnitude(R(j,j)) > (_prec * NormVecIn[j]) ) {
	//
	// Normalize qj to make it into a unit vector.
	//
	ScalarType rjj = one / R(j,j);
	MVT::MvAddMv( rjj, *qj, zero, *qj, *qj );
	cols->push_back( j );  
      }
      else {
	// 
	if (_om->isVerbosityAndPrint( OrthoDetails )){
	  *_os << "Rank deficiency at column index: " << j << endl;
	}
	//
	// Don't normalize qj, enter one on diagonal of R,
	// and zeros in the row to the right of the diagonal -- this
	// requires updating the indices of the dependent columns
	//
	R(j,j) = one;
	dep_idx.push_back( j );
	num_dep++;
      }	
      
    } // for (j=0; j<nb; j++)
    //
    // Return true if we could not create any independent direction vectors (failure).
    //
    if (cols->size() == 0) return true;
    return false;
    //
  } // end QRFactorDef
  
  
  template<class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::CheckOrthogonality(MV& P1, MV& P2) 
  {
    //
    // This routine computes P2^T * P1
    //
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    int i, k;
    //
    int numvecs1 = MVT::GetNumberVecs( P1 );
    int numvecs2 = MVT::GetNumberVecs( P2 );
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> PtP(numvecs2, numvecs1);
    MVT::MvTransMv( one, P2, P1, PtP);
    //
    ScalarType* ptr = PtP.values();
    ScalarType column_sum;
    //
    for (k=0; k<numvecs1; k++) {
      column_sum = zero;
      for (i=0; i<numvecs2; i++) {
	column_sum += ptr[i];
      }
      *_os << " P2^T*P1 for column "
	   << k << " is  " << Teuchos::ScalarTraits<ScalarType>::magnitude(column_sum) << endl;
      ptr += numvecs2;
    }
    *_os << endl;
    //
  } // end check_orthog
  //
  
  template<class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::CheckAOrthogonality(MV& P1, MV& P2) 
  {
    //
    // This routine computes P2^T * A * P1
    // Checks the orthogonality wrt A between any two blocks of multivectors with the same length
    //
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    int i, k;
    //
    int numvecs1 = MVT::GetNumberVecs( P1 );
    int numvecs2 = MVT::GetNumberVecs( P2 );
    //
    RefCountPtr<MV> AP = MVT::Clone( P1, numvecs1 );
    {
      Teuchos::TimeMonitor OpTimer(*_timerOp);
      _lp->ApplyOp(P1, *AP);
    }
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> PAP(numvecs2, numvecs1);
    MVT::MvTransMv( one, P2, *AP, PAP);
    //
    ScalarType* ptr = PAP.values();
    ScalarType column_sum;
    //
    for (k=0; k<numvecs1; k++) {
      column_sum = zero;
      for (i=0; i<numvecs2; i++) {
	column_sum += ptr[i];
    }
      *_os << " P2^T*A*P1 for column "
	   << k << " is  " << Teuchos::ScalarTraits<ScalarType>::magnitude(column_sum) << endl;
      ptr += numvecs2;
    }
    *_os << " " << endl;
    //
  } // end check_orthog
  
  
  template<class ScalarType, class MV, class OP>
  void 
  BlockCG<ScalarType,MV,OP>::PrintCGIterInfo( const std::vector<int> &ind )
  {
    //
    int i;
    int indsz = ind.size();
    *_os << "# of independent direction vectors: " << indsz << endl;    
    *_os << " Independent indices: " << endl;
    for (i=0; i<indsz; i++){
      *_os << ind[i] << " ";
    }
    *_os << endl << endl;
    //
  } // end Print_CGiter_info

} // end namespace Belos

#endif
// End of file BelosBlockCG.hpp


