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

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"

/*!	
  \class Belos::BlockGmres
  
  \brief This class implements the Restarted Block GMRES algorithm
  for solving real nonsymmetric linear systems of equations AX = B,
  where B is a matrix containing one or more right-hand sides, and 
  X is the matrix of corresponding solutions.
  
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
  template <class TYPE>
  class BlockGmres : public IterativeSolver<TYPE> { 
  public:
    //@{ \name Constructor/Destructor.
    //! %Belos::BlockGmres constructor.
    BlockGmres(LinearProblemManager<TYPE>& lp, 
	       StatusTest<TYPE>& stest,
	       const int numrhs, 
	       const TYPE tol=1.0e-6, 
	       const int length=25, 
	       const int block=1, 
	       bool vb = false);
    
    //! %Belos::BlockGmres destructor.
    virtual ~BlockGmres();
    //@}
    
    //@{ \name Accessor methods
    
    //! Get the iteration count for the current block of linear systems.
    int GetNumIters() const { return( _totaliter ); }
    
    //! Get the restart count of the iteration method for the current block of linear systems.
    int GetNumRestarts() const { return( _restartiter ); }
    
    //! Get the current number of linear system being solved for.
    /*! Since the block size is independent of the number of right-hand sides, 
      it is important to know how many linear systems
      are being solved for when the status is checked.  This is informative for residual
      checks because the entire block of residuals may not be of interest.  Thus, this 
      number can be anywhere between 1 and the block size for the solver.
    */
    int GetNumToSolve() const { return( _num_to_solve ); };
    
    //! Get the current number of the right-hand side block being solved for.
    /*! Since the block size is independent of the number of right-hand sides for
      some solvers (GMRES, CG, etc.), it is important to know which right-hand side
      block is being solved for.  That may mean you need to update the information
      about the norms of your initial residual vector for weighting purposes.  This
      information can keep you from querying the solver for information that rarely
      changes.
    */
    int GetRHSBlockNum() const { return( _rhs_iter ); };
    
    //! Get the solvers native residuals for the current block of linear systems.
    /*! For GMRES this is not the same as the actual residual of the linear system.
        \note If the blocksize is one, then the only available norm for the native
	residual is the TwoNorm, any other norm will return undefined.
     */
    ReturnType GetNativeResidNorms(TYPE *normvec, NormType norm_type) const;

    //! Get the true residuals for the current block of linear systems.
    /*! For GMRES this will force the solver to compute a current residual for its linear 
      systems, the current solution is not stored. <b> This is an expensive computation, 
      so a convergence test using these residuals should be secondary to using the native 
      residuals. </b>
    */
    MultiVec<TYPE>* GetTrueResidVecs();

    //! Get a constant reference to the current linear problem.  
    /*! This may include a current solution, if the solver has recently restarted or completed.
     */
    const LinearProblemManager<TYPE>& GetLinearProblem() const { return( _lp ); }

    /*! \brief Get information whether the solution contained in the linear problem
      is current.  In the case of GMRES, this solution is only updated during a restart
      or when the solver completes.
    */
    bool IsSolutionCurrent() { return (_iter == 0); };
    
    //@} 

    //@{ \name Solver application method.
    
    /*! \brief This method uses the iterative method to compute approximate
      solutions to the original problem.  This method can return unconverged if the
      maximum number of iterations is reached, or numerical breakdown is observed.
    */
    void Solve(bool);
    //@}
    
    //@{ \name Set methods.
    
    //! This method sets the number of allowable restarts.
    void SetRestart( const int restart ) { _restart = (restart > 0 ? restart : 0 ); } 
    //@}
    
    //@{ \name Output methods.
    
    /*! \brief This method allows the user to set the solver's level of visual output
      during computations.
    */
    void SetDebugLevel( const int debuglevel ) { _debuglevel = debuglevel; }
    
    //! This method requests that the solver print out its current residuals.
    void PrintResids(bool) const;

    //@}

  private:

    //! Method for setting the basis dependency tolerances for extending the Krylov basis.
    void SetGmresBlkTols();

    /*! Method for augmenting the solution / right-hand side block when the number of linear problems
      being solved for is less than the blocksize
    */
    void SetUpBlocks();

    /*! Method for extracting out the solution from the solver when the number of linear problems
      being solved for is less than the blocksize.
    */
    void ExtractCurSolnBlock();

    //! Method for performing the block Krylov decomposition.
    bool BlockReduction(bool&, bool);

    //! Method for orthogonalization of one block.
    bool QRFactorAug(MultiVec<TYPE>&, Teuchos::SerialDenseMatrix<int,TYPE>&,
		     bool, bool);

    //! Method for block orthogonalization when a dependency has not been detected in the Krylov basis.
    bool BlkOrth(MultiVec<TYPE>&, bool);

    //! Method for block orthogonalization when a dependency has been detected in the Krylov basis.
    bool BlkOrthSing(MultiVec<TYPE>&, bool);

    //! Method for checking the orthogonality of the Krylov basis.
    void CheckKrylovOrth(const int, bool);

    //! Reference to the linear problem being solver for with the solver. [passed in by user]
    LinearProblemManager<TYPE>& _lp; 

    //! Reference to the status test, which provides the stopping criteria for the solver. [passed in by user]
    StatusTest<TYPE>& _stest; 

    //! Pointers to the right-hand side and solution multivecs held by the linear problem manager.
    MultiVec<TYPE> *_rhs, *_solutions; 

    //! Pointers to the Krylov basis constructed by the solver.
    MultiVec<TYPE> *_basisvecs;

    //! Pointers to the current right-hand side and solution multivecs being solved for.
    MultiVec<TYPE> *_cur_block_rhs, *_cur_block_sol;

    //! Dense matrices for holding the upper Hessenberg matrix (H) of the Arnoldi factorization 
    Teuchos::SerialDenseMatrix<int,TYPE> _hessmatrix;

    //! Dense vector for holding the right-hand side of the least squares problem.
    Teuchos::SerialDenseMatrix<int,TYPE> _z;

    const int _length, _blocksize, _numrhs;
    const TYPE _residual_tolerance;
    TYPE *_residerrors; 
    int _num_to_solve;
    int _rhs_iter, _restartiter, _totaliter, _iter;
    int _debuglevel, _restart;
    TYPE _blkerror, _dep_tol, _blk_tol, _sing_tol;
  };
  //
  // Implementation
  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class TYPE>
  BlockGmres<TYPE>::BlockGmres(LinearProblemManager<TYPE>& lp, 
			       StatusTest<TYPE>& stest,
			       const int numrhs, 
			       const TYPE tol, 
			       const int length, 
			       const int block, 
			       bool vb) : 
    _lp(lp),
    _stest(stest),
    _rhs(_lp.GetRHS()), 
    _solutions(_lp.GetLHS()),
    _basisvecs(0),     
    _cur_block_rhs(0),
    _cur_block_sol(0),
    _length(length), 
    _blocksize(block), 
    _numrhs(numrhs), 
    _residual_tolerance(tol), 
    _residerrors(0), 
    _num_to_solve(block),
    _rhs_iter(0), 
    _restartiter(0), 
    _totaliter(0),
    _iter(0),
    _debuglevel(0), 
    _restart(0),
    _blkerror(1.0)
  {
    //
    // Initial check that input information is valid.
    //
    assert(_length > 0); assert(_blocksize > 0); assert(_numrhs > 0);
    //
    // Make room for the Arnoldi vectors and F.
    //
    _basisvecs = _rhs->Clone((_length+1)*_blocksize); assert(_basisvecs!=NULL);
    //
    // Create the rectangular Hessenberg matrix and right-hand side of least squares problem.
    //
    _hessmatrix.shape((_length+1)*_blocksize, _length*_blocksize);
    _z.shape((_length+1)*_blocksize, _blocksize); 
    _residerrors = new TYPE[_numrhs + _blocksize]; assert(_residerrors!=NULL);
    //
    // Set up the block orthogonality tolerances
    //
    SetGmresBlkTols();	
  }
    
  template <class TYPE>
  BlockGmres<TYPE>::~BlockGmres() 
  {
    if (_basisvecs) delete _basisvecs;
    if (_residerrors) delete [] _residerrors;
  }
  
  template <class TYPE>
  void BlockGmres<TYPE>::SetGmresBlkTols() 
  {
    const TYPE two = 2.0;
    TYPE eps;
    char precision = 'P';
    Teuchos::LAPACK<int,TYPE> lapack;
    eps = lapack.LAMCH(precision);
    _dep_tol = 1/sqrt(two);
    _blk_tol = 10*sqrt(eps);
    _sing_tol = 10 * eps;
  }
  
  template <class TYPE>
  ReturnType BlockGmres<TYPE>::GetNativeResidNorms(TYPE *normvec, NormType norm_type) const 
  {
    if (normvec) {
      if (norm_type == TwoNorm) {
	Teuchos::BLAS<int,TYPE> blas;
	for (int j=0; j<_blocksize; j++) {
	  normvec[j] = blas.NRM2( _blocksize, &_z(_iter*_blocksize, j ), 1);
	}
	return Ok;
      } else
	return Undefined;
    }
    return Error;
  }
  
  template <class TYPE>
  MultiVec<TYPE>* BlockGmres<TYPE>::GetTrueResidVecs()
  {    
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    int i, m = _iter*_blocksize;
    Teuchos::BLAS<int,TYPE> blas;
    int *index = new int[m]; assert(index!=NULL);
    for ( i=0; i<m; i++ ) {   
      index[i] = i;
    }
    MultiVec<TYPE> * Vjp1 = _basisvecs->CloneView(index, m); assert(Vjp1!=NULL);
    MultiVec<TYPE> * RV = _basisvecs->Clone(_blocksize); assert(RV!=NULL);
    MultiVec<TYPE> * cur_sol_copy = _cur_block_sol->CloneCopy(); assert(cur_sol_copy!=NULL);
    //
    //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
    //
    Teuchos::SerialDenseMatrix<int,TYPE> y( Teuchos::Copy, _z, m, _blocksize );
    //
    //  Solve the least squares problem and compute current solutions.
    //
    blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	       Teuchos::NON_UNIT_DIAG, m, _blocksize, one,  
	       _hessmatrix.values(), _hessmatrix.stride(), y.values(), y.stride() );
    
    cur_sol_copy->MvTimesMatAddMv( one, *Vjp1, y, one );
    _lp.ComputeResVec( RV, cur_sol_copy, _cur_block_rhs );
    
    if (Vjp1) delete Vjp1;
    if (cur_sol_copy) delete cur_sol_copy;
    if (index) delete [] index;

    return RV;
  }
  

  template<class TYPE>
  void BlockGmres<TYPE>::SetUpBlocks ()  
  {
    //
    // We only need to copy the solutions and right-hand sides if the linear systems of
    // interest are less than the block size.
    //    
    if (_num_to_solve < _blocksize ) {
      int i;
      int *index = new int[_num_to_solve]; assert(index!=NULL);
      for ( i=0; i<_num_to_solve; i++) { index[i] = _rhs_iter*_blocksize+i; }
      const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
      //
      // Logic to handle the number of righthand sides solved for at this iteration,
      // when it is less than the blocksize of the solver.  In this case, the right-hand 
      // side and solution block have to be augmented with some random vectors.
      //
      // Fill up the right-hand side block with random vectors, then place
      // the remaining (unsolved) right hand sides into the initial portion
      // of the right-hand side block.
      //	
      int *index2 = new int[_num_to_solve]; assert(index!=NULL);
      _cur_block_rhs->MvRandom();
      //
      MultiVec<TYPE> *tptr = _rhs->CloneView(index,_num_to_solve); assert(tptr!=NULL);
      for (i=0; i<_num_to_solve; i++) {
	index2[i] = i;
      }
      _cur_block_rhs->SetBlock( *tptr, index2, _num_to_solve);
      //
      // Now deal with solution block, augment with zero vectors.
      //
      _cur_block_sol->MvInit( zero );
      MultiVec<TYPE> *tptr2 = _solutions->CloneView(index,_num_to_solve); assert(tptr2!=NULL);
      _cur_block_sol->SetBlock( *tptr2, index2, _num_to_solve);
      //
      // Delete the temporary views
      //
      delete tptr; tptr = 0;
      delete tptr2; tptr2 = 0;
      delete [] index2; index2=0;
      delete [] index; index=0;
    }
    //
  }
  

  template <class TYPE>
  void BlockGmres<TYPE>::ExtractCurSolnBlock() 
  {
    //
    // We only need to copy the solutions back in if the linear systems of
    // interest are less than the block size.
    //
    if (_num_to_solve < _blocksize) {
      int i;
      //
      int * index = new int[_num_to_solve]; assert(index!=NULL);
      MultiVec<TYPE> *tptr=0;
      //
      // Get a view of the current solutions and correction vector.
      //
      for (i=0; i<_num_to_solve; i++) { 
	index[i] = i;	
      }
      tptr = _cur_block_sol->CloneView(index,_num_to_solve); assert(tptr!=NULL);
      //
      // Copy the correction vector to the solution vector.
      //
      for (i=0; i< _num_to_solve; i++) { 
	index[i] = _rhs_iter*_blocksize+i; 
      }
      _solutions->SetBlock( *tptr, index, _num_to_solve);
      //
      // Clean up.
      //
      delete tptr;
      delete [] index; index=0;
    }
  }    
  
  template <class TYPE>   
  void BlockGmres<TYPE>::PrintResids(bool vb) const 
  {
    //
    int i;
    //
    if (vb) {
      cout << "--------------------------------------" << endl;
      for (i=0; i<_numrhs; i++){
	cout << "_residerrors[" << i << "] = " << _residerrors[i] << endl;
      }
      cout << endl;
    }
    //
  } // end PrintResids
  
  
  template <class TYPE>
  void BlockGmres<TYPE>::Solve (bool vb) 
  {
    int i,j, maxidx;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    TYPE sigma, mu, vscale, maxelem;
    MultiVec<TYPE> *U_vec=0, *F_vec=0;
    TYPE *beta = new TYPE[(_length+1)*_blocksize];
    int *index = new int[ (_length+1)*_blocksize ]; assert(index!=NULL);
    bool dep_flg = false, exit_flg = false, brkflg = false;
    Teuchos::LAPACK<int, TYPE> lapack;
    Teuchos::BLAS<int, TYPE> blas;
    //	
    // Each pass through the solver solves for _blocksize right-hand sides.
    // max_rhs_iters is the number of passes through the
    // solver required to solve for all right-hand sides	
    //
    int max_rhs_iters = (_numrhs+_blocksize-1) / _blocksize;
    for (i=0; i < (_length+1)*_blocksize; i++) { index[i] = i; }
    //
    //  Start executable statements. 
    //
    for ( _rhs_iter=0; _rhs_iter < max_rhs_iters; _rhs_iter++ ) {
      //
      if (_debuglevel > 3 && vb) {
	cout << "_rhs_iter: " << _rhs_iter <<  endl
	     <<  endl;
      }
      brkflg = false;
      //
      // Compute the the remaining number of right-hand sides to be solved.
      //
      if ( _numrhs - (_rhs_iter * _blocksize) < _blocksize ) {
	_num_to_solve = _numrhs - (_rhs_iter * _blocksize);
      }
      //
      // Reset the iteration counter for this block of right-hand sides.
      //
      _totaliter = 0;
      //
      // Put the current initial guesses and right-hand sides into current blocks
      //
      if ( _num_to_solve < _blocksize ) {
	_cur_block_sol = _solutions->Clone(_blocksize);
	_cur_block_rhs = _solutions->Clone(_blocksize);
	SetUpBlocks();
      } else {
	_cur_block_sol = _solutions->CloneView( index+(_rhs_iter*_blocksize), _blocksize);
	_cur_block_rhs = _rhs->CloneView( index+(_rhs_iter*_blocksize), _blocksize);
      }
      //
      _blkerror = one;
      TYPE init_norm = one;
      //
      for (_restartiter=0; _restartiter < _restart+1 && _blkerror > _residual_tolerance ; _restartiter++) {
	//
	if (brkflg){
	  break;
	}
	//		
	// Associate the initial block of _basisvecs with U_vec.
	//
	U_vec = _basisvecs->CloneView(index, _blocksize);
	assert(U_vec!=NULL);
	//
	// Compute current residual and place into 1st block
	//
	_lp.ComputeResVec( U_vec, _cur_block_sol, _cur_block_rhs );
	//
	dep_flg = false; exit_flg = false;
	//
	// Re-initialize RHS of the least squares system and create a view.
	//
	_z.putScalar();
	Teuchos::SerialDenseMatrix<int,TYPE> G10(Teuchos::View, _z, _blocksize, _blocksize);
	exit_flg = QRFactorAug( *U_vec, G10, true, vb );
	//
	if (exit_flg){
	  if (vb){
	    cout << "Exiting Block GMRES" << endl;
	    cout << "  RHS pass# " << _rhs_iter+1
		 << "  Restart iteration# " << _restartiter
		 << "  Iteration# " << _iter << endl;
	    cout << "  Reason: Failed to compute initial block of orthonormal basis vectors"
		 << endl << endl;
	  }
	  if (U_vec) {delete U_vec; U_vec = 0;}
	  break;
	}
	//
	TYPE norm_G10 = G10.normFrobenius();  
	if (_restartiter == 0) {
	  init_norm = norm_G10;
	}
	//
	// The block error used here is an average residual error over
	// the current block. This is set to one to start with since all 
	// initial residuals have norm one
	//
	_blkerror = one; 	
	//
	for (_iter=0; _iter<_length && _blkerror > _residual_tolerance; _iter++, ++_totaliter) {
	  //
	  // Compute a length _length block Arnoldi Reduction
	  //    (one step at a time)
	  //
	  //dep_flg = true;
	  exit_flg = false;
	  exit_flg = BlockReduction(dep_flg, vb);
	  //
	  if (exit_flg){
	    brkflg = true; // set flag so we can also break out of the restart loop
	    break;
	  }
	  if (_stest.CheckStatus(this) != Converged) {
	    //_stest.Print(cout);
	  }
	  //
	  // QR factorization of Least-Squares system with Householder reflectors
	  //
	  for (j=0; j<_blocksize; j++) {
	    //
	    // Apply previous Householder reflectors to new block of Hessenberg matrix
	    //
	    for (i=0; i<_iter*_blocksize+j; i++) {
	      sigma = blas.DOT( _blocksize, &_hessmatrix(i+1,i), 1, &_hessmatrix(i+1,_iter*_blocksize+j), 1);
	      sigma += _hessmatrix(i,_iter*_blocksize+j);
	      sigma *= beta[i];
	      blas.AXPY(_blocksize, -sigma, &_hessmatrix(i+1,i), 1, &_hessmatrix(i+1,_iter*_blocksize+j), 1);
	      _hessmatrix(i,_iter*_blocksize+j) -= sigma;
	    }
	    //
	    // Compute new Householder reflector
	    //
	    maxidx = blas.IAMAX( _blocksize+1, &_hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j), 1 );
	    maxelem = _hessmatrix(_iter*_blocksize+j+maxidx-1,_iter*_blocksize+j);
	    for (i=0; i<_blocksize+1; i++) 
	      _hessmatrix(_iter*_blocksize+j+i,_iter*_blocksize+j) /= maxelem;
	    sigma = blas.DOT( _blocksize, &_hessmatrix(_iter*_blocksize+j+1,_iter*_blocksize+j), 1, 
			      &_hessmatrix(_iter*_blocksize+j+1,_iter*_blocksize+j), 1 );
	    if (sigma == zero) {
	      beta[_iter*_blocksize + j] = zero;
	    } else {
	      mu = sqrt(_hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j)*_hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j)+sigma);
	      if ( _hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j) < zero ) {
		vscale = _hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j) - mu;
	      } else {
		vscale = -sigma / (_hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j) + mu);
	      }
	      beta[_iter*_blocksize+j] = 2.0*vscale*vscale/(sigma + vscale*vscale);
	      _hessmatrix(_iter*_blocksize+j,_iter*_blocksize+j) = maxelem*mu;
	      for (i=0; i<_blocksize; i++)
		_hessmatrix(_iter*_blocksize+j+1+i,_iter*_blocksize+j) /= vscale;
	    }
	    //
	    // Apply new Householder reflector to rhs
	    //
	    for (i=0; i<_blocksize; i++) {
	      sigma = blas.DOT( _blocksize, &_hessmatrix(_iter*_blocksize+j+1,_iter*_blocksize+j), 1, &_z(_iter*_blocksize+j+1,i), 1);
	      sigma += _z(_iter*_blocksize+j,i);
	      sigma *= beta[_iter*_blocksize+j];
	      blas.AXPY(_blocksize, -sigma, &_hessmatrix(_iter*_blocksize+j+1,_iter*_blocksize+j), 1, &_z(_iter*_blocksize+j+1,i), 1);
	      _z(_iter*_blocksize+j,i) -= sigma;
	    }
	  }
	  //
	  // Compute the residuals and the block error
	  //
	  _blkerror = zero;
	  for (j=0; j<_blocksize; j++ ) {
	    _residerrors[_rhs_iter*_blocksize+j] = blas.NRM2( _blocksize, &_z((_iter+1)*_blocksize, j ), 1);
	    //
	    if (norm_G10) _residerrors[_rhs_iter*_blocksize+j] /= init_norm; 
	    //
	    _blkerror += _residerrors[_rhs_iter*_blocksize+j];
	  }
	  _blkerror = _blkerror / _blocksize;
	  //
	  // Print out residuals
	  //
	  if (_debuglevel>0 && vb) {
	    cout << " " << endl;
	    cout << "------------------------------------------------------------------------" << endl;
	    cout << "Computed GMRES Residual norms -- " << endl;
	    cout << "  RHS pass# " << _rhs_iter+1 
		 << "  Restart iteration# " << _restartiter 
		 << "  Iteration# " << _iter << endl;
	    for (j=0; j<_blocksize; j++) {
	      cout << "  _residerrors[" << _rhs_iter*_blocksize+j << "] = " << 
		_residerrors[_rhs_iter*_blocksize+j] << endl;
	    }
	    cout << " " << endl;
	  }
	  //
	} // end for (_iter=0;...
	//
	// Update the solutions by solving the triangular system to get the Krylov weights.
	//
	blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, _iter*_blocksize, _blocksize, one,
		   _hessmatrix.values(), _hessmatrix.stride(), _z.values(), _z.stride() ); 
	// Create view into current basis vectors.
	MultiVec<TYPE> * Vjp1 = _basisvecs->CloneView(index, _iter*_blocksize);
	Teuchos::SerialDenseMatrix<int,TYPE> rhs_view(Teuchos::View, _z, _iter*_blocksize, _blocksize);
	_cur_block_sol->MvTimesMatAddMv( one, *Vjp1, rhs_view, one );
	delete Vjp1;
	//
	// Insert the current block of solutions into _solutions and let the
	// linear problem know we updated to solution, so we can
	// continue if we have more right-hand sides to solve for
	//	  
	ExtractCurSolnBlock();
	_lp.SolutionUpdated();
	//
	// Print out solver status
	// 
	if (vb) {
	  cout << "  RHS pass# " << _rhs_iter+1 
	       << "  Restart iteration# " << _restartiter 
	       << "  Iteration# " << _iter << endl;
	  if (exit_flg) {
	    cout << " Exiting Block GMRES --- " << endl;
	    cout << "  Reason: Failed to compute new block of orthonormal basis vectors" << endl;
	    cout << "  ***Solution from previous step will be returned***"<< endl<< endl;
	  }
	  if (_restartiter == _restart && _blkerror > _residual_tolerance) {
	    cout << " Exiting Block GMRES --- " << endl;
	    cout << "  Reason: maximum number of iterations has been reached"
		 << endl << endl;
	  }
	  if (_blkerror < _residual_tolerance) {
	    cout << " Exiting Block GMRES --- " << endl;
	    cout << "  Reason: Block GMRES has converged" << endl << endl;
	  }
	} 
	if (U_vec) {delete U_vec; U_vec = 0;}
	//
      } // end for (_restartiter=0;...
      //
      // **************Free heap space**************
      //
      if (_cur_block_sol) {delete _cur_block_sol; _cur_block_sol=0;}
      if (_cur_block_rhs) {delete _cur_block_rhs; _cur_block_rhs=0;}
      //
    } // end for (_rhs_iter =0;...)
    //
    if (index) {delete [] index; index=0;}
    if (beta) {delete [] beta; beta=0; }
    //
  } // end Solve()
  
    
  template<class TYPE>
  bool BlockGmres<TYPE>::BlockReduction ( bool& dep_flg, bool vb ) 
  {
    //
    int i;	
    int *index = new int[_blocksize]; assert(index!=NULL);
    MultiVec<TYPE> *AU_vec = _basisvecs->Clone( _blocksize );
    //
    // Associate the j-th block of _basisvecs with U_vec.
    //
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = _iter*_blocksize+i;
    }
    MultiVec<TYPE>* U_vec = _basisvecs->CloneView(index, _blocksize);
    assert(U_vec!=NULL);
    //
    _lp.Apply( *U_vec, *AU_vec ); 
    //
    bool dep = false;
    if (!dep_flg){
      dep = BlkOrth(*AU_vec, vb);
      if (dep) {
	dep_flg = true;
      }
    }
    // If any dependencies have been detected during this step of
    // Block Reduction, or any previous steps (within the construction
    // of the current Krylov subspaces), block orthogonalization is 
    // implemented with a variant of A. Ruhe's approach.
    //
    bool flg = false;
    if (dep_flg){
      flg = BlkOrthSing(*AU_vec, vb);
    }
    //
    delete U_vec; delete [] index;
    delete AU_vec;
    //
    return flg;
    //
  } // end BlockReduction()
  
  
  template<class TYPE>
  bool BlockGmres<TYPE>::BlkOrth( MultiVec<TYPE>& VecIn, bool vb) 
  {
    //
    // Orthogonalization is first done between the new block of 
    // vectors and all previous blocks, then the vectors within the
    // new block are orthogonalized.
    //
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    const int max_num_orth = 2;
    int i, k, row_offset, col_offset;
    int * index = new int[(_length+1)*_blocksize]; assert(index!=NULL);
    TYPE * norm1 = new TYPE[_blocksize]; assert(norm1!=NULL);
    TYPE * norm2 = new TYPE[_blocksize]; assert(norm2!=NULL);
    //
    // Initialize index vector.
    //
    for (i=0; i<(_iter+2)*_blocksize; i++) { index[i] = i; }
    //
    // Associate (j+1)-st block of ArnoldiVecs with F_vec.
    //
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index+(_iter+1)*_blocksize, _blocksize);
    assert(F_vec!=NULL);
    //
    // Copy preconditioned AU_vec into (j+1)st block of _basisvecs
    //
    F_vec->MvAddMv( one, VecIn, zero, VecIn);
    //
    // Zero out the full block column of the Hessenberg matrix 
    // even though we're only going to set the coefficients in 
    // rows [0:(j+1)*_blocksize-1]
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_blocksize; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix(i,_iter*_blocksize+k) = zero;
      }
    }
    //
    // Grab all previous Arnoldi vectors
    //
    int num_prev = (_iter+1)*_blocksize;
    MultiVec<TYPE>* V_prev = _basisvecs->CloneView(index, num_prev);
    assert(V_prev!=NULL);
    //
    // Create a matrix to store the product trans(V_prev)*F_vec
    //
    Teuchos::SerialDenseMatrix<int,TYPE> dense_mat( num_prev, _blocksize );
    //
    F_vec->MvNorm(norm1);
    //
    // Perform two steps of block classical Gram-Schmidt so that
    // F_vec is orthogonal to the columns of V_prev.
    //
    for ( int num_orth=0; num_orth<max_num_orth; num_orth++ ) {
      //
      // Compute trans(V_prev)*F_vec and store in the j'th diagonal
      // block of the Hessenberg matrix
      //
      F_vec->MvTransMv (one, *V_prev, dense_mat);
      //
      // Update the orthogonalization coefficients for the j-th block
      // column of the Hessenberg matrix.
      //
      for ( k=0; k<_blocksize; k++ ) {
	for ( i=0; i<num_prev; i++ ) {
	  _hessmatrix(i,_iter*_blocksize+k) += dense_mat(i,k);
	}
      }
      //
      // F_vec <- F_vec - V(0:(j+1)*block-1,:) * H(0:num_prev-1,j:num_prev-1)
      //
      F_vec->MvTimesMatAddMv( -one, *V_prev, dense_mat, one );
    } // end for num_orth=0;...)
      //
    F_vec->MvNorm(norm2);
    //
    // Check to make sure the new block of Arnoldi vectors are 
    // not dependent on previous Arnoldi vectors
    //
    bool flg = false; // This will get set true if dependencies are detected
    //
    for (i=0; i<_blocksize; i++){
      if (norm2[i] < norm1[i] * _blk_tol) {
	flg = true;
	if (_debuglevel > 3 && vb){
	  cout << "Col " << num_prev+i << " is dependent on previous "
	       << "Arnoldi vectors in V_prev" << endl;
	  cout << endl;
	}
      }
    } // end for (i=0;...)
      //
    if (_debuglevel>2) {
      if (vb){
	cout << endl;
	cout << "Checking Orthogonality after BlkOrth()"
	     << " Iteration: " << _iter << endl;
      }
      CheckKrylovOrth(_iter, vb);
    }
    //
    // If dependencies have not already been detected, compute
    // the QR factorization of the next block. Otherwise,
    // this block of Arnoldi vectors will be re-computed via and 
    // implementation of A. Ruhe's block Arnoldi.
    //
    if (!flg) {
      //
      // Compute the QR factorization of F_vec
      //
      row_offset = (_iter+1)*_blocksize; col_offset = _iter*_blocksize;
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _blocksize, _blocksize,
							  row_offset, col_offset);
      flg = QRFactorAug( *F_vec, sub_block_hess, false, vb );
    }
    //
    delete F_vec;
    delete V_prev;
    delete [] index;
    delete [] norm1;
    delete [] norm2;
    //
    return flg;
    //
  }  // end BlkOrth()
  
  
  template<class TYPE>
  bool BlockGmres<TYPE>::BlkOrthSing( MultiVec<TYPE>& VecIn, bool vb) 
  {
    //
    // This is a variant of A. Ruhe's block Arnoldi
    // The orthogonalization of the vectors AU_vec is done
    // one at a time. If a dependency is detected, a random
    // vector is added and orthogonalized against all previous
    // Arnoldi vectors.
    // 
    const int IntOne = 1;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::SerialDenseVector<int,TYPE> dense_vec;
    int i, k, num_orth;
    int * index = new int[(_length+1)*_blocksize]; assert(index!=NULL);
    int * index2 = new int[IntOne];
    TYPE nm1[IntOne];
    TYPE nm2[IntOne];
    //
    // Initialize index vector.
    //
    for ( i=0; i<(_length+1)*_blocksize; i++ ) { index[i] = i; }
    //
    // Associate (j+1)-st block of ArnoldiVecs with F_vec.
    //
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index+(_iter+1)*_blocksize, _blocksize);
    assert(F_vec!=NULL);
    //
    // Copy preconditioned AU_vec into (j+1)st block of _basisvecs
    //
    F_vec->MvAddMv( one, VecIn, zero, VecIn);
    //
    // Zero out the full block column of the Hessenberg matrix 
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_blocksize; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix(i, _iter*_blocksize+k) = zero;
      }
    }
    //
    MultiVec<TYPE> *q_vec=0, *Q_vec=0, *tptr=0;
    tptr = F_vec->Clone(IntOne); assert(tptr!=NULL);
    //
    // Start a loop to orthogonalize each of the _blocksize
    // columns of F_vec against all previous _basisvecs
    //
    bool flg = false;
    //
    for (int num_prev = (_iter+1)*_blocksize; num_prev < (_iter+2)*_blocksize; num_prev++) {
      // Initialize dense vector.
      dense_vec.size(num_prev);
      //
      // Grab the next column of _basisvecs
      //
      index2[0] = num_prev;
      q_vec = _basisvecs->CloneView(index2, IntOne); assert(q_vec!=NULL);
      //
      // Grab all previous columns of _basisvecs
      //
      Q_vec = _basisvecs->CloneView(index, num_prev); assert(Q_vec!=NULL);
      //
      // Do one step of classical Gram-Schmidt orthogonalization
      // with a 2nd correction step if needed.
      //
      bool dep = false;
      q_vec->MvNorm(nm1);
      //
      // Compute trans(Q_vec)*q_vec
      //
      q_vec->MvTransMv(one, *Q_vec, dense_vec);
      //
      // Sum results [0:num_prev-1] into column (num_prev-_blocksize)
      // of the Hessenberg matrix
      //
      for (k=0; k<num_prev; k++){
	_hessmatrix(k,num_prev-_blocksize) += dense_vec(k);
      }
      // Compute q_vec<- q_vec - Q_vec * dense_vec
      //
      q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
      //
      q_vec->MvNorm(nm2);
      //
      if (nm2[0] < nm1[0] * _dep_tol) {
	// 
	// Repeat process with newly computed q_vec
	//
	// Compute trans(Q_vec)*q_vec
	//
	q_vec->MvTransMv(one, *Q_vec, dense_vec);
	//
	// Sum results [0:num_prev-1] into column (num_prev-_blocksize)
	// of the Hessenberg matrix
	//
	for (k=0; k<num_prev; k++){
	  _hessmatrix(k,num_prev-_blocksize) += dense_vec(k);
	}
	// Compute q_vec<- q_vec - Q_vec * dense_vec
	//
	q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
	//
	q_vec->MvNorm(nm2);
      }
      //
      // Check for linear dependence
      //
      if (nm2[0] < nm1[0] * _sing_tol) {
	dep = true;
      }
      if (!dep){
	//
	// Normalize the new q_vec
	//
	TYPE rjj = one/nm2[0];
	q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
	//
	// Enter norm of q_vec to the [(j+1)*_blocksize + iter] row
	// in the [(j*_blocksize + iter] column of the Hessenberg matrix
	// 
	_hessmatrix( num_prev, num_prev-_blocksize ) = nm2[0];
      }
      else { 
	//
	if (_debuglevel > 3 && vb) {
	  cout << "Column " << num_prev << " of _basisvecs is dependent" << endl;
	  cout << endl;
	}
	//
	// Create a random vector and orthogonalize it against all 
	// previous cols of _basisvecs
	// We could try adding a random unit vector instead -- not 
	// sure if this would make any difference.
	//
	tptr->MvRandom();
	tptr->MvNorm(nm1);
	//
	// This code  is automatically doing 2 steps of orthogonalization
	// after adding a random vector. We could do one step of
	// orthogonalization with a correction step if needed.
	//
	for (num_orth=0; num_orth<2; num_orth++){
	  tptr->MvTransMv(one, *Q_vec, dense_vec);
	  // Note that we don't change the entries of the
	  // Hessenberg matrix when we orthogonalize a 
	  // random vector
	  tptr->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
	}
	//
	tptr->MvNorm(nm2);
	//
	if (nm2[0] >= nm1[0] * _sing_tol){ 
	  // Copy vector into the current column of _basisvecs
	  q_vec->MvAddMv( one, *tptr, zero, *tptr );
	  q_vec->MvNorm(nm2);
	  // Normalize the new q_vec
	  //
	  TYPE rjj = one/nm2[0];
	  q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
	  //
	  // Enter a zero in the [(j+1)*_blocksize + iter] row in the
	  // [(j*_blocksize + iter] column of the Hessenberg matrix
	  //
	  _hessmatrix( num_prev, num_prev-_blocksize ) = zero;
	}
	else {
	  // Can't produce a new orthonormal basis vector
	  // Return a flag so we can exit this pass of block GMRES
	  flg = true;
	  // Clean up memory
	  delete [] index;
	  delete q_vec; q_vec = 0;
	  delete Q_vec; Q_vec = 0;
	  delete tptr; tptr = 0;
	  delete F_vec; F_vec = 0;
	  return flg;
	}
	//
      } // end else 
	//
    } // end for (iter=0;...)
      //
    if (_debuglevel > 2){
      if (vb){
	cout << endl;
	cout << "Checking Orthogonality after BlkOrthSing()"
	     << " Iteration: " << _iter << endl;
      }
      CheckKrylovOrth(_iter, vb);
    }
    //
    //	free heap space
    //
    delete [] index;
    delete q_vec; q_vec=0;
    delete Q_vec; Q_vec=0;
    delete tptr; tptr=0;
    delete F_vec;
    //
    return flg;
    //
  } // end BlkOrthSing()
  
  
  
  template<class TYPE>
  bool BlockGmres<TYPE>::QRFactorAug(MultiVec<TYPE>& VecIn, 
				     Teuchos::SerialDenseMatrix<int,TYPE>& FouierR, 
				     bool blkone, 
				     bool vb) 
  {
    int i,j,k;
    int nb = VecIn.GetNumberVecs(); assert (nb == _blocksize);
    int *index = new int[nb]; assert(index!=NULL);
    const int IntOne = 1;
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    bool addvec = false;
    bool flg = false;
    //
    TYPE norm1[IntOne];
    TYPE norm2[IntOne];
    Teuchos::SerialDenseVector<int,TYPE> rj; 
    MultiVec<TYPE> *qj = 0, *Qj = 0, *tptr = 0;
    tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
    //
    // Zero out the array that will contain the Fourier coefficients.
    //
    for ( j=0; j<nb; j++ ) {
      for ( i=0; i<nb; i++ ) {
	FouierR(i,j) = zero;
      }
      index[j] = j;
    }
    //
    // Start the loop to orthogonalize the nb columns of VecIn.
    //
    for ( j=0; j<nb; j++ ) {
      //
      flg = false;
      //
      // Grab the j-th column of VecIn (the first column is indexed to 
      // be the zero-th one).
      //
      qj = VecIn.CloneView(index+j, IntOne); assert(qj!=NULL);
      //
      // If we are beyond the 1st column, orthogonalize against the previous
      // vectors in the current block
      //
      if ( j ) {
	//
	// Grab the first j columns of VecIn (that are now an orthogonal
	// basis for first j columns of the entering VecIn).
	//
	rj.size(j);
	Qj = VecIn.CloneView(index, j);
	qj->MvNorm(norm1);
	//
	// Do one step of classical Gram-Schmidt orthogonalization
	// with a second correction step if needed
	//
	// Determine the Fouier coefficients for orthogonalizing column
	// j of VecIn against columns 0:j-1 of VecIn. In other words,
	// result = trans(Qj)*qj.
	//
	qj->MvTransMv( one, *Qj, rj );
	//
	// Sum results[0:j-1] into column j of R.
	//
	for ( k=0; k<j; k++ ) {
	  FouierR(k,j) += rj(k);
	}
	//
	// Compute qj <- qj - Qj * rj.
	//
	qj->MvTimesMatAddMv(-one, *Qj, rj, one);
	//
	qj->MvNorm(norm2);
	//
	if (norm2[0] < norm1[0] * _dep_tol){
	  //
	  // Repeat process with newly computed qj
	  //
	  qj->MvTransMv( one, *Qj, rj );
	  //
	  // Sum results[0:j-1] into column j of R.
	  //
	  for ( k=0; k<j; k++ ) {
	    FouierR(k,j) += rj(k);
	  }
	  //
	  // Compute qj <- qj - Qj * rj.
	  //
	  qj->MvTimesMatAddMv(-one, *Qj, rj, one);
	  //
	  qj->MvNorm(norm2);
	}
	//
	// Check for dependencies
	//
	if (!blkone) {
	  // This is not the 1st block. A looser tolerance is used to 
	  // determine dependencies. If a dependency is detected, a flag
	  // is set so we can back out this routine and out of BlkOrth. 
	  // The routine BlkOrthSing is used to construct the new block 
	  // of orthonormal basis vectors one at a time. If a dependency
	  // is detected within this routine, a random vector is added 
	  // and orthogonalized against all previous basis vectors.
	  // 
	  //
	  if (norm2[0] < norm1[0] * _blk_tol) {
	    if (_debuglevel > 3 && vb) {
	      cout << "Column " << j << " of current block is dependent" << endl;
	    }
	    flg = true;  
	    delete qj; delete Qj; delete tptr;
	    delete [] index;
	    return flg;
	  }
	}
	else {
	  // This is the 1st block of basis vectors.
	  // Use a tighter tolerance to determine dependencies, because
	  // if a dependency is detected we will be adding a random
	  // vector and orthogonalizing it against previous vectors
	  // in the 1st block
	  //
	  if (norm2[0] < norm1[0] * _sing_tol) {
	    // The 1st block of vectors are dependent
	    // Add a random vector and orthogonalize it against
	    // previous vectors in block.
	    //
	    addvec = true;
	    Teuchos::SerialDenseVector<int,TYPE> tj(j);
	    //
	    tptr->MvRandom();
	    tptr->MvNorm(norm1);
	    //
	    int num_orth;
	    for (num_orth=0; num_orth<2; num_orth++){
	      tptr->MvTransMv(one, *Qj, tj);
	      tptr->MvTimesMatAddMv(-one, *Qj, tj, one);
	    }
	    tptr->MvNorm(norm2); 
	    //
	    if (norm2[0] >= norm1[0] * _sing_tol){
	      // Copy vector into current column of _basisvecs
	      qj->MvAddMv(one, *tptr, zero, *tptr);
	    }
	    else {
	      flg = true;
	      delete qj; delete Qj; delete tptr;
	      delete [] index;
	      return flg;
	    } 
	  } 
	} // end else
      } // end if (j)
      //
      // If we have not exited, compute the norm of column j of
      // VecIn (qj), then normalize qj to make it into a unit vector
      //
      TYPE normq[IntOne];
      qj->MvNorm(normq);
      //
      TYPE rjj = one / normq[0];
      qj->MvAddMv ( rjj, *qj, zero, *qj );
      //
      if (addvec){
	// We've added a random vector, so
	// enter a zero in j'th diagonal element of R
	FouierR(j,j) = zero;
      }
      else {
	FouierR(j,j) = normq[0];
      }
      delete qj; delete Qj;
      //
    } // end for (j=0; j<nb; j++)
      //
    delete [] index;
    delete tptr; tptr = 0;
    return flg;
    //
  } // end QRFactorAug()
  
  
  template<class TYPE>
  void BlockGmres<TYPE>::CheckKrylovOrth( const int j, bool vb ) 
  {
    int i,k,m=(j+1)*_blocksize;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::SerialDenseMatrix<int,TYPE> VTV; VTV.shape(m,m);
    int *index = new int[m];
    
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = m+i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _blocksize);
    assert(F_vec!=NULL);
    
    TYPE *ptr_norms = new double[m];
    TYPE sum = zero;
    
    F_vec->MvNorm(ptr_norms);
    for ( i=0; i<_blocksize; i++ ) {
      sum += ptr_norms[i];
    }
    
    for ( i=0; i<m; i++ ) {
      index[i] = i;
    }
    MultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
    assert(Vj!=NULL);
    Vj->MvTransMv(one, *Vj, VTV);
    TYPE column_sum;
    //
    if (vb){
      cout << " " <<  endl;
      cout << "********Block Arnoldi iteration******** " << j <<  endl;
      cout << " " <<  endl;
    }
    //
    for (k=0; k<m; k++) {
      column_sum = zero;
      for (i=0; i<m; i++) {
	if (i==k) {
	  VTV(i,i) -= one;
	}
	column_sum += VTV(i,k);
      }
      if (vb){
	cout <<  " V^T*V-I " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) <<  endl;
      }
    }
    if (vb) {cout << " " <<  endl;}
    
    Teuchos::SerialDenseMatrix<int,TYPE> E; E.shape(m,_blocksize);
    
    F_vec->MvTransMv(one, *Vj, E);
    
    for (k=0;k<_blocksize;k++) {
      column_sum = zero;
      for (i=0; i<m; i++) {
	column_sum += E(i,k);
      }
      if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
      if (vb){
	cout << " Orthogonality with F " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) <<  endl;
      }
    }
    if (vb) {cout << " " <<  endl;}
    //
    delete F_vec;
    delete Vj;
    delete [] index;
    delete [] ptr_norms;
    //
  } // end CheckKrylovOrth
    
} // end namespace Belos

#endif
// End of file BelosBlockGmres.hpp


