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
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosMultiVecTraits.hpp"

/*!	
  \class Belos::BlockGmres
  
  \brief This class implements the Restarted Block GMRES algorithm
  for solving real nonsymmetric linear systems of equations AX = B,
  where B is a matrix containing one or more right-hand sides, and 
  X is the matrix of corresponding solutions.
  
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
  template <class TYPE, class OP, class MV>
  class BlockGmres : public IterativeSolver<TYPE,OP,MV> { 
  public:
    //@{ \name Constructor/Destructor.
    //! %Belos::BlockGmres constructor.
    BlockGmres(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp, 
	       const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest,
               const RefCountPtr<OutputManager<TYPE> > &om,
	       const int length=25);
    
    //! %Belos::BlockGmres destructor.
    virtual ~BlockGmres();
    //@}
    
    //@{ \name Accessor methods
    
    //! Get the iteration count for the current block of linear systems.
    int GetNumIters() const { return( _totaliter ); }
    
    //! Get the restart count of the iteration method for the current block of linear systems.
    int GetNumRestarts() const { return( _restartiter ); }
    
    //! Get the solvers native residuals for the current block of linear systems.
    /*! For GMRES this is not the same as the actual residual of the linear system and the
	residual is not in MultiVec form, so the normvec will be populated with the residual norm.
	
	\param  normvec  [out] Array (length <tt>blocksize</tt>) of computed norm values
	                 for current residuals.
     */
    RefCountPtr<const MV> GetNativeResiduals( TYPE *normvec ) const;

    //! Get the true residuals for the current block of linear systems.
    /*! For GMRES this will force the solver to compute a current residual for its linear 
      systems, the current solution is not stored. <b> This is an expensive computation, 
      so a convergence test using these residuals should be secondary to using the native 
      residuals. </b>
    */
    RefCountPtr<MV> GetCurrentSoln();

    //! Get a constant reference to the current linear problem.  
    /*! This may include a current solution, if the solver has recently restarted or completed.
     */
    LinearProblemManager<TYPE,OP,MV>& GetLinearProblem() const { return( *_lp ); }

    //@} 

    //@{ \name Solver application method.
    
    /*! \brief This method uses the iterative method to compute approximate
      solutions to the original problem.  This method can return unconverged if the
      maximum number of iterations is reached, or numerical breakdown is observed.
    */
    void Solve();
    //@}
    
  private:

    //! Method for setting the basis dependency tolerances for extending the Krylov basis.
    void SetGmresBlkTols();

    //! Method for performing the block Krylov decomposition.
    bool BlockReduction(bool&);

    //! Method for orthogonalization of one block.
    bool QRFactorAug(MV&, Teuchos::SerialDenseMatrix<int,TYPE>&,
		     bool);

    //! Method for block orthogonalization when a dependency has not been detected in the Krylov basis.
    bool BlkOrth(MV&);

    //! Method for block orthogonalization when a dependency has been detected in the Krylov basis.
    bool BlkOrthSing(MV&);

    //! Method for checking the orthogonality of the Krylov basis.
    void CheckKrylovOrth(const int);

    //! Reference to the linear problem being solver for with the solver. [passed in by user]
    RefCountPtr<LinearProblemManager<TYPE,OP,MV> > _lp;

    //! Reference to the status test, which provides the stopping criteria for the solver. [passed in by user]
    RefCountPtr<StatusTest<TYPE,OP,MV> > _stest; 

    //! Reference to the output manager for this linear solver. [passed in by user]
    RefCountPtr<OutputManager<TYPE> > _om;

    //! Pointers to the Krylov basis constructed by the solver.
    RefCountPtr<MV> _basisvecs;

    //! Pointers to the current right-hand side and solution multivecs being solved for.
    RefCountPtr<MV> _cur_block_rhs, _cur_block_sol;

    //! Dense matrices for holding the upper Hessenberg matrix (H) of the Arnoldi factorization 
    Teuchos::SerialDenseMatrix<int,TYPE> _hessmatrix;

    //! Dense vector for holding the right-hand side of the least squares problem.
    Teuchos::SerialDenseMatrix<int,TYPE> _z;

    //! The output stream for sending solver information.
    ostream *_os;

    const int _length;
    int _blocksize;
    int _restartiter, _totaliter, _iter;
    TYPE _dep_tol, _blk_tol, _sing_tol;

    typedef MultiVecTraits<TYPE,MV>  MVT;
  };
  //
  // Implementation
  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class TYPE, class OP, class MV>
  BlockGmres<TYPE,OP,MV>::BlockGmres(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp, 
			       const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest,
			       const RefCountPtr<OutputManager<TYPE> >&om,
			       const int length) : 
    _lp(lp),
    _stest(stest),
    _om(om),
    _os(&om->GetOStream()),
    _length(length), 
    _blocksize(0), 
    _restartiter(0), 
    _totaliter(0),
    _iter(0)
  {
    //
    // Set up the block orthogonality tolerances
    //
    SetGmresBlkTols();	
  }
    
  template <class TYPE, class OP, class MV>
  BlockGmres<TYPE,OP,MV>::~BlockGmres() 
  {}
  
  template <class TYPE, class OP, class MV>
  void BlockGmres<TYPE,OP,MV>::SetGmresBlkTols() 
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
  
  template <class TYPE, class OP, class MV>
  RefCountPtr<const MV> BlockGmres<TYPE,OP,MV>::GetNativeResiduals( TYPE *normvec ) const 
  {
    //
    // If this is the first iteration for a new right-hand side return the
    // residual for the current block rhs and solution.
    //
    if (_totaliter == 0) {
      RefCountPtr<MV> temp_res = MVT::Clone(*_cur_block_rhs,_blocksize);
      _lp->ComputeResVec( &*temp_res, &*_cur_block_sol, &*_cur_block_rhs );
      MVT::MvNorm( *temp_res, normvec, TwoNorm );
    } else {
      if (normvec) {
        Teuchos::BLAS<int,TYPE> blas;
	for (int j=0; j<_blocksize; j++)
	  normvec[j] = blas.NRM2( _blocksize, &_z(_iter*_blocksize, j ), 1);
      }
    }
    return null;
  }
  
  template <class TYPE, class OP, class MV>
  RefCountPtr<MV> BlockGmres<TYPE,OP,MV>::GetCurrentSoln()
  {    
    //
    // If this is the first iteration of the Arnoldi factorization, return the current solution.
    // It has either been updated recently, if there was a restart, or we haven't computed anything yet.
    //
    RefCountPtr<MV> cur_sol_copy = MVT::CloneCopy(*_cur_block_sol);
    if (_iter==0) { 
        return cur_sol_copy;
    } else {
      const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
      int i, m = _iter*_blocksize;
      Teuchos::BLAS<int,TYPE> blas;
      int *index = new int[m]; assert(index!=NULL);
      for ( i=0; i<m; i++ ) {   
        index[i] = i;
      }
      RefCountPtr<const MV> Vjp1 = MVT::CloneView( MVT::c(*_basisvecs), index, m );
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
    
      MVT::MvTimesMatAddMv( one, *Vjp1, y, one, *cur_sol_copy );
    
      if (index) delete [] index;
    }
    return cur_sol_copy;
  }
    
  template <class TYPE, class OP, class MV>
  void BlockGmres<TYPE,OP,MV>::Solve () 
  {
    int i,j, maxidx;
    TYPE *beta=0;
    int *index=0;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    TYPE sigma, mu, vscale, maxelem;
    RefCountPtr<MV> U_vec;
    bool dep_flg = false, exit_flg = false;
    Teuchos::LAPACK<int, TYPE> lapack;
    Teuchos::BLAS<int, TYPE> blas;
    //
    // Obtain the first block linear system form the linear problem manager.
    //
    _cur_block_sol = _lp->GetCurrLHSVec();
    _cur_block_rhs = _lp->GetCurrRHSVec();
    _blocksize = _lp->GetBlockSize();
    //
    //  Start executable statements. 
    //
    while ( _cur_block_rhs.get() && _cur_block_sol.get() ) {
      //
      if (_om->doOutput( 0 )) {
        *_os << endl;
        *_os << "===================================================" << endl;
        *_os << "Solving linear system(s):  " << _lp->GetRHSIndex() << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve() << endl;
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
      //
      // Create the rectangular Hessenberg matrix and right-hand side of least squares problem.
      //
      _hessmatrix.shape((_length+1)*_blocksize, _length*_blocksize);
      _z.shape((_length+1)*_blocksize, _blocksize); 
      //
      beta = new TYPE[(_length+1)*_blocksize]; assert(beta!=NULL);
      index = new int[ (_length+1)*_blocksize ]; assert(index!=NULL);
      for (i=0; i < (_length+1)*_blocksize; i++) { index[i] = i; }
      //
      for (_restartiter=0; _stest->CheckStatus(this) == Unconverged && !exit_flg; _restartiter++) {
	//
	// Associate the initial block of _basisvecs with U_vec.
	//
	U_vec = MVT::CloneView(*_basisvecs, index, _blocksize);
	//
	// Compute current residual and place into 1st block
	//
	_lp->ComputeResVec( &*U_vec, &*_cur_block_sol, &*_cur_block_rhs );
	//
	dep_flg = false; exit_flg = false;
	//
	// Re-initialize RHS of the least squares system and create a view.
	//
	_z.putScalar();
	Teuchos::SerialDenseMatrix<int,TYPE> G10(Teuchos::View, _z, _blocksize, _blocksize);
	exit_flg = QRFactorAug( *U_vec, G10, true );
	//
	if (exit_flg){
	  if (_om->doOutput( 0 )){
	    *_os << "Exiting Block GMRES" << endl;
	    *_os << "  Restart iteration# " << _restartiter
		 << "  Iteration# " << _iter << endl;
	    *_os << "  Reason: Failed to compute initial block of orthonormal basis vectors"
		 << endl << endl;
	  }
	  if (U_vec.get()) {U_vec = null;}
	}
	//
	for (_iter=0; _iter<_length && _stest->CheckStatus(this) == Unconverged && !exit_flg; _iter++, ++_totaliter) {
	  //
	  // Compute a length _length block Arnoldi Reduction (one step at a time),
	  // the exit_flg indicates if we cannot extend the Arnoldi Reduction.
          // If exit_flg is true, then we need to leave this loop and compute the latest solution.
	  //
	  //dep_flg = true;
	  exit_flg = BlockReduction(dep_flg);
	  if (exit_flg){ 
	    break;
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
	} // end for (_iter=0;...
	//
	// Update the solutions by solving the triangular system to get the Krylov weights.
	//
        if (_iter) {
	  // Make a copy of _z since it may be used in the convergence test to compute native residuals.
	  Teuchos::SerialDenseMatrix<int,TYPE> _z_copy( Teuchos::Copy,_z, _iter*_blocksize, _blocksize );	
	  blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, _iter*_blocksize, _blocksize, one,
		   _hessmatrix.values(), _hessmatrix.stride(), _z_copy.values(), _z_copy.stride() ); 
	  // Create view into current basis vectors.
	  RefCountPtr<const MV> Vjp1 = MVT::CloneView(*_basisvecs, index, _iter*_blocksize);
          RefCountPtr<MV> solnUpdate = MVT::Clone(*_cur_block_sol,_blocksize );
	  MVT::MvTimesMatAddMv( one, *Vjp1, _z_copy, zero, *solnUpdate );
	  //
	  // Update the solution held by the linear problem.
	  //	  
	  _lp->SolutionUpdated( &*solnUpdate );
        }
	if (_om->doOutput( 0 )) {
	  if (exit_flg) {
	    *_os << " Exiting Block GMRES --- " << endl;
	    *_os << "  Reason: Failed to compute new block of orthonormal basis vectors" << endl;
	    *_os << "  ***Solution from previous step will be returned***"<< endl<< endl;
	  }
	}
	if (U_vec.get()) {U_vec = null;}
	//
	// Break out of this loop before the _restartiter is incremented if we are finished.
	//
        if ( _stest->GetStatus() != Unconverged || exit_flg ) { break; }
        //
      } // end for (_restartiter=0;...
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
      // **************Free heap space**************
      //
      _basisvecs = null;
      if (index) {delete [] index; index=0;}
      if (beta) {delete [] beta; beta=0; }
      //
      // Print out solver status
      // 
      if (_om->doOutput( 0 )) {
	*_os << endl;
	_stest->Print(*_os); 
	*_os << endl;
	*_os << "===================================================" << endl;
      }
    } // end while( _cur_block_sol && _cur_block_rhs )
    //
  } // end Solve()
  
    
  template<class TYPE, class OP, class MV>
  bool BlockGmres<TYPE,OP,MV>::BlockReduction ( bool& dep_flg ) 
  {
    //
    int i;	
    int *index = new int[_blocksize]; assert(index!=NULL);
    RefCountPtr<MV> AU_vec = MVT::Clone(*_basisvecs,_blocksize);
    //
    // Associate the j-th block of _basisvecs with U_vec.
    //
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = _iter*_blocksize+i;
    }
    RefCountPtr<MV> U_vec = MVT::CloneView(*_basisvecs, index, _blocksize);
    //
    _lp->Apply( *U_vec, *AU_vec ); 
    //
    bool dep = false;
    if (!dep_flg){
      dep = BlkOrth(*AU_vec);
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
      flg = BlkOrthSing(*AU_vec);
    }
    //
    delete [] index;
    //
    return flg;
    //
  } // end BlockReduction()
  
  
  template<class TYPE, class OP, class MV>
  bool BlockGmres<TYPE,OP,MV>::BlkOrth( MV& VecIn ) 
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
    RefCountPtr<MV> F_vec = MVT::CloneView(*_basisvecs, index+(_iter+1)*_blocksize, _blocksize);
    //
    // Copy preconditioned AU_vec into (j+1)st block of _basisvecs
    //
    MVT::MvAddMv( one, VecIn, zero, VecIn, *F_vec );
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
    RefCountPtr<MV> V_prev = MVT::CloneView(*_basisvecs, index, num_prev);
    //
    // Create a matrix to store the product trans(V_prev)*F_vec
    //
    Teuchos::SerialDenseMatrix<int,TYPE> dense_mat( num_prev, _blocksize );
    //
    MVT::MvNorm(*F_vec,norm1);
    //
    // Perform two steps of block classical Gram-Schmidt so that
    // F_vec is orthogonal to the columns of V_prev.
    //
    for ( int num_orth=0; num_orth<max_num_orth; num_orth++ ) {
      //
      // Compute trans(V_prev)*F_vec and store in the j'th diagonal
      // block of the Hessenberg matrix
      //
      MVT::MvTransMv(*F_vec,one, *V_prev, dense_mat);
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
      MVT::MvTimesMatAddMv( -one, *V_prev, dense_mat, one, *F_vec );
    } // end for num_orth=0;...)
      //
    MVT::MvNorm(*F_vec,norm2);
    //
    // Check to make sure the new block of Arnoldi vectors are 
    // not dependent on previous Arnoldi vectors
    //
    bool flg = false; // This will get set true if dependencies are detected
    //
    for (i=0; i<_blocksize; i++){
      if (norm2[i] < norm1[i] * _blk_tol) {
	flg = true;
	if (_om->doOutput( 3 )){
	  *_os << "Col " << num_prev+i << " is dependent on previous "
	       << "Arnoldi vectors in V_prev" << endl;
	  *_os << endl;
	}
      }
    } // end for (i=0;...)
      //
    if (_om->doOutput( 2 )) {
      *_os << "Checking Orthogonality after BlkOrth()"
	     << " Iteration: " << _iter << endl;
      CheckKrylovOrth(_iter);
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
      flg = QRFactorAug( *F_vec, sub_block_hess, false );
    }
    //
    delete [] index;
    delete [] norm1;
    delete [] norm2;
    //
    return flg;
    //
  }  // end BlkOrth()
  
  
  template<class TYPE, class OP, class MV>
  bool BlockGmres<TYPE,OP,MV>::BlkOrthSing( MV& VecIn ) 
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
    RefCountPtr<MV> F_vec = MVT::CloneView(*_basisvecs, index+(_iter+1)*_blocksize, _blocksize);
    //
    // Copy preconditioned AU_vec into (j+1)st block of _basisvecs
    //
    MVT::MvAddMv( one, VecIn, zero, VecIn, *F_vec );
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
    RefCountPtr<const MV> Q_vec;
    RefCountPtr<MV> q_vec, tptr;
    tptr = MVT::Clone(*F_vec,IntOne);
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
      q_vec = MVT::CloneView(*_basisvecs, index2, IntOne);
      //
      // Grab all previous columns of _basisvecs
      //
      Q_vec = MVT::CloneView(MVT::c(*_basisvecs), index, num_prev);
      //
      // Do one step of classical Gram-Schmidt orthogonalization
      // with a 2nd correction step if needed.
      //
      bool dep = false;
      MVT::MvNorm(*q_vec,nm1);
      //
      // Compute trans(Q_vec)*q_vec
      //
      MVT::MvTransMv(*q_vec, one, *Q_vec, dense_vec);
      //
      // Sum results [0:num_prev-1] into column (num_prev-_blocksize)
      // of the Hessenberg matrix
      //
      for (k=0; k<num_prev; k++){
	_hessmatrix(k,num_prev-_blocksize) += dense_vec(k);
      }
      // Compute q_vec<- q_vec - Q_vec * dense_vec
      //
      MVT::MvTimesMatAddMv( -one, *Q_vec, dense_vec, one, *q_vec );
      //
      MVT::MvNorm(*q_vec,nm2);
      //
      if (nm2[0] < nm1[0] * _dep_tol) {
	// 
	// Repeat process with newly computed q_vec
	//
	// Compute trans(Q_vec)*q_vec
	//
	MVT::MvTransMv(*q_vec, one, *Q_vec, dense_vec);
	//
	// Sum results [0:num_prev-1] into column (num_prev-_blocksize)
	// of the Hessenberg matrix
	//
	for (k=0; k<num_prev; k++){
	  _hessmatrix(k,num_prev-_blocksize) += dense_vec(k);
	}
	// Compute q_vec<- q_vec - Q_vec * dense_vec
	//
	MVT::MvTimesMatAddMv(-one, *Q_vec, dense_vec, one, *q_vec);
	//
	MVT::MvNorm(*q_vec,nm2);
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
	MVT::MvAddMv( rjj, *q_vec, zero, *q_vec, *q_vec ); // RAB: Warning, aliasing of arguments!
	//
	// Enter norm of q_vec to the [(j+1)*_blocksize + iter] row
	// in the [(j*_blocksize + iter] column of the Hessenberg matrix
	// 
	_hessmatrix( num_prev, num_prev-_blocksize ) = nm2[0];
      }
      else { 
	//
	if (_om->doOutput( 3 )) {
	  *_os << "Column " << num_prev << " of _basisvecs is dependent" << endl;
	  *_os << endl;
	}
	//
	// Create a random vector and orthogonalize it against all 
	// previous cols of _basisvecs
	// We could try adding a random unit vector instead -- not 
	// sure if this would make any difference.
	//
	MVT::MvRandom(*tptr);
	MVT::MvNorm(*tptr,nm1);
	//
	// This code  is automatically doing 2 steps of orthogonalization
	// after adding a random vector. We could do one step of
	// orthogonalization with a correction step if needed.
	//
	for (num_orth=0; num_orth<2; num_orth++){
	  MVT::MvTransMv(*tptr, one, *Q_vec, dense_vec);
	  // Note that we don't change the entries of the
	  // Hessenberg matrix when we orthogonalize a 
	  // random vector
	  MVT::MvTimesMatAddMv(-one, *Q_vec, dense_vec, one, *tptr);
	}
	//
	MVT::MvNorm(*tptr,nm2);
	//
	if (nm2[0] >= nm1[0] * _sing_tol){ 
	  // Copy vector into the current column of _basisvecs
	  MVT::MvAddMv( one, *tptr, zero, *tptr, *q_vec );
	  MVT::MvNorm(*q_vec,nm2);
	  // Normalize the new q_vec
	  //
	  TYPE rjj = one/nm2[0];
	  MVT::MvAddMv( rjj, *q_vec, zero, *q_vec, *q_vec ); // RAB: Warning, aliasing of arguments!
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
	  return flg;
	}
	//
      } // end else 
	//
    } // end for (iter=0;...)
      //
    if (_om->doOutput( 2 )){
      	*_os << endl;
	*_os << "Checking Orthogonality after BlkOrthSing()"
	     << " Iteration: " << _iter << endl;
      CheckKrylovOrth(_iter);
    }
    //
    //	free heap space
    //
    delete [] index;
    //
    return flg;
    //
  } // end BlkOrthSing()
  

  template<class TYPE, class OP, class MV>
  bool BlockGmres<TYPE,OP,MV>::QRFactorAug(MV& VecIn, 
				     Teuchos::SerialDenseMatrix<int,TYPE>& FouierR, 
				     bool blkone) 
  {
    int i,j,k;
    int nb = MVT::GetNumberVecs(VecIn); assert (nb == _blocksize);
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
    RefCountPtr<const MV> Qj;
    RefCountPtr<MV> qj, tptr;
    tptr = MVT::Clone(*_basisvecs,IntOne);
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
      qj = MVT::CloneView( VecIn, index+j, IntOne );
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
	Qj = MVT::CloneView( MVT::c(VecIn), index, j );
	MVT::MvNorm(*qj,norm1);
	//
	// Do one step of classical Gram-Schmidt orthogonalization
	// with a second correction step if needed
	//
	// Determine the Fouier coefficients for orthogonalizing column
	// j of VecIn against columns 0:j-1 of VecIn. In other words,
	// result = trans(Qj)*qj.
	//
	MVT::MvTransMv( *qj, one, *Qj, rj );
	//
	// Sum results[0:j-1] into column j of R.
	//
	for ( k=0; k<j; k++ ) {
	  FouierR(k,j) += rj(k);
	}
	//
	// Compute qj <- qj - Qj * rj.
	//
	MVT::MvTimesMatAddMv(-one, *Qj, rj, one, *qj);
	//
	MVT::MvNorm(*qj,norm2);
	//
	if (norm2[0] < norm1[0] * _dep_tol){
	  //
	  // Repeat process with newly computed qj
	  //
	  MVT::MvTransMv( *qj, one, *Qj, rj );
	  //
	  // Sum results[0:j-1] into column j of R.
	  //
	  for ( k=0; k<j; k++ ) {
	    FouierR(k,j) += rj(k);
	  }
	  //
	  // Compute qj <- qj - Qj * rj.
	  //
	  MVT::MvTimesMatAddMv(-one, *Qj, rj, one, *qj);
	  //
	  MVT::MvNorm(*qj,norm2);
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
	    if (_om->doOutput( 3 )) {
	      *_os << "Column " << j << " of current block is dependent" << endl;
	    }
	    flg = true;  
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
	    MVT::MvRandom(*tptr);
	    MVT::MvNorm(*tptr,norm1);
	    //
	    int num_orth;
	    for (num_orth=0; num_orth<2; num_orth++){
	      MVT::MvTransMv(*tptr, one, *Qj, tj);
	      MVT::MvTimesMatAddMv(-one, *Qj, tj, one, *tptr);
	    }
	    MVT::MvNorm(*tptr,norm2); 
	    //
	    if (norm2[0] >= norm1[0] * _sing_tol){
	      // Copy vector into current column of _basisvecs
	      MVT::MvAddMv( one, *tptr, zero, *tptr, *qj );
	    }
	    else {
	      flg = true;
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
      MVT::MvNorm(*qj,normq);
      //
      TYPE rjj = one / normq[0];
      MVT::MvAddMv( rjj, *qj, zero, *qj, *qj ); // RAB: Warning aliasing of arguments!
      //
      if (addvec){
	// We've added a random vector, so
	// enter a zero in j'th diagonal element of R
	FouierR(j,j) = zero;
      }
      else {
	FouierR(j,j) = normq[0];
      }
      //
    } // end for (j=0; j<nb; j++)
      //
    delete [] index;
    return flg;
    //
  } // end QRFactorAug()
  
  
  template<class TYPE, class OP, class MV>
  void BlockGmres<TYPE,OP,MV>::CheckKrylovOrth( const int j )
  {
    int i,k,m=(j+1)*_blocksize;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::SerialDenseMatrix<int,TYPE> VTV; VTV.shape(m,m);
    int *index = new int[m];
    
    for ( i=0; i<_blocksize; i++ ) {
      index[i] = m+i;
    }
    RefCountPtr<MV> F_vec = MVT::CloneView(*_basisvecs, index, _blocksize);
    
    TYPE *ptr_norms = new double[m];
    TYPE sum = zero;
    
    MVT::MvNorm(*F_vec,ptr_norms);
    for ( i=0; i<_blocksize; i++ ) {
      sum += ptr_norms[i];
    }
    
    for ( i=0; i<m; i++ ) {
      index[i] = i;
    }
    RefCountPtr<MV> Vj = MVT::CloneView(*_basisvecs,index, m);
    MVT::MvTransMv(*Vj, one, *Vj, VTV);
    TYPE column_sum;
    //
    *_os << " " <<  endl;
    *_os << "********Block Arnoldi iteration******** " << j <<  endl;
    *_os << " " <<  endl;
    //
    for (k=0; k<m; k++) {
      column_sum = zero;
      for (i=0; i<m; i++) {
	if (i==k) {
	  VTV(i,i) -= one;
	}
	column_sum += VTV(i,k);
      }
      *_os <<  " V^T*V-I " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) <<  endl;
    }
    *_os << " " <<  endl;
    
    Teuchos::SerialDenseMatrix<int,TYPE> E; E.shape(m,_blocksize);
    
    MVT::MvTransMv(*F_vec, one, *Vj, E);
    
    for (k=0;k<_blocksize;k++) {
      column_sum = zero;
      for (i=0; i<m; i++) {
	column_sum += E(i,k);
      }
      if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
      *_os << " Orthogonality with F " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) <<  endl;
    }
    *_os << " " <<  endl;
    //
    delete [] index;
    delete [] ptr_norms;
    //
  } // end CheckKrylovOrth
    
} // end namespace Belos

#endif
// End of file BelosBlockGmres.hpp


