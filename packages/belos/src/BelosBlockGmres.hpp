//
// File BelosBlockGmres.hpp
//
// Beta version of code put in repository on 4/01/03
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
#ifndef BLOCK_GMRES_HPP
#define BLOCK_GMRES_HPP

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "BelosConfigDefs.hpp"


/*!	\class Belos::BlockGmres

	\brief This class implements the Restarted Block GMRES algorithm
	for solving real nonsymmetric linear systems of equations AX = B,
	where B is a matrix containing one or more right-hand sides, and 
	X is the matrix of corresponding solutions.

	\author Teri Barth
*/

namespace Belos {

template <class TYPE>
class BlockGmres { 
public:
	//@{ \name Constructor/Destructor.
	//! %Belos::BlockGmres constructor.
	BlockGmres(Anasazi::Matrix<TYPE> & mat, Anasazi::Precondition<TYPE> &precond,
		Anasazi::MultiVec<TYPE>& rhs, 
		const int numrhs, const TYPE tol=1.0e-6, const int maxits=25, 
		const int block=1, bool vb = false);

	//! %Belos::BlockGmres destructor.
	virtual ~BlockGmres();
	//@}

	//@{ \name Solver application method.

	/*! \brief This method uses the iterative method to compute approximate
	solutions to the original problem.  This method can return unconverged if the
	maximum number of iterations is reached, or numerical breakdown is observed.
	*/
	void Solve(bool);
	//@}

	//@{ \name Solution return methods.

	//! This method puts the current solutions in %soln.
	void GetSolutions(Anasazi::MultiVec<TYPE>& soln);
	
	//! This method computes the true residuals for the current solutions. 
	void TrueResiduals(bool);
	//@}

	//@{ \name Set methods.

	//! This method sets the initial guess to be %iguess.
	void SetInitGuess(Anasazi::MultiVec<TYPE>& iguess);

	//! This method sets the number of allowable restarts.
	void SetRestart(const int);
	//@}

	//@{ \name Output methods.

	/*! \brief This method allows the user to set the solver's level of visual output
		during computations.
	*/
	void SetDebugLevel(const int);

	//! This method requests that the solver print out its current residuals.
	void PrintResids(bool)const;
	//@}
private:
	void SetGmresBlkTols();
	void SetUpBlocks(Anasazi::MultiVec<TYPE>&, Anasazi::MultiVec<TYPE>&, int);
	void ExtractCurSolnBlock(Anasazi::MultiVec<TYPE>&, int);
	bool BlockReduction(bool&, bool);
	bool QRFactorAug(Anasazi::MultiVec<TYPE>&, Anasazi::DenseMatrix<TYPE>&,
		                 bool, bool);
	bool BlkOrth(Anasazi::MultiVec<TYPE>&, bool);
	bool BlkOrthSing(Anasazi::MultiVec<TYPE>&, bool);
	void CheckGmresOrth(const int, bool);
	void CheckBlkArnRed(const int, bool);
	void CheckGmresResids(Anasazi::MultiVec<TYPE> &, Anasazi::MultiVec<TYPE> &,
		Anasazi::DenseMatrix<TYPE> &, bool) const;
	Anasazi::Matrix<TYPE> &_amat; // must be passed in by the user
	Anasazi::Precondition<TYPE> &_precond; // must be passed in by the user
	Anasazi::MultiVec<TYPE> &_rhs, *_basisvecs, *_solutions;
	Anasazi::DenseMatrix<TYPE>* _hessmatrix;
	const int _maxits, _blocksize, _numrhs;
	const TYPE _residual_tolerance;
	TYPE *_residerrors, *_trueresids;
	int _rhs_iter, _restartiter, _iter;
	bool _startblock;
	int _debuglevel, _restart;
	TYPE _dep_tol, _blk_tol, _sing_tol, _blkerror;
};
//
// Implementation
//
// Note: I should define a copy constructor and overload = because of the use of new
//
template <class TYPE>
BlockGmres<TYPE>::BlockGmres(Anasazi::Matrix<TYPE> & mat, Anasazi::Precondition<TYPE> &precond, 
							 Anasazi::MultiVec<TYPE>& rhs, 
							 const int numrhs, const TYPE tol, const int maxits, 
							 const int blksz, bool vb) : 
							_amat(mat), _precond(precond),
							_rhs(rhs), _basisvecs(0), _solutions(0),
							_hessmatrix(0),
							_maxits(maxits), _blocksize(blksz), _numrhs(numrhs), _restart(1),
							_residual_tolerance(tol), _residerrors(0), _trueresids(0),
							_rhs_iter(0), _restartiter(0), _iter(0),
							_startblock(false), _debuglevel(0), _dep_tol(0.75),
							_blk_tol(5.0e-7), _sing_tol(5.0e-14), _blkerror(1.0){
	//	cout << "ctor:BlockGmres " << this << endl;
	assert(_maxits>=0); assert(_blocksize>=0); assert(_numrhs>=0);
	//
	// Make room for the Arnoldi vectors and F.
	//
	_basisvecs = _rhs.Clone((_maxits+1)*_blocksize); assert(_basisvecs!=NULL);
	if (_maxits*_blocksize && _basisvecs) {
		//
		// Create the rectangular Hessenberg matrix
		//
		_hessmatrix = new Anasazi::DenseMatrix<TYPE>((_maxits+1)*_blocksize, _maxits*_blocksize); 
		assert(_hessmatrix!=NULL);
		//_residerrors = new TYPE[ _blocksize > _numrhs ? _blocksize : _numrhs ]; assert(_residerrors!=NULL);
		_residerrors = new TYPE[_numrhs + _blocksize]; assert(_residerrors!=NULL);
	}
	else {
		 cout << "BlockGmres:ctor " << _maxits << _blocksize << _basisvecs <<  endl;
		exit(-1);
	}
	//
	// Set up the block orthogonality tolerances
	//
	SetGmresBlkTols();	
}

template <class TYPE>
BlockGmres<TYPE>::~BlockGmres() {
	//	cout << "dtor:BlockGmres " << this << endl;
	if (_basisvecs) delete _basisvecs;
	if (_hessmatrix) delete _hessmatrix;
	if (_solutions) delete _solutions;
	if (_residerrors) delete [] _residerrors;
	if (_trueresids) delete [] _trueresids;
}

template <class TYPE>
void BlockGmres<TYPE>::SetRestart(const int in) {
	if (in > 0 ) {
		//
		// Set the number of restarts
		//
		_restart = in;
	}
}

template <class TYPE>
void BlockGmres<TYPE>::SetDebugLevel(const int in) {
	_debuglevel=in;
}
 
  
template <class TYPE>
void BlockGmres<TYPE>::SetGmresBlkTols() {
	const TYPE two = 2.0;
	TYPE eps;
	char precision = 'P';
	Teuchos::LAPACK<int,TYPE> lapack;
	eps = lapack.LAMCH(precision);
	_blk_tol = 10*sqrt(eps);
	_sing_tol = 10 * eps;
	_dep_tol = 1/sqrt(two);
}

template<class TYPE>
void BlockGmres<TYPE>::SetUpBlocks (Anasazi::MultiVec<TYPE>& sol_block,  
				     Anasazi::MultiVec<TYPE>& rhs_block,  
				     int num_to_solve) {
	//
	int i;
	int *index = new int[_blocksize + _numrhs]; assert(index!=NULL);
	const TYPE one=1.0;
	const TYPE zero=0.0;
	//
	// Logic to handle the number of righthand sides solved
	// for at this iteration.
	//
	if (num_to_solve >= _blocksize) {
		//
		// Easy case: The number of right-hand sides left to solve for is >= the 
		// size of a block. Solve for the next _blocksize of these right-hand sides
		// at this iteration.
	        //
		// Put the next _blocksize of the right-hand sides the rhs_block 
		//
		for ( i=0;i<_blocksize; i++ ) {
			index[i] = _rhs_iter*_blocksize + i;
		}
		Anasazi::MultiVec<TYPE> *tptr = _rhs.CloneView(index,_blocksize); assert(tptr!=NULL);
		rhs_block.MvAddMv(one, *tptr, zero, *tptr);
		//
		// Put the next _blocksize of the initial guesses
		// into the sol_block
		//
		Anasazi::MultiVec<TYPE> *tptr2 = _solutions->CloneView(index,_blocksize); assert(tptr2!=NULL);
		sol_block.MvAddMv(one, *tptr2, zero, *tptr2);
		delete tptr; tptr = 0;
		delete tptr2; tptr2 = 0;
	}
	else {
		// More involved. This is the case where the number of right-hand
		// sides left to solve for at this iteration is less than the block size.
		//
		// Fill up the right-hand side block with random vectors, then place
		// the remaining (unsolved) right hand sides into the initial portion
		// of the right-hand side block.
		//
		int *index2 = new int[num_to_solve]; assert(index2!=NULL);
		rhs_block.MvRandom();
		//
		for ( i=0; i<num_to_solve; i++ ) {
			index[i] = _rhs_iter*_blocksize + i;
		}
		Anasazi::MultiVec<TYPE> *tptr = _rhs.CloneView(index,num_to_solve); assert(tptr!=NULL);
		for (i=0; i<num_to_solve; i++) {
			index2[i] = i;
		}
		rhs_block.SetBlock( *tptr, index2, num_to_solve);
		//
		// Fill up the sol_block with zero vectors, then
		// place the remaining (unsolved) initial guesses into the initial portion
		// of the sol_block.
		//
		sol_block.MvInit( 0.0 );
		//
		Anasazi::MultiVec<TYPE> *tptr2 = _solutions->CloneView(index,num_to_solve); assert(tptr2!=NULL);
		sol_block.SetBlock( *tptr2, index2, num_to_solve);
		//
		// Delete the temporary views
		//
		delete tptr; tptr = 0;
		delete tptr2; tptr2 = 0;
		delete [] index2; index2=0;
	}
	delete [] index; index=0;
	//
}
//
template <class TYPE>
void BlockGmres<TYPE>::ExtractCurSolnBlock(Anasazi::MultiVec<TYPE>& sol_block,
								  int num_to_solve) {
	//
	int i;
	const TYPE one = 1.0;
	const TYPE zero = 0.0;
	int * index = new int[_blocksize + _numrhs]; assert(index!=NULL);
	Anasazi::MultiVec<TYPE> *tptr=0, *tptr2=0;
	//
	if (num_to_solve >= _blocksize) {
		 for (i=0; i<_blocksize; i++) {
			 index[i] = _rhs_iter * _blocksize + i;
		 }
		 tptr = _solutions->CloneView(index,_blocksize); assert(tptr!=NULL);
		 tptr->MvAddMv(one, sol_block, zero, sol_block);
	 }
	 else {
		 for (i=0; i<num_to_solve; i++) {
			 index[i] = _rhs_iter * _blocksize + i;
		 }
		 tptr = _solutions->CloneView(index,num_to_solve); assert(tptr!=NULL);
		 //
		 for (i=0; i<num_to_solve; i++) {
			 index[i] = i;
		 }
		 tptr2 = sol_block.CloneView(index,num_to_solve); assert(tptr2!=NULL);
		 //


		 tptr->MvAddMv(one, *tptr2, zero, *tptr2);
	 }
	 //
	 if (tptr) {
		 delete tptr; tptr=0;
	 }
	 if (tptr2) {
		 delete tptr2; tptr2=0;
	 }
	 delete [] index; index=0;
	//
}
//


template <class TYPE>
void BlockGmres<TYPE>::SetInitGuess(Anasazi::MultiVec<TYPE>& iguess) {
//  This will set the initial guess to the input vector.  If it has less
//  columns than the number of right hand sides, then the rest of _solutions
//  will be filled up with random vectors.
        if (!_startblock) {
                int i, numvecs = iguess.GetNumberVecs();
                int* index = new int[ _numrhs ];
                _solutions = _rhs.Clone(_numrhs); assert(_solutions!=NULL);
                if (numvecs < _numrhs) {
                        for (i=0; i<numvecs; i++) {
                                index[i] = i;
                        }
                        _solutions->SetBlock( iguess, index, numvecs );
                        for (i=numvecs; i<_numrhs; i++) {
                                index[i-numvecs] = i;
                        }       
                        Anasazi::MultiVec<TYPE>* U_vec = _solutions->CloneView( index, _numrhs-numvecs );
                        assert(U_vec!=NULL);
                        U_vec->MvRandom();
                        delete U_vec;
                }
                else {
                        for (i=0; i<_numrhs; i++) {
                                index[i] = i;
                        }
                        _solutions->SetBlock( iguess, index, _numrhs );
                }
                _startblock = true;
		delete [] index; index = 0;
        }
}


template <class TYPE>
void BlockGmres<TYPE>::TrueResiduals (bool vb) {
	Anasazi::MultiVec<TYPE> * AX = _solutions->Clone(_numrhs); assert(AX!=NULL);
	_amat.ApplyMatrix(*_solutions, *AX);
	assert(AX->GetNumberVecs()==_rhs.GetNumberVecs());
	
	TYPE* norms_AX=new TYPE[AX->GetNumberVecs()]; assert(norms_AX!=NULL);
	AX->MvNorm(norms_AX);
	int i;
	TYPE* norms_rhs=new TYPE[_rhs.GetNumberVecs()]; assert(norms_rhs!=NULL);
	_rhs.MvNorm(norms_rhs);
	
	_trueresids = new TYPE[AX->GetNumberVecs()]; assert(_trueresids!=NULL);
	const TYPE one=1.0;
	AX->MvAddMv(one, _rhs, -one, *AX);	
	AX->MvNorm(_trueresids);
    	if(vb){
		cout << "--------------------------------------" << endl;
      		for (i=0; i<AX->GetNumberVecs(); i++){
		  cout <<"Unscaled "<< i <<"-th true residual "<< _trueresids[i] << endl;
	  	}
	  	cout << endl << endl;
	}
	TYPE scale;
	for (i=0; i<AX->GetNumberVecs(); i++ ) {
		scale = AX->GetVecLength();
		scale = sqrt(scale)*(norms_AX[i] + norms_rhs[i])/AX->GetNumberVecs();
		if(vb){
		  if (scale) {
			  cout << "Scaled "<< i << "-th true residual " << _trueresids[i]/scale << endl;
		  }
		  else {
			  if (vb){
			     cout << " scale is zero " << endl;
			     cout << "Scaled "<< i << "-th true residual " << _trueresids[i] << endl;
			  }
		  }
		}
	}
	if (vb) { cout << endl << endl; }
    	//
	if (AX) delete AX;
	if (norms_AX) delete [] norms_AX;
	if (norms_rhs) delete [] norms_rhs;
}


template <class TYPE>   
void BlockGmres<TYPE>::PrintResids(bool vb)const {
	//
	int i;
	//
	if (vb) {
		cout << "--------------------------------------" << endl;
		for (i=0; i<_numrhs; i++){
			cout << "_residerrors[" << i << "] = "
				 << _residerrors[i] << endl;
		}
		cout << endl;
	}
	//
} // end PrintResids



template <class TYPE>
void BlockGmres<TYPE>::GetSolutions(Anasazi::MultiVec<TYPE>& soln) {
	int i, numvecs = soln.GetNumberVecs();
	if (numvecs > _numrhs) {
		numvecs = _numrhs;
	}
	int* index = new int[ numvecs ];
	for (i=0; i<numvecs; i++) {
		index[i] = i;
	}
	soln.SetBlock( *_solutions, index, numvecs );

	delete [] index;
}


template <class TYPE>
void BlockGmres<TYPE>::Solve (bool vb) {
	int i,j, maxidx, numrhs_to_solve;
	const int izero=0;
	const TYPE one=1.0;
	const TYPE zero=0.0;
	TYPE sigma, mu, vscale, maxelem;
	Anasazi::MultiVec<TYPE> *cur_block_sol=0, *cur_block_rhs=0;
	Anasazi::MultiVec<TYPE> *U_vec=0, *F_vec=0;
	int *index = new int[ (_maxits+1)*_blocksize ]; assert(index!=NULL);
	bool dep_flg = false, exit_flg = false, brkflg = false;
	Teuchos::LAPACK<int, TYPE> lapack;
	Teuchos::BLAS<int, TYPE> blas;
	//	
	Anasazi::DenseMatrix<TYPE> rhs((_maxits+1)*_blocksize,_blocksize);
	TYPE *ptr_rhs = rhs.getarray();
	TYPE *beta = new TYPE[(_maxits+1)*_blocksize];
	int ldrhs = rhs.getld();
	//
	// Each pass through the solver solves for _blocksize right-hand sides.
	// max_rhs_iters is the number of passes through the
	// solver required to solve for all right-hand sides	
	//
	int max_rhs_iters = (_numrhs+_blocksize-1) / _blocksize;
	// 
	// If not provided, set initial guesses to AX = B to random vectors
	//
	if (!_startblock) {
	  _solutions = _rhs.Clone(_numrhs); assert(_solutions!=NULL);
	  _solutions->MvRandom();
	  _startblock = true;
	}
	cur_block_sol = _solutions->Clone(_blocksize); assert(cur_block_sol!=NULL);
	cur_block_rhs = _solutions->Clone(_blocksize); assert(cur_block_rhs!=NULL);
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
		if ( _blocksize < _numrhs ) {
			numrhs_to_solve = _numrhs - (_rhs_iter * _blocksize);
		} else {
			numrhs_to_solve = _numrhs;
		}
		//
		// Put the current initial guesses and right-hand sides into current blocks
		//
        	SetUpBlocks(*cur_block_sol, *cur_block_rhs, numrhs_to_solve);
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
		  for ( i=0; i<_blocksize; i++ ) {
		    index[i] = i;
		  }
		  U_vec = _basisvecs->CloneView(index, _blocksize);
		  assert(U_vec!=NULL);
		  //
		  // Copy current solution (initial guess) into 1st block
		  //
		  U_vec->MvAddMv(one, *cur_block_sol, zero, *U_vec);
		  //
		  // Compute the initial residuals then store them in 1st
		  // block of _basisvecs
		  //
		  for ( i=0; i<_blocksize; i++ ) {
		    index[i] = _blocksize+i;
		  }
		  F_vec = _basisvecs->CloneView(index, _blocksize);
		  assert(F_vec!=NULL);
		  //
		  // F_vec <- A*U_vec
		  //
		  _amat.ApplyMatrix( *U_vec,*F_vec );
		  //
		  // U_vec <- cur_block_rhs - A*U_vec
		  //
		  U_vec->MvAddMv(one, *cur_block_rhs, -one, *F_vec);
		  //
		  // Apply the preconditioner
		  //
		  Anasazi::MultiVec<TYPE>* AU_vec = U_vec->CloneCopy(); assert(AU_vec!=NULL);
		  _precond.ApplyPrecondition( *AU_vec, *U_vec );
		  //
		  if (AU_vec) { 
		    delete AU_vec; AU_vec=0;
		  }
		  dep_flg = false; exit_flg = false;
		  //
		  // Re-initialize RHS of the least squares system and create a view.
		  //
		  rhs.init();
		  Anasazi::DenseMatrix<TYPE> G10( rhs, izero, izero, _blocksize, _blocksize);
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
		    if (F_vec) {delete F_vec; F_vec = 0;}
		    break;
		  }
		  //
		  TYPE norm_G10 = G10.getfronorm();
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
		  for (_iter=0; _iter<_maxits && _blkerror > _residual_tolerance; _iter++) {
		    //
		    // Compute a length _maxits block Arnoldi Reduction
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
		    //
		    // Create a view into the rectangular matrix
		    //
		    Anasazi::DenseMatrix<TYPE> Hj(*_hessmatrix, 0, 0, 
						  (_iter+2)*_blocksize, (_iter+1)*_blocksize);
		    int ldhj = Hj.getld();
		    TYPE *ptr_hj = Hj.getarray();
		    //
		    // QR factorization of Least-Squares system with Householder reflectors
		    //
		    for (j=0; j<_blocksize; j++) {
		      //
		      // Apply previous Householder reflectors to new block of Hessenberg matrix
		      //
		      for (i=0; i<_iter*_blocksize+j; i++) {
			sigma = blas.DOT( _blocksize, ptr_hj+i*(ldhj+1)+1, 1, ptr_hj+(_iter*_blocksize+j)*ldhj+i+1, 1);
			sigma += ptr_hj[(_iter*_blocksize+j)*ldhj+i];
			sigma *= beta[i];
			blas.AXPY(_blocksize, -sigma, ptr_hj+i*(ldhj+1)+1, 1, ptr_hj+(_iter*_blocksize+j)*ldhj+i+1, 1);
			ptr_hj[(_iter*_blocksize+j)*ldhj+i] -= sigma;
		      }
		      //
		      // Compute new Householder reflector
		      //
		      maxidx = blas.IAMAX( _blocksize+1, ptr_hj+(_iter*_blocksize+j)*(ldhj+1), 1 );
		      maxelem = ptr_hj[(_iter*_blocksize+j)*(ldhj+1)+maxidx-1];
		      for (i=0; i<_blocksize+1; i++) 
			ptr_hj[(_iter*_blocksize+j)*(ldhj+1)+i] /= maxelem;
		      sigma = blas.DOT( _blocksize, ptr_hj+(_iter*_blocksize+j)*(ldhj+1)+1, 1, ptr_hj+(_iter*_blocksize+j)*(ldhj+1)+1, 1 );
		      if (sigma == zero) {
			beta[_iter*_blocksize + j] = zero;
		      } else {
			mu = sqrt(ptr_hj[(_iter*_blocksize+j)*(ldhj+1)]*ptr_hj[(_iter*_blocksize+j)*(ldhj+1)]+sigma);
			if ( ptr_hj[(_iter*_blocksize+j)*(ldhj+1)] < zero ) {
			  vscale = ptr_hj[(_iter*_blocksize+j)*(ldhj+1)] - mu;
			} else {
			  vscale = -sigma / (ptr_hj[(_iter*_blocksize+j)*(ldhj+1)] + mu);
			}
			beta[_iter*_blocksize+j] = 2.0*vscale*vscale/(sigma + vscale*vscale);
			ptr_hj[(_iter*_blocksize+j)*(ldhj+1)] = maxelem*mu;
			for (i=0; i<_blocksize; i++)
			  ptr_hj[(_iter*_blocksize+j)*(ldhj+1)+1+i] /= vscale;
		      }
		      //
		      // Apply new Householder reflector to rhs
		      //
		      for (i=0; i<_blocksize; i++) {
			sigma = blas.DOT( _blocksize, ptr_hj+(_iter*_blocksize+j)*(ldhj+1)+1, 1, ptr_rhs+(i*ldrhs)+(_iter*_blocksize)+j+1, 1);
			sigma += ptr_rhs[(i*ldrhs)+(_iter*_blocksize)+j];
			sigma *= beta[_iter*_blocksize+j];
			blas.AXPY(_blocksize, -sigma, ptr_hj+(_iter*_blocksize+j)*(ldhj+1)+1, 1, ptr_rhs+(i*ldrhs)+(_iter*_blocksize)+j+1, 1 );
			ptr_rhs[(i*ldrhs)+(_iter*_blocksize)+j] -= sigma;
		      }
		    }
		    //
		    // Compute the residuals and the block error
		    //
		    _blkerror = zero;
		    for (j=0; j<_blocksize; j++ ) {
		      _residerrors[_rhs_iter*_blocksize+j] = blas.NRM2( _blocksize, ptr_rhs+j*ldrhs+(_iter+1)*_blocksize, 1);
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
		    // Compute the true residuals and print them out
		    //
		    if (_debuglevel>1) CheckGmresResids (*cur_block_sol, *cur_block_rhs, rhs, vb);
		    //
		  } // end for (_iter=0;...
		  //
		  // Update the solutions by solving the triangular system to get the Krylov weights.
		  //
		  blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
			     Teuchos::NON_UNIT_DIAG, _iter*_blocksize, _blocksize, one,
			     _hessmatrix->getarray(), _hessmatrix->getld(), ptr_rhs, ldrhs ); 
		  for ( i=0; i<_iter*_blocksize; i++ ) {
		    index[i] = i;
		  }
		  Anasazi::MultiVec<TYPE> * Vjp1 = _basisvecs->CloneView(index, _iter*_blocksize);
		  Anasazi::DenseMatrix<TYPE> rhs_view(rhs, izero, izero, _iter*_blocksize, _blocksize);
		  cur_block_sol->MvTimesMatAddMv( one, *Vjp1, rhs_view, one );
		  delete Vjp1;
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
		  if (F_vec) {delete F_vec; F_vec = 0;}
		  //
		  //
		} // end for (_restartiter=0;...
		//
		//
		// Insert the current block of solutions into _solutions so we can
		// continue if we have more right-hand sides to solve for
		//
        	ExtractCurSolnBlock(*cur_block_sol, numrhs_to_solve);
		//
		//**************Free heap space**************
		//
	} // end for (_rhs_iter =0;...)
	//
	if (index) {delete [] index; index=0;}
	if (cur_block_sol) {delete cur_block_sol; cur_block_sol=0;}
	if (cur_block_rhs) {delete cur_block_rhs; cur_block_rhs=0;}
	//
	//
} // end Solve()


template<class TYPE>
bool BlockGmres<TYPE>::BlockReduction ( bool& dep_flg, bool vb ) {
  //
  int i;	
  int *index = new int[_blocksize]; assert(index!=NULL);
  //
  // Associate the j-th block of _basisvecs with U_vec.
  //
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = _iter*_blocksize+i;
  }
  Anasazi::MultiVec<TYPE>* U_vec = _basisvecs->CloneView(index, _blocksize);
  assert(U_vec!=NULL);
  //
  Anasazi::MultiVec<TYPE>* Temp_vec = _basisvecs->Clone(_blocksize); assert(Temp_vec!=NULL);
  //
  //  Compute Temp_vec = A * U_vec
  //
  _amat.ApplyMatrix( *U_vec, *Temp_vec ); 
  //
  // Apply the preconditioner and store result in AU_vec
  //
  Anasazi::MultiVec<TYPE>* AU_vec = _basisvecs->Clone(_blocksize); assert(AU_vec!=NULL);
  _precond.ApplyPrecondition( *Temp_vec, *AU_vec );
  //
  delete Temp_vec; delete U_vec; delete [] index;
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
  delete AU_vec;
  return flg;
  //
} // end BlockReduction()


template<class TYPE>
bool BlockGmres<TYPE>::BlkOrth( Anasazi::MultiVec<TYPE>& VecIn, bool vb) {
  //
  // Orthogonalization is first done between the new block of 
  // vectors and all previous blocks, then the vectors within the
  // new block are orthogonalized.
  //
  const TYPE one = 1.0;
  const TYPE zero = 0.0;
  const int max_num_orth = 2;
  int i, k, row_offset, col_offset;
  int * index = new int[_blocksize * _maxits]; assert(index!=NULL);
  TYPE * norm1 = new TYPE[_blocksize]; assert(norm1!=NULL);
  TYPE * norm2 = new TYPE[_blocksize]; assert(norm2!=NULL);
  //
  // Associate (j+1)-st block of ArnoldiVecs with F_vec.
  //
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = (_iter+1)*_blocksize+i;
  }
  Anasazi::MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _blocksize);
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
  int ldh = _hessmatrix->getld();
  int n_row = _hessmatrix->getrows();
  int n_col = _hessmatrix->getcols();
  //
  TYPE* ptr_hess = _hessmatrix->getarray();
  for ( k=0; k<_blocksize; k++ ) {
    for ( i=0; i<n_row ; i++ ) {
      ptr_hess[_iter*ldh*_blocksize + k*ldh + i] = zero;
    }
  }
  //
  // Grab all previous Arnoldi vectors
  //
  int num_prev = (_iter+1)*_blocksize;
  for (i=0; i<num_prev; i++){
    index[i] = i;
  }
  Anasazi::MultiVec<TYPE>* V_prev = _basisvecs->CloneView(index,num_prev);
  assert(V_prev!=NULL);
  //
  // Create a matrix to store the product trans(V_prev)*F_vec
  //
  Anasazi::DenseMatrix<TYPE> dense_mat(num_prev, _blocksize );
  TYPE* ptr_dense = dense_mat.getarray();
  int ld_dense = dense_mat.getld();
  int num_orth;
  //
  F_vec->MvNorm(norm1);
  //
  // Perform two steps of block classical Gram-Schmidt so that
  // F_vec is orthogonal to the columns of V_prev.
  //
  for ( num_orth=0; num_orth<max_num_orth; num_orth++ ) {
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
	ptr_hess[_iter*ldh*_blocksize + k*ldh + i] += 
	  ptr_dense[k*ld_dense + i];
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
    CheckGmresOrth(_iter, vb);
    if (_debuglevel>3){
      // CheckBlkArnRed(_iter, vb);
    }
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
    Anasazi::DenseMatrix<TYPE> sub_block_hess(*_hessmatrix, row_offset, col_offset,
					      _blocksize, _blocksize);
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
bool BlockGmres<TYPE>::BlkOrthSing( Anasazi::MultiVec<TYPE>& VecIn, bool vb) {
  //
  // This is a variant of A. Ruhe's block Arnoldi
  // The orthogonalization of the vectors AU_vec is done
  // one at a time. If a dependency is detected, a random
  // vector is added and orthogonalized against all previous
  // Arnoldi vectors.
  // 
  const int IntOne = 1;
  const TYPE one = 1.0;
  const TYPE zero = 0.0;
  int i, k;
  int * index = new int[_blocksize * _maxits +_blocksize]; assert(index!=NULL);
  TYPE nm1[IntOne];
  TYPE nm2[IntOne];
  //
  // Associate (j+1)-st block of ArnoldiVecs with F_vec.
  //
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = (_iter+1)*_blocksize+i;
  }
  Anasazi::MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _blocksize);
  assert(F_vec!=NULL);
  //
  // Copy preconditioned AU_vec into (j+1)st block of _basisvecs
  //
  F_vec->MvAddMv( one, VecIn, zero, VecIn);
  //
  // Zero out the full block column of the Hessenberg matrix 
  //
  int ldh = _hessmatrix->getld();
  int n_row = _hessmatrix->getrows();
  int n_col = _hessmatrix->getcols();
  //
  TYPE* ptr_hess = _hessmatrix->getarray();
  for ( k=0; k<_blocksize; k++ ) {
    for ( i=0; i<n_row ; i++ ) {
      ptr_hess[_iter*ldh*_blocksize + k*ldh + i] = zero;
    }
  }
  //
  Anasazi::MultiVec<TYPE> *q_vec=0, *Q_vec=0, *tptr=0;
  tptr = F_vec->Clone(IntOne); assert(tptr!=NULL);
  //
  // Start a loop to orthogonalize each of the _blocksize
  // columns of F_vec against all previous _basisvecs
  //
  int iter, num_prev;
  bool flg = false;
  //
  for (iter=0; iter<_blocksize; iter++){
    num_prev = (_iter+1)*_blocksize + iter; // number of previous _basisvecs
    //
    // Grab the next column of _basisvecs
    //
    index[0] = num_prev;
    q_vec = _basisvecs->CloneView(index, IntOne); assert(q_vec!=NULL);
    //
    // Grab all previous columns of _basisvecs
    //
    for (i=0; i<num_prev; i++){
      index[i] = i;
    }
    Q_vec = _basisvecs->CloneView(index, num_prev); assert(Q_vec!=NULL);
    //
    // Create matrix to store product trans(Q_vec)*q_vec
    //
    Anasazi::DenseMatrix<TYPE> dense_mat(num_prev, IntOne);
    TYPE* ptr_dense = dense_mat.getarray();
    //
    // Do one step of classical Gram-Schmidt orthogonalization
    // with a 2nd correction step if needed.
    //
    bool dep = false;
    q_vec->MvNorm(nm1);
    //
    // Compute trans(Q_vec)*q_vec
    //
    q_vec->MvTransMv(one, *Q_vec, dense_mat);
    //
    // Sum results [0:num_prev-1] into column (num_prev-_blocksize)
    // of the Hessenberg matrix
    //
    for (k=0; k<num_prev; k++){
      ptr_hess[(_iter*_blocksize + iter)*ldh +k] += ptr_dense[k];
    }
    // Compute q_vec<- q_vec - Q_vec * dense_mat
    //
    q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
    //
    q_vec->MvNorm(nm2);
    //
    if (nm2[0] < nm1[0] * _dep_tol) {
      // 
      // Repeat process with newly computed q_vec
      //
      // Compute trans(Q_vec)*q_vec
      //
      q_vec->MvTransMv(one, *Q_vec, dense_mat);
      //
      // Sum results [0:num_prev-1] into column (num_prev-_blocksize)
      // of the Hessenberg matrix
      //
      for (k=0; k<num_prev; k++){
	ptr_hess[(_iter*_blocksize + iter)*ldh +k] += ptr_dense[k];
      }
      // Compute q_vec<- q_vec - Q_vec * dense_mat
      //
      q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
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
      ptr_hess[(_iter*_blocksize+iter)*ldh + (_iter+1)*_blocksize+iter] = nm2[0];
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
      int num_orth;
      //
      // This code  is automatically doing 2 steps of orthogonalization
      // after adding a random vector. We could do one step of
      // orthogonalization with a correction step if needed.
      //
      for (num_orth=0; num_orth<2; num_orth++){
	tptr->MvTransMv(one, *Q_vec, dense_mat);
	// Note that we don't change the entries of the
	// Hessenberg matrix when we orthogonalize a 
	// random vector
	tptr->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
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
	ptr_hess[(_iter*_blocksize+iter)*ldh + (_iter+1)*_blocksize+iter] = zero;
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
    CheckGmresOrth(_iter, vb);
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
bool BlockGmres<TYPE>::QRFactorAug(Anasazi::MultiVec<TYPE>& VecIn, 
		Anasazi::DenseMatrix<TYPE>& FouierR, bool blkone, bool vb) {
	int i,j,k;
	int nb = VecIn.GetNumberVecs(); assert (nb == _blocksize);
	int ldR = FouierR.getld();
	int *index = new int[nb]; assert(index!=NULL);
	const int IntOne=1;
	const TYPE zero=0.0;
	const TYPE one=1.0;
	bool addvec = false;
	bool flg = false;
	//
	TYPE * R = FouierR.getarray();
	TYPE norm1[IntOne];
	TYPE norm2[IntOne];
	Anasazi::MultiVec<TYPE> *qj = 0, *Qj = 0, *tptr = 0;
	tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
	//
	// Zero out the array that will contain the Fourier coefficients.
	//
	for ( j=0; j<nb; j++ ) {
	  for ( i=0; i<nb; i++ ) {
	    R[j*ldR+i] = zero;
	  }
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
	  index[0] = j;
	  qj = VecIn.CloneView(index, IntOne); assert(qj!=NULL);
	  //
	  // If we are beyond the 1st column, orthogonalize against the previous
	  // vectors in the current block
	  //
	  if ( j ) {
	    for ( i=0; i<j; i++ ) {
	      index[i] = i;
	    }
	    //
	    // Grab the first j columns of VecIn (that are now an orthogonal
	    // basis for first j columns of the entering VecIn).
	    //
	    Qj = VecIn.CloneView(index, j);
	    Anasazi::DenseMatrix<TYPE> rj(j,1);
	    TYPE * result = rj.getarray();
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
	      R[j*ldR+k] += result[k];
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
		R[j*ldR+k] += result[k];
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
		Anasazi::DenseMatrix<TYPE> tj(j,1);
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
	    *(R+j*ldR+j) = zero;
	  }
	  else {
	    *(R+j*ldR+j) = normq[0];
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
void BlockGmres<TYPE>::CheckGmresOrth( const int j, bool vb ) {
  int i,k,m=(j+1)*_blocksize;
  int *index = new int[m];
  
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = m+i;
  }
  Anasazi::MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _blocksize);
  assert(F_vec!=NULL);
  
  TYPE *ptr_norms = new double[m];
  TYPE sum=0.0;
  
  F_vec->MvNorm(ptr_norms);
  for ( i=0; i<_blocksize; i++ ) {
    sum += ptr_norms[i];
  }
  
  for ( i=0; i<m; i++ ) {
    index[i] = i;
  }
  Anasazi::MultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
  assert(Vj!=NULL);
  const TYPE one=1.0;
  const TYPE zero=0.0;
  Anasazi::DenseMatrix<TYPE> VTV(m,m);
  Vj->MvTransMv(one,*Vj,VTV);
  TYPE* ptr=VTV.getarray();
  TYPE column_sum;
  //
  if (vb){
    cout << " " <<  endl;
    cout << "********Block Arnoldi iteration******** " << j <<  endl;
    cout << " " <<  endl;
  }
  //
  for (k=0; k<m; k++) {
    column_sum=zero;
    for (i=0; i<m; i++) {
      if (i==k) {
	ptr[i] -= one;
      }
      column_sum += ptr[i];
    }
    if (vb){
      cout <<  " V^T*V-I " << "for column " << k << " is " << fabs(column_sum) <<  endl;
    }
    ptr += m;
  }
  if (vb) {cout << " " <<  endl;}
  
  Anasazi::DenseMatrix<TYPE> E(m,_blocksize);
  
  F_vec->MvTransMv(one,*Vj,E);
  TYPE* ptr_Ej=E.getarray();
  
  for (k=0;k<_blocksize;k++) {
    column_sum=zero;
    for (i=0; i<m; i++) {
      column_sum += ptr_Ej[i];
    }
    ptr_Ej += m;
    if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
    if (vb){
      cout << " Orthogonality with F " << "for column " << k << " is " << fabs(column_sum) <<  endl;
    }
  }
  if (vb) {cout << " " <<  endl;}
  //
  delete F_vec;
  delete Vj;
  delete [] index;
  delete [] ptr_norms;
  //
} // end CheckGmresOrth
  

template<class TYPE>
void BlockGmres<TYPE>::CheckBlkArnRed( const int j, bool vb ) {
  int i,m=(j+1)*_blocksize;
  const TYPE one=1.0;
  const TYPE zero=0.0;
  TYPE * ptr_norms = new TYPE[m];
  int *index = new int[m];
  //
  for ( i=0; i<m; i++ ) {
    index[i] = i;
  }
  Anasazi::MultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
  assert(Vj!=NULL);
  //
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = m+i;
  }
  Anasazi::MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _blocksize);
  assert(F_vec!=NULL);
  //
  Anasazi::MultiVec<TYPE>* AVj = _basisvecs->Clone(m); assert(AVj!=NULL);
  _amat.ApplyMatrix(*Vj,*AVj);
  //
  // Apply the preconditioner
  //
  Anasazi::MultiVec<TYPE>* MAVj = _basisvecs->Clone(m); assert(AVj_vec!=NULL);
  _precond.ApplyPrecondition( *AVj, *MAVj );

  Anasazi::DenseMatrix<TYPE> Hj(*_hessmatrix, 0, 0, m, m);
  MAVj->MvTimesMatAddMv(-one, *Vj, Hj, one);
  for ( i=0; i<_blocksize; i++ ) {
    index[i] = j*_blocksize+i;
  }
  Anasazi::MultiVec<TYPE>* Fj = MAVj->CloneView(index, _blocksize);
  Fj->MvAddMv(-one, *F_vec, one, *Fj);	
  MAVj->MvNorm(ptr_norms);
  //
  if (vb) {
    cout << " " <<  endl;
    cout << "********Block Arnoldi iteration******** " << j <<  endl;
    cout << " " <<  endl;
    for ( i=0; i<m; i++ ) {
      cout << " Arnoldi relation " << "for column " << i << " is " << fabs(ptr_norms[i]) <<  endl;
    }
    cout << " " <<  endl;
  }
  
  delete F_vec;
  delete Fj;
  delete MAVj, AVj;
  delete Vj;
  delete [] index;
  delete [] ptr_norms;
  //
} // end CheckBlkArnRed
  
  
template<class TYPE>
void BlockGmres<TYPE>::CheckGmresResids(Anasazi::MultiVec<TYPE> & x, Anasazi::MultiVec<TYPE> & b,
					Anasazi::DenseMatrix<TYPE> & z, bool vb) const {
        const TYPE one=1.0, zero=0.0;
        int i, m = (_iter+1)*_blocksize;
	Teuchos::BLAS<int,TYPE> blas;
        int *index = new int[m]; assert(index!=NULL);
        for ( i=0; i<m; i++ ) {   
                index[i] = i;
        }
        Anasazi::MultiVec<TYPE> * Vjp1 = _basisvecs->CloneView(index, m); assert(Vjp1!=NULL);
        Anasazi::MultiVec<TYPE> * Ax_copy = _basisvecs->Clone(_blocksize); assert(Ax_copy!=NULL);
        Anasazi::MultiVec<TYPE> * x_copy = x.CloneCopy(); assert(x_copy!=NULL);
	//
	//  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
	//
        Anasazi::DenseMatrix<TYPE> y_temp( z, 0, 0, m, _blocksize );
	Anasazi::DenseMatrix<TYPE> y(y_temp);
	//
	//  Solve the least squares problem and compute current solutions.
	//
  	blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		   Teuchos::NON_UNIT_DIAG, m, _blocksize, one,  
		   _hessmatrix->getarray(), _hessmatrix->getld(), y.getarray(), y.getld() );

        x_copy->MvTimesMatAddMv( one, *Vjp1, y, one );
        _amat.ApplyMatrix( *x_copy, *Ax_copy );
        x_copy->MvAddMv( one, b, -one, *Ax_copy );
	//
	//  Compute un-normalized residuals for current solutions.
	//
        TYPE *residuals = new TYPE[_blocksize]; assert(residuals!=NULL);
        x_copy->MvNorm(residuals);

        if (vb){                                                 
            cout << " " <<  endl;
            cout << "*****True GMRES Residuals***** " <<  endl;
           for (i=0; i<_blocksize; i++) {
                    cout << "True " << _rhs_iter * _blocksize + i << "-th residual "
                                  << residuals[i] <<  endl;
           }
            cout << " " <<  endl;
        }       
         
        if (Vjp1) delete Vjp1;
        if (x_copy) delete x_copy;
        if (Ax_copy) delete Ax_copy;
        if (index) delete [] index;
        if (residuals) delete [] residuals;
}
  
} // end namespace Belos
#endif
// End of file BelosBlockGmres.hpp


