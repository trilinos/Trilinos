//
// File BelosBlockCG.hpp
//
// Beta version of code put in repository on 3/21/03
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
#ifndef BLOCK_CG_HPP
#define BLOCK_CG_HPP

#include "Epetra_LAPACK.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

//
// BlockCG base class
//
template <class TYPE>
class BlockCG { 
public:
	BlockCG(AnasaziMatrix<TYPE> & mat, 
		AnasaziPrecondition<TYPE> &precond, AnasaziMultiVec<TYPE>& rhs, 
		const int numrhs, const TYPE tol=1.0e-6, const int maxits=25, 
		const int block=1, bool=false);
	virtual ~BlockCG();
	void GetSolutions(TYPE [], const int);
	void SetInitGuess(const TYPE [], const int, const int);
	void SetDebugLevel(const int);
	void SetCGBlkTols();
    void SetUpBlocks(AnasaziMultiVec<TYPE>&, AnasaziMultiVec<TYPE>&, int);
	void ExtractCurSolnBlock(AnasaziMultiVec<TYPE>&, int);
    void Solve(bool);
	void QRFactorDef(AnasaziMultiVec<TYPE>&, AnasaziDenseMatrix<TYPE>&, bool&,int,
		                 int[],int&,bool);
    void TrueResiduals(bool);
	void PrintResids(bool)const;
private:
    void CheckCGOrth(AnasaziMultiVec<TYPE>&, AnasaziMultiVec<TYPE>&, bool);
	void PrintCGIterInfo(int[], int, int[], int, int[], int);
	void CheckCGResids(AnasaziMultiVec<TYPE>&, AnasaziMultiVec<TYPE>&, bool)const;
	AnasaziMatrix<TYPE> &_amat; // must be passed in by the user
	AnasaziPrecondition<TYPE> &_precond; // must be passed in by user
	AnasaziMultiVec<TYPE> &_rhs, *_basisvecs, *_solutions, *_residvecs;
	const int _maxits, _blocksize, _numrhs;
	const TYPE _residual_tolerance;
	TYPE *_trueresids, *_residerrors;
	bool _startblock;
	int _debuglevel;
	int _rhs_iter, _iter;
	TYPE _prec, _dep_tol, _blkerror;
};

//
// Implementation
//

template <class TYPE>
BlockCG<TYPE>::BlockCG(AnasaziMatrix<TYPE> & mat, 
					         AnasaziPrecondition<TYPE> &precond, AnasaziMultiVec<TYPE>& rhs, 
							 const int numrhs, const TYPE tol, const int maxits, 
							 const int block, bool vb) : 
							_amat(mat), _precond(precond), 
							_rhs(rhs), _basisvecs(0),
							_solutions(0),_residvecs(0),
							_maxits(maxits), _blocksize(block), _numrhs(numrhs),
							_residual_tolerance(tol), _trueresids(0), _residerrors(0), 
							_startblock(false), _debuglevel(0), 
							_rhs_iter(0), _iter(0),
							_prec(5.0e-15), _dep_tol(0.75), _blkerror(1.0) { 
	//if (vb) cout << "ctor:BlockCG " << this << endl << endl;
	assert(_maxits>=0); assert(_blocksize>=0); assert(_numrhs>=0);
	//
	// Make room for the direction and residual vectors
	// We save 2 blocks of these vectors
	//
    _basisvecs = _rhs.Clone(2*_blocksize); assert(_basisvecs);
	_residvecs = _rhs.Clone(2*_blocksize); assert(_residvecs);
    //
    if (2*_blocksize && _basisvecs) {
		//
		_residerrors = new TYPE[_numrhs + _blocksize]; assert(_residerrors);
		//if (vb) cout << "BlockCG:ctor " << _maxits << _blocksize <<  _basisvecs << endl;
	}
	else {
		//if (vb) cout << "BlockCG:ctor " << _maxits << _blocksize << _basisvecs << endl;
		exit(-1);
	}
}


template <class TYPE>
BlockCG<TYPE>::~BlockCG() {
	//	cout << "dtor:BlockCG " << this << endl;
	if (_basisvecs) delete _basisvecs;
    if (_residvecs) delete _residvecs;
	if (_solutions) delete _solutions;
	if (_trueresids) delete [] _trueresids;
	if (_residerrors) delete [] _residerrors;
}



template <class TYPE>
void BlockCG<TYPE>::SetDebugLevel(const int in) {
	_debuglevel=in;
}


template <class TYPE>
void BlockCG<TYPE>::SetCGBlkTols() {
	const TYPE two = 2.0;
	TYPE eps;
	char precision = 'P';
	Epetra_LAPACK lapack;
	lapack.LAMCH(precision, eps);
	_prec = eps;
	_dep_tol = 1/sqrt(two);
}


template <class TYPE>
void BlockCG<TYPE>::SetInitGuess(const TYPE x[], const int cols, const int ldx) {
	if (_startblock==false && cols==_numrhs) {
		//
		// Set _solutions to the guesses
		//
		_solutions = _rhs.Clone(cols);
		assert(_solutions);
		_solutions->SetVecValues( x, ldx );
		_startblock = true;
	}
}
//

template<class TYPE>
void BlockCG<TYPE>::SetUpBlocks (AnasaziMultiVec<TYPE>& sol_block,  
				     AnasaziMultiVec<TYPE>& rhs_block,  
				     int num_to_solve) {
	//
	int i,j;
	int *index = new int[_blocksize + _numrhs]; assert(index);
	AnasaziMultiVec<TYPE> *tptr=0, *tptr2=0;
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
		// Put the next _blocksize of the right-hand sides into  
		// the rhs_block 
		//
		for ( i=0;i<_blocksize; i++ ) {
			index[i] = _rhs_iter*_blocksize + i;
		}
		tptr = _rhs.CloneView(index,_blocksize); assert(tptr);
		rhs_block.MvAddMv(one, *tptr, zero, *tptr);
		//
		// Put the next _blocksize of the initial guesses
		// into the sol_block
		//
		tptr = _solutions->CloneView(index,_blocksize); assert(tptr);
		sol_block.MvAddMv(one, *tptr, zero, *tptr);
		if (tptr) {
			delete tptr; tptr = 0;
		}
	}
	else {
		// More involved. This is the case where the number of right-hand
		// sides left to solve for at this iteration is less than the block size.
		//
		// Fill up the right-hand side block with random vectors, then place
		// the remaining (unsolved) right hand sides into the initial portion
		// of the right-hand side block.
		//
		rhs_block.MvRandom();
	
		for ( i=0; i<num_to_solve; i++ ) {
			index[i] = _rhs_iter*_blocksize + i;
		}
		tptr = _rhs.CloneView(index,num_to_solve); assert(tptr);
		//
		for (i=0; i<num_to_solve; i++) {
			index[i] = i;
		}
		tptr2 = rhs_block.CloneView(index, num_to_solve); assert(tptr2);
		tptr2->MvAddMv(one, *tptr, zero, *tptr);
		//
        // Fill up the sol_block and the with random vectors, then
		// place the remaining (unsolved) initial guesses into the initial portion
		// of the sol_block.
		//
		//sol_block.MvRandom();
		//
		// Fill up sol_block with zero vectors.
		int numrows = sol_block.GetVecLength();
		double * array = new double[_blocksize * numrows]; 
        for (j=0;j<_blocksize;j++){
			for (i=0; i<numrows; i++){
				array[i+j*numrows] = 0.0;
			}
		}
		sol_block.SetVecValues(array,numrows);
		//
		//
		//
		for ( i=0; i<num_to_solve; i++ ) {
			index[i] = _rhs_iter*_blocksize + i;
		}
        tptr = _solutions->CloneView(index,num_to_solve); assert(tptr);
		//
		for (i=0; i<num_to_solve; i++) {
			index[i] = i;
		}
		tptr2 = sol_block.CloneView(index, num_to_solve); assert(tptr2);
		tptr2->MvAddMv(one, *tptr, zero, *tptr);
		//
		// Delete the temporary views
		//
		if (tptr) {
			delete tptr; tptr = 0;
		}
		if (tptr2) {
			delete tptr2; tptr2 = 0;
		}
		delete [] array; array=0;
	}
	delete [] index; index=0;
	//
}
//


template <class TYPE>
void BlockCG<TYPE>::ExtractCurSolnBlock(AnasaziMultiVec<TYPE>& sol_block,
										  int num_to_solve) {
	//
	int i;
	const TYPE one = 1.0;
	const TYPE zero = 0.0;
	int * index = new int[_blocksize + _numrhs]; assert(index);
	AnasaziMultiVec<TYPE> *tptr=0, *tptr2=0;
	//
	if (num_to_solve >= _blocksize) {
		 for (i=0; i<_blocksize; i++) {
			 index[i] = _rhs_iter * _blocksize + i;
		 }
		 tptr = _solutions->CloneView(index,_blocksize); assert(tptr);
		 tptr->MvAddMv(one, sol_block, zero, sol_block);
	 }
	 else {
		 for (i=0; i<num_to_solve; i++) {
			 index[i] = _rhs_iter * _blocksize + i;
		 }
		 tptr = _solutions->CloneView(index,num_to_solve); assert(tptr);
		 //
		 for (i=0; i<num_to_solve; i++) {
			 index[i] = i;
		 }
		 tptr2 = sol_block.CloneView(index,num_to_solve); assert(tptr2);
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
void BlockCG<TYPE>::TrueResiduals (bool vb) {
	AnasaziMultiVec<TYPE> * AX = _solutions->Clone(_numrhs); assert(AX);
	//
	_amat.ApplyMatrix(*_solutions, *AX);
	assert(AX->GetNumberVecs()==_rhs.GetNumberVecs());
	
	TYPE* norms_AX=new TYPE[AX->GetNumberVecs()]; assert(norms_AX);
	AX->MvNorm(norms_AX);
	int i;
	//
	TYPE* norms_rhs=new TYPE[_rhs.GetNumberVecs()]; assert(norms_rhs);
	_rhs.MvNorm(norms_rhs);
	
	_trueresids = new TYPE[AX->GetNumberVecs()]; assert(_trueresids);
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
	cout << endl << endl;
	//
	if (AX) delete AX;
	if (norms_AX) delete [] norms_AX;
	if (norms_rhs) delete [] norms_rhs;
}


template <class TYPE>
void BlockCG<TYPE>::PrintResids(bool vb)const {
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
void BlockCG<TYPE>::GetSolutions(TYPE x[], const int ldx) {
	if (_solutions) {
		_solutions->GetVecValues(x,ldx);
	}
	else {
		assert(_solutions);
	}
}



template <class TYPE>
void BlockCG<TYPE>::Solve (bool vb) {
    //
	int i,j,k,num_ind;
	bool pflg, exit_flg = false;
	const TYPE one=1.0;
	const TYPE zero=0.0;
	AnasaziMultiVec<TYPE> *cur_block_sol=0, *cur_block_rhs=0;
	AnasaziMultiVec<TYPE> *R_prev=0, *R_new=0, *P_prev=0, *P_new=0, *AP_prev=0;
	AnasaziMultiVec<TYPE> *temp_block=0, *PC_resid;
	AnasaziMultiVec<TYPE> *precond_resid=0, *cur_sol=0;
	TYPE * ptr_T1 = 0;
	TYPE * ptr_T2 = 0;
	TYPE * cur_resid_norms=0;
	TYPE * init_resid_norms=0; 
	int *index1 = new int[_numrhs + _blocksize]; assert(index1);
    int *index2 = new int[_numrhs + _blocksize]; assert(index2);
    //
	//********************************************************************************
	//
	// max_rhs_iters is the number of times that we need to iterate in order to
	// solve for all right-hand sides (_numrhs), _blocksize at a time.
	//
	int max_rhs_iters = (_numrhs+_blocksize-1) / _blocksize;
	//
	//  Start executable statements. 
    //
	for ( _rhs_iter=0; _rhs_iter < max_rhs_iters; _rhs_iter++ ) {
		//
		if (vb && _debuglevel > 2) {
		   cout << endl;
		   cout << "===================================================" << endl;
		   cout << "RHS pass: " << _rhs_iter << endl;
		   cout << endl;
		}	
		//
		int cur_blksz = _blocksize;
	    int ind_blksz = _blocksize;
		int prev_ind_blksz = _blocksize;
		int num_conv = 0;
		//
		// If not provided, and if _rhs_iter = 0, 
		// set initial guesses to AX = B to random vectors
		//
		if (!_startblock) {
			_solutions = _rhs.Clone(_numrhs); assert(_solutions);
			_solutions->MvRandom();
			_startblock = true;
		}
		
		// Make room for current blocks of solutions and right-hand sides
		//
		 cur_block_sol = _solutions->Clone(_blocksize); assert(cur_block_sol);
		 cur_block_rhs = _solutions->Clone(_blocksize); assert(cur_block_rhs);
		 //
		 // Compute the number of right-hand sides remaining to be solved 
	     //
		 int numrhs_to_solve;
	     if ( _blocksize < _numrhs ) {
		     numrhs_to_solve = _numrhs - (_rhs_iter * _blocksize);
		 }
	     else {
		     numrhs_to_solve = _numrhs;
		 }
		 //
		 // Put the current initial guesses and right-hand sides into current blocks
		 //
         SetUpBlocks(*cur_block_sol, *cur_block_rhs, numrhs_to_solve);
         //
		 // Make additional space needed during iteration
		 //
         temp_block = _solutions->Clone(_blocksize); assert(temp_block);
		 PC_resid = _solutions->Clone(_blocksize); assert(PC_resid);
         //
		 // Additional initialization
		 //
		 int *ind_idx = new int[_blocksize]; assert(ind_idx);
		 int *cur_idx = new int[_blocksize]; assert(cur_idx);
		 int *conv_idx = new int[_blocksize]; assert(conv_idx);
		 int *cols = new int[_blocksize]; assert(cols);
		 //
		 for (i=0;i<_blocksize;i++){
			 ind_idx[i] = i; cur_idx[i] = i; conv_idx[i] = 0;
			 cols[i] = i;
		 }
		//
		// ************ Compute the initial residuals ********************************
		//
		// Associate the first block of _basisvecs with P_prev and the
        // first block of _residvecs with R_prev
		//
		for ( i=0; i<_blocksize; i++ ) {
			index1[i] = i;
		}
		P_prev = _basisvecs->CloneView(index1, _blocksize); assert(P_prev);
        R_prev = _residvecs->CloneView(index1, _blocksize); assert(R_prev);
		AP_prev = temp_block->CloneView(index1, _blocksize); assert(AP_prev);
		cur_sol = cur_block_sol->CloneView(index1, _blocksize); assert(cur_sol);
        //
        // Store initial guesses to AX = B in 1st block of _basisvecs
        //         P_prev = one*cur_block_sol + zero*P_prev
		//
		P_prev->MvAddMv(one, *cur_block_sol, zero, *P_prev);
		//
		// Multiply by A and store in AP_prev
		//       AP_prev = A*P_prev
		//
		_amat.ApplyMatrix( *P_prev,*AP_prev );
		//
		// Compute initial residual block and store in 1st block of _residvecs
		//     R_prev = cur_block_rhs - A*P_prev
		//
		R_prev->MvAddMv(one, *cur_block_rhs, -one, *AP_prev);
		//
		//*******Compute and save the initial residual norms*******
		//
		init_resid_norms = new TYPE[_blocksize]; assert(init_resid_norms);
        R_prev->MvNorm(init_resid_norms);
		//
        // Update indices of current (independent) blocks.
		// If a residual is too small, it will be dropped from
		// the current block, thus, from future computations
		//
		k = 0; j = 0;
		for (i=0; i<_blocksize; i++){
			_residerrors[_rhs_iter*_blocksize +i] = init_resid_norms[i];
			if (init_resid_norms[i] > _prec){
				_residerrors[_rhs_iter*_blocksize +i] = one;
				cur_idx[k] = i; ind_idx[k] = i;
				k = k+1;
			}
			else {
				conv_idx[j] = i;
				j = j+1;
			}
		}
		cur_blksz = k; ind_blksz = k; num_conv = j;
		//
        //***** If _debuglevel > 0, print out information****************************
        //
        if (vb && _debuglevel > 0) {
			cout << endl;
            cout << " CG Initial Residual Norms" << endl; 
			for (i=0; i<_blocksize; i++){
				cout << _residerrors[_rhs_iter*_blocksize + i] << " ";
			}
            cout << endl << endl;
		}
	    //
		//
		if (cur_blksz > 0) { 
			_blkerror = 1.0;
			// All initial residuals have not converged -- continue Block CG	
			//Compute the initial block of direciton vectors
			//
		    // Associate current blocks of residuals, directions, and precond residuals
		    // with R_prev, P_prev, and precond_resid
		    //
		    for (i=0; i< cur_blksz; i++){
			    index1[i] = cur_idx[i];
			}
		    R_prev = _residvecs->CloneView(index1,cur_blksz);
		    P_prev = _basisvecs->CloneView(index1,cur_blksz);
		    precond_resid = PC_resid->CloneView(index1,cur_blksz);
		    //
		    // Compute the preconditioned initial residual, store in precond_resid
		    //
	    	_precond.ApplyPrecondition(*R_prev, *precond_resid);
		    //
		    //**************Compute initial direction vectors************************
			// Initially, they are set to the preconditioned residuals
		    //
		    P_prev->MvAddMv(one, *precond_resid, zero, *precond_resid);
		    //
		    // Compute an orthonormal block of initial direction vectors,
			// and check for dependencies, adjusting indices of independent
			// vectors if needed
		    //
            AnasaziDenseMatrix<TYPE> G(cur_blksz, cur_blksz);
		    num_ind = 0; pflg = false;
		    QRFactorDef(*P_prev, G, pflg, cur_blksz, cols, num_ind, vb);
            //
	 	    if (pflg) {
			    // The initial block of direction vectors are linearly dependent
			    if (vb && _debuglevel > 2) {
			       cout << " Initial direction vectors are dependent" << endl;
			       cout << " Adjusting blocks and indices for iteration" << endl;
				}
			    //
			    ind_blksz = num_ind;
			    for (i=0; i< ind_blksz; i++){
				    ind_idx[i] = ind_idx[cols[i]];
				}	
			}  // end if (pflg)
			//
			if (ind_blksz == 0){
				if (vb) {
					cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
					cout << " Reason: No independent initial direction vectors" << endl;
				}
				exit_flg = true;
			}

		}  // end if (cur_blksz > 0)
	    //		
	    else {  // all initial residuals have converged
			if (vb) {
				cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				cout << "Reason: All initial residuals have converged" << endl;
			}
			exit_flg = true;
		}
        // ***************************************************************************
        // ************************Main CG Loop***************************************
		// ***************************************************************************
        // 
		int new_blk = 2;
        if (vb && _debuglevel > 2) cout << "Entering main CG loop" << endl << endl;
	    //
	    for (_iter=0; _iter<_maxits; _iter++) {
			//
			if (exit_flg) break;
		    // 
		   //*******Compute the new blocks of iterates and residuals******************
		   //
		   // Get views of the previous blocks of residuals, direction vectors, etc.
		   //
		   for (i=0; i< ind_blksz; i++) {
			   index1[i] = ind_idx[i];
			   index2[i] = _blocksize + ind_idx[i];
		   } 
		   if (new_blk == 2){
             P_prev = _basisvecs->CloneView(index1,ind_blksz);
		   }
		   else {
			 P_prev = _basisvecs->CloneView(index2,ind_blksz);
		   }
		   AP_prev = temp_block->CloneView(index1,ind_blksz);
		   //
		   for (i=0; i < cur_blksz; i++){
			   index1[i] = cur_idx[i];
               index2[i] = _blocksize + cur_idx[i];
		   }
		   if (new_blk == 2){
		      R_prev = _residvecs->CloneView(index1,cur_blksz);
			  R_new = _residvecs->CloneView(index2,cur_blksz);
		   }
		   else {
			  R_prev = _residvecs->CloneView(index2,cur_blksz);
			  R_new = _residvecs->CloneView(index1,cur_blksz);
		   }
		   cur_sol = cur_block_sol->CloneView(index1,cur_blksz); 
           //
		   // Compute the coefficient matrix alpha
		   //
           // P_prev^T * A * P_prev * alpha = P_prev^T * R_prev
		   // 1) Compute P_prev^T * A * P_prev = T2 and P_prev^T * R_prev = T1
           // 2) Compute the Cholesky Factorization of T2
           // 3) Back and Forward Solves for alpha
           //
		   AnasaziDenseMatrix<TYPE> alpha(ind_blksz,cur_blksz);
           AnasaziDenseMatrix<TYPE> T1(ind_blksz,cur_blksz);
           AnasaziDenseMatrix<TYPE> T2(ind_blksz,ind_blksz);
		   ptr_T1 = 0; ptr_T2 = 0;
		   char UPLO = 'U';
           int ii = 0;
           int * info = &ii;
           int numrhs = _numrhs;
           Epetra_LAPACK lapack;
		   //
		   // 1)
		   _amat.ApplyMatrix(*P_prev, *AP_prev);
		   R_prev->MvTransMv(one, *P_prev, T1);
           AP_prev->MvTransMv(one, *P_prev, T2);
		   // 2)
		   ptr_T1 = T1.getarray(); assert(ptr_T1);
           ptr_T2 = T2.getarray(); assert(ptr_T2);
           lapack.POTRF(UPLO, ind_blksz, ptr_T2, ind_blksz, info);
		   if (*info != 0) {
			   if(vb){
				   cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				   cout << " Reason: Cannot compute coefficient matrix alpha" << endl;
				   cout << " P_prev'* A*P_prev is singular" << endl;
				   cout << " Solution will be updated upon exiting loop" << endl;
			   }
			   break;
		   }
		   // 3)
		   lapack.POTRS(UPLO, ind_blksz, cur_blksz, ptr_T2, ind_blksz, ptr_T1, ind_blksz, info);
		   // Note: solution returned in ptr_T1
		   if (*info != 0) {
			   if(vb){
				   cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				   cout << " Reason: Cannot compute coefficient matrix alpha" << endl;
				   cout << " Solution will be updated upon exiting loop" << endl;
			   }
			   break;
		   }
           //
           alpha.setvalues(ptr_T1,ind_blksz);
           //
           // Update the solution: cur_sol = one*cur_sol + one*P_prev*alpha
		   //
           cur_sol->MvTimesMatAddMv(one, *P_prev, alpha, one);
           //
           // Update the residual vectors: R_new = R_prev - A*P_prev*alpha
           //
           R_new->MvAddMv(one, *R_prev, zero, *R_prev);
           R_new->MvTimesMatAddMv(-one, *AP_prev, alpha, one);
		   //
		   // ****Compute the Current Relative Residual Norms and the Block Error****
           //
		   cur_resid_norms = new TYPE[cur_blksz]; assert(cur_resid_norms);
           R_new->MvNorm(cur_resid_norms);
		   for (i=0; i<cur_blksz; i++){
			   cur_resid_norms[i] = cur_resid_norms[i]/init_resid_norms[cur_idx[i]];
		   }
		   // Update _residerrors
		   for (i=0; i<cur_blksz; i++){
			   _residerrors[_rhs_iter*_blocksize + cur_idx[i]] = cur_resid_norms[i];
		   }
		   //
		   if (cur_resid_norms) {
			   delete [] cur_resid_norms; cur_resid_norms=0;
		   }  
		   prev_ind_blksz = ind_blksz; // Save old ind_blksz of P_prev
		   //
		   // Update the number of residuals in current block, their indices,
		   // and the block error
		   //
		   k=0;
		   _blkerror = _residerrors[_rhs_iter*_blocksize +cur_idx[0]];
		   for (i=0; i<cur_blksz; i++){
			   if (_residerrors[_rhs_iter*_blocksize + cur_idx[i]] > _blkerror){ // get new max error
				   _blkerror = _residerrors[_rhs_iter*_blocksize + cur_idx[i]];
			   }
			   if (_residerrors[_rhs_iter*_blocksize + cur_idx[i]] > _prec){ // get new cur_idx
				   cur_idx[k] = cur_idx[i]; k = k+1;
			   }
		   }
		   cur_blksz = k;
		   //
		   // Update the number of current residuals that correspond
		   // to linearly independent direction vectors. Note that
		   // ind_idx are a subset of cur_idx.
		   //
		   k = 0;
		   for (i=0; i< ind_blksz; i++){
               if (_residerrors[_rhs_iter*_blocksize +ind_idx[i]] > _prec){
				   ind_idx[k] = ind_idx[i]; k = k+1;
			   }
		   }
		   ind_blksz = k;
           //
		   // Update the number of converged residuals and their indices 
		   //
		   k = 0;
		   for (i=0; i<_blocksize; i++){
               if (_residerrors[_rhs_iter*_blocksize + i] < _residual_tolerance){
				   conv_idx[k] = i; k = k+1;
			   }
		   }
		   num_conv = k;
		   //
		   //****************Print iteration information*****************************
		   //
		   if (vb && _debuglevel > 0) {
			   PrintCGIterInfo(cur_idx,cur_blksz,
				                      ind_idx,ind_blksz,conv_idx,num_conv);
		   }
		   //
		   if (_debuglevel > 1) {
			   CheckCGResids(*cur_block_sol, *cur_block_rhs, vb);
		   }
		   //
		   //****************Test for convergence*************************************
		   //
		   if (_blkerror <= _residual_tolerance) {
			   if (vb) {
				   cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				   cout << " Reason: Block CG has converged" << endl; 
			   }
			   break;
		   }
		   //
		   if (ind_blksz <= 0){
				if (vb) {
					cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				   cout << " Reason: No more independent direction vectors" << endl;
				   cout << " Solution will be updated upon exiting loop" << endl;
				}
			    break;
			}
		   //
		   // **************Compute the new block of direction vectors****************
		   //
           // Get views of the new blocks of independent direction vectors and
		   // the corresponding residuals. Note: ind_idx are a subset of cur_idx.
		   //
		   for (i=0; i<ind_blksz; i++){
			   index1[i] = ind_idx[i]; index2[i] = _blocksize + ind_idx[i];
		   }
		   if (new_blk == 2) {
		      R_new = _residvecs->CloneView(index2,ind_blksz);
		      P_new = _basisvecs->CloneView(index2,ind_blksz);
		   }
		   else {
              R_new = _residvecs->CloneView(index1,ind_blksz);
		      P_new = _basisvecs->CloneView(index1,ind_blksz);
		   }
		   precond_resid = PC_resid->CloneView(index1,ind_blksz);
		   //
           // Compute preconditioned residuals
		   //
		   _precond.ApplyPrecondition(*R_new, *precond_resid);
		   // 
		   // Compute coefficient matrix beta
		   //
           // P_prev^T A * P_prev * beta = P_prev^T A * precond_resid
           // 1) Compute P_prev^T A * P_prev = T2 and P_prev^T * A * precond_resid = T3
		   //                                     or (A*P_prev)^T * precond_resid (A SPD)
           // 2) Compute the Cholesky Factorization of T2
           // 3) Back and Forward Solves for beta
           //
           AnasaziDenseMatrix<TYPE> T3(prev_ind_blksz,ind_blksz);
           AnasaziDenseMatrix<TYPE> beta(prev_ind_blksz,ind_blksz);
		   ptr_T1 = 0;
		   //
		   // 1 & 2)  Note: we already have computed T2 and its Cholesky
		   //         factorization during computation of alpha
		   precond_resid->MvTransMv(-one, *AP_prev, T3);
		   // 3)
           ptr_T1 = T3.getarray(); assert(ptr_T1);  
           lapack.POTRS(UPLO, prev_ind_blksz, ind_blksz, ptr_T2, prev_ind_blksz, ptr_T1, prev_ind_blksz, info);
		   // Note: Solution returned in ptr_T1
		   if (*info != 0) {
			   if (vb) {
				   cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
			      cout << "Reason: Cannot compute coefficient matrix beta" << endl;
			      cout << "Solution will be updated upon exiting loop" << endl;
			   }
			   break;
		   }
           //
           beta.setvalues(ptr_T1,prev_ind_blksz);
           //
           // Compute: P_new = precond_resid + P_prev * beta
           //
           P_new->MvAddMv(one, *precond_resid, zero, *R_new);
           P_new->MvTimesMatAddMv(one, *P_prev, beta, one);
           //
           // Check A-orthogonality of new and previous blocks of direction vectors
           //
		   if (_debuglevel > 2) {
			   if(vb){
		         cout << "Orthogonality check" << endl;
			   }
		       CheckCGOrth(*P_prev, *P_new, vb);   
		   }
		   //
		   // Compute orthonormal block of direction vectors,
		   // and check for dependencies, adjusting indices of
		   // independent vectors if needed
		   //
		   AnasaziDenseMatrix<TYPE> G(ind_blksz,ind_blksz);
		   num_ind = 0; pflg = false;
		   QRFactorDef(*P_new, G, pflg, ind_blksz, cols, num_ind, vb);
		   //
		   if (pflg) {
			  ind_blksz = num_ind;
			  // The new block of direction vectors are linearly dependent
			  if (vb && _debuglevel > 2) {
			     cout << "The new block of direction vectors are dependent " << endl;
                 cout << "# independent direction vectors: " << ind_blksz << endl;
			     cout << "independent indices: " << endl;
		         for (i=0; i<ind_blksz ; i++) {
			        cout << cols[i];
				 }
		         cout << endl << endl;
			  }
			  if (ind_blksz == 0){
				  if (vb) {
					  cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
				      cout << "Reason: No more linearly independent direction vectors" << endl;
				      cout << " Solution will be updated upon exiting loop" << endl;
				  }
			      break;
			  }
			  for (i=0; i<ind_blksz; i++){
				  ind_idx[i] = ind_idx[cols[i]];
			  }
			  //
           }  // end of if (pflg)
		   //
		   // Check A-orthogonality after orthonormalization
		   //
		   if (_debuglevel > 2) {
			   if(vb){
				  cout << "Orthogonality check after orthonormalization " << endl; 
			   }
				  CheckCGOrth(*P_prev, *P_new, vb);
		   }
		   // 
		   // *****Update index of new blocks*****************************************
		   //
           if (new_blk == 2){
			   new_blk = 1;
		   }
		   else {
			   new_blk = 2;
		   }
		   //
     } // end of the main CG loop -- for(_iter = 0;...)
	 //*******************************************************************************
	 //
	 // Iteration has converged, maxits has been reached,
	 // or direction vector block is linearly dependent
	 //
	 if (_iter >= _maxits && _blkerror > _residual_tolerance){
		 if (vb){
			 cout << " Exiting Block CG iteration -- RHS pass# " << _rhs_iter 
						<< " -- Iteration# " << _iter << endl;
			 cout << "Reason: maximum number of iterations has been reached" << endl;
		 }
	 }
	 //
	 if(vb){
	   cout << " Iteration has stopped at step: " << _iter << endl;
	   cout << endl;
	 }
	 //
	 // Insert the current block of solutions into _solutions so we can
	 // continue if we have more right-hand sides to solve for
	 //
	 ExtractCurSolnBlock(*cur_block_sol, numrhs_to_solve);	 
     //
     // ****************Free heap space***********************************************
	 //
	 if (P_prev) {delete P_prev; P_prev=0;}
	 if (P_new) { delete P_new; P_new=0;}
     if (R_prev) { delete R_prev; R_prev=0;}
     if (R_new) { delete R_new; R_new=0;}
	 if (AP_prev) { delete AP_prev; AP_prev=0;}
	 if (PC_resid) {delete PC_resid; PC_resid=0;}
	 if (precond_resid) {delete precond_resid; precond_resid=0;}
	 if (temp_block) { delete temp_block; temp_block=0;}
     if (cur_block_sol) { delete cur_block_sol; cur_block_sol=0;}
	 if (cur_sol) { delete cur_sol; cur_sol=0;}
	 if (cur_block_rhs) { delete cur_block_rhs; cur_block_rhs=0;}
	 if (ind_idx) { delete [] ind_idx; ind_idx=0;}
	 if (cur_idx) { delete [] cur_idx; cur_idx=0;}
	 if (cols) { delete [] cols; cols=0;}
	 if (init_resid_norms) { delete [] init_resid_norms; init_resid_norms=0;}
	 //
  } // end if ( _rhs_iter = 0;... )
  //**********************************************************************************
  //
  if (index1) { delete [] index1; index1=0;}
  if (index2) { delete [] index2; index2=0;}
  //
} // end CGSolve()
//


template<class TYPE>
void BlockCG<TYPE>::QRFactorDef (AnasaziMultiVec<TYPE>& VecIn, 
				     AnasaziDenseMatrix<TYPE>& FouierR, bool &flg, int blksz,
					 int cols[], int &num, bool vb) {
	int i,j,k;
	int nb = VecIn.GetNumberVecs(); assert (nb == blksz);
	int ldR = FouierR.getld();
	int *index = new int[nb]; assert(index);
	int *dep_idx = new int[nb]; assert(dep_idx);
	int num_dep = 0;
	TYPE * R = FouierR.getarray();
	AnasaziMultiVec<TYPE> *qj = 0, *Qj = 0;
	const int IntOne=1;
	const int IntZero=0;
	const TYPE zero=0.0;
	const TYPE one=1.0;
	TYPE * NormVecIn = new TYPE[nb]; assert(NormVecIn);
	//
	//
    // Zero out the array that will contain the Fourier coefficients.
	//
	for ( j=0; j<nb; j++ ) {
		for ( i=0; i<nb; i++ ) {
			R[j*ldR+i] = zero;
		}
	}
	//
	// Compute the Frobenius norm of VecIn -- this will be used to
	// determine rank deficiency of VecIn
	//
	VecIn.MvNorm(NormVecIn);
	TYPE FroNorm = 0.0;
	for (j=0;j<nb;j++) {
		FroNorm += NormVecIn[j] * NormVecIn[j];
	}
	FroNorm = sqrt(FroNorm);
	num = 0; flg = false;
	//
    // Start the loop to orthogonalize the nb columns of VecIn.
	//
	for ( j=0; j<nb; j++ ) {
		//
        // Grab the j-th column of VecIn (the first column is indexed to 
        // be the zero-th one).
		//
		index[0] = j;
		qj = VecIn.CloneView(index, IntOne); assert(qj);
		if ( j ) {
			int num_orth;
			for ( i=0; i<j; i++ ) {
				index[i] = i;
			}
			//
			// Grab the first j columns of VecIn (that are now an orthogonal
			// basis for first j columns of the entering VecIn).
			//
			Qj = VecIn.CloneView(index, j);
			AnasaziDenseMatrix<TYPE> rj(j,1);
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
				qj->MvTransMv( one, *Qj, rj );
				//
				// Sum result[0:j-1] into column j of R.
				//
				TYPE* result = rj.getarray();
				//
				for (k=0; k<num_dep; k++) {
					result[dep_idx[k]] = zero;
				}
				//
				for ( k=0; k<j; k++ ) {
					R[j*ldR+k] += result[k];
				}
				//
				//   Compute qj <- qj - Qj * rj.
				//
				qj->MvTimesMatAddMv(-one, *Qj, rj, one);
			}
		}
		//
		// Compute the norm of column j of VecIn (=qj).
		//
		qj->MvNorm ( R+j*ldR+j );
	    //
		if ( NormVecIn[j] > _prec && *(R+j*ldR+j) > (_prec * NormVecIn[j]) ) {
			//
			// Normalize qj to make it into a unit vector.
			//
			TYPE rjj = one / *(R+j*ldR+j);
			qj->MvAddMv ( rjj, *qj, zero, *qj );
			cols[num] = j;
			num = num + 1;
		}
		else {
			// 
			if (vb && _debuglevel > 2){
			    cout << "Rank deficiency at column index: " << j << endl;
			}
			flg = true;
			//
			// Don't normalize qj, enter one on diagonal of R,
			// and zeros in the row to the right of the diagonal -- this
			// requires updating the indices of the dependent columns
			//
			R[j*ldR+j] = one;
			dep_idx[num_dep] = j;
			num_dep = num_dep + 1;
		}	
		delete qj; qj = 0;
		delete Qj; Qj = 0;
	}
	assert ((num + num_dep) == blksz);
	delete [] index;
	delete [] dep_idx;
	if (NormVecIn) delete [] NormVecIn;
	//
} // end QRFactorDef


template<class TYPE>
void BlockCG<TYPE>::CheckCGOrth(AnasaziMultiVec<TYPE>& P1, AnasaziMultiVec<TYPE>& P2,
								bool vb) {
  //
  // This routine computes P2^T * A * P1
  // Checks the orthogonality wrt A between any two blocks of multivectors with the same length
  //
  const TYPE one = 1.0;
  const TYPE zero = 0.0;
  int i, k;
  int veclen1 = P1.GetVecLength();
  int veclen2 = P2.GetVecLength();
  assert(veclen1 == veclen2);
  //
  int numvecs1 = P1.GetNumberVecs();
  int numvecs2 = P2.GetNumberVecs();
  //
  AnasaziMultiVec<TYPE>* AP = P1.CloneCopy();
  assert(AP);
  _amat.ApplyMatrix(P1, *AP);
  //
  AnasaziDenseMatrix<TYPE> PAP(numvecs2, numvecs1);
  AP->MvTransMv(one, P2, PAP);
  //
  TYPE* ptr = PAP.getarray();
  TYPE column_sum;
  //
  for (k=0; k<numvecs1; k++) {
	  column_sum = zero;
	  for (i=0; i<numvecs2; i++) {
		  column_sum += ptr[i];
	  }
	  if (vb) {
	     cout << " P2^T*A*P1 " << " for column "
		      << k << " is  " << fabs(column_sum) << endl;
	  }
	  ptr += numvecs2;
  }
  if (vb) {
     cout << " " << endl;
  }
  //
  //PAP.DisplayMat();

  if(AP) {
    delete AP; AP=0;
  }  
//  
} // end check_orthog
//


template<class TYPE>
void BlockCG<TYPE>::PrintCGIterInfo(int cur[], const int cursz, 
									  int ind[], const int indsz,
									  int conv[], const int convsz) {
	//
	int i;
	if (_debuglevel > 0){
		cout << "--------------------------------------------------" << endl;
		cout << endl;
        cout << " CG Residual Norms -- Iteration# " << _iter 
			 << "  RHS pass# " << _rhs_iter << endl;
		for (i=0;i<_blocksize;i++){
			cout << "_residerrors[" << _rhs_iter*_blocksize+i << "] = " 
				<< _residerrors[_rhs_iter*_blocksize + i] << endl;
		}
		cout << endl;
		cout << "Maximum Current Residual Error:  " << _blkerror << endl;
		cout << endl;
	}  // end if (_debuglevel > 0)
	if (_debuglevel > 2){	
		cout << "# of current residuals in block: " << cursz << endl;
		cout << "# of converged residuals in block: " << convsz << endl;
        cout << "# of independent direction vectors: " << indsz << endl;    
		cout << " Current indices: " << endl;
		for (i=0; i<cursz; i++){
			cout << cur[i] << " ";
		}
        cout << endl;
        cout << " Converged indices: " << endl;
		for (i=0; i<convsz; i++){
			cout << conv[i] << " ";
		}
        cout << endl;
        cout << " Independent indices: " << endl;
		for (i=0; i<indsz; i++){
			cout << ind[i] << " ";
		}
        cout << endl << endl;
	} // end if (_debuglevel > 1)
	//
} // end Print_CGiter_info


template<class TYPE>
void BlockCG<TYPE>::CheckCGResids(AnasaziMultiVec<TYPE>& X, AnasaziMultiVec<TYPE>& B, bool vb) const {
	//
	int i;
	const TYPE one = 1.0;
	//
	AnasaziMultiVec<TYPE> *AX = X.CloneCopy(); assert(AX);
	AnasaziMultiVec<TYPE> *TR = X.CloneCopy(); assert(TR);
	//
	_amat.ApplyMatrix(X, *AX);
	TR->MvAddMv(one, B, -one, *AX);
	//
	TYPE *resids = new TYPE[_blocksize]; assert(resids);
	TR->MvNorm(resids);
	//
	if (vb){
		cout << endl;
		cout << "*****True CG Residual Norms*****" << endl;
		for (i=0; i<_blocksize; i++){
			cout << "True " << _rhs_iter*_blocksize+i << "-th residual  "
				 << resids[i] << endl;
		}
		cout << endl;
	}
	//
	if (AX) delete AX;
	if (TR) delete TR;
	if (resids) delete [] resids;
	//
} // end CheckCGResids()

//
//
#endif
// End of file BelosBlockCG.hpp
