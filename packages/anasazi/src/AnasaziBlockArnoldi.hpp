// File BlockArnoldi.hpp
#ifndef BLOCK_ARNOLDI_HPP
#define BLOCK_ARNOLDI_HPP

#include "AnasaziLAPACK.hpp"
#include "AnasaziBLAS.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziCommon.hpp"

/*!	\class Anasazi::BlockArnoldi

	\brief This class implements the Implicitly Restarted Block Arnoldi Method,
	an iterative method for solving eigenvalue problems.

	\author Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

template <class TYPE>
class BlockArnoldi { 
public:
	//@{ \name Constructor/Destructor.
	//! %Anasazi::BlockArnoldi constructor.
	BlockArnoldi( AnasaziMatrix<TYPE> & mat, AnasaziMultiVec<TYPE>& vec,
		const TYPE tol=1.0e-6, const int nev=5, const int length=25, 
		const int block=1, const string which="LM", const int step=25, 
		const int restarts=0 );
	//! %Anasazi::BlockArnoldi destructor.
	virtual ~BlockArnoldi();
	//@}

	//@{ \name Solver application methods.

	/*! \brief This method performs a given number of steps of the Block Arnoldi
		Method, returning upon completion or convergence.
	*/
	void iterate( const int steps=1 );

	/*! \brief This method uses iterate to compute approximate solutions to the
		original problem.  It may return without converging if it has taken the
		maximum number of iterations or numerical breakdown is observed.
	*/
	void solve();
	//@}

	//@{ \name Solution return methods.

	//! This method puts the real part of the computed eigenvectors in %evecs.
	void getEvecs(AnasaziMultiVec<TYPE>& evecs); 

	//! This method puts the imaginary part of the computed eigenvectors in %ievecs.
	void getiEvecs(AnasaziMultiVec<TYPE>& ievecs);

	//! This method returns the real part of the computed eigenvalues.
	TYPE * getEvals();

	//! This method returns the imaginary part of the computed eigenvalues.
	TYPE * getiEvals();

	//! This method returns the residuals for the computed eigenpairs.
	TYPE * getResiduals();
	//@}
	
	//@{ \name Problem attribute method.

	/*! \brief This method allows the user to inform the solver of a problem's symmetry.
		Some computational work can be avoided by setting this properly in the
		symmetric case.
	*/
	void setSymmetric( const bool );
	//@}

	//@{ \name Output methods.

	/*! \brief This method allows the user to set the solver's level of visual output
		during computations.
	*/
	void setDebugLevel( const int );

	//! This method requests that the solver print out its current status to screen.
	void currentStatus();
	//@}
private:
	void QRFactorization( AnasaziMultiVec<TYPE>&, AnasaziDenseMatrix<TYPE>& );
	void BlockReduction();
	void BlkOrth( const int j );
	void BlkOrthSing( const int j );
	void ComputeResiduals( const bool );
	void ComputeEvecs();
	void Restart();
	void Sort( const bool );
	void SetInitBlock();
	void SetBlkTols();
	void CheckBlkArnRed( const int );
	AnasaziMatrix<TYPE> &_amat; // must be passed in by the user
	AnasaziMultiVec<TYPE> &_ivec; // must be passed in by the user
	AnasaziMultiVec<TYPE> *_basisvecs, *_evecr, *_eveci;
	AnasaziDenseMatrix<TYPE>* _hessmatrix;
	const int _nev, _length, _block, _restarts, _step;
	const TYPE _residual_tolerance;
	TYPE *_ritzresiduals, *_actualresiduals, *_evalr, *_evali;
	int _restartiter, _iter, _jstart, _jend, _nevblock, _debuglevel, _nconv, _defblock;
	bool _initialguess, _issym, _isdecompcurrent, _isevecscurrent, exit_flg, dep_flg;
	TYPE _schurerror, _dep_tol, _blk_tol, _sing_tol, _def_tol;
	string _which;
	int *_order;
};
//
// Implementation
//
// Note: I should define a copy constructor and overload = because of the use of new
//
template <class TYPE>
BlockArnoldi<TYPE>::BlockArnoldi(AnasaziMatrix<TYPE> & mat, AnasaziMultiVec<TYPE> & vec, 
				const TYPE tol, const int nev, const int length, const int block,
				const string which, const int step, const int restarts) : 
				_amat(mat), _ivec(vec), _basisvecs(0), _evecr(0), _eveci(0), 
				_hessmatrix(0), _nev(nev), _length(length), _block(block), 
				_restarts(restarts), _residual_tolerance(tol), _step(step),
				_ritzresiduals(0), _evalr(0), _evali(0), _restartiter(0), 
				_iter(0), _jstart(0), _jend(0), _which(which),
				_initialguess(true), _debuglevel(0), _nevblock(0), _defblock(0),
				_issym(false), _nconv(0), _schurerror(1.0), dep_flg(false),
				_dep_tol(0), _sing_tol(0), _blk_tol(0), _def_tol(0), 
				exit_flg(false), _isevecscurrent(false), _isdecompcurrent(false) {
	//	cout << "ctor:BlockArnoldi " << this << endl;
	//
	// Determine _nevblock : how many blocks it will take to contain the _nev eigenvalues/vectors
	//
	_nevblock = _nev/_block;
	if (_nev%_block) { 
		// Another block is needed to contain the _nev eigenvalues
		_nevblock++;
	}
	assert(_length>=0); assert(_block>=0);
	//
	// Make room for the Arnoldi vectors and F.
	//
	_basisvecs = _ivec.Clone((_length+1)*_block); assert(_basisvecs);
	//
	// Make room for the eigenvectors
	//
	_evecr = _ivec.Clone(_nev); assert(_evecr);
	_eveci = _ivec.Clone(_nev); assert(_eveci);
	if (_length*_block && _basisvecs && _evecr && _eveci) {
		//
		// Create the rectangular Hessenberg matrix
		//
		_hessmatrix = new AnasaziDenseMatrix<TYPE>((_length+1)*_block, _length*_block); 
		assert(_hessmatrix);
		//
		// Create the vectors for eigenvalues and their residual errors and
		// initialize them.
		//
		_evalr = new TYPE[ _block*_length ]; assert(_evalr);  
		_evali = new TYPE[ _block*_length ]; assert(_evali);  
		_ritzresiduals = new TYPE[ _block*_length ]; assert(_ritzresiduals);
		_actualresiduals = new TYPE[ _block*_length ]; assert(_actualresiduals);
		_order = new int[ _block*_length ]; assert(_order);
		const TYPE one = 1.0, zero = 0.0;
		for (int i=0; i< _block*_length; i++) {
			_evalr[i] = zero; _evali[i] = zero;
			_ritzresiduals[i] = one;
			_actualresiduals[i] = one;
		}			
		//
		//  Set the tolerances for block orthogonality
		//
		SetBlkTols();  
	}
	else {
		cout << "BlockArnoldi:ctor " << _length << _block << _basisvecs << endl;
		exit(-1);
	}
}

template <class TYPE>
void BlockArnoldi<TYPE>::SetBlkTols() {
        const TYPE two = 2.0;
        TYPE eps;
        char precision = 'P';
        AnasaziLAPACK lapack;
        lapack.LAMCH(precision, eps);
        _blk_tol = 10*sqrt(eps);
        _sing_tol = 10 * eps;
        _dep_tol = 1/sqrt(two);
	_def_tol = eps;
}

template <class TYPE>
BlockArnoldi<TYPE>::~BlockArnoldi() {
	//	cout << "dtor:BlockArnoldi " << this << endl;
	if (_basisvecs) delete _basisvecs;
	if (_evecr) delete _evecr;
	if (_eveci) delete _eveci;
	if (_hessmatrix) delete _hessmatrix;
	if (_ritzresiduals) delete [] _ritzresiduals;
	if (_evalr) delete [] _evalr;
	if (_evali) delete [] _evali;
	if (_order) delete [] _order;
}

template <class TYPE>
void BlockArnoldi<TYPE>::getEvecs(AnasaziMultiVec<TYPE>& evecs) {
	//
	//  Compute the current eigenvectors if they are not current.
	//
	if (!_isevecscurrent) { ComputeEvecs(); }
        int i, numvecs = evecs.GetNumberVecs(); 
        if (numvecs > _nev) {
                numvecs = _nev;
        } 
        int* index = new int[ numvecs ];
        for (i=0; i<numvecs; i++) {
                index[i] = i;
        }
        evecs.SetBlock( *_evecr, index, numvecs );
                 
        delete [] index;
}
 
template <class TYPE>
void BlockArnoldi<TYPE>::getiEvecs(AnasaziMultiVec<TYPE>& ievecs) {
	//
	//  Compute the current eigenvectors if they are not current.
	//
	if (!_isevecscurrent) { ComputeEvecs(); }
        int i, numvecs = ievecs.GetNumberVecs(); 
        if (numvecs > _nev) {
                numvecs = _nev;
        } 
        int* index = new int[ numvecs ];
        for (i=0; i<numvecs; i++) {
                index[i] = i;
        }
        ievecs.SetBlock( *_eveci, index, numvecs );
                 
        delete [] index;
}

template <class TYPE>
TYPE * BlockArnoldi<TYPE>::getEvals() {
	int i;
	TYPE *temp_evals = new TYPE[ _nev ];
	for (i=0; i<_nev; i++) {
		temp_evals[i] = _evalr[i];
	}
	return temp_evals;		
}

template <class TYPE>
TYPE * BlockArnoldi<TYPE>::getiEvals() {
	int i;
	TYPE *temp_evals = new TYPE[ _nev ];
	for (i=0; i<_nev; i++) {
		temp_evals[i] = _evali[i];
	}
	return temp_evals;		
}

template <class TYPE>
TYPE * BlockArnoldi<TYPE>::getResiduals() {
	int i;
	TYPE *temp_resids = new TYPE[ _nev ];
	for (i=0; i<_nev; i++) {
		temp_resids[i] = _ritzresiduals[i];
	}
	return temp_resids;		
}

template <class TYPE>
void BlockArnoldi<TYPE>::setDebugLevel( const int level ) {
	_debuglevel = level;
}

template <class TYPE>
void BlockArnoldi<TYPE>::setSymmetric( const bool sym ) {
	_issym = sym;
}

template <class TYPE>
void BlockArnoldi<TYPE>::currentStatus() {
	int i;
	cout<<" "<<endl;
	cout<<"********************CURRENT STATUS********************"<<endl;
	cout<<"Iteration\t"<<_iter<<" of\t"<< _length+_restarts*(_length-_nevblock)<<endl;
	cout<<"Restart \t"<<_restartiter<<" of\t"<< _restarts<<endl;
	cout<<"Requested Eigenvalues : "<<_nev<<endl;
	cout<<"Requested Ordering : "<<_which<<endl;
	cout<<"Residual Tolerance : "<<_residual_tolerance<<endl;	
	cout<<"Error for the partial Schur decomposition is : "<< _schurerror <<endl;
	//
	//  Determine status of solver and output information correctly.
	//
	if ( _schurerror < _residual_tolerance ) {
		cout<<"------------------------------------------------------"<<endl;
		cout<<"Computed Eigenvalues: "<<endl;
	} else {
		if (exit_flg && _iter != _length+_restarts*(_length-_nevblock)) {
			cout<<"ERROR: Complete orthogonal basis could not be computed"<<endl;
		}
		cout<<"------------------------------------------------------"<<endl;
		cout<<"Current Eigenvalue Estimates: "<<endl;
	}
	//
	//  Print out current computed eigenvalues.  If we don't have all the requested
	//  eigenvalues yet, print out the ones we have.
	//
	int _nevtemp = _nev;
	if (_jstart < _nevblock) { _nevtemp = _jstart*_block; }
	//
	if (_issym) {
		cout<<"Eigenvalue \t Ritz Residual"<<endl;
		cout<<"------------------------------------------------------"<<endl;
		for (i=0; i<_nevtemp; i++) {
			cout.width(10);
			cout<<_evalr[i]<<"\t"<<_ritzresiduals[i]<<endl;
		}
		cout<<"------------------------------------------------------"<<endl;
		} else {
		cout<<"Real Part \t Imag Part \t Ritz Residual"<<endl;
		cout<<"------------------------------------------------------"<<endl;
		for (i=0; i<_nevtemp; i++) {
			cout.width(10);
			cout<<_evalr[i]<<"\t"<<_evali[i]<<"\t"<<_ritzresiduals[i]<<endl;
		}
		cout<<"------------------------------------------------------"<<endl;
		cout<<" "<<endl;
	}
	cout<<"******************************************************"<<endl;
}	

template <class TYPE>
void BlockArnoldi<TYPE>::SetInitBlock() {
	int i,j;
	int *index = new int[ _block ]; assert(index);

	// This method will set the first block of _basisvecs to the initial guess,
	// if one is given, else it will fill the block with random vectors.

	if (_initialguess) {
		const int cols = _ivec.GetNumberVecs();
		const int rows = _ivec.GetVecLength();
		if (cols < _block) {

			// Copy the given vectors in the first positions in the block
			// and fill the rest with random vectors.
			for (i=0; i<cols; i++) {
				index[i] = i;
			}
			_basisvecs->SetBlock( _ivec, index, cols );			

			// Initialize the rest of the block with random vectors
			for (i=cols; i<_block; i++) {
				index[i-cols] = i;
			}
			AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block-cols);
			assert(U_vec);
			U_vec->MvRandom();
			delete U_vec;
		}
		else {
			// Copy the first _block of the given vectors into the first _block
			// of _basisvecs, any additional vectors will be ignored.

			for (i=0; i<_block; i++) {
				index[i] = i;
			}
			_basisvecs->SetBlock( _ivec, index, _block );
		}
	}
	else {
		// No initial guess is given, so initialize block with random vectors

		for (i=0; i<_block; i++) {
         		index[i] = i;
		}
 		AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block);
		assert(U_vec);
		U_vec->MvRandom();
		delete U_vec;		
	}

	// Clean up
	delete [] index;
}

template <class TYPE>
void BlockArnoldi<TYPE>::iterate(const int steps) {
	int i,j,temp;
	int tempsteps = steps;
	const int izero=0;
	const TYPE one=1.0;
	const TYPE zero=0.0;
	//
	// If this is the first steps of Block Arnoldi, initialize the first block of _basisvecs
	//
	if (!_iter) {
		SetInitBlock();
		int *index = new int[ (_length+1)*_block ]; assert(index);
		for ( i=0; i<_block; i++ ) {
			index[i] = i;
		}
		AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block);
		assert(U_vec);
		AnasaziDenseMatrix<TYPE> G10(_block,_block);
		QRFactorization( *U_vec, G10 );
		delete U_vec;
		delete [] index;
	}				
	//
	// Leave the iteration method now if the orthogonal subspace can't be extended.
	//
	if (exit_flg) { return; }	
	//			
	// Now we go the number of steps requested by the user.  This may cause
	// a restart or hit the number of maximum iterations (restarts).  
	//
	while(tempsteps > 0 && _restartiter <= _restarts && !exit_flg) {
		// If we don't need to restart, just get it over with and return.
		if (_jstart+tempsteps < _length) {
			_jend = _jstart+tempsteps;
			_iter += tempsteps;
			tempsteps = 0;
			BlockReduction();
			if (exit_flg) { break; } // We need to leave before we move the pointer
			_jstart = _jend;  // Move the pointer
			ComputeResiduals( false );		
			_isdecompcurrent = false;
		}
		// Finish off this factorization and restart.
		else {  
			_jend = _length;
			temp = _length-_jstart;
			_iter += temp;
			tempsteps -= temp;
			BlockReduction();
			if (exit_flg) { break; } // We need to leave before we move the pointer
			_jstart = _length; // Move the pointer
			//
			//  Compute the Schur factorization and prepare for a restart.
			//
			ComputeResiduals( true );  
			Restart();  
			_isdecompcurrent = true;
			_restartiter++;
		}
	}
	//
	// Compute the current eigenvalue estimates before returning.
	//
	//cout<<"Upper Hessenberg matrix as of iteration :"<<_iter<<endl<<endl;
	//_hessmatrix->print();

	// Output current information if necessary
	if (_debuglevel > 0) {
		currentStatus();
	}
}

template <class TYPE>
void BlockArnoldi<TYPE>::solve () {
	int rem_iters = _length+_restarts*(_length-_nevblock)-_iter;
	//
	// Right now the solver will just go the remaining iterations, but this design will allow
	// for checking of the residuals every so many iterations, independent of restarts.
	//
	while (_schurerror > _residual_tolerance && _iter < _length+_restarts*(_length-_nevblock) && !exit_flg) {
		iterate( _step );
	}
}

template<class TYPE>
void BlockArnoldi<TYPE>::BlockReduction () {
	int i,j;
	
	int *index = new int[ _block ]; assert(index);
	
	for ( j = _jstart; j < _jend; j++ ) {
		//
		// Associate the j-th block of _basisvecs with U_vec.
		//
		for ( i=0; i<_block; i++ ) {
			index[i] = j*_block+i;
		}
		AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index, _block);
		assert(U_vec);
		//
		// Associate (j+1)-st block of ArnoldiVecs with F_vec.
		//
		for ( i=0; i<_block; i++ ) {
			index[i] = (j+1)*_block+i;
		}
		AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
		assert(F_vec);
		//
		//  Compute F_vec = OP * U_vec
		//
		_amat.ApplyMatrix( *U_vec, *F_vec ); 
		//
		// Use previous dependency information to decide which orthogonalization
		// method to use for the new block.  The global variable dep_flg tells us
		// if we've had problems with orthogonality before.  If no problems have
		// been detected before we will use standard block orthogonalization.
		// However, if there are problems with this, we will use a more stringent
		// approach.
		//
		if (!dep_flg) {
			BlkOrth(j);
		}
		//
		// If any block dependency was detected previously, then the more stringent
		// orthogonalization will be used.  If this method can't resolve the
		// dependency, then the exit_flg will be set indicating that we can't proceed
		// any further.
		//			
		if (dep_flg) {
			BlkOrthSing(j);
		}
		//
		delete U_vec, F_vec;
		//
		// If we cannot go any further with the factorization, then we need to exit
		// this method.
		//
		if (exit_flg) { return; }
	}
	delete [] index;
} // end BlockReduction()


template<class TYPE>
void BlockArnoldi<TYPE>::BlkOrth( const int j ) {
        //
        // Orthogonalization is first done between the new block of
        // vectors and all previous blocks, then the vectors within the
        // new block are orthogonalized.
        //
        const int IntOne = 1;
        const TYPE one = 1.0;
        const TYPE zero = 0.0;
        const int max_num_orth = 2;
        int i, k, row_offset, col_offset;
        int * index = new int[ (_length+1)*_block ]; assert(index);
        TYPE * norm1 = new TYPE[_block]; assert(norm1);
        TYPE * norm2 = new TYPE[_block]; assert(norm2);
        //
        // Associate (j+1)-st block of ArnoldiVecs with F_vec.
        //
        for ( i=0; i<_block; i++ ) {
                        index[i] = (j+1)*_block+i;
        }
        AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
        assert(F_vec);
        //
        // Zero out the full block column of the Hessenberg matrix
        // even though we're only going to set the coefficients in
        // rows [0:(j+1)*_block-1]
        //
        int ldh = _hessmatrix->getld();
        int n_row = _hessmatrix->getrows();
        int n_col = _hessmatrix->getcols();
        //
        TYPE* ptr_hess = _hessmatrix->getarray();
        for ( k=0; k<_block; k++ ) {
                for ( i=0; i<n_row ; i++ ) {
                                ptr_hess[j*ldh*_block + k*ldh + i] = zero;
                }
        }
        //
        // Grab all previous Arnoldi vectors
        //
        int num_prev = (j+1)*_block;
        for (i=0; i<num_prev; i++){
                index[i] = i;
        }
        AnasaziMultiVec<TYPE>* V_prev = _basisvecs->CloneView(index,num_prev);
        assert(V_prev);
        //
        // Create a matrix to store the product trans(V_prev)*F_vec
        //
        AnasaziDenseMatrix<TYPE> dense_mat(num_prev, _block );
        TYPE* ptr_dense = dense_mat.getarray();
        int ld_dense = dense_mat.getld();
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
                for ( k=0; k<_block; k++ ) {
                        for ( i=0; i<num_prev; i++ ) {
                                ptr_hess[j*ldh*_block + k*ldh + i] +=
                                        ptr_dense[k*ld_dense + i];
                        }
                }
                //
                // F_vec <- F_vec - V(0:(j+1)*block-1,:) * H(0:num_prev-1,j:num_prev-1)
                //
                F_vec->MvTimesMatAddMv( -one, *V_prev, dense_mat, one );
		//
        } // end for num_orth=0;...)
        //
        F_vec->MvNorm(norm2);
        //
        // Check to make sure the new block of Arnoldi vectors are
        // not dependent on previous Arnoldi vectors
        //
        for (i=0; i<_block; i++){
                if (norm2[i] < norm1[i] * _blk_tol) {
                        dep_flg = true;
                        if (_debuglevel > 2 ){
                           cout << "Col " << num_prev+i << " is dependent on previous "
                                    << "Arnoldi vectors in V_prev" << endl;
                           cout << endl;
                        }
                }
        } // end for (i=0;...)
        //
        if (_debuglevel>2) {
                CheckBlkArnRed(j);
        }
        //
        // If dependencies have not already been detected, compute
        // the QR factorization of the next block. Otherwise,
        // this block of Arnoldi vectors will be re-computed via and
        // implementation of A. Ruhe's block Arnoldi.
        //
        if (!dep_flg) {
                //
                // Compute the QR factorization of F_vec
                //
                row_offset = (j+1)*_block; col_offset = j*_block;
                AnasaziDenseMatrix<TYPE> sub_block_hess(*_hessmatrix, row_offset, col_offset,
                        _block, _block);
                QRFactorization( *F_vec, sub_block_hess );
        }
        //
        delete F_vec, V_prev;
        delete [] index;
        delete [] norm1;
        delete [] norm2;
        //
}  // end BlkOrth()


template<class TYPE>
void BlockArnoldi<TYPE>::BlkOrthSing( const int j ) {
        //
        // This is a variant of A. Ruhe's block Arnoldi
        // The orthogonalization of the vectors F_vec is done
        // one at a time. If a dependency is detected, a random
        // vector is added and orthogonalized against all previous
        // Arnoldi vectors.
        //
        const int IntOne = 1;
        const TYPE one = 1.0;
        const TYPE zero = 0.0;
        int i, k, iter, num_prev;
        int * index = new int[ (_length+1)*_block ]; assert(index);
        TYPE norm1[IntOne];
        TYPE norm2[IntOne];
        //
        // Associate (j+1)-st block of ArnoldiVecs with F_vec.
        //
        for ( i=0; i<_block; i++ ) {
        	index[i] = (j+1)*_block+i;
        }
        AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
        assert(F_vec);
        //
        // Zero out the full block column of the Hessenberg matrix
        //
        int ldh = _hessmatrix->getld();
        int n_row = _hessmatrix->getrows();
        int n_col = _hessmatrix->getcols();
        //
        TYPE* ptr_hess = _hessmatrix->getarray();
        for ( k=0; k<_block; k++ ) {
                for ( i=0; i<n_row ; i++ ) {
                        ptr_hess[j*ldh*_block + k*ldh + i] = zero;
                }
        }
        //
        AnasaziMultiVec<TYPE> *q_vec=0, *Q_vec=0, *tptr=0;
        tptr = F_vec->Clone(IntOne); assert(tptr);
        //
        // Start a loop to orthogonalize each of the _block
        // columns of F_vec against all previous _basisvecs
        //
        for (int iter=0; iter<_block; iter++){
                num_prev = (j+1)*_block + iter; // number of previous _basisvecs
                //
                // Grab the next column of _basisvecs
                //
                index[0] = num_prev;
                q_vec = _basisvecs->CloneView(index, IntOne); assert(q_vec);
                //
                // Grab all previous columns of _basisvecs
                //
                for (i=0; i<num_prev; i++){
                        index[i] = i;
                }
                Q_vec = _basisvecs->CloneView(index, num_prev); assert(Q_vec);
                //
                // Create matrix to store product trans(Q_vec)*q_vec
                //
                AnasaziDenseMatrix<TYPE> dense_mat(num_prev, IntOne);
                TYPE* ptr_dense = dense_mat.getarray();
                //
                // Do one step of classical Gram-Schmidt orthogonalization
                // with a 2nd correction step if needed.
                //
                q_vec->MvNorm(norm1);
                //
                // Compute trans(Q_vec)*q_vec
                //
                q_vec->MvTransMv(one, *Q_vec, dense_mat);
                //
                // Sum results [0:num_prev-1] into column (num_prev-_block)
                // of the Hessenberg matrix
                //
                for (k=0; k<num_prev; k++){
                        ptr_hess[(j*_block + iter)*ldh +k] += ptr_dense[k];
                }
		//
                // Compute q_vec<- q_vec - Q_vec * dense_mat
                //
                q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
                //
                q_vec->MvNorm(norm2);
                //
                if (norm2[0] < norm1[0] * _dep_tol) {
                        //
                        // Repeat process with newly computed q_vec
                        //
                    	// Compute trans(Q_vec)*q_vec
                    	//
                    	q_vec->MvTransMv(one, *Q_vec, dense_mat);
                    	//
                    	// Sum results [0:num_prev-1] into column (num_prev-_block)
                    	// of the Hessenberg matrix
                    	//
                    	for (k=0; k<num_prev; k++){
                            	ptr_hess[(j*_block + iter)*ldh +k] += ptr_dense[k];
                        }
			//
                    	// Compute q_vec<- q_vec - Q_vec * dense_mat
                    	//
                    	q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
                    	//
                    	q_vec->MvNorm(norm2);
                }
                //
                // Check for linear dependence
                //
                if (norm2[0] < norm1[0] * _sing_tol) {
                        if (_debuglevel > 2) {
                           cout << "Column " << num_prev << " of _basisvecs is dependent" 
				<< endl<<endl;
                        }
                        //
                        // Create a random vector and orthogonalize it against all
                        // previous cols of _basisvecs
                        // We could try adding a random unit vector instead -- not
                        // sure if this would make any difference.
                        //
                        tptr->MvRandom();
                        tptr->MvNorm(norm1);
                        //
                        // This code  is automatically doing 2 steps of orthogonalization
                        // after adding a random vector. We could do one step of
                        // orthogonalization with a correction step if needed.
                        //
                        for (int num_orth=0; num_orth<2; num_orth++){
                                tptr->MvTransMv(one, *Q_vec, dense_mat);
                                // Note that we don't change the entries of the
                                // Hessenberg matrix when we orthogonalize a
                                // random vector
                                tptr->MvTimesMatAddMv(-one, *Q_vec, dense_mat, one);
                        }
                        //
                        tptr->MvNorm(norm2);
                    	//
                        if (norm2[0] > norm1[0] * _sing_tol){
                                // Copy vector into the current column of _basisvecs
                		q_vec->MvAddMv( one, *tptr, zero, *tptr );
                                q_vec->MvNorm(norm2);
                		//
				// Normalize the new q_vec
                        	//
                        	TYPE rjj = one/norm2[0];
                        	q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
                        	//
                        	// Enter a zero in the [(j+1)*_block + iter] row in the
                        	// [(j*_block + iter] column of the Hessenberg matrix
                        	//
                		ptr_hess[(j*_block+iter)*ldh + (j+1)*_block+iter] = zero;
                        }
                        else {
                                // Can't produce a new orthonormal basis vector
                                // Clean up and exit this block Arnoldi factorization!
                                exit_flg = true;
				delete [] index;
        			delete q_vec; q_vec=0;
        			delete Q_vec; Q_vec=0;
        			delete tptr; tptr=0;
        			delete F_vec;
				return;
                        }
                }
                else {
                        //
                    	// Normalize the new q_vec
                    	//
                    	TYPE rjj = one/norm2[0];
                    	q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
                        //
                        // Enter norm of q_vec to the [(j+1)*_block + iter] row
                    	// in the [(j*_block + iter] column of the Hessenberg matrix
                        //
                        ptr_hess[(j*_block+iter)*ldh + (j+1)*_block+iter] = norm2[0];
                } // end else ...
        } // end for (i=0;...)
        //
        if (_debuglevel > 2){
                cout << "Checking Orthogonality after BlkOrthSing()"
                            << " Iteration: " << j << endl<<endl;
                CheckBlkArnRed(j);
        }
        //
        //      free heap space
        //
        delete [] index;
        delete q_vec; q_vec=0;
        delete Q_vec; Q_vec=0;
        delete tptr; tptr=0;
        delete F_vec;
} // end BlkOrthSing()

template<class TYPE>
void BlockArnoldi<TYPE>::QRFactorization (AnasaziMultiVec<TYPE>& VecIn, 
						AnasaziDenseMatrix<TYPE>& FouierR) {
	int i,j,k;
	int nb = VecIn.GetNumberVecs(); assert (nb == _block);
	int ldR = FouierR.getld();
	int *index = new int[nb]; assert(index);
	const int IntOne=1;
	const int IntZero=0;
	const TYPE zero=0.0;
	const TYPE one=1.0;
	bool addvec = false, flg = false;
	//
	TYPE * R = FouierR.getarray();
	TYPE norm1[IntOne];
	TYPE norm2[IntOne];
	AnasaziMultiVec<TYPE> *qj = 0, *Qj = 0, *tptr = 0;
	tptr = _basisvecs->Clone(IntOne); assert(tptr);
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
		qj = VecIn.CloneView(index, IntOne); assert(qj);
		//
		// If we are beyong the 1st column, orthogonalize against the previous
		// vectors in the current block.
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
			AnasaziDenseMatrix<TYPE> rj(j,1);
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
			if (_iter) {
        			// This is not the 1st block. A looser tolerance is used to
        			// determine dependencies. If a dependency is detected, a flag
        			// is set so we can back out this method and out of BlkOrth.
        			// The method BlkOrthSing is used to construct the new block
        			// of orthonormal basis vectors one at a time. If a dependency
        			// is detected within this method, a random vector is added
        			// and orthogonalized against all previous basis vectors.
        			//
        			if (norm2[0] < norm1[0] * _blk_tol) {
            				if (_debuglevel > 2) {
               					cout << "Column " << j << " of current block is dependent"<<endl;
	                		}
        	    			dep_flg = true;
            				delete qj; delete Qj; delete tptr;
                			delete [] index;
                			return;
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
                                 	AnasaziDenseMatrix<TYPE> tj(j,1);
                                 	//
                                 	tptr->MvRandom();
                                 	tptr->MvNorm(norm1);
                                     	//
                                 	for (int num_orth=0; num_orth<2; num_orth++){
                                         	tptr->MvTransMv(one, *Qj, tj);
                                         	tptr->MvTimesMatAddMv(-one, *Qj, tj, one);
                                     	}
                                 	tptr->MvNorm(norm2);
                                 	//
                                 	if (norm2[0] > norm1[0] * _sing_tol){
                                         	// Copy vector into current column of _basisvecs
                                         	qj->MvAddMv(one, *tptr, zero, *tptr);
                                     	}
					else {
                                         	exit_flg = true;
                                         	delete qj; delete Qj; delete tptr;
                                             	delete [] index;
                                             	return;
                                     	}
                             	}
			} // if (_iter) ...
		} // if (j) ...
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
      	} // for (j=0; j<nb; j++) ...
	//
	delete [] index;	
}

template<class TYPE>
void BlockArnoldi<TYPE>::ComputeResiduals( bool apply ) {
	int i=0,j=0;
	int m = _jstart*_block, n=_jstart*_block, info=0;
	int mm1 = (_jstart-1)*_block;
	const TYPE one = 1.0;
	const TYPE zero = 0.0;
	AnasaziDenseMatrix<TYPE>* Hj;	
	//
	// If we are going to restart then we can overwrite the 
	// hessenberg matrix with the Schur factorization, else we
	// will just use a copy of it.
	//
	if (apply) {
		Hj = new AnasaziDenseMatrix<TYPE>(*_hessmatrix, i, j, m, n);		
	} else {	
		// Create a view into the current hessenberg matrix and
		// make a copy.
		AnasaziDenseMatrix<TYPE> Hj_temp(*_hessmatrix, i, j, m, n);
		Hj = new AnasaziDenseMatrix<TYPE>(Hj_temp);
	}
	//
	//---------------------------------------------------
	// Compute the current eigenvalue estimates
	// ---> Use driver GEES to first reduce to upper Hessenberg 
	// 	form and then compute Schur form, outputting eigenvalues
	//---------------------------------------------------
	//
	AnasaziDenseMatrix<TYPE> Q(n,n);
	int sdim = 0; 
	int lwork = 4*n;
	int ldhj = Hj->getld();
	int ldq = Q.getld();
	int *select = new int[ n ];
	int *bwork = new int[ n ];
	char * jobvs = "V";
	char * sort = "N";
	TYPE *ptr_hj = Hj->getarray();
	TYPE *ptr_q = Q.getarray();
	TYPE *work = new TYPE[lwork]; assert(work);
	AnasaziLAPACK lapack;
	lapack.GEES( *jobvs, *sort, select, n, ptr_hj, ldhj, sdim,_evalr,
			_evali, ptr_q, ldq, work, lwork, bwork, &info );
	assert(info==0);
	//
	// Sort the eigenvalues, this also sorts the _order vector so we know
	// which ones we want. 
	//
//	cout<<"Before sorting (Hj):"<<endl;
//	Hj->print();
//	cout<<"Before sorting (Q):"<<endl;
//	Q.print();

	Sort( true );
	//
	// Reorder real Schur factorization, remember to add one to the indices for the
	// fortran call and determine offset.  The offset is necessary since the TREXC
	// method reorders in a nonsymmetric fashion, thus we use the reordering in
	// a stack-like fashion.
	//
	AnasaziBLAS blas;

	int _nevtemp, _nevtemp2;
	if (_jstart < _nevblock) {
		_nevtemp = n; _nevtemp2 = n;
	} else {
		_nevtemp = _nevblock*_block; _nevtemp2 = _nev;
	}		
	char * compq = "V";
	int *offset = new int[ _nevtemp ]; assert(offset);
	for (i=0; i<_nevtemp; i++) {
		offset[i] = 0;
		for (j=i; j<_nevtemp; j++) {
			if (_order[j] > _order[i]) { offset[i]++; }
		}
	}
	for (i=_nevtemp-1; i>=0; i--) {
		lapack.TREXC( *compq, n, ptr_hj, ldhj, ptr_q, ldq, _order[i]+1+offset[i], 
				1, work, &info );
		assert(info==0);
//		cout<<"Moving "<<_order[i]+1+offset[i]<<" to "<<1<<endl;

//		for (j=0; j<n; j++) {
//	  		cout<<j<<"\t"<<ptr_hj[j*ldhj + j]<<endl;
//		}
	}
//	cout<<"After sorting and reordering (Hj):"<<endl;
//	Hj->print();
//	cout<<"After sorting and reordering(Q):"<<endl;
//	Q.print();
	//
	// Check the residual error for the Krylov-Schur decomposition.
	// The residual for the Schur decomposition A(VQ) = (VQ)T + FB_m^TQ
	// where HQ = QT is || FB_m^TQ || <= || H_{m+1,m} || || B_m^TQ ||.
	//
	// We are only interested in the partial Krylov-Schur decomposition corresponding
	// to the _nev eigenvalues of interest or the _nevblock*_block number of
	// eigenvalues we're keeping.
	//
	//  Calculate the B matrix for the Krylov-Schur basis F_vec*B^T
	//
	AnasaziDenseMatrix<TYPE> sub_block_hess(*_hessmatrix, m, mm1, _block, _block);
	TYPE *ptr_sbh = sub_block_hess.getarray();
	int ld_sbh = sub_block_hess.getld();
	//
	AnasaziDenseMatrix<TYPE> sub_block_q( Q, mm1, 0, _block, _nevtemp );
	TYPE *ptr_sbq = sub_block_q.getarray();
	int ld_sbq = sub_block_q.getld();
	//
	AnasaziDenseMatrix<TYPE> sub_block_b( _block, _nevtemp );
	TYPE *ptr_sbb = sub_block_b.getarray();
	char* trans = "N";	
	blas.GEMM( *trans, *trans, _block, _nevtemp, _block, one, ptr_sbh, ld_sbh, 
		ptr_sbq, ld_sbq, zero, ptr_sbb, _block );
	AnasaziDenseMatrix<TYPE> sub_block_b2(sub_block_b, 0, 0, _block, _nevtemp2);
	AnasaziDenseMatrix<TYPE> sub_block_h(*_hessmatrix,0,0, _nevtemp, _nevtemp );
	//
	//  Compute approximate ritzresiduals for each eigenvalue
	//		
	TYPE _scalefactor = sub_block_h.getfronorm();
	_schurerror = sub_block_b2.getfronorm()/_scalefactor;
	//
	// ------------>  NOT SURE IF RITZRESIDUALS CAN BE UPDATED AFTER DEFLATION!
	//
        //for (i=0; i<_nevtemp ; i++) {
        for (i=_defblock*_block; i<_nevtemp ; i++) {
        	AnasaziDenseMatrix<TYPE> s(sub_block_b,0,i,_block,1);
                TYPE *ptr_s=s.getarray();
                _ritzresiduals[i] = blas.NRM2(_block, ptr_s)/_scalefactor;
        }   
	//
	//  We are going to restart, so update the Krylov-Schur decomposition.
	//
	if (apply) {	
		//
		// Update the Krylov basis.  Take into account that deflated blocks
		// need not be updated.
		//	
		int *index = new int[ n ]; assert(index);
		for (i = 0; i < n; i++ ) {
			index[i] = i;
		}
		AnasaziDenseMatrix<TYPE> Qnev(Q, n, _nevtemp);
		AnasaziMultiVec<TYPE>* basistemp = _basisvecs->CloneView( index, _nevtemp );
		AnasaziMultiVec<TYPE>* basistemp2 = _basisvecs->CloneCopy( index, n );
		basistemp->MvTimesMatAddMv ( one, *basistemp2, Qnev, zero );
		//
		// Update the Krylov-Schur form (quasi-triangular matrix).
		//
		AnasaziDenseMatrix<TYPE> Hjp1(*_hessmatrix,_nevtemp,0,_block,_nevtemp);
		TYPE* ptr_hjp1 = Hjp1.getarray();
		int ld_hjp1 = Hjp1.getld();
		for (i=0; i<_block; i++) {
		    for (j=0; j<_nevtemp; j++) {
			ptr_hjp1[j*ld_hjp1 + i] = ptr_sbb[j*_block + i];
		    }
		}
		delete basistemp, basistemp2;
		delete [] index;
	}			
	delete [] work; 
	delete [] bwork, select;
	delete [] offset;
}


template<class TYPE>
void BlockArnoldi<TYPE>::ComputeEvecs() {
	int i=0,j=0;
	int m = _jstart*_block, n=_jstart*_block, info=0;
	const TYPE one = 1.0;
	const TYPE zero = 0.0;
	//
	//  If the Krylov-Schur decomposition is not current, compute residuals
	//  like we are going to restart to update decomposition.
	//
	if (!_isdecompcurrent) { ComputeResiduals(true); }
	//
	//------------------------------------------------------------------------
	//  Now the eigenvalues and Krylov-Schur decomposition are current.
	//		A*V = V*T + F*B^T
	//  ----->  Compute the eigenvectors of T which were already sorted by
	//		ComputeResiduals.
	//------------------------------------------------------------------------
	//
	//  Create view into Schur matrix, then copy it so it doesn't get overwritten.
	//
	AnasaziDenseMatrix<TYPE> Hj_temp(*_hessmatrix, i, j, m, n);
	AnasaziDenseMatrix<TYPE> Hj(Hj_temp);
        //
	//  Now compute the eigenvectors of the Schur form
	//
        AnasaziDenseMatrix<TYPE> Q(n,n);
	char * side = "R";
	char * howmny = "A";
        int *select = new int[ n ];
        int lwork = 4*n;
        int ldhj = Hj.getld();
        int ldq = Q.getld();
        TYPE *ptr_hj = Hj.getarray();
        TYPE *ptr_q = Q.getarray();   
        TYPE *work = new TYPE[lwork]; assert(work); 
        AnasaziLAPACK lapack;
	int mm, ldvl = 1;
	TYPE *vl = new TYPE[ ldvl ];
	lapack.TREVC( *side, *howmny, select, n, ptr_hj, ldhj, vl, ldvl,
			ptr_q, ldq, n, &mm, work, &info );
	assert(info==0);
	//
	//  Convert back to approximate eigenvectors of the operator.
	//
	int * index = new int [ n ];
	for (i=0; i<n; i++) {
		index[i] = i;
	}
	AnasaziMultiVec<TYPE>* evecstemp = _basisvecs->Clone( n );
	AnasaziMultiVec<TYPE>* basistemp = _basisvecs->CloneView( index, n );
	evecstemp->MvTimesMatAddMv( one, *basistemp, Q, zero );
	//
	// Sort the eigenvectors, if symmetric problem don't worry about complex
	// eigenvectors.
	//
	if (_issym) {
		_evecr->SetBlock( *evecstemp, index, _nev );
	} else {  // Right now only the real part of the eigenvector is being set!
		AnasaziMultiVec<TYPE>* evecstempi = _basisvecs->Clone( n );
		evecstempi->MvAddMv( -one, *evecstemp, zero, *evecstemp );
		int conjprs=0;
		int * _orderr = new int [ _nev ];
		int * _orderi = new int [ _nev ];
		i = 0;
		while ( i<_nev ) {	
			if (_evali[i] != zero) {
				_orderr[i] = i;
				_orderr[i+1] = i;										
				i = i+2; conjprs++;
			} else {
				_orderr[i] = i;
				i++;
			}
		}
		cout<< "There are "<< conjprs <<" conjugate pairs of eigenvalues"
			<<endl;
		AnasaziMultiVec<TYPE>* evecstemp2 = evecstemp->CloneView( _orderr, _nev );
		_evecr->SetBlock( *evecstemp2, index, _nev );
		delete evecstempi, evecstemp2;
	}

	_isevecscurrent = true;
	delete [] index;
	delete evecstemp, basistemp;
}

template<class TYPE>
void BlockArnoldi<TYPE>::Restart() {
	//  This method assumes the ComputeResiduals has been called before it
	//  to compute the Schur vectors and residuals.  This information is used to 
	//  restart the factorization.
	//
	int i,j, defcnt;
	int _nevtemp = _nevblock*_block;
	int *index = new int[ _nevtemp ];
	//
	//  Move the F_vec block to the _jstart+1 position.	
	//
	for (i=0; i<_block; i++) {
		index[i] = _jstart*_block + i;
	}
	AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneCopy(index, _block);
	for (i=0; i<_block; i++) {
		index[i] = _nevtemp + i;
	}
	_basisvecs->SetBlock( *F_vec, index, _block);
	//
	//  Check for blocks to deflate
	//
	i = _defblock;
	while ( i<_nevblock ) {
		defcnt = 0;
		for (j=0; j<_block; j++) {
			if (_ritzresiduals[i*_block+j] < _def_tol ) { defcnt++; }
		}
		if (defcnt == _block) {
			_defblock++;
		}
		i++;
	}		
	//
	//  If there are blocks to deflate, we need to set the subdiagonal entries to zero
	//
	if (_defblock > 0) {
		if (_debuglevel > 2) {
			cout<<"Deflating blocks with eigenvalue residuals below : "<<_def_tol<<endl;
			cout<<"Number of blocks being deflated : "<<_defblock<<endl;
		}
		TYPE zero = 0.0;
		AnasaziDenseMatrix<TYPE> Hj_temp(*_hessmatrix, _nevtemp, 0, _block, _defblock*_block);
		TYPE *ptr_hj = Hj_temp.getarray();
		int ld_hj = Hj_temp.getld();
		for (i=0; i<_block; i++) {
		    for (j=0; j<_defblock*_block; j++) {
			ptr_hj[j*ld_hj + i] = zero;
		    }
		}
	}
	//
	//  Reset the pointer.
	//
	_jstart = _nevblock; 
	//
	//  Clean up
	//
	delete F_vec;
	delete [] index;
}

template<class TYPE>
void BlockArnoldi<TYPE>::CheckBlkArnRed( const int j ) {
        int i,k,m=(j+1)*_block;
        int *index = new int[m];
                
        for ( i=0; i<_block; i++ ) {
                index[i] = m+i;
        }
        AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
        assert(F_vec);
                        
        TYPE *ptr_norms = new TYPE[m];
        TYPE sum=0.0;
                        
        F_vec->MvNorm(ptr_norms);
        for ( i=0; i<_block; i++ ) {
                sum += ptr_norms[i];
        }
                
        for ( i=0; i<m; i++ ) {
                index[i] = i;  
        }
        AnasaziMultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
        assert(Vj);   
        cout << " " << endl;
        cout << "********Block Arnoldi iteration******** " << j << endl;
        cout << " " << endl;
                        
        const TYPE one=1.0;
        const TYPE zero=0.0;
        AnasaziDenseMatrix<TYPE> VTV(m,m);
        Vj->MvTransMv(one,*Vj,VTV);
        TYPE* ptr=VTV.getarray();
        TYPE column_sum;
        
        for (k=0; k<m; k++) {
                column_sum=zero;
                for (i=0; i<m; i++) {
                        if (i==k) {
                                ptr[i] -= one;
                        }
                        column_sum += ptr[i];
                }
                cout <<  " V^T*V-I " << "for column " << k << " is " << fabs(column_sum) << endl;
                ptr += m;
        }
        cout << " " << endl;
        
        AnasaziDenseMatrix<TYPE> E(m,_block);
        
        F_vec->MvTransMv(one,*Vj,E);
        TYPE* ptr_Ej=E.getarray();
                        
        for (k=0;k<_block;k++) {
                column_sum=zero;
                for (i=0; i<m; i++) {
                        column_sum += ptr_Ej[i];
                }
                ptr_Ej += m;
                if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
                cout << " Orthogonality with F " << "for column " << k << " is " << fabs(column_sum) << endl;
}
        cout << " " << endl;
                 
        AnasaziMultiVec<TYPE>* AVj = _basisvecs->Clone(m); assert(AVj);
        _amat.ApplyMatrix(*Vj,*AVj);
                        
        int row_offset=0;
        int col_offset=0;
        AnasaziDenseMatrix<TYPE> Hj(*_hessmatrix, row_offset, col_offset, m, m);
        AVj->MvTimesMatAddMv(-one, *Vj, Hj, one);
        for ( i=0; i<_block; i++ ) {  
                index[i] = j*_block+i;
        }

        AnasaziMultiVec<TYPE>* Fj = AVj->CloneView(index, _block);
        Fj->MvAddMv(-one, *F_vec, one, *Fj);
	
        AVj->MvNorm(ptr_norms);
        
        for ( i=0; i<m; i++ ) { 
                cout << " Arnoldi relation " << "for column " << i << " is " << fabs(ptr_norms[i]) << endl;        
	}
        cout << " " << endl;
                
        delete F_vec;
        delete Fj;
        delete AVj;
        delete Vj;
        delete [] index;
        delete [] ptr_norms;
}


template<class TYPE>
void BlockArnoldi<TYPE>::Sort( const bool apply ) {
	int i, j, tempord;
	const int n = _jstart*_block;
	TYPE temp, tempr, tempi;
	AnasaziLAPACK lapack;
	//
	// Reset the index
	//		
	for (i=0; i < n; i++) {
		_order[i] = i;
	}
	//
	// These methods use an insertion sort method to circument recursive calls.
	//---------------------------------------------------------------
	// Sort eigenvalues in increasing order of magnitude
	//---------------------------------------------------------------
	if (!_which.compare("SM")) {
		if (_issym) {  // The eigenvalues are real
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; 
				tempord = _order[j];
				temp = abs(_evalr[j]);
				for (i=j-1; i>=0 && abs(_evalr[i])>temp; --i) {
					_evalr[i+1]=_evalr[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _order[i+1] = tempord;	
			}
		}
		else {  // The eigenvalues may be complex
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; tempi = _evali[j]; 
				tempord = _order[j];
				temp=lapack.LAPY2(_evalr[j],_evali[j]);
				for (i=j-1; i>=0 && lapack.LAPY2(_evalr[i],_evali[i])>temp; --i) {
					_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
			}	
		}
	}
	//---------------------------------------------------------------
	// Sort eigenvalues in increasing order of real part
	//---------------------------------------------------------------
	if (!_which.compare("SR")) {
		if (_issym) {  // The eigenvalues are real
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; 
				tempord = _order[j];
					for (i=j-1; i>=0 && _evalr[i]>tempr; --i) {
						_evalr[i+1]=_evalr[i];
						_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _order[i+1] = tempord;	
			}
		}
		else {  // The eigenvalues may be complex
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; tempi = _evali[j]; 
				tempord = _order[j];
				for (i=j-1; i>=0 && _evalr[i]>tempr; --i) {
					_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
			}	
		}
	}
	//---------------------------------------------------------------
	// Sort eigenvalues in increasing order of imaginary part
	//---------------------------------------------------------------
	if (!_which.compare("SI")) {
		for (j=1; j < n; ++j) {
			tempr = _evalr[j]; tempi = _evali[j]; 
			tempord = _order[j];
			for (i=j-1; i>=0 && _evali[i]>tempi; --i) {
				_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
				_order[i+1]=_order[i];
			}
			_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
		}
	}
	//---------------------------------------------------------------
	// Sort eigenvalues in decreasing order of magnitude
	//---------------------------------------------------------------
	if (!_which.compare("LM")) {
		if (_issym) {  // The eigenvalues are real
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; 
				tempord = _order[j];
				temp = abs(_evalr[j]);
				for (i=j-1; i>=0 && abs(_evalr[i])<temp; --i) {
					_evalr[i+1]=_evalr[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _order[i+1] = tempord;	
			}
		}
		else {  // The eigenvalues may be complex
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; tempi = _evali[j]; 
				tempord = _order[j];
				temp=lapack.LAPY2(_evalr[j],_evali[j]);
				for (i=j-1; i>=0 && lapack.LAPY2(_evalr[i],_evali[i])<temp; --i) {
					_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
			}	
		}
	}
	//---------------------------------------------------------------
	// Sort eigenvalues in decreasing order of real part
	//---------------------------------------------------------------
	if (!_which.compare("LR")) {
		if (_issym) {  // The eigenvalues are real
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; 
				tempord = _order[j];
					for (i=j-1; i>=0 && _evalr[i]<tempr; --i) {
						_evalr[i+1]=_evalr[i];
						_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _order[i+1] = tempord;	
			}
		}
		else {  // The eigenvalues may be complex
			for (j=1; j < n; ++j) {
				tempr = _evalr[j]; tempi = _evali[j]; 
				tempord = _order[j];
				for (i=j-1; i>=0 && _evalr[i]<tempr; --i) {
					_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
					_order[i+1]=_order[i];
				}
				_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
			}	
		}
	}
	//---------------------------------------------------------------
	// Sort eigenvalues in decreasing order of imaginary part
	//---------------------------------------------------------------
	if (!_which.compare("LI")) {
		for (j=1; j < n; ++j) {
			tempr = _evalr[j]; tempi = _evali[j]; 
			tempord = _order[j];
			for (i=j-1; i>=0 && _evali[i]<tempi; --i) {
				_evalr[i+1]=_evalr[i]; _evali[i+1]=_evali[i];
				_order[i+1]=_order[i];
			}
			_evalr[i+1] = tempr; _evali[i+1] = tempi; _order[i+1] = tempord;	
		}
	}
}

} // End of namespace Anasazi
#endif
// End of file BlockArnoldi.hpp

