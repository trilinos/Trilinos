// File BlockArnoldi.hpp
#ifndef BLOCK_ARNOLDI_HPP
#define BLOCK_ARNOLDI_HPP

#include "AnasaziLAPACK.hpp"
#include "AnasaziBLAS.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziCommon.hpp"

//using namespace std;

// 
// BlockArnoldi base class
//
template <class TYPE>
class BlockArnoldi { 
public:
	BlockArnoldi( AnasaziMatrix<TYPE> & mat, AnasaziMultiVec<TYPE>& vec,
		const TYPE tol=1.0e-6, const int nev=5, const int length=25, 
		const int block=1, const string which="LM", const int restarts=0 );
	virtual ~BlockArnoldi();
	void iterate( int steps=1 );
	void solve();
	AnasaziMultiVec<TYPE>* getEvecs(); 
	AnasaziMultiVec<TYPE>* getiEvecs();
	TYPE * getEvals();
	TYPE * getiEvals();
	TYPE * getResiduals();
	void setDebugLevel( const int );
	void setSymmetric( const bool );
	void currentStatus();
private:
	void QRFactorization( AnasaziMultiVec<TYPE>&, AnasaziDenseMatrix<TYPE>& );
	void BlockReduction();
	void ComputeResiduals( const bool );
	void ComputeEvecs();
	void Restart();
	void Sort( const bool );
	void SetInitBlock();
	void Check_Block_Arn_Red( const int );
	AnasaziMatrix<TYPE> &_amat; // must be passed in by the user
	AnasaziMultiVec<TYPE> &_ivec; // must be passed in by the user
	AnasaziMultiVec<TYPE> *_basisvecs, *_evecr, *_eveci;
	AnasaziDenseMatrix<TYPE>* _hessmatrix;
	AnasaziDenseMatrix<TYPE>* _B; // needed for block Krylov decomposition
	const int _nev, _length, _block, _restarts;
	const TYPE _residual_tolerance;
	TYPE *_ritzresiduals, *_evalr, *_evali;
	int _restartiter, _iter, _jstart, _jend, _nevblock, _debuglevel, _nconv;
	bool _initialguess, _qrfact, _issym;
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
				const string which, const int restarts ) : 
				_amat(mat), _ivec(vec), _basisvecs(0), _evecr(0), _eveci(0), 
				_hessmatrix(0), _nev(nev), _length(length), _block(block), 
				_restarts(restarts), _residual_tolerance(tol), _ritzresiduals(0), 
				_evalr(0), _evali(0), _restartiter(0), _B(0), 
				_iter(0), _jstart(0), _jend(0), _which(which),
				_initialguess(false), _debuglevel(0), _nevblock(0),
				_qrfact(false), _issym(false), _nconv(0) {
	//	std::cout << "ctor:BlockArnoldi " << this << std::endl;
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
		// Create the vectors for eigenvalues and their residual errors
		//
		_evalr = new TYPE[ _block*_length ]; assert(_evalr);  
		_evali = new TYPE[ _block*_length ]; assert(_evali);  
		_ritzresiduals = new TYPE[ _block*_length ]; assert(_ritzresiduals);
		_order = new int[ _block*_length ]; assert(_order);
		//
		// Create the block matrix B need for Krylov decomposition, initialize it to I_b. 
		//
		_B = new AnasaziDenseMatrix<TYPE>(_block, _block); assert(_B);
		int i, _block2 = _block*_block;
		TYPE * array = new TYPE [ _block2 ];
		TYPE one = 1.0, zero = 0.0;
		for (i=0; i<_block2; i++) {
			array[i] = zero;
		}
		for (i=0; i<_block; i++) {
			array[i*_block + i] = one;				
		}
		_B->setvalues( array, _block );
		delete [] array;
	}
	else {
		std::cout << "BlockArnoldi:ctor " << _length << _block << _basisvecs << std::endl;
		exit(-1);
	}
}

template <class TYPE>
BlockArnoldi<TYPE>::~BlockArnoldi() {
	//	std::cout << "dtor:BlockArnoldi " << this << std::endl;
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
AnasaziMultiVec<TYPE> * BlockArnoldi<TYPE>::getEvecs() {
	return _evecr->CloneCopy();
}
 
template <class TYPE>
AnasaziMultiVec<TYPE> * BlockArnoldi<TYPE>::getiEvecs() {
	return _eveci->CloneCopy();
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
	std::cout<<" "<<std::endl;
	std::cout<<"********************CURRENT STATUS********************"<<std::endl;
	std::cout<<"Iteration\t"<<_iter<<" of\t"<< _length*(_restarts+1)<<std::endl;
	std::cout<<"Restart \t"<<_restartiter<<" of\t"<< _restarts<<std::endl;
	std::cout<<"Converged Eigenvalues : "<<_nconv<<std::endl;
	std::cout<<"Requested Eigenvalues : "<<_nev<<std::endl;
	std::cout<<"Requested Ordering : "<<_which<<std::endl;
	std::cout<<"Residual Tolerance : "<<_residual_tolerance<<std::endl;	
	std::cout<<"------------------------------------------------------"<<std::endl;
	std::cout<<"Current Eigenvalue Estimates: "<<std::endl;
	if (_issym) {
		std::cout<<"Eigenvalue \t Residual"<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
		for (i=0; i<_nev; i++) {
			std::cout.width(10);
			std::cout<<_evalr[i]<<"\t"<<_ritzresiduals[i]<<std::endl;
		}
		std::cout<<"------------------------------------------------------"<<std::endl;
	} else {
		std::cout<<"Real Part \t Imag Part \t Residual"<<std::endl;
		std::cout<<"------------------------------------------------------"<<std::endl;
		for (i=0; i<_nev; i++) {
			std::cout.width(10);
			std::cout<<_evalr[i]<<"\t"<<_evali[i]<<"\t"<<_ritzresiduals[i]<<std::endl;
		}
		std::cout<<"------------------------------------------------------"<<std::endl;
	}
	std::cout<<" "<<std::endl;
}	

template <class TYPE>
void BlockArnoldi<TYPE>::SetInitBlock() {
	int i,j;
	int *index = new int[ _block ]; assert(index);

	// This routine will set the first block of _basisvecs to the initial guess,
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
	delete [] index; index=0;
}


template <class TYPE>
void BlockArnoldi<TYPE>::iterate(int steps) {
	int i,j,loops,temp;
	const int izero=0;
	const TYPE one=1.0;
	const TYPE zero=0.0;
	int *index = new int[ (_length+1)*_block ]; assert(index);

	// If this is the first steps of Block Arnoldi, initialize the first block of _basisvecs

	if (!_iter) {
		SetInitBlock();
		for ( i=0; i<_block; i++ ) {
			index[i] = i;
		}
		AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block);
		assert(U_vec);
		AnasaziDenseMatrix<TYPE> G10(_block,_block);
		QRFactorization( *U_vec, G10 );
		_qrfact = true;
		delete U_vec;
	}				
				
	// Now we go the number of steps requested by the user.  This may cause
	// a restart or hit the number of maximum iterations (restarts).  Currently this
	// doesn't look at the residuals.
	//
	// Lets see how many loops we'll have to do, take into account that we may not
	// be in the beginning of the factorization.
	//

	if (_jstart+steps < _length) {
		// If we don't need to restart, just get it over with.
		_jend = _jstart+steps;
		_iter += steps;
		BlockReduction();
		_jstart = _jend;  // Move the pointer
	}
	else {  
		// Finish off this factorization and restart;
		_jend = _length;
		temp = _length-_jstart;
		_iter += temp;
		steps -= temp;
		BlockReduction();
		_jstart = _length; // Move the pointer
		//
		// Now restart and do the necessary loops required by the number of steps given.
		//
		while (steps > 0 && _restartiter <  _restarts) {
			ComputeResiduals( true );
			Restart();  // Assuming _jstart is set in Restart()

			// Output current information if necessary
			if (_debuglevel > 0) {
				currentStatus();
			}
			_restartiter++;

			//  Reset _jend if we don't have many steps left
			if (steps < (_length-_jstart)) {
				_jend = _jstart+steps;
			}
			BlockReduction();
			temp = _jend-_jstart;
			_iter += temp;
			steps -= temp;
			_jstart = _jend; // Move the pointer
		}
	}
	//
	// Compute the current eigenvalue estimates before returning.
	//
	ComputeResiduals( false );		

	// Output current information if necessary
	if (_debuglevel > 0) {
		currentStatus();
	}
}

template <class TYPE>
void BlockArnoldi<TYPE>::solve () {
	int rem_iters = _length*(_restarts+1)-_iter;

	// Right now the solver will just go the remaining iterations, but this design will allow
	// for checking of the residuals every so many iterations, independent of restarts.
	iterate( rem_iters );
}

template<class TYPE>
void BlockArnoldi<TYPE>::BlockReduction () {
	int i,j,k,row_offset,col_offset;
	int *index=0;
	const int max_num_orth=2;
	const int IntOne=1;
	const int IntZero=0;
	const TYPE zero=0.0;
	const TYPE one=1.0;
	
	index = new int[ _block*_length ]; assert(index);
	
	int ldh = _hessmatrix->getld();
	int n_row = _hessmatrix->getrows();
	int n_col = _hessmatrix->getcols();

	for ( j = _jstart; j < _jend; j++ ) {
		//
		// Associate the j-th block of _basisvecs with U_vec.
		//
		for ( i=0; i<_block; i++ ) {
			index[i] = j*_block+i;
		}
		AnasaziMultiVec<TYPE>* U_vec = _basisvecs->CloneView(index, _block);
		assert(U_vec);
		if (!_qrfact) {
			//
			// Compute the QR factorization of U_vec
			//
			if (j==0) {
				row_offset = 0; col_offset = 0;
			}
			else {
				row_offset = j*_block; col_offset = (j-1)*_block;
			}
			AnasaziDenseMatrix<TYPE> sub_block_hess(*_hessmatrix, row_offset, col_offset,
				_block, _block);
			QRFactorization( *U_vec, sub_block_hess );
		}
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
		//  Compute trans(V_j)*F_vec and store in the j-th diagonal
		//  block of Hess_matrix.
		//
		AnasaziDenseMatrix<TYPE> dense_mat((j+1)*_block, _block );
		for ( i=0; i<(j+1)*_block; i++ ) {
			index[i] = i;
		}
		AnasaziMultiVec<TYPE>* Vj = _basisvecs->CloneView(index, (j+1)*_block);
		assert(Vj);
		//
		// Zero out the full block column of Hess_matrix even though we're only
		// going to set the coefficients in rows 0:(j+1)*block-1
		//
		TYPE* ptr_hess = _hessmatrix->getarray();
		for ( k=0; k<_block; k++ ) {
			for ( i=0; i<n_row ; i++ ) {
				ptr_hess[j*ldh*_block + k*ldh + i] = zero;
			}
		}
		//
		// Perform two steps of block classical Gram-Schmidt so that
		// F_vec is orthogonal to the columns of Vj.
		//
		TYPE* ptr_dense = dense_mat.getarray();
		int ld_dense = dense_mat.getld();
		int num_orth;
		for ( num_orth=0; num_orth<max_num_orth; num_orth++ ) {
			F_vec->MvTransMv (one, *Vj, dense_mat);
			//
			// Update the orthogonalization coefficients for the j-th block
			// column of Hess_matrix.
			//
			for ( k=0; k<_block; k++ ) {
				for ( i=0; i<(j+1)*_block; i++ ) {
					ptr_hess[j*ldh*_block + k*ldh + i] += 
						ptr_dense[k*ld_dense + i];
				}
			}
			//
			// F_vec <- F_vec - V(0:(j+1)*block-1,:) * H(0:(j+1)*block-1,j:(j+1)*block-1)
			//
			F_vec->MvTimesMatAddMv( -one, *Vj, dense_mat, one );
		}
		if (_debuglevel>2) {
                        Check_Block_Arn_Red(j);  
                }
		//
		//	Compute the QR factorization of the next block
		//
		if (_qrfact) {
			//
			// Compute the QR factorization of U_vec
			//
			row_offset = (j+1)*_block; col_offset = j*_block;
			AnasaziDenseMatrix<TYPE> sub_block_hess(*_hessmatrix, row_offset, col_offset,
				_block, _block);
			QRFactorization( *F_vec, sub_block_hess );
		}
		//
		//	free heap for located allocated memory
		//
		delete U_vec;
		delete F_vec;
		delete Vj;
	}
	delete [] index;
}

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
	const TYPE eps = 1.e-16; // replace this with call to LAPACK lamch
	TYPE * R = FouierR.getarray();
	TYPE * NormVecIn = new TYPE[nb]; assert(NormVecIn);
	AnasaziMultiVec<TYPE> *qj = 0, *Qj = 0;
	//
    // Zero out the array that will contain the Fourier coefficients.
	//
	for ( j=0; j<nb; j++ ) {
		for ( i=0; i<nb; i++ ) {
			R[j*ldR+i] = zero;
		}
	}
	//
	//  Get the norms of the columns of VecIn; this will be used to
	//  determine whether there is rank deficiency in VecIn.
	//
	VecIn.MvNorm(NormVecIn);
	//
	// Compute the Forbenuis norm of VecIn
	//
	TYPE FroNorm=0.0;
	for ( j=0; j<nb; j++ ) {
		FroNorm += NormVecIn[j]*NormVecIn[j];
	}
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
		if ( *(R+j*ldR+j) > eps*FroNorm ) {
			//
			// Normalize qj to make it into a unit vector.
			//
			TYPE rjj = one / *(R+j*ldR+j);
			qj->MvAddMv ( rjj, *qj, zero, *qj );
		}
		else {
			// 
			// If qj is the zero vector, get ready to exit this function.
			//
			std::cout << " ***rank deficient block*** " << j << std::endl;
		}	
		delete qj;
		delete Qj;
	}
	delete [] index;
	delete [] NormVecIn;
}

template<class TYPE>
void BlockArnoldi<TYPE>::ComputeResiduals( bool apply ) {
	int i=0,j=0;
	int m = _jstart*_block, n=_jstart*_block, info=0;
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
	int sdim; 
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
	Sort( true );
	//
	// Reorder real Schur factorization, remember to add one to the indices for the
	// fortran call and determine offset.  The offset is necessary since the TREXC
	// routine reorders in a nonsymmetric fashion, thus we use the reordering in
	// a stack-like fashion.  Only necessary if we are restarting.
	//

	int _nevtemp;
	if (_jstart < _nevblock) {
		_nevtemp = n;
	} else {
		_nevtemp = _nevblock*_block;
	}		
	int *index = new int[ n ]; assert(index);
	for (i = 0; i < n; i++ ) {
		index[i] = i;
	}
	AnasaziBLAS blas;

	if (apply) {
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
//			std::cout<<"Moving "<<_order[i]+1+offset[i]<<" to "<<1<<std::endl;
//			std::cout<<"After sorting and reordering"<<std::endl;
//			for (j=0; j<n; j++) {
//	  			std::cout<<j<<"\t"<<ptr_hj[j*ldhj + j]<<std::endl;
//			}
		}

		AnasaziDenseMatrix<TYPE> Qnev(Q,n,_nevtemp);
		AnasaziMultiVec<TYPE>* basistemp = _basisvecs->CloneView( index, _nevtemp );
		AnasaziMultiVec<TYPE>* basistemp2 = _basisvecs->CloneCopy( index, n );
		basistemp->MvTimesMatAddMv ( one, *basistemp2, Qnev, zero );
		//
		// Check the residual errors
		//
	//==========================================================
	//	need a matrix matrix multiply here, but not sure what form B should be in !
	//	AnasaziDenseMatrix<TYPE> s(Q, (_jstart-1)*_block, 0, _block, n)
	//==========================================================
		for (i=0; i<n ; i++) {
 	      		AnasaziDenseMatrix<TYPE> s(Q,(_jstart-1)*_block,i,_block,1);
        		TYPE *ptr_s=s.getarray();
        		_ritzresiduals[i] = blas.NRM2(_block, ptr_s);
		}
		delete basistemp, basistemp2;
		delete [] offset;
	}
	else {
		//
		// Check the residual errors, the Schur form was never reordered, so use
		// the _order vector from the sorting routine.
		//
	//==========================================================
	//	need a matrix matrix multiply here, but not sure what form B should be in !
	//	AnasaziDenseMatrix<TYPE> s(Q, (_jstart-1)*_block, 0, _block, n)
	//==========================================================
		for (i=0; i<n ; i++) {
 	      		AnasaziDenseMatrix<TYPE> s(Q,(_jstart-1)*_block,_order[i],_block,1);
        		TYPE *ptr_s=s.getarray();
        		_ritzresiduals[i] = blas.NRM2(_block, ptr_s);
		}
	}		
	//---------------------------------------------------------------------
	//  Check convergence before returning to iterate/solve routine
	//---------------------------------------------------------------------
	_nconv = 0;
	std::cout<<"Checking residuals for tolerance : "<<_residual_tolerance<<std::endl;
	for (i=0; i<_nev; i++) {
		std::cout<<"Eigenvalue "<<i<<" : "<<_ritzresiduals[i]<<std::endl;
		if ( _ritzresiduals[i] < _residual_tolerance ) {
			_nconv++;		
		}
	}			
	std::cout<<"Converged eigenvalues : "<<_nconv<<std::endl<<std::endl;

	delete [] work; 
	delete [] index;
}


template<class TYPE>
void BlockArnoldi<TYPE>::ComputeEvecs() {
	int i=0,j=0;
	int m = _jstart*_block, n=_jstart*_block, info=0;
	const TYPE one = 1.0;
	const TYPE zero = 0.0;
	//
	//  Create view into upper Hessenberg matrix, then copy it so it
	//  doesn't get overwritten.
	//
	AnasaziDenseMatrix<TYPE> Hj_temp(*_hessenberg, i, j, m, n);
	AnasaziDenseMatrix<TYPE> Hj(Hj_temp);
        //      
        //---------------------------------------------------
        // Compute the current eigenvalue estimates
        // ---> Use driver GEES to first reduce to upper Hessenberg
        //      form and then compute Schur form, outputting eigenvalues
        //---------------------------------------------------
        //
        AnasaziDenseMatrix<TYPE> Q(n,n);
        int sdim;
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
        Sort( true );
	//
	// Now compute the eigenvectors of the Schur form
	//
	char * side = "R";
	char * howmny = "A";
	int mm, ldvl = 1;
	TYPE *vl = new TYPE[ ldvl ];
	lapack.TREVC( *side, *howmny, select, n, ptr_hj, ldhj, vl, ldvl,
			ptr_q, ldq, n, mm, work, &info );
	assert(info==0);
	//
	// Convert back to approximate eigenvectors of the matrix A, then sort
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
		AnasaziMultiVec<TYPE>* evecstemp2 = evecstemp->CloneView( _order, _nev );
		_evecr->SetBlock( *evecstemp2, index, _nev );
		delete evecstemp2;
	} else {  // Right now only the real part of the eigenvector is being set!
		AnasaziMultiVec<TYPE>* evecstempi = _basisvecs->Clone( n );
		evecstempi->MvAddMv( -one, evecstemp, zero, evecstemp );
		int conjprs=0;
		int * _orderr = new int [ _nev ];
		int * indexi = new int [ _nev ];
		i = 0;
		while ( i<_nev ) {	
			ord_i = _order[i];
			if (_evali[i] != zero) {
				_orderr[i] = ord_i;
				_orderr[i+1] = ord_i;										
				i = i+2;
			} else {
				_orderr[i] = ord_i;
				i++;
			}
		}
		AnasaziMultiVec<TYPE>* evecstemp2 = evecstemp->CloneView( _orderr, _nev );
		_evecr->SetBlock( *evecstemp2, index, _nev );
		delete evecstempi, evecstemp2;
	}

	delete [] index;
	delete evecstemp, basistemp;
}

template<class TYPE>
void BlockArnoldi<TYPE>::Restart() {
	//  This routine assumes the ComputeResiduals has been called before it
	//  to compute the Schur vectors and residuals.  This information is used to 
	//  restart the factorization.

	int i;
	int *index = new int[ _block ];
	for (i=0; i<_block; i++) {
		index[i] = _jstart*_block + i;
	}
	// Move the F_vec block to the _jstart+1 position before resetting the pointer.	

	AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneCopy(index, _block);
	for (i=0; i<_block; i++) {
		index[i] = _nevblock*_block + i;
	}
	_basisvecs->SetBlock( *F_vec, index, _block);
	_jstart = _nevblock; // or is it +1?
}

template<class TYPE>
void BlockArnoldi<TYPE>::Check_Block_Arn_Red( const int j ) {
        int i,k,m=(j+1)*_block;
        int *index = new int[m];
                
        for ( i=0; i<_block; i++ ) {
                index[i] = m+i;
        }
        AnasaziMultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
        assert(F_vec);
                        
        TYPE *ptr_norms = new double[m];
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
        std::cout << " " << std::endl;
        std::cout << "********Block Arnoldi iteration******** " << j << std::endl;
        std::cout << " " << std::endl;
                        
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
                std::cout <<  " V^T*V-I " << "for column " << k << " is " << fabs(column_sum) << std::endl;
                ptr += m;
        }
        std::cout << " " << std::endl;
        
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
                std::cout << " Orthogonality with F " << "for column " << k << " is " << fabs(column_sum) << std::endl;
}
        std::cout << " " << std::endl;
                 
        AnasaziMultiVec<TYPE>* AVj = _basisvecs->Clone(m); assert(_basisvecs);
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
                std::cout << " Arnoldi relation " << "for column " << i << " is " << fabs(ptr_norms[i]) << std::endl;        
	}
        std::cout << " " << std::endl;
                
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
	// These routines use an insertion sort routine to circument recursive calls.
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

#endif
// End of file BlockArnoldi.hpp

