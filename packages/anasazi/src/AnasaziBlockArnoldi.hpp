// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef BLOCK_ARNOLDI_HPP
#define BLOCK_ARNOLDI_HPP

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigenproblem.hpp"

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
    BlockArnoldi( Eigenproblem<TYPE>& problem, 
		  const TYPE tol=1.0e-6,
		  const int nev=5, 
		  const int length=25, 
		  const int block=1,
		  const string which="LM", 
		  const int step=25, 
		  const int restarts=0 
		  );
    
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
    void getEvecs(MultiVec<TYPE>& evecs); 
    
    //! This method puts the imaginary part of the computed eigenvectors in %ievecs.
    void getiEvecs(MultiVec<TYPE>& ievecs);
    
    //! This method returns the real part of the \c nev computed eigenvalues.
    TYPE * getEvals();
    
    //! This method returns a requested number of computed eigenvalues.
    /*! The input \c num can be greater than \c nev, but the method will only return
      the number of eigenvalues it has approximations to.  On exit, \c num will be the
      number of eigenvalue approximations that were returned 
    */
    TYPE * getEvals( int& num );
    
    //! This method returns the imaginary part of the computed eigenvalues.
    TYPE * getiEvals();
    
    //! This method returns a requested number of computed eigenvalues.
    /*! The input \c num can be greater than \c nev, but the method will only return
      the number of eigenvalues it has approximations to.  On exit, \c num will be the
      number of eigenvalue approximations that were returned 
    */
    TYPE * getiEvals( int& num );
    
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
    void QRFactorization( MultiVec<TYPE>&, Teuchos::SerialDenseMatrix<int,TYPE>& );
    void SortSchurForm( Teuchos::SerialDenseMatrix<int,TYPE>& H, Teuchos::SerialDenseMatrix<int,TYPE>& Q );
    void BlockReduction();
    void BlkOrth( MultiVec<TYPE>& Vec_in, const int j );
    void BlkOrthSing( MultiVec<TYPE>& Vec_in, const int j );
    void ComputeResiduals( const bool );
    void ComputeEvecs();
    void Restart();
    void SortEvals();
    void SetInitBlock();
    void SetBlkTols();
    void CheckBlkArnRed( const int j );
    void CheckSchurVecs( const int j ); 
    Eigenproblem<TYPE> &_problem; // must be passed in by the user
    MultiVec<TYPE> *_basisvecs, *_evecr, *_eveci;
    Teuchos::SerialDenseMatrix<int,TYPE> _hessmatrix;
    const int _nev, _length, _block, _restarts, _step;
    const TYPE _residual_tolerance;
    string _which;
    TYPE *_ritzresiduals, *_actualresiduals, *_evalr, *_evali;
    int *_order;
    int _restartiter, _iter, _jstart, _jend, _nevblock, _debuglevel, _defblock;
    int _offset, _maxoffset;
    bool _initialguess, _issym, _isdecompcurrent, _isevecscurrent, _exit_flg, _dep_flg;
    TYPE _schurerror, _scalefactor, _dep_tol, _blk_tol, _sing_tol, _def_tol;
  };
  //
  // Implementation
  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class TYPE>
  BlockArnoldi<TYPE>::BlockArnoldi(Eigenproblem<TYPE> & problem, 
				   const TYPE tol, 
				   const int nev, 
				   const int length, 
				   const int block,
				   const string which, 
				   const int step, 
				   const int restarts
				   ): 
    _problem(problem), 
    _basisvecs(0), 
    _evecr(0), 
    _eveci(0), 
    _hessmatrix(),
    _nev(nev), 
    _length(length), 
    _block(block), 
    _restarts(restarts),
    _step(step),
    _residual_tolerance(tol),
    _which(which),
    _ritzresiduals(0), 
    _actualresiduals(0),
    _evalr(0), 
    _evali(0), 
    _order(0),
    _restartiter(0), 
    _iter(0), 
    _jstart(0), 
    _jend(0), 
    _nevblock(0),
    _debuglevel(0),
    _defblock(0),
    _offset(0),
    _maxoffset(0),
    _initialguess(true), 
    _issym(false),
    _isdecompcurrent(false),
    _isevecscurrent(false),
    _exit_flg(false),
    _dep_flg(false),
    _schurerror(1.0), 
    _scalefactor(1.0),
    _dep_tol(1.0), 
    _blk_tol(1.0),
    _sing_tol(1.0),
    _def_tol(1.0)
    {     
    //
    // Determine _nevblock : how many blocks it will take to contain the _nev eigenvalues/vectors
    // NOTE: An additional block is kept if _nev is a perfect multiple of _block because of the
    // potential presence of complex eigenvalue pairs.  Additional blocks can be retained, up to
    // _maxoffset if the block ends with one eigenvalue of a complex conjugate pair.
    //
    _nevblock = _nev/_block + 1;
    _maxoffset = (_length-_nevblock)/2;
    //
    // Retrieve the initial vector from the Anasazi::Eigenproblem.
    //
    MultiVec<TYPE>* ivec = _problem.GetInitVec();
    assert(ivec!=NULL);
    
    assert(_length>0); assert(_block>0); assert(_step>0);
    //
    // Make room for the Arnoldi vectors and F.
    //
    _basisvecs = ivec->Clone((_length+1)*_block); assert(_basisvecs!=NULL);
    //
    // Make room for the eigenvectors
    //
    _evecr = ivec->Clone(_nev+1); assert(_evecr!=NULL);
    _eveci = ivec->Clone(_nev+1); assert(_eveci!=NULL);
    //
    // Create the rectangular Hessenberg matrix
    //
    _hessmatrix.shape((_length+1)*_block, _length*_block); 
    //
    // Create the vectors for eigenvalues and their residual errors and
    // initialize them.
    //
    _evalr = new TYPE[ _block*_length ]; assert(_evalr!=NULL);  
    _evali = new TYPE[ _block*_length ]; assert(_evali!=NULL);  
    _ritzresiduals = new TYPE[ _block*_length ]; assert(_ritzresiduals!=NULL);
    _actualresiduals = new TYPE[ _block*_length ]; assert(_actualresiduals!=NULL);
    _order = new int[ _block*_length ]; assert(_order!=NULL);
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
  
  template <class TYPE>
  void BlockArnoldi<TYPE>::SetBlkTols() {
    const TYPE two = 2.0;
    TYPE eps;
    char precision = 'P';
    Teuchos::LAPACK<int,TYPE> lapack;
    eps = lapack.LAMCH(precision);
    _blk_tol = 10*sqrt(eps);
    _sing_tol = 10 * eps;
    _dep_tol = 1/sqrt(two);
    _def_tol = eps;
  }
  
  template <class TYPE>
  BlockArnoldi<TYPE>::~BlockArnoldi() 
  {
    if (_basisvecs) delete _basisvecs;
    if (_evecr) delete _evecr;
    if (_eveci) delete _eveci;
    if (_ritzresiduals) delete [] _ritzresiduals;
    if (_actualresiduals) delete [] _actualresiduals;
    if (_evalr) delete [] _evalr;
    if (_evali) delete [] _evali;
    if (_order) delete [] _order;
  }
  
  template <class TYPE>
  void BlockArnoldi<TYPE>::getEvecs(MultiVec<TYPE>& evecs) {
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
  void BlockArnoldi<TYPE>::getiEvecs(MultiVec<TYPE>& ievecs) {
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
  TYPE * BlockArnoldi<TYPE>::getEvals( int& num ) {
    int i;
    
    // Correct the value of num if it's greater than the number of eigenvalue
    // approximations available.  If there was a restart recently, then there
    // may be more eigenvalue approximations than _jstart would lead you to
    // believe.
    switch ( _restarts ) {
    case 0 :
      if ( num > _jstart*_block ) {
	num = _jstart*_block;
      }
      break;
    default :
      if ( _jstart==_nevblock && num > _length*_block ) {
	num = _length*_block;
      }
      else if ( _jstart!=_nevblock && num > _jstart*_block ) {
		    num = _jstart*_block;
      }
      break;
    }
    
    // Now copy the eigenvalues.
    TYPE *temp_evals = new TYPE[ num ];
    for (i=0; i<num; i++) {
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
  TYPE * BlockArnoldi<TYPE>::getiEvals( int& num ) {
    int i;
    // Correct the value of num if it's greater than the number of eigenvalue
    // approximations available.  If there was a restart recently, then there
    // may be more eigenvalue approximations than _jstart would lead you to
    // believe.
    switch ( _restarts ) {
    case 0 :
      if ( num > _jstart*_block ) {
	num = _jstart*_block;
      }
      break;
    default :
      if ( _jstart==_nevblock && num > _length*_block ) {
	num = _length*_block;
      }
      else if ( _jstart!=_nevblock && num > _jstart*_block ) {
	num = _jstart*_block;
      }
      break;
    }
    
    // Now copy the eigenvalues.
    TYPE *temp_evals = new TYPE[ num ];
    for (i=0; i<num; i++) {
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
    cout<<"Iterations :\t"<<_iter<<endl;
    
    if (_restartiter > _restarts) 
      cout<<"Restarts :\t"<<_restartiter-1<<" of\t"<< _restarts<<endl;
    else
      cout<<"Restarts :\t"<<_restartiter<<" of\t"<< _restarts<<endl;
    
    cout<<"Block Size :\t"<<_block<<endl;
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
	  if (_exit_flg && _iter != _length+_restarts*(_length-_nevblock)) {
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
      cout<<"Eigenvalue\tRitz Residual"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for (i=0; i<_nevtemp; i++) {
	cout.width(10);
	cout<<_evalr[i]<<"\t"<<_ritzresiduals[i]<<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
    } else {
      cout<<"Real Part\tImag Part\tRitz Residual"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for (i=0; i<_nevtemp; i++) {
	cout.width(10);
	cout<<_evalr[i]<<"\t"<<_evali[i]<<"\t\t"<<_ritzresiduals[i]<<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
      cout<<" "<<endl;
    }
    cout<<"******************************************************"<<endl;
  }	

  template <class TYPE>
  void BlockArnoldi<TYPE>::SetInitBlock() {
    int i;
    int *index = new int[ _block ]; assert(index!=NULL);
    
    // This method will set the first block of _basisvecs to the initial guess,
    // if one is given, else it will fill the block with random vectors.
    
    if (_initialguess) {
      MultiVec<TYPE>* ivec = _problem.GetInitVec();
      assert(ivec!=NULL);
      const int cols = ivec->GetNumberVecs();
      if (cols < _block) {
	
	// Copy the given vectors in the first positions in the block
	// and fill the rest with random vectors.
	for (i=0; i<cols; i++) {
	  index[i] = i;
	}
	_basisvecs->SetBlock( *ivec, index, cols );			
	
	// Initialize the rest of the block with random vectors
	for (i=cols; i<_block; i++) {
	  index[i-cols] = i;
	}
	MultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block-cols);
	assert(U_vec!=NULL);
	U_vec->MvRandom();
	delete U_vec;
      }
      else {
	// Copy the first _block of the given vectors into the first _block
	// of _basisvecs, any additional vectors will be ignored.
	
	for (i=0; i<_block; i++) {
	  index[i] = i;
	}
	_basisvecs->SetBlock( *ivec, index, _block );
      }
    }
    else {
      // No initial guess is given, so initialize block with random vectors
      
      for (i=0; i<_block; i++) {
	index[i] = i;
      }
      MultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block);
      assert(U_vec!=NULL);
      U_vec->MvRandom();
      delete U_vec;		
    }
    
    // Clean up
    delete [] index;
  }
  
  template <class TYPE>
  void BlockArnoldi<TYPE>::solve () {
    //int rem_iters = _length+_restarts*(_length-_nevblock)-_iter;
    //
    // Right now the solver will just go the remaining iterations, but this design will allow
    // for checking of the residuals every so many iterations, independent of restarts.
    //
    while (_schurerror > _residual_tolerance && _restartiter <= _restarts && !_exit_flg) {
      iterate( _step );
    }
  }
  
  template <class TYPE>
  void BlockArnoldi<TYPE>::iterate(const int steps) {
    int i,temp;
    int tempsteps = steps;
    //
    // If this is the first steps of Block Arnoldi, initialize the first block of _basisvecs
    //
    if (!_iter) {
      SetInitBlock();
      int *index = new int[ (_length+1)*_block ]; assert(index!=NULL);
      for ( i=0; i<_block; i++ ) {
	index[i] = i;
      }
      MultiVec<TYPE>* U_vec = _basisvecs->CloneView(index,_block);
      assert(U_vec!=NULL);
      Teuchos::SerialDenseMatrix<int,TYPE> G10(_block,_block);
      QRFactorization( *U_vec, G10 );
      delete U_vec;
      delete [] index;
    }				
    //
    // Leave the iteration method now if the orthogonal subspace can't be extended.
    //
    if (_exit_flg) { return; }	
    //			
    // Now we go the number of steps requested by the user.  This may cause
    // a restart or hit the number of maximum iterations (restarts).  
    //
    while(tempsteps > 0 && _restartiter <= _restarts && !_exit_flg) {
      _isevecscurrent = false;
      // If we don't need to restart, just get it over with and return.
      if (_jstart+tempsteps < _length) {
	_jend = _jstart+tempsteps;
	_iter += tempsteps;
	tempsteps = 0;
	BlockReduction();
	if (_exit_flg && (_block > 1 )) { break; } 
	// We need to leave before we move the pointer if the blocksize is > 1,
	// otherwise it's a lucky breakdown.
	_jstart = _jend;  // Move the pointer
	ComputeResiduals( false );		
	_isdecompcurrent = false;
	// Output current information if necessary
	if (_debuglevel > 0) {
	  currentStatus();
	}
      }
      // Finish off this factorization and restart.
      else {  
	_jend = _length;
	temp = _length-_jstart;
	_iter += temp;
	tempsteps -= temp;
	BlockReduction();
	// We need to leave before we move the pointer if the blocksize is > 1,
	// otherwise it's a lucky breakdown.
	if (_exit_flg && (_block > 1)) { break; } 
	_jstart = _length; // Move the pointer
	//
	//  Compute the Schur factorization and prepare for a restart.  Don't
	//  compute restart if at end of iterations.
	//
	if (_restartiter < _restarts) {
	  ComputeResiduals( true );  
	  Restart();  
	  _isdecompcurrent = true;
	  _restartiter++;
	} else {
	  ComputeResiduals( false );
	  _restartiter++;
	  _isdecompcurrent = false;
	}
	// Output current information if necessary
	if (_debuglevel > 0) {
	  currentStatus();
	}
      }
    }
    //
    // Compute the current eigenvalue estimates before returning.
    //
    //cout<<"Upper Hessenberg matrix as of iteration :"<<_iter<<endl<<endl;       
    //cout<<_hessmatrix<<endl;
    
  }
  
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::BlockReduction () {
    int i,j;
    ReturnType ret;
    int *index = new int[ _block ]; assert(index!=NULL);
    MultiVec<TYPE> *U_vec=0, *F_vec=0;
    
    for ( j = _jstart; j < _jend; j++ ) {
      //
      // Associate the j-th block of _basisvecs with U_vec.
      //
      for ( i=0; i<_block; i++ ) {
	index[i] = j*_block+i;
      }
      U_vec = _basisvecs->CloneView(index, _block);
      assert(U_vec!=NULL);
      //
      // Associate (j+1)-st block of ArnoldiVecs with F_vec.
      //
      //for ( i=0; i<_block; i++ ) {
      //	index[i] = (j+1)*_block+i;
      //}
      F_vec = _basisvecs->Clone(_block);
      //F_vec = _basisvecs->CloneView(index, _block);
      //assert(F_vec!=NULL);
      //
      //  Compute F_vec = OP * U_vec
      //
      ret =_problem.ApplyOp( *U_vec, *F_vec ); 
      //
      // Use previous dependency information to decide which orthogonalization
      // method to use for the new block.  The global variable _dep_flg tells us
      // if we've had problems with orthogonality before.  If no problems have
      // been detected before we will use standard block orthogonalization.
      // However, if there are problems with this, we will use a more stringent
      // approach.
      //
      if (!_dep_flg) {
	BlkOrth( *F_vec, j );
      }
      //
      // If any block dependency was detected previously, then the more stringent
      // orthogonalization will be used.  If this method can't resolve the
      // dependency, then the _exit_flg will be set indicating that we can't proceed
      // any further.
      //			
      if (_dep_flg) {
	BlkOrthSing( *F_vec, j );
      }
      //
      delete U_vec; U_vec=0;
      delete F_vec; F_vec=0;
      //
      // If we cannot go any further with the factorization, then we need to exit
      // this method.
      //
      if (_exit_flg) { return; }
    }
    delete [] index;
  } // end BlockReduction()
  
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::BlkOrth( MultiVec<TYPE>& Vec_in, const int j ) {
    //
    // Orthogonalization is first done between the new block of
    // vectors and all previous blocks, then the vectors within the
    // new block are orthogonalized.
    //
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    const int max_num_orth = 2;
    int i, k, row_offset, col_offset;
    int * index = new int[ (_length+1)*_block ]; assert(index!=NULL);
    TYPE * norm1 = new TYPE[_block]; assert(norm1!=NULL);
    TYPE * norm2 = new TYPE[_block]; assert(norm2!=NULL);
    ReturnType ret; 
    //
    // Associate (j+1)-st block of ArnoldiVecs with F_vec.
    //
    for ( i=0; i<_block; i++ ) {
      index[i] = (j+1)*_block+i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
    assert(F_vec!=NULL);
    F_vec->MvAddMv(one, Vec_in, zero, Vec_in);
    //
    // Zero out the full block column of the Hessenberg matrix
    // even though we're only going to set the coefficients in
    // rows [0:(j+1)*_block-1]
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_block; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix( i, j*_block+k ) = zero;
      }
    }
    //
    // Grab all previous Arnoldi vectors
    //
    int num_prev = (j+1)*_block;
    for (i=0; i<num_prev; i++){
      index[i] = i;
    }
    MultiVec<TYPE>* V_prev = _basisvecs->CloneView(index,num_prev);
    assert(V_prev!=NULL);
    //
    // Create a matrix to store the product trans(V_prev)*B*F_vec
    //
    Teuchos::SerialDenseMatrix<int,TYPE> dense_mat(num_prev, _block );
    //
    F_vec->MvNorm(norm1);
    //
    // Perform two steps of block classical Gram-Schmidt so that
    // F_vec is B-orthogonal to the columns of V_prev.
    //
    for ( int num_orth=0; num_orth<max_num_orth; num_orth++ ) {
      //
      // Compute trans(V_prev)*B*F_vec and store in the j'th diagonal
      // block of the Hessenberg matrix
      //
      ret = _problem.BInProd( one, *V_prev, *F_vec, dense_mat );
      //
      // Update the orthogonalization coefficients for the j-th block
      // column of the Hessenberg matrix.
      //
      for ( k=0; k<_block; k++ ) {
	for ( i=0; i<num_prev; i++ ) {
	  _hessmatrix( i, j*_block+k ) += dense_mat(i,k);
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
	_dep_flg = true;
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
    if (!_dep_flg) {
      //
      // Compute the QR factorization of F_vec
      //
      row_offset = (j+1)*_block; col_offset = j*_block;
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _block, _block, 
							  row_offset, col_offset);
      
      QRFactorization( *F_vec, sub_block_hess );
    }
    //
    delete F_vec;
    delete V_prev;
    delete [] index;
    delete [] norm1;
    delete [] norm2;
    //
  }  // end BlkOrth()
  
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::BlkOrthSing( MultiVec<TYPE>& Vec_in, const int j ) {
    //
    // This is a variant of A. Ruhe's block Arnoldi
    // The orthogonalization of the vectors F_vec is done
    // one at a time. If a dependency is detected, a random
    // vector is added and orthogonalized against all previous
    // Arnoldi vectors.
    //
    const int IntOne = 1;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    int i, k, num_prev;
    int * index = new int[ (_length+1)*_block ]; assert(index!=NULL);
    Teuchos::SerialDenseVector<int,TYPE> dense_vec;
    TYPE norm1[IntOne];
    TYPE norm2[IntOne];
    ReturnType ret;
    //
    // Place the candidate vectors Vec_in into the (j+1)-st block of ArnoldiVecs.
    //
    for ( i=0; i<_block; i++ ) {
      index[i] = (j+1)*_block+i;
    }
    _basisvecs->SetBlock( Vec_in, index, _block ); 
    //
    // Zero out the full block column of the Hessenberg matrix
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_block; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix(i, j*_block+k) = zero;
      }
    }
    //
    MultiVec<TYPE> *q_vec=0, *Q_vec=0, *tptr=0;
    tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
    //
    // Start a loop to orthogonalize each of the _block
    // columns of the (j+1)-st block of _basisvecs against all 
    // the others.
    //
    for (int iter=0; iter<_block; iter++){
      num_prev = (j+1)*_block + iter; // number of previous _basisvecs
      dense_vec.size(num_prev);
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
      // Create matrix to store product trans(Q_vec)*B*q_vec
      //
      // Do one step of classical Gram-Schmidt B-orthogonalization
      // with a 2nd correction step if needed.
      //
      q_vec->MvNorm(norm1);
      //
      // Compute trans(Q_vec)*B*q_vec
      //
      ret = _problem.BInProd( one, *Q_vec, *q_vec, dense_vec );
      //
      // Sum results [0:num_prev-1] into column (num_prev-_block)
      // of the Hessenberg matrix
      //
      for (k=0; k<num_prev; k++){
	_hessmatrix(k, j*_block+iter) += dense_vec(k);
      }
      //
      // Compute q_vec<- q_vec - Q_vec * dense_vec
      //
      q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
      //
      q_vec->MvNorm(norm2);
      //
      if (norm2[0] < norm1[0] * _dep_tol) {
	//
	// Repeat process with newly computed q_vec
	//
	// Compute trans(Q_vec)*q_vec
	//
	ret = _problem.BInProd( one, *Q_vec, *q_vec, dense_vec );
	//
	// Sum results [0:num_prev-1] into column (num_prev-_block)
	// of the Hessenberg matrix
	//
	for (k=0; k<num_prev; k++){
	  _hessmatrix(k, j*_block+iter) += dense_vec(k);
	}
	//
	// Compute q_vec<- q_vec - Q_vec * dense_vec
	//
	q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
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
	// This code  is automatically doing 2 steps of B-orthogonalization
	// after adding a random vector. We could do one step of
	// orthogonalization with a correction step if needed.
	//
	for (int num_orth=0; num_orth<2; num_orth++){
	  ret = _problem.BInProd( one, *Q_vec, *tptr, dense_vec );
	  // Note that we don't change the entries of the
	  // Hessenberg matrix when we orthogonalize a
	  // random vector
	  tptr->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
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
	  _hessmatrix((j+1)*_block+iter, j*_block+iter) = zero;
	}
	else {
	  // Can't produce a new orthonormal basis vector
	  // Clean up and exit this block Arnoldi factorization!
	  _exit_flg = true;
	  delete [] index;
	  delete q_vec; q_vec=0;
	  delete Q_vec; Q_vec=0;
	  delete tptr; tptr=0;
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
	_hessmatrix((j+1)*_block+iter, j*_block+iter) = norm2[0];
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
  } // end BlkOrthSing()
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::QRFactorization (MultiVec<TYPE>& VecIn, 
					    Teuchos::SerialDenseMatrix<int,TYPE>& FouierR) {
    int i,j,k;
    int nb = VecIn.GetNumberVecs(); assert (nb == _block);
    int *index = new int[nb]; assert(index!=NULL);
    const int IntOne=1;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    bool addvec = false, flg = false;
    ReturnType ret;
    //
    TYPE norm1[IntOne];
    TYPE norm2[IntOne];
    MultiVec<TYPE> *qj = 0, *Qj = 0, *tptr = 0;
    tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
    //
    // Zero out the array that will contain the Fourier coefficients.
    //
    for ( j=0; j<nb; j++ ) {
      for ( i=0; i<nb; i++ ) {
	FouierR(i,j) = zero;
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
	Teuchos::SerialDenseVector<int,TYPE> rj(j);
	_problem.BMvNorm( *qj, norm1 );
	//
	// Do one step of classical Gram-Schmidt orthogonalization
	// with a second correction step if needed
	//
	// Determine the Fouier coefficients for B-orthogonalizing column
	// j of VecIn against columns 0:j-1 of VecIn. In other words,
	// result = trans(Qj)*B*qj.
	//
	ret = _problem.BInProd( one, *Qj, *qj, rj );
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
	_problem.BMvNorm( *qj, norm2 );			
	//
	if (norm2[0] < norm1[0] * _dep_tol){
	  //
	  // Repeat process with newly computed qj
	  //
	  ret = _problem.BInProd( one, *Qj, *qj, rj );
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
	  _problem.BMvNorm( *qj, norm2 );
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
	    _dep_flg = true;
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
	    Teuchos::SerialDenseVector<int,TYPE> tj(j);
	    //
	    tptr->MvRandom();
	    _problem.BMvNorm( *tptr, norm1 );
	    //
	    for (int num_orth=0; num_orth<2; num_orth++){
	      ret = _problem.BInProd( one, *Qj, *tptr, tj );
	      tptr->MvTimesMatAddMv(-one, *Qj, tj, one);
	    }
	    _problem.BMvNorm( *tptr, norm2 );
	    //
	    if (norm2[0] > norm1[0] * _sing_tol){
	      // Copy vector into current column of _basisvecs
	      qj->MvAddMv(one, *tptr, zero, *tptr);
	    }
	    else {
	      _exit_flg = true;
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
      _problem.BMvNorm( *qj, normq );
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
    } // for (j=0; j<nb; j++) ...
    //
    delete tptr;
    delete [] index;	
  }
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::ComputeResiduals( const bool apply ) {
    int i=0,j=0;
    int m = _jstart*_block, n=_jstart*_block;
    int mm1 = (_jstart-1)*_block;
    int _nevtemp, _nevtemp2, numimag;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::SerialDenseMatrix<int,TYPE> Q(n,n);
    Teuchos::LAPACK<int,TYPE> lapack;
    Teuchos::BLAS<int,TYPE> blas;
    //
    // If we are going to restart then we can overwrite the 
    // hessenberg matrix with the Schur factorization, else we
    // will just use a copy of it.  The Schur vectors will be in Q
    // on return.
    //
    if (apply) {
      Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::View, _hessmatrix, m, n);		
      SortSchurForm( Hj, Q );
      //
      // Determine new offset depending upon placement of conjugate pairs.	
      //
      _offset = _maxoffset;
      for (i=0; i<_maxoffset; i++) {
	numimag = 0;
	for (j=0; j<(_nevblock+i)*_block; j++) { 
	  if (_evali[j]!=zero) { numimag++; }; 
	}
	if (!(numimag % 2)) { _offset = i; break; }
      }
    } else {	
      //
      // Create a view into the current hessenberg matrix and make a copy.
      //
      Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::Copy, _hessmatrix, m, n);
      SortSchurForm( Hj, Q );
    }
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
    if (_jstart < _nevblock+_offset) {
      _nevtemp = n; _nevtemp2 = n;
    } else {
      _nevtemp = (_nevblock+_offset)*_block; _nevtemp2 = _nev;
    }
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _block, _block, m, mm1);
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_q(Teuchos::View, Q, _block, _nevtemp, mm1 );
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_b( _block, _nevtemp );
    blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _block, _nevtemp, _block, one, sub_block_hess.values(), sub_block_hess.stride(), 
	       sub_block_q.values(), sub_block_q.stride(), zero, sub_block_b.values(), _block );
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_b2(Teuchos::View, sub_block_b, _block, _nevtemp2);
    //
    //  Compute approximate ritzresiduals for each eigenvalue and the new scaling factor which
    //  will provide an approximate 2-norm to scale.
    //
    TYPE tempsf; 
    _scalefactor = lapack.LAPY2(_evalr[0],_evali[0]);
    for (i=1; i<n; i++) {
      tempsf = lapack.LAPY2(_evalr[i],_evali[i]);
      if (tempsf > _scalefactor) _scalefactor = tempsf;
    }
    _scalefactor = sqrt(_scalefactor);
    _schurerror = sub_block_b2.normFrobenius()/_scalefactor;
    //
    // ------------>  NOT SURE IF RITZRESIDUALS CAN BE UPDATED AFTER DEFLATION!
    //
    //for (i=0; i<_nevtemp ; i++) {
    for (i=_defblock*_block; i<_nevtemp ; i++) {
      Teuchos::SerialDenseMatrix<int,TYPE> s(Teuchos::View, sub_block_b, _block, 1, 0, i);
      _ritzresiduals[i] = blas.NRM2(_block, s.values(), 1)/_scalefactor;
    }   
    //
    //  We are going to restart, so update the Krylov-Schur decomposition.
    //
    if (apply) {	
      //
      // Update the Krylov basis.  Take into account that deflated blocks
      // need not be updated.
      //	
      int *index = new int[ n ]; assert(index!=NULL);
      for (i = 0; i < n; i++ ) {
	index[i] = i;
      }
      Teuchos::SerialDenseMatrix<int,TYPE> Qnev(Teuchos::View, Q, n, _nevtemp);
      MultiVec<TYPE>* basistemp = _basisvecs->CloneView( index, _nevtemp );
      MultiVec<TYPE>* basistemp2 = _basisvecs->CloneCopy( index, n );
      basistemp->MvTimesMatAddMv ( one, *basistemp2, Qnev, zero );
      //
      // Update the Krylov-Schur form (quasi-triangular matrix).
      //
      Teuchos::SerialDenseMatrix<int,TYPE> Hjp1(Teuchos::View, _hessmatrix,_block,_nevtemp, _nevtemp );
      for (i=0; i<_block; i++) {
	for (j=0; j<_nevtemp; j++) {
	  Hjp1(i, j) = sub_block_b(i, j);
	}
      }      
      //
      // Determine whether we need to continue with the computations.
      //
      if (_schurerror < _residual_tolerance )
	{
	  _exit_flg = true;
	}
      
      delete basistemp;
      delete basistemp2;
      delete [] index;
    }			
  }
  
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::ComputeEvecs() {
    int i=0,j=0,k=0;
    int n=_jstart*_block, info=0;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::LAPACK<int,TYPE> lapack;
    Teuchos::BLAS<int,TYPE> blas;
    MultiVec<TYPE>* basistemp=0;
    Teuchos::SerialDenseMatrix<int,TYPE> Q(n,n);
    int * index = new int [ n ]; assert(index!=NULL);
    //
    //  Set the index array.
    //
    for (k=0; k<n; k++) { index[k] = k; }
    //
    //  By some chance the number of eigenvalues requested may be more than
    //  the size of the factorization.  To prevent a seg fault in creating the
    //  eigenvectors, determine if this is the current situation.
    //
    int curr_nev = _nev;
    if ( n < _nev ) { curr_nev = n; }
    //
    //  Get a view into the current Hessenberg matrix.
    //
    Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::Copy, _hessmatrix, n, n, i, j);
    //
    //  If the Krylov-Schur decomposition is not current, compute residuals
    //  like we are going to restart to update decomposition.
    //
    if (!_isdecompcurrent) { 
      //
      // Update the Krylov basis using the Schur form.  
      // Take into account that deflated blocks need not be updated.
      //	
      SortSchurForm( Hj, Q );
      basistemp = _basisvecs->Clone( n );
      MultiVec<TYPE>* basistemp2 = _basisvecs->CloneView( index, n );
      basistemp->MvTimesMatAddMv ( one, *basistemp2, Q, zero );
      delete basistemp2;
    } else {
      //
      // We can aquire the Ritz vectors from the current decomposition.
      //
      basistemp = _basisvecs->CloneCopy( index, n );
    }
    //
    // Check the Schur form.
    //
    if (_debuglevel > 2 )
      CheckSchurVecs( _jstart );
    //
    //  If the operator is symmetric, then the Ritz vectors are the eigenvectors.
    //  So, copy the Ritz vectors.  Else, we need to compute the eigenvectors of the
    //  Schur form to compute the eigenvectors of the non-symmetric operator.
    //
    if (_issym) {
      _evecr->SetBlock( *basistemp, index, curr_nev );
    } else {  
      //
      //  Now compute the eigenvectors of the Schur form
      //  Reset the dense matrix and compute the eigenvalues of the Schur form.
      //
      int lwork = 4*n;
      TYPE *work = new TYPE[lwork]; assert(work!=NULL);
      int *select = new int[ n ];	  
      char * side = "R";
      char * howmny = "A";
      int mm, ldvl = 1;
      TYPE *vl = new TYPE[ ldvl ];
      lapack.TREVC( *side, *howmny, select, n, Hj.values(), Hj.stride(), vl, ldvl,
		    Q.values(), Q.stride(), n, &mm, work, &info );
      assert(info==0);
      delete [] work;
      delete [] select;
      delete [] vl;
      //
      //  Convert back to approximate eigenvectors of the operator and compute their norms.
      //
      TYPE* evecnrm = new double[ n ];
      MultiVec<TYPE>* evecstemp = _basisvecs->Clone( n );
      evecstemp->MvTimesMatAddMv( one, *basistemp, Q, zero );
      evecstemp->MvNorm( evecnrm );
      //
      // Sort the eigenvectors.
      //
      // Initialize imaginary part of eigenvector to zero, so we won't have to deal
      // with tracking it when complex eigenvalues are present.
      _eveci->MvInit(zero);
      
      int conjprs=0;
      int * indexi = new int [ curr_nev+1 ];
      MultiVec<TYPE> *evecstempr, *evecr1;
      TYPE t_evecnrm;
      i = 0;
      while ( i < curr_nev ) {	
	if (_evali[i] != zero) {
	  t_evecnrm = one/lapack.LAPY2(evecnrm[i],evecnrm[i+1]);
	  // Copy the real part of the eigenvector.  Scale by square-root of 2 to normalize the vector.
	  evecstempr = evecstemp->CloneView( index+i, 1 );
	  evecr1 = _evecr->CloneView( index+i, 1 );
	  evecr1->MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr );
	  delete evecr1; evecr1=0;
	  evecr1 = _evecr->CloneView( index+i+1, 1 );
	  evecr1->MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr );
	  // Note where imaginary part of eigenvector is.
	  indexi[conjprs] = i+1;
	  
	  // Increment counters.
	  conjprs++;
	  i = i+2; 				
	} else {
	  // Copy the real part of the eigenvector, scale to be norm one.
	  // We don't have to do anything for the imaginary
	  // part since we initialized the vectors to zero.
	  evecstempr = evecstemp->CloneView( index+i, 1 );
	  evecr1 = _evecr->CloneView( index+i, 1 );
	  evecr1->MvAddMv( one/evecnrm[i], *evecstempr, zero, *evecstempr );
	  // Increment counter.
	  i++;			
	}
	delete evecr1, evecr1=0;
	delete evecstempr; evecstempr=0;
      }
      // Set the imaginary part of the eigenvectors if conjugate pairs exist.
      // If the last eigenvector has a split conjugate pair, don't set negative imaginary
      // part.
      if (conjprs) {	
	MultiVec<TYPE>  *evecstempi=0, *eveci1=0;
	//
	// There is storage for an extra eigenvector.  
	// So, when the last eigenvalues is the first of a conjugate pair, that eigenvector will be computed.
	//
	for (i=0; i<conjprs; i++) {
	  t_evecnrm = one/lapack.LAPY2(evecnrm[indexi[i]],evecnrm[indexi[i]-1]);
	  evecstempi = evecstemp->CloneView( indexi+i, 1 ); 
	  eveci1 = _eveci->CloneView( indexi+i, 1 );
	  eveci1->MvAddMv( t_evecnrm*Teuchos::ScalarTraits<TYPE>::magnitude(_evali[indexi[i]])/_evali[indexi[i]],
			   *evecstempi, zero, *evecstempi );
	  delete eveci1; eveci1=0;
	  // Change index and set non-conjugate part of imag eigenvector.
	  indexi[i]--;
	  eveci1 = _eveci->CloneView( indexi+i, 1 );
	  eveci1->MvAddMv( t_evecnrm*Teuchos::ScalarTraits<TYPE>::magnitude(_evali[indexi[i]])/_evali[indexi[i]],
			   *evecstempi, zero, *evecstempi );
	  delete eveci1; eveci1=0;
	  delete evecstempi; evecstempi=0;
	}	      
      }
      
      // Clean up.
      delete evecstemp; 
      delete [] indexi; 
      delete [] evecnrm;
    }
    
    _isevecscurrent = true;
    delete [] index;
    delete basistemp;
  }
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::SortSchurForm( Teuchos::SerialDenseMatrix<int,TYPE>& H, Teuchos::SerialDenseMatrix<int,TYPE>& Q ) {
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::LAPACK<int,TYPE> lapack; 
    int i, j, info=0;
    int n = H.numRows(), ldh = H.stride(), ldq = Q.stride(); 
    TYPE* ptr_h = H.values();
    TYPE* ptr_q = Q.values();
    //
    //  If the operator is symmetric, analyze the block tridiagonal matrix
    //  and enforce symmetry.
    //
    if (_issym) {
      if (_restartiter > 0 && _restarts!=0) {
	//
	// The method has been restarted, so more caution must be used in
	// imposing symmetry.
	//
	for(j=_nevblock*_block; j<n; j++) {
	  for(i=0; i<j; i++) {
	    H( i, j ) = H( j, i );
	  }
	}
      } else {
	//
	// We haven't restarted, so just enforce symmetry throughout Hj.
	//
	for( j=0; j<n; j++ ) {
	  for( i=0; i<j; i++ ) {
	    H( i, j ) = H( j, i );
	  }
	}
      }
    }
    //
    //---------------------------------------------------
    // Compute the current eigenvalue estimates
    // ---> Use driver GEES to first reduce to upper Hessenberg 
    // 	form and then compute Schur form, outputting eigenvalues
    //---------------------------------------------------
    //
    int lwork = 4*n;
    TYPE *work = new TYPE[lwork]; assert(work!=NULL);
    int *select = new int[ n ];
    int sdim = 0; 
    int *bwork = new int[ n ];
    char * jobvs = "V";
    char * sort = "N";
    lapack.GEES( *jobvs, *sort, select, n, ptr_h, ldh, &sdim,_evalr,
		 _evali, ptr_q, ldq, work, lwork, bwork, &info );
    assert(info==0);
    //
    // Sort the eigenvalues, this also sorts the _order vector so we know
    // which ones we want. 
    //
    //cout<<"Before sorting the Schur form (H):"<<endl;
    //H.print(cout);	  
    SortEvals();
    //
    // Reorder real Schur factorization, remember to add one to the indices for the
    // fortran call and determine offset.  The offset is necessary since the TREXC
    // method reorders in a nonsymmetric fashion, thus we use the reordering in
    // a stack-like fashion.  Also take into account conjugate pairs, which may mess
    // up the reordering, since the pair is moved if one of the pair is moved.
    //
    int _nevtemp = 0;
    char * compq = "V";
    int *offset2 = new int[ n ]; assert(offset2!=NULL);
    int *_order2 = new int[ n ]; assert(_order2!=NULL);
    i = 0; 
    while (i < n) {
      if (_evali[i] != zero) {
	offset2[_nevtemp] = 0;
	for (j=i; j<n; j++) {
	  if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
	}
	_order2[_nevtemp] = _order[i];
	i = i+2;
      } else {
	offset2[_nevtemp] = 0;
	for (j=i; j<n; j++) {
	  if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
	}
	_order2[_nevtemp] = _order[i];
	i++;
      }
      _nevtemp++;
    }
    for (i=_nevtemp-1; i>=0; i--) {
      lapack.TREXC( *compq, n, ptr_h, ldh, ptr_q, ldq, _order2[i]+1+offset2[i], 
		    1, work, &info );
      assert(info==0);
    }
    //cout<<"After sorting and reordering the Schur form (H):"<<endl;
    //H.print(cout); 
    //
    // Determine largest off diagonal element of Schur matrix for symmetric case.
    //
    TYPE _maxsymmelem = zero;
    if (_issym) {
      for(j=0; j<n; j++){
	for(i=0; i<j; i++) {
	  if(Teuchos::ScalarTraits<TYPE>::magnitude(H(i, j))>_maxsymmelem) { _maxsymmelem = H(i, j); }
	}
      }
    }
    delete [] work; 
    delete [] bwork;
    delete [] select;
    delete [] offset2;
    delete [] _order2;
  }
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::Restart() {
    //  This method assumes the ComputeResiduals has been called before it
    //  to compute the Schur vectors and residuals.  This information is used to 
    //  restart the factorization.
    //
    int i,j;
    int _nevtemp = (_nevblock+_offset)*_block;
    int *index = new int[ _nevtemp ];
    //
    //  Move the F_vec block to the _jstart+1 position.	
    //
    for (i=0; i<_block; i++) {
      index[i] = _jstart*_block + i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneCopy(index, _block);
    for (i=0; i<_block; i++) {
      index[i] = _nevtemp + i;
    }
    _basisvecs->SetBlock( *F_vec, index, _block);
    //
    //  Check for blocks to deflate
    //  DEFLATION IS NOT READY RIGHT NOW!!!!!!!!
    //	int defcnt;
    //	i = _defblock;
    //	while ( i<_nevblock ) {
    //		defcnt = 0;
    //		for (j=0; j<_block; j++) {
    //			if (_ritzresiduals[i*_block+j] < _def_tol ) { defcnt++; }
    //		}
    //		if (defcnt == _block) {
    //			_defblock++;
    //		}
    //		i++;
    //	}		
    //
    //  If there are blocks to deflate, we need to set the subdiagonal entries to zero
    //
    if (_defblock > 0) {
      if (_debuglevel > 2) {
	cout<<"Deflating blocks with eigenvalue residuals below : "<<_def_tol<<endl;
	cout<<"Number of blocks being deflated : "<<_defblock<<endl;
      }
      TYPE zero = 0.0;
      Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::View, _hessmatrix, _block, _defblock*_block, _nevtemp, 0);
      for (i=0; i<_block; i++) {
	for (j=0; j<_defblock*_block; j++) {
	  Hj( i, j ) = zero;
	}
      }
    }
    //
    //  Reset the pointer.
    //
    _jstart = _nevblock+_offset; 
    //
    //  Clean up
    //
    delete F_vec;
    delete [] index;
  }
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::SortEvals() {
    int i, j, tempord;
    const int n = _jstart*_block;
    TYPE temp, tempr, tempi;
    Teuchos::LAPACK<int,TYPE> lapack;
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
	  temp = _evalr[j]*_evalr[j];
	  for (i=j-1; i>=0 && (_evalr[i]*_evalr[i])>temp; --i) {
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
	  temp = _evalr[j]*_evalr[j];
	  for (i=j-1; i>=0 && (_evalr[i]*_evalr[i])<temp; --i) {
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
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::CheckSchurVecs( const int j ) {
    //
    // Check the difference between the projection of A with the Schur vectors and the Schur matrix.
    // 
    int i, n = j*_block;
    int* index = new int[ n ];
    for( i=0; i<n; i++ ) { index[i] = i; } 
    TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    Teuchos::SerialDenseMatrix<int,TYPE> Hj( Teuchos::View, _hessmatrix, n, n );
    Teuchos::SerialDenseMatrix<int,TYPE> SchurProj( n, n );
    MultiVec<TYPE>* Z = _basisvecs->CloneView( index, n );
    MultiVec<TYPE>* basistemp = _basisvecs->Clone( n );
    _problem.ApplyOp( *Z, *basistemp );
    basistemp->MvTransMv( one, *Z, SchurProj );
    SchurProj.scale( -one );
    SchurProj += Hj;
    cout<< "Error in Schur Projection ( || (VQ)^T*A*(VQ) - S || ) at restart " << _restartiter+1 << " is "<< SchurProj.normFrobenius()/_scalefactor<<" (should be small)"<<endl;
  }
  
  
  template<class TYPE>
  void BlockArnoldi<TYPE>::CheckBlkArnRed( const int j ) {
    int i,k,m=(j+1)*_block;
    int *index = new int[m];
    ReturnType ret;       
    
    for ( i=0; i<_block; i++ ) {
      index[i] = m+i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
    assert(F_vec!=NULL);
    
    TYPE *ptr_norms = new TYPE[m];
    TYPE sum=0.0;
    
    F_vec->MvNorm(ptr_norms);
    for ( i=0; i<_block; i++ ) {
      sum += ptr_norms[i];
    }
    
    for ( i=0; i<m; i++ ) {
      index[i] = i;  
    }
    MultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
    assert(Vj!=NULL);   
    cout << " " << endl;
    cout << "********Block Arnoldi iteration******** " << j << endl;
    cout << " " << endl;
    
    const TYPE one=1.0;
    const TYPE zero=0.0;
    Teuchos::SerialDenseMatrix<int,TYPE> VTV(m,m);
    ret = _problem.BInProd( one, *Vj, *Vj, VTV );
    if (ret != Ok) { }
    TYPE* ptr=VTV.values();
    TYPE column_sum;
    
    for (k=0; k<m; k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
	if (i==k) {
	  ptr[i] -= one;
	}
	column_sum += ptr[i];
      }
      cout <<  " V^T*B*V-I " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
      ptr += m;
    }
    cout << " " << endl;
    
    Teuchos::SerialDenseMatrix<int,TYPE> E(m,_block);
    
    ret = _problem.BInProd( one, *Vj, *F_vec, E );
    if (ret != Ok) { }
    TYPE* ptr_Ej=E.values();
    
    for (k=0;k<_block;k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
	column_sum += ptr_Ej[i];
      }
      ptr_Ej += m;
      if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
      cout << " B-Orthogonality with F " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
    }
    cout << " " << endl;
                 
    MultiVec<TYPE>* AVj = _basisvecs->Clone(m); assert(AVj!=NULL);
    ret = _problem.ApplyOp(*Vj,*AVj);
    Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::View, _hessmatrix, m, m);
    AVj->MvTimesMatAddMv(-one, *Vj, Hj, one);
    for ( i=0; i<_block; i++ ) {  
      index[i] = j*_block+i;
    }
    
    MultiVec<TYPE>* Fj = AVj->CloneView(index, _block);
    Fj->MvAddMv(-one, *F_vec, one, *Fj);
    
    AVj->MvNorm(ptr_norms);
    
    for ( i=0; i<m; i++ ) { 
      cout << " Arnoldi relation " << "for column " << i << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(ptr_norms[i])/_scalefactor << endl;  
    }
    cout << " " << endl;
    
    delete F_vec;
    delete Fj;
    delete AVj;
    delete Vj;
    delete [] index;
    delete [] ptr_norms;
  }
  
} // End of namespace Anasazi
#endif
// End of file BlockArnoldi.hpp



