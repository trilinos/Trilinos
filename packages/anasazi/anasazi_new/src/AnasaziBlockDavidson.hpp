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

#ifndef ANASAZI_BLOCK_DAVIDSON_HPP
#define ANASAZI_BLOCK_DAVIDSON_HPP

#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

/*!	\class Anasazi::BlockDavidson

	\brief This class implements the Restarted Block Krylov Schur Method,
	an iterative method for solving eigenvalue problems.

	This method is a block version of the method presented by G.W. Stewart 
	in "A Krylov-Schur Algorithm for Large Eigenproblems", 
	SIAM J. Matrix Anal. Appl., Vol 28, No. 8, pp. 601-614.

	\author Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {
  
  template <class STYPE, class MV, class OP>
  class BlockDavidson : public Eigensolver<STYPE,MV,OP> { 
  public:
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockDavidson constructor.
    BlockDavidson( const Teuchos::RefCountPtr<Eigenproblem<STYPE,MV,OP> > &problem, 
		      const Teuchos::RefCountPtr<SortManager<STYPE,MV,OP> > &sm,
		      const Teuchos::RefCountPtr<OutputManager<STYPE> > &om,
		      const STYPE tol=1.0e-6,
		      const int length=25, 
		      const int step=25, 
		      const int maxIter=300 
		      );
    
    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockDavidson() {};
    //@}
    
    //@{ \name Solver application methods.
    
    /*! \brief This method uses iterate to compute approximate solutions to the
      original problem.  It may return without converging if it has taken the
      maximum number of iterations or numerical breakdown is observed.
    */
    void solve();
    //@}
    
    //@{ \name Solver status methods.
    
    //! This method returns the computed Ritz values.
    Teuchos::RefCountPtr<const std::vector<STYPE> > GetRitzValues() const { return(_ritzvalues); };
    
    //! This method returns the Ritz residuals for the computed eigenpairs.
    Teuchos::RefCountPtr<const std::vector<STYPE> > GetRitzResiduals() const { return(_ritzresiduals); };

    //! Get the current iteration count.
    int GetNumIters() const { return(_iter); };
    
    //! Get the current restart count of the iteration method.
    /*! Some eigensolvers can perform restarts (i.e.Arnoldi) to reduce memory
      and orthogonalization costs.  For other eigensolvers that don't
      perform restarts (i.e. LOBPCG), this is not a valid stopping criteria.
    */
    int GetNumRestarts() const { return(_restarts); };
    
    //! Get the solvers native residuals for the current eigenpairs. 
    /*! This is not be the same as the true residuals for most solvers. Sometimes the native
      residuals are not in multivector form, so the norm type is solver dependent.  
      
      \note
      <ol>
      <li> If the native residual is in multivector form then a non-null pointer will be
      returned, else the normvec will be populated with the current residual norms. 
      <li> If the native residual is returned in multivector form, the memory is managed
      by the calling routine.
      </ol>
    */
    Teuchos::RefCountPtr<const MV> GetNativeResiduals( STYPE* normvec ) const { return Teuchos::null; };
    
    /*! \brief Get a constant reference to the current linear problem, 
      which may include a current solution.
    */
    Eigenproblem<STYPE,MV,OP>& GetEigenproblem() const { return(*_problem); };
    
    //@}
    
    //@{ \name Output methods.
    
    //! This method requests that the solver print out its current status to screen.
    void currentStatus();
    //@}
  private:
    //
    // Internal methods
    //
    void accuracyCheck(const MV *X, const MV *MX, const MV *Q) const; 
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    Teuchos::RefCountPtr<Eigenproblem<STYPE,MV,OP> > _problem; 
    Teuchos::RefCountPtr<SortManager<STYPE,MV,OP> > _sm; 
    Teuchos::RefCountPtr<OutputManager<STYPE> > _om; 
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    Teuchos::RefCountPtr<OP> _BOp;
    Teuchos::RefCountPtr<OP> _Prec;
    Teuchos::RefCountPtr<MV> _evecs;
    const int _nev, _blockSize;  
    STYPE* _evals;
    //
    // Internal data.
    //
    const int _numBlocks, _restarts, _step, _maxIter;
    const STYPE _residual_tolerance;
    int _restartiter, _iter, _dimSearch;
    int _lwork;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<STYPE,MV,OP> _MSUtils;
    //
    // Internal storage for eigensolver
    //
    Teuchos::RefCountPtr<MV> _X, _MX, _KX;
    Teuchos::SerialDenseMatrix<int,STYPE> _K, _S;
    std::vector<STYPE> _theta, _normR, _work;

    typedef MultiVecTraits<STYPE,MV> MVT;
    typedef OperatorTraits<STYPE,MV,OP> OPT;
  };
  //
  // Implementation
  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class STYPE, class MV, class OP>
  BlockDavidson<STYPE,MV,OP>::BlockDavidson(const Teuchos::RefCountPtr<Eigenproblem<STYPE,MV,OP> > &problem, 
					   const Teuchos::RefCountPtr<SortManager<STYPE,MV,OP> > &sm,
					   const Teuchos::RefCountPtr<OutputManager<STYPE> > &om,
					   const STYPE tol, 
					   const int numBlocks, 
					   const int step, 
					   const int maxIter
					   ): 
    _problem(problem), 
    _sm(sm),
    _om(om),
    _Op(_problem->GetOperator()),
    _BOp(_problem->GetB()),
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()), 
    _nev(problem->GetNEV()), 
    _blockSize(problem->GetBlockSize()), 
    _evals(problem->GetEvals()), 
    _numBlocks(numBlocks), 
    _restarts(25),
    _step(step),
    _maxIter(maxIter),
    _residual_tolerance(tol),
    _restartiter(0), 
    _iter(0), 
    _dimSearch(0),
    _MSUtils(om)
  {     
    //
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    //
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    //
    // Define local block vectors
    //
    // MX = Working vectors (storing M*X if M is specified, else pointing to X)
    // KX = Working vectors (storing K*X)
    //
    // R = Residuals
    //
    _dimSearch = _blockSize*_numBlocks;
    _X = MVT::Clone( *iVec, _dimSearch + _blockSize );
    _KX = MVT::Clone( *iVec, _blockSize );
    _R = MVT::Clone( *iVec, _blockSize );
    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of _X.]
    //
    if (_BOp->get())
      _MX = MVT::Clone( *iVec, _blockSize );
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
	index[i] = i;
      _MX = MVT::CloneView( *_X, &index[0], _blockSize );
    }
    //
    // Initialize the workspace.
    //
    _X->MvRandom();
    //
    // Define dense local matrices and arrays
    //
    // theta = Storage for local eigenvalues     (size: dimSearch)
    // normR = Storage for the norm of residuals (size: blockSize)
    //
    // KK = Local stiffness matrix               (size: dimSearch x dimSearch)
    //
    // S = Local eigenvectors                    (size: dimSearch x dimSearch)
    //
    // tmpKK = Local workspace                   (size: dimSearch x dimSearch)
    //
    _theta.size( _dimSearch );
    _normR.size( _blockSize );
    //    
    _KK.size( _dimSearch, _dimSearch );
    _S.size( _dimSearch, _dimSearch );
        
  }
  
  template <class STYPE, class MV, class OP>
  void BlockDavidson<STYPE,MV,OP>::solve () 
  {
    int i, j;
    int info, nb;
    int bStart = 0, offset = 0;
    int knownEV = 0;
    int nFound = _blockSize;
    bool reStart = false;
    bool criticalExit = false;
    Teuchos::RefCountPtr<MV> Xcurrent, Xprev, Xtotal, Xnext;
    STYPE one = Teuchos::ScalarTraits<STYPE>::one();
    STYPE zero = Teuchos::ScalarTraits<STYPE>::zero();
    Teuchos::BLAS<int,STYPE> blas;
    Teuchos::LAPACK<int,STYPE> lapack;
    //
    // Determine the maximum number of blocks for this factorization.
    // ( NOTE:  This will be the _numBlocks since we don't know about already converged vectors here )
    int maxBlock = _numBlocks;
    //
    // Compute the required workspace for this factorization.
    //
    _lwork = lapack.ILAENV( 1, "geqrf", "", maxBlock*_blockSize, maxBlock*_blockSize );
    _lwork *= _blockSize;
    _work.size( _lwork );
    //
    while ( _iter <= _maxIter ) {
      
      for (nb = bStart; nb < maxBlock ; nb++) 
      {
	// 
	// Increment iteration counter
	//
	_iter++;
	if ( _iter > _maxIter )
	  break;
	//
	// Get a view of the current vectors being worked on.
	// The offset is dependent on the number of iterations been taken so far (nb)
	// and the number of known eigenvalues (knownEV).
	//
	int localSize = nb*_blockSize;
	std::vector<int> index( _blockSize );
	for (i=0; i < _blockSize; i++)
	  index[i] = localSize + knownEV + i;
	//
	Xcurrent = MVT::CloneView( *_X, &index[0], _blockSize );
	//
	index.resize( localSize );
	for (i=0; i < knownEV + localSize; i++)
	  index[i] = i;

	Xprev = MVT::CloneView( *_X, &index[0], knownEV + localSize );
	//
	// Apply the mass matrix.
	//	
	if (_BOp->get())
	  OPT::Apply( *_BOp, *Xcurrent, *_MX );
	//
	// Orthonormalize Xcurrent again the known eigenvectors and previous vectors.
	//
	if (nb == bStart) {
	  if (nFound > 0) {
	    if (knownEV == 0) {
	      info = _MSUtils.massOrthonormalize( Xcurrent, MX, _BOp->get(), Xprev, nFound, 2 );
	    }
	    else {
	      info = _MSUtils.massOrthonormalize( Xcurrent, MX, _BOp->get(), Xprev, nFound, 0 );
	    }
	  }
	  nFound = 0;
	} 
	else {
	  info = _MSUtils.massOrthonormalize( Xcurrent, MX, _BOp->get(), Xprev, _blockSize, 0 );
	}
	//
	// Exit the code if there has been a problem.
	//
	if (info < 0) 
	  return;
	//
	// Check orthogonality of X ( if required )
	//
	if ( 0 ) {
	  if (localSize > 0)
	    accuracyCheck( &Xcurrent, &MX, &Xprev );
	  else
	    accuracyCheck( &Xcurrent, &MX, Teuchos::null );
	}
	//
	// Apply the stiffness matrix.
	//
	OPT::Apply( *_Op, *Xcurrent, *_KX );
	//
	// Update the local stiffness matrix ( Xtotal^T * K * Xcurrent where Xtotal = [Xprev Xcurrent] )
	// Note:  Only the upper half of the matrix is stored in KK
	//
	index.resize( localSize + _blockSize );
	for (i=0; i < localSize + _blockSize; i++)
	  index[i] = knownEV + i;
	Xtotal = MVT::CloneView( *_X, &index[0], localSize + _blockSize );
	Teuchos::SerialDenseMatrix<int,STYPE> subKK( View, _KK, 0, localSize+_blockSize, _blockSize, 0, localSize );
	MVT::MvTransMv( one, *Xtotal, *_KX, KKsub );
	//
	// Perform spectral decomposition
	//
	info = _MSUtils.directSolver(localSize+_blockSize, _KK, 0, &_S, &theta, localSize+_blockSize, 10);
	//
	// Exit the code if there has been a problem.
	//
	if (info != 0) {
	  //
	  // Detect critical failure
	  // 
	  if (info<0) {
	    criticalExit = true;
	    break;
	  }
	  //
	  // Restart:  Spectral decomposition failed.
	  // Reinitialize all counters and randomize the starting block
	  //
	  reStart = true;
	  _restartIter++;
	  index.resize( _blockSize );
	  for (i=0; i<_blockSize; i++)
	    index[i] = knownEV + i;
	  Teuchos::RefCountPtr<MV> Xinit = MVT::CloneView( *_X, &index[0], _blockSize );
	  MVT::MvRandom( *Xinit );
	  nFound = _blockSize;
	  bStart = 0;
	  break;
	}
	//
	// Update the search space :
	// KX = Xtotal * S where S is the eigenvectors of the projected problem.
	//
	Teuchos::SerialDenseMatrix<int,STYPE> subS( View, _S, localSize+_blockSize, _blockSize );
	MVT::MvTimesMatAddMv( one, *_Xtotal, subS, zero, *_KX );
	//
	// Apply the mass matrix for the next block
	// 
	if (_BOp->get())
	  OPT::Apply( *_BOp, *_KX, *_MX );
	//
	// Apply the stiffness matrix for the next block
	//
	OPT::Apply( *_Op, *_KX, *_R );
	//
	// Compute the residual :
	// R = KX - diag(theta)*MX
	// 
	Teuchos::SerialDenseMatrix<int,STYPE> D(_blockSize, _blockSize);
	for (i=0; i<_blockSize; i++ )
	  D(i,i) = -theta[i];
	//
	if (_BOp->get()) {
	  MVT::MvTimesMatAddMv( one, *_MX, D, one, *_R );
	}
	else {
	  MVT::MvTimesMatAddMv( one, *_KX, D, one, *_R );
	}
	_problem->MvNorm( *_R, &_normR[0] );
	//
	// Scale the norms of residuals with the eigenvalues and check for converged eigenvectors.
	//
	nFound = 0;
	for (j=0; j<_blockSize; j++) {
	  // Scale eigenvalues if _theta is non-zero.
	  if ( _theta[j] != zero )
	    _normR[j] /= _theta[j];
	  // Check for convergence
	  if (_normR[j] < _residual_tolerance)
	    nFound ++;	  
	}
	//
	// Exit the loop to treat the converged eigenvectors
	//
	if (nFound > 0) {
	  nb += 1;
	  offSet = 0;
	  break;
	}
	//
	// Apply the preconditioner on the residuals
	//
	if (maxBlock == 1) {
	  if (_Prec->get()) {
	    OPT::Apply( *_Prec, *_R, *Xcurrent );
	  }
	  else
	    MVT::MvAddMv( one, *_R, zero, *_R, *Xcurrent );
	  //
	  // Update the preconditioned residual 
	  //
	  MVT::MvAddMv( one, *_KX, -one, *Xcurrent, *Xcurrent );
	  break;
	} // if (maxBlock == 1)

	if (nb == maxBlock - 1) {
	  nb++;
	  break;
	}
	//
	// Prepare next block in factorization
	//
	index.resize( _blockSize );
	for( i=0; i<_blockSize; i++) 
	  index[i] = knownEV + localSize + _blockSize + i;
	Xnext = MVT::CloneView( *_X, &index[0], _blockSize );
	if (_Prec) {
	  OPT::Apply( *_Prec, *_R, *Xnext );
	}
	else 
	  MVT::MvAddMv( one, *_R, zero, *_R, *Xnext );
	
      } // for (nb = bStart; nv < maxBlock; nb++)

      //
      // Check if there is any reason we need to skip the rest of the while loop.
      //
      if (_iter > _maxIter)
	break;

      if (reStart == true) {
	reStart = false;
	continue;
      }

      if (criticalExit == true)
	break;

      //
      // Store the final converged eigenvectors
      //
      if (knownEV + nFound >= _nev) {
	std::vector<int> index2(1);
	for (j=0; j<_blockSize; j++) {
	  if (_normR[j] < _residual_tolerance) {
	    index2[0] = j;
	    Teuchos::RefCountPtr<MV> tmp_KX = MVT::CloneView( *_KX, &index[0], 1 );
	    index2[0] = knownEV;
	    MVT::SetBlock( *tmp_KX, &index2[0], 1, *_X );
	    _evals[knownEV] = theta[j];
	    knownEV++;
	  }
	}
	break;
      } // if (knownEV + nFound >= _nev)      
      //
      // Store the converged eigenvectors and define the restarting block
      // in the particular case of 1 block ( maxBlock == 1 ).
      //
      if (maxBlock == 1) {
	if (nFound > 0) {
	  std::vector<int> index2(1);
	  int tmp_ptr = knownEV + nFound;
	  nFound = 0;
	  for (j=0; j<_blockSize; j++) {
	    //
	    // Get a view of the current prospective eigenvector.
	    //
	    index2[0] = j;
	    Teuchos::RefCountPtr<MV> tmp_KX = MVT::CloneView( *_KX, &index[0], 1 );
	    if (_normR[j] < _residual_tolerance) {
	      index2[0] = knownEV;
	      MVT::SetBlock( *tmp_KX, &index2[0], 1, *_X );	      
	      _evals[knownEV] = theta[j];
	      knownEV++;
	      nFound++;	      
	    }
	    else {
	      index2[0] = tmp_ptr + (j-nFound);
	      MVT::SetBlock( *tmp_KX, &index2[0], 1, *_X );
	    }
	  } // for (j=0; j<_blockSize; j++)
	  //
	  index.resize( nFound );
	  for (i=0; i<nFound; i++) 
	    index[i] = knownEV + _blockSize - nFound + i;
	  Xnext = MVT::CloneView( *_X, &index[0], nFound );
	  MVT::MvRandom( *Xnext );
	}
	else {
	  nFound = _blockSize;
	}
	continue;
      } // if (maxBlock==1)
      //
      // Define the restarting block when maxBlock > 1
      //
      if (nFound > 0) {
	int firstIndex = _blockSize;
	//
	// Find the first index where the unconverged eigenvectors start.
	//
	for (j=0; j<_blockSize; j++) {
	  if (_normR[j] >= _residual_tolerance) {
	    firstIndex = j;
	    break;
	  }
	}
	//
	// If the firstIndex comes before the number found, then we need
	// to move the converged eigenvectors to the front of the spectral
	// transformation.
	//
	while (firstIndex < nFound) {
	  for (j=firstIndex; j<_blockSize; j++) {
	    if (_normR[j] < _residual_tolerance) {
	      //
	      // Swap and j-th and firstIndex-th position
	      //
	      blas.SWAP(nb*_blockSize, _S.values() + j*_dimSearch, 1, _S.values() + firstIndex*_dimSearch, 1 );
	      blas.SWAP(1, &_theta[0] + j, 1, &_theta[0] + firstIndex, 1 );
	      blas.SWAP(1, &_normR[0] + j, 1, &_normR[0] + firstIndex, 1 );
	      break;
	    }
	  } 
	  for (j=0; j<_blockSize; j++) {
	    if (_normR[j] >= _residual_tolerance) {
	      firstIndex = j;
	      break;
	    }
	  }
	} // while (firstIndex < nFound)
	//
	// Copy the converged eigenvalues from "theta" to "evals".
	//
	blas.COPY( nFound, &_theta[0], 1, _evals + knownEV, 1 ); 

      } // if (nFound > 0)

      //
      // Define the restarting size
      //
      bStart = ((nb - offSet) > 2 ) ? (nb - offSet)/2 : 0;
      //
      // Define the restarting space and local stiffness matrix
      //
      _KK.PutScalar( zero );
      for (j=0; j<bStart*_blockSize; j++)
	KK(j,j) = _theta[j + nFound];
      //
      // Form the restarting space
      //
      int oldCol = nb*_blockSize;
      int newCol = nFound + (bStart+1)*_blockSize;
      newCol = (newCol > oldCol) ? oldCol : newCol;
      lapack.GEQRF(oldCol, newCol, _S.values(), _S.stride(), &_theta[0], &_work[0], _lwork, &info);
      lapack.ORGQR(oldCol, newCol, newCol, _S.values(), _S,stride(), &_theta[0], &_work[0], _lwork, &info);      
      index.resize( oldCol );
      for (i=0; i<oldCol; i++)
	index[i] = knownEV + i;
      Teuchos::RefCountPtr<MV> oldX = MVT::CloneView( *_X, &index[0], oldCol );
      index.resize( newCol );
      for (i=0; i<newCol; i++)
	index[i] = knownEV + i; 
      Teuchos::RefCountPtr<MV> newX = MVT::CloneView( *_X, &index[0], newCol );
      Teuchos::RefCountPtr<MV> temp_newX = MVT::Clone( *_X, newCol );
      Teuchos::SerialDenseMatrix<int,STYPE> _Sview( View, _S, oldCol, newCol );
      MVT::MvTimesMatAddMv( one, *oldX, _Sview, zero, *temp_newX );
      MVT::MvAddMv( one, *temp_newX, zero, *temp_newX, *newX ); 
      
      if (nFound == 0)
	offSet++;
      
      knownEV += nFound;
      maxBlock = (dimSearch/_blockSize) - (knownEV/_blockSize);
      //
      // Put random vectors if the Rayleigh Ritz vectors are not enough
      // 
      newCol = nFound + (bStart+1)*_blockSize;
      if (newCol > oldCol) {
	index.resize( nFound );
	for (i=0; i<nFound; i++)
	  index[i] = knownEV + _blockSize - nFound + i;
	Xnext = MVT::CloneView( *_X, &index[0], nFound );
	MVT::MvRandom( *Xnext );
	continue;
      }
      //
      // Reset nFound counter.
      // 
      nFound = 0;
      //
    } // while (_iter < _maxIter)
    //
    // Sort the eigenvectors
    //
    if ((info==0) && (knownEV > 0))
      _MSUtils.sortScalars_Vectors(knownEV, _evals, _evecs);
 
  } // end solve()


  template <class STYPE, class MV, class OP>
  void BlockDavidson<STYPE,MV,OP>::accuracyCheck(const MV *X, const MV *MX, const MV *Q) const 
    {
      
      cout.precision(2);
      cout.setf(ios::scientific, ios::floatfield);
      double tmp;
      
      int myPid = MyComm.MyPID();
      
      if (X) {
	if (M) {
	  if (MX) {
	    tmp = _MSUtils.errorEquality(X, MX, M);
	    if (myPid == 0)
	      cout << " >> Difference between MX and M*X = " << tmp << endl;
	  }
	  tmp = _MSUtils.errorOrthonormality(X, M);
	  if (myPid == 0)
	    cout << " >> Error in X^T M X - I = " << tmp << endl;
	}
	else {
	  tmp = _MSUtils.errorOrthonormality(X, 0);
	  if (myPid == 0)
	    cout << " >> Error in X^T X - I = " << tmp << endl;
	}
      }
      
      if (Q == 0)
	return;
      
      if (M) {
	tmp = _MSUtils.errorOrthonormality(Q, M);
	if (myPid == 0)
	  cout << " >> Error in Q^T M Q - I = " << tmp << endl;
	if (X) {
	  tmp = _MSUtils.errorOrthogonality(Q, X, M);
	  if (myPid == 0)
	    cout << " >> Orthogonality Q^T M X up to " << tmp << endl;
	}
      }
      else {
	tmp = _MSUtils.errorOrthonormality(Q, 0);
	if (myPid == 0)
	  cout << " >> Error in Q^T Q - I = " << tmp << endl;
	if (X) {
	  tmp = _MSUtils.errorOrthogonality(Q, X, 0);
	  if (myPid == 0)
	    cout << " >> Orthogonality Q^T X up to " << tmp << endl;
	}
      }
      
    }
  
    } // End of namespace Anasazi

#endif

// End of file AnasaziBlockDavidson.hpp



