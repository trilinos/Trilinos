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

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"
#include "BelosStatusTest.hpp"

/*!	\class Belos::BlockCG

	\brief This class implements the preconditioned Conjugate Gradient algorithm for
	solving real symmetric positive definite linear systems of equations
	AX = B, where B is a matrix containing one or more right-hand sides.

	\author Teri Barth and Heidi Thornquist
*/

namespace Belos {

template <class TYPE>
class BlockCG : public IterativeSolver<TYPE> { 
public:
  //@{ \name Constructor/Destructor.
  //! %Belos::BlockCG constructor.
  BlockCG(LinearProblemManager<TYPE>& lp,
	  StatusTest<TYPE>& stest,
	  OutputManager<TYPE>& om);
  
  //! %BlockCG destructor.
  virtual ~BlockCG();
  //@}
  
  //@{ \name Accessor methods
  
  //! Get the iteration count for the current block of linear systems.
  int GetNumIters() const { return( _iter ); }
  
  //! Get the restart count of the iteration method for the current block of linear systems [not valid for CG].
  int GetNumRestarts() const { return(0); }
   
  //! Get the solvers native residuals for the current block of linear systems.
  /*! 
   */
  ReturnType GetNativeResidNorms(TYPE *normvec, NormType norm_type) const;
  
  //! Get the actual residual vectors for the current block of linear systems.
  /*! This may force the solver to compute a current residual for its linear
  	systems.  For CG, this method is only useful when the blocksize is larger
	than the current number of linear systems being solved for.  Otherwise,
	the solution held in the linear problem manager is current.
  */
  MultiVec<TYPE>* GetCurrentSoln() { return _cur_block_sol; };
  
  //! Get a constant reference to the current linear problem.  
  /*! This may include a current solution, if the solver has recently restarted or completed.
   */
  LinearProblemManager<TYPE>& GetLinearProblem() const { return( _lp ); }

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
  bool QRFactorDef(MultiVec<TYPE>&, Teuchos::SerialDenseMatrix<int,TYPE>&,
		   int[], int&);
  void CheckCGOrth(MultiVec<TYPE>&, MultiVec<TYPE>&);
  void PrintCGIterInfo(int[], const int);

  //! Linear problem manager [ must be passed in by the user ]
  LinearProblemManager<TYPE>& _lp; 

  //! Status test [ must be passed in by the user ]
  StatusTest<TYPE>& _stest; 

  //! Output manager [ must be passed in by the user ]
  OutputManager<TYPE>& _om;

  //! Pointer to current linear systems block of solution vectors [obtained from linear problem manager]
  MultiVec<TYPE> *_cur_block_sol;

  //! Pointer to current linear systems block of right-hand sides [obtained from linear problem manager]
  MultiVec<TYPE> *_cur_block_rhs; 

  //! Vector of the current residual norms.
  TYPE * _cur_resid_norms; 

  //! Vector of the initial residual norms.
  TYPE * _init_resid_norms;

  //! Current blocksize and iteration number.
  int _blocksize, _iter;

  //! Numerical breakdown tolerances.
  TYPE _prec, _dep_tol;

  //! Output stream.
  ostream& _os;
};

//
// Implementation
//

template <class TYPE>
BlockCG<TYPE>::BlockCG(LinearProblemManager<TYPE>& lp,
		       StatusTest<TYPE>& stest,
		       OutputManager<TYPE>& om) : 
  _lp(lp), 
  _stest(stest),
  _om(om),
  _cur_block_rhs(0),
  _cur_block_sol(0),
  _blocksize(0), 
  _cur_resid_norms(0), 
  _init_resid_norms(0), 
  _iter(0),
  _prec(5.0e-15), 
  _dep_tol(0.75),
  _os(om.GetOStream())
{ 
  //
  // Set the block orthogonality tolerances
  //
  SetCGBlkTols();
}

template <class TYPE>
BlockCG<TYPE>::~BlockCG() 
{
 if (_cur_resid_norms) delete _cur_resid_norms;
 if (_init_resid_norms) delete _init_resid_norms;
}

template <class TYPE>
void BlockCG<TYPE>::SetCGBlkTols() 
{
  const TYPE two = 2.0;
  TYPE eps;
  char precision = 'P';
  Teuchos::LAPACK<int,TYPE> lapack;
  eps = lapack.LAMCH(precision);
  _prec = eps;
  _dep_tol = 1/sqrt(two);
}

template <class TYPE>
ReturnType BlockCG<TYPE>::GetNativeResidNorms(TYPE *normvec, NormType norm_type) const 
{
  if (normvec) {
    for (int i=0; i<_blocksize; i++) {
      normvec[i] = _cur_resid_norms[i];
    }
    return Ok;
  }
  return Error;
}

template <class TYPE>
void BlockCG<TYPE>::Solve () 
{
  //
  int i, j, k, info, num_ind;
  int ind_blksz, prev_ind_blksz;
  bool exit_flg = false;
  char UPLO = 'U';
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  Teuchos::LAPACK<int,TYPE> lapack;
  MultiVec<TYPE> *_basisvecs=0, *_residvecs=0;
  MultiVec<TYPE> *R_prev=0, *R_new=0, *P_prev=0, *P_new=0;
  MultiVec<TYPE> *AP_prev = 0, *temp_blk=0;
  //
  // Retrieve the first linear system to be solved.
  //
  _cur_block_sol = _lp.GetCurrLHSVec();
  _cur_block_rhs = _lp.GetCurrRHSVec();
  //
  //  Start executable statements. 
  //
  while (_cur_block_sol && _cur_block_rhs ) {
    //
    if (_om.doOutput( 0 )) {
      _os << endl;
      _os << "===================================================" << endl;
      _os << "Solving linear system(s):  " << _lp.GetRHSIndex() << " through " << _lp.GetRHSIndex()+_lp.GetNumToSolve() << endl;
      _os << endl;
    }	
    //
    // Get the blocksize for this set of linear systems.
    //
    _blocksize = _lp.GetBlockSize();
    //
    // Make additional space needed during iteration
    //
    int *index2 = new int[ _blocksize ]; assert(index2!=NULL);
    int *index = new int[ _blocksize ]; assert(index!=NULL);
    int *ind_idx = new int[_blocksize]; assert(ind_idx!=NULL);
    int *cols = new int[_blocksize]; assert(cols!=NULL);
    //
    // Make room for the direction, residual, and operator applied to direction vectors.
    // We save 3 blocks of these vectors.  Also create 2 vectors of _blocksize to store
    // The residual norms used to determine valid direction/residual vectors.
    //
    _basisvecs = _cur_block_sol->Clone(2*_blocksize); assert(_basisvecs!=NULL);
    _residvecs = _cur_block_sol->Clone(2*_blocksize); assert(_residvecs!=NULL);
    temp_blk = _cur_block_sol->Clone(_blocksize); assert(temp_blk!=NULL);
    _cur_resid_norms = new TYPE[_blocksize]; assert(_cur_resid_norms!=NULL);
    _init_resid_norms= new TYPE[_blocksize]; assert(_init_resid_norms!=NULL);
    //
    ind_blksz = _blocksize;
    prev_ind_blksz = _blocksize;
    Teuchos::SerialDenseMatrix<int,TYPE> alpha( _blocksize, _blocksize );
    Teuchos::SerialDenseMatrix<int,TYPE> beta( _blocksize, _blocksize );
    Teuchos::SerialDenseMatrix<int,TYPE> T2( _blocksize, _blocksize );
    //
    for (i=0;i<_blocksize;i++){
      index[i] = i;
      ind_idx[i] = i; 
      cols[i] = i;
    }
    //
    // ************ Compute the initial residuals ********************************
    //
    // Associate the first block of _basisvecs with P_prev and the
    // first block of _residvecs with R_prev
    //
    P_prev = _basisvecs->CloneView(ind_idx, _blocksize); assert(P_prev!=NULL);
    R_prev = _residvecs->CloneView(ind_idx, _blocksize); assert(R_prev!=NULL);
    AP_prev = temp_blk->CloneView(ind_idx, _blocksize); assert(AP_prev!=NULL);
    //
    // Store initial guesses to AX = B in 1st block of _basisvecs
    //         P_prev = one*cur_block_sol + zero*P_prev
    //
    P_prev->MvAddMv(one, *_cur_block_sol, zero, *P_prev);
    //
    // Multiply by A and store in AP_prev
    //       AP_prev = A*P_prev
    //
    _lp.ApplyOp( *P_prev, *AP_prev );
    //
    // Compute initial residual block and store in 1st block of _residvecs
    //     R_prev = cur_block_rhs - A*P_prev
    //
    R_prev->MvAddMv(one, *_cur_block_rhs, -one, *AP_prev);
    //
    //-------Compute and save the initial residual norms----------
    //
    R_prev->MvNorm(_init_resid_norms);
    //
    // Update indices of current (independent) blocks.
    // If a residual is too small, it will be dropped from
    // the current block, thus, from future computations
    //
    k = 0; j = 0;
    for (i=0; i<_blocksize; i++){
      _cur_resid_norms[i] = _init_resid_norms[i];
      if (_init_resid_norms[i] > _prec){
	ind_idx[k] = i;
	k = k+1;
      }
    }
    ind_blksz = k; 
    //
    if (ind_blksz > 0) { 
      //
      // All initial residuals have not converged -- continue Block CG	
      // Compute the initial block of direciton vectors
      //
      // Associate current blocks of residuals, directions, and solution block
      // with R_prev, P_prev, and cur_sol
      //
      delete R_prev; R_prev=0;
      delete P_prev, P_prev=0;
      R_prev = _residvecs->CloneView( ind_idx, ind_blksz );
      P_prev = _basisvecs->CloneView( ind_idx, ind_blksz );
      //
      //----------------Compute initial direction vectors--------------------------
      // Initially, they are set to the preconditioned residuals
      //
      if (_lp.ApplyLeftPrec( *R_prev, *P_prev ) != Ok ) { P_prev->MvAddMv( one , *R_prev, zero, *R_prev); }
      //
      // Compute an orthonormal block of initial direction vectors,
      // and check for dependencies, adjusting indices of independent
      // vectors if needed
      //
      Teuchos::SerialDenseMatrix<int,TYPE> G(ind_blksz, ind_blksz);
      num_ind = 0; exit_flg = false;
      exit_flg = QRFactorDef(*P_prev, G, cols, num_ind);
      //
      if ( exit_flg ) {
	if (_om.doOutput( 0 )) {
	  _os << " Exiting Block CG iteration " << endl; 
	  _os << " Reason: No independent initial direction vectors" << endl;
	}
      }		
      if (num_ind < ind_blksz) {
	// The initial block of direction vectors are linearly dependent
	if (_om.doOutput( 0 )) {
	  _os << " Initial direction vectors are dependent" << endl;
	  _os << " Adjusting blocks and indices for iteration" << endl;
	}
	//
	ind_blksz = num_ind;
	for (i=0; i< ind_blksz; i++){
	  ind_idx[i] = ind_idx[cols[i]];
	}	
      }  // end if (num < ind_blksz)
    }  // end if (ind_blksz > 0)
    //		
    else {  // all initial residuals have converged
      if (_om.doOutput( 0 )) {
	_os << " Exiting Block CG iteration " 
	     << " -- Iteration# " << _iter << endl;
	_os << "Reason: All initial residuals have converged" << endl;
      }
      exit_flg = true;
    }

    // ***************************************************************************
    // ************************Main CG Loop***************************************
    // ***************************************************************************
    // 
    if (_om.doOutput( 2 )) _os << "Entering main CG loop" << endl << endl;
    //
    int new_blk = 1;
    for (_iter=0; _stest.CheckStatus(this) == Unconverged && !exit_flg; _iter++) 
{
      //
      //  Clean up before computing another iterate.
      //  Current information is kept for the StatusTest.
      //
      delete P_prev; P_prev=0;
      delete R_prev; R_prev=0;
      if (P_new) delete P_new; P_new=0;
      if (R_new) delete R_new; R_new=0;
      //
      //----------------Compute the new blocks of iterates and residuals------------------
      //
      // Get views of the previous blocks of residuals, direction vectors, etc.
      //
      if (new_blk){
	P_prev = _basisvecs->CloneView( ind_idx, ind_blksz );
      }
      else {
        for (i=0; i< ind_blksz; i++) {
	  index2[i] = _blocksize + ind_idx[i];
        } 
	P_prev = _basisvecs->CloneView(index2,ind_blksz);
      }
      //
      for (i=0; i < _blocksize; i++){
	index2[i] = _blocksize + i;
      }
      if (new_blk){
	R_prev = _residvecs->CloneView( index, _blocksize );
	R_new = _residvecs->CloneView( index2, _blocksize );
      }
      else {
	R_prev = _residvecs->CloneView( index2, _blocksize );
	R_new = _residvecs->CloneView( index, _blocksize );
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
	delete AP_prev; 
	AP_prev = temp_blk->CloneView( ind_idx, ind_blksz ); 
	alpha.reshape( ind_blksz, _blocksize );
	T2.reshape( ind_blksz, ind_blksz );
      }
      _lp.ApplyOp( *P_prev, *AP_prev );
      R_prev->MvTransMv(one, *P_prev, alpha);   
      AP_prev->MvTransMv(one, *P_prev, T2);
      //
      // 2)
      lapack.POTRF(UPLO, ind_blksz, T2.values(), ind_blksz, &info);
      if (info != 0) {
	if(_om.doOutput( 0 )){
	  _os << " Exiting Block CG iteration "
	       << " -- Iteration# " << _iter << endl;
	  _os << " Reason: Cannot compute coefficient matrix alpha" << endl;
	  _os << " P_prev'* A*P_prev is singular" << endl;
	  _os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      // 3)
      lapack.POTRS(UPLO, ind_blksz, _blocksize, T2.values(), ind_blksz, alpha.values(), ind_blksz, &info);
      // Note: solution returned in alpha
      if (info != 0) {
	if(_om.doOutput( 0 )){
	  _os << " Exiting Block CG iteration "
	       << " -- Iteration# " << _iter << endl;
	  _os << " Reason: Cannot compute coefficient matrix alpha" << endl;
	  _os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Update the solution: cur_sol = one*cur_sol + one*P_prev*alpha
      // 
      _cur_block_sol->MvTimesMatAddMv(one, *P_prev, alpha, one);
      _lp.SolutionUpdated();	// Inform the linear problem that the solution was updated.
      //
      // Update the residual vectors: R_new = R_prev - A*P_prev*alpha
      //
      R_new->MvAddMv(one, *R_prev, zero, *R_prev);
      R_new->MvTimesMatAddMv(-one, *AP_prev, alpha, one);
      //
      // ****Compute the Current Relative Residual Norms and the Block Error****
      //
      R_new->MvNorm(_cur_resid_norms);
      //
      prev_ind_blksz = ind_blksz; // Save old ind_blksz of P_prev
      //
      // Update the number of current residuals that correspond
      // to linearly independent direction vectors. Note that
      // ind_idx are a subset of cur_idx.
      //
      k = 0;
      for (i=0; i< ind_blksz; i++){
	if (_cur_resid_norms[ ind_idx[i] ] / _init_resid_norms[ ind_idx[i] ] > _prec){
	  ind_idx[k] = ind_idx[i]; k = k+1;
	}
      }
      ind_blksz = k;
      //
      // ****************Print iteration information*****************************
      //
      if (_om.doOutput( 2 )) {
	PrintCGIterInfo( ind_idx, ind_blksz );
      }
      //
      // ****************Test for breakdown*************************************
      //
      if (ind_blksz <= 0){
	if (_om.doOutput( 0 )) {
	  _os << " Exiting Block CG iteration " 
	       << " -- Iteration# " << _iter << endl;
	  _os << " Reason: No more independent direction vectors" << endl;
	  _os << " Solution will be updated upon exiting loop" << endl;
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
	index2[i] = _blocksize + ind_idx[i];
      }
      delete R_new; R_new=0;
      
      if (new_blk) {
	R_new = _residvecs->CloneView( index2, ind_blksz );
	P_new = _basisvecs->CloneView( index2, ind_blksz );
      }
      else {
	R_new = _residvecs->CloneView( ind_idx, ind_blksz );
	P_new = _basisvecs->CloneView( ind_idx, ind_blksz );
      }
      //
      // Put the current preconditioned initial residual into P_new since P_new = precond_resid + P_prev * beta
      //
      if (_lp.ApplyLeftPrec( *R_new, *P_new ) != Ok ) { P_new->MvAddMv( one, *R_new, zero, *R_new ); }
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
      P_new->MvTransMv(-one, *AP_prev, beta);
      // 3)
      lapack.POTRS(UPLO, prev_ind_blksz, ind_blksz, T2.values(), prev_ind_blksz, beta.values(), prev_ind_blksz, &info);
      // Note: Solution returned in beta
      if (info != 0) {
	if (_om.doOutput( 0 )) {
	  _os << " Exiting Block CG iteration " 
	       << " -- Iteration# " << _iter << endl;
	  _os << "Reason: Cannot compute coefficient matrix beta" << endl;
	  _os << "Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Compute: P_new = precond_resid + P_prev * beta
      // ( remember P_new already = precond_resid )
      P_new->MvTimesMatAddMv(one, *P_prev, beta, one);
      //
      // Check A-orthogonality of new and previous blocks of direction vectors
      //
      if (_om.doOutput( 2 )) {
	_os << "Orthogonality check" << endl;
	CheckCGOrth(*P_prev, *P_new);   
      }
      //
      // Compute orthonormal block of direction vectors,
      // and check for dependencies, adjusting indices of
      // independent vectors if needed
      //
      Teuchos::SerialDenseMatrix<int,TYPE> G(ind_blksz,ind_blksz);
      exit_flg = QRFactorDef(*P_new, G, cols, num_ind);
      //
      // Check if the orthogonalization has failed.
      //
      if ( exit_flg ) {
        ind_blksz = num_ind;
	if (_om.doOutput( 0 )) {
	  _os  << " Exiting Block CG iteration "  
		<< " -- Iteration# " << _iter << endl;
	  _os << "Reason: No more linearly independent direction vectors" << endl;
	  _os << " Solution will be updated upon exiting loop" << endl;
	}
	break;
      }
      //
      // Check if the new block of direction vectors are linearly dependent
      //
      if (num_ind < ind_blksz) {
	if (_om.doOutput( 2 )) {
	  _os << "The new block of direction vectors are dependent " << endl;
	  _os << "# independent direction vectors: " << num_ind << endl;
	  _os << " Independent indices: " << endl;
	  for (i=0; i<num_ind ; i++)
	    _os << cols[i] << " ";
	  _os << endl << endl;
	}
	ind_blksz = num_ind;
        for (i=0; i<ind_blksz; i++)
	  ind_idx[i] = ind_idx[cols[i]];
      }
      //
      // Check A-orthogonality after orthonormalization
      //
      if (_om.doOutput( 2 )) {
	_os << "Orthogonality check after orthonormalization " << endl; 
	CheckCGOrth(*P_prev, *P_new);
      }
      // 
      // *****Update index of new blocks*****************************************
      //
      if (new_blk)
	new_blk--;
      else
	new_blk++;
      //
    } // end of the main CG loop -- for(_iter = 0;...)
    // *******************************************************************************
    //
    // Inform the linear problem manager that we are done with the current block of linear systems.
    //
    _lp.SetCurrLSVec();
    _lp.SolutionUpdated();
    //
    // Get the next block of linear systems, if it returns the null pointer we are done.
    //
    _cur_block_sol = _lp.GetCurrLHSVec();
    _cur_block_rhs = _lp.GetCurrRHSVec();
    //
    // Print out solver status.
    //
    if (_om.doOutput( 0 )) {
      _stest.Print(_os);
    }
    //
    // **************Free heap space**************
    //   
    delete temp_blk; temp_blk = 0;
    delete _basisvecs; _basisvecs = 0;
    delete _residvecs; _residvecs = 0;
    delete [] index; index=0;
    delete [] index2; index2=0;
    delete [] ind_idx; ind_idx = 0;
    delete [] cols; cols = 0;
    delete [] _cur_resid_norms; _cur_resid_norms = 0;
    delete [] _init_resid_norms; _init_resid_norms = 0;
    if (AP_prev) { delete AP_prev; AP_prev=0; }
    if (P_prev) { delete P_prev; P_prev=0;}
    if (P_new) { delete P_new; P_new=0; }
    if (R_prev) { delete R_prev; R_prev=0;}
    if (R_new) { delete R_new; R_new=0;}
    //
  } // end while ( _cur_block_sol && _cur_block_rhs )
  // **********************************************************************************
  //
} // end CGSolve()
//

template<class TYPE>
bool BlockCG<TYPE>::QRFactorDef (MultiVec<TYPE>& VecIn, 
				 Teuchos::SerialDenseMatrix<int,TYPE>& FouierR, 
				 int cols[], int &num) 
{
  int i, j, k;
  int num_orth, num_dep = 0;
  int nb = VecIn.GetNumberVecs();
  int *index = new int[ nb ]; assert(index!=NULL);
  int *dep_idx = new int[ nb ]; assert(dep_idx!=NULL);
  MultiVec<TYPE> *qj = 0, *Qj = 0;
  Teuchos::SerialDenseVector<int,TYPE> rj;
  const int IntOne = Teuchos::OrdinalTraits<int>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  TYPE * NormVecIn = new TYPE[nb]; assert(NormVecIn!=NULL);
  //
  // Set the index vector.
  //
  for ( i=0; i<nb; i++ ) { index[i] = i; }
  //
  // Zero out the array that will contain the Fourier coefficients.
  //
  FouierR.putScalar();
  //
  // Compute the Frobenius norm of VecIn -- this will be used to
  // determine rank deficiency of VecIn
  //
  VecIn.MvNorm(NormVecIn);
  TYPE FroNorm = zero;
  for (j=0; j<nb; j++) {
    FroNorm += NormVecIn[j] * NormVecIn[j];
  }
  FroNorm = sqrt(FroNorm);
  num = 0; 
  //
  // Start the loop to orthogonalize the nb columns of VecIn.
  //
  for ( j=0; j<nb; j++ ) {
    //
    // Grab the j-th column of VecIn (the first column is indexed to 
    // be the zero-th one).
    //
    qj = VecIn.CloneView(index+j, IntOne); assert(qj!=NULL);
    if ( j ) {
      //
      // Grab the first j columns of VecIn (that are now an orthogonal
      // basis for first j columns of the entering VecIn).
      //
      Qj = VecIn.CloneView(index, j);
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
	qj->MvTransMv( one, *Qj, rj );
	//
	// Sum result[0:j-1] into column j of R.
	//
	for (k=0; k<num_dep; k++) {
	  rj[dep_idx[k]] = zero;
	}
	//
	for ( k=0; k<j; k++ ) {
	  FouierR(k,j) += rj[k];
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
    qj->MvNorm ( &FouierR(j,j) );
    //
    if ( NormVecIn[j] > _prec && FouierR(j,j) > (_prec * NormVecIn[j]) ) {
      //
      // Normalize qj to make it into a unit vector.
      //
      TYPE rjj = one / FouierR(j,j);
      qj->MvAddMv ( rjj, *qj, zero, *qj );
      cols[num] = j;
      num++;
    }
    else {
      // 
      if (_om.doOutput( 2 )){
	_os << "Rank deficiency at column index: " << j << endl;
      }
      //
      // Don't normalize qj, enter one on diagonal of R,
      // and zeros in the row to the right of the diagonal -- this
      // requires updating the indices of the dependent columns
      //
      FouierR(j,j) = one;
      dep_idx[num_dep] = j;
      num_dep++;
    }	
    delete qj; qj = 0;
    delete Qj; Qj = 0;
  }
  delete [] index;
  delete [] dep_idx;
  delete [] NormVecIn;
  //
  // Return true if we could not create any independent direction vectors (failure).
  //
  if (!num) return true;
  return false;
  //
} // end QRFactorDef


template<class TYPE>
void BlockCG<TYPE>::CheckCGOrth(MultiVec<TYPE>& P1, MultiVec<TYPE>& P2) 
{
  //
  // This routine computes P2^T * A * P1
  // Checks the orthogonality wrt A between any two blocks of multivectors with the same length
  //
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  int i, k;
  //
  int numvecs1 = P1.GetNumberVecs();
  int numvecs2 = P2.GetNumberVecs();
  //
  MultiVec<TYPE>* AP = P1.CloneCopy();
  assert(AP!=NULL);
  _lp.ApplyOp(P1, *AP);
  //
  Teuchos::SerialDenseMatrix<int,TYPE> PAP(numvecs2, numvecs1);
  AP->MvTransMv(one, P2, PAP);
  //
  TYPE* ptr = PAP.values();
  TYPE column_sum;
  //
  for (k=0; k<numvecs1; k++) {
    column_sum = zero;
    for (i=0; i<numvecs2; i++) {
      column_sum += ptr[i];
    }
    _os << " P2^T*A*P1 for column "
	 << k << " is  " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
    ptr += numvecs2;
  }
  _os << " " << endl;
  //
  if(AP) { delete AP; AP=0; }  
  //  
} // end check_orthog
//


template<class TYPE>
void BlockCG<TYPE>::PrintCGIterInfo( int ind[], const int indsz )
{
  //
  int i;
  _os << "# of independent direction vectors: " << indsz << endl;    
  _os << " Independent indices: " << endl;
  for (i=0; i<indsz; i++){
    _os << ind[i] << " ";
  }
  _os << endl << endl;
  //
} // end Print_CGiter_info
//
//
} // end namespace Belos

#endif
// End of file BelosBlockCG.hpp


