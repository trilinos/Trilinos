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
#include "BelosStatusTest.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosCG.hpp"

/*!	\class Belos::BlockCG

	\brief This class implements the preconditioned Conjugate Gradient algorithm for
	solving real symmetric positive definite linear systems of equations
	AX = B, where B is a matrix containing one or more right-hand sides.

	\author Teri Barth and Heidi Thornquist
*/

namespace Belos {

template <class TYPE, class OP, class MV>
class BlockCG : public IterativeSolver<TYPE,OP,MV> { 
public:
  //@{ \name Constructor/Destructor.
  //! %Belos::BlockCG constructor.
  BlockCG(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp,
	  const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest,
	  const RefCountPtr<OutputManager<TYPE> > &om);
  
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
      \note The memory for the residual MultiVec must be handled by the calling routine.
   */
  RefCountPtr<const MV> GetNativeResiduals( TYPE *normvec ) const;
  
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
  LinearProblemManager<TYPE,OP,MV>& GetLinearProblem() const { return( *_lp ); }

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
  bool QRFactorDef(MV&, Teuchos::SerialDenseMatrix<int,TYPE>&,
		   int[], int&);
  void CheckCGOrth(MV&, MV&);
  void PrintCGIterInfo(int[], const int);
  void BlockIteration();

  //! Linear problem manager [ must be passed in by the user ]
  RefCountPtr<LinearProblemManager<TYPE,OP,MV> > _lp; 

  //! Status test [ must be passed in by the user ]
  RefCountPtr<StatusTest<TYPE,OP,MV> > _stest; 

  //! Output manager [ must be passed in by the user ]
  RefCountPtr<OutputManager<TYPE> > _om;

  //! Pointer to current linear systems block of solution vectors [obtained from linear problem manager]
  RefCountPtr<MV> _cur_block_sol;

  //! Pointer to current linear systems block of right-hand sides [obtained from linear problem manager]
  RefCountPtr<MV> _cur_block_rhs; 

  //! Pointer to block of the current residual vectors.
  RefCountPtr<MV> _residvecs;

  //! Output stream.
  ostream* _os;

  //! Current blocksize, iteration number, and basis pointer.
  int _blocksize, _iter, _new_blk;

  //! Numerical breakdown tolerances.
  TYPE _prec, _dep_tol;

  typedef MultiVecTraits<TYPE,MV> MVT;
};

//
// Implementation
//

template <class TYPE, class OP, class MV>
BlockCG<TYPE,OP,MV>::BlockCG(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp,
		       const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest,
		       const RefCountPtr<OutputManager<TYPE> > &om) : 
  _lp(lp), 
  _stest(stest),
  _om(om),
  _os(&om->GetOStream()),
  _blocksize(0), 
  _iter(0),
  _new_blk(1),
  _prec(5.0e-15), 
  _dep_tol(0.75)
{ 
  //
  // Set the block orthogonality tolerances
  //
  SetCGBlkTols();
}

template <class TYPE, class OP, class MV>
BlockCG<TYPE,OP,MV>::~BlockCG() 
{
}

template <class TYPE, class OP, class MV>
void BlockCG<TYPE,OP,MV>::SetCGBlkTols() 
{
  const TYPE two = 2.0;
  TYPE eps;
  char precision = 'P';
  Teuchos::LAPACK<int,TYPE> lapack;
  eps = lapack.LAMCH(precision);
  _prec = eps;
  _dep_tol = 1/sqrt(two);
}

template <class TYPE, class OP, class MV>
RefCountPtr<const MV> BlockCG<TYPE,OP,MV>::GetNativeResiduals( TYPE *normvec ) const 
{
  int i;
  int* index = new int[ _blocksize ];
  if (_new_blk)
    for (i=0; i<_blocksize; i++) { index[i] = i; }
  else
    for (i=0; i<_blocksize; i++) { index[i] = _blocksize + i; }
  RefCountPtr<MV> ResidMV = MVT::CloneView( *_residvecs, index, _blocksize );
  delete [] index;
  return ResidMV;
}

template <class TYPE, class OP, class MV>
void BlockCG<TYPE,OP,MV>::Solve () 
{
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
    _blocksize = _lp->GetBlockSize();
    if (_blocksize == 1 ) {
      //
      // Create single vector CG solver for this linear system.
      //
      Belos::CG<TYPE,OP,MV> CGSolver(_lp, _stest, _om);
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
} // end Solve()
//


template <class TYPE, class OP, class MV>
void BlockCG<TYPE,OP,MV>::BlockIteration ( ) 
{
  //
  int i, j, k, info, num_ind;
  int ind_blksz, prev_ind_blksz;
  bool exit_flg = false;
  char UPLO = 'U';
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  Teuchos::LAPACK<int,TYPE> lapack;
  RefCountPtr<MV> P_prev, R_prev, AP_prev, P_new, R_new;
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
  RefCountPtr<MV> _basisvecs = MVT::Clone( *_cur_block_sol, 2*_blocksize ); 
  RefCountPtr<MV> temp_blk = MVT::Clone( *_cur_block_sol, _blocksize );
  TYPE* _cur_resid_norms = new TYPE[_blocksize]; assert(_cur_resid_norms!=NULL);
  TYPE* _init_resid_norms= new TYPE[_blocksize]; assert(_init_resid_norms!=NULL);
  _residvecs = MVT::Clone( *_cur_block_sol, 2*_blocksize ); 
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
  if (_om->doOutput( 0 )) {
    *_os << endl;
    *_os << "===================================================" << endl;
    *_os << "Solving linear system(s):  " << _lp->GetRHSIndex() << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve() << endl;
    *_os << endl;
  }	
  //
  //
  // ************ Compute the initial residuals ********************************
  //
  // Associate the first block of _basisvecs with P_prev and the
  // first block of _residvecs with R_prev
  //
  P_prev = MVT::CloneView( *_basisvecs, ind_idx, _blocksize);
  R_prev = MVT::CloneView( *_residvecs, ind_idx, _blocksize);
  AP_prev = MVT::CloneView( *temp_blk, ind_idx, _blocksize);
  //
  // Store initial guesses to AX = B in 1st block of _basisvecs
  //         P_prev = one*cur_block_sol + zero*P_prev
  //
  MVT::MvAddMv(one, *_cur_block_sol, zero, *_cur_block_sol, *P_prev);
  //
  // Multiply by A and store in AP_prev
  //       AP_prev = A*P_prev
  //
  _lp->ApplyOp( *P_prev, *AP_prev );
  //
  // Compute initial residual block and store in 1st block of _residvecs
  //     R_prev = cur_block_rhs - A*P_prev
  //
  MVT::MvAddMv(one, *_cur_block_rhs, -one, *AP_prev, *R_prev );
  //
  //-------Compute and save the initial residual norms----------
  //
  MVT::MvNorm(*R_prev, _init_resid_norms);
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
    R_prev = MVT::CloneView( *_residvecs, ind_idx, ind_blksz );
    P_prev = MVT::CloneView( *_basisvecs, ind_idx, ind_blksz );
    //
    //----------------Compute initial direction vectors--------------------------
    // Initially, they are set to the preconditioned residuals
    //
    if (_lp->ApplyLeftPrec( *R_prev, *P_prev ) != Ok ) { MVT::MvAddMv( one , *R_prev, zero, *R_prev, *P_prev); }
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
      if (_om->doOutput( 0 )) {
	*_os << " Exiting Block CG iteration " << endl; 
	*_os << " Reason: No independent initial direction vectors" << endl;
      }
    }		
    if (num_ind < ind_blksz) {
      // The initial block of direction vectors are linearly dependent
      if (_om->doOutput( 0 )) {
	*_os << " Initial direction vectors are dependent" << endl;
	*_os << " Adjusting blocks and indices for iteration" << endl;
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
    if (_om->doOutput( 0 )) {
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
  if (_om->doOutput( 2 )) *_os << "Entering main CG loop" << endl << endl;
  //
  _new_blk = 1;
  for (_iter=0; _stest->CheckStatus(this) == Unconverged && !exit_flg; _iter++) {
    //
    //----------------Compute the new blocks of iterates and residuals------------------
    //
    // Get views of the previous blocks of residuals, direction vectors, etc.
    //
    if (_new_blk){
      P_prev = MVT::CloneView( *_basisvecs, ind_idx, ind_blksz );
    }
    else {
      for (i=0; i< ind_blksz; i++) {
	index2[i] = _blocksize + ind_idx[i];
      } 
      P_prev = MVT::CloneView( *_basisvecs, index2, ind_blksz);
    }
    //
    for (i=0; i < _blocksize; i++){
      index2[i] = _blocksize + i;
    }
    if (_new_blk){
      R_prev = MVT::CloneView( *_residvecs, index, _blocksize );
      R_new = MVT::CloneView( *_residvecs, index2, _blocksize );
    }
    else {
      R_prev = MVT::CloneView( *_residvecs, index2, _blocksize );
      R_new = MVT::CloneView( *_residvecs, index, _blocksize );
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
      AP_prev = MVT::CloneView( *temp_blk, ind_idx, ind_blksz ); 
      alpha.reshape( ind_blksz, _blocksize );
      T2.reshape( ind_blksz, ind_blksz );
    }
    _lp->ApplyOp( *P_prev, *AP_prev );
    MVT::MvTransMv( *R_prev, one, *P_prev, alpha);   
    MVT::MvTransMv( *AP_prev, one, *P_prev, T2);
    //
    // 2)
    lapack.POTRF(UPLO, ind_blksz, T2.values(), ind_blksz, &info);
    if (info != 0) {
      if(_om->doOutput( 0 )){
	*_os << " Exiting Block CG iteration "
	    << " -- Iteration# " << _iter << endl;
	*_os << " Reason: Cannot compute coefficient matrix alpha" << endl;
	*_os << " P_prev'* A*P_prev is singular" << endl;
	*_os << " Solution will be updated upon exiting loop" << endl;
      }
      break;
    }
    // 3)
    lapack.POTRS(UPLO, ind_blksz, _blocksize, T2.values(), ind_blksz, alpha.values(), ind_blksz, &info);
    // Note: solution returned in alpha
    if (info != 0) {
      if(_om->doOutput( 0 )){
	*_os << " Exiting Block CG iteration "
	    << " -- Iteration# " << _iter << endl;
	*_os << " Reason: Cannot compute coefficient matrix alpha" << endl;
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
    MVT::MvNorm( *R_new, _cur_resid_norms);
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
    if (_om->doOutput( 2 )) {
      PrintCGIterInfo( ind_idx, ind_blksz );
    }
    //
    // ****************Test for breakdown*************************************
    //
    if (ind_blksz <= 0){
      if (_om->doOutput( 0 )) {
	*_os << " Exiting Block CG iteration " 
	    << " -- Iteration# " << _iter << endl;
	*_os << " Reason: No more independent direction vectors" << endl;
	*_os << " Solution will be updated upon exiting loop" << endl;
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
    
    if (_new_blk) {
      R_new = MVT::CloneView( *_residvecs, index2, ind_blksz );
      P_new = MVT::CloneView( *_basisvecs, index2, ind_blksz );
    }
    else {
      R_new = MVT::CloneView( *_residvecs, ind_idx, ind_blksz );
      P_new = MVT::CloneView( *_basisvecs, ind_idx, ind_blksz );
    }
    //
    // Put the current preconditioned initial residual into P_new since P_new = precond_resid + P_prev * beta
    //
    if (_lp->ApplyLeftPrec( *R_new, *P_new ) != Ok ) { MVT::MvAddMv( one, *R_new, zero, *R_new, *P_new ); }
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
    MVT::MvTransMv(*P_new, -one, *AP_prev, beta);
    // 3)
    lapack.POTRS(UPLO, prev_ind_blksz, ind_blksz, T2.values(), prev_ind_blksz, beta.values(), prev_ind_blksz, &info);
    // Note: Solution returned in beta
    if (info != 0) {
      if (_om->doOutput( 0 )) {
	*_os << " Exiting Block CG iteration " 
	    << " -- Iteration# " << _iter << endl;
	*_os << "Reason: Cannot compute coefficient matrix beta" << endl;
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
    if (_om->doOutput( 2 )) {
      *_os << "Orthogonality check" << endl;
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
      if (_om->doOutput( 0 )) {
	*_os  << " Exiting Block CG iteration "  
	     << " -- Iteration# " << _iter << endl;
	*_os << "Reason: No more linearly independent direction vectors" << endl;
	*_os << " Solution will be updated upon exiting loop" << endl;
      }
      break;
    }
    //
    // Check if the new block of direction vectors are linearly dependent
    //
    if (num_ind < ind_blksz) {
      if (_om->doOutput( 2 )) {
	*_os << "The new block of direction vectors are dependent " << endl;
	*_os << "# independent direction vectors: " << num_ind << endl;
	*_os << " Independent indices: " << endl;
	for (i=0; i<num_ind ; i++)
	  *_os << cols[i] << " ";
	*_os << endl << endl;
      }
      ind_blksz = num_ind;
      for (i=0; i<ind_blksz; i++)
	ind_idx[i] = ind_idx[cols[i]];
    }
    //
    // Check A-orthogonality after orthonormalization
    //
    if (_om->doOutput( 2 )) {
      *_os << "Orthogonality check after orthonormalization " << endl; 
      CheckCGOrth(*P_prev, *P_new);
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
  if (_om->doOutput( 0 )) {
    _stest->Print(*_os);
  }  
  // *******************************************************************************
  // **************Free heap space**************
  //   
  delete [] index; index=0;
  delete [] index2; index2=0;
  delete [] ind_idx; ind_idx = 0;
  delete [] cols; cols = 0;
  delete [] _cur_resid_norms; _cur_resid_norms = 0;
  delete [] _init_resid_norms; _init_resid_norms = 0;
  //
} // end BlockIteration()
//

template<class TYPE, class OP, class MV>
bool BlockCG<TYPE,OP,MV>::QRFactorDef (MV& VecIn, 
				 Teuchos::SerialDenseMatrix<int,TYPE>& FouierR, 
				 int cols[], int &num) 
{
  int i, j, k;
  int num_orth, num_dep = 0;
  int nb = MVT::GetNumberVecs( VecIn );
  int *index = new int[ nb ]; assert(index!=NULL);
  int *dep_idx = new int[ nb ]; assert(dep_idx!=NULL);
  RefCountPtr<MV> qj, Qj;
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
  MVT::MvNorm( VecIn, NormVecIn );
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
    qj = MVT::CloneView( VecIn, index+j, IntOne);
    if ( j ) {
      //
      // Grab the first j columns of VecIn (that are now an orthogonal
      // basis for first j columns of the entering VecIn).
      //
      Qj = MVT::CloneView( VecIn, index, j);
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
	MVT::MvTransMv( *qj, one, *Qj, rj );
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
	MVT::MvTimesMatAddMv(-one, *Qj, rj, one, *qj);
      }
    }
    //
    // Compute the norm of column j of VecIn (=qj).
    //
    MVT::MvNorm( *qj, &FouierR(j,j) );
    //
    if ( NormVecIn[j] > _prec && FouierR(j,j) > (_prec * NormVecIn[j]) ) {
      //
      // Normalize qj to make it into a unit vector.
      //
      TYPE rjj = one / FouierR(j,j);
      MVT::MvAddMv( rjj, *qj, zero, *qj, *qj );
      cols[num] = j;
      num++;
    }
    else {
      // 
      if (_om->doOutput( 2 )){
	*_os << "Rank deficiency at column index: " << j << endl;
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


template<class TYPE, class OP, class MV>
void BlockCG<TYPE,OP,MV>::CheckCGOrth(MV& P1, MV& P2) 
{
  //
  // This routine computes P2^T * A * P1
  // Checks the orthogonality wrt A between any two blocks of multivectors with the same length
  //
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  int i, k;
  //
  int numvecs1 = MVT::GetNumberVecs( P1 );
  int numvecs2 = MVT::GetNumberVecs( P2 );
  //
  RefCountPtr<MV> AP = MVT::CloneCopy( P1 );
  _lp->ApplyOp(P1, *AP);
  //
  Teuchos::SerialDenseMatrix<int,TYPE> PAP(numvecs2, numvecs1);
  MVT::MvTransMv(*AP, one, P2, PAP);
  //
  TYPE* ptr = PAP.values();
  TYPE column_sum;
  //
  for (k=0; k<numvecs1; k++) {
    column_sum = zero;
    for (i=0; i<numvecs2; i++) {
      column_sum += ptr[i];
    }
    *_os << " P2^T*A*P1 for column "
	 << k << " is  " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
    ptr += numvecs2;
  }
  *_os << " " << endl;
  //
} // end check_orthog
//


template<class TYPE, class OP, class MV>
void BlockCG<TYPE,OP,MV>::PrintCGIterInfo( int ind[], const int indsz )
{
  //
  int i;
  *_os << "# of independent direction vectors: " << indsz << endl;    
  *_os << " Independent indices: " << endl;
  for (i=0; i<indsz; i++){
    *_os << ind[i] << " ";
  }
  *_os << endl << endl;
  //
} // end Print_CGiter_info
//
//
} // end namespace Belos

#endif
// End of file BelosBlockCG.hpp


