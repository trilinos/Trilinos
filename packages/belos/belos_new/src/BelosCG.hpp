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
// This file contains an implementation of the CG algorithm
// for solving real symmetric positive definite linear systems of 
// equations Ax = b, where b is a single-vector and x is the corresponding solution.
// This implementation uses the preconditioned inner-product to retain symmetry of
// the linear system.
//
#ifndef BELOS_CG_HPP
#define BELOS_CG_HPP

/*!
  \file BelosCG.hpp

  \brief Belos concrete class for solving symmetric positive definite linear systems with the
  preconditioned Conjugate Gradient (CG) method.
*/

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOperator.hpp"
#include "BelosStatusTest.hpp"
#include "BelosMultiVecTraits.hpp"

/*!	\class Belos::CG

	\brief This class implements the preconditioned Conjugate Gradient algorithm for
	solving real symmetric positive definite linear systems of equations
	Ax = b, where b is the right-hand side vector and x is the corresponding solution.

	\author Heidi Thornquist
*/

namespace Belos {

template <class TYPE, class OP, class MV>
class CG : public IterativeSolver<TYPE,OP,MV> { 
public:
  //@{ \name Constructor/Destructor.
  //! %Belos::CG constructor.
  CG(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp, 
     const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest, 
     const RefCountPtr<OutputManager<TYPE> > &om);
  
  //! %CG destructor.
  virtual ~CG();
  //@}
  
  //@{ \name Accessor methods
  
  //! Get the iteration count for the current linear system.
  int GetNumIters() const { return( _iter ); }
  
  //! Get the restart count of the iteration method for the current linear system [not valid for CG].
  int GetNumRestarts() const { return(0); }
   
  //! Get the solvers native residual for the current linear system.
  /*! 
      \note The memory for the residual MultiVec must be handled by the calling routine.
   */
  RefCountPtr<const MV> GetNativeResiduals( TYPE *normvec ) const;
  
  //! Get the actual residual vector for the current linear system.
  /*! This may force the solver to compute a current residual for its linear
  	system.  For CG, this method is not useful since the linear problem
	manager always has the current solution.
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
  RefCountPtr<MV> _residvec;

  //! Output stream.
  ostream* _os;

  //! Current blocksize, iteration number, and basis pointer.
  int _iter;

  //! Numerical breakdown tolerances.
  TYPE _prec, _dep_tol;

  typedef MultiVecTraits<TYPE,MV> MVT;
};

//
// Implementation
//

template <class TYPE, class OP, class MV>
CG<TYPE,OP,MV>::CG(const RefCountPtr<LinearProblemManager<TYPE,OP,MV> > &lp,
		   const RefCountPtr<StatusTest<TYPE,OP,MV> > &stest,
		   const RefCountPtr<OutputManager<TYPE> > &om) : 
  _lp(lp), 
  _stest(stest),
  _om(om),
  _os(&om->GetOStream()),
  _iter(0),
  _prec(5.0e-15), 
  _dep_tol(0.75)
{ 
  //
  // Set the block orthogonality tolerances
  //
  SetCGBlkTols();
}

template <class TYPE, class OP, class MV>
CG<TYPE,OP,MV>::~CG() 
{
}

template <class TYPE, class OP, class MV>
void CG<TYPE,OP,MV>::SetCGBlkTols() 
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
RefCountPtr<const MV> CG<TYPE,OP,MV>::GetNativeResiduals( TYPE *normvec ) const 
{
  int i;
  int* index = new int[ 1 ];
  index[ 0 ] = 0;
  RefCountPtr<MV> ResidMV = MVT::CloneView( *_residvec, index, 1 );
  delete [] index;
  return ResidMV;
}

template <class TYPE, class OP, class MV>
void CG<TYPE,OP,MV>::Solve () 
{
  //
  int i, j, k, info, num_ind;
  int ind_blksz, prev_ind_blksz;
  bool exit_flg = false;
  char UPLO = 'U';
  const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
  const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
  Teuchos::LAPACK<int,TYPE> lapack;
  Teuchos::SerialDenseMatrix<int,TYPE> alpha( 1, 1 );
  Teuchos::SerialDenseMatrix<int,TYPE> beta( 1, 1 );
  Teuchos::SerialDenseMatrix<int,TYPE> temp_sdm( 1, 1 );
  RefCountPtr<MV> _p, _Ap, _z;
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
    // Only continue if the linear system is single-vector.
    //
    if ( _lp->GetBlockSize() > 1 ) return;
    //
    if (_om->doOutput( 0 )) {
      *_os << endl;
      *_os << "===================================================" << endl;
      *_os << "Solving linear system(s):  " << _lp->GetRHSIndex() << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve() << endl;
      *_os << endl;
    }	
    //
    _residvec = MVT::Clone( *_cur_block_sol, 1 );
    _p = MVT::Clone( *_cur_block_sol, 1 );
    _Ap = MVT::Clone( *_cur_block_sol, 1 ); 
    _z = MVT::Clone( *_cur_block_sol, 1 ); 
    //
    // ************ Compute the initial residual ********************************
    //
    // p0 = r0 = cur_block_rhs - A * cur_block_sol 
    //
    MVT::MvAddMv(one, *_cur_block_sol, zero, *_cur_block_sol, *_p);
    //
    // Multiply the current solution by A and store in _Ap
    //       _Ap = A*_p 
    //
    _lp->ApplyOp( *_p, *_Ap );
    //
    // Compute initial residual and store in _residvec
    //     _residvec = cur_block_rhs - _Ap
    //
    MVT::MvAddMv(one, *_cur_block_rhs, -one, *_Ap, *_residvec);
    MVT::MvAddMv( one, *_residvec, zero, *_residvec, *_p );
    //
    //----------------Compute initial direction vectors--------------------------
    // Initially, they are set to the preconditioned residuals
    //
    if (_lp->ApplyLeftPrec( *_residvec, *_z ) != Ok ) { MVT::MvAddMv( one , *_residvec, zero, *_residvec, *_z); }
    //
    // ***************************************************************************
    // ************************Main CG Loop***************************************
    // ***************************************************************************
    // 
    if (_om->doOutput( 2 )) *_os << "Entering main CG loop" << endl << endl;
    //
    for (_iter=0; _stest->CheckStatus(this) == Unconverged && !exit_flg; _iter++) 
    {
      //
      // Multiply the current direction vector by A and store in _Ap
      //       _Ap = A*_p 
      //
      _lp->ApplyOp( *_p, *_Ap );  
      //
      // Compute alpha := <_residvec, _z> / <_p, _Ap >
      //
      MVT::MvTransMv( *_p, one, *_Ap, temp_sdm );
      MVT::MvTransMv( *_residvec, one/temp_sdm(0,0), *_z, alpha );
      //
      // Check that alpha is a positive number!
      //
      if ( alpha(0,0) <= zero ) {
	if (_om->doOutput( 0 )) {
	  *_os << " Exiting CG iteration " << endl;
	  *_os << " Reason: Non-positive value for p^T*A*p ("<< alpha(0,0) <<") !!! "<< endl;
	}
	break; // Get out from this solve.
      }
      //
      // Update the solution vector x := x + alpha * _p
      //
      MVT::MvAddMv( one, *_cur_block_sol, alpha(0,0), *_p, *_cur_block_sol );
      _lp->SolutionUpdated();
      //
      // Compute the denominator of beta before residual is updated [ old <_residvec, _z> ]
      //
      MVT::MvTransMv( *_residvec, one, *_z, temp_sdm );
      //
      // Compute the new residual _residvec := _residvec - alpha * _Ap
      //
      MVT::MvAddMv( one, *_residvec, -alpha(0,0), *_Ap, *_residvec );
      //
      // Compute beta := [ new <_residvec, _z> ] / [ old <_residvec, _z> ], and the new direction vector
      //
      if (_lp->ApplyLeftPrec( *_residvec, *_z ) != Ok ) { MVT::MvAddMv( one, *_residvec, zero, *_residvec, *_z); }
      //
      MVT::MvTransMv( *_residvec, one/temp_sdm(0,0), *_z, beta );
      //
      MVT::MvAddMv( one, *_residvec, beta(0,0), *_p, *_p );
      //
    } // end of the main CG loop -- for(_iter = 0;...)
    // *******************************************************************************
    //
    // Inform the linear problem manager that we are done with the current block of linear systems.
    //
    _lp->SetCurrLSVec();
    //
    // Get the next block of linear systems, if it returns the null pointer we are done.
    //
    _cur_block_sol = _lp->GetCurrLHSVec();
    _cur_block_rhs = _lp->GetCurrRHSVec();
    //
    // Print out solver status.
    //
    if (_om->doOutput( 0 )) {
      _stest->Print(*_os);
    }
    //
  } // end while ( _cur_block_sol && _cur_block_rhs )
  // **********************************************************************************
  //
} // end CGSolve()
//
} // namespace Belos
//
#endif // BELOS_CG_HPP
//
// End of file BelosCG.hpp


