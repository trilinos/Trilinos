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
#include "BelosMultiVec.hpp"
#include "BelosStatusTest.hpp"

/*!	\class Belos::CG

	\brief This class implements the preconditioned Conjugate Gradient algorithm for
	solving real symmetric positive definite linear systems of equations
	Ax = b, where b is the right-hand side vector and x is the corresponding solution.

	\author Heidi Thornquist
*/

namespace Belos {

template <class TYPE>
class CG : public IterativeSolver<TYPE> { 
public:
  //@{ \name Constructor/Destructor.
  //! %Belos::CG constructor.
  CG(LinearProblemManager<TYPE>& lp, StatusTest<TYPE>& stest, OutputManager<TYPE>& om);
  
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
  MultiVec<TYPE>* GetNativeResiduals( TYPE *normvec ) const;
  
  //! Get the actual residual vector for the current linear system.
  /*! This may force the solver to compute a current residual for its linear
  	system.  For CG, this method is not useful since the linear problem
	manager always has the current solution.
  */
  MultiVec<TYPE>* GetCurrentSoln() { return _cur_block_sol->CloneCopy(); };
  
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

  //! Pointer to block of the current residual vectors.
  MultiVec<TYPE> *_residvec;

  //! Current blocksize, iteration number, and basis pointer.
  int _iter;

  //! Numerical breakdown tolerances.
  TYPE _prec, _dep_tol;

  //! Output stream.
  ostream& _os;
};

//
// Implementation
//

template <class TYPE>
CG<TYPE>::CG(LinearProblemManager<TYPE>& lp,
		       StatusTest<TYPE>& stest,
		       OutputManager<TYPE>& om) : 
  _lp(lp), 
  _stest(stest),
  _om(om),
  _cur_block_rhs(0),
  _cur_block_sol(0),
  _residvec(0),
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
CG<TYPE>::~CG() 
{
 if (_residvec) delete _residvec; 
}

template <class TYPE>
void CG<TYPE>::SetCGBlkTols() 
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
MultiVec<TYPE>* CG<TYPE>::GetNativeResiduals( TYPE *normvec ) const 
{
  int i;
  int* index = new int[ 1 ];
  index[ 0 ] = 0;
  MultiVec<TYPE>* ResidMV = _residvec->CloneView( index, 1 );
  delete [] index;
  return ResidMV;
}

template <class TYPE>
void CG<TYPE>::Solve () 
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
  MultiVec<TYPE> *_p=0, *_Ap=0, *_z=0;
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
    // Only continue if the linear system is single-vector.
    //
    if ( _lp.GetBlockSize() > 1 ) return;
    //
    if (_om.doOutput( 0 )) {
      _os << endl;
      _os << "===================================================" << endl;
      _os << "Solving linear system(s):  " << _lp.GetRHSIndex() << " through " << _lp.GetRHSIndex()+_lp.GetNumToSolve() << endl;
      _os << endl;
    }	
    //
    _residvec = _cur_block_sol->Clone( 1 ); assert(_residvec!=NULL);
    _p = _cur_block_sol->Clone( 1 ); assert(_p!=NULL);
    _Ap = _cur_block_sol->Clone( 1 ); assert(_Ap!=NULL);
    _z = _cur_block_sol->Clone( 1 ); assert(_z!=NULL);
    //
    // ************ Compute the initial residual ********************************
    //
    // p0 = r0 = cur_block_rhs - A * cur_block_sol 
    //
    _p->MvAddMv(one, *_cur_block_sol, zero, *_cur_block_sol);
    //
    // Multiply the current solution by A and store in _Ap
    //       _Ap = A*_p 
    //
    _lp.ApplyOp( *_p, *_Ap );
    //
    // Compute initial residual and store in _residvec
    //     _residvec = cur_block_rhs - _Ap
    //
    _residvec->MvAddMv(one, *_cur_block_rhs, -one, *_Ap);
    _p->MvAddMv( one, *_residvec, zero, *_residvec );
    //
    //----------------Compute initial direction vectors--------------------------
    // Initially, they are set to the preconditioned residuals
    //
    if (_lp.ApplyLeftPrec( *_residvec, *_z ) != Ok ) { _z->MvAddMv( one , *_residvec, zero, *_residvec); }
    //
    // ***************************************************************************
    // ************************Main CG Loop***************************************
    // ***************************************************************************
    // 
    if (_om.doOutput( 2 )) _os << "Entering main CG loop" << endl << endl;
    //
    for (_iter=0; _stest.CheckStatus(this) == Unconverged && !exit_flg; _iter++) 
    {
      //
      // Multiply the current direction vector by A and store in _Ap
      //       _Ap = A*_p 
      //
      _lp.ApplyOp( *_p, *_Ap );  
      //
      // Compute alpha := <_residvec, _z> / <_p, _Ap >
      //
      _p->MvTransMv( one, *_Ap, temp_sdm );      
      _residvec->MvTransMv( one/temp_sdm(0,0), *_z, alpha );
      //
      // Check that alpha is a positive number!
      //
      if ( alpha(0,0) <= zero ) {
	if (_om.doOutput( 0 )) {
	  _os << " Exiting CG iteration " << endl;
	  _os << " Reason: Non-positive value for p^T*A*p ("<< alpha(0,0) <<") !!! "<< endl;
	}
	break; // Get out from this solve.
      }
      //
      // Update the solution vector x := x + alpha * _p
      //
      _cur_block_sol->MvAddMv( one, *_cur_block_sol, alpha(0,0), *_p );
      _lp.SolutionUpdated();
      //
      // Compute the denominator of beta before residual is updated [ old <_residvec, _z> ]
      //
      _residvec->MvTransMv( one, *_z, temp_sdm );
      //
      // Compute the new residual _residvec := _residvec - alpha * _Ap
      //
      _residvec->MvAddMv( one, *_residvec, -alpha(0,0), *_Ap );
      //
      // Compute beta := [ new <_residvec, _z> ] / [ old <_residvec, _z> ], and the new direction vector
      //
      if (_lp.ApplyLeftPrec( *_residvec, *_z ) != Ok ) { _z->MvAddMv( one, *_residvec, zero, *_residvec ); }
      //
      _residvec->MvTransMv( one/temp_sdm(0,0), *_z, beta );
      //
      _p->MvAddMv( one, *_residvec, beta(0,0), *_p );
      //
    } // end of the main CG loop -- for(_iter = 0;...)
    // *******************************************************************************
    //
    // Inform the linear problem manager that we are done with the current block of linear systems.
    //
    _lp.SetCurrLSVec();
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
    if (_residvec) { delete _residvec; _residvec=0; }
    if (_p) { delete _p; _p=0; }
    if (_Ap) { delete _Ap; _Ap=0;}
    if (_z) { delete _z; _z=0; }
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


