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
// This file contains an implementation of the TFQMR algorithm
// for solving non-Hermitian linear systems of equations Ax = b, 
// where b is a single-vector and x is the corresponding solution.
//
// The implementation is a slight modification on the TFQMR algorithm
// found in Saad's "Iterative Methods for Sparse Linear Systems".
//

#ifndef BELOS_TFQMR_HPP
#define BELOS_TFQMR_HPP

/*!
  \file BelosTFQMR.hpp

  \brief Belos concrete class for solving non-Hermitian linear systems with the
  preconditioned tranpose-free QMR (TFQMR) method.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOperator.hpp"
#include "BelosStatusTest.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	\class Belos::TFQMR

	\brief This class implements the preconditioned transpose-free QMR algorithm for
	solving non-Hermitian linear systems of equations Ax = b, where b is the right-hand 
	side vector and x is the corresponding solution.

	\author Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType, class MV, class OP>
  class TFQMR : public IterativeSolver<ScalarType,MV,OP> { 
  public:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    
    //@{ \name Constructor/Destructor.

    //! %Belos::TFQMR constructor.
    TFQMR(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp, 
	  const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest, 
	  const RefCountPtr<OutputManager<ScalarType> > &om,
	  const RefCountPtr<ParameterList> &pl = Teuchos::null
	  );
    
    //! %Belos::TFQMR destructor.
    virtual ~TFQMR() {};
    //@}
    
    //@{ \name Accessor methods
    
    //! Get the iteration count for the current linear system.
    int GetNumIters() const { return( _iter ); }
    
    //! Get the restart count of the iteration method for the current linear system [not valid for TFQMR].
    int GetNumRestarts() const { return(0); }
    
    //! Get the solvers native residual for the current linear system.
    /*! 
      \note The memory for the residual MultiVec must be handled by the calling routine.
    */
    RefCountPtr<const MV> GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const;
    
    //! Get the actual residual vector for the current linear system.
    /*! This may force the solver to compute a current residual for its linear
      system.  For TFQMR, this method is not useful since the linear problem
      manager always has the current solution.
    */
    RefCountPtr<MV> GetCurrentSoln() { return MVT::CloneCopy( *_cur_block_sol ); };
  
    //! Get a constant reference to the current linear problem.  
    /*! This may include a current solution, if the solver has recently restarted or completed.
     */
    RefCountPtr<LinearProblem<ScalarType,MV,OP> > GetLinearProblem() const { return( _lp ); }
    
    RefCountPtr<StatusTest<ScalarType,MV,OP> > GetStatusTest() const { return( _stest ); }

    //@} 
    
    //@{ \name Solver application method.
    
    /*! \brief This method uses the iterative method to compute approximate solutions
      to the original problem.  This method can return unconverged if the maximum number
      of iterations is reached, or numerical breakdown is observed.
    */
    void Solve();
    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{
    
    /** \brief Method to return description of the block GMRES solver */
    std::string description() const;
    
    //@}
    
  private:
    
    //! Linear problem manager [ must be passed in by the user ]
    RefCountPtr<LinearProblem<ScalarType,MV,OP> > _lp; 
    
    //! Status test [ must be passed in by the user ]
    RefCountPtr<StatusTest<ScalarType,MV,OP> > _stest; 
    
    //! Output manager [ must be passed in by the user ]
    RefCountPtr<OutputManager<ScalarType> > _om;
    
    //! Parameter list containing information for configuring the linear solver. [ must be passed in by the user ]
    RefCountPtr<ParameterList> _pl;     

    //! Pointer to current linear systems block of solution vectors [obtained from linear problem manager]
    RefCountPtr<MV> _cur_block_sol;
    
    //! Pointer to current linear systems block of right-hand sides [obtained from linear problem manager]
    RefCountPtr<MV> _cur_block_rhs; 
    
    //! Pointer to block of the current residual vectors.
    RefCountPtr<MV> _residvec;
    
    //! Output stream.
    RefCountPtr<ostream> _os;
    
    //! Current iteration number.
    int _iter;

    //! Current residual norm estimate.
    std::vector<MagnitudeType> _tau;

    //! Restart timers each time Solve() is called.
    bool _restartTimers;

    //! Internal timers
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerTotal;
  };
  
  //
  // Implementation
  //
  
  template <class ScalarType, class MV, class OP>
  TFQMR<ScalarType,MV,OP>::TFQMR(const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp,
				 const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest,
				 const RefCountPtr<OutputManager<ScalarType> > &om,
				 const RefCountPtr<ParameterList> &pl 
				 ) : 
    _lp(lp), 
    _stest(stest),
    _om(om),
    _pl(pl),
    _os(om->GetOStream()),
    _iter(0),
    _tau(1),
    _restartTimers(true),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time"))
  { 
  }
  
  template <class ScalarType, class MV, class OP>
  RefCountPtr<const MV> 
  TFQMR<ScalarType,MV,OP>::GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
  {
    MagnitudeType one = Teuchos::ScalarTraits<MagnitudeType>::one();
    if (normvec)
      (*normvec)[0] = Teuchos::ScalarTraits<MagnitudeType>::squareroot( _iter + one )*_tau[0];

    return Teuchos::null;
  }
  
  template <class ScalarType, class MV, class OP>
  void 
  TFQMR<ScalarType,MV,OP>::Solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);
    
    if ( _restartTimers ) {
      _timerOp->reset();
    }
    
    // Get the output stream from the OutputManager
    _os = _om->GetOStream();
    //
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType STzero = Teuchos::ScalarTraits<ScalarType>::zero();
    const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();
    bool exit_flg = false;
    Teuchos::SerialDenseMatrix<int,ScalarType> _alpha( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> _rho( 1, 1 ), _rho_old( 1, 1 );
    RefCountPtr<MV> _v, _w, _u, _Au, _d, _rtilde;
    ScalarType _eta = STzero, _beta = STzero;
    std::vector<MagnitudeType> _cs(1,MTzero), _theta(1,MTzero);
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
      if (_om->isVerbosityAndPrint( IterationDetails )) {
	*_os << endl;
	*_os << "===================================================" << endl;
	*_os << "Solving linear system(s):  " << _lp->GetRHSIndex() << " through " << _lp->GetRHSIndex()+_lp->GetNumToSolve() << endl;
	*_os << endl;
      }	
      //
      _d = MVT::Clone( *_cur_block_sol, 1 );
      _v = MVT::Clone( *_cur_block_sol, 1 );
      MVT::MvInit( *_d );
      //
      // ************ Compute the initial residual ********************************
      //
      // w0 = y0 = r0 = cur_block_rhs - A * cur_block_sol 
      //
      _residvec = MVT::Clone( *_cur_block_sol, 1 );
      _lp->ComputeResVec( &*_residvec, &*_cur_block_sol, &*_cur_block_rhs );
      //
      _w = MVT::CloneCopy( *_residvec ); 
      _u = MVT::CloneCopy( *_residvec ); 
      _rtilde = MVT::CloneCopy( *_residvec ); 
      //
      // Multiply the current residual by Op and store in _v
      //       _v = _Op*_residvec 
      //
      {
	Teuchos::TimeMonitor OpTimer(*_timerOp);
	_lp->Apply( *_residvec, *_v );
      }
      _Au = MVT::CloneCopy( *_v ); 
      //
      // Compute initial scalars: theta, eta, tau, rho_old
      //
      MVT::MvNorm( *_residvec, &_tau );                         // tau = ||r_0||
      MVT::MvTransMv( one, *_residvec, *_rtilde, _rho_old );    // rho = (r_0, r_tilde)
      //
      // ***************************************************************************
      // ************************Main TFQMR Loop***************************************
      // ***************************************************************************
      // 
      for (_iter=0; _stest->CheckStatus(this) == Unconverged && !exit_flg; _iter++) 
	{
	  //
	  //--------------------------------------------------------
	  // Compute the new alpha if we need to
	  //--------------------------------------------------------
	  //
	  if (_iter%2 == 0) {
	    MVT::MvTransMv( one, *_v, *_rtilde, _alpha );      //   alpha = rho / (v, r_tilde) 
	    _alpha(0,0) = _rho_old(0,0)/_alpha(0,0);
	  }
	  //
	  //--------------------------------------------------------
	  // Update d.
	  //   d = u + (theta^2/alpha)eta*d
	  //--------------------------------------------------------
	  //
	  MVT::MvAddMv( one, *_u, (_theta[0]*_theta[0]/_alpha(0,0))*_eta, *_d, *_d );
	  //
	  //--------------------------------------------------------
	  // Update w.
	  //   w = w - alpha*Au
	  //--------------------------------------------------------
	  //
	  MVT::MvAddMv( one, *_w, -_alpha(0,0), *_Au, *_w );
	  //
	  //--------------------------------------------------------
	  // Update u if we need to.
	  //   u = u - alpha*v
	  //   
	  // Note: This is usually computed with alpha (above), but we're trying be memory efficient.
	  //--------------------------------------------------------
	  //
	  if (_iter%2 == 0) {
	    MVT::MvAddMv( one, *_u, -_alpha(0,0), *_v, *_u );
	  }
	  //
	  //--------------------------------------------------------
	  // Compute the new theta, c, eta, tau; i.e. the update to the least squares solution.
	  //--------------------------------------------------------
	  //
	  MVT::MvNorm( *_w, &_theta );     // theta = ||w|| / tau
	  _theta[0] /= _tau[0];

	  // cs = sqrt( 1 + theta^2 )
	  _cs[0] = Teuchos::ScalarTraits<ScalarType>::squareroot(one + _theta[0]*_theta[0]);

	  _tau[0] *= _theta[0]*_cs[0];     // tau = tau * theta * cs
	  _eta = _cs[0]*_cs[0]*_alpha(0,0);     // eta = cs^2 * alpha

	  //
	  //--------------------------------------------------------
	  // Update the solution.
	  //--------------------------------------------------------
	  //
	  _lp->SolutionUpdated( &*_d, _eta );
	  //
	  if (_iter%2) {
	    //
	    //--------------------------------------------------------
	    // Compute the new rho, beta if we need to.
	    //--------------------------------------------------------
	    //
	    MVT::MvTransMv( one, *_w, *_rtilde, _rho );       // rho = ( w, r_tilde )
	    _beta = _rho(0,0)/_rho_old(0,0);                  // beta = rho / rho_old
	    _rho_old(0,0) = _rho(0,0);                        // rho_old = rho
	    //
	    //--------------------------------------------------------
	    // Update u, v, and Au if we need to.
	    // Note: We are updating v in two stages to be memory efficient
	    //--------------------------------------------------------
	    //
	    MVT::MvAddMv( one, *_w, _beta, *_u, *_u );       // u = w + beta*u

	    // First stage of v update.
	    MVT::MvAddMv( one, *_Au, _beta, *_v, *_v );      // v = Au + beta*v 

	    // Update Au.
	    {
	      Teuchos::TimeMonitor OpTimer(*_timerOp);
	      _lp->Apply( *_u, *_Au );                       // Au = A*u
	    }

	    // Second stage of v update.
	    MVT::MvAddMv( one, *_Au, _beta, *_v, *_v );      // v = Au + beta*v
	  }

	} // end of the main TFQMR loop -- for(_iter = 0;...)
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
      if (_om->isVerbosityAndPrint( FinalSummary )) {
	_stest->Print(*_os);
      }
      //
    } // end while ( _cur_block_sol && _cur_block_rhs )
    // **********************************************************************************
    //
    // Print timing details 

    // Stop timer.
    _timerTotal->stop();

    // Reset format that will be used to print the summary
    Teuchos::TimeMonitor::format().setPageWidth(54);

    if (_om->isVerbosity( Belos::TimingDetails )) {
      if (_om->doPrint())
        *_os <<"********************TIMING DETAILS********************"<<endl;
      Teuchos::TimeMonitor::summarize( *_os );
      if (_om->doPrint())
        *_os <<"******************************************************"<<endl;
    }
  } // end TFQMRSolve()


  // Overridden from Teuchos::Describable
  template<class ScalarType, class MV, class OP>
  std::string 
  TFQMR<ScalarType,MV,OP>::description() const
  {
    std::ostringstream oss;
    oss << "Belos::TFQMR<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
    oss << "{";
    oss << "}";
    return oss.str();
  }

  //
} // namespace Belos
//
#endif // BELOS_TFQMR_HPP
//
// End of file BelosTFQMR.hpp


