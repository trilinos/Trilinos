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

#ifndef BELOS_STATUS_TEST_RESNORM_H
#define BELOS_STATUS_TEST_RESNORM_H

/*!
  \file BelosStatusTestResNorm.hpp
  \brief Belos::StatusTest for specifying a residual norm stopping criteria.
*/

#include "BelosStatusTest.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosMultiVec.hpp"
#include <vector>

/*! 
  \class Belos::StatusTestResNorm
  \brief An implementation of StatusTest using a family of residual norms.

  StatusTestResNorm is an implementation of StatusTest that allows a user to construct
  one of a family of residual tests for use as a status/convergence test for Belos.  
  The form of the test is
   \f[
   \frac{\|r_i\|}{\sigma_i} \le \tau
   \f]
   where 
   <ul>
   <li> \f$r_i\f$ is the i-th residual vector, implicitly or explicitly computed (determined by enum ResType),
   <li> \f$\|r_i\|\f$ is the i-th residual norm determined by the enum NormType  (1-norm, 2-norm or inf-norm), 
   <li> \f$\sigma_i\f$ is the i-th scale factor that can be passed in as a precomputed number of the templated type, 
   or can be selected from by the enum ScaleType (norm of RHS, norm of initial residual).
   <li> \f$\tau\f$ is the tolerance that is passed in as a number of the templated type to the constructor.  
   The value of \f$\tau\f$ can be reset using the ResetTolerance() method.
   </ul>

*/

namespace Belos {

template <class TYPE>
class StatusTestResNorm: public StatusTest<TYPE> {

 public:

  //@{ \name Enums.
  /*! 
    \brief Select how the residual vector is produced.
  */
  enum ResType {Implicit, /*!< Use the residual vector produced by the iterative solver. */
		Explicit  /*!< Explicitly compute the residual vector r = b - A*x using the 
			    linear problem. */
  };
  /*! 

    \brief Select the scale type.
  */
  enum ScaleType {NormOfRHS,     /*!< Use the norm of the right-hand-side. */
		  NormOfInitRes, /*!< Use the initial residual vector (exlicitly computed). */
		  None,          /*!< Use unscaled residual. */
		  UserProvided   /*!< User provides an explicit value that the norm of 
				   the residual will be divided by. */
  };
  //@}

  //@{ \name Constructors/destructors.
  //! Constructor
  /*! The constructor takes a single argument specifying the tolerance (\f$\tau\f$).  
    If none of the form definition methods are called, we use \f$\|r\|_2/\|r^{(0)}\|_2 \le \tau\f$ 
    as the stopping criterion, where \f$\|r\|_2\f$ uses the least costly form of the 2-norm of 
    residual available from the iterative method and \f$\|r^{(0)}\|_2\f$ is the corresponding norm 
    of the initial residual.  The least costly form of the 2-norm depends on the chosen iterative 
    method.  Most Krylov methods produce the preconditioned residual vector in a form that would be 
    exact in infinite precision arithmetic.  This vector may be different from the true residual 
    either because left scaling or preconditioning was used, or because round-off error has 
    introduced significant error, or both.
  */
  StatusTestResNorm( TYPE Tolerance );

  //! Destructor
  virtual ~StatusTestResNorm() {};
  //@}

  //@{ \name Form and parameter definition methods.

  //! Define form of the residual, its norm and optional weighting vector.
  /*! This method defines the form of \f$\|r\|\f$.  We specify:
    <ul>
    <li> Whether the residual vector should be explicitly computed, or taken from the iterative method.
    <li> The norm to be used on the residual (this may be different than the norm used in 
    DefineScaleForm()).
    </ul>
  */
  int DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm);
  
  //! Define form of the scaling, its norm, its optional weighting vector, or, alternatively, define an explicit value.
  /*! This method defines the form of how the residual is scaled (if at all).  It operates in two modes:
    <ol>
    <li> User-provided scaling value:
    <ul> 
    <li> Set argument TypeOfScaling to UserProvided.
    <li> Set ScaleValue to a non-zero value that the residual norm will be divided by.
    <li> TypeOfNorm argument will be ignored.
    <li> Sample use:  Define ScaleValue = \f$\|A\|_{\infty}\f$ where \f$ A \f$ is the matrix 
    of the linear problem.
    </ul>
    
    <li> Use a supported Scaling Form:
    <ul>
    <li> Define TypeOfScaling to be the norm of the right hand side, the initial residual vector, 
    or to none.
    <li> Define norm to be used on the scaling vector (this may be different than the norm used 
    in DefineResForm()).
    </ul>
    </ol>
  */
  int DefineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, TYPE ScaleValue = 1.0);

  //! Reset the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual, 
    we find that the initial tolerance was too tight or too lax.
  */
  int ResetTolerance(TYPE Tolerance) {tolerance_ = Tolerance; return(0);};
  //@}

  //@{ \name Status methods
  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  
    Depending on how the residual test is constructed this method will return 
    the appropriate status type.

    \return StatusType: Unconverged, Converged or Failed.
  */
  StatusType CheckStatus(IterativeSolver<TYPE>* iSolver);

  //! Return the result of the most recent CheckStatus call.
  StatusType GetStatus() const {return(status_);};
  //@}

  //@{ \name Attribute methods

  //! Indicates if residual vector is required by this convergence test.
  /*! The value returned by this method will depend on several factors.  
    Once an StatusTestResNorm object is constructed and the DefineResForm 
    and DefineScaleForm methods are optionally called, this method can tested.  
    For most Krylov solvers, there is no extra cost to providing the residual vector.  
    However, GMRES and Transpose-free QMR will need to explicitly compute this vector 
    if ResidualVectorRequired() returns true, so this is an extra cost for
    these two iterative methods.
  */
  bool ResidualVectorRequired() const { return(resvecrequired_); };

  //@}

  //@{ \name Print methods

  //! Output formatted description of stopping test to output stream.
  ostream& Print(ostream& os, int indent = 0) const;
  //@}

  //@{ \name Methods to access data members.

  //! Returns the value of the tolerance, \f$ \tau \f$, set in the constructor.
  TYPE GetTolerance() const {return(tolerance_);};
  
  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  TYPE* GetTestValue() const {return(testvector_);};

  //! Returns the residual norm value, \f$ \|r\| \f$, computed in most recent call to CheckStatus.
  TYPE* GetResNormValue() const {return(resvector_);};

  //! Returns the scaled norm value, \f$ \sigma \f$.
  TYPE* GetScaledNormValue() const {return(scalevector_);};

  //@}

 protected:

 private:

  //@{ \name Private data members.
  
  //! Tolerance used to determine convergence
  TYPE tolerance_;
 
  //! Type of residual to use (explicit or implicit)
  ResType restype_;
  
  //! Type of norm to use on residual (OneNorm, TwoNorm, or InfNorm).
  NormType resnormtype_;
  
  //! Type of scaling to use (Norm of RHS, Norm of Initial Residual, None or User provided)
  ScaleType scaletype_;
  
  //! Type of norm to use on the scaling (OneNorm, TwoNorm, or InfNorm)
  NormType scalenormtype_;

  //! Scaling value.
  TYPE scalevalue_;

  //! Scaling vector.
  TYPE* scalevector_;
  
  //! Residual norm vector.
  TYPE* resvector_;

  //! Test vector = resvector_ / scalevector_
  TYPE* testvector_;
  
  //! Status
  StatusType status_;
  
  //! The index of the right-hand side vector at the beginning of the block being solved for.
  int cur_rhs_num_;

  //! The current size of the block being solved for.
  int cur_blksz_;

  //! The total number of right-hand sides being solved for.
  int numrhs_;

  //! Is residual vector required?
  bool resvecrequired_;

  //! Is this the first time CheckStatus is called?
  bool firstcallCheckStatus_;

  //! Is this the first time DefineResForm is called?
  bool firstcallDefineResForm_;

  //! Is this the first time DefineScaleForm is called?
  bool firstcallDefineScaleForm_;

  //@}

};


  template <class TYPE>
  StatusTestResNorm<TYPE>::StatusTestResNorm( TYPE Tolerance )
    : tolerance_(Tolerance),
      restype_(Implicit),
      resnormtype_(TwoNorm),	
      scaletype_(NormOfInitRes),
      scalenormtype_(TwoNorm),
      scalevalue_(1.0),
      scalevector_(0),
      resvector_(0),
      testvector_(0),
      status_(Unchecked),
      cur_rhs_num_(0),
      cur_blksz_(0),
      numrhs_(0),
      resvecrequired_(false),
      firstcallCheckStatus_(true),
      firstcallDefineResForm_(true),
      firstcallDefineScaleForm_(true)
  {
    // This constructor will compute the residual ||r_i||/||r0_i|| <= tolerance using the 2-norm of
    // the implicit residual vector.
  }

  template <class TYPE>
  int StatusTestResNorm<TYPE>::DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm )
  {    
    assert( firstcallDefineResNorm_ );
    firstcallDefineResForm_ = false;
    
    restype_ = TypeOfResidual;
    resnormtype_ = TypeOfNorm;
    
    // These conditions force the residual vector to be computed
    if (restype_==Explicit)
      resvecrequired_ = true;
    
    return(0);
  }

  template <class TYPE> 
  int StatusTestResNorm<TYPE>::DefineScaleForm(ScaleType TypeOfScaling, NormType TypeOfNorm,
					       TYPE ScaleValue )
  {
    
    assert( firstcallDefineScaleForm_ );
    firstcallDefineScaleForm_ = false;
    
    scaletype_ = TypeOfScaling;
    scalenormtype_ = TypeOfNorm;
    scalevalue_ = ScaleValue;
    
    return(0);
  }

  template <class TYPE>
  StatusType StatusTestResNorm<TYPE>::CheckStatus( IterativeSolver<TYPE>* iSolver )
  {
  int i;
  ReturnType ret;
  LinearProblemManager<TYPE>& lp = iSolver->GetLinearProblem();
  // Compute scaling term (done once for each block that's being solved)
  if (firstcallCheckStatus_) {
    //
    // Get some current solver information.
    //
    firstcallCheckStatus_ = false;
    cur_rhs_num_ = lp.GetRHSIndex();
    cur_blksz_ = lp.GetNumToSolve();
    //
    if (scaletype_== NormOfRHS) {
      MultiVec<TYPE>* rhs = lp.GetRHS();
      numrhs_ = rhs->GetNumberVecs();
      scalevector_ = new TYPE[ numrhs_ ];
      resvector_ = new TYPE[ numrhs_ + cur_blksz_ ];
      testvector_ = new TYPE[ numrhs_ ];
      rhs->MvNorm( scalevector_, scalenormtype_ );
    }
    else if (scaletype_==NormOfInitRes) {
      MultiVec<TYPE>* init_res = lp.GetInitResVec();
      numrhs_ = init_res->GetNumberVecs();
      scalevector_ = new TYPE[ numrhs_ + cur_blksz_ ];
      resvector_ = new TYPE[ numrhs_ + cur_blksz_ ];
      testvector_ = new TYPE[ numrhs_ ];
      init_res->MvNorm( scalevector_, scalenormtype_ );
    }
    // Initialize the testvector.
    for (i=0; i<numrhs_; i++) { testvector_[i] = 1.0; }

    // Return an error if the scaling is zero.
    if (scalevalue_==0.0) {
      status_ = Failed;
      return(status_);
    }
  }
  //
  // This section computes the norm of the residual vector
  //
  if ( cur_rhs_num_ != lp.GetRHSIndex() || cur_blksz_ != lp.GetNumToSolve() ) {
    //
    // We have moved on to the next rhs block
    //
    cur_rhs_num_ = lp.GetRHSIndex();
    cur_blksz_ = lp.GetNumToSolve();
  }
  if (restype_==Implicit) {
    //
    // Get the native residual norms from the solver for this block of right-hand sides.
    // If the residual is returned in multivector form, use the resnormtype to compute the residual norms.
    // Otherwise the native residual is assumed to be stored in the resvector_.
    //
    MultiVec<TYPE>* residMV = iSolver->GetNativeResiduals( resvector_ + cur_rhs_num_ );     
    if ( residMV != NULL ) { 
  	residMV->MvNorm( resvector_ + cur_rhs_num_, resnormtype_ );    
	delete residMV;
    } 
  }
  else if (restype_==Explicit) {
    //
    // Request the true residual for this block of right-hand sides.
    // See if the linear problem manager has been updated before
    // asking for the true residual from the solver.
    // (NOTE:  This residual may be preconditioned)
    //
    MultiVec<TYPE>* cur_res=0;
    if ( lp.IsSolutionUpdated() ) {
      cur_res = lp.GetCurrResVec();
      ret = cur_res->MvNorm( resvector_ + cur_rhs_num_, resnormtype_ );
      if ( ret != Ok ) {
        status_ = Failed;
        return(status_);
      }
    } else {
      MultiVec<TYPE>* cur_soln = iSolver->GetCurrentSoln();
      cur_res = cur_soln->Clone( cur_blksz_ );
      ret = lp.ComputeResVec( cur_res, cur_soln, lp.GetCurrRHSVec() );
      ret = cur_res->MvNorm( resvector_ + cur_rhs_num_, resnormtype_ );
      if ( ret != Ok ) {
        status_ = Failed;
        return(status_);
      }
      delete cur_soln;
    }
    delete cur_res;
  }
  //
  // Compute the new linear system residuals for testing.
  // (if any of them don't meet the tolerance or are NaN, then we exit with that status)
  //
  status_ = Converged; // This will be set to unconverged or NaN.
  if ( scalevector_ ) {
    for (i = cur_rhs_num_; i < (cur_rhs_num_ + cur_blksz_); i++) {
      testvector_[ i ] = resvector_[ i ] / scalevector_[ i ] / scalevalue_;
      if (testvector_[ i ] > tolerance_)
	status_ = Unconverged;
      else if (testvector_[ i ] < tolerance_) { 
	// do nothing.
      } else {
	status_ = NaN;            
	return(status_); // Return immediately if we detect a NaN.
      }
    } 
  }
  else {
    for (i = cur_rhs_num_; i < (cur_rhs_num_ + cur_blksz_); i++) {
      testvector_[ i ] = resvector_[ i ] / scalevalue_;
      if (testvector_[ i ] > tolerance_)
	status_ = Unconverged;
      else if (testvector_[ i ] < tolerance_) { 
	// do nothing.
      } else {
	status_ = NaN;            
	return(status_); // Return immediately if we detect a NaN.
      }
    } 
  }	
  return status_;
}


template <class TYPE>
ostream& StatusTestResNorm<TYPE>::Print(ostream& os, int indent) const
{
  for (int j = 0; j < indent; j ++)
    os << ' ';
  PrintStatus(os, status_);
  os << "(";
  os << ((resnormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
  os << ((restype_==Explicit) ? " Exp" : " Imp");
  os << " Res Vec) ";
  if (scaletype_!=None)
    os << "/ ";
  if (scaletype_==UserProvided)
    os << " (User Scale)";
  else {
    os << "(";
    os << ((scalenormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
    if (scaletype_==NormOfInitRes)
      os << " Res0";
    else
      os << " RHS ";
    os << ")";
  }
  if (status_==Unchecked)
    os << " Unchecked ( tol = " << tolerance_ << " ) "<<endl;
  else {
    os << endl;
    for ( int i=0; i<numrhs_; i++ ) {
      for (int j = 0; j < indent + 13; j ++)
    	os << ' ';
      os << "residual [ " << i << " ] = " << testvector_[ i ];
      os << ((testvector_[i]<tolerance_) ? " < " : (testvector_[i]==tolerance_) ? " == " : (testvector_[i]>tolerance_) ? " > " : " "  ) << tolerance_ << endl;
    }
  }
  os << endl;
  return os;
}

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_RESNORM_H */
