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
#include "BelosMultiVec.hpp"
#include <vector>

/*! 
  \class Belos::StatusTestResNorm
  \brief An implementation of StatusTest using a family of residual norms.

  StatusTestResNorm is an implementation of StatusTest that allows a user to construct
  one of a family of residual tests for use as a status/convergence test for Belos.  
  The form of the test is
   \f[
   \frac{\|r\|}{\sigma} \le \tau
   \f]
   where 
   <ul>
   <li> \f$r\f$ is the residual vector, implicitly or explicitly computed (determined by enum ResType),
   <li> \f$\|r\|\f$ is the residual norm determined by the enum NormType  (1-norm, 2-norm or inf-norm), 
   <li> \f$\sigma\f$ is the scale factor that can be passed in as a precomputed double precision number, 
   or can be selected from by the enum ScaleType (norm of RHS, norm of initial residual).
   <li> \f$\tau\f$ is the tolerance that is passed in as a double precision number to the constructor.  
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
  bool ResidualVectorRequired() const { return(_resvecrequired); };

  //@}

  //@{ \name Print methods

  //! Output formatted description of stopping test to output stream.
  ostream& Print(ostream& os, int indent = 0) const;
  //@}

  //@{ \name Methods to access data members.

  //! Returns the value of the tolerance, \f$ \tau \f$, set in the constructor.
  TYPE GetTolerance() const {return(tolerance_);};
  
  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  TYPE* GetTestValue() const {return(testvalue_);};

  //! Returns the residual norm value, \f$ \|r\| \f$, computed in most recent call to CheckStatus.
  TYPE* GetResNormValue() const {return(resvalue_);};

  //! Returns the scaled norm value, \f$ \sigma \f$.
  TYPE* GetScaledNormValue() const {return(scalevalue_);};

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
  NormType scalenormtype_

  //! Scale value.
  TYPE* scalevalue_;
  
  //! Residual norm value.
  TYPE* resvalue_;

  //! Test value = resvalue_ / scalevalue_
  TYPE* testvalue_;
  
  //! Status
  StatusType status_;
  
  //! The current right-hand side block being solved for.
  int rhs_block_;

  //! The current size of the block being solved for.
  int cur_blksz_;

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
      resvalue_(0),
      init_resvalue_(0),
      testvalue_(0),
      status_(Unchecked),
      rhs_block_(0),
      cur_blksz_(0),
      resvecrequired_(false),
      firstcallCheckStatus_(true),
      firstcallDefineResForm_(true),
      firstcallDefineScaleForm_(true)
  {
    // This constructor will compute the residual ||r||/||r0|| <= tolerance using the 2-norm of
    // the implicit residual vector.
  }

  template <class TYPE>
  int StatusTestResNorm<TYPE>::DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm)
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
    // This section computes the norm of the residual vector
    if (restype_==Implicit) {
      
    }
    else if (restype_==Explicit) {
      curresvecexplicit_ = true;
      if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
      // Compute explicit residual
      operator_.Apply(lhs_, *localresvector_);
      localresvector_->Update(1.0, rhs_, -1.0); // localresvector_ = rhs_ - operator_* lhs_
    }
    resvalue_ = ComputeNorm(*localresvector_, resnormtype_);
  }
  else {
    curresvecexplicit_ = false;
    if (resweights_!=0) { // Check if we should scale the vector
      if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
      // localresvector_ = resweights_ * localresvector_
      localresvector_->Multiply(1.0, *resweights_, *crv, 0.0);
      resvalue_ = ComputeNorm(*localresvector_, resnormtype_);
    }
    else
      resvalue_ = ComputeNorm(*crv, resnormtype_);
  }
  // Compute scaling term (done once)
  if (firstcallCheckStatus_) {
    if (scaletype_==NormOfRHS) {
      if (scaleweights_!=0) { // Check if we should scale the vector
        if (localresvector_==0) localresvector_ = new Epetra_Vector(rhs_.Map());
        // localresvector = scaleweights_ * rhs_
        localresvector_->Multiply(1.0, *scaleweights_, rhs_, 0.0);
        scalevalue_ = ComputeNorm(*localresvector_, resnormtype_);
      }
      else {
        scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
      }
    }
    else if (scaletype_==NormOfInitRes) {
      if (restype_==Implicit && scalenormtype_==TwoNorm && CurrentResNormEst!=-1.0)
        scalevalue_ = CurrentResNormEst;
      else {
        if (scaleweights_!=0) { // Check if we should scale the vector
          if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
          // weightedrhs = scaleweights_ * initial residual
          localresvector_->Multiply(1.0, *scaleweights_, *crv, 0.0);
          scalevalue_ = ComputeNorm(*localresvector_, resnormtype_);
        }
        else {
          scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
        }
      }
    }
    if (scalevalue_==0.0) {
      status_ = Failed;
      return(status_);
    }
  }
  
  testvalue_ = resvalue_/scalevalue_;
  if (testvalue_>tolerance_)
    status_ = Unconverged;
  else if (testvalue_<=tolerance_)
    status_ = Converged;
  else
  status_ = NaN;
  
  firstcallCheckStatus_ = false;
  return status_;
}

template <class TYPE>
ostream& StatusTestResNorm<TYPE>::Print(ostream& os, int indent) const
{
  for (int j = 0; j < indent; j ++)
    os << ' ';
  PrintStatus(os, status_);
  os << "(";
  if (resweights_!=0) os << "Weighted ";
  os << ((resnormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
  os << ((curresvecexplicit_) ? " Exp" : " Imp");
  os << " Res Vec) ";
  if (scaletype_!=None)
    os << "/";
  if (scaletype_==UserProvided)
    os << " (User Scale)";
  else {
    os << "(";
    if (scaleweights_!=0) os << "Weighted ";
    os << ((scalenormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
    if (scaletype_==NormOfInitRes)
      os << " Res0";
    else
      os << " RHS ";
    os << ")";
  }
  if (status_==Unchecked)
    os << " Unchecked << ";
  else {
    os << " = " << testvalue_;
    os << ((testvalue_<tolerance_) ? " < " : (testvalue_==tolerance_) ? " = " : (testvalue_>tolerance_) ? "$  }
  os << tolerance_;
  os << endl;

  return os;
}

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_RESNORM_H */
