
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

#ifndef BELOS_LINEAR_PROBLEM_MANAGER_HPP
#define BELOS_LINEAR_PROBLEM_MANAGER_HPP

/*! \file BelosLinearProblemManager.hpp
    \brief Class which describes the linear problem to be solved by the iterative solver.
*/
#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"

/*! \class Belos::LinearProblemManager  
  \brief The Belos::LinearProblemManager class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
  The general information is being held as either Belos::Operator or
  Belos::MultiVec objects.
*/

namespace Belos {

template<class TYPE>
class LinearProblemManager {
    
  public:

  //@{ \name Constructors/Destructor.
  //!  Default Constructor.
  /*! Creates an empty Belos::LinearProblemManager instance. The operator A, left-hand-side X
    and right-hand-side B must be set using the SetOperator(), SetLHS() and SetRHS()
    methods respectively.
  */
  LinearProblemManager(void);
  
  //! Unpreconditioned linear system constructor.
  /*! Creates an unpreconditioned LinearProblemManager instance with the 
    Belos::Operator (\c A), initial guess (\c X), and right hand side (\c B). 
    Preconditioners can be set using the SetLeftPrec() and SetRightPrec() methods, and
    scaling can also be set using the SetLeftScale() and SetRightScale() methods.
  */
  LinearProblemManager(Operator<TYPE> * A, 
		       MultiVec<TYPE> * X, 
		       MultiVec<TYPE> * B);
  
  //! Copy Constructor.
  /*! Makes copy of an existing LinearProblemManager instance.
  */
  LinearProblemManager(const LinearProblemManager<TYPE>& Problem);

  //! Destructor.
  /*! Completely deletes a LinearProblemManager object.  
  */
  virtual ~LinearProblemManager(void);
  //@}
  
  //@{ \name Set methods

  //! Set Operator A of linear problem AX = B.
  /*! Sets a pointer to an Operator.  No copy of the operator is made.
  */
  void SetOperator(Operator<TYPE> * A) { A_ = A; };

  //! Set left-hand-side X of linear problem AX = B.
  /*! Sets a pointer to a MultiVec.  No copy of the object is made.
  */
  void SetLHS(MultiVec<TYPE> * X);

  //! Set right-hand-side B of linear problem AX = B.
  /*! Sets a pointer to a MultiVec.  No copy of the object is made.
  */
  void SetRHS(MultiVec<TYPE> * B) { B_ = B; };

  //! Set left preconditioning operator (\c LP) of linear problem AX = B.
  /*! Sets a pointer to an Operator.  No copy of the operator is made.
  */
  void SetLeftPrec(Operator<TYPE> * LP) {  LP_ = LP; Left_Prec_ = true; };

  //! Set right preconditioning operator (\c RP) of linear problem AX = B.
  /*! Sets a pointer to an Operator.  No copy of the operator is made.
  */
  void SetRightPrec(Operator<TYPE> * RP) { RP_ = RP; Right_Prec_ = true; 
};

  //! Inform the linear problem that the operator is symmetric.
  /*! This knowledge may allow the operator to take advantage of the linear problem symmetry.
    However, this should not be set to true if the preconditioner is not symmetric, or symmetrically
    applied.
  */
  void AssertSymmetric(){ operatorSymmetric_ = true; };

  //! Inform the linear problem that the solution has been updated.
  /*! Next time GetCurrResVecs is called, a new residual will be computed.  This keeps the
    linear problem from having to recompute the residual vector everytime it's asked for if
    the solution hasn't been updated.
  */
  void SolutionUpdated(){ solutionUpdated_ = true; };

  //@}
  
  //@{ \name Accessor methods

  //! Get a pointer to the operator A.
  Operator<TYPE> * GetOperator() const { return(A_); };

  //! Get a pointer to the left-hand-side X.
  MultiVec<TYPE> * GetLHS() const { return(X_); };

  //! Get a pointer to the right-hand-side B.
  MultiVec<TYPE> * GetRHS() const { return(B_); };

  //! Get a pointer to the initial residual vector.
  /*! \note This may be the preconditioned residual, if the linear problem is left-preconditioned.
   */
  MultiVec<TYPE> * GetInitResVec();

  //! Get a pointer to the current residual vector.
  /*! \note This may be the preconditioned residual, if the linear problem is left-preconditioned.
   */
  MultiVec<TYPE> * GetCurrResVec();

  //! Get a pointer to the left preconditioning operator.
  Operator<TYPE> * GetLeftPrec() const { return(LP_); };

  //! Get a pointer to the right preconditioning operator.
  Operator<TYPE> * GetRightPrec() const { return(RP_); };

  //! Get operator symmetry bool.
  bool IsOperatorSymmetric() const { return(operatorSymmetric_); };

  //@}

  //@{ \name Apply / Compute methods
  
  //! Apply the composite operator of this linear problem to \c x, returning \c y.
  /*! This application is the composition of the left/right preconditioner and operator.
    Most Krylov methods will use this application method within their code.
  */
  ReturnType Apply( const MultiVec<TYPE>& x, MultiVec<TYPE>& y );

  //! Apply ONLY the operator to \c x, returning \c y.
  /*! This application is only of the linear problem operator, no preconditioners are applied.
    Flexible variants of Krylov methods will use this application method within their code.
  */
  ReturnType ApplyOp( const MultiVec<TYPE>& x, MultiVec<TYPE>& y );
  
  //! Apply ONLY the left preconditioner to \c x, returning \c y.  
  /*! This application is only of the left preconditioner, which may be required for flexible variants
    of Krylov methods.
    \note This will return Undefined if the left preconditioner is not defined for this operator.
   */  
  ReturnType ApplyLeftPrec( const MultiVec<TYPE>& x, MultiVec<TYPE>& y );
  
  //! Apply ONLY the right preconditioner to \c x, returning \c y.
  /*! This application is only of the right preconditioner, which may be required for flexible variants
    of Krylov methods.
    \note This will return Undefined if the right preconditioner is not defined for this operator.
   */
  ReturnType ApplyRightPrec( const MultiVec<TYPE>& x, MultiVec<TYPE>& y );

  //! Compute a residual \c R for this operator given a solution \c X, and right-hand side \c B.
  /*! This method will compute the residual for the current linear system if \c X and \c B are null pointers.
    The result will be returned into R.  Otherwise <tt>R = OP(A)X - B</tt> will be computed and returned.
    \note This residual will be a preconditioned residual if the system has a left preconditioner.
  */
  ReturnType ComputeResVec( MultiVec<TYPE>* R, MultiVec<TYPE>* X = 0, MultiVec<TYPE>* B = 0 );
  //@}

 private:

  //! Operator of linear system. 
  Operator<TYPE> * A_;

  //! Solution vector of linear system.
  MultiVec<TYPE> * X_;

  //! Right-hand side of linear system.
  MultiVec<TYPE> * B_;

  //! Current residual of the linear system.
  MultiVec<TYPE> * R_;

  //! Initial residual of the linear system.
  MultiVec<TYPE> * R0_;

  //! Left preconditioning operator of linear system
  Operator<TYPE> * LP_;  

  //! Right preconditioning operator of linear system
  Operator<TYPE> * RP_;

  bool Left_Prec_;
  bool Right_Prec_;
  bool Left_Scale_;
  bool Right_Scale_;
  bool operatorSymmetric_;
  bool solutionUpdated_;    
  bool initresidsComputed_;
};

  //--------------------------------------------
  //  Constructor Implementations
  //--------------------------------------------

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(void) : 
  A_(0), 
  X_(0), 
  B_(0),
  R_(0),
  R0_(0),
  LP_(0),
  RP_(0),
  Left_Prec_(false),
  Right_Prec_(false),
  Left_Scale_(false),
  Right_Scale_(false),
  operatorSymmetric_(false),
  solutionUpdated_(false),
  initresidsComputed_(false)
{
}

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(Operator<TYPE> * A, 
						 MultiVec<TYPE> * X, 
						 MultiVec<TYPE> * B):
  A_(A),
  X_(X),
  B_(B),
  R_(0),
  R0_(0),
  LP_(0),
  RP_(0),
  Left_Prec_(false),
  Right_Prec_(false),
  Left_Scale_(false),
  Right_Scale_(false),
  operatorSymmetric_(false),
  solutionUpdated_(true),
  initresidsComputed_(false)
{
  R_ = X_->Clone( X_->GetNumberVecs() );
  R0_ = X_->Clone( X_->GetNumberVecs() );
}

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(const LinearProblemManager<TYPE>& Problem) :
  A_(Problem.A_),
  X_(Problem.X_),
  B_(Problem.B_),
  R_(Problem.R_),
  R0_(Problem.R0_),
  LP_(Problem.LP_),
  RP_(Problem.RP_),
  Left_Prec_(Problem.Left_Prec_),
  Right_Prec_(Problem.Right_Prec_),
  Left_Scale_(Problem.Left_Scale_),
  Right_Scale_(Problem.Right_Scale_),
  operatorSymmetric_(Problem.operatorSymmetric_),
  solutionUpdated_(Problem.solutionUpdated_),
  initresidsComputed_(Problem.initresidsComputed_)
{
}

template<class TYPE>
LinearProblemManager<TYPE>::~LinearProblemManager(void)
{
  if (R_) delete R_; R_ = 0;
  if (R0_) delete R0_; R0_ = 0;
}

template<class TYPE>
void LinearProblemManager<TYPE>::SetLHS(MultiVec<TYPE> * X)
{
  X_ = X; 
  solutionUpdated_ = true; 
  R0_ = X_->Clone( X_->GetNumberVecs() ); 
  R_ = X_->Clone( X_->GetNumberVecs() );
}

template<class TYPE>
MultiVec<TYPE>* LinearProblemManager<TYPE>::GetInitResVec() 
{
  // Compute the initial residual if it hasn't been computed
  // and all the components of the linear system are there.
  // The left preconditioner will be applied if it exists, resulting
  // in a preconditioned residual.
  if (!initresidsComputed_ && A_ && X_ && B_) 
    {
      ComputeResVec( R0_ );
      initresidsComputed_ = true;
    }
  return (R0_);
}

template<class TYPE>
MultiVec<TYPE>* LinearProblemManager<TYPE>::GetCurrResVec() 
{
  // Compute the current residual of the linear system.
  // This should be used if the solution has been updated.
  // In the case of GMRES, this solution is only updated at restarts,
  // or when the solver exists.
  //
  // The left preconditioner will be applied if it exists, resulting
  // in a preconditioned residual.
  //
  if (solutionUpdated_ && A_ && X_ && B_) 
    {
      ComputeResVec( R_ );
      solutionUpdated_ = false;
    }
  return (R_);
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::Apply( const MultiVec<TYPE>& x, MultiVec<TYPE>& y )
{
  MultiVec<TYPE>* ytemp = y.Clone(y.GetNumberVecs());
  //
  // No preconditioning.
  // 
  if (!Left_Prec_ && !Right_Prec_){ A_->Apply( x, y );}
  //
  // Preconditioning is being done on both sides
  //
  else if( Left_Prec_ && Right_Prec_ ) 
    {
      RP_->Apply( x, y );   
      A_->Apply( y, *ytemp );
      LP_->Apply( *ytemp, y );
    }
  //
  // Preconditioning is only being done on the left side
  //
  else if( Left_Prec_ ) 
    {
      A_->Apply( x, *ytemp );
      LP_->Apply( *ytemp, y );
    }
  //
  // Preconditioning is only being done on the right side
  //
  else 
    {
      RP_->Apply( x, *ytemp );
      A_->Apply( *ytemp, y );
    }  
  delete ytemp; ytemp = 0;
  return Ok;
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::ApplyOp( const MultiVec<TYPE>& x, MultiVec<TYPE>& y )
{
  if (A_)
    return ( A_->Apply(x, y) );   
  else
    return Undefined;
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::ApplyLeftPrec( const MultiVec<TYPE>& x, MultiVec<TYPE>& y )
{
  if (Left_Prec_)
    return (LP_->Apply(x, y) );
  else 
    return Undefined;
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::ApplyRightPrec( const MultiVec<TYPE>& x, MultiVec<TYPE>& y )
{
  if (Right_Prec_)
    return (RP_->Apply(x, y) );
  else
    return Undefined;
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::ComputeResVec( MultiVec<TYPE>* R, 
						      MultiVec<TYPE>* X, 
						      MultiVec<TYPE>* B )
{
  if (X && B) // The entries are specified, so compute the residual of Op(A)X = B
    {
      if (Left_Prec_)
	{
	  MultiVec<TYPE>* R_temp = X->Clone( X->GetNumberVecs() );
	  A_->Apply( *X, *R_temp );
	  R_temp->MvAddMv( -1.0, *R_temp, 1.0, *B );
	  LP_->Apply( *R_temp, *R );
	  delete R_temp;
	}
      else 
	{
	  A_->Apply( *X, *R );
	  R->MvAddMv( -1.0, *R, 1.0, *B );
	}
    }
  else { 
    // One of the entries is not specified, so just use the linear system information we have.
    // Later we may want to check to see which multivec is not specified, and use what is specified.
    if (Left_Prec_)
	{
	  MultiVec<TYPE>* R_temp = X_->Clone( X_->GetNumberVecs() );
	  A_->Apply( *X_, *R_temp );
	  R_temp->MvAddMv( -1.0, *R_temp, 1.0, *B_ );
	  LP_->Apply( *R_temp, *R );
	  delete R_temp;
	}
      else 
	{
	  A_->Apply( *X_, *R );
	  R->MvAddMv( -1.0, *R, 1.0, *B_ );
	}
  }    
  return Ok;
}

} // end Belos namespace

#endif /* BELOS_LINEAR_PROBLEM_MANAGER_HPP */

