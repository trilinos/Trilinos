
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
  void SetRightPrec(Operator<TYPE> * RP) { RP_ = RP; Right_Prec_ = true; };

  //! Set the blocksize of the linear problem.  This information is used to set up the linear problem for block solvers.
  void SetBlockSize(int blocksize) { blocksize_ = blocksize; };

  //! Inform the linear problem that the solver is finished with the current linear system.
  /*! \note This method is to be <b> only </b> used by the solver to inform the linear problem manager that it's
	finished with this block of linear systems.  The next time the Curr(RHS/LHS)Vec() is called, the next
	linear system will be returned.  Computing the next linear system isn't done in this method in case the 
	blocksize is changed.
  */
  void SetCurrLSVec();

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
  void SolutionUpdated( MultiVec<TYPE>* SolnUpdate = 0 );

  //@}
  
  //@{ \name Accessor methods

  //! Get a pointer to the operator A.
  Operator<TYPE> * GetOperator() const { return(A_); };

  //! Get a pointer to the left-hand side X.
  MultiVec<TYPE> * GetLHS() const { return(X_); };

  //! Get a pointer to the right-hand side B.
  MultiVec<TYPE> * GetRHS() const { return(B_); };

  //! Get a pointer to the initial residual vector.
  /*! \note This may be the preconditioned residual, if the linear problem is left-preconditioned.
   */
  MultiVec<TYPE> * GetInitResVec();

  //! Get a pointer to the current residual vector.
  /*! \note <ul>
	    <li> This method computes the true residual of the current linear system
		 obtained using GetCurr[RHS/LHS]Vec().
	    <li> If the solution hasn't been updated in the LinearProblemManager and
                 a current solution has been computed by the solver (like GMRES), it can
                 be passed into this method to update the residual.
 	    </ul>
   */
  MultiVec<TYPE> * GetCurrResVec( MultiVec<TYPE>* CurrSoln = 0 );

  //! Get a pointer to the current left-hand side (solution) of the linear system.
  /*! This method is called by the solver or any method that is interested in the current linear system
	being solved for.  
	<ol>
	<li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
	<li> If there is no linear system to solve, this method will return a NULL pointer
	</ol>
  */
  MultiVec<TYPE> * GetCurrLHSVec();

  //! Get a pointer to the current right-hand side of the linear system.
  /*! This method is called by the solver of any method that is interested in the current linear system
	being solved for.  
	<ol>
	<li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
	<li> If there is no linear system to solve, this method will return a NULL pointer
	</ol>
  */	
  MultiVec<TYPE> * GetCurrRHSVec();
 
  //! Get a pointer to the left preconditioning operator.
  Operator<TYPE> * GetLeftPrec() const { return(LP_); };

  //! Get a pointer to the right preconditioning operator.
  Operator<TYPE> * GetRightPrec() const { return(RP_); };

  //! Get the current blocksize of the linear problem manager.
  int GetBlockSize() const { return( blocksize_ ); };
  
  //! Get the current number of linear system being solved for.
  /*! Since the block size is independent of the number of right-hand sides, 
    it is important to know how many linear systems
    are being solved for when the status is checked.  This is informative for residual
    checks because the entire block of residuals may not be of interest.  Thus, this 
    number can be anywhere between 1 and the blocksize of the linear system.
  */
  int GetNumToSolve() const { return( num_to_solve_ ); };

  //! Get the index of the first vector in the current right-hand side block being solved for.
  /*! Since the block size is independent of the number of right-hand sides for
    some solvers (GMRES, CG, etc.), it is important to know which right-hand sides
    are being solved for.  That may mean you need to update the information
    about the norms of your initial residual vector for weighting purposes.  This
    information can keep you from querying the solver for information that rarely
    changes.
  */
  int GetRHSIndex() const { return( rhs_index_ ); };
  
  //! Get the current status of the solution.
  /*! This only means that the current linear system being solved for ( obtained by GetCurr<LHS/RHS>Vec() )
      has been updated by the solver.  This will be true every iteration for solvers like CG, but not
      true until restarts for GMRES.
  */
  bool IsSolutionUpdated() const { return(solutionUpdated_); };

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

  //! Private method for populating the next block linear system.
  void SetUpBlocks();

  //! Operator of linear system. 
  Operator<TYPE> * A_;

  //! Solution vector of linear system.
  MultiVec<TYPE> * X_;

  //! Current solution vector of the linear system.
  MultiVec<TYPE> * CurX_;

  //! Right-hand side of linear system.
  MultiVec<TYPE> * B_;

  //! Current right-hand side of the linear system.
  MultiVec<TYPE> * CurB_;

  //! Current residual of the linear system.
  MultiVec<TYPE> * R_;

  //! Initial residual of the linear system.
  MultiVec<TYPE> * R0_;

  //! Left preconditioning operator of linear system
  Operator<TYPE> * LP_;  

  //! Right preconditioning operator of linear system
  Operator<TYPE> * RP_;

  //! Block size of linear system.
  int blocksize_;

  //! Number of linear systems that are currently being solver for ( <= blocksize_ )
  int num_to_solve_;

  //! Index of current block of right-hand sides being solver for ( RHS[:, rhs_index_:(rhs_index_+_blocksize)] ).
  int rhs_index_;

  //! Booleans to keep track of linear problem attributes/status.
  bool Left_Prec_;
  bool Right_Prec_;
  bool Left_Scale_;
  bool Right_Scale_;
  bool operatorSymmetric_;
  bool solutionUpdated_;    
  bool solutionFinal_;
  bool initresidsComputed_;
};

  //--------------------------------------------
  //  Constructor Implementations
  //--------------------------------------------

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(void) : 
  A_(0), 
  X_(0), 
  CurX_(0),
  B_(0),
  CurB_(0),
  R_(0),
  R0_(0),
  LP_(0),
  RP_(0),
  blocksize_(1),
  num_to_solve_(0),
  rhs_index_(0),  
  Left_Prec_(false),
  Right_Prec_(false),
  Left_Scale_(false),
  Right_Scale_(false),
  operatorSymmetric_(false),
  solutionUpdated_(false),
  solutionFinal_(true),
  initresidsComputed_(false)
{
}

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(Operator<TYPE> * A, 
						 MultiVec<TYPE> * X, 
						 MultiVec<TYPE> * B):
  A_(A),
  X_(X),
  CurX_(0),
  B_(B),
  CurB_(0),
  R_(0),
  R0_(0),
  LP_(0),
  RP_(0),
  blocksize_(1),
  num_to_solve_(1),
  rhs_index_(0),
  Left_Prec_(false),
  Right_Prec_(false),
  Left_Scale_(false),
  Right_Scale_(false),
  operatorSymmetric_(false),
  solutionUpdated_(false),
  solutionFinal_(true),
  initresidsComputed_(false)
{
  R0_ = X_->Clone( X_->GetNumberVecs() );
}

template<class TYPE>
LinearProblemManager<TYPE>::LinearProblemManager(const LinearProblemManager<TYPE>& Problem) :
  A_(Problem.A_),
  X_(Problem.X_),
  CurX_(Problem.CurX_),
  B_(Problem.B_),
  CurB_(Problem.CurB_),
  R_(Problem.R_),
  R0_(Problem.R0_),
  LP_(Problem.LP_),
  RP_(Problem.RP_),
  blocksize_(Problem.blocksize_),
  num_to_solve_(Problem.num_to_solve_),
  rhs_index_(Problem.rhs_index_),
  Left_Prec_(Problem.Left_Prec_),
  Right_Prec_(Problem.Right_Prec_),
  Left_Scale_(Problem.Left_Scale_),
  Right_Scale_(Problem.Right_Scale_),
  operatorSymmetric_(Problem.operatorSymmetric_),
  solutionUpdated_(Problem.solutionUpdated_),
  solutionFinal_(Problem.solutionFinal_),
  initresidsComputed_(Problem.initresidsComputed_)
{
}

template<class TYPE>
LinearProblemManager<TYPE>::~LinearProblemManager(void)
{
  if (R_) delete R_; R_ = 0;
  if (R0_) delete R0_; R0_ = 0;
  if (CurX_) delete CurX_; CurX_ = 0;
  if (CurB_) delete CurB_; CurB_ = 0;
}

template<class TYPE>
void LinearProblemManager<TYPE>::SetUpBlocks()
{
  // Compute the new block linear system.
  // ( first clean up old linear system )
  if (CurB_) delete CurB_; CurB_ = 0;
  if (CurX_) delete CurX_; CurX_ = 0;
  if (R_) delete R_; R_ = 0;
  //
  // Determine how many linear systems are left to solve for and populate LHS and RHS vector.
  // If the number of linear systems left are less than the current blocksize, then
  // we create a multivector and copy the left over LHS and RHS vectors into them.
  // The rest of the multivector is populated with random vectors (RHS) or zero vectors (LHS).
  //
  num_to_solve_ = X_->GetNumberVecs() - rhs_index_;
  //
  // Return the NULL pointer if we don't have any more systems to solve for.
  if ( num_to_solve_ <= 0 ) { return; }  
  //
  int i;
  int* index = new int[blocksize_];
  for ( i=0; i<blocksize_; i++ ) { index[i] = rhs_index_ + i; }
  //
  if ( num_to_solve_ < blocksize_ ) 
  {
    int* index2 = new int[ num_to_solve_ ];
    for (i=0; i<num_to_solve_; i++) {
      index2[i] = i;
    }
    //
    // First create multivectors of blocksize and fill the RHS with random vectors LHS with zero vectors.
    CurX_ = X_->Clone(blocksize_); assert(CurX_!=NULL); 
    CurX_->MvInit();
    CurB_ = B_->Clone(blocksize_); assert(CurB_!=NULL);
    CurB_->MvRandom();
    R_ = X_->Clone(blocksize_); assert(R_!=NULL);
    //
    MultiVec<TYPE> *tptr = B_->CloneView(index,num_to_solve_); assert(tptr!=NULL);
    CurB_->SetBlock( *tptr, index2, num_to_solve_);
    //
    MultiVec<TYPE> *tptr2 = X_->CloneView(index,num_to_solve_); assert(tptr2!=NULL);
    CurX_->SetBlock( *tptr2, index2, num_to_solve_);
    //
    // Clean up.
    //
    delete tptr; tptr = 0;
    delete tptr2; tptr2 = 0;
    delete [] index2; index2=0;
    //
  } else { 
    //
    // If the number of linear systems left are more than or equal to the current blocksize, then
    // we create a view into the LHS and RHS.
    //
    num_to_solve_ = blocksize_;
    CurX_ = X_->CloneView( index, num_to_solve_);
    CurB_ = B_->CloneView( index, num_to_solve_);
    R_ = X_->Clone( num_to_solve_ );
    //
  }
  //
  // Compute the current residual.
  // 
  if (R_) {
    A_->Apply( *CurX_, *R_ );
    R_->MvAddMv( 1.0, *CurB_, -1.0, *R_ );
    solutionUpdated_ = false;
  }
  delete [] index; index=0;
}

template<class TYPE>
void LinearProblemManager<TYPE>::SetLHS(MultiVec<TYPE> * X)
{
  X_ = X; 
  R0_ = X_->Clone( X_->GetNumberVecs() ); 
}

template<class TYPE>
void LinearProblemManager<TYPE>::SetCurrLSVec() 
{ 
  int i;
  //
  // We only need to copy the solutions back if the linear systems of
  // interest are less than the block size.
  //
  if (num_to_solve_ < blocksize_) {
    //
    int * index = new int[num_to_solve_]; assert(index!=NULL);
    MultiVec<TYPE> *tptr=0;
    //
    // Get a view of the current solutions and correction vector.
    //
    for (i=0; i<num_to_solve_; i++) { 
      index[i] = i;	
    }
    tptr = CurX_->CloneView( index, num_to_solve_ ); assert(tptr!=NULL);
    //
    // Copy the correction vector to the solution vector.
    //
    for (i=0; i<num_to_solve_; i++) { 
      index[i] = rhs_index_ + i; 
    }
    X_->SetBlock( *tptr, index, num_to_solve_);
    //
    // Clean up.
    //
    delete tptr;
    delete [] index; index=0;
  }
  //
  // Get the linear problem ready to determine the next linear system.
  //
  solutionFinal_ = true; 
  rhs_index_ += num_to_solve_; 
}

template<class TYPE>
void LinearProblemManager<TYPE>::SolutionUpdated( MultiVec<TYPE>* SolnUpdate )
{ 
  if (SolnUpdate) {
    if (Right_Prec_) {
      //
      // Apply the right preconditioner before computing the current solution.
      MultiVec<TYPE>* TrueUpdate = SolnUpdate->Clone( SolnUpdate->GetNumberVecs() );
      RP_->Apply( *SolnUpdate, *TrueUpdate ); 
      CurX_->MvAddMv( 1.0, *CurX_, 1.0, *TrueUpdate ); 
      delete TrueUpdate; TrueUpdate = 0;
    } else {
      CurX_->MvAddMv( 1.0, *CurX_, 1.0, *SolnUpdate ); 
    }
  }
  solutionUpdated_ = true; 
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
      A_->Apply( *X_, *R0_ );
      R0_->MvAddMv( 1.0, *B_, -1.0, *R0_ );
      initresidsComputed_ = true;
    }
  return (R0_);
}

template<class TYPE>
MultiVec<TYPE>* LinearProblemManager<TYPE>::GetCurrResVec( MultiVec<TYPE>* CurrSoln ) 
{
  // Compute the residual of the current linear system.
  // This should be used if the solution has been updated.
  // Alternatively, if the current solution has been computed by GMRES
  // this can be passed in and the current residual will be updated using
  // it.
  //
  if (solutionUpdated_) 
    {
      A_->Apply( *GetCurrLHSVec(), *R_ );
      R_->MvAddMv( 1.0, *GetCurrRHSVec(), -1.0, *R_ ); 
      solutionUpdated_ = false;
    }
  else if (CurrSoln) 
    {
      A_->Apply( *CurrSoln, *R_ );
      R_->MvAddMv( 1.0, *GetCurrRHSVec(), -1.0, *R_ ); 
    }
  return (R_);
}

template<class TYPE>
MultiVec<TYPE> * LinearProblemManager<TYPE>::GetCurrLHSVec()
{
  if (solutionFinal_) {
    solutionFinal_ = false;	// make sure we don't populate the current linear system again.
    SetUpBlocks();
  }
  return CurX_; 
}

template<class TYPE>
MultiVec<TYPE> * LinearProblemManager<TYPE>::GetCurrRHSVec()
{
  if (solutionFinal_) {
    solutionFinal_ = false;	// make sure we don't populate the current linear system again.
    SetUpBlocks();
  }
  return CurB_;
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


