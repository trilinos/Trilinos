
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
  LinearProblemManager(const RefCountPtr<const Operator<TYPE> > &A, 
		       const RefCountPtr<MultiVec<TYPE> > &X, 
		       const RefCountPtr<const MultiVec<TYPE> > &B
		       );
  
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
  void SetOperator(const RefCountPtr<const Operator<TYPE> > &A) { A_ = A; };

  //! Set left-hand-side X of linear problem AX = B.
  /*! Sets a pointer to a MultiVec.  No copy of the object is made.
  */
  void SetLHS(const RefCountPtr<MultiVec<TYPE> > &X);

  //! Set right-hand-side B of linear problem AX = B.
  /*! Sets a pointer to a MultiVec.  No copy of the object is made.
  */
  void SetRHS(const RefCountPtr<const MultiVec<TYPE> > &B) { B_ = B; };

  //! Set left preconditioning operator (\c LP) of linear problem AX = B.
  /*! Sets a pointer to an Operator.  No copy of the operator is made.
  */
  void SetLeftPrec(const RefCountPtr<const Operator<TYPE> > &LP) {  LP_ = LP; Left_Prec_ = true; };

  //! Set right preconditioning operator (\c RP) of linear problem AX = B.
  /*! Sets a pointer to an Operator.  No copy of the operator is made.
  */
  void SetRightPrec(const RefCountPtr<const Operator<TYPE> > &RP) { RP_ = RP; Right_Prec_ = true; };

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
  void SolutionUpdated( const MultiVec<TYPE>* SolnUpdate = 0 );

  //@}
  
  //@{ \name Reset method
  
  //! Reset the linear problem manager.
  /*! This is useful for solving the linear system with another right-hand side.  
    The internal flags will be set as if the linear system manager was just initialized.
   */
  void Reset( const RefCountPtr<MultiVec<TYPE> > &newX = null, const RefCountPtr<const MultiVec<TYPE> > &newB = null );
  //@}

  //@{ \name Accessor methods

  //! Get a pointer to the operator A.
  RefCountPtr<const Operator<TYPE> > GetOperator() const { return(A_); };

  //! Get a pointer to the left-hand side X.
  RefCountPtr<MultiVec<TYPE> > GetLHS() const { return(X_); };

  //! Get a pointer to the right-hand side B.
  RefCountPtr<const MultiVec<TYPE> > GetRHS() const { return(B_); };

  //! Get a pointer to the initial residual vector.
  /*! \note This may be the preconditioned residual, if the linear problem is left-preconditioned.
   */
  const MultiVec<TYPE>& GetInitResVec();

  //! Get a pointer to the current residual vector.
  /*!
    
  \param  CurrSoln  [in] If non-null, then this is the LHS that is used to compute
                    the current residual.  If null, then GetCurrLHSVec() is used.

  Note, the current residual is always computed with respect to GetCurrRHSVec().		    

  \note <ul>
   <li> This method computes the true residual of the current linear system
	 with respect to GetCurrRHSVec() and GetCurrLHSVec() if CurrSoln==NULL
         or with respect to *CurrSoln if CurrSoln!=NULL.  
   <li> If the solution hasn't been updated in the LinearProblemManager and
        a current solution has been computed by the solver (like GMRES), it can
        be passed into this method to compute the residual.
   </ul>
   */
  const MultiVec<TYPE>& GetCurrResVec( const MultiVec<TYPE>* CurrSoln = 0 );

  //! Get a pointer to the current left-hand side (solution) of the linear system.
  /*! This method is called by the solver or any method that is interested in the current linear system
	being solved for.  
	<ol>
	<li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
	<li> If there is no linear system to solve, this method will return a NULL pointer
	</ol>
  */
  RefCountPtr<MultiVec<TYPE> > GetCurrLHSVec();

  //! Get a pointer to the current right-hand side of the linear system.
  /*! This method is called by the solver of any method that is interested in the current linear system
	being solved for.  
	<ol>
	<li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
	<li> If there is no linear system to solve, this method will return a NULL pointer
	</ol>
  */	
  RefCountPtr<MultiVec<TYPE> > GetCurrRHSVec();
 
  //! Get a pointer to the left preconditioning operator.
  RefCountPtr<const Operator<TYPE> > GetLeftPrec() const { return(LP_); };

  //! Get a pointer to the right preconditioning operator.
  RefCountPtr<const Operator<TYPE> > GetRightPrec() const { return(RP_); };

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
    
    Precondition:<ul>
    <li><tt>GetOperator().get()!=NULL</tt>
    </ul>
  */
  ReturnType Apply( const MultiVec<TYPE>& x, MultiVec<TYPE>& y );

  //! Apply ONLY the operator to \c x, returning \c y.
  /*! This application is only of the linear problem operator, no preconditioners are applied.
    Flexible variants of Krylov methods will use this application method within their code.
    
    Precondition:<ul>
    <li><tt>GetOperator().get()!=NULL</tt>
    </ul>
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
  ReturnType ComputeResVec( MultiVec<TYPE>* R, const MultiVec<TYPE>* X = 0, const MultiVec<TYPE>* B = 0 );

  //@}

 private:

  //! Private method for populating the next block linear system.
  void SetUpBlocks();

  //! Operator of linear system. 
  RefCountPtr<const Operator<TYPE> > A_;

  //! Solution vector of linear system.
  RefCountPtr<MultiVec<TYPE> > X_;

  //! Current solution vector of the linear system.
  RefCountPtr<MultiVec<TYPE> > CurX_;

  //! Right-hand side of linear system.
  RefCountPtr<const MultiVec<TYPE> > B_;

  //! Current right-hand side of the linear system.
  RefCountPtr<MultiVec<TYPE> > CurB_;

  //! Current residual of the linear system.
  RefCountPtr<MultiVec<TYPE> > R_;

  //! Initial residual of the linear system.
  RefCountPtr<MultiVec<TYPE> > R0_;

  //! Left preconditioning operator of linear system
  RefCountPtr<const Operator<TYPE> > LP_;  

  //! Right preconditioning operator of linear system
  RefCountPtr<const Operator<TYPE> > RP_;

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
LinearProblemManager<TYPE>::LinearProblemManager(const RefCountPtr<const Operator<TYPE> > &A, 
						 const RefCountPtr<MultiVec<TYPE> > &X, 
						 const RefCountPtr<const MultiVec<TYPE> > &B
						 ) :
  A_(A),
  X_(X),
  B_(B),
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
  R0_ = rcp( X_->Clone( X_->GetNumberVecs() ) );
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
{}

template<class TYPE>
void LinearProblemManager<TYPE>::SetUpBlocks()
{
  // Compute the new block linear system.
  // ( first clean up old linear system )
  if (CurB_.get()) CurB_ = null;
  if (CurX_.get()) CurX_ = null;
  if (R_.get()) R_ = null;
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
    std::vector<int> index2(num_to_solve_);
    //int* index2 = new int[ num_to_solve_ ];
    for (i=0; i<num_to_solve_; i++) {
      index2[i] = i;
    }
    //
    // First create multivectors of blocksize and fill the RHS with random vectors LHS with zero vectors.
    CurX_ = rcp( X_->Clone(blocksize_) );
    CurX_->MvInit();
    CurB_ = rcp( const_cast<MultiVec<TYPE>&>(*B_).Clone(blocksize_) );
    CurB_->MvRandom();
    R_ = rcp( X_->Clone(blocksize_) );
    //
    RefCountPtr<MultiVec<TYPE> > tptr = rcp( const_cast<MultiVec<TYPE>&>(*B_).CloneView(index,num_to_solve_) );
    CurB_->SetBlock( *tptr, &index2[0], num_to_solve_);
    //
    RefCountPtr<MultiVec<TYPE> > tptr2 = rcp( X_->CloneView(index,num_to_solve_) );
    CurX_->SetBlock( *tptr2, &index2[0], num_to_solve_);
  } else { 
    //
    // If the number of linear systems left are more than or equal to the current blocksize, then
    // we create a view into the LHS and RHS.
    //
    num_to_solve_ = blocksize_;
    CurX_ = rcp( X_->CloneView( index, num_to_solve_) );
    CurB_ = rcp( const_cast<MultiVec<TYPE>&>(*B_).CloneView( index, num_to_solve_) );
    R_ = rcp( X_->Clone( num_to_solve_ ) );
    //
  }
  //
  // Compute the current residual.
  // 
  if (R_.get()) {
    A_->Apply( *CurX_, *R_ );
    R_->MvAddMv( 1.0, *CurB_, -1.0, *R_ );
    solutionUpdated_ = false;
  }
  delete [] index; index=0;
}

template<class TYPE>
void LinearProblemManager<TYPE>::SetLHS(const RefCountPtr<MultiVec<TYPE> > &X)
{
  X_ = X; 
  R0_ = rcp( X_->Clone( X_->GetNumberVecs() ) ); 
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
    RefCountPtr<MultiVec<TYPE> > tptr;
    //
    // Get a view of the current solutions and correction vector.
    //
    for (i=0; i<num_to_solve_; i++) { 
      index[i] = i;	
    }
    tptr = rcp( CurX_->CloneView( index, num_to_solve_ ) );
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
    delete [] index; index=0;
  }
  //
  // Get the linear problem ready to determine the next linear system.
  //
  solutionFinal_ = true; 
  rhs_index_ += num_to_solve_; 
}

template<class TYPE>
void LinearProblemManager<TYPE>::SolutionUpdated( const MultiVec<TYPE>* SolnUpdate )
{ 
  if (SolnUpdate) {
    if (Right_Prec_) {
      //
      // Apply the right preconditioner before computing the current solution.
      RefCountPtr<MultiVec<TYPE> > TrueUpdate = rcp( const_cast<MultiVec<TYPE>*>(SolnUpdate)->Clone( SolnUpdate->GetNumberVecs() ) );
      RP_->Apply( *SolnUpdate, *TrueUpdate ); 
      CurX_->MvAddMv( 1.0, *CurX_, 1.0, *TrueUpdate ); 
    } else {
      CurX_->MvAddMv( 1.0, *CurX_, 1.0, const_cast<MultiVec<TYPE>&>(*SolnUpdate) ); 
    }
  }
  solutionUpdated_ = true; 
}

template<class TYPE>
void LinearProblemManager<TYPE>::Reset( const RefCountPtr<MultiVec<TYPE> > &newX, const RefCountPtr<const MultiVec<TYPE> > &newB )
{
  solutionUpdated_ = false;
  solutionFinal_ = true;
  initresidsComputed_ = false;
  rhs_index_ = 0;

  X_ = newX;
  B_ = newB;
  GetInitResVec();
}

template<class TYPE>
const MultiVec<TYPE>& LinearProblemManager<TYPE>::GetInitResVec() 
{
  // Compute the initial residual if it hasn't been computed
  // and all the components of the linear system are there.
  // The left preconditioner will be applied if it exists, resulting
  // in a preconditioned residual.
  if (!initresidsComputed_ && A_.get() && X_.get() && B_.get()) 
    {
      if (R0_.get()) R0_ = null;
      R0_ = rcp( X_->Clone( X_->GetNumberVecs() ) );
      A_->Apply( *X_, *R0_ );
      R0_->MvAddMv( 1.0, const_cast<MultiVec<TYPE>&>(*B_), -1.0, *R0_ );
      initresidsComputed_ = true;
    }
  return (*R0_);
}

template<class TYPE>
const MultiVec<TYPE>& LinearProblemManager<TYPE>::GetCurrResVec( const MultiVec<TYPE>* CurrSoln ) 
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
  return (*R_);
}

template<class TYPE>
RefCountPtr<MultiVec<TYPE> > LinearProblemManager<TYPE>::GetCurrLHSVec()
{
  if (solutionFinal_) {
    solutionFinal_ = false;	// make sure we don't populate the current linear system again.
    SetUpBlocks();
  }
  return CurX_; 
}

template<class TYPE>
RefCountPtr<MultiVec<TYPE> > LinearProblemManager<TYPE>::GetCurrRHSVec()
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
  RefCountPtr<MultiVec<TYPE> > ytemp = rcp(y.Clone(y.GetNumberVecs()));
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
  return Ok;
}

template<class TYPE>
ReturnType LinearProblemManager<TYPE>::ApplyOp( const MultiVec<TYPE>& x, MultiVec<TYPE>& y )
{
  if (A_.get())
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
ReturnType LinearProblemManager<TYPE>::ComputeResVec( MultiVec<TYPE>* R, const MultiVec<TYPE>* X, const MultiVec<TYPE>* B )
{
  if (X && B) // The entries are specified, so compute the residual of Op(A)X = B
    {
      if (Left_Prec_)
	{
	  RefCountPtr<MultiVec<TYPE> > R_temp = rcp( const_cast<MultiVec<TYPE>*>(X)->Clone( X->GetNumberVecs() ) );
	  A_->Apply( *X, *R_temp );
	  R_temp->MvAddMv( -1.0, *R_temp, 1.0, const_cast<MultiVec<TYPE>&>(*B) );
	  LP_->Apply( *R_temp, *R );
	}
      else 
	{
	  A_->Apply( *X, *R );
	  R->MvAddMv( -1.0, *R, 1.0, const_cast<MultiVec<TYPE>&>(*B) );
	}
    }
  else { 
    // One of the entries is not specified, so just use the linear system information we have.
    // Later we may want to check to see which multivec is not specified, and use what is specified.
    if (Left_Prec_)
	{
	  RefCountPtr<MultiVec<TYPE> > R_temp = rcp( X_->Clone( X_->GetNumberVecs() ) );
	  A_->Apply( *X_, *R_temp );
	  R_temp->MvAddMv( -1.0, *R_temp, 1.0, const_cast<MultiVec<TYPE>&>(*B_) );
	  LP_->Apply( *R_temp, *R );
	}
      else 
	{
	  A_->Apply( *X_, *R );
	  R->MvAddMv( -1.0, *R, 1.0, const_cast<MultiVec<TYPE>&>(*B_) );
	}
  }    
  return Ok;
}

} // end Belos namespace

#endif /* BELOS_LINEAR_PROBLEM_MANAGER_HPP */


