
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

#ifndef BELOS_LINEAR_PROBLEM_HPP
#define BELOS_LINEAR_PROBLEM_HPP

/*! \file BelosLinearProblem.hpp
    \brief Class which describes the linear problem to be solved by the iterative solver.
*/

#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::null;
using Teuchos::rcp_const_cast;
using Teuchos::ParameterList;

/*! \class Belos::LinearProblem  
  \brief The Belos::LinearProblem class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
*/

namespace Belos {

  //! @name LinearProblem Exceptions
  //@{
  
  /** \brief Exception thrown to signal error with the Belos::LinearProblem object.
   */
  class LinearProblemError : public BelosError
  {public: LinearProblemError(const std::string& what_arg) : BelosError(what_arg) {}};
  
  //@}
  

  template <class ScalarType, class MV, class OP>
  class LinearProblem {
   
  public:
    
    //! @name Constructors/Destructor
    //@{ 
    //!  Default Constructor.
    /*! Creates an empty Belos::LinearProblem instance. The operator A, left-hand-side X
      and right-hand-side B must be set using the setOperator(), setLHS() and setRHS()
      methods respectively.
    */
    LinearProblem(void);
    
    //! Unpreconditioned linear system constructor.
    /*! Creates an unpreconditioned LinearProblem instance with the 
      Belos::Operator (\c A), initial guess (\c X), and right hand side (\c B). 
      Preconditioners can be set using the setLeftPrec() and setRightPrec() methods, and
      scaling can also be set using the setLeftScale() and setRightScale() methods.
    */
    LinearProblem(const RCP<const OP> &A, 
		  const RCP<MV> &X, 
		  const RCP<const MV> &B
		  );
    
    //! Copy Constructor.
    /*! Makes copy of an existing LinearProblem instance.
     */
    LinearProblem(const LinearProblem<ScalarType,MV,OP>& Problem);
    
    //! Destructor.
    /*! Completely deletes a LinearProblem object.  
     */
    virtual ~LinearProblem(void);
    //@}
    
    //! @name Set methods
    //@{ 
    
    //! Set Operator A of linear problem AX = B.
    /*! Sets a pointer to an Operator.  No copy of the operator is made.
     */
    void setOperator(const RCP<const OP> &A) { A_ = A; isSet_=false; }
    
    //! Set left-hand-side X of linear problem AX = B.
    /*! Sets a pointer to a MultiVec.  No copy of the object is made.
     */
    void setLHS(const RCP<MV> &X) { X_ = X; isSet_=false; }
    
    //! Set right-hand-side B of linear problem AX = B.
    /*! Sets a pointer to a MultiVec.  No copy of the object is made.
     */
    void setRHS(const RCP<const MV> &B) { B_ = B; isSet_=false; }
    
    //! Set left preconditioning operator (\c LP) of linear problem AX = B.
    /*! Sets a pointer to an Operator.  No copy of the operator is made.
     */
    void setLeftPrec(const RCP<const OP> &LP) {  LP_ = LP; }
    
    //! Set right preconditioning operator (\c RP) of linear problem AX = B.
    /*! Sets a pointer to an Operator.  No copy of the operator is made.
     */
    void setRightPrec(const RCP<const OP> &RP) { RP_ = RP; }
    
    //! Inform the linear problem that the solver is finished with the current linear system.
    /*! \note This method is to be <b> only </b> used by the solver to inform the linear problem that it's
      finished with this block of linear systems.  The next time the Curr(RHS/LHS)Vec() is called, the next
      linear system will be returned.  Computing the next linear system isn't done in this method in case the 
      blocksize is changed.
    */
    void setCurrLS();
    
    //! Inform the linear problem of the linear systems that need to be solved next.
    /*! Any calls to get the current RHS/LHS vectors after this method is called will return the new
      linear systems indicated by \c index.  The length of \c index is assumed to be the blocksize and entries
      of \c index must be between 0 and the number of vectors in the RHS/LHS multivector.  An entry of the
      \c index std::vector can also be -1, which means this column of the linear system is augmented using a random
      std::vector.
    */
    void setLSIndex(std::vector<int>& index); 
    
    //! Inform the linear problem that the operator is Hermitian.
    /*! This knowledge may allow the operator to take advantage of the linear problem symmetry.
      However, this should not be set to true if the preconditioner is not Hermitian, or symmetrically
      applied.
    */
    void setHermitian(){ isHermitian_ = true; }
   
    //! Set the label prefix used by the timers in this object.  The default is "Belos".
    /*! \note The timers are created during the first call to setProblem().  Any calls to this method to change 
        the label after that will not change the label used in the timer.
    */ 
    void setLabel(const std::string& label) { label_ = label; }

    //! Compute the new solution to the linear system given the /c update.
    /*! \note If \c updateLP is true, then the next time GetCurrResVecs is called, a new residual will be computed.  
      This keeps the linear problem from having to recompute the residual vector everytime it's asked for if
      the solution hasn't been updated.  If \c updateLP is false, the new solution is computed without actually 
      updating the linear problem.
    */
    RCP<MV> updateSolution( const RCP<MV>& update = null,
				    bool updateLP = false,
                                    ScalarType scale = Teuchos::ScalarTraits<ScalarType>::one() );    

    //! Compute the new solution to the linear system given the /c update without updating the linear problem.
    RCP<MV> updateSolution( const RCP<MV>& update = null,
                                    ScalarType scale = Teuchos::ScalarTraits<ScalarType>::one() ) const
    { return const_cast<LinearProblem<ScalarType,MV,OP> *>(this)->updateSolution( update, false, scale ); }

    //@}
    
    //! @name Set / Reset method
    //@{ 
    
    //! Setup the linear problem manager.
    /*! This is useful for solving the linear system with another right-hand side or getting
      the linear problem prepared to solve the linear system that was already passed in.  
      The internal flags will be set as if the linear system manager was just initialized 
      and the initial residual will be computed.
    */
    bool setProblem( const RCP<MV> &newX = null, const RCP<const MV> &newB = null );

    //@}
    
    //! @name Accessor methods
    //@{ 
    
    //! Get a pointer to the operator A.
    RCP<const OP> getOperator() const { return(A_); }
    
    //! Get a pointer to the left-hand side X.
    RCP<MV> getLHS() const { return(X_); }
    
    //! Get a pointer to the right-hand side B.
    RCP<const MV> getRHS() const { return(B_); }
    
    //! Get a pointer to the initial residual vector.
    /*! \note This is the preconditioned residual if the linear system is preconditioned on the left.
     */
    RCP<const MV> getInitResVec() const { return(R0_); }
    
    //! Get a pointer to the preconditioned initial residual vector.
    /*! \note This is the unpreconditioned residual.
     */
    RCP<const MV> getInitPrecResVec() const { return(PR0_); }
    
    //! Get a pointer to the current left-hand side (solution) of the linear system.
    /*! This method is called by the solver or any method that is interested in the current linear system
      being solved for.  
      <ol>
      <li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
      <li> If there is no linear system to solve, this method will return a NULL pointer
      </ol>
    */
    RCP<MV> getCurrLHSVec();
    
    //! Get a pointer to the current right-hand side of the linear system.
    /*! This method is called by the solver of any method that is interested in the current linear system
      being solved for.  
      <ol>
      <li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
      <li> If there is no linear system to solve, this method will return a NULL pointer
      </ol>
    */	
    RCP<MV> getCurrRHSVec();
    
    //! Get a pointer to the left preconditioning operator.
    RCP<const OP> getLeftPrec() const { return(LP_); };
    
    //! Get a pointer to the right preconditioning operator.
    RCP<const OP> getRightPrec() const { return(RP_); };
    
    //! Get the 0-based index std::vector indicating the current linear systems being solved for.
    /*! Since the block size is independent of the number of right-hand sides for
      some solvers (GMRES, CG, etc.), it is important to know which linear systems
      are being solved for.  That may mean you need to update the information
      about the norms of your initial residual std::vector for weighting purposes.  This
      information can keep you from querying the solver for information that rarely
      changes.
      \note The length of the index std::vector is the number of right-hand sides being solved for.
            If an entry of this std::vector is -1 then that linear system is an augmented linear
	    system and doesn't need to be considered for convergence.
      \note The std::vector returned from this method is valid if isProblemSet() returns true.
    */
    const std::vector<int> getLSIndex() const { return(rhsIndex_); }

    //! Get the number of linear systems that have been set with this LinearProblem object.
    /* This can be used by status test classes to determine if the solver manager has advanced 
       and is solving another linear system.
    */
    int getLSNumber() const { return(lsNum_); }

    //@}
    
    //! @name State methods
    //@{ 
    
    //! Get the current status of the solution.
    /*! This only means that the current linear system being solved for ( obtained by getCurr<LHS/RHS>Vec() )
      has been updated by the solver.  This will be true every iteration for solvers like CG, but not
      true until restarts for GMRES.
    */
    bool isSolutionUpdated() const { return(solutionUpdated_); }

    //! If the problem has been set, this will return true.
    bool isProblemSet() const { return(isSet_); }
    
    //! Get the current symmetry of the operator.
    bool isHermitian() const { return(isHermitian_); }
    
    //! Get information on whether the linear system is being preconditioned on the left.
    bool isLeftPrec() const { return(LP_!=null); }

    //! Get information on whether the linear system is being preconditioned on the right.
    bool isRightPrec() const { return(RP_!=null); }
 
    //@}
    
    //! @name Apply / Compute methods
    //@{ 
    
    //! Apply the composite operator of this linear problem to \c x, returning \c y.
    /*! This application is the composition of the left/right preconditioner and operator.
      Most Krylov methods will use this application method within their code.
      
      Precondition:<ul>
      <li><tt>getOperator().get()!=NULL</tt>
      </ul>
    */
    void apply( const MV& x, MV& y ) const;
    
    //! Apply ONLY the operator to \c x, returning \c y.
    /*! This application is only of the linear problem operator, no preconditioners are applied.
      Flexible variants of Krylov methods will use this application method within their code.
      
      Precondition:<ul>
      <li><tt>getOperator().get()!=NULL</tt>
      </ul>
    */
    void applyOp( const MV& x, MV& y ) const;
    
    //! Apply ONLY the left preconditioner to \c x, returning \c y.  
    /*! This application is only of the left preconditioner, which may be required for flexible variants
      of Krylov methods.
      \note This will return Undefined if the left preconditioner is not defined for this operator.
    */  
    void applyLeftPrec( const MV& x, MV& y ) const;
    
    //! Apply ONLY the right preconditioner to \c x, returning \c y.
    /*! This application is only of the right preconditioner, which may be required for flexible variants
      of Krylov methods.
      \note This will return Undefined if the right preconditioner is not defined for this operator.
    */
    void applyRightPrec( const MV& x, MV& y ) const;
    
    //! Compute a residual \c R for this operator given a solution \c X, and right-hand side \c B.
    /*! This method will compute the residual for the current linear system if \c X and \c B are null pointers.
      The result will be returned into R.  Otherwise <tt>R = OP(A)X - B</tt> will be computed and returned.
      \note This residual will not be preconditioned if the system has a left preconditioner.
    */
    void computeCurrResVec( MV* R , const MV* X = 0, const MV* B = 0 ) const;

    //! Compute a residual \c R for this operator given a solution \c X, and right-hand side \c B.
    /*! This method will compute the residual for the current linear system if \c X and \c B are null pointers.
      The result will be returned into R.  Otherwise <tt>R = OP(A)X - B</tt> will be computed and returned.
      \note This residual will be preconditioned if the system has a left preconditioner.
    */
    void computeCurrPrecResVec( MV* R, const MV* X = 0, const MV* B = 0 ) const;
    
    //@}
    
  private:
    
    //! Operator of linear system. 
    RCP<const OP> A_;
    
    //! Solution std::vector of linear system.
    RCP<MV> X_;
    
    //! Current solution std::vector of the linear system.
    RCP<MV> curX_;
    
    //! Right-hand side of linear system.
    RCP<const MV> B_;
    
    //! Current right-hand side of the linear system.
    RCP<MV> curB_;
    
    //! Initial residual of the linear system.
    RCP<MV> R0_;
   
    //! Preconditioned initial residual of the linear system.
    RCP<MV> PR0_;
 
    //! Left preconditioning operator of linear system
    RCP<const OP> LP_;  
    
    //! Right preconditioning operator of linear system
    RCP<const OP> RP_;
    
    //! Timers
    mutable Teuchos::RCP<Teuchos::Time> timerOp_, timerPrec_;

    //! Current block size of linear system.
    int blocksize_;

    //! Number of linear systems that are currently being solver for ( <= blocksize_ )
    int num2Solve_;
    
    //! Indices of current linear systems being solver for.
    std::vector<int> rhsIndex_;    

    //! Number of linear systems that have been loaded in this linear problem object.
    int lsNum_;

    //! Booleans to keep track of linear problem attributes/status.
    bool Left_Scale_;
    bool Right_Scale_;
    bool isSet_;
    bool isHermitian_;
    bool solutionUpdated_;    
   
    //! Linear problem label that prefixes the timer labels.
    std::string label_;
 
    typedef MultiVecTraits<ScalarType,MV>  MVT;
    typedef OperatorTraits<ScalarType,MV,OP>  OPT;
  };
  
  //--------------------------------------------
  //  Constructor Implementations
  //--------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  LinearProblem<ScalarType,MV,OP>::LinearProblem(void) : 
    blocksize_(0),
    num2Solve_(0),
    rhsIndex_(0),
    lsNum_(0),
    Left_Scale_(false),
    Right_Scale_(false),
    isSet_(false),
    isHermitian_(false),
    solutionUpdated_(false),
    label_("Belos")
  {
  }
  
  template <class ScalarType, class MV, class OP>
  LinearProblem<ScalarType,MV,OP>::LinearProblem(const RCP<const OP> &A, 
						 const RCP<MV> &X, 
						 const RCP<const MV> &B
						 ) :
    A_(A),
    X_(X),
    B_(B),
    blocksize_(0),
    num2Solve_(0),
    rhsIndex_(0),
    lsNum_(0),
    Left_Scale_(false),
    Right_Scale_(false),
    isSet_(false),
    isHermitian_(false),
    solutionUpdated_(false),
    label_("Belos")
  {
  }
  
  template <class ScalarType, class MV, class OP>
  LinearProblem<ScalarType,MV,OP>::LinearProblem(const LinearProblem<ScalarType,MV,OP>& Problem) :
    A_(Problem.A_),
    X_(Problem.X_),
    curX_(Problem.curX_),
    B_(Problem.B_),
    curB_(Problem.curB_),
    R0_(Problem.R0_),
    PR0_(Problem.PR0_),
    LP_(Problem.LP_),
    RP_(Problem.RP_),
    timerOp_(Problem.timerOp_),
    timerPrec_(Problem.timerPrec_),
    blocksize_(Problem.blocksize_),
    num2Solve_(Problem.num2Solve_),
    rhsIndex_(Problem.rhsIndex_),
    lsNum_(Problem.lsNum_),
    Left_Scale_(Problem.Left_Scale_),
    Right_Scale_(Problem.Right_Scale_),
    isSet_(Problem.isSet_),
    isHermitian_(Problem.isHermitian_),
    solutionUpdated_(Problem.solutionUpdated_),
    label_(Problem.label_)
  {
  }
  
  template <class ScalarType, class MV, class OP>
  LinearProblem<ScalarType,MV,OP>::~LinearProblem(void)
  {}
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::setLSIndex(std::vector<int>& index)
  {
    // Set new linear systems using the indices in index.
    rhsIndex_ = index;
    
    // Compute the new block linear system.
    // ( first clean up old linear system )
    curB_ = null;
    curX_ = null;
   
    // Create indices for the new linear system.
    int validIdx = 0, ivalidIdx = 0;
    blocksize_ = rhsIndex_.size();
    std::vector<int> vldIndex( blocksize_ );
    std::vector<int> newIndex( blocksize_ );
    std::vector<int> iIndex( blocksize_ );
    for (int i=0; i<blocksize_; ++i) {
      if (rhsIndex_[i] > -1) {
        vldIndex[validIdx] = rhsIndex_[i];
        newIndex[validIdx] = i;
        validIdx++;
      }
      else {
        iIndex[ivalidIdx] = i;
        ivalidIdx++;
      }
    }
    vldIndex.resize(validIdx);
    newIndex.resize(validIdx);   
    iIndex.resize(ivalidIdx);
    num2Solve_ = validIdx;

    // Create the new linear system using index
    if (num2Solve_ != blocksize_) {
      newIndex.resize(num2Solve_);
      vldIndex.resize(num2Solve_);
      //
      // First create multivectors of blocksize.
      // Fill the RHS with random vectors LHS with zero vectors.
      curX_ = MVT::Clone( *X_, blocksize_ );
      MVT::MvInit(*curX_);
      curB_ = MVT::Clone( *B_, blocksize_ );
      MVT::MvRandom(*curB_);
      //
      // Now put in the part of B into curB 
      RCP<const MV> tptr = MVT::CloneView( *B_, vldIndex );
      MVT::SetBlock( *tptr, newIndex, *curB_ );
      //
      // Now put in the part of X into curX
      tptr = MVT::CloneView( *X_, vldIndex );
      MVT::SetBlock( *tptr, newIndex, *curX_ );
      //
      solutionUpdated_ = false;
    }
    else {
      curX_ = MVT::CloneView( *X_, rhsIndex_ );
      curB_ = rcp_const_cast<MV>(MVT::CloneView( *B_, rhsIndex_ ));
    }
    //
    // Increment the number of linear systems that have been loaded into this object.
    //
    lsNum_++;
  }


  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::setCurrLS() 
  { 
    //
    // We only need to copy the solutions back if the linear systems of
    // interest are less than the block size.
    //
    if (num2Solve_ < blocksize_) {
      //
      // Get a view of the current solutions and correction std::vector.
      //
      int validIdx = 0;
      std::vector<int> newIndex( num2Solve_ );
      std::vector<int> vldIndex( num2Solve_ );
      for (int i=0; i<blocksize_; ++i) {
        if ( rhsIndex_[i] > -1 ) { 
          vldIndex[validIdx] = rhsIndex_[i];
	  newIndex[validIdx] = i;
          validIdx++;
        }	
      }
      RCP<MV> tptr = MVT::CloneView( *curX_, newIndex );
      MVT::SetBlock( *tptr, vldIndex, *X_ );
    }
    //
    // Clear the current vectors of this linear system so that any future calls
    // to get the vectors for this system return null pointers.
    //
    curX_ = null;
    curB_ = null;
    rhsIndex_.resize(0);
  }
  

  template <class ScalarType, class MV, class OP>
  RCP<MV> LinearProblem<ScalarType,MV,OP>::updateSolution( const RCP<MV>& update, 
							   bool updateLP,
							   ScalarType scale )
  { 
    RCP<MV> newSoln;
    if (update != null) {
      if (updateLP == true) {
	if (RP_!=null) {
	  //
	  // Apply the right preconditioner before computing the current solution.
	  RCP<MV> TrueUpdate = MVT::Clone( *update, MVT::GetNumberVecs( *update ) );
	  {
	    Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	    OPT::Apply( *RP_, *update, *TrueUpdate ); 
	  }
	  MVT::MvAddMv( 1.0, *curX_, scale, *TrueUpdate, *curX_ ); 
	} 
	else {
	  MVT::MvAddMv( 1.0, *curX_, scale, *update, *curX_ ); 
	}
	solutionUpdated_ = true; 
	newSoln = curX_;
      }
      else {
	newSoln = MVT::Clone( *update, MVT::GetNumberVecs( *update ) );
	if (RP_!=null) {
	  //
	  // Apply the right preconditioner before computing the current solution.
	  RCP<MV> trueUpdate = MVT::Clone( *update, MVT::GetNumberVecs( *update ) );
	  {
	    Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	    OPT::Apply( *RP_, *update, *trueUpdate ); 
	  }
	  MVT::MvAddMv( 1.0, *curX_, scale, *trueUpdate, *newSoln ); 
	} 
	else {
	  MVT::MvAddMv( 1.0, *curX_, scale, *update, *newSoln ); 
	}
      }
    }
    else {
      newSoln = curX_;
    }
    return newSoln;
  }
  

  template <class ScalarType, class MV, class OP>
  bool LinearProblem<ScalarType,MV,OP>::setProblem( const RCP<MV> &newX, const RCP<const MV> &newB )
  {
    // Create timers if the haven't been created yet.
    if (timerOp_ == Teuchos::null) {
      std::string opLabel = label_ + ": Operation Op*x";
      timerOp_ = Teuchos::TimeMonitor::getNewTimer( opLabel );
    }
    if (timerPrec_ == Teuchos::null) {
      std::string precLabel = label_ + ": Operation Prec*x";
      timerPrec_ = Teuchos::TimeMonitor::getNewTimer( precLabel );
    }

    // Set the linear system using the arguments newX and newB
    if (newX != null)
      X_ = newX;
    if (newB != null)
      B_ = newB;

    // Invalidate the current linear system indices and multivectors.
    rhsIndex_.resize(0);
    curX_ = null;
    curB_ = null;

    // Check the validity of the linear problem object.
    // If no operator A exists, then throw an std::exception.
    if (A_ == null || X_ == null || B_ == null) {
      isSet_ = false;
      return isSet_;
    }

    // Initialize the state booleans
    solutionUpdated_ = false;
    
    // Compute the initial residuals.
    if (R0_==null || MVT::GetNumberVecs( *R0_ )!=MVT::GetNumberVecs( *X_ )) {
      R0_ = MVT::Clone( *X_, MVT::GetNumberVecs( *X_ ) );
    }
    computeCurrResVec( &*R0_, &*X_, &*B_ );

    if (LP_!=null) {
      if (PR0_==null || MVT::GetNumberVecs( *PR0_ )!=MVT::GetNumberVecs( *X_ )) {
        PR0_ = MVT::Clone( *X_, MVT::GetNumberVecs( *X_ ) );
      }
      {
        Teuchos::TimeMonitor PrecTimer(*timerPrec_);
        OPT::Apply( *LP_, *R0_, *PR0_ );
      }
    } 
    else {
      PR0_ = R0_;
    }    

    // The problem has been set and is ready for use.
    isSet_ = true;
    
    // Return isSet.
    return isSet_;
  }
  
  template <class ScalarType, class MV, class OP>
  RCP<MV> LinearProblem<ScalarType,MV,OP>::getCurrLHSVec()
  {
    if (isSet_) {
      return curX_;
    }
    else {
      return Teuchos::null;
    }
  }
  
  template <class ScalarType, class MV, class OP>
  RCP<MV> LinearProblem<ScalarType,MV,OP>::getCurrRHSVec()
  {
    if (isSet_) {
      return curB_;
    }
    else {
      return Teuchos::null;
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::apply( const MV& x, MV& y ) const
  {
    RCP<MV> ytemp = MVT::Clone( y, MVT::GetNumberVecs( y ) );
    bool leftPrec = LP_!=null;
    bool rightPrec = RP_!=null;
    //
    // No preconditioning.
    // 
    if (!leftPrec && !rightPrec){ 
      Teuchos::TimeMonitor OpTimer(*timerOp_);
      OPT::Apply( *A_, x, y );
    }
    //
    // Preconditioning is being done on both sides
    //
    else if( leftPrec && rightPrec ) 
      {
        {
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	  OPT::Apply( *RP_, x, y );   
        }
        {
          Teuchos::TimeMonitor OpTimer(*timerOp_);
	  OPT::Apply( *A_, y, *ytemp );
        }
        {
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	  OPT::Apply( *LP_, *ytemp, y );
        }
      }
    //
    // Preconditioning is only being done on the left side
    //
    else if( leftPrec ) 
      {
        {
          Teuchos::TimeMonitor OpTimer(*timerOp_);
	  OPT::Apply( *A_, x, *ytemp );
        }
        {
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	  OPT::Apply( *LP_, *ytemp, y );
        }
      }
    //
    // Preconditioning is only being done on the right side
    //
    else 
      {
        {
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	  OPT::Apply( *RP_, x, *ytemp );
        }
        {
          Teuchos::TimeMonitor OpTimer(*timerOp_);
      	  OPT::Apply( *A_, *ytemp, y );
        }
      }  
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyOp( const MV& x, MV& y ) const {
    if (A_.get()) {
      Teuchos::TimeMonitor OpTimer(*timerOp_);
      OPT::Apply( *A_,x, y);   
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x, 
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyLeftPrec( const MV& x, MV& y ) const {
    if (LP_!=null) {
      Teuchos::TimeMonitor PrecTimer(*timerPrec_);
      return ( OPT::Apply( *LP_,x, y) );
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x, 
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyRightPrec( const MV& x, MV& y ) const {
    if (RP_!=null) {
      Teuchos::TimeMonitor PrecTimer(*timerPrec_);
      return ( OPT::Apply( *RP_,x, y) );
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x, 
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::computeCurrPrecResVec( MV* R, const MV* X, const MV* B ) const {

    if (R) {
      if (X && B) // The entries are specified, so compute the residual of Op(A)X = B
	{
	  if (LP_!=null)
	    {
	      RCP<MV> R_temp = MVT::Clone( *X, MVT::GetNumberVecs( *X ) );
              {
                Teuchos::TimeMonitor OpTimer(*timerOp_);
	        OPT::Apply( *A_, *X, *R_temp );
              }
	      MVT::MvAddMv( -1.0, *R_temp, 1.0, *B, *R_temp );
              {
                Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	        OPT::Apply( *LP_, *R_temp, *R );
              }
	    }
	  else 
	    {
              {
                Teuchos::TimeMonitor OpTimer(*timerOp_);
	        OPT::Apply( *A_, *X, *R );
              }
	      MVT::MvAddMv( -1.0, *R, 1.0, *B, *R );
	    }
	}
      else { 
	// The solution and right-hand side may not be specified, check and use which ones exist.
	RCP<const MV> localB, localX;
	if (B)
	  localB = rcp( B, false );
	else
	  localB = curB_;
	
	if (X)
	  localX = rcp( X, false );
	else
	  localX = curX_;
	
	if (LP_!=null)
	  {
	    RCP<MV> R_temp = MVT::Clone( *localX, MVT::GetNumberVecs( *localX ) );
            {
              Teuchos::TimeMonitor OpTimer(*timerOp_);
	      OPT::Apply( *A_, *localX, *R_temp );
            }
	    MVT::MvAddMv( -1.0, *R_temp, 1.0, *localB, *R_temp );
            {
              Teuchos::TimeMonitor PrecTimer(*timerPrec_);
	      OPT::Apply( *LP_, *R_temp, *R );
            }
	  }
	else 
	  {
            {
              Teuchos::TimeMonitor OpTimer(*timerOp_);
  	      OPT::Apply( *A_, *localX, *R );
            }
	    MVT::MvAddMv( -1.0, *R, 1.0, *localB, *R );
	  }
      }    
    } 
  }
  
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::computeCurrResVec( MV* R, const MV* X, const MV* B ) const {

    if (R) {
      if (X && B) // The entries are specified, so compute the residual of Op(A)X = B
	{
          {
            Teuchos::TimeMonitor OpTimer(*timerOp_);
	    OPT::Apply( *A_, *X, *R );
          }
	  MVT::MvAddMv( -1.0, *R, 1.0, *B, *R );
	}
      else { 
	// The solution and right-hand side may not be specified, check and use which ones exist.
	RCP<const MV> localB, localX;
	if (B)
	  localB = rcp( B, false );
	else
	  localB = curB_;
	
	if (X)
	  localX = rcp( X, false );
	else
	  localX = curX_;
	  
        {
          Teuchos::TimeMonitor OpTimer(*timerOp_);
	  OPT::Apply( *A_, *localX, *R );
        }
	MVT::MvAddMv( -1.0, *R, 1.0, *localB, *R );
      }    
    }
  }
  
} // end Belos namespace

#endif /* BELOS_LINEAR_PROBLEM_HPP */


