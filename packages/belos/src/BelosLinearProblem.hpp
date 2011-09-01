
//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOS_LINEAR_PROBLEM_HPP
#define BELOS_LINEAR_PROBLEM_HPP

/*! \file BelosLinearProblem.hpp
    \brief Class which describes the linear problem to be solved by the iterative solver.
*/

#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

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
    LinearProblem(const Teuchos::RCP<const OP> &A, 
		  const Teuchos::RCP<MV> &X, 
		  const Teuchos::RCP<const MV> &B
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
    void setOperator(const Teuchos::RCP<const OP> &A) { A_ = A; isSet_=false; }
    
    //! Set left-hand-side X of linear problem AX = B.
    /*! Sets a pointer to a MultiVec.  No copy of the object is made.
     */
    void setLHS(const Teuchos::RCP<MV> &X) { X_ = X; isSet_=false; }
    
    //! Set right-hand-side B of linear problem AX = B.
    /*! Sets a pointer to a MultiVec.  No copy of the object is made.
     */
    void setRHS(const Teuchos::RCP<const MV> &B) { B_ = B; isSet_=false; }
    
    //! Set left preconditioning operator (\c LP) of linear problem AX = B.
    /*! Sets a pointer to an Operator.  No copy of the operator is made.
     */
    void setLeftPrec(const Teuchos::RCP<const OP> &LP) {  LP_ = LP; }
    
    //! Set right preconditioning operator (\c RP) of linear problem AX = B.
    /*! Sets a pointer to an Operator.  No copy of the operator is made.
     */
    void setRightPrec(const Teuchos::RCP<const OP> &RP) { RP_ = RP; }
    
    //! Inform the linear problem that the solver is finished with the current linear system.
    /*! \note This method is <b> only </b> to be used by the solver to inform the linear problem that it is
      finished with this block of linear systems.  The next time the Curr(RHS/LHS)Vec() is called, the next
      linear system will be returned.  Computing the next linear system isn't done in this method in case the 
      blocksize is changed.
    */
    void setCurrLS();

    /// \brief Tell the linear problem which linear system(s) need to be solved next.
    ///
    /// Any calls to get the current RHS/LHS vectors after this method
    /// is called will return the new linear system(s) indicated by \c
    /// index.  The length of \c index is assumed to be the blocksize.
    /// Entries of \c index must be between 0 and the number of
    /// vectors in the RHS/LHS multivector.  An entry of \c index may
    /// also be -1, which means this column of the linear system is
    /// augmented using a random vector.
    void setLSIndex(const std::vector<int>& index); 
    
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
    Teuchos::RCP<MV> updateSolution( const Teuchos::RCP<MV>& update = Teuchos::null,
				    bool updateLP = false,
                                    ScalarType scale = Teuchos::ScalarTraits<ScalarType>::one() );    

    //! Compute the new solution to the linear system given the /c update without updating the linear problem.
    Teuchos::RCP<MV> updateSolution( const Teuchos::RCP<MV>& update = Teuchos::null,
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
    bool setProblem( const Teuchos::RCP<MV> &newX = Teuchos::null, const Teuchos::RCP<const MV> &newB = Teuchos::null );

    //@}
    
    //! @name Accessor methods
    //@{ 
    
    //! Get a pointer to the operator A.
    Teuchos::RCP<const OP> getOperator() const { return(A_); }
    
    //! Get a pointer to the left-hand side X.
    Teuchos::RCP<MV> getLHS() const { return(X_); }
    
    //! Get a pointer to the right-hand side B.
    Teuchos::RCP<const MV> getRHS() const { return(B_); }
    
    //! Get a pointer to the initial residual vector.
    /*! \note This is the unpreconditioned residual.
     */
    Teuchos::RCP<const MV> getInitResVec() const { return(R0_); }
    
    //! Get a pointer to the preconditioned initial residual vector.
    /*! \note This is the preconditioned residual if the linear system is preconditioned on the left.
     */
    Teuchos::RCP<const MV> getInitPrecResVec() const { return(PR0_); }
    
    //! Get a pointer to the current left-hand side (solution) of the linear system.
    /*! This method is called by the solver or any method that is interested in the current linear system
      being solved for.  
      <ol>
      <li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
      <li> If there is no linear system to solve, this method will return a NULL pointer
      </ol>
    */
    Teuchos::RCP<MV> getCurrLHSVec();
    
    //! Get a pointer to the current right-hand side of the linear system.
    /*! This method is called by the solver of any method that is interested in the current linear system
      being solved for.  
      <ol>
      <li> If the solution has been updated by the solver, then this vector is current ( see SolutionUpdated() ).
      <li> If there is no linear system to solve, this method will return a NULL pointer
      </ol>
    */	
    Teuchos::RCP<const MV> getCurrRHSVec();
    
    //! Get a pointer to the left preconditioning operator.
    Teuchos::RCP<const OP> getLeftPrec() const { return(LP_); };
    
    //! Get a pointer to the right preconditioning operator.
    Teuchos::RCP<const OP> getRightPrec() const { return(RP_); };
    
    //! Get the 0-based index vector indicating the current linear systems being solved for.
    /*! Since the block size is independent of the number of right-hand sides for
      some solvers (GMRES, CG, etc.), it is important to know which linear systems
      are being solved for.  That may mean you need to update the information
      about the norms of your initial residual vector for weighting purposes.  This
      information can keep you from querying the solver for information that rarely
      changes.
      \note The length of the returned index vector is the number of right-hand sides being solved for.
            If an entry of the index vector is -1, then that linear system is an augmented linear
	    system and doesn't need to be considered for convergence.
      \note The index vector returned from this method is valid if isProblemSet() returns true.
    */
    const std::vector<int> getLSIndex() const { return(rhsIndex_); }

    //! Get the number of linear systems that have been set with this LinearProblem object.
    /* This can be used by status test classes to determine if the solver manager has advanced 
       and is solving another linear system.
    */
    int getLSNumber() const { return(lsNum_); }

    /*! \brief Return the timers for this object.
     *
     * The timers are ordered as follows:
     *   - time spent applying operator
     *   - time spent applying preconditioner
     */
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
      return tuple(timerOp_,timerPrec_);
    }


    //@}
    
    //! @name State methods
    //@{ 
    
    /// \brief Has the current approximate solution been updated?
    ///
    /// This only means that the current linear system for which the
    /// solver is solving (as obtained by getCurr{LHS, RHS}Vec()) has
    /// been updated by the solver.  This will be true every iteration
    /// for solvers like CG, but not true for solvers like GMRES until
    /// the solver restarts.
    bool isSolutionUpdated() const { return(solutionUpdated_); }

    //! Whether the problem has been set.
    bool isProblemSet() const { return(isSet_); }
    
    /// Whether the operator A is symmetric (in real arithmetic, or
    /// Hermitian in complex arithmetic).
    bool isHermitian() const { return(isHermitian_); }
    
    //! Whether the linear system is being preconditioned on the left.
    bool isLeftPrec() const { return(LP_!=Teuchos::null); }

    //! Whether the linear system is being preconditioned on the right.
    bool isRightPrec() const { return(RP_!=Teuchos::null); }
 
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
    Teuchos::RCP<const OP> A_;
    
    //! Solution vector of linear system.
    Teuchos::RCP<MV> X_;
    
    //! Current solution vector of the linear system.
    Teuchos::RCP<MV> curX_;
    
    //! Right-hand side of linear system.
    Teuchos::RCP<const MV> B_;
    
    //! Current right-hand side of the linear system.
    Teuchos::RCP<const MV> curB_;
    
    //! Initial residual of the linear system.
    Teuchos::RCP<MV> R0_;
   
    //! Preconditioned initial residual of the linear system.
    Teuchos::RCP<MV> PR0_;
 
    //! Left preconditioning operator of linear system
    Teuchos::RCP<const OP> LP_;  
    
    //! Right preconditioning operator of linear system
    Teuchos::RCP<const OP> RP_;
    
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
  LinearProblem<ScalarType,MV,OP>::LinearProblem(const Teuchos::RCP<const OP> &A, 
						 const Teuchos::RCP<MV> &X, 
						 const Teuchos::RCP<const MV> &B
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
  void LinearProblem<ScalarType,MV,OP>::setLSIndex(const std::vector<int>& index)
  {
    // Set new linear systems using the indices in index.
    rhsIndex_ = index;
    
    // Compute the new block linear system.
    // ( first clean up old linear system )
    curB_ = Teuchos::null;
    curX_ = Teuchos::null;
   
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
      Teuchos::RCP<MV> tmpCurB = MVT::Clone( *B_, blocksize_ );
      MVT::MvRandom(*tmpCurB);
      //
      // Now put in the part of B into curB 
      Teuchos::RCP<const MV> tptr = MVT::CloneView( *B_, vldIndex );
      MVT::SetBlock( *tptr, newIndex, *tmpCurB );
      curB_ = tmpCurB;
      //
      // Now put in the part of X into curX
      tptr = MVT::CloneView( *X_, vldIndex );
      MVT::SetBlock( *tptr, newIndex, *curX_ );
      //
      solutionUpdated_ = false;
    }
    else {
      curX_ = MVT::CloneViewNonConst( *X_, rhsIndex_ );
      curB_ = MVT::CloneView( *B_, rhsIndex_ );
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
      Teuchos::RCP<MV> tptr = MVT::CloneViewNonConst( *curX_, newIndex );
      MVT::SetBlock( *tptr, vldIndex, *X_ );
    }
    //
    // Clear the current vectors of this linear system so that any future calls
    // to get the vectors for this system return null pointers.
    //
    curX_ = Teuchos::null;
    curB_ = Teuchos::null;
    rhsIndex_.resize(0);
  }
  

  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<MV> LinearProblem<ScalarType,MV,OP>::updateSolution( const Teuchos::RCP<MV>& update, 
							   bool updateLP,
							   ScalarType scale )
  { 
    Teuchos::RCP<MV> newSoln;
    if (update != Teuchos::null) {
      if (updateLP == true) {
	if (RP_!=Teuchos::null) {
	  //
	  // Apply the right preconditioner before computing the current solution.
	  Teuchos::RCP<MV> TrueUpdate = MVT::Clone( *update, MVT::GetNumberVecs( *update ) );
	  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
	    Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
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
	if (RP_!=Teuchos::null) {
	  //
	  // Apply the right preconditioner before computing the current solution.
	  Teuchos::RCP<MV> trueUpdate = MVT::Clone( *update, MVT::GetNumberVecs( *update ) );
	  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
	    Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
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
  bool LinearProblem<ScalarType,MV,OP>::setProblem( const Teuchos::RCP<MV> &newX, const Teuchos::RCP<const MV> &newB )
  {
    // Create timers if the haven't been created yet.
    if (timerOp_ == Teuchos::null) {
      std::string opLabel = label_ + ": Operation Op*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerOp_ = Teuchos::TimeMonitor::getNewTimer( opLabel );
#endif
    }
    if (timerPrec_ == Teuchos::null) {
      std::string precLabel = label_ + ": Operation Prec*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerPrec_ = Teuchos::TimeMonitor::getNewTimer( precLabel );
#endif
    }

    // Set the linear system using the arguments newX and newB
    if (newX != Teuchos::null)
      X_ = newX;
    if (newB != Teuchos::null)
      B_ = newB;

    // Invalidate the current linear system indices and multivectors.
    rhsIndex_.resize(0);
    curX_ = Teuchos::null;
    curB_ = Teuchos::null;

    // If we didn't set a matrix A, a left-hand side X, or a
    // right-hand side B, then we didn't set the problem.
    if (A_ == Teuchos::null || X_ == Teuchos::null || B_ == Teuchos::null) {
      isSet_ = false;
      return isSet_;
    }

    // Reset whether the solution has been updated.  (We're just
    // setting the problem now, so of course we haven't updated the
    // solution yet.)
    solutionUpdated_ = false;
    
    // Compute the initial residuals.
    if (R0_==Teuchos::null || MVT::GetNumberVecs( *R0_ )!=MVT::GetNumberVecs( *B_ )) {
      R0_ = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );
    }
    computeCurrResVec( &*R0_, &*X_, &*B_ );

    if (LP_!=Teuchos::null) {
      if (PR0_==Teuchos::null || MVT::GetNumberVecs( *PR0_ )!=MVT::GetNumberVecs( *B_ )) {
        PR0_ = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );
      }
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
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
  Teuchos::RCP<MV> LinearProblem<ScalarType,MV,OP>::getCurrLHSVec()
  {
    if (isSet_) {
      return curX_;
    }
    else {
      return Teuchos::null;
    }
  }
  
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> LinearProblem<ScalarType,MV,OP>::getCurrRHSVec()
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
    Teuchos::RCP<MV> ytemp = MVT::Clone( y, MVT::GetNumberVecs( y ) );
    bool leftPrec = LP_!=Teuchos::null;
    bool rightPrec = RP_!=Teuchos::null;
    //
    // No preconditioning.
    // 
    if (!leftPrec && !rightPrec){ 
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
      OPT::Apply( *A_, x, y );
    }
    //
    // Preconditioning is being done on both sides
    //
    else if( leftPrec && rightPrec ) 
      {
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	  OPT::Apply( *RP_, x, y );   
        }
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	  OPT::Apply( *A_, y, *ytemp );
        }
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	  OPT::Apply( *LP_, *ytemp, y );
        }
      }
    //
    // Preconditioning is only being done on the left side
    //
    else if( leftPrec ) 
      {
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	  OPT::Apply( *A_, x, *ytemp );
        }
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	  OPT::Apply( *LP_, *ytemp, y );
        }
      }
    //
    // Preconditioning is only being done on the right side
    //
    else 
      {
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	  OPT::Apply( *RP_, x, *ytemp );
        }
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
      	  OPT::Apply( *A_, *ytemp, y );
        }
      }  
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyOp( const MV& x, MV& y ) const {
    if (A_.get()) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
      OPT::Apply( *A_,x, y);   
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x, 
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyLeftPrec( const MV& x, MV& y ) const {
    if (LP_!=Teuchos::null) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
      return ( OPT::Apply( *LP_,x, y) );
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x, 
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::applyRightPrec( const MV& x, MV& y ) const {
    if (RP_!=Teuchos::null) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
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
	  if (LP_!=Teuchos::null)
	    {
	      Teuchos::RCP<MV> R_temp = MVT::Clone( *B, MVT::GetNumberVecs( *B ) );
              {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	        OPT::Apply( *A_, *X, *R_temp );
              }
	      MVT::MvAddMv( -1.0, *R_temp, 1.0, *B, *R_temp );
              {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	        OPT::Apply( *LP_, *R_temp, *R );
              }
	    }
	  else 
	    {
              {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	        OPT::Apply( *A_, *X, *R );
              }
	      MVT::MvAddMv( -1.0, *R, 1.0, *B, *R );
	    }
	}
      else { 
	// The solution and right-hand side may not be specified, check and use which ones exist.
	Teuchos::RCP<const MV> localB, localX;
	if (B)
	  localB = Teuchos::rcp( B, false );
	else
	  localB = curB_;
	
	if (X)
	  localX = Teuchos::rcp( X, false );
	else
	  localX = curX_;
	
	if (LP_!=Teuchos::null)
	  {
	    Teuchos::RCP<MV> R_temp = MVT::Clone( *localB, MVT::GetNumberVecs( *localB ) );
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	      OPT::Apply( *A_, *localX, *R_temp );
            }
	    MVT::MvAddMv( -1.0, *R_temp, 1.0, *localB, *R_temp );
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
	      OPT::Apply( *LP_, *R_temp, *R );
            }
	  }
	else 
	  {
            {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
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
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	    OPT::Apply( *A_, *X, *R );
          }
	  MVT::MvAddMv( -1.0, *R, 1.0, *B, *R );
	}
      else { 
	// The solution and right-hand side may not be specified, check and use which ones exist.
	Teuchos::RCP<const MV> localB, localX;
	if (B)
	  localB = Teuchos::rcp( B, false );
	else
	  localB = curB_;
	
	if (X)
	  localX = Teuchos::rcp( X, false );
	else
	  localX = curX_;
	  
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
	  OPT::Apply( *A_, *localX, *R );
        }
	MVT::MvAddMv( -1.0, *R, 1.0, *localB, *R );
      }    
    }
  }
  
} // end Belos namespace

#endif /* BELOS_LINEAR_PROBLEM_HPP */


