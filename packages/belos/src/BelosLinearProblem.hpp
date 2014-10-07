
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

/// \file BelosLinearProblem.hpp 
/// \brief Class which describes the linear problem to be solved by
///   the iterative solver.
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Belos {

  //! @name LinearProblem Exceptions
  //@{
  
  /// \class LinearProblemError
  /// \brief Exception thrown to signal error with the Belos::LinearProblem object.
  class LinearProblemError : public BelosError {
  public: 
    LinearProblemError (const std::string& what_arg) : 
      BelosError(what_arg) {}
  };
  
  //@}

  /// \class LinearProblem  
  /// \brief A linear system to solve, and its associated information.
  ///
  /// This class encapsulates the general information needed for
  /// solving a linear system of equations using an iterative method.
  ///
  /// \tparam ScalarType The type of the entries in the matrix and
  ///   vectors.
  /// \tparam MV The (multi)vector type.
  /// \tparam OP The operator type.  (Operators are functions that
  ///   take a multivector as input and compute a multivector as
  ///   output.)
  template <class ScalarType, class MV, class OP>
  class LinearProblem {
  public:
    
    //! @name Constructors/Destructor
    //@{ 

    /// \brief Default constructor.
    ///
    /// Creates an empty Belos::LinearProblem instance.  The operator
    /// A, left-hand-side X and right-hand-side B must be set using
    /// the \c setOperator(), \c setLHS(), and \c setRHS() methods
    /// respectively.
    LinearProblem (void);
    
    /// \brief Unpreconditioned linear system constructor.
    ///
    /// Creates an unpreconditioned LinearProblem instance with the
    /// operator (\c A), initial guess (\c X), and right hand side (\c
    /// B).  Preconditioners can be set using the \c setLeftPrec() and
    /// \c setRightPrec() methods, and scaling can also be set using
    /// the \c setLeftScale() and \c setRightScale() methods.
    LinearProblem (const Teuchos::RCP<const OP> &A, 
		   const Teuchos::RCP<MV> &X, 
		   const Teuchos::RCP<const MV> &B);
    
    /// \brief Copy constructor.
    ///
    /// Makes a copy of an existing LinearProblem instance.
    LinearProblem (const LinearProblem<ScalarType,MV,OP>& Problem);
    
    //! Destructor (declared virtual for memory safety of derived classes).
    virtual ~LinearProblem (void);

    //@}
    
    //! @name Set methods
    //@{ 
    
    /// \brief Set the operator A of the linear problem \f$AX = B\f$.
    /// 
    /// The operator is set by pointer; no copy of the operator is made.
    void setOperator (const Teuchos::RCP<const OP> &A) { 
      A_ = A; 
      isSet_=false; 
    }
    
    /// \brief Set left-hand-side X of linear problem \f$AX = B\f$.
    ///
    /// Setting the "left-hand side" sets the starting vector (also
    /// called "initial guess") of an iterative method.  The
    /// multivector is set by pointer; no copy of the object is made.
    /// Belos' solvers will modify this multivector in place.
    void setLHS (const Teuchos::RCP<MV> &X) { 
      X_ = X; 
      isSet_=false; 
    }
    
    /// \brief Set right-hand-side B of linear problem \f$AX = B\f$.
    ///
    /// The multivector is set by pointer; no copy of the object is
    /// made.
    void setRHS (const Teuchos::RCP<const MV> &B) { 
      B_ = B; 
      isSet_=false; 
    }
    
    /// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
    ///
    /// The operator is set by pointer; no copy of the operator is made.
    void setLeftPrec(const Teuchos::RCP<const OP> &LP) {  LP_ = LP; }

    /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
    ///
    /// The operator is set by pointer; no copy of the operator is made.
    void setRightPrec(const Teuchos::RCP<const OP> &RP) { RP_ = RP; }

    /// Tell the linear problem that the solver is finished with the current linear system.
    ///
    /// \note This method is <b>only</b> to be used by the solver to
    ///   inform the linear problem that it is finished with the
    ///   current block of linear systems.  The next time that
    ///   Curr{RHS, LHS}Vec() is called, the next linear system will
    ///   be returned.  Computing the next linear system isn't done in
    ///   this method in case the block size is changed.
    void setCurrLS ();

    /// \brief Tell the linear problem which linear system(s) need to be solved next.
    ///
    /// Any calls to get the current RHS/LHS vectors after this method
    /// is called will return the new linear system(s) indicated by \c
    /// index.  The length of \c index is assumed to be the blocksize.
    /// Entries of \c index must be between 0 and the number of
    /// vectors in the RHS/LHS multivector.  An entry of \c index may
    /// also be -1, which means this column of the linear system is
    /// augmented using a random vector.
    void setLSIndex (const std::vector<int>& index); 
    
    /// \brief Tell the linear problem that the (preconditioned) operator is Hermitian.
    ///
    /// This knowledge may allow the operator to take advantage of the
    /// linear problem's symmetry.  However, this method should not be
    /// called if the preconditioned operator is not Hermitian (or
    /// symmetric in real arithmetic).
    ///
    /// We make no attempt to detect the symmetry of the operators, so
    /// we cannot check whether this method has been called
    /// incorrectly.
    void setHermitian( bool isSym = true ) { isHermitian_ = isSym; }
   
    /// \brief Set the label prefix used by the timers in this object.  
    ///
    /// The default label prefix is "Belos".  The timers are created
    /// during the call to \c setProblem().  If they have already been
    /// created and this label is different than the current one, then
    /// this method will generate a new timer.
    void setLabel (const std::string& label);

    /// \brief Compute the new solution to the linear system using the
    ///   given update vector.
    ///
    /// Let \f$\delta\f$ be the update vector, \f$\alpha\f$ the scale
    /// factor, and \f$x\f$ the current solution.  If there is a right
    /// preconditioner \f$M_R^{-1}\f$, then we compute the new
    /// solution as \f$x + \alpha M_R^{-1} \delta\f$.  Otherwise, if
    /// there is no right preconditioner, we compute the new solution
    /// as \f$x + \alpha \delta\f$.
    ///
    /// This method always returns the new solution.  If updateLP is
    /// false, it computes the new solution as a deep copy, without
    /// modifying the internally stored current solution.  If updateLP
    /// is true, it computes the new solution in place, and returns a
    /// pointer to the internally stored solution.
    ///
    /// \param update [in/out] The solution update vector.  If null,
    ///   this method returns a pointer to the new solution.
    ///
    /// \param updateLP [in] This is ignored if the update vector is
    ///   null.  Otherwise, if updateLP is true, the following things
    ///   happen: (a) this LinearProblem's stored solution is updated
    ///   in place, and (b) the next time \c GetCurrResVecs() is
    ///   called, a new residual will be computed.  If updateLP is
    ///   false, then the new solution is computed and returned as a
    ///   copy, without modifying this LinearProblem's stored current
    ///   solution.
    ///
    /// \param scale [in] The factor \f$\alpha\f$ by which to multiply
    ///   the solution update vector when computing the update.  This
    ///   is ignored if the update vector is null.
    ///
    /// \return A pointer to the new solution.  This is freshly
    ///   allocated if updateLP is false, otherwise it is a view of
    ///   the LinearProblem's stored current solution.
    ///
    Teuchos::RCP<MV> 
    updateSolution (const Teuchos::RCP<MV>& update = Teuchos::null,
		    bool updateLP = false,
		    ScalarType scale = Teuchos::ScalarTraits<ScalarType>::one());    

    /// \brief Compute the new solution to the linear system using the
    ///   given update vector.
    ///
    /// This method does the same thing as calling the three-argument
    /// version of updateSolution() with updateLP = false.  It does
    /// not update the linear problem or change the linear problem's
    /// state in any way.
    ///
    /// \param update [in/out] The solution update vector.  If null,
    ///   this method returns a pointer to the new solution.
    ///
    /// \param scale [in] The factor \f$\alpha\f$ by which to multiply
    ///   the solution update vector when computing the update.  This
    ///   is ignored if the update vector is null.
    ///
    /// \return A pointer to the new solution. 
    ///
    Teuchos::RCP<MV> updateSolution( const Teuchos::RCP<MV>& update = Teuchos::null,
                                    ScalarType scale = Teuchos::ScalarTraits<ScalarType>::one() ) const
    { return const_cast<LinearProblem<ScalarType,MV,OP> *>(this)->updateSolution( update, false, scale ); }

    //@}
    
    //! @name Set / Reset method
    //@{ 
    
    /// \brief Set up the linear problem manager.
    ///
    /// Call this method if you want to solve the linear system with a
    /// different left- or right-hand side, or if you want to prepare
    /// the linear problem to solve the linear system that was already
    /// passed in.  (In the latter case, call this method with the
    /// default arguments.)  The internal flags will be set as if the
    /// linear system manager was just initialized, and the initial
    /// residual will be computed.
    ///
    /// Many of Belos' solvers require that this method has been
    /// called on the linear problem, before they can solve it.
    ///
    /// \param newX [in/out] If you want to solve the linear system
    ///   with a different left-hand side, pass it in here.
    ///   Otherwise, set this to null (the default value).
    ///
    /// \param newB [in] If you want to solve the linear system with a
    ///   different right-hand side, pass it in here.  Otherwise, set
    ///   this to null (the default value).
    ///
    /// \return Whether the linear problem was successfully set up.
    ///   Successful setup requires at least that the matrix operator
    ///   A, the left-hand side X, and the right-hand side B all be
    ///   non-null.
    bool 
    setProblem (const Teuchos::RCP<MV> &newX = Teuchos::null, 
		const Teuchos::RCP<const MV> &newB = Teuchos::null);

    //@}
    
    //! @name Accessor methods
    //@{ 
    
    //! A pointer to the (unpreconditioned) operator A.
    Teuchos::RCP<const OP> getOperator() const { return(A_); }
    
    //! A pointer to the left-hand side X.
    Teuchos::RCP<MV> getLHS() const { return(X_); }
    
    //! A pointer to the right-hand side B.
    Teuchos::RCP<const MV> getRHS() const { return(B_); }
    
    //! A pointer to the initial unpreconditioned residual vector.
    Teuchos::RCP<const MV> getInitResVec() const { return(R0_); }
    
    /// \brief A pointer to the preconditioned initial residual vector.
    ///
    /// \note This is the preconditioned residual if the linear system
    ///   is left preconditioned.
    Teuchos::RCP<const MV> getInitPrecResVec() const { return(PR0_); }
    
    /// \brief Get a pointer to the current left-hand side (solution) of the linear system.
    ///
    /// This method is called by the solver or any method that is
    /// interested in the current linear system being solved.
    ///   - If the solution has been updated by the solver, then this
    ///     vector is current ( see \c isSolutionUpdated() ).
    ///   - If there is no linear system to solve, this method returns
    ///     a null pointer.
    ///
    /// This method is <i>not</i> the same thing as \c getLHS().  The
    /// \c getLHS() method just returns a pointer to the original
    /// left-hand side vector.  This method only returns a valid
    /// vector if the current subset of right-hand side(s) to solve
    /// has been set (via the \c setLSIndex() method).
    Teuchos::RCP<MV> getCurrLHSVec();
    
    /// \brief Get a pointer to the current right-hand side of the linear system.
    ///
    /// This method is called by the solver or any method that is
    /// interested in the current linear system being solved.
    ///   - If the solution has been updated by the solver, then this
    ///     vector is current ( see \c isSolutionUpdated() ).
    ///   - If there is no linear system to solve, this method returns
    ///     a null pointer.
    ///
    /// This method is <i>not</i> the same thing as \c getRHS().  The
    /// \c getRHS() method just returns a pointer to the original
    /// right-hand side vector.  This method only returns a valid
    /// vector if the current subset of right-hand side(s) to solve
    /// has been set (via the \c setLSIndex() method).
    Teuchos::RCP<const MV> getCurrRHSVec();
    
    //! Get a pointer to the left preconditioner.
    Teuchos::RCP<const OP> getLeftPrec() const { return(LP_); };
    
    //! Get a pointer to the right preconditioner.
    Teuchos::RCP<const OP> getRightPrec() const { return(RP_); };
    
    /// \brief (Zero-based) indices of the linear system(s) currently being solved.
    ///
    /// Since the block size is independent of the number of
    /// right-hand sides for some solvers (GMRES, CG, etc.), it is
    /// important to know which linear systems are being solved.  That
    /// may mean you need to update the information about the norms of
    /// your initial residual vector for weighting purposes.  This
    /// information can help you avoid querying the solver for
    /// information that rarely changes.
    ///
    /// \note The length of the returned index vector is the number of
    ///   right-hand sides currently being solved.  If an entry of the
    ///   index vector is -1, then the corresponding linear system is
    ///   an augmented linear system and doesn't need to be considered
    ///   for convergence.
    /// 
    /// \note The index vector returned from this method can only be
    ///   nonempty if \c setLSIndex() has been called with a nonempty
    ///   index vector argument, or if this linear problem was
    ///   constructed via the copy constructor of a linear problem
    ///   with a nonempty index vector.
    const std::vector<int> getLSIndex() const { return(rhsIndex_); }

    /// \brief The number of linear systems that have been set.
    ///
    /// This can be used by status test classes to determine if the
    /// solver manager has advanced and is solving another linear
    /// system.  This is incremented by one every time that \c
    /// setLSIndex() completes successfully.
    int getLSNumber() const { return(lsNum_); }

    /*! \brief The timers for this object.
     *
     * The timers are ordered as follows:
     *   - time spent applying operator
     *   - time spent applying preconditioner
     */
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
      return Teuchos::tuple(timerOp_,timerPrec_);
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

    /// \brief Whether the (preconditioned) operator is Hermitian.
    ///
    /// If preconditioner(s) are defined and this method returns true,
    /// then the entire preconditioned operator is Hermitian (or
    /// symmetric in real arithmetic).
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
    
    /// \brief Apply ONLY the operator to \c x, returning \c y.
    ///
    /// This method only applies the linear problem's operator,
    /// without any preconditioners that may have been defined.
    /// Flexible variants of Krylov methods will use this method.  If
    /// no operator has been defined, this method just copies x into
    /// y.
    void applyOp( const MV& x, MV& y ) const;
    
    /// \brief Apply ONLY the left preconditioner to \c x, returning \c y.
    ///
    /// This method only applies the left preconditioner.  This may be
    /// required for flexible variants of Krylov methods.  If no left
    /// preconditioner has been defined, this method just copies x
    /// into y.
    void applyLeftPrec( const MV& x, MV& y ) const;

    /// \brief Apply ONLY the right preconditioner to \c x, returning \c y.
    ///
    /// This method only applies the right preconditioner.  This may
    /// be required for flexible variants of Krylov methods.  If no
    /// right preconditioner has been defined, this method just copies
    /// x into y.
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

    //! @name Booleans to keep track of linear problem attributes and status.
    //@{ 

    //! Is there a left scaling?
    bool Left_Scale_;

    //! Is there a right scaling?
    bool Right_Scale_;

    //! Has the linear problem to solve been set?
    bool isSet_;

    /// Whether the operator A is symmetric (in real arithmetic, or
    /// Hermitian in complex arithmetic).
    bool isHermitian_;

    //! Has the current approximate solution been updated?
    bool solutionUpdated_;    

    //@}
   
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
  Teuchos::RCP<MV> 
  LinearProblem<ScalarType,MV,OP>::
  updateSolution (const Teuchos::RCP<MV>& update, 
		  bool updateLP,
		  ScalarType scale)
  { 
    using Teuchos::RCP;
    using Teuchos::null;

    RCP<MV> newSoln;
    if (update.is_null())
      { // The caller didn't supply an update vector, so we assume
	// that the current solution curX_ is unchanged, and return a
	// pointer to it.
	newSoln = curX_;
      }
    else // the update vector is NOT null.
      { 
	if (updateLP) // The caller wants us to update the linear problem.
	  { 
	    if (RP_.is_null()) 
	      { // There is no right preconditioner.
		// curX_ := curX_ + scale * update.
		MVT::MvAddMv( 1.0, *curX_, scale, *update, *curX_ ); 
	      }
	    else 
	      { // There is indeed a right preconditioner, so apply it
		// before computing the new solution.
		RCP<MV> rightPrecUpdate = 
		  MVT::Clone (*update, MVT::GetNumberVecs (*update));
		{
#ifdef BELOS_TEUCHOS_TIME_MONITOR
		  Teuchos::TimeMonitor PrecTimer (*timerPrec_);
#endif
		  OPT::Apply( *RP_, *update, *rightPrecUpdate ); 
		}
		// curX_ := curX_ + scale * rightPrecUpdate.
		MVT::MvAddMv( 1.0, *curX_, scale, *rightPrecUpdate, *curX_ ); 
	      } 
	    solutionUpdated_ = true; 
	    newSoln = curX_;
	  }
	else 
	  { // Rather than updating our stored current solution curX_,
	    // we make a copy and compute the new solution in the
	    // copy, without modifying curX_.
	    newSoln = MVT::Clone (*update, MVT::GetNumberVecs (*update));
	    if (RP_.is_null())
	      { // There is no right preconditioner.
		// newSoln := curX_ + scale * update.
		MVT::MvAddMv( 1.0, *curX_, scale, *update, *newSoln ); 
	      }
	    else 
	      { // There is indeed a right preconditioner, so apply it
		// before computing the new solution.
		RCP<MV> rightPrecUpdate = 
		  MVT::Clone (*update, MVT::GetNumberVecs (*update));
		{
#ifdef BELOS_TEUCHOS_TIME_MONITOR
		  Teuchos::TimeMonitor PrecTimer(*timerPrec_);
#endif
		  OPT::Apply( *RP_, *update, *rightPrecUpdate ); 
		}
		// newSoln := curX_ + scale * rightPrecUpdate.
		MVT::MvAddMv( 1.0, *curX_, scale, *rightPrecUpdate, *newSoln ); 
	      } 
	  }
      }
    return newSoln;
  }
  
  template <class ScalarType, class MV, class OP>
  void LinearProblem<ScalarType,MV,OP>::setLabel(const std::string& label)
  {
    if (label != label_) {
      label_ = label;
      // Create new timers if they have already been created.
      if (timerOp_ != Teuchos::null) {
        std::string opLabel = label_ + ": Operation Op*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        timerOp_ = Teuchos::TimeMonitor::getNewCounter( opLabel );
#endif
      }
      if (timerPrec_ != Teuchos::null) {
        std::string precLabel = label_ + ": Operation Prec*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        timerPrec_ = Teuchos::TimeMonitor::getNewCounter( precLabel );
#endif
      }
    }
  }

  template <class ScalarType, class MV, class OP>
  bool 
  LinearProblem<ScalarType,MV,OP>::
  setProblem (const Teuchos::RCP<MV> &newX, 
	      const Teuchos::RCP<const MV> &newB)
  {
    // Create timers if the haven't been created yet.
    if (timerOp_ == Teuchos::null) {
      std::string opLabel = label_ + ": Operation Op*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerOp_ = Teuchos::TimeMonitor::getNewCounter( opLabel );
#endif
    }
    if (timerPrec_ == Teuchos::null) {
      std::string precLabel = label_ + ": Operation Prec*x";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerPrec_ = Teuchos::TimeMonitor::getNewCounter( precLabel );
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
    using Teuchos::null;
    using Teuchos::RCP;

    const bool leftPrec = LP_ != null;
    const bool rightPrec = RP_ != null;

    // We only need a temporary vector for intermediate results if
    // there is a left or right preconditioner.  We really should just
    // keep this temporary vector around instead of allocating it each
    // time.
    RCP<MV> ytemp = (leftPrec || rightPrec) ? MVT::Clone (y, MVT::GetNumberVecs (y)) : null;

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


