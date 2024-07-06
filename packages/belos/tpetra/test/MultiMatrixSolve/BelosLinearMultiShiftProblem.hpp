
// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_LINEAR_MULTIPROBLEM_HPP
#define BELOS_LINEAR_MULTIPROBLEM_HPP

/// \file BelosLinearMultiShiftProblem.hpp
/// \brief Class which describes the linear problem to be solved by
///   the iterative solver.
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Belos {

  //@}

  /// \class LinearMultiShiftProblem
  /// \brief A linear system to solve multiple matrix and right-hand sides.
  ///
  /// This class encapsulates the general information needed for
  /// solving a set of linear systems of equations using an iterative method.
  /// A_i*x_i = b_i for i=1,...,N, where A_i is A+i*I
  ///
  /// \tparam ScalarType The type of the entries in the matrix and
  ///   vectors.
  /// \tparam MV The (multi)vector type.
  /// \tparam OP The operator type.  (Operators are functions that
  ///   take a multivector as input and compute a multivector as
  ///   output.)
  template <class ScalarType, class MV, class OP>
  class LinearMultiShiftProblem : public LinearProblem<ScalarType, MV, OP> {
  public:

    //! @name Constructors/Destructor
    //@{

    /// \brief Default constructor.
    ///
    /// Creates an empty Belos::LinearMultiShiftProblem instance.  The operator
    /// A, left-hand-side X and right-hand-side B must be set using
    /// the \c setOperator(), \c setLHS(), and \c setRHS() methods
    /// respectively.
    LinearMultiShiftProblem (void);
    /// \brief Unpreconditioned linear system constructor.
    ///
    /// Creates an unpreconditioned LinearMultiShiftProblem instance with the
    /// operator (\c A), initial guess (\c X), and right hand side (\c
    /// B).  Preconditioners can be set using the \c setLeftPrec() and
    /// \c setRightPrec() methods.
    LinearMultiShiftProblem (const Teuchos::RCP<const OP> &A,
		   const Teuchos::RCP<MV> &X,
		   const Teuchos::RCP<const MV> &B);

    //! Destructor (declared virtual for memory safety of derived classes).
    virtual ~LinearMultiShiftProblem (void) {}

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


    void
    setShift( bool addShift )
    {
      addShift_ = addShift;
    }

    //@}

    //! @name Apply / Compute methods
    //@{

    /// \brief Apply ONLY the operator to \c x, returning \c y.
    ///
    /// This method only applies the linear problem's operator,
    /// without any preconditioners that may have been defined.
    /// Flexible variants of Krylov methods will use this method.  If
    /// no operator has been defined, this method just copies x into
    /// y.
    void applyOp( const MV& x, MV& y ) const;

    //@}

  private:
    bool addShift_;

    using LinearProblem<ScalarType, MV, OP>::A_;
    using LinearProblem<ScalarType, MV, OP>::X_;
    using LinearProblem<ScalarType, MV, OP>::curX_;
    using LinearProblem<ScalarType, MV, OP>::B_;
    using LinearProblem<ScalarType, MV, OP>::curB_;
    using LinearProblem<ScalarType, MV, OP>::R0_;
    using LinearProblem<ScalarType, MV, OP>::PR0_;
    using LinearProblem<ScalarType, MV, OP>::R0_user_;
    using LinearProblem<ScalarType, MV, OP>::PR0_user_;
    using LinearProblem<ScalarType, MV, OP>::LP_;
    using LinearProblem<ScalarType, MV, OP>::RP_;
    using LinearProblem<ScalarType, MV, OP>::timerOp_;
    using LinearProblem<ScalarType, MV, OP>::timerPrec_;
    using LinearProblem<ScalarType, MV, OP>::blocksize_;
    using LinearProblem<ScalarType, MV, OP>::num2Solve_;
    using LinearProblem<ScalarType, MV, OP>::rhsIndex_;
    using LinearProblem<ScalarType, MV, OP>::lsNum_;
    using LinearProblem<ScalarType, MV, OP>::isSet_;
    using LinearProblem<ScalarType, MV, OP>::solutionUpdated_;
    using LinearProblem<ScalarType, MV, OP>::label_;

    typedef MultiVecTraits<ScalarType,MV>  MVT;
    typedef OperatorTraits<ScalarType,MV,OP>  OPT;
  };

  //--------------------------------------------
  //  Constructor Implementations
  //--------------------------------------------

  template <class ScalarType, class MV, class OP>
  LinearMultiShiftProblem<ScalarType,MV,OP>::LinearMultiShiftProblem(void)
  : LinearProblem<ScalarType,MV,OP>(),
    addShift_(true)
  {}

  template <class ScalarType, class MV, class OP>
  LinearMultiShiftProblem<ScalarType,MV,OP>::LinearMultiShiftProblem(const Teuchos::RCP<const OP> &A, 
						 const Teuchos::RCP<MV> &X,
						 const Teuchos::RCP<const MV> &B
						 )
  : LinearProblem<ScalarType,MV,OP>(A, X, B),
    addShift_(true)
  {
    int numVecs = MVT::GetNumberVecs( *B );
    rhsIndex_.resize(numVecs);
    for (int i=0; i<numVecs; ++i)
      rhsIndex_[i] = i;
  }

  template <class ScalarType, class MV, class OP>
  bool
  LinearMultiShiftProblem<ScalarType,MV,OP>::
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
    curX_ = Teuchos::null;
    curB_ = Teuchos::null;

    // Temporarily set the rhsIndex_ to all the vectors of the linear system.
    // This is necessary for the applyOp to perform correctly in computing the 
    // initial residual.  It will be reset at the end of this method.
    int numVecs = MVT::GetNumberVecs( *B_ );
    rhsIndex_.resize(numVecs);
    for (int i=0; i<numVecs; ++i)
      rhsIndex_[i] = i;

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
    if(Teuchos::is_null(R0_user_)) {
      if (R0_==Teuchos::null || MVT::GetNumberVecs( *R0_ )!=MVT::GetNumberVecs( *B_ )) {
        R0_ = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );
      }
      this->computeCurrResVec( &*R0_, &*X_, &*B_ );

      if (LP_!=Teuchos::null) {
        if (PR0_==Teuchos::null || (PR0_==R0_) || (MVT::GetNumberVecs(*PR0_)!=MVT::GetNumberVecs(*B_))) {
          PR0_ = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );
        }
        this->applyLeftPrec( *R0_, *PR0_ );
      }
      else {
        PR0_ = R0_;
      }
    }
    else { // User specified the residuals
      // If the user did not specify the right sized residual, create one and set R0_user_ to null.
      if (MVT::GetNumberVecs( *R0_user_ )!=MVT::GetNumberVecs( *B_ )) {
        Teuchos::RCP<MV> helper = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );
        this->computeCurrResVec( &*helper, &*X_, &*B_ );
        R0_user_ = Teuchos::null;
        R0_ = helper;
      }

      if (LP_!=Teuchos::null) {
        // If the user provided preconditioned residual is the wrong size or pointing at
        // the wrong object, create one and set the PR0_user_ to null.
        if (PR0_user_==Teuchos::null || (PR0_user_==R0_) || (PR0_user_==R0_user_)
          || (MVT::GetNumberVecs(*PR0_user_)!=MVT::GetNumberVecs(*B_))) {
          Teuchos::RCP<MV> helper = MVT::Clone( *B_, MVT::GetNumberVecs( *B_ ) );

          // Get the initial residual from getInitResVec because R0_user_ may be null from above.
          this->applyLeftPrec( *(this->getInitResVec()), *helper );
          PR0_user_ = Teuchos::null;
          PR0_ = helper;
        }
      }
      else {
        // The preconditioned initial residual vector is the residual vector.
        // NOTE:  R0_user_ could be null if incompatible.
        if (R0_user_!=Teuchos::null)
        {
          PR0_user_ = R0_user_;
        }
        else
        {
          PR0_user_ = Teuchos::null;
          PR0_ = R0_;
        }
      }
    }

    // Resize the rhsIndex_ back to zero length to correspond with null pointers for curX_ and curB_.
    rhsIndex_.resize(0);

    // The problem has been set and is ready for use.
    isSet_ = true;
    // Return isSet.
    return isSet_;
  }

  template <class ScalarType, class MV, class OP>
  void LinearMultiShiftProblem<ScalarType,MV,OP>::applyOp( const MV& x, MV& y ) const {

    if (A_.get())
    {
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor OpTimer(*timerOp_);
#endif
      OPT::Apply( *A_,x, y);
      }
      // Now add in shift based on rhsIndex_
      if (addShift_)
      {
        std::vector<int> newIdx( 1 );
        for (unsigned int i=0; i<rhsIndex_.size(); ++i)
        {
          newIdx[0] = i;
          Teuchos::RCP<const MV> tmp_x = MVT::CloneView( x, newIdx );
          Teuchos::RCP<MV> tmp_y = MVT::CloneViewNonConst( y, newIdx );
          MVT::MvAddMv( rhsIndex_[i]*Teuchos::ScalarTraits<ScalarType>::one(), *tmp_x,
                        Teuchos::ScalarTraits<ScalarType>::one(), *tmp_y, *tmp_y );
        }
      }
    }
    else {
      MVT::MvAddMv( Teuchos::ScalarTraits<ScalarType>::one(), x,
		    Teuchos::ScalarTraits<ScalarType>::zero(), x, y );
    }
  }

} // end Belos namespace

#endif /* BELOS_LINEAR_MULTIPROBLEM_HPP */


