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

#ifndef BELOS_GCR_ITER_HPP
#define BELOS_GCR_ITER_HPP

/*! \file BelosGCRIter.hpp
    \brief Belos concrete class for performing the pseudo-block GCR iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosCGIteration.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!
  \class Belos::GCRIter

  \brief This class implements the pseudo-block GCR iteration, where the basic GCR
  algorithm is performed on all of the linear systems sequentially. The algorithm
  is based on the GCR algorithm described in the "Iterative Methods for Sparse
  Linear Systems" (by Yousef Saad)

  \ingroup belos_solver_framework

  \author Vinh Dang
*/

namespace Belos {

  //! @name GCRIteration Structures
  //@{

  /** \brief Structure to contain pointers to GCRIteration state variables.
   *
   * This struct is utilized by GCRIteration::initialize() and GCRIteration::getState().
   */
  template <class ScalarType, class MV>
  struct GCRIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> R;

    //std::vector<ScalarType> rho_old, alpha, omega;

    GCRIterationState() : R(Teuchos::null)
    {}
  };

  template<class ScalarType, class MV, class OP>
  class GCRIter : virtual public Iteration<ScalarType,MV,OP> {

  public:

    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

    //! @name Constructors/Destructor
    //@{

    /*! \brief %GCRIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver.
     */
    GCRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
             const Teuchos::RCP<OutputManager<ScalarType> > &printer,
             const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
             Teuchos::ParameterList &params );

    //! Destructor.
    virtual ~GCRIter() {};
    //@}


    //! @name Solver methods
    //@{

    /*! \brief This method performs GCR iterations on each linear system until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * std::exception is thrown).
     *
     * iterate() will first determine whether the solver is initialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs GCR iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller.
     *
     * The status test is queried at the beginning of the iteration.
     *
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %GCRIter contains a certain amount of state, consisting of the current
     * direction vectors and residuals.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules given
     * below.
     *
     * \post
     * <li>isInitialized() == \c true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any necessary component of the
     * state not given to initialize() will be generated.
     *
     * \note For any pointer in \c newstate which directly points to the multivectors in
     * the solver, the data is not copied.
     */
    void initializeGCR(GCRIterationState<ScalarType,MV>& newstate);

    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      GCRIterationState<ScalarType,MV> empty;
      initializeGCR(empty);
    }

    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A GCRIterationState object containing const pointers to the current
     * solver state.
     */
    GCRIterationState<ScalarType,MV> getState() const {
      GCRIterationState<ScalarType,MV> state;
      state.R = R_;
      return state;
    }

    //@}


    //! @name Status methods
    //@{

    //! \brief Get the current iteration count.
    int getNumIters() const { return iter_; }

    //! \brief Reset the iteration count.
    void resetNumIters( int iter = 0 ) { iter_ = iter; }

    //! Get the norms of the residuals native to the solver.
    //! \return A std::vector of length blockSize containing the native residuals.
    // amk TODO: are the residuals actually being set?  What is a native residual?
    Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> * /* norms */ ) const { return R_; }

    //! Get the current update to the linear system.
    /*! \note This method returns a null pointer because the linear problem is current.
    */
    // amk TODO: what is this supposed to be doing?
    Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

    //! Has breakdown been detected in any linear system.
    bool breakdownDetected() { return breakdown_; }

    //@}

    //! @name Accessor methods
    //@{

    //! Get a constant reference to the linear problem.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }

    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
                         "Belos::GCRIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! Get the maximum number of vectors in the Krylov subspace used by the iterative solver in solving this linear problem.
    int getNumKrylovVecs() const { return numKrylovVecs_; }
    
    //! \brief Set the maximum number of vectors in the Krylov subspace used by the iterative solver.
    void setNumKrylovVecs(int numKrylovVecs);

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }

    //@}

  private:

    void axpy(const ScalarType alpha, const MV & A,
              const std::vector<ScalarType> beta, const MV& B, MV& mv, bool minus=false);

    void deep_copy(const MV& Src, const std::vector<int>& colIndex,
                   const std::vector<Teuchos::RCP<MV> >& Dst);

    //
    // Classes inputed through constructor that define the linear problem to be solved.
    //
    const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >  lp_;
    const Teuchos::RCP<OutputManager<ScalarType> >        om_;
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >     stest_;

    //
    // Algorithmic parameters
    //
    // numRHS_ is the current number of linear systems being solved.
    int numRHS_;
    // numKrylovVecs_ is the number of vectors in the Krylov subspace.
    int numKrylovVecs_;

    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;

    // Breakdown has been observed for at least one of the linear systems
    bool breakdown_;

    // Current number of iterations performed.
    int iter_;

    //
    // State Storage
    //
    // Residual and temporary multivecs
    Teuchos::RCP<MV> R_, AxR_;
    //
    std::vector<ScalarType> pone_;
    std::vector<int> curIndex_, newIndex_;
    //
    std::vector<Teuchos::RCP<MV> > U_, C_;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  GCRIter<ScalarType,MV,OP>::GCRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                     const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                     const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                     Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numRHS_(0),
    numKrylovVecs_(0),
    initialized_(false),
    breakdown_(false),
    iter_(0)
  {
    // Get the maximum number of vectors allowed in the Krylov subspace
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Num KrylovVecs"), std::invalid_argument,
                       "Belos::GCRIter::constructor: mandatory parameter 'Num KrylovVecs' is not specified.");
    int nKrylovVecs = Teuchos::getParameter<int>(params, "Num KrylovVecs");

    setNumKrylovVecs( nKrylovVecs );
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the maximum number of vectors allowed in the Krylov subspace and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void GCRIter<ScalarType,MV,OP>::setNumKrylovVecs (int numKrylovVecs)
  { 
    TEUCHOS_TEST_FOR_EXCEPTION(numKrylovVecs <= 0, std::invalid_argument, "Belos::GCRIter::setNumKrylovVecs was passed a non-positive argument.");

    numKrylovVecs_ = numKrylovVecs;

    initialized_ = false;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void GCRIter<ScalarType,MV,OP>::initializeGCR(GCRIterationState<ScalarType,MV>& newstate)
  {
    // Check if there is any multivector to clone from.
    Teuchos::RCP<const MV> lhsMV = lp_->getCurrLHSVec();
    Teuchos::RCP<const MV> rhsMV = lp_->getCurrRHSVec();
    TEUCHOS_TEST_FOR_EXCEPTION((lhsMV==Teuchos::null && rhsMV==Teuchos::null),std::invalid_argument,
                       "Belos::GCRIter::initialize(): Cannot initialize state storage!");

    // Get the multivector that is not null.
    Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );

    // Get the number of right-hand sides we're solving for now.
    int numRHS = MVT::GetNumberVecs(*tmp);
    numRHS_ = numRHS;

    // Initialize the state storage
    if (Teuchos::is_null(R_) || MVT::GetNumberVecs(*R_)!=numRHS_) {
      R_   = MVT::Clone( *tmp, numRHS_ );
      AxR_ = MVT::Clone( *tmp, numRHS_ );

      MVT::MvInit(*R_);
      MVT::MvInit(*AxR_);

      pone_.resize(numRHS_);
      curIndex_.resize(numRHS_);
      newIndex_.resize(numRHS_);
    }
    U_.resize(numRHS_);
    C_.resize(numRHS_);
    for (int i=0; i<numRHS_; ++i) {
      if (Teuchos::is_null(U_[i])) {
        U_[i] = MVT::Clone(*tmp, numKrylovVecs_);
        MVT::MvInit(*U_[i]);
      }
      if (Teuchos::is_null(C_[i])) {
        C_[i] = MVT::Clone(*tmp, numKrylovVecs_);
        MVT::MvInit(*C_[i]);
      }
    }

    // Reset breakdown to false before initializing iteration
    breakdown_ = false;

    // NOTE:  In GCRIter R_, the initial residual, is required!!!
    //
    std::string errstr("Belos::GCRIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for one and zero.
    const ScalarType one = SCT::one();
    const int zero = Teuchos::ScalarTraits<int>::zero();

    //if (!Teuchos::is_null(newstate.R)) {
	//
    //  TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
    //                      std::invalid_argument, errstr );
    //  TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != numRHS_,
    //                      std::invalid_argument, errstr );
	//
    //  // Copy residual vectors from newstate into R
    //  if (newstate.R != R_) {
    //    // Assigned by the new state
    //    MVT::Assign(*newstate.R, *R_);
    //  }
    //  else {
    //    // Computed
    //    lp_->computeCurrResVec(R_.get());
    //  }
	//
    //}
    //else {
    //  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(newstate.R),std::invalid_argument,
    //                     "Belos::GCRIter::initialize(): GCRStateIterState does not have initial residual.");
    //}

    // Set pone to 1
    pone_.assign(numRHS_,one);

    // Set curIndex_ to 0
    curIndex_.assign(numRHS_,zero);

    // Set newIndex_ to 0
    newIndex_.assign(numRHS_,zero);

    // Pr(:,:) = Prec*RHS(:,:)
    Teuchos::RCP<const MV> Pr_ = lp_->getInitPrecResVec();

    // R = Pr - R;
    MVT::Assign(*Pr_, *R_); // *R_ was already initialized with zeros

    // The solver is initialized
    initialized_ = true;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void GCRIter<ScalarType,MV,OP>::iterate()
  {
    using Teuchos::RCP;

    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Create convenience variable for one.
    const ScalarType one = SCT::one();
    //
    std::vector<MagnitudeType> nrmval(1);
    std::vector<ScalarType> gcrAlpha(1);
    std::vector<ScalarType> gcrBeta(1);

    // Get the current solution std::vector.
    Teuchos::RCP<MV> X = lp_->getCurrLHSVec();

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed && !breakdown_) {

      // Increment the iteration
      iter_++;

      curIndex_ = newIndex_;

      // U[0:numRHS_-1](:,curIndex_(0:numRHS_-1))]=R(:,0:numRHS_-1);
      deep_copy(*R_, curIndex_, U_);

      // AxR_ = A*R_;
      lp_->applyOp(*R_,*AxR_);

      for(int i=0; i<numRHS_; i++) {
        std::vector<int> index1(1);
        std::vector<int> index2(1);
        std::vector<int> index3(1);
        index1[0] = i;
        index2[0] = curIndex_[i]; 
        Teuchos::RCP<const MV> AxR_sub = MVT::CloneView(*AxR_, index1);
        Teuchos::RCP<MV> C_sub = MVT::CloneViewNonConst(*C_[i], index2);
        Teuchos::RCP<MV> U_sub = MVT::CloneViewNonConst(*U_[i], index2);

        // C[i](:,curIndex_[i]) = Prec*AxR_(:,i)
        if(lp_->isLeftPrec()) {
          lp_->applyLeftPrec(*AxR_sub, *C_sub);
        }
        /*else if(lp_->isRightPrec()) {
          lp_->applyRightPrec(*AxR_sub, *C_sub);
        }*/
        else { // no preconditioner
          MVT::Assign(*AxR_sub, *C_sub);
        }

        MagnitudeType betaMin = 1.0e10;

        for (int j=0; j<std::min(numKrylovVecs_,iter_); j++) {
          if (j != curIndex_[i]) {
            //gcrBeta = dotHerm(C[i](:,j),C[i](:,curIndex_[i])
            index3[0] = j;
            Teuchos::RCP<MV> C_prev_sub = MVT::CloneViewNonConst(*C_[i], index3);
            Teuchos::RCP<MV> U_prev_sub = MVT::CloneViewNonConst(*U_[i], index3);

            MVT::MvDot(*C_sub, *C_prev_sub, gcrBeta);
            if (SCT::magnitude(gcrBeta[0]) < betaMin) {
              betaMin = SCT::magnitude(gcrBeta[0]);
              newIndex_[i] = j;
            }

            //C[i](:,curIndex_[i])=C[i](:,curIndex_[i])-gcrBeta*C[i](:,j);
            //U[i](:,curIndex_[i])=U[i](:,curIndex_[i])-gcrBeta*U[i](:,j);
            axpy(one, *C_sub, gcrBeta, *C_prev_sub, *C_sub, true);
            axpy(one, *U_sub, gcrBeta, *U_prev_sub, *U_sub, true);
          }
        }

        // update solution and residual
        // gcrAlpha = one/norm2(C[i](:,curIndex_[i]));
        MVT::MvNorm(*C_sub, nrmval);
        gcrAlpha[0] = one/nrmval[0];

        // U[i](:,curIndex_[i])=gcrAlpha*U[i](:,curIndex_[i]);
        // C[i](:,curIndex_[i])=gcrAlpha*C[i](:,curIndex_[i]);
        MVT::MvScale(*U_sub, gcrAlpha[0]);
        MVT::MvScale(*C_sub, gcrAlpha[0]);

        // gcrAlpha = dotHerm(C[i](:,curIndex_[i]), R(:,i));
        Teuchos::RCP<MV> R_sub = MVT::CloneViewNonConst(*R_, index1);
        MVT::MvDot(*R_sub, *C_sub, gcrAlpha);

        // X(:,i) = X(:,i) + gcrAlpha*U[i](:,curIndex_[i]);
        Teuchos::RCP<MV> X_sub = MVT::CloneViewNonConst(*X, index1);
        axpy(one, *X_sub, gcrAlpha, *U_sub, *X_sub);

        // R(:,i) = R(:,i) - gcrAlpha*C[i](:,curIndex_[i]);
        axpy(one, *R_sub, gcrAlpha, *C_sub, *R_sub, true);
      } // end for(int i=0; i<numRHS_; i++)

    } // end while (sTest_->checkStatus(this) != Passed)
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void GCRIter<ScalarType,MV,OP>::axpy(const ScalarType alpha, const MV & A,
                                       const std::vector<ScalarType> beta, const MV& B, MV& mv, bool minus)
  {
    Teuchos::RCP<const MV> A1, B1;
    Teuchos::RCP<MV> mv1;
    std::vector<int> index(1);

    for(int i=0; i<numRHS_; i++) {
      index[0] = i;
      A1 = MVT::CloneView(A,index);
      B1 = MVT::CloneView(B,index);
      mv1 = MVT::CloneViewNonConst(mv,index);
      if(minus) {
        MVT::MvAddMv(alpha,*A1,-beta[i],*B1,*mv1);
      }
      else {
        MVT::MvAddMv(alpha,*A1,beta[i],*B1,*mv1);
      }
    }
  }

  template <class ScalarType, class MV, class OP>
  void GCRIter<ScalarType,MV,OP>::deep_copy(const MV& Src,
                                            const std::vector<int>& colIndex,
                                            const std::vector<Teuchos::RCP<MV> >& Dst)
  {
    Teuchos::RCP<const MV> subSrc;
    std::vector<int> index1(1);
    std::vector<int> index2(1);

    for(int i=0; i<numRHS_; i++) {
      //Dst[i](:,colIndex[i]) = Src(:,i);
      index1[0] = i;  
      index2[0] = colIndex[i];
      subSrc   = MVT::CloneView(Src, index1);
      MVT::SetBlock(*subSrc, index2, *Dst[i]);
    }
  }

} // end Belos namespace

#endif /* BELOS_GCR_ITER_HPP */
