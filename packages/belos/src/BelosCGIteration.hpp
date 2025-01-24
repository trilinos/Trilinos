// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_CG_ITERATION_HPP
#define BELOS_CG_ITERATION_HPP

/*! \file BelosCGIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a conjugate gradient linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"
#include "Teuchos_Assert.hpp"

namespace Belos {

enum CGIterType {
  CG,
  CGSingleReduce,
  PseudoBlockCG,
  BlockCG
};

  //! @name CGIteration Structures
  //@{

  /** \brief Structure to contain pointers to CGIteration state variables.
   *
   * This struct is utilized by CGIteration::initialize() and CGIteration::getState().
   */
  template <class ScalarType, class MV>
  class CGIterationState {

  public:
    CGIterationState()
    : S(Teuchos::null), R(Teuchos::null), Z(Teuchos::null), P(Teuchos::null),
          AP(Teuchos::null), isInitialized_(false) {}

    CGIterationState(Teuchos::RCP<const MV> tmp) {
      initialize(tmp);
    }

    bool isInitialized() const {
      return isInitialized_;
    }

    void initialize(Teuchos::RCP<const MV> tmp, CGIterType iterType, int numRHS_) {
      using MVT = MultiVecTraits<ScalarType, MV>;
      if (iterType == CG) {
        TEUCHOS_ASSERT(numRHS_ == 1);

        // S = (R, Z)
        // This allows to compute the inner products (R, S) = ((R, R), (R, Z)) using a single reduction.
        S = MVT::Clone( *tmp, 2 );
        std::vector<int> index(1,0);
        index[0] = 0;
        R = MVT::CloneViewNonConst( *S, index );
        index[0] = 1;
        Z = MVT::CloneViewNonConst( *S, index );

        P = MVT::Clone( *tmp, 1 );
        AP = MVT::Clone(*tmp, 1);
      } else if (iterType == CGSingleReduce) {
        TEUCHOS_ASSERT(numRHS_ == 1);

        // W = (AZ, R, Z)
        W = MVT::Clone( *tmp, 3 );
        std::vector<int> index2(2,0);
        std::vector<int> index(1,0);

        // S = (AZ, R)
        index2[0] = 0;
        index2[1] = 1;
        S = MVT::CloneViewNonConst( *W, index2 );

        // U = (AZ, Z)
        index2[0] = 0;
        index2[1] = 2;
        U = MVT::CloneViewNonConst( *W, index2 );

        index[0] = 1;
        R = MVT::CloneViewNonConst( *W, index );
        index[0] = 0;
        AZ = MVT::CloneViewNonConst( *W, index );
        index[0] = 2;
        Z = MVT::CloneViewNonConst( *W, index );

        // T = (R, Z)
        index2[0] = 1;
        index2[1] = 2;
        T = MVT::CloneViewNonConst( *W, index2 );

        // V = (AP, P)
        V = MVT::Clone( *tmp, 2 );
        index[0] = 0;
        AP = MVT::CloneViewNonConst( *V, index );
        index[0] = 1;
        P = MVT::CloneViewNonConst( *V, index );
      } else if ((iterType == PseudoBlockCG) || (iterType == BlockCG)) {
        R = MVT::Clone( *tmp, numRHS_ );
        Z = MVT::Clone( *tmp, numRHS_ );
        P = MVT::Clone( *tmp, numRHS_ );
        AP = MVT::Clone(*tmp, numRHS_ );
      }
      numRHS = numRHS_;
      isInitialized_ = true;
    }

    bool matches(Teuchos::RCP<const MV> tmp, int numRHS_=1) const {
      using MVT = MultiVecTraits<ScalarType, MV>;
      return (isInitialized() && !S.is_null() && !R.is_null() &&
              !Z.is_null() && !P.is_null() && !AP.is_null() && (numRHS == numRHS_) && (MVT::GetGlobalLength(*tmp) == MVT::GetGlobalLength(*R)));
    }

    Teuchos::RCP<MV> W;
    Teuchos::RCP<MV> V;
    Teuchos::RCP<MV> U;
    Teuchos::RCP<MV> S;
    Teuchos::RCP<MV> T;

    /*! \brief The current residual. */
    Teuchos::RCP<MV> R;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<MV> Z;

    /*! \brief The current decent direction vector */
    Teuchos::RCP<MV> P;

    /*! \brief The matrix A applied to current decent direction vector */
    Teuchos::RCP<MV> AP;

    Teuchos::RCP<MV> AZ;

  private:

    bool isInitialized_;
    int numRHS;

};

  //! @name CGIteration Exceptions
  //@{ 
  
  /** \brief CGIterationInitFailure is thrown when the CGIteration object is unable to
   * generate an initial iterate in the CGIteration::initialize() routine. 
   *
   * This std::exception is thrown from the CGIteration::initialize() method, which is
   * called by the user or from the CGIteration::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown, 
   * CGIteration::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class CGIterationInitFailure : public BelosError {public:
    CGIterationInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGIterateFailure is thrown when the CGIteration object is unable to
   * compute the next iterate in the CGIteration::iterate() routine. 
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterateFailure : public BelosError {public:
    CGIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGPositiveDefiniteFailure is thrown when the the CG 'alpha = p^H*A*P' value
   * is less than zero, indicating a breakdown in the computation due to roundoff or
   * a non-positive-definite matrix.
   *
   */
  class CGPositiveDefiniteFailure : public CGIterateFailure {public:
    CGPositiveDefiniteFailure(const std::string& what_arg) : CGIterateFailure(what_arg)
    {}};

  /** \brief CGIterationOrthoFailure is thrown when the CGIteration object is unable to
   * compute independent direction vectors in the CGIteration::iterate() routine. 
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterationOrthoFailure : public BelosError {public:
    CGIterationOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGIterationLAPACKFailure is thrown when a nonzero return value is passed back
   * from an LAPACK routine.
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterationLAPACKFailure : public BelosError {public:
    CGIterationLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


template<class ScalarType, class MV, class OP>
class CGIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %CGIteration contains a certain amount of state, consisting of the current 
   * residual, preconditioned residual, and decent direction.
   *
   * initialize() gives the user the opportunity to manually set these,
   * although only the current unpreconditioned residual is required.
   *
   * \post 
   * <li>isInitialized() == \c true (see post-conditions of isInitialize())
   *
   * \note For any pointer in \c newstate which directly points to the multivectors in 
   * the solver, the data is not copied.
   */
  virtual void initializeCG(CGIterationState<ScalarType,MV>& newstate, Teuchos::RCP<MV> R_0) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A CGIterationState object containing const pointers to the current solver state.
   */
  virtual CGIterationState<ScalarType,MV> getState() const = 0;
  //@}


  //! Sets whether or not to store the diagonal for condition estimation
  virtual void setDoCondEst(bool val) = 0;

  //! Gets the diagonal for condition estimation
  virtual Teuchos::ArrayView<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getDiag() = 0;

  //! Gets the off-diagonal for condition estimation 
  virtual Teuchos::ArrayView<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getOffDiag() = 0;

};

} // end Belos namespace

#endif /* BELOS_CG_ITERATION_HPP */
