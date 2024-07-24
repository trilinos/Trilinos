// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_SmootherPreconditionerFactory_hpp__
#define __Teko_SmootherPreconditionerFactory_hpp__

// Teko includes
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_ImplicitLinearOp.hpp"
#include "Teko_RequestHandlerContainer.hpp"

namespace Teko {

/** This class applies an operator as a smoother.
 * The main idea being that a residual is used and
 * corrected to get a Gauss-Seidel like method.
 */
class SmootherLinearOp : public ImplicitLinearOp, public RequestHandlerContainer {
 public:
  SmootherLinearOp(const LinearOp &A, const LinearOp &invM, unsigned int applications,
                   bool useDestAsInitialGuess = false);
  SmootherLinearOp(const LinearOp &A, const LinearOp &invM, unsigned int applications,
                   unsigned int block);

  /** @brief Range space of this operator */
  virtual VectorSpace range() const { return invM_->range(); }

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const { return invM_->domain(); }

  /** @brief Perform a matrix vector multiply with this implicitly
   * defined blocked operator.
   *
   * The <code>apply</code> function takes one vector as input
   * and applies a linear operator. The result
   * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
   * \f$ y = \alpha M x + \beta y \f$
   *
   * @param[in]     x
   * @param[in,out] y
   * @param[in]     alpha (default=1)
   * @param[in]     beta  (default=0)
   */
  virtual void implicitApply(const MultiVector &x, MultiVector &y, const double alpha = 1.0,
                             const double beta = 0.0) const;

  //! Set the request handler with pointers to the appropriate callbacks
  virtual void setRequestHandler(const Teuchos::RCP<RequestHandler> &rh);

  //! Get the request handler with pointers to the appropriate callbacks
  virtual Teuchos::RCP<RequestHandler> getRequestHandler() const;

 private:
  // enum describing initial guess type
  typedef enum {
    Unspecified,
    RequestInitialGuess,
    DestAsInitialGuess,
    NoInitialGuess
  } InitialGuessType;

  LinearOp A_;                         // forward operator
  LinearOp invM_;                      // preconditioner
  unsigned int applications_;          // how much smoothing is required
  InitialGuessType initialGuessType_;  // type of initial guess to use

  // for producing an initial guess
  Teuchos::RCP<RequestMesg> requestMesg_;

  // used by RequestHandlerContainer interface
  Teuchos::RCP<RequestHandler> requestHandler_;

  SmootherLinearOp();                          // hide me
  SmootherLinearOp(const SmootherLinearOp &);  // hide me
};

LinearOp buildSmootherLinearOp(const LinearOp &A, const LinearOp &invM, unsigned int applications,
                               bool useDestAsInitialGuess = false);
LinearOp buildSmootherLinearOp(const LinearOp &A, const LinearOp &invM, unsigned int applications,
                               unsigned int initialGuessBlock);

/**
 */
class SmootherPreconditionerFactory : public virtual Teko::PreconditionerFactory {
 public:
  //! Default constructor, for use with the AutoClone class.
  SmootherPreconditionerFactory();

  /** \brief Function that is called to build the preconditioner
   *        for the linear operator that is passed in.
   *
   * This function builds a preconditioner based on the passed
   * in LinearOp.
   *
   * \param[in] lo    Source linear operator that is to be preconditioned.
   * \param[in] state An object associated with this operator to store
   *                  the preconditioner state.
   *
   * \returns The preconditioner as a linear operator (i.e. to perform
   *           a matrix-vector operation simply call "apply").
   */
  virtual LinearOp buildPreconditionerOperator(LinearOp &lo, PreconditionerState &state) const;

  //! @name Methods for construction from a parameter list entry
  //@{

  /** \brief This function builds the internals of the preconditioner factory
   *        from a parameter list.
   *
   * This function builds the internals of the preconditioner factory
   * from a parameter list. Furthermore, it allows a preconditioner factory
   * developer to easily add a factory to the build system. This function
   * is required for building a preconditioner from a parameter list.
   *
   * \param[in] settings Parmaeter list to use as the internal settings
   */
  virtual void initializeFromParameterList(const Teuchos::ParameterList &settings);

  //@}

 private:
  // enum describing initial guess type
  typedef enum {
    Unspecified,
    RequestInitialGuess,
    DestAsInitialGuess,
    NoInitialGuess
  } InitialGuessType;

  // parameters specifying behavior of smoother operator
  unsigned int sweepCount_;
  InitialGuessType initialGuessType_;
  unsigned int initialGuessBlock_;

  // prectionditioner to use as residual correction in smoother
  Teuchos::RCP<Teko::InverseFactory> precFactory_;
};

}  // end namespace Teko

#endif
