// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_STATUS_TEST_RESNORM_H
#define BELOS_STATUS_TEST_RESNORM_H

/*!
  \file BelosStatusTestResNorm.hpp
  \brief Belos::StatusTest abstract class for specifying a residual norm stopping criteria.
*/

#include "BelosStatusTest.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"

/*! 
  \class Belos::StatusTestResNorm
  \brief An abstract class of StatusTest for stopping criteria using residual norms.
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestResNorm: public StatusTest<ScalarType,MV,OP> {

 public:

  // Convenience typedefs
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef MultiVecTraits<ScalarType,MV>  MVT;

  //! @name Form and parameter definition methods.
  //@{ 

  //! Set the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual, 
    we find that the initial tolerance was too tight or too lax.
  */
  virtual int setTolerance(MagnitudeType tolerance) = 0;

  //! Sets the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  virtual int setQuorum(int quorum) = 0;

  //! Set whether the only maximum residual norm is displayed when the print() method is called
  virtual int setShowMaxResNormOnly(bool showMaxResNormOnly) = 0;

  //! Define the form of the scaling for the residual.
  virtual int defineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, MagnitudeType ScaleValue = Teuchos::ScalarTraits<MagnitudeType>::one()) = 0;

  //@}

  //! @name Methods to access data members.
  //@{ 

  //! Returns the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  virtual int getQuorum() const = 0;

  //! Returns whether the only maximum residual norm is displayed when the print() method is called
  virtual bool getShowMaxResNormOnly() = 0;

  //! Returns the std::vector containing the indices of the residuals that passed the test.
  virtual std::vector<int> convIndices() = 0;

  //! Returns the value of the tolerance, \f$ \tau \f$, set in the constructor.
  virtual MagnitudeType getTolerance() const = 0;
 
  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  virtual const std::vector<MagnitudeType>* getTestValue() const = 0;
 
  //! Returns the current solution estimate that was computed for the most recent residual test.
  //! \note This method will return a null pointer if no vector was computed.
  virtual Teuchos::RCP<MV> getSolution() = 0;

  //! Returns a boolean indicating a loss of accuracy has been detected in computing the residual.
  virtual bool getLOADetected() const = 0;

  //@}

};

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_RESNORM_H */
