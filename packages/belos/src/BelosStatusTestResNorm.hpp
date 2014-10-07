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
