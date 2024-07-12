// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_ABSTRACTSTRATEGY_H
#define LOCA_MULTICONTINUATION_ABSTRACTSTRATEGY_H

#include <vector>

#include "Teuchos_RCP.hpp"

#include "LOCA_Extended_MultiAbstractGroup.H" // base class
#include "LOCA_Abstract_Iterator.H" // for StepStatus

// forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedVector;
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    //! Abstract interface class for continuation strategies
    /*!
     * AbstractStrategy defines an abstract interface for continuation
     * strategies.  This interface is used by the LOCA::Stepper to manipulate
     * continuation groups in a consistent manner.  It defines a number
     * pure virtual methods that all continuation groups must implement.
     */
    class AbstractStrategy :
      public virtual LOCA::Extended::MultiAbstractGroup {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source) = 0;

      //! Returns number of parameters
      virtual int getNumParams() const = 0;

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.
       */
      virtual void
      preProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus) = 0;

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.
       */
      virtual void
      postProcessContinuationStep(
             LOCA::Abstract::Iterator::StepStatus stepStatus) = 0;

      //! Compute predictor directions
      virtual NOX::Abstract::Group::ReturnType
      computePredictor() = 0;

      //! Is Predictor valid
      virtual bool isPredictor() const = 0;

      //! Scales tangent to predictor
      virtual void scaleTangent() = 0;

      //! Sets tangent to predictor
      /*!
       * This is required by MF which takes the tangent space,
       * orthogonalizes it, and then sets it back in the group.
       */
      virtual void
      setPredictorTangentDirection(
                const LOCA::MultiContinuation::ExtendedVector& v,
                int i) = 0;

      //! Returns tangent to predictor
      virtual const LOCA::MultiContinuation::ExtendedMultiVector&
      getPredictorTangent() const = 0;

      //! Returns scaled tangent to predictor
      virtual const LOCA::MultiContinuation::ExtendedMultiVector&
      getScaledPredictorTangent() const = 0;

      //! Set the previous solution vector y
      virtual void setPrevX(const NOX::Abstract::Vector& y) = 0;

      //! Gets the previous solution vector
      virtual const LOCA::MultiContinuation::ExtendedVector&
      getPrevX() const = 0;

      //! Set step size for continuation constraint equation \em i
      virtual void setStepSize(double deltaS, int i = 0) = 0;

      //! Get step size for continuation constraint equation \em i
      virtual double getStepSize(int i = 0) const = 0;

      //! Sets the value for continuation parameter \em i
      virtual void setContinuationParameter(double val, int i = 0) = 0;

      //! Returns the value for continuation parameter \em i
      virtual double getContinuationParameter(int i = 0) const = 0;

      //! Get the continuation parameter id for parameter \em i
      virtual int getContinuationParameterID(int i = 0) const = 0;

      //! Get the continuation parameter ids
      virtual const std::vector<int>& getContinuationParameterIDs() const = 0;

      //! Get the continuation parameter id for parameter \em i
      virtual std::string getContinuationParameterName(int i = 0) const = 0;

      //! Returns step size scale factor for constraint equation \em i
      virtual double getStepSizeScaleFactor(int i = 0) const = 0;

      //! Prints the group
      virtual void printSolution() const = 0;

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
                   const NOX::Abstract::Vector& x,
                   const NOX::Abstract::Vector& y) const = 0;

      //! Returns dimension of project to draw array
      virtual int projectToDrawDimension() const = 0;

      //! Fills the project to draw array
      virtual void projectToDraw(
                const LOCA::MultiContinuation::ExtendedVector& x,
                double *px) const = 0;

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    }; // Class AbstractStrategy

  } // Namespace MultiContinuation

} // Namespace LOCA

#endif // LOCA_MULTICONTINUATION_ABSTRACTSTRATEGY_H
