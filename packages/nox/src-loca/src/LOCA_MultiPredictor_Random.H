// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTIPREDICTOR_RANDOM_H
#define LOCA_MULTIPREDICTOR_RANDOM_H

#include "LOCA_MultiPredictor_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA{
  class GlobalData;
}

namespace LOCA {

  namespace MultiPredictor {

    //! %Random predictor strategy
    /*!
     * This class computes the predictor direction where the solution
     * component is filled with random values and the parameter component
     * equal to 1.  Each componenet of the solution vector \f$v_i\f$ of the
     * predictor is given by \f$v_i = \epsilon r_i x_i\f$ where
     * \f$r_i\f$ is a random value between -1 and 1, \f$x_i\f$ is the
     * corresponding component of the solution vector, and \f$\epsilon\f$
     * is a parameter.
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "Epsilon" - \f$\epsilon\f$ as defined above (Default 1.0e-3)
     * </ul>
     */
    class Random : public LOCA::MultiPredictor::AbstractStrategy {

    public:

      //! Constructor.
      /*!
       * \param global_data [in] Global data object
       * \param predParams [in] Predictor parameters as described above.
       */
      Random(const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const Teuchos::RCP<Teuchos::ParameterList>& predParams);

      //! Destructor
      virtual ~Random();

      //! Copy constructor
      Random(const Random& source, NOX::CopyType type = NOX::DeepCopy);

      //! Assignment operator
      virtual LOCA::MultiPredictor::AbstractStrategy&
      operator=(const LOCA::MultiPredictor::AbstractStrategy& source);

      //! Clone function
      virtual Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      /*!
       * \brief Compute the predictor given the current and
       * previous solution vectors.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      /*!
       * This method actually implements the predictor computation described
       * above
       */
      virtual NOX::Abstract::Group::ReturnType
      compute(bool baseOnSecant, const std::vector<double>& stepSize,
          LOCA::MultiContinuation::ExtendedGroup& grp,
          const LOCA::MultiContinuation::ExtendedVector& prevXVec,
          const LOCA::MultiContinuation::ExtendedVector& xVec);

      //! Evaluate predictor with step size \c stepSize.
      /*!
       * This method computes
       * \c result[i] = \c xVec[i] + \c stepSize[i] * v[i]
       * for each i, where \c v[i] is the ith predictor direction.
       */
      virtual NOX::Abstract::Group::ReturnType
      evaluate(const std::vector<double>& stepSize,
           const LOCA::MultiContinuation::ExtendedVector& xVec,
           LOCA::MultiContinuation::ExtendedMultiVector& result) const;

      //! Compute tangent to predictor and store in \c tangent.
      virtual NOX::Abstract::Group::ReturnType
      computeTangent(LOCA::MultiContinuation::ExtendedMultiVector& tangent);

      //! Is the tangent vector for this predictor scalable
      /*!
       * For the random predictor, this always returns false.
       */
      virtual bool isTangentScalable() const;

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Stores predictor vector
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> predictor;

      //! Stores secant vector for setting orientation
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> secant;

      //! Flag indicating whether vectors have been initialized
      bool initialized;

      //! Size of random value
      double epsilon;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
