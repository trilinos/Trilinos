// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BIFURCATION_PITCHFORKBORD_STATUSTEST_NULLVECTORNORMWRMS_H
#define LOCA_BIFURCATION_PITCHFORKBORD_STATUSTEST_NULLVECTORNORMWRMS_H

#include "NOX_StatusTest_Generic.H"    // base class

namespace LOCA {

  namespace Bifurcation {

    namespace PitchforkBord {

      namespace StatusTest {

    /*!
     * \brief A convergence test based on the weighted root-mean-square
     * norm of the update to the null vector component for pitchfork
     * location.
     */
    /*!
     * Let \f$n\f$ be the approximation to the null vector for pitchfork
     * tracking (see LOCA::Bifurcation::PitchforkBord::ExtendedGroup).
     * This convergence test defines convergence for the null vector when
     * the following is true
     * \f[
     \sqrt{\frac{1}{N}\sum_{i=1}^N\left(\frac{|n_i-(n_0)_i|}{\epsilon_r|n_i| + \epsilon_a}\right)} < \tau
     * \f]
     * where \f$n_0\f$ is the previous approximation to the null vector,
     * \f$N\f$ is the length of \f$n\f$, \f$\epsilon_r\f$
     * is the relative tolerance, \f$\epsilon_a\f$ is the absolute
     * tolerance, and \f$\tau\f$ is an overall scale factor (typically
     * \f$\tau = 1\f$).
     *
     * Note that this status test deals only with the null vector component
     * of the pitchfork equations.  This status test should be combined
     * with other status tests for the solution and parameter components
     * (using NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to
     * build a composite status test for the entire system.
     *
     * Also note that if the group returned by the getSolutionGroup()
     * method of the solver supplied in checkStatus() is not a
     * pitchfork group (i.e., not derived from
     * LOCA::Bifurcation::PitchforkBord::ExtendedGroup), checkStatus()
     * returns
     * NOX::StatusTest::Converged.  This allows the status test to
     * be used in situations other than pitchfork tracking, e.g.,
     * steady-state solves, without raising error conditions.
     */
    class NullVectorNormWRMS : public NOX::StatusTest::Generic {

    public:

      //! Constructor
      /*!
       * \em rtol is the relative tolerance \f$\epsilon_r\f$, \em atol
       * is the absolute tolerance \f$\epsilon_a\f$, and \em tol is the
       * overall scale factor \f$\tau\f$ defined above.
       */
      NullVectorNormWRMS(double rtol, double atol, double tol);

      //! Destructor.
      virtual ~NullVectorNormWRMS();

      //! Evaluates convergence criteria specified above
      virtual NOX::StatusTest::StatusType
      checkStatus(const NOX::Solver::Generic& problem);

      //! Returns status as defined above
      virtual NOX::StatusTest::StatusType
      getStatus() const;

      //! Prints current status
      virtual std::ostream& print(ostream& stream, int indent = 0) const;


      /* @name Accessor Functions
       * Used to query current values of variables in the status test.
       */
      //@{

      //! Returns the value of weighted parameter update norm.
      double getNullVectorNormWRMS() const;

      //! Returns the realative tolerance set in the constructor.
      double getRTOL() const;

      //! Returns the absolute tolerance set in the constructor.
      double getATOL() const;

      //! Returns the tolerance set in the constructor
      double getTOL() const;

      //@}

    protected:

      //! Relative tolerance for convergence.
      double rtol;

      //! Absolute tolerance for convergence.
      double atol;

      //! Tolerance for convergence
      double tol;

      //! Current value of the parameter update norm
      double normWRMS;

      //! %Status of entire %status test
      NOX::StatusTest::StatusType status;
    };
      } // namespace StatusTest
    } // namespace PitchforkBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
