// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_ANASAZIOPERATOR_CAYLEY_H
#define LOCA_ANASAZIOPERATOR_CAYLEY_H

#include "Teuchos_RCP.hpp"

#include "LOCA_AnasaziOperator_AbstractStrategy.H"  // base class
#include "LOCA_TimeDependent_AbstractGroup.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace AnasaziOperator {

    /*! \brief Anasazi operator for computing generalized eigenvalues
     * using Cayley transformations.
     */
    /*!
     * This class implements the LOCA::AnasaziOperator::AbstractStrategy
     * interface for computing generalized eigenvalues \f$\lambda\f$ and
     * eigenvectors \f$z\f$ of the system
     *\f[
     *    J z = \lambda M z
     *\f]
     * where \f$J\f$ is the Jacobian matrix and \f$M\f$ is the mass matrix.
     * The  eigenvalues are computed using a Cayley transformation, i.e.
     * solving
     * \f[
     *     (J - \sigma M) z = (J - \mu M) r
     * \f]
     * where \f$\sigma\f$ is the Cayley pole and \f$\mu\f$ is the Cayley
     * zero.
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "Cayley Pole" - \f$\sigma\f$ as defined above (Default 0.0)
     * <li> "Cayley Zero" - \f$\mu\f$ as defined above (Default 0.0)
     * </ul>
     *
     * Also the \em grp argument to the constructor must
     * be a child of LOCA::TimeDependent::AbstractGroup for the shift-invert
     * operations.
     */
    class Cayley : public LOCA::AnasaziOperator::AbstractStrategy {

    public:

      //! Constructor.
      /*!
       * Argument \em grp must be of type LOCA::TimeDependent::AbstractGroup.
       * See class description for a list of eigenParams.
       */
      Cayley(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
    const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
    const Teuchos::RCP<LOCA::TimeDependent::AbstractGroup>& grp);

      //! Destructor
      virtual ~Cayley();

      /*!
       * \brief Return name of this operator
       */
      virtual const std::string& label() const;

      /*!
       * \brief Apply the operator.
       *
       * Applies the inverse of the shifted operator, i.e., solves
       * \f[
       *      (J-\omega I)z = M r
       * \f]
       * for \f$z\f$, where \f$r = \mbox{input}\f$ and \f$z = \mbox{output}\f$.
       */
      virtual void
      apply(const NOX::Abstract::MultiVector& input,
        NOX::Abstract::MultiVector& output) const;

      /*!
       * \brief PreProcess the random seed vector
       *
       * Performs one backward Euler iteration on the random
       * initial seed vector, to satisfy contraints
       */
      virtual void
      preProcessSeedVector(NOX::Abstract::MultiVector& ivec);

      /*!
       * \brief Transform eigenvalue.
       *
       * Transforms the given eigenvalue to the eigenvalue of the
       * Jacobian-mass matrix system by shifting and inverting it.
       */
      virtual void
      transformEigenvalue(double& ev_r, double& ev_i) const;

      /*!
       * \brief Compute Rayleigh quotient
       *
       * Computes the Rayleigh quotient \f$z^T J z / z^T M z\f$ for
       * the eigenvector \f$z\f$.
       */
      virtual NOX::Abstract::Group::ReturnType
      rayleighQuotient(NOX::Abstract::Vector& evec_r,
               NOX::Abstract::Vector& evec_i,
               double& rq_r, double& rq_i) const;

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Name of this operator
      std::string myLabel;

      //! Stores parameters relating to the operator
      mutable Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! Stores linear solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Stores group representing Jacobian and Mass matrix
      Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> grp;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::MultiVector> tmp_r;

      //! Stores a temporary vector for computing Rayleigh quotients
      mutable Teuchos::RCP<NOX::Abstract::MultiVector> tmp_i;

      //! Stores Cayley pole \f$\sigma\f$
      double sigma;

      //! Stores Cayley zero \f$\mu\f$
      double mu;

    };
  } // namespace AnasaziOperator
} // namespace LOCA

#endif
