// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PITCHFORK_MOORESPENCE_ABSTRACTGROUP_H
#define LOCA_PITCHFORK_MOORESPENCE_ABSTRACTGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"    // base class

namespace LOCA {

  //! Groups and vectors for pitchfork bifurcations
  namespace Pitchfork {

    /*!
     * \brief Groups and vectors for locating pitchfork bifurcations
     * using the Moore-Spence formulation.
     */
    namespace MooreSpence {

      /*!
       * \brief Interface to underlying groups for pitchfork calculations
       * using the Moore-Spence formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate pitchforks using the bordering algorithms for the
       * Moore-Spence pitchfork formulation (see
       * LOCA::Pitchfork::MooreSpence::ExtendedGroup for a description of the
       * governing equations).
       *
       * This class is derived from the
       * LOCA::TurningPoint::MooreSpence::AbstractGroup and declares a single
       * virtual method, innerProduct(), to compute the inner product of a
       * vector with the asymmetry vector.  It has a default implementation
       * given by the dot product, but should be overloaded for any
       * problem that has a different definition for the inner product.
       */
      class AbstractGroup :
    public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup {

      public:

    //! Default constructor.
    AbstractGroup() {}

    //! Destructor
    virtual ~AbstractGroup() {}

    //! Compute the inner product of \c a and \c b.
    /*!
     * The default implementation is given by the dot product of
     * \c a and \c b.
     */
    virtual double innerProduct(const NOX::Abstract::Vector& a,
                    const NOX::Abstract::Vector& b) const {
      return a.innerProduct(b);
    }

    //! Compute the inner product of \c a and \c b.
    /*!
     * The default implementation is given by the dot product of
     * \c a and \c b.
     */
    virtual void innerProduct(
                const NOX::Abstract::MultiVector& a,
                const NOX::Abstract::MultiVector& b,
                NOX::Abstract::MultiVector::DenseMatrix& c) const {
      b.multiply(1.0, a, c);
    }

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MooreSpence

  } // namespace Pitchfork

} // namespace LOCA

#endif
