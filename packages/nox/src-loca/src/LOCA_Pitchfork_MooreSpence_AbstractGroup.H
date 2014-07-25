// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
