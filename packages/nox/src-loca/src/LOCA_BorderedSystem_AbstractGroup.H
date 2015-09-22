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

#ifndef LOCA_BORDEREDSYSTEM_ABSTRACTGROUP_H
#define LOCA_BORDEREDSYSTEM_ABSTRACTGROUP_H

#include "Teuchos_RCP.hpp"
#include "NOX_Abstract_Group.H"       // base class
#include "NOX_Abstract_MultiVector.H" // for DenseMatrix

namespace LOCA {

  namespace BorderedSystem {

    /*!
     * \brief An interface for groups that are bordered systems
     */
    /*!
     * This class provides an interface for groups whose Jacobian is of
     * the form
     * \f[
     *     \begin{bmatrix}
     *          J & A    \\
     *        B^T & C
     *     \end{bmatrix}
     * \f]
     * where \f$A\f$ and \f$B\f$ are multivectors and \f$C\f$ is a dense
     * matrix.  It provides methods for determining the width of the bordered
     * rows/columns and for extracting these components.  It is intended
     * to be a recusive interface, in that if the group representing
     * \f$J\f$ also has this form, the extracted rows/columns should be
     * the combined rows/columns of this group and the underlying group
     * (and so on).
     */
    class AbstractGroup : public virtual NOX::Abstract::Group {

    public:

      //! Constructor
      AbstractGroup() {}

      //! Destructor.
      virtual ~AbstractGroup() {}

      /*!
       * @name Pure virtual methods
       * These methods must be defined by any concrete implementation
       */
      //@{

      //! Return the total width of the bordered rows/columns
      virtual int getBorderedWidth() const = 0;

      //! Get bottom-level unbordered group
      virtual Teuchos::RCP<const NOX::Abstract::Group>
      getUnborderedGroup() const = 0;

      //! Indicates whether combined A block is zero
      virtual bool isCombinedAZero() const = 0;

      //! Indicates whether combined B block is zero
      virtual bool isCombinedBZero() const = 0;

      //! Indicates whether combined C block is zero
      virtual bool isCombinedCZero() const = 0;

      /*!
       * Given the vector \c v, extract the underlying solution component
       * corresponding to the unbordered group.
       */
      virtual void
      extractSolutionComponent(const NOX::Abstract::MultiVector& v,
                   NOX::Abstract::MultiVector& v_x) const = 0;

      /*!
       * Given the vector \c v, extract the parameter components of all
       * of the nested subvectors in \c v down to the solution component
       * for the unbordered group.
       */
      virtual void
      extractParameterComponent(
               bool use_transpose,
               const NOX::Abstract::MultiVector& v,
                       NOX::Abstract::MultiVector::DenseMatrix& v_p) const = 0;

      /*!
       * Given the solution component \c v_x and combined parameter
       * components \c v_p, distribute these components through the nested
       * sub-vectors in \c v.
       */
      virtual void
      loadNestedComponents(const NOX::Abstract::MultiVector& v_x,
               const NOX::Abstract::MultiVector::DenseMatrix& v_p,
               NOX::Abstract::MultiVector& v) const = 0;

      //! Fill the combined A block as described above
      virtual void fillA(NOX::Abstract::MultiVector& A) const = 0;

      //! Fill the combined B block as described above
      virtual void fillB(NOX::Abstract::MultiVector& B) const = 0;

      //! Fill the combined C block as described above
      virtual void fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const = 0;

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    }; // Class AbstractGroup

  } // Namespace BorderedSystem

} // Namespace LOCA

#endif // LOCA_BORDEREDSYSTEM_ABSTRACTGROUP_H
