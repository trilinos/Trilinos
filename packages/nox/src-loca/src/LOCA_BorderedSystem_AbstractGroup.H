// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

#include "Teuchos_RefCountPtr.hpp"
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
      virtual Teuchos::RefCountPtr<const NOX::Abstract::Group>
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
