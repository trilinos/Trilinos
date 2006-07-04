// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H
#define LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H

#include "NOX_Abstract_Group.H"           // base class

namespace LOCA {

  namespace Abstract {
    
    /*! 
     * \brief Abstract group interface class for solving the transpose
     * of the Jacobian.
     */
    /*!
     * This interface, derived from NOX::Abstract::Group, provides the
     * additional interface for solving the transpose of the Jacobian.
     */
    class TransposeSolveGroup : public virtual NOX::Abstract::Group {

    public:

      //! Constructor
      TransposeSolveGroup() {}

      //! Destructor
      virtual ~TransposeSolveGroup() {}

      //! Solve Jacobian-tranpose system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverse(Teuchos::ParameterList& params, 
				    const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const = 0;

      //! Solve Jacobian-tranpose system with multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverseMultiVector(
				Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input, 
				NOX::Abstract::MultiVector& result) const = 0;

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    }; // Class TransposeSolveGroup

  } // Namespace Abstract

} // Namespace LOCA

#endif // LOCA_ABSTRACT_TRANSPOSESOLVEGROUP_H
