// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_BELOS_JACOBIANOPERATOR_H
#define NOX_BELOS_JACOBIANOPERATOR_H

#include "BelosOperator.hpp"            // base class
#include "NOX_Abstract_Group.H"         // data member

namespace NOX { 

  namespace Belos { 

    /*! 
     * \brief Implementation of Belos::Operator for the Jacobian 
     * matrix-vector product defined by a NOX::Abstract::Group
     */
    /*!
     * This class implements the Belos::Operator interface representing a 
     * Jacobian matrix using the NOX::Abstract::Group.  It implements the 
     * operator Apply() method via applyJacobianMultiVector or 
     * applyJacobianTransposeMultiVector.
     */
    class JacobianOperator : public ::Belos::Operator<double> {

      public:

      //! Constructor.
      JacobianOperator(NOX::Abstract::Group& g);

      //! Destructor.
      virtual ~JacobianOperator();

      //! Apply operator
      /*!
       * Calls underlying group's applyJacobianMultiVector or 
       * applyJacobianTransposeMultiVector depending on \em trans
       */
      virtual ::Belos::ReturnType 
      Apply(const ::Belos::MultiVec<double>& x, 
	    ::Belos::MultiVec<double>& y, 
	    ::Belos::ETrans trans=::Belos::NOTRANS) const;

      //! ApplyInverse operator
      /*!
       * Returns Belos::Undefined
       */
      virtual ::Belos::ReturnType 
      ApplyInverse(const ::Belos::MultiVec<double>& x, 
		   ::Belos::MultiVec<double>& y, 
		   ::Belos::ETrans trans=::Belos::NOTRANS) const;

    protected:

      /*! 
       * \brief Converts NOX::Abstract::Group::ReturnType to the corresponding
       * Belos type.
       */
      ::Belos::ReturnType 
      noxReturnTypeToBelos(NOX::Abstract::Group::ReturnType noxStatus) const;

    protected:

      //! Underlying group
      NOX::Abstract::Group& grp;  
  
    };
  } // namespace Belos
} // namespace LOCA

#endif
