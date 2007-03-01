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

#ifndef LOCA_HOMOTOPY_ABSTRACTGROUP_H
#define LOCA_HOMOTOPY_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA { 

  namespace Homotopy {

    /*!  
     * \brief Interface to underlying groups for homotopy calculations
     */
    /*!
     * This abstract class provides an interface for a homotopy technique
     * for solving nonlinear equations.  See LOCA::Homotopy::Group for a 
     * description of the technique used.  This class provides a single
     * pure virtual method, augmentJacobianForHomotopy(), which scales the
     * diagonal of the Jacobian by a constant times the identity matrix.
     */
    class AbstractGroup : 
      public virtual LOCA::MultiContinuation::AbstractGroup {

    public:
  
      //! Default constructor.
      AbstractGroup() {}
 
      //! Destructor
      virtual ~AbstractGroup() {}

      /*! 
       * @name Pure virtual methods 
       * These methods must be defined by any concrete implementation
       */
      //@{

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$aJ+bI\f$ where \f$I\f$ is 
       * the identity matrix
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double a, double b) = 0;

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace Homotopy
} // namespace LOCA

#endif
