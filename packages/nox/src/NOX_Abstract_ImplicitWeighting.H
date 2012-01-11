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

#ifndef NOX_ABSTRACT_IMPLICIT_WEIGHTING_H
#define NOX_ABSTRACT_IMPLICIT_WEIGHTING_H

#include "NOX_Common.H"  // for NOX_Config.h

namespace NOX {
namespace Abstract {

/*! 
  \brief A pure virtual interface for enabling/disabling any user defined implicit weighting of the concrete object.

  This class is used to toggle on/off any implicit weighting that
  might be handled internal to a particular user defined class.  The
  classical use case is to disable an implicitly defined weighting on
  the inner product and norm computations of the abstract vector class
  so that classes like the NormWRMS can apply their own weighting.
  This prevents a "double weight" effect being applied to specific
  status tests.

  The typical use case is the get the current weighting flag, then
  disable the weighting, perform the computations on the unweighted
  object and then set the flag back to its original value.  By storing
  the original value, one can return it to its prior state regardless
  if an implicit weighting is defined.

  \code
  NOX::ImplicitWeighting& v;

  bool saved_status = v.getImplicitWeighting();
  v.setImplicitWeighting(false);
  .
  <perform unweighted computaiton>
  .
  v.setImplicitWeighting(saved_status);
  \endcode

  CAUTION: This is a power user feature and should only be used in
  concert with specialized NOX::Abstract::Vector implementations.

  \author Roger Pawlowski (SNL 1444)
*/
class ImplicitWeighting {

public:
  
  virtual ~ImplicitWeighting() {}

  /** \brief Return true if implicit weighting is currently enabled

      CAUTION: This is a power user feature and should only be used in
      concert with specialized NOX::Abstract::Vector implementations.
  */
  virtual bool getImplicitWeighting() const = 0;

  /** \brief Set to true to enable implicit weighting, false disables

      CAUTION: This is a power user feature and should only be used in
      concert with specialized NOX::Abstract::Vector implementations.
  */
  virtual void setImplicitWeighting(bool do_implicit_weighting) = 0;

};

} // namespace Abstract
} // namespace NOX

#endif
