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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_ABSTRACT_VECTOR_H
#define LOCA_ABSTRACT_VECTOR_H

#include "NOX_Abstract_Vector.H"  // Base class

//! Nonlinear solvers package namespace
namespace LOCA {

//! %NOX abstract interface for vector and group
namespace Abstract {

/*! 
  \brief %LOCA's pure abstract vector interface for vectors that are
  used by the NLS solver.  Derived from the NOX::Abstract::Vector class.

  The documentation throughout assumes that \f$x\f$ denotes "this"
  vector, \f$x_i\f$ is the i-th component (or entry) of \f$x\f$, and
  \f$n\f$ is the length of \f$x\f$.

  \author Roger Pawlowski (SNL 9233)
*/

class Vector : public NOX::Abstract::Vector {

public:
  
  //! %Vector constructor.
  Vector() {};

  //! %Vector destructor.
  virtual ~Vector() {};

  //! Copies source vector y into "this".
  /*! \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& operator=(const NOX::Abstract::Vector& y) = 0;

  //! Copies source vector y into "this".
  /*! \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& operator=(const Vector& y) = 0;

  /*! \brief Create a new %Vector of the same underlying type by
    cloning "this", and return a pointer to the new vector.  */
  /*! If type is \c DeepCopy, then we need to create an exact replica
    of "this". Otherwise, if type is \c ShapeCopy, we need only
    replicate the shape of "this". Returns \c NULL if clone is not
    supported. */
  virtual Vector* clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

  //! Fill the vector with random scalars between zero and one.
  virtual Vector& random() = 0;

}; // class Vector
} // namespace Abstract
} // namespace LOCA

#endif
