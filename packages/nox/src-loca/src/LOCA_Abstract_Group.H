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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "NOX_Abstract_Group.H"    // base class

#include "LOCA_Abstract_Vector.H"  // for NOX::CopyType
#include "NOX_Common.H"            // for string 
#include "LOCA_DerivUtils.H"

namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace LOCA {
  class ParameterVector;
}

namespace LOCA { 
namespace Abstract { 

/*!  \brief %NOX pure abstract interface to a "group"; i.e., a
  solution vector and the corresponding F-vector, Jacobian matrix,
  gradient vector, and Newton vector.
*/
/*! 
  This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.

  Recall that our goal is to solve
  \f[
  F(x) = \left[ \begin{array}{c} F_1(x) \\ F_2(x) \\ \vdots \\ F_n(x) \\ \end{array} \right] = 0.
  \f]
  We may also consider the equivalent optimization problem
  \f[
  \min f(x) \equiv \frac{1}{2} F(x)^T F(x).
  \f]

  The \b Jacobian
  is denoted by \f$J\f$, and defined by
  \f[
  J_{ij} = \frac{\partial F_i}{\partial x_j} (x).
  \f]
  The \b Newton \b direction is the solution of
  \f[
  Js = -F.
  \f]
  The \b gradient (of \f$f\f$) is defined as \f$g \equiv J^T F\f$.
  
*/

class Group : public virtual NOX::Abstract::Group {

public:
  
  //! Constructor.
  Group(const DerivUtils& deriv = DerivUtils());
 
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();
  
  /*! 
    \brief Copies the values of all vectors and any other data in
    the source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual Group& operator=(const Group& source);
     
  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is \c DeepCopy, then we need to create an exact replica of
    "this". Otherwise, if type is \c ShapeCopy, we need only replicate
    the shape of "this". Returns \c NULL if clone is not supported.
  */
//  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Set the parameter list in the group to p (pVector = p).  
  virtual bool setParams(const ParameterVector& p) = 0;

  //! Compute a new parameter vector where this.pVector = grp.pVector + step*d.
  virtual bool computeParams(const ParameterVector& oldParams, 
			     const ParameterVector& direction, 
			     double step) = 0;

  //! Return a const reference to the ParameterVector owned by the group. 
  virtual const ParameterVector& getParams() const = 0;

  //! Compute tangent direction
  virtual bool computeTangent(NOX::Parameter::List& params, int paramID) = 0;

  //! Get tangent direction vector, needed for predictor
  virtual const NOX::Abstract::Vector& getTangent() const = 0;

  virtual bool computeDfDp(int paramID, NOX::Abstract::Vector& result);

  virtual bool print() const = 0;
protected:
  DerivUtils* derivPtr;
};
} // namespace Abstract
} // namespace LOCA

#endif
