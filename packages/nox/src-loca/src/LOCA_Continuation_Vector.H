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

#ifndef LOCA_CONTINUATION_VECTOR_H
#define LOCA_CONTINUATION_VECTOR_H

#include "LOCA_MultiVector.H"  // Base class

//! Nonlinear solvers package namespace
namespace LOCA {

//! %NOX  bifurcation vector and group
namespace Continuation {

class Vector : public LOCA::MultiVector {

public:
  
  //! %Vector constructor.
  Vector(const NOX::Abstract::Vector& xVec, double param);

  //! %Copy constructor.
  Vector(const Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //! %Vector destructor.
  virtual ~Vector();

  virtual Vector& operator=(const Vector& y);

  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Sets the %ArcLengthVector by setting its two components 
  virtual void setVec(const NOX::Abstract::Vector& xVec, double param);

  //! Returns the solution vector component of extended vector
  virtual const NOX::Abstract::Vector& getXVec() const;

  //! Get %Arclength parameter
  virtual double getParam() const;

  //! Returns the solution vector component of extended vector
  virtual NOX::Abstract::Vector& getXVec();

  //! Get %Arclength parameter
  virtual double& getParam();

}; // class Vector
} // namespace Continuation
} // namespace LOCA

#endif
