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

#ifndef LOCA_BIFURCATION_TPBORDVECTOR_H
#define LOCA_BIFURCATION_TPBORDVECTOR_H

#include "LOCA_MultiVector.H"  // Base class

//! Nonlinear solvers package namespace
namespace LOCA {

//! %NOX  bifurcation vector and group
namespace Bifurcation {

class ArcLengthVector : public LOCA::MultiVector {

public:
  
  //! %Vector constructor.
  ArcLengthVector(const NOX::Abstract::Vector& xVec,
   	          double arcParam);

  //! %Copy constructor.
  ArcLengthVector(const ArcLengthVector& source, NOX::CopyType type = NOX::DeepCopy);

  //! %Vector destructor.
  virtual ~ArcLengthVector();

  virtual LOCA::Bifurcation::ArcLengthVector& operator=(const LOCA::Bifurcation::ArcLengthVector& y);

  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Sets the %ArcLengthVector by setting its two components 
  virtual void setVec(const NOX::Abstract::Vector& xVec,
                      double arcPar);

  //! Returns the solution vector component of extended vector
  virtual const NOX::Abstract::Vector& getXVec() const;

  //! Get %Arclength parameter
  virtual double getArcParam() const;

  //! Returns the solution vector component of extended vector
  virtual NOX::Abstract::Vector& getXVec();

  //! Get %Arclength parameter
  virtual double& getArcParam();

}; // class Vector
} // namespace Bifurcation
} // namespace LOCA

#endif
