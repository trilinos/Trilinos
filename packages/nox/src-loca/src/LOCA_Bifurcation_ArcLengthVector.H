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

#include "NOX_Abstract_Vector.H"  // Base class

//! Nonlinear solvers package namespace
namespace LOCA {

//! %NOX  bifurcation vector and group
namespace Bifurcation {

/*! 
  \brief %LOCA's pure abstract vector interface for vectors that are
  used by the NLS solver.  Derived from the NOX::Abstract::Vector class.

  The documentation throughout assumes that \f$x\f$ denotes "this"
  vector, \f$x_i\f$ is the i-th component (or entry) of \f$x\f$, and
  \f$n\f$ is the length of \f$x\f$.

  \author Roger Pawlowski (SNL 9233)
*/

class ArcLengthVector : public NOX::Abstract::Vector {

public:
  
  //! %Vector constructor.
  ArcLengthVector(const NOX::Abstract::Vector& xVec,
   	          double arcParam);

  //! %Copy constructor.
  ArcLengthVector(const ArcLengthVector& source);

  //! %Vector destructor.
  virtual ~ArcLengthVector();

  virtual NOX::Abstract::Vector& init(double gamma);

  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

  virtual LOCA::Bifurcation::ArcLengthVector& operator=(const LOCA::Bifurcation::ArcLengthVector& y);

  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector& scale(double gamma);

  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

  virtual NOX::Abstract::Vector& update(double alpha,
                                        const NOX::Abstract::Vector& a,
                                        double gamma = 0.0);

  virtual NOX::Abstract::Vector& update(double alpha,
                         const NOX::Abstract::Vector& a, double beta,
                         const NOX::Abstract::Vector& b, double gamma = 0.0);

  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

  virtual double norm(NormType type = TwoNorm) const;

  virtual double norm(const NOX::Abstract::Vector& weights) const;

  virtual double dot(const NOX::Abstract::Vector& y) const;
  
  virtual int length() const;

  virtual bool print() const;

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

protected:

  //! Solution vector
  NOX::Abstract::Vector* xVectorPtr;

  //! %Arclength parameter that's part of the augmented solution vector
  double arcParam;

}; // class Vector
} // namespace Bifurcation
} // namespace LOCA

#endif
