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

#ifndef LOCA_MULTIVECTOR_H
#define LOCA_MULTIVECTOR_H

#include "NOX_Common.H"  // for <iostream>, <vector>, and <string>
#include "NOX_Abstract_Vector.H"  // Base class

//! Nonlinear solvers package namespace
namespace LOCA {

class MultiVector : public NOX::Abstract::Vector {

public:
  
  //! %MultiVector constructor.
  MultiVector(const vector< NOX::Abstract::Vector* > vv,
	      const vector< double > vs);

  //! %Copy constructor.
  MultiVector(const MultiVector& source, NOX::CopyType type = NOX::DeepCopy);

  //! %Vector destructor.
  virtual ~MultiVector();

  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
  virtual MultiVector& operator=(const MultiVector& y);

  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

  virtual NOX::Abstract::Vector& init(double gamma);

  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

  virtual NOX::Abstract::Vector& scale(double gamma);

  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

  virtual NOX::Abstract::Vector& update(double alpha,
                                        const NOX::Abstract::Vector& a,
                                        double gamma = 0.0);

  virtual NOX::Abstract::Vector& update(double alpha,
                         const NOX::Abstract::Vector& a, double beta,
                         const NOX::Abstract::Vector& b, double gamma = 0.0);

  virtual double norm(NormType type = TwoNorm) const;

  virtual double norm(const NOX::Abstract::Vector& weights) const;

  virtual double dot(const NOX::Abstract::Vector& y) const;
  
  virtual int length() const;

  virtual void print() const;

  //! Sets the ith vector
  virtual void setVector(int i, const NOX::Abstract::Vector& v);

  //! Sets the ith scalar
  virtual void setScalar(int i, double s);

  //! Returns const reference to the ith vector
  virtual const NOX::Abstract::Vector& getVector(int i) const;

  //! Returns reference to the ith vector
  virtual NOX::Abstract::Vector& getVector(int i);

  //! Returns copy of the ith scalar
  virtual double getScalar(int i) const;

  //! Returns reference to the ith scalar
  virtual double& getScalar(int i);

protected:

  //! %MultiVector constructor.
  MultiVector(int nvecs, int nscalars);

protected:

  //! Array of vector points
  vector< NOX::Abstract::Vector* > vectorPtrs;

  //! Array of doubles
  vector< double > scalars;

}; // class MultiVector
} // namespace LOCA

#endif
