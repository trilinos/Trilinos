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

#ifndef LOCA_PARAMETER_VECTOR_H
#define LOCA_PARAMETER_VECTOR_H

#include "NOX_Common.H"  // for <iostream>, <vector>, and <string>

namespace LOCA {

/*! 
  \brief %LOCA's container for holding a set of parameters that are
  used by the LOCA continuation routines.

  \author Roger Pawlowski (SNL 9233)
*/

class ParameterVector {

public:
  
  //! Construct a vector of length n. Can optionally be initialized to a scalar value.
  ParameterVector();

  //! Copy constructor
  ParameterVector(const ParameterVector& source);

  //! Clone
  ParameterVector* clone() const;

  //! Destructor
  ~ParameterVector();

  //! Adds a parameter to the list.  Returns the index value assigned to the parameter.
  int addParameter(string label, double value = 0.0);

  //! Initialize the vector.  Returns true if successful.
  bool init(double value);

  //! Scales the entire vector by value.  Returns true if successful.
  bool scale(double value);

  //! Scales the vactor with another vector (element-wise multiply).  Returns true if successful.
  bool scale(const ParameterVector& p);

  //! Updates the parameter vector: this = alpha * alphaVector + b * this.  Returns true if successful.
  bool update(double alpha, const ParameterVector& alphaVector, double b);

  //! Set a parameter vector equal to another. This will erase any old parameters in the object. 
  ParameterVector& operator=(const ParameterVector& y);

  //! Return a reference to the i-th element.  Will throw an error if index is out of range.
  double& operator[] (unsigned int i);

  //! Return a reference to the i-th element (no l-value type).  Will throw an error if index is out of range.
  const double& operator[] (unsigned int i) const;

  //! Set the value of the parameter with index i.  Will throw an error if index is out of range.
  void setValue(unsigned int i, double value);

  //! Set the value of the parameter with the corresponding label.  Will throw an error if "label" is not valid.
  void setValue(string label, double value);

  //! Returns the value of the parameter with index i.  Will throw an error if index is out of range.
  double getValue(unsigned int i) const;

  //! Returns the value of the parameter with the corresponding label.  Will throw an error if "label" is not valid.
  double getValue(string label) const;

  //! Returns the index of the parameter with the corresponding label.  Returns a -1 if "label" is not found.
  int getIndex(string label) const;

  //! Returns a pointer to a C-style array of the parameter values.
  double* getDoubleArrayPointer();

  //! Returns true if the parameter string "label" corresponds to a parameter label in the object.
  bool isParameter(string label) const;

  //! Returns the label of the parameter with index i.
  string getLabel(unsigned int i) const;

  //! Returns the length of parameter vector.
  int length() const;

  //! Prints the vector to cout.
  void print(ostream& stream) const;

protected:

  //! Parameter values
  vector<double> x;

  //! Parameter names
  vector<string> l;

};
} //namespace LOCA

//! Overload the operator<< to print ParameterVector objects.
ostream& operator<<(ostream& stream, const LOCA::ParameterVector& p);

#endif
