// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#ifndef LOCA_PARAMETER_VECTOR_H
#define LOCA_PARAMETER_VECTOR_H

#include "NOX_Common.H"  // for <iostream>, <vector>, and <string>
#include <string>
#include <vector>
#include <ostream>

// Forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

/*! 
  \brief %LOCA's container for holding a set of parameters that are
  used by the LOCA continuation routines.

  \author Roger Pawlowski (SNL 9233)
*/

class ParameterVector {

public:
  
  //! Constructor.
  ParameterVector();

  //! Copy constructor
  ParameterVector(const ParameterVector& source);

  //! Clone
  ParameterVector* clone() const;

  //! Destructor
  ~ParameterVector();

  //! Adds a parameter to the list.  Returns the index value assigned to the parameter.
  int addParameter(std::string label, double value = 0.0);

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
  void setValue(std::string label, double value);

  //! Returns the value of the parameter with index i.  Will throw an error if index is out of range.
  double getValue(unsigned int i) const;

  //! Returns the value of the parameter with the corresponding label.  Will throw an error if "label" is not valid.
  double getValue(std::string label) const;

  //! Returns the index of the parameter with the corresponding label.  Returns a -1 if "label" is not found.
  int getIndex(std::string label) const;

  //! Returns a pointer to a C-style array of the parameter values.
  double* getDoubleArrayPointer();

  //! Returns true if the parameter std::string "label" corresponds to a parameter label in the object.
  bool isParameter(std::string label) const;

  //! Returns the label of the parameter with index i.
  std::string getLabel(unsigned int i) const;

  //! Returns the length of parameter vector.
  int length() const;

  //! Prints the vector to cout.
  void print(std::ostream& stream) const;

  //! Accessor to get the underlying stl vector with all parameter values
  const std::vector<double>& getValuesVector() const;

  //! Accessor to get the underlying stl vector with all parameter names
  const std::vector<std::string>& getNamesVector() const;

protected:

  //! Parameter values
  std::vector<double> x;

  //! Parameter names
  std::vector<std::string> l;

};

//! Overload the operator<< to print ParameterVector objects.
std::ostream& operator<<(std::ostream& stream, const LOCA::ParameterVector& p);

} //namespace LOCA

#endif
