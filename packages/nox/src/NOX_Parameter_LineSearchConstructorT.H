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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_PARAMETER_LINESEARCHCONSTRUCTORT_H
#define NOX_PARAMETER_LINESEARCHCONSTRUCTORT_H

#include "NOX_Common.H"		// class data element (string)
#include "NOX_Parameter_LineSearchConstructor.H"

namespace NOX {

namespace Parameter {

/*! 
  \brief Templated instantiation of NOX::Parameter::LineSearchConstructor.

  This templated class provides an instantiation of
  NOX::Parameter::LineSearchConstructor for any line search class derived
  from NOX::LineSearch::Generic.

  The name of the line search can optionally be included as an argument
  to the constructor.
*/
template <class T> class LineSearchConstructorT : public LineSearchConstructor
{
  
public:
  
  //! Default Constructor. Sets the name to "Templated LineSearch Constructor".
  LineSearchConstructorT() : name("Templated LineSearch Constructor") {};

  //! Constructor with the line search name.
  LineSearchConstructorT(const string& n) : name(n) {};
  
  //! Destructor
  ~LineSearchConstructorT() {};

  //derived
  LineSearchConstructor* clone() const { return new LineSearchConstructorT<T>(name); };

  //derived
  const string& getType() const { return name; };

  //derived
  ostream& print(ostream& stream, int indent = 0) const {return stream;} ;

  //! Constructs and linesearch of type T and returns its pointer.
  NOX::LineSearch::Generic* newLineSearch(const NOX::Utils& u, NOX::Parameter::List& p) const { return new T(u,p); };

private:

  //! Descriptive name of the linesearch that this object can generate.
  const string name;

};

} // namespace Parameter
} // namespace NOX

#endif
