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

#ifndef NOX_PARAMETER_DIRECTIONCONSTRUCTORT_H
#define NOX_PARAMETER_DIRECTIONCONSTRUCTORT_H

#include "NOX_Common.H"		// class data element (string)
#include "NOX_Parameter_DirectionConstructor.H"

namespace NOX {

namespace Parameter {

/*! 
  \brief Templated instantiation of NOX::Parameter::DirectionConstructor.

  This templated class provides an instantiation of
  NOX::Parameter::DirectionConstructor for any direction class derived
  from NOX::Direction::Generic.

  The name of the direction can optionally be included as an argument
  to the constructor.
*/
template <class T> class DirectionConstructorT : public DirectionConstructor
{
  
public:
  
  //! Default Constructor. Sets the name to "Templated Direction Constructor".
  DirectionConstructorT() : name("Templated Direction Constructor") {};

  //! Constructor with the direction name.
  DirectionConstructorT(const string& n) : name(n) {};
  
  //! Destructor
  ~DirectionConstructorT() {};

  //derived
  DirectionConstructor* clone() const { return new DirectionConstructorT<T>(name); };

  //derived
  const string& getType() const { return name; };

  //derived
  ostream& print(ostream& stream, int indent = 0) const {return stream;} ;

  //! Constructs and direction of type T and returns its pointer.
  NOX::Direction::Generic* newDirection(const NOX::Utils& u, NOX::Parameter::List& p) const { return new T(u,p); };

private:

  //! Descriptive name of the direction that this object can generate.
  const string name;

};

} // namespace Parameter
} // namespace NOX

#endif
