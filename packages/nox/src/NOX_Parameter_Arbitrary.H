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

#ifndef NOX_PARAMETER_ARBITRARY_H
#define NOX_PARAMETER_ARBITRARY_H

#include "NOX_Common.H"		// class data element (string)

namespace NOX {

namespace Parameter {

/*! 
  \brief Support for Arbitrary parameters. 

  Any object which derives from this class may be passed as a
  parameter via a NOX::Parameter::List.

  For example, we include a very basic instantiation to support the
  passing of a vector<double> object.

\code

class DoubleVectorParameter : public NOX::Parameter::Arbitrary {

public:

DoubleVectorParameter(const vector<double>& x) { v = x; } ;
~DoubleVectorParameter() {};
Arbitrary* clone() const { return new DoubleVectorParameter(v); };
const string& getType() const { return "Double Vector"; };
const vector<double>& getDoubleVector() const {return v;};

private:

vector<double> v;

};

\endcode

This paramater would be stored in a parameter list as follows.

\code

vector<double> x;

// .... Fill in x here ...

NOX::Parameter::List params;
params.setParameter("My List of Doubles", DoubleVectorParameter(x));

\endcode

This parameter could later be accessed as follows.

\code

function foo (NOX::Parameter::List& params)
{
  const NOX::Parameter::Arbitrary& a = params.getArbitraryParameter("My List of Doubles");
  const DoubleVectorParameter& d = dynamic_cast<const DoubleVectorParameter&>(a);
  vector<double> x = d.getDoubleVector();

  // ... Use x here ...
}

\endcode

This is a very simplistic example of what one could do with an
Arbitrary parameter. For example, even in this case we are likely
doing far too much copying of the vectors. Ideally, we would use a
shared_ptr object to store the object and not actually copy it every
time we did a clone.

*/
class Arbitrary {

public:
  
  //! Default Constructor
  Arbitrary() {};
  
  //! Destructor
  virtual ~Arbitrary() {};

  //! Clone a exact replica of yourself and pass back a pointer.
  /*!  It's left to the implementer to decide precisely how clone
    should function. The main requirement is that the cloned copy
    should still be valid even if the original is destroyed and even
    if multiple cloned copies exist.
   */
  virtual Arbitrary* clone() const = 0;

  //! Get a short descriptive string describing Arbitrary object
  virtual const string& getType() const = 0;

  //! Print out detailed information describing the Arbitrary object
  /*!
    Each line of output should be indented by the number of spaces specified by \c indent.
  */
  virtual ostream& print(ostream& stream, int indent = 0) const {return stream;} ;

};

} // namespace Parameter
} // namespace NOX

#endif
