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

#ifndef NOX_PARAMETER_LIST_H
#define NOX_PARAMETER_LIST_H

#include "NOX_Common.H"		 // class data element (string)
#include "NOX_Parameter_Entry.H" // class data element 

namespace NOX {
namespace Parameter {

//! Manipulating lists of parameters.
class List {

  //! Parameter container typedef
  typedef map<string, Entry> ParamContainer;

  //! Parameter container const iterator typedef
  typedef ParamContainer::const_iterator PCConstIterator;

  //! Parameter container iterator typedef
  typedef ParamContainer::iterator PCIterator;

public:

  //! Constructor
  List();

  //! Copy Constructor
  List(const List& source);

  //! Copy
  List& operator=(const List& source);

  //! Deconstructor
  ~List();

  //! List %unused parameters
  void unused() const;

  //! Creates and empty sublist and returns a reference to the
  //! sublist. If the list already exists, returns reference to that
  //! sublist. If the name exists but is not a sublist, throws an error.
  List& sublist(const string& name);

  /** @name Setting Parameters 
   *
   * Sets different types of parameters. The type depends on the
   * second entry. Be sure to use static_cast<type>() when the type is
   * ambiguous. Both char* and string map to are stored as strings
   * internally. Sets the parameter as "unused". */
  //@{
  void setParameter(const string& name, bool value);
  void setParameter(const string& name, int value);
  void setParameter(const string& name, double value);
  void setParameter(const string& name, const char* value);
  void setParameter(const string& name, const string& value);
  //@}

  /** @name Getting Parameters 
   *
   * Get different types of parameters. The type depends on the second
   * entry. Returns the nominal value if that parameter has not been *
   * specified. Be sure to use static_cast<type>() when the type is *
   * ambiguous. Both char* and string map return string values.  Sets
   * the parameters as "used". */
  //@{
  bool getParameter(const string& name, bool nominal) const;
  int getParameter(const string& name, int nominal) const;
  double getParameter(const string& name, double nominal) const;
  const string& getParameter(const string& name, const char* nominal) const;
  const string& getParameter(const string& name, const string& nominal) const;
  //@}

  //! Return true if a parameter with this name exists.
  bool isParameter(const string& name) const;


  /** @name Is Parameter Equal Value?
   *
   * Returns true if the specified parameter exists AND is equal to
   * the specified value.  (Not value for List
   * parameters.) */
  //@{
  bool isParameterEqual(const string& name, bool value) const;
  bool isParameterEqual(const string& name, int value) const;
  bool isParameterEqual(const string& name, double value) const;
  bool isParameterEqual(const string& name, const char* value) const;
  bool isParameterEqual(const string& name, const string& value) const;
  //@}

  //! Printing 
  ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Check to see if "l" or any of its sublists is "this"
  bool isRecursive(const List& l) const;


private:

  //! Parameter list
  ParamContainer params;
 
  //! Used to create a string when the getParameter is called with a
  //! char* nominal value. A new string is created for each such
  //! argument. The whole group of strings is destroyed when this object
  //! is destroyed. This is really annoying, but I don't know a better 
  //! way.
  mutable vector<string> tmpstrings;
};
}
}

#endif


