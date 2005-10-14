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

#ifndef NOX_PARAMETER_LIST_H
#define NOX_PARAMETER_LIST_H

#include "NOX_Common.H"		 // class data element (string)
#include "NOX_Parameter_Entry.H" // class data element 

#include "NOX_Parameter_AnyPtr.H"
#include "NOX_Parameter_Rcp.H"

namespace NOX {
namespace Parameter {

//! Manipulating lists of parameters.
class List {

  //! Parameter container typedef
  typedef map<string, Entry> Map;

  //! Parameter container iterator typedef
  typedef Map::iterator Iterator;

public:

  /* Made this typedef public so that we can use iterators in
   * converting to Teuchos parameter list.  Just like the
   * Teuchos paramter list -- CRJ 21 Sept 2005 */
  //! Parameter container const iterator typedef
  typedef Map::const_iterator ConstIterator;

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

  //! Returns a const reference to the sublist
  /*! 
    If the list does not already exist, throws an error. If the name
    exists but is not a sublist, throws an error.
  */
  const List& sublist(const string& name) const;

  /** @name Setting Parameters 
   
    Sets different types of parameters. The type depends on the second
    entry. Be sure to use static_cast<type>() when the type is
    ambiguous. Both char* and string map to are stored as strings
    internally. Sets the parameter as "unused".
  */
  //@{

  template <typename T> 
  void setParameter(const string& name, T* value);

  template <typename T> 
  void setParameter(const string& name, const T* value);

  template <typename T> 
  void setParameter(const string& name, const Teuchos::RefCountPtr<T>& value);

  void setParameter(const string& name, bool value);
  void setParameter(const string& name, int value);
  void setParameter(const string& name, double value);
  void setParameter(const string& name, const char* value);
  void setParameter(const string& name, char *value);
  void setParameter(const string& name, const string& value);
  void setParameter(const string& name, 
		    const NOX::Parameter::Arbitrary& value);
  //@}

  /** @name Getting Parameters 
   
    Get different types of parameters. The type depends on the second
    entry. Returns the nominal value if that parameter has not been
    specified. The non-const version adds the (name, nominal) pair to
    the list if it's not already specified. Be sure to use
    static_cast<type>() when the type is ambiguous. Both char* and
    string map return string values. Sets the parameters as "used".
  */
  //@{
  bool getParameter(const string& name, bool nominal);
  int getParameter(const string& name, int nominal);
  double getParameter(const string& name, double nominal);
  const string& getParameter(const string& name, const char* nominal);
  const string& getParameter(const string& name, const string& nominal);
  const NOX::Parameter::Arbitrary& 
  getParameter(const string& name, const NOX::Parameter::Arbitrary& nominal);

  bool getParameter(const string& name, bool nominal) const;
  int getParameter(const string& name, int nominal) const;
  double getParameter(const string& name, double nominal) const;
  const string& getParameter(const string& name, const char* nominal) const;
  const string& getParameter(const string& name, const string& nominal) const;
  const NOX::Parameter::Arbitrary& 
  getParameter(const string& name, 
	       const NOX::Parameter::Arbitrary& nominal) const;

  //@}

  /** @name Getting Arbitrary Parameters Without Nominal Value

  Special commands for accessing parameters without specifying a
  nominal value.  Will throw an error if the parameter does not exist
  or is of the wrong type. Recommend called isParameterArbitrary()
  first.
  */
  //@{
  const NOX::Parameter::Arbitrary& 
  getArbitraryParameter(const string& name) const;

  template <typename T> 
  T* getAnyPtrParameter(const string& name) const;

  template <typename T> 
  const T* getAnyConstPtrParameter(const string& name) const;

  template <typename T>
  Teuchos::RefCountPtr<T> getRcpParameter(const string& name) const;

  //@}

  //! Return true if a parameter with this name exists.
  bool isParameter(const string& name) const;

  /** @name Is Parameter of the Specified Type?
   
    Returns true if the specified parameter exists AND is of the
    specified type.
  */
  //@{
  bool isParameterBool(const string& name) const;
  bool isParameterInt(const string& name) const;
  bool isParameterDouble(const string& name) const;
  bool isParameterString(const string& name) const;
  bool isParameterSublist(const string& name) const;
  bool isParameterArbitrary(const string& name) const;
  template <typename T> bool isParameterAnyPtr(const string& name) const;
  template <typename T> bool isParameterRcp(const string& name) const;
  //@}

  /** @name Is Parameter Equal Value?
   
    Returns true if the specified parameter exists AND is equal to the
    specified value.  (Not valid for Arbitrary or List parameters.)
  */
  //@{
  bool isParameterEqual(const string& name, bool value) const;
  bool isParameterEqual(const string& name, int value) const;
  bool isParameterEqual(const string& name, double value) const;
  bool isParameterEqual(const string& name, const char* value) const;
  bool isParameterEqual(const string& name, const string& value) const;
  //@}

  //! Printing 
  ostream& print(ostream& stream, int indent = 0) const;

  /* Added these methods to facilitate conversion to Teuchos parameter lists
   * -- crj 28 Sept 2005 */
  //@{ \name Read-only access to the iterator
  //! An iterator pointing to the first entry
  ConstIterator begin() const ;

  //! An iterator pointing beyond the last entry
  ConstIterator end() const ;


public:

  //! Check to see if "l" or any of its sublists is "this"
  bool isRecursive(const List& l) const;

  //! Access to name (i.e., returns i->first)
  const string& name(ConstIterator i) const;

  //! Access to Entry (i.e., returns i->second)
  Entry& entry(Iterator i);

  //! Access to Entry (i.e., returns i->second)
  const Entry& entry(ConstIterator i) const;

private:

  //! Parameter list
  Map params;
 
  //! Used to create a string when the getParameter is called with a
  //! char* nominal value. A new string is created for each such
  //! argument. The whole group of strings is destroyed when this object
  //! is destroyed. This is really annoying, but I don't know a better 
  //! way.
  mutable vector<string> tmpstrings;
};
}
}

//
// Template definitions here
//

template <typename T> 
void NOX::Parameter::List::setParameter(const string& name, T* value)
{
  setParameter(name, NOX::Parameter::AnyPtr(value));
}

template <typename T> 
void NOX::Parameter::List::setParameter(const string& name, const T* value)
{
  setParameter(name, NOX::Parameter::AnyPtr(value));
}

template <typename T> 
void NOX::Parameter::List::setParameter(const string& name, 
					const Teuchos::RefCountPtr<T>& value)
{
  setParameter(name, NOX::Parameter::Rcp(value));
}

template <typename T> 
T* NOX::Parameter::List::getAnyPtrParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isArbitrary())) {
    const NOX::Parameter::Arbitrary& arb = entry(i).getArbitraryValue();
    const NOX::Parameter::AnyPtr* anyPtr = 
      dynamic_cast<const NOX::Parameter::AnyPtr*>(&arb);
    if (anyPtr != NULL)
      return (*anyPtr).template getPointer<T>();
    else {
      cerr << "NOX::Parameter::List::getAnyPtrParameter - "
	   << "parameter " << name << "is not a pointer" << endl;
      throw "NOX Error";
    }
  }

  cerr << "NOX::Parameter::List::getAnyPtrParameter - no such parameter" 
       << endl;
  throw "NOX Error";
}

template <typename T> 
const T* NOX::Parameter::List::getAnyConstPtrParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isArbitrary())) {
    const NOX::Parameter::Arbitrary& arb = entry(i).getArbitraryValue();
    const NOX::Parameter::AnyPtr* anyPtr = 
      dynamic_cast<const NOX::Parameter::AnyPtr*>(&arb);
    if (anyPtr != NULL)
      return (*anyPtr).template getConstPointer<T>();
    else {
      cerr << "NOX::Parameter::List::getAnyPtrParameter - "
	   << "parameter " << name << "is not a pointer" << endl;
      throw "NOX Error";
    }
  }

  cerr << "NOX::Parameter::List::getAnyPtrParameter - no such parameter" 
       << endl;
  throw "NOX Error";
}

template <typename T> 
Teuchos::RefCountPtr<T> 
NOX::Parameter::List::getRcpParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isArbitrary())) {
    const NOX::Parameter::Arbitrary& arb = entry(i).getArbitraryValue();
    const NOX::Parameter::Rcp* rcpPtr = 
      dynamic_cast<const NOX::Parameter::Rcp*>(&arb);
    if (rcpPtr != NULL)
      return (*rcpPtr).template getRcp<T>();
    else {
      cerr << "NOX::Parameter::List::getRcpParameter - "
	   << "parameter " << name << "is not a ref-count pointer" << endl;
      throw "NOX Error";
    }
  }

  cerr << "NOX::Parameter::List::getRcpParameter - no such parameter" 
       << endl;
  throw "NOX Error";
}

template <typename T> 
bool NOX::Parameter::List::isParameterAnyPtr(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end()) {
    if (entry(i).isArbitrary()) {
      const NOX::Parameter::Arbitrary& arb = entry(i).getArbitraryValue();
      const NOX::Parameter::AnyPtr* anyPtr = 
	dynamic_cast<const NOX::Parameter::AnyPtr*>(&arb);
      if (anyPtr != NULL)
	return (*anyPtr).template isType<T>();
    }
  }

  return false;
}

template <typename T> 
bool NOX::Parameter::List::isParameterRcp(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end()) {
    if (entry(i).isArbitrary()) {
      const NOX::Parameter::Arbitrary& arb = entry(i).getArbitraryValue();
      const NOX::Parameter::Rcp* rcpPtr = 
	dynamic_cast<const NOX::Parameter::Rcp*>(&arb);
      if (rcpPtr != NULL)
	return (*rcpPtr).template isType<T>();
    }
  }

  return false;
}

#endif


