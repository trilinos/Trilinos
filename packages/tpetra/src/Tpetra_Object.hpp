// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_OBJECT_HPP
#define TPETRA_OBJECT_HPP

#include "Tpetra_ConfigDefs.hpp" // for iostream and string

namespace Tpetra
{

//! Tpetra::Object:  The base Tpetra class.
/*! The Object class provides capabilities common to all Tpetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
*/

class Object
{
  public:
  //@{ \name Constructors/destructor.
  //! Object Constructor.
  /*! Object is the primary base class in Tpetra.  All Tpetra class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Object(int tracebackModeIn = -1) 
		: label_("Tpetra::Object")
	{
		tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
	};

  //! Object Constructor.
  /*! Creates an Object with the given label.
  */
  Object(std::string const& label, int tracebackModeIn = -1)
		: label_(label)
	{
		tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
	};

  //! Object Copy Constructor.
  /*! Makes an exact copy of an existing Object instance.
  */
  Object(Object const& obj)
		: label_(obj.label_)
	{
	};

  //! Object Destructor.
  /*! Completely deletes an Object object.  
  */
  virtual ~Object()
	{};
	
  //@}
  
  //@{ \name Attribute set/get methods.

  //! Object Label definition using string.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void setLabel(std::string const label)
	{ 
		label_ = label;
  }

  //! Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual std::string label() const
	{
		return(label_);
	};

  //! Set the value of the Object error traceback report mode.
  /*! Sets the integer error traceback behavior.  
      TracebackMode controls whether or not traceback information is printed when run time 
      integer errors are detected:

      <= 0 - No information report

       = 1 - Fatal (negative) values are reported

      >= 2 - All values (except zero) reported.

      Default is set to 1.
  */
  static void setTracebackMode(int tracebackModeValue)
	{
		if(tracebackModeValue < 0)
			tracebackModeValue = 0;
		Object tempObject(tracebackModeValue);
	};
	
  //! Get the value of the Object error report mode.
  static int getTracebackMode()
	{
		int temp = Object::tracebackMode;
		if(temp == -1)
			temp = Tpetra_DefaultTracebackMode;
		return(temp);
	};
	
  //@}

  //@{ \name Miscellaneous

  //! Print object to an output stream
  //! Print method
  virtual void print(ostream& os) const
	{
		  // os << label_; // No need to print label, since ostream does it already
	};

	//! Error reporting method.
	virtual int reportError(std::string const message, int errorCode) const {
		cerr << endl << "Error in Tpetra Object with label: " << label_ << endl 
				 << "Error Message:  " << message << "  Error Code:  " << errorCode << endl;
		return(errorCode);
	}
  
// tracebackMode controls how much traceback information is printed when run time 
// integer errors are detected:
// = 0 - No information report
// = 1 - Fatal (negative) values are reported
// = 2 - All values (except zero) reported.

// Default is set to 2.  Can be set to different value using setTracebackMode() method in
// Object class
  static int tracebackMode;
	
private:
	
  std::string label_;

	//! Assignment operator (declared but not defined, do not use)
	Object& operator = (Object const& Source);
	
}; // class Object

// Set TracebackMode value to default
int Object::tracebackMode(-1);

inline ostream& operator<<(ostream& os, Tpetra::Object const& Obj)
{
  os << Obj.label() << endl;
  Obj.print(os); 
  return(os);
}

} // namespace Tpetra


#endif // TPETRA_OBJECT_HPP
