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

#ifndef NOX_PARAMETER_ANYPTR_H
#define NOX_PARAMETER_ANYPTR_H

#include "NOX_Common.H" // for NOX_Config.h

#ifdef HAVE_NOX_ANY

#include "NOX_Parameter_Arbitrary.H"		// base class
#include "Teuchos_any.hpp"                      // for "any" class

namespace NOX {

  namespace Parameter {

    //! Implementation of NOX::Parameter::Arbitrary to hold pointers
    /*!
     * %AnyPtr is an implementation of NOX::Parameter::Arbitrary to allow
     * passing of pointers through parameter lists.  It uses Teuchos::any
     * to store the pointer in a type-safe manner.
     */
    class AnyPtr : public NOX::Parameter::Arbitrary {

    public:
  
      //! Constructor, non-const reference
      template <class T> AnyPtr(T* ptr) 
	: data(ptr), 
	  type_string("Generic AnyPtr Container") {}

      //! Constructor, const reference
      template <class T> AnyPtr(const T* ptr) 
	: data(ptr), 
	  type_string("Generic AnyPtr Container") {}

      AnyPtr(const Teuchos::any& otherData)
        : data(otherData), type_string("Generic AnyPtr Container") {}

      //! Copy constructor
      AnyPtr(const AnyPtr& src)
	: data(src.data) , type_string(src.type_string) {}
  
      //! Destructor
      virtual ~AnyPtr() {}

      //! Clone a exact replica of yourself and pass back a pointer.
      virtual Arbitrary* clone() const {
	return new AnyPtr(*this); }

      //! Get a short descriptive string describing the %Arbitrary object
      virtual const string& getType() const {
	return type_string; }

      //! Print out detailed information describing the %Arbitrary object
      virtual ostream& print(ostream& stream, int indent = 0) const {
	for (int i=0; i<indent; i++)
	  stream << " ";
	return stream << data << endl; } 

      //! Return whether stored data type is T
      template <class T> bool isType() const {
	return (data.type() == typeid(T));
      }

      //! Return stored pointer
      template <class T> T* getPointer() const {
	return Teuchos::any_cast<T*>(data);
      }

      //! Return stored pointer
      template <class T> const T* getConstPointer() const {
	return Teuchos::any_cast<T*>(data);
      }

    protected:
      
      //! Teuchos any container
      Teuchos::any data;

      //! Type name for Arbitrary
      string type_string;
   
    }; // class AnyPtr
  } // namespace Parameter
} // namespace NOX

#endif // HAVE_NOX_ANY

#endif // NOX_PARAMETER_ANYPTR_H
