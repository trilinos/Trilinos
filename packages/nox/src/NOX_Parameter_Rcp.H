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

#ifndef NOX_PARAMETER_RCP_H
#define NOX_PARAMETER_RCP_H

#include "NOX_Common.H" // for NOX_Config.h

#ifdef HAVE_NOX_ANY

#include "NOX_Parameter_Arbitrary.H"		// base class
#include "Teuchos_any.hpp"                      // for "any" class
#include "Teuchos_RefCountPtr.hpp"              // for ref-count pointer

namespace NOX {

  namespace Parameter {

    //! Implementation of NOX::Parameter::Arbitrary to ref-count pointers
    /*!
     * %Rcp is an implementation of NOX::Parameter::Arbitrary to allow
     * passing of ref-count pointers through parameter lists.  It uses 
     * Teuchos::any to store the ref-count pointer in a type-safe manner.
     */
    class Rcp : public NOX::Parameter::Arbitrary {

    public:
  
      //! Constructor, non-const reference
      template <typename T> Rcp(const Teuchos::RefCountPtr<T>& ptr) 
	: data(ptr), 
	  type_string("Generic Ref-Count Pointer Container") {}

      //! Constructor taking Teuchos::any
      Rcp(const Teuchos::any& otherData)
        : data(otherData), type_string("Generic Ref-Count Pointer Container") 
      {}

      //! Copy constructor
      Rcp(const Rcp& src)
	: data(src.data) , type_string(src.type_string) {}
  
      //! Destructor
      virtual ~Rcp() {}

      //! Clone a exact replica of yourself and pass back a pointer.
      virtual Arbitrary* clone() const {
	return new Rcp(*this); }

      //! Get a short descriptive string describing the %Arbitrary object
      virtual const string& getType() const {
	return type_string; }

      //! Print out detailed information describing the %Arbitrary object
      virtual ostream& print(ostream& stream, int indent = 0) const {
	for (int i=0; i<indent; i++)
	  stream << " ";
	return stream << data << endl; } 

      //! Return whether stored data type is RefCountPtr<T>
      template <typename T> bool isType() const {
	return (data.type() == typeid(Teuchos::RefCountPtr<T>));
      }

      //! Return stored pointer
      template <typename T> Teuchos::RefCountPtr<T> getRcp() const {
	return Teuchos::any_cast< Teuchos::RefCountPtr<T> >(data);
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
