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

#ifndef LOCA_PARAMETER_LIBRARY_H
#define LOCA_PARAMETER_LIBRARY_H

#include "NOX_Common.H"  // for <map>

namespace LOCA {
  
  namespace Parameter {

    // Forward declarations:

    // Abstract data type for each entry in library
    class AbstractEntry;  

    // Data type for each entry templated on value type
    template <class ValueType> class Entry;

    /*! 
     * \brief Class to provide a centralized library for setting/retrieving 
     * numerical parameter values in application codes.
     */
    /*!
     * This class provides a mechanism for setting and retrieving arbitrary
     * numerical parameter values throughout an application code.  Parameters 
     * can be
     * material properties, coefficients in source functions, etc.  The 
     * purpose of this class is to allow external libraries to set and 
     * retrieve parameters values to perform, for example, numerical 
     * continuation and optimization.
     *
     * This class in currently under development and is far from complete.
     */
    class Library {

      //! Map of entries for a parameter name
      typedef map<string, AbstractEntry*> ValueTypeMap;

      //! Const iterator for ValueTypeMap
      typedef ValueTypeMap::const_iterator ValueTypeMapConstIterator;

      //! Iterator for ValueTypeMap
      typedef ValueTypeMap::iterator ValueTypeMapIterator;

      //! Map of maps for all parameter names
      typedef map<string, ValueTypeMap*> ParameterMap;

      //! Const iterator for ParameterMap
      typedef ParameterMap::const_iterator ParameterMapConstIterator;

      //! Iterator for ParameterMap
      typedef ParameterMap::iterator ParameterMapIterator;

    public:
  
      //! Default constructor
      Library() {}

      //! Destructor
      ~Library();

      //! Set parameter given by \em name to value \em value
      template <class ValueType>
      void setParameterValue(const string& name, const ValueType& value);

      //! Add a new parameter to library using the default setting mechanism
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class ObjectType, class ValueType>
      bool addParameterEntry(const string& name, ObjectType& object, 
			     ValueType ObjectType::* object_val_ptr);

      //! Add a new parameter to library using functor setting mechanism
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class FunctorType, class ValueType>
      bool addParameterEntry(const string& name, FunctorType* fctr);

      //! Add a new parameter using custom entry
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class ValueType>
      bool addParameterEntry(const string& name, 
			     Entry<ValueType>* entry);

    protected:

      // Make copy constructor, assignment operator protected to prohibit 
      // copying
  
      //! Copy constructor
      Library(const Library& l);

      //! Assignment operator
      Library& operator = (const Library& l);

      /*! 
       * \brief Finds the set of entries corresponding to parameter \em name 
       * and sets the iterator \em it to point to it.  
       */
      ParameterMapIterator getEntryMapIterator(const string& name);

      /*!
       * \brief Finds the entry corresponding to entry \em valueTypeString for 
       * parameter corresponding to \em paramIterator and sets the iterator
       * \em valueIterator to point to it.  Returns false if there
       * is no entry for type given by \em valueTypeString.
       */
      ValueTypeMapIterator getEntryIterator(
				   const string& valueTypeString, 
				   const ParameterMapIterator& paramIterator);

      /*! 
       * \brief Gets the entry corresponding to parameter name \em name and 
       * type \em ValueType
       */
      template <class ValueType>
      Entry<ValueType>& getEntry(const string& name);

      //! Returns a string representation of type \em ValueType
      template <class ValueType> string getTypeName();

    protected:

      //! Library of Entries
      ParameterMap library;
  
    };
  }
}

// Include template definitions
#include "LOCA_Parameter_LibraryT.H"

#endif
