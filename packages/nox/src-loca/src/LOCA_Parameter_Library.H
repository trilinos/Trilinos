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

#ifndef LOCA_PARAMETER_LIBRARY_H
#define LOCA_PARAMETER_LIBRARY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Common.H"
#include <map>

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
      typedef std::map<std::string, AbstractEntry*> ValueTypeMap;

      //! Const iterator for ValueTypeMap
      typedef ValueTypeMap::const_iterator ValueTypeMapConstIterator;

      //! Iterator for ValueTypeMap
      typedef ValueTypeMap::iterator ValueTypeMapIterator;

      //! Map of maps for all parameter names
      typedef std::map<std::string, ValueTypeMap*> ParameterMap;

      //! Const iterator for ParameterMap
      typedef ParameterMap::const_iterator ParameterMapConstIterator;

      //! Iterator for ParameterMap
      typedef ParameterMap::iterator ParameterMapIterator;

    public:
  
      //! Default constructor
      Library() : library() {}

      //! Destructor
      ~Library();

      //! Set parameter given by \em name to value \em value
      template <class ValueType>
      void setValue(const std::string& name, const ValueType& value);

      //! Get parameter given by \em name
      template <class ValueType>
      ValueType getValue(const std::string& name) const;

      //! Add a new parameter to library using the default setting mechanism
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class ObjectType, class ValueType>
      bool addParameterEntry(const std::string& name, ObjectType& object, 
			     ValueType ObjectType::* object_val_ptr);

      //! Add a new parameter to library using functor setting mechanism
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class FunctorType, class ValueType>
      bool addParameterEntry(const std::string& name, FunctorType* fctr);

      //! Add a new parameter using custom entry
      /*!
       * Returns true if successful in adding entry to library, false 
       * otherwise.
       */
      template <class ValueType>
      bool addParameterEntry(const std::string& name, 
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
      ParameterMapIterator getEntryMapIterator(const std::string& name);

      /*! 
       * \brief Finds the set of entries corresponding to parameter \em name 
       * and sets the iterator \em it to point to it.  
       */
      ParameterMapConstIterator getEntryMapIterator(const std::string& name) const;

      /*!
       * \brief Finds the entry corresponding to entry \em valueTypeString for 
       * parameter corresponding to \em paramIterator and sets the iterator
       * \em valueIterator to point to it.  Returns false if there
       * is no entry for type given by \em valueTypeString.
       */
      ValueTypeMapIterator getEntryIterator(
				   const std::string& valueTypeString, 
				   const ParameterMapIterator& paramIterator);

      /*!
       * \brief Finds the entry corresponding to entry \em valueTypeString for 
       * parameter corresponding to \em paramIterator and sets the iterator
       * \em valueIterator to point to it.  Returns false if there
       * is no entry for type given by \em valueTypeString.
       */
      ValueTypeMapConstIterator getEntryIterator(
			 const std::string& valueTypeString, 
			 const ParameterMapConstIterator& paramIterator) const;

      /*! 
       * \brief Gets the entry corresponding to parameter name \em name and 
       * type \em ValueType
       */
      template <class ValueType>
      Entry<ValueType>& getEntry(const std::string& name);

      /*! 
       * \brief Gets the entry corresponding to parameter name \em name and 
       * type \em ValueType
       */
      template <class ValueType>
      const Entry<ValueType>& getEntry(const std::string& name) const;

      //! Returns a std::string representation of type \em ValueType
      template <class ValueType> std::string getTypeName() const;

    protected:

      //! Library of Entries
      ParameterMap library;
  
    };
  }
}

// Include template definitions
#include "LOCA_Parameter_LibraryT.H"

#endif
