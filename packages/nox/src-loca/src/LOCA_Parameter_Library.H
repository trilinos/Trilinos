// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
