// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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


#ifndef TEUCHOS_PARAMETER_ENTRY_H
#define TEUCHOS_PARAMETER_ENTRY_H

/*! \file Teuchos_ParameterEntry.hpp
    \brief Object held as the "value" in the Teuchos::ParameterList map.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_any.hpp"

/*! \class Teuchos::ParameterEntry
    \brief This object is held as the "value" in the Teuchos::ParameterList map.

    This structure holds a \c Teuchos::any value and information on the status of this
    parameter (isUsed, isDefault, etc.).  The type of parameter is chosen through the
    templated Set/Get methods.
*/

namespace Teuchos {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
class ParameterList; // another parameter type (forward declaration)
#endif

class ParameterEntry {

public:
  //@{ \name Constructors/Destructor

  //! Default Constructor
  ParameterEntry();
  
  //! Copy constructor
  ParameterEntry(const ParameterEntry& source);

  //! Templated constructor
  template<typename T>
  ParameterEntry(T value, bool isDefault = false)
  : val_(value),
    isUsed_(false),
    isDefault_(isDefault)
  {}

  //! Destructor
  ~ParameterEntry() {};

  //@}

  //@{ \name Set Methods

  //! Replace the current parameter entry with \c source.
  ParameterEntry& operator=(const ParameterEntry& source);

  /*! \brief Templated set method that uses the input value type to determine the type of parameter.  
      
      \note <ul>
	    <li> Invalidates any previous values stored by this object although it doesn't necessarily erase them.  
            <li> Resets 'isUsed' functionality.  
	    </ul>
  */
  template<typename T>
  void setValue(T value, bool isDefault = false)
  {
    val_ = value;
    isDefault_ = isDefault;
  }

  //! Create a parameter entry that is an empty list.
  ParameterList& setList(bool isDefault = false);

  //@}

  //@{ \name Get Methods 
   
  /*! \brief Templated get method that uses the input pointer type to determine the type of parameter to return.  

      \note This method will cast the value to the type requested.  If that type is incorrect, 
	    an exception will be thrown by the any_cast.
  */
  template<typename T>
  T& getValue(T *ptr) const
  {
    isUsed_ = true;
    return const_cast<T&>(Teuchos::any_cast<T>( val_ ));
  }

  //@}

  //@{ \name Attribute Methods  
  
  //! Return whether or not the value has been used; i.e., whether or not the value has been retrieved via a get function.
  bool isUsed() const { return isUsed_; }

  //! Return whether or not the value itself is a list.
  bool isList() const { return isList_; }
  //@}

  //@{ \name I/O Methods

  /*! \brief Output a non-list parameter to the given output stream.  

      The parameter is followed by "[default]" if it is the default value given through a 
      Set method.  Otherwise, if the parameter was unused (not accessed through a Get method), 
      it will be followed by "[unused]".  This function is called by the "ostream& operator<<". 
  */
  ostream& leftshift(ostream& os) const;
  //@}

private:

  //! Reset the entry
  void reset();
  
  //! Templated Datatype
  any val_;

  //! List flag
  bool isList_;

  //! Has this parameter been accessed by a "get" function?
  mutable bool isUsed_;

  //! Was this parameter a default value assigned by a "get" function?
  mutable bool isDefault_;

};

/*! \relates ParameterEntry 
    \brief A templated helper function for returning the value of type \c T held in the ParameterEntry object,
    where the type \c T can be specified in the call.  This is an easier way to call the getValue method
    in the ParameterEntry class, since the user does not have to pass in a pointer of type \c T.
*/
template<typename T>
T& getValue( const ParameterEntry &entry )
{
  return entry.getValue((T*)NULL);
}

/*! \relates ParameterEntry 
    \brief Output stream operator for handling the printing of parameter entries.  
*/
inline ostream& operator<<(ostream& os, const ParameterEntry& e) 
{ 
  return e.leftshift(os);
}

} // namespace Teuchos


#endif
