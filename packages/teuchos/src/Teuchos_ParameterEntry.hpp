
#ifndef TEUCHOS_PARAMETER_ENTRY_H
#define TEUCHOS_PARAMETER_ENTRY_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_any.hpp"

/*! \file Teuchos_ParameterEntry.hpp
    \brief Object held as the "value" in the Teuchos::ParameterList map.
*/

namespace Teuchos {

class ParameterList; // another parameter type (forward declaration)

/*! \class ParameterEntry
    \brief This object is held as the "value" in the Teuchos::ParameterList map.
    This structure holds a \c Teuchos::any value and information on the status of this
    parameter (isUsed, isDefault, etc.).  The type of parameter is chosen through the
    templated Set/Get methods.
*/

class ParameterEntry {

public:
  //@{ \name Constructor/Destructor methods

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

  //@{ \name Set methods.

  //! Replace the current parameter entry with \c source.
  ParameterEntry& operator=(const ParameterEntry& source);

  //! Templated set method that uses the input value type to determine the
  //! type of parameter.  Invalidates any previous values stored by this object
  //! although it doesn't necessarily erase them.  Resets 'isUsed' functionality.
  
  template<typename T>
  void setValue(T value, bool isDefault = false)
  {
    val_ = value;
    isDefault_ = isDefault;
  }

  //! Create a parameter entry that is an empty list.
  ParameterList& setList(bool isDefault = false);

  //@}

  //@{ \name Get methods. 
   
  //! Templated get method that uses the input pointer type to determine the
  //! type of parameter to return.  The returned value is nonsense if we do not
  //! request the correct type of value. 

  template<typename T>
  T& getValue(T *ptr) const
  {
    isUsed_ = true;
    return const_cast<T&>(Teuchos::any_cast<T>( val_ ));
  }

  //@}

  //@{ \name Is methods.  
  
  //! Return whether or not the value is used; i.e., whether or not
  //! the value has been retrieved via a get function.

  bool isUsed() const { return isUsed_; };

  //! Return whether or not the value itself is a list.
  bool isList() const { return isList_; };
  //@}

  //@{ \name I/O Methods

  //! Output a non-list parameter to the given output stream.  The parameter is followed by
  //! "[default]" if it is the default value given through a Set method.  Otherwise,
  //! if the parameter was unused (not accessed through a Get method), it will be 
  //! followed by "[unused]".  This function is called by the "ostream& operator<<". 

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
    A templated helper function for returning the value of type \c T held in the ParameterEntry object,
    where the type \c T can be specified in the call.  This is an easier way to call the getValue method
    in the ParameterEntry class, since the user does not have to pass in a pointer of type \c T.
*/
template<typename T>
T& getValue( const ParameterEntry &entry )
{
  return entry.getValue((T*)NULL);
}

/*! \relates ParameterEntry 
    Output stream operator for handling the printing of parameter entries.  
*/
inline ostream& operator<<(ostream& os, const ParameterEntry& e) 
{ 
  return e.leftshift(os);
}

} // namespace Teuchos


#endif
