
#ifndef TEUCHOS_PARAMETER_ENTRY_H
#define TEUCHOS_PARAMETER_ENTRY_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_any.hpp"

namespace Teuchos {

//! %Teuchos %Parameter support.

class ParameterList; // anonther parameter type (forward declaration)

//! Manipulating single parameters, including sublists.
class ParameterEntry {

public:
  
  //! Default Constructor
  ParameterEntry();
  
  //! Copy constructor
  ParameterEntry(const ParameterEntry& source);

  //! Copy
  ParameterEntry& operator=(const ParameterEntry& source);

  //! Templated constructor
  template<typename T>
  ParameterEntry(T value, bool isDefault = false)
  : val_(value),
    isUsed_(false),
    isDefault_(isDefault)
  {}

  //! Destructor
  ~ParameterEntry() {};

  /** @name ParameterList parameters
   *
   * Function for creating a parameter that is an empty list.  */
  //@{
  ParameterList& setList(bool isDefault = false);
  //@}

  /** @name Set functions. 
   
    The input value type determines the type of parameter
    stored. Invalidates any previous values stored by this object,
    although it doesn't necessarily erase them. Resets 'isused'
    functionality. 
  */
  //@{ 
  template<typename T>
  void setValue(T value, bool isDefault = false)
  {
    val_ = value;
    isDefault_ = isDefault;
  }
  //@}

  //! Return whether or not the value is used; i.e., whether or not
  //! the value has been retrieved via a get function.
  bool isList() const { return isList_; };
  bool isUsed() const { return isUsed_; };

  /** @name Get functions. 
   
    Returns value of parameter. The value is nonsense if we do not
    request the correct type of value. We cannot name all of these
    functions the same since the language does not allow us to
    overload functions based solely on return value. 
  */
  //@{

  template<typename T>
  T& getValue(T *ptr) const
  {
    isUsed_ = true;
    return const_cast<T&>(Teuchos::any_cast<T>( val_ ));
  }

  //@}

  //! Output the parameter to the given os. 
  /*! 
    Formats the output as "<type,value>", except in the case of a
    list which just outputs "\<sublist\>". If the parameter has not yet
    been set, it outputs "\<NONE\>". This is the function called by the
    ostream operator<<. 
  */
  ostream& leftshift(ostream& os) const;

private:

  //! Reset the entry
  void reset();
  
  //! Data
  any val_;

  //! List flag
  bool isList_;

  //! Has this parameter been accessed by a "get" function?
  mutable bool isUsed_;

  //! Was this parameter a default value assigned by a "get" function?
  mutable bool isDefault_;

};

template<typename T>
T& getValue( const ParameterEntry &entry )
{
  return entry.getValue((T*)NULL);
}

//! Output the parameter. Relies of leftshift operator defined in the class.
inline ostream& operator<<(ostream& os, const ParameterEntry& e) 
{ 
  return e.leftshift(os);
}

} // namespace Teuchos


#endif
