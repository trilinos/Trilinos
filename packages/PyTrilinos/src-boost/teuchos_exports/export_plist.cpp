
// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

#include <Teuchos_ParameterList.hpp>


typedef Teuchos::ParameterList plist;
void expose_plist()
{
    
    //   ParameterList();
    //   ParameterList(const std::string &name);
    //   ParameterList(const ParameterList& source);
    //void ( matrix::*fillComplete1 )( void ) = &matrix::fillComplete;
    
    class_<plist>( "ParameterList",
    
        "Parameters can be added and retreived with the templated \"get\" and \"set\"   \n"
        "functions.  These parameters can any data type which uses value sementics      \n"
        "(e.g. double, float, int, *double, *float, *int, ...) which includes other     \n"
        "parameter lists, allowing for a hierarchy of parameter lists.  These           \n"
        "parameters can also be pointers to vectors or functions.                       \n" , 
        init<>()
        )


    //  void setName( const std::string &name );
        .def("setName", &plist::setName, args("name"),
                "Set the name of this list.")

    //   ParameterList& setParameters(const ParameterList& source);
        .def("setParameters", &plist::setParameters, args("source"),
                
                "Set the parameters in source.  \n"
                "\n"
                "Note, this function will set the parameters and sublists from \n"
                "source into this but will not result in parameters \n"
                "being removed from this.  Parameters in this with the \n"
                "same names as those in source will be overwritten. \n",
                return_internal_reference<>()    )


    //   ParameterList& setParametersNotAlreadySet(const ParameterList& source);
        .def("setParametersNotAlreadySet", &plist::setParametersNotAlreadySet, args("source"),
             
             "Set the parameters in source that are not already set in \n"
             "this. \n"
             "\n"
             "Note, this function will set the parameters and sublists from \n"
             "source into this but will not result in parameters \n"
             "being removed from this or in parameters already set in \n"
             "this being overrided.  Parameters in this with the\n"
             "same names as those in source will not be overwritten. \n" ,
             return_internal_reference<>() )
    
    //   /*!  Sets different types of parameters. The type depends on the second entry.  
    //     
    //     \note <ul>
    //     <li> Use static_cast<T>() when the type is ambiguous. 
    //     <li> Both char* and string map to are stored as strings internally. 
    //     <li> Sets the parameter as "unused".
    //     </ul>
    //   */
    //   template<typename T>
    //   void set(
    //     std::string const& name, T const& value, std::string const& docString = ""
    //     ,RefCountPtr<const ParameterEntryValidator> const& validator = null
    //     );
    
    // 
    //   /*!  Template specialization for the case when a user sets the parameter with a character
    //     string in parenthesis.
    //   */
    //   void set(
    //     std::string const& name, char value[], std::string const& docString = ""
    //     ,RefCountPtr<const ParameterEntryValidator> const& validator = null
    //     );
    // 
    //   /*!  Template specialization for the case when a user sets the parameter with a character
    //     string in parenthesis.
    //   */
    //   void set(
    //     std::string const& name, const char value[], string const& docString = ""
    //     ,RefCountPtr<const ParameterEntryValidator> const& validator = null
    //     );
    // 
    //   /*!  Template specialization for the case when a user sets the parameter with a ParameterList.
    // /
    //   void set(
    //     std::string const& name, ParameterList const& value, string const& docString = ""
    //     );
    // 
    //   /*!  Set a parameter directly as a ParameterEntry. 
    //  \note This is required to preserve the isDefault value when reading back
    //  from XML. KL 7 August 2004 
    // /
    //   void setEntry(const string& name, const ParameterEntry& entry);
    // 
    //   
    //   
    //   @name Get Functions 
    //   
    //   
    //   /*!  Retrieves parameter \c name of type \c T from list, if it exists, else the \c def_value is
    //     used to enter a new parameter into the list.
    //     
    //     \note <ul> 
    //     <li> Use the static_cast<T>() when the type is ambiguous.
    //     <li> Both char* and string map to are stored as strings internally. 
    //     <li> Sets the parameter as "used".
    //     <li> Exception is thrown if \c name exists, but is not of type \c T.
    //     </ul>
    //   */
    //   template<typename T>
    //   T& get(const string& name, T def_value);
    // 
    //   /*!  Template specialization of get, where the nominal value is a character string in parenthesis.
    //     Both char* and string are stored as strings and return string values.
    //   */
    //   std::string& get(const string& name, char def_value[]);
    //   
    //   /*!  Template specialization of get, where the nominal value is a character string in parenthesis.
    //     Both char* and string are stored as strings and return string values.
    //   */
    //   std::string& get(const string& name, const char def_value[]);
    //   
    //   /*!  Retrieves parameter \c name of type \c T from a list, an
    //     Exceptions::InvalidParameter exception is thrown if this
    //     parameter doesn't exist (Exceptions::InvalidParameterName) or is
    //     the wrong type (Exceptions::InvalidParameterName).  \note The
    //     syntax for calling this method is:  list.template get<int>( "Iters" )
    //     
    //   */
    //   template<typename T>
    //   T& get(const string& name);
    //   
    //   /*!  Retrieves parameter \c name of type \c T from a constant list, an
    //     Exceptions::InvalidParameter exception is thrown if this
    //     parameter doesn't exist (Exceptions::InvalidParameterName) or is
    //     the wrong type (Exceptions::InvalidParameterName).  \note The
    //     syntax for calling this method is:  list.template get<int>( "Iters" )
    //     
    //   */
    //   template<typename T>
    //   const T& get(const string& name) const;  
    //   
    //   /*!  Retrieves the pointer for parameter \c name of type \c T from a
    //     list.  A null pointer is returned if this parameter doesn't exist or is
    //     the wrong type.  \note The syntax for calling this method is: 
    //     list.template getPtr<int>( "Iters" ) 
    //   */
    //   template<typename T>
    //   T* getPtr(const string& name);
    //   
    //   /*!  Retrieves the pointer for parameter \c name of type \c T from a
    //     constant list.  A null pointer is returned if this parameter doesn't exist
    //     or is the wrong type.  \note The syntax for calling this method is: 
    //     list.template getPtr<int>( "Iters" ) 
    //   */
    //   template<typename T>
    //   const T* getPtr(const string& name) const;  
    //   
    //   /*!  Retrieves an entry with the name name.
    // 
    //  Throws Exceptions::InvalidParameterName if this parameter does
    //  not exist.
    // /
    //   ParameterEntry& getEntry(const string& name);  
    //   
    //   /*!  Retrieves a const entry with the name name.
    // 
    //  Throws Exceptions::InvalidParameterName if this parameter does
    //  not exist.
    // /
    //   const ParameterEntry& getEntry(const string& name) const;  
    //   
    //   /*!  Retrieves the pointer for an entry with the name name if
    //   it exists. */
    //   ParameterEntry* getEntryPtr(const string& name);  
    //   
    //   /*!  Retrieves the pointer for a constant entry with the name name if
    //   it exists. */
    //   const ParameterEntry* getEntryPtr(const string& name) const;  
    // 
    //   
    // 
    //   @name Parameter removal functions
    //   
    //  
    //  Remove a parameter (does not depend on the type).  
    // 
    //  \param  name
    //            [in] The name of the parameter to remove
    //  \param  throwIfNotExists
    //            [in] If true then if the parameter with
    //            the name name an exception will be thrown!
    // 
    //  Returns true of the parameter was removed, and false if
    //  the parameter was not removed (return value possible only if
    //  throwIfExists==false).
    // /
    //   bool remove(
    //     std::string const& name, bool throwIfNotExists = true
    //     );
    // 
    //   
    //   
    //   @name Sublist Functions 
    //   
    // 
    //   /*!  Creates an empty sublist and returns a reference to the sublist
    //   \c name. If the list already exists, returns reference to that
    //   sublist. If the name exists but is not a sublist, an exception is
    //   thrown.
    // /
    //   ParameterList& sublist(
    //     const string& name, bool mustAlreadyExist = false
    //     ,const string& docString = ""
    //     );
    //   
    //   /*!  Return a const reference to an existing sublist \c name.  If the
    //   list does not already exist or the name exists but is not a sublist, an
    //   exception is thrown.
    // /
    //   const ParameterList& sublist(const string& name) const;
    // 
    //   
    //   
    //   @name Attribute Functions 
    //   
    // 
    //   /*!  Query the name of this parameter list. */
    //   const std::string& name() const;
    // 
    //   /*!  Query the existence of a parameter.  \return "true" if a
    //     parameter with this \c name exists, else "false".
    //   */
    //   bool isParameter(const string& name) const;
    //   
    //   /*!  Query the existence of a parameter and whether it is a parameter
    //     list.  \return "true" if a parameter with this \c name exists and is
    //     itself a parameter list, else "false".
    //   */
    //   bool isSublist(const string& name) const;
    //   
    //   /*!  Query the existence and type of a parameter.  \return "true" is a
    //     parameter with this \c name exists and is of type \c T, else "false".
    //     \note The syntax for calling this method is:  list.template
    //     isType<int>( "Iters" ) 
    //   */
    //   template<typename T>
    //   bool isType(const string& name) const;
    // 
    // #ifndef DOXYGEN_SHOULD_SKIP_THIS  
    //   /*!  Query the existence and type of a parameter.
    //     \return "true" is a parameter with this \c name exists and is of type \c T, else "false".
    //     \note <b>It is not recommended that this method be used directly!</b>  
    //     Please use either the helper function <b>isParameterType</b> or non-nominal <b>isType</b> method. 
    //   */
    //   template<typename T>
    //   bool isType(const string& name, T* ptr) const;
    // #endif
    // 
    //   
    //   
    //   @name I/O Functions 
    //   
    // 
    //   /*!  Printing method for parameter lists which takes an print options
    //   object.*/
    //   ostream& print(ostream& os, const PrintOptions &printOptions ) const;
    // 
    //   /*!  Printing method for parameter lists.  Indenting is used to indicate
    //     parameter list hierarchies. */
    //   ostream& print(ostream& os, int indent = 0, bool showTypes = false, bool showFlags = true ) const;
    //   
    //   Print out unused parameters in the ParameterList.
    //   void unused(ostream& os) const;
    // 
    //   Create a single formated string of all of the zero-level parameters in this list
    //   std::string currentParametersString() const;
    // 
    //   
    // 
    //   @name Read-only access to the iterator 
    //   
    // 
    //   An iterator pointing to the first entry
    //   ConstIterator begin() const ;
    // 
    //   An iterator pointing beyond the last entry
    //   ConstIterator end() const ;
    // 
    //   Access to ParameterEntry (i.e., returns i->second)
    //   const ParameterEntry& entry(ConstIterator i) const;
    //   
    //   Access to name (i.e., returns i->first)
    //   const string& name(ConstIterator i) const;
    // 
    //   
    // 
    //   @name Validation Functions 
    //   
    // 
    //  Validate the parameters is this list given valid selections in
    //  the input list.
    // 
    //  \param  validParamList
    //               [in] This is the list that the parameters and sublist in this
    //               are compared against.
    //  \param  depth
    //               [in] Determines the number of levels of depth that the validation will
    //               recurse into.  A value of dpeth=0 means that only the top level
    //               parameters and sublists will be checked.  Default: depth = large number.
    //  \param  validateUsed
    //               [in] Determines if parameters that have been used are checked
    //               against those in validParamList.  Default:
    //               validateDefaults = VALIDATE_DEFAULTS_ENABLED.
    //  \param  validateDefaults
    //               [in] Determines if parameters set at their default values using get(name,defaultVal)
    //               are checked against those in validParamList.  Default:
    //               validateDefaults = VALIDATE_DEFAULTS_ENABLED.
    // 
    //  If a parameter in this is not found in validParamList
    //  then an exception of type Exceptions::InvalidParameterName will
    //  be thrown which will contain an excellent error message returned by
    //  excpt.what().  If the parameter exists but has the wrong type,
    //  then an exception type Exceptions::InvalidParameterType will be
    //  thrown.  If the parameter exists and has the right type, but the value is
    //  not valid then an exception type
    //  Exceptions::InvalidParameterValue will be thrown.
    // 
    //  A breath-first search is performed to validate all of the parameters in
    //  one sublist before moving into nested subslist.
    // /
    //   void validateParameters(
    //     ParameterList const& validParamList,
    //     int const depth = 1000,
    //     EValidateUsed const validateUsed = VALIDATE_USED_ENABLED,
    //     EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED
    //     ) const;
    // 
    //  Validate the parameters is this list given valid selections in
    //  the input list and set defaults for those not set.
    // 
    //  \param  validParamList
    //               [in] This is the list that the parameters and sublist in this
    //               are compared against.
    //  \param  depth
    //               [in] Determines the number of levels of depth that the validation will
    //               recurse into.  A value of dpeth=0 means that only the top level
    //               parameters and sublists will be checked.  Default: depth = large number.
    // 
    //  If a parameter in this is not found in validParamList
    //  then an exception of type Exceptions::InvalidParameterName will
    //  be thrown which will contain an excellent error message returned by
    //  excpt.what().  If the parameter exists but has the wrong type,
    //  then an exception type Exceptions::InvalidParameterType will be
    //  thrown.  If the parameter exists and has the right type, but the value is
    //  not valid then an exception type
    //  Exceptions::InvalidParameterValue will be thrown.  If a
    //  parameter in validParamList does not exist in this,
    //  then it will be set at its default value as determined by
    //  validParamList.
    // 
    //  A breath-first search is performed to validate all of the parameters in
    //  one sublist before moving into nested subslist.
    // /
    //   void validateParametersAndSetDefaults(
    //     ParameterList const& validParamList,
    //     int const depth = 1000
    //     );
    ;
}