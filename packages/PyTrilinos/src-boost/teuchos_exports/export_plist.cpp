
// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

//#include <Teuchos_ParameterList.hpp>
#include "plist_inject.hpp"
#include "python_ParameterEntry.hpp"

typedef Teuchos::ParameterList plist;

string get_name( plist const& self ){ return self.name(); }

// ####################################################################################
// ## Added Class Functions
// ####################################################################################

// PyObject* asDict( plist& self )
// { return python_plist_tools::parameterListToNewPyDict( self, raiseError); }

// object items( plist& self )
// { return object(); }


//####################################################################################    //####################################################################################
void expose_plist()
{
    
    enum_<Teuchos::EValidateUsed>("EValidateUsed")
      .value("VALIDATE_USED_ENABLED",
            Teuchos::VALIDATE_USED_ENABLED )
      .value("VALIDATE_USED_DISABLED",
            Teuchos::VALIDATE_USED_DISABLED )
    ;

    enum_<Teuchos::EValidateDefaults>("EValidateDefaults")
        .value( "VALIDATE_DEFAULTS_ENABLED",
            Teuchos::VALIDATE_DEFAULTS_ENABLED )
        .value( "VALIDATE_DEFAULTS_DISABLED",
            Teuchos::VALIDATE_DEFAULTS_DISABLED )
    ;
    
    
    Teuchos::ParameterEntry& (plist::*getEntry1)(const string&) = &plist::getEntry;
    plist& (plist::*sublist1)(const string&, bool,const string&) = &plist::sublist;
    // ostream& (plist::*print1)(ostream&, int, bool, bool ) const = &plist::print;
    
    //   ParameterList();
    class_<plist>( "ParameterList",
    
        "Parameters can be added and retreived with the templated \"get\" and \"set\"   \n"
        "functions.  These parameters can any data type which uses value sementics      \n"
        "(e.g. double, float, int, ...) which includes other                            \n"
        "parameter lists, allowing for a hierarchy of parameter lists.  These           \n"
        "parameters can also be pointers to vectors or functions.                       \n" , 
        init<>()
        )
    //   ParameterList(const ParameterList& source);            
        .def( init<plist>() )

    //   ParameterList(const std::string &name);
        .def( init<string>() )
        
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
             return_internal_reference<>() 
        )
        
        
        .def("set", &python_plist_tools::setPythonParameter,
                args("self","name","value"),
                "plist.get(name,value) -> bool " )
                
        .def("get", &python_plist_tools::getPythonParameter,
                args("self","name"),
                "plist.get(name) -> object " )

        // ParameterEntry& getEntry(const string& name);
        .def("getEntry", getEntry1,
                args("name"),
                
                "plist.getEntry( name ) -> ParameterEntry \n\n" 
                "Retrieves an entry with the name name.  \n"
                "       \n "
                "Throws Exceptions::InvalidParameterName if this parameter does \n"
                "not exist. \n",
                
                return_internal_reference<>()
                )
                
    //   bool remove(
    //     std::string const& name, bool throwIfNotExists = true
    //     );
        .def("remove",&plist::remove, ( args("name"),args("throwIfNotExists")=true ),
                "plist.remove( name, throwIfNotExists=True ) -> bool \n\n"
                "Returns true of the parameter was removed, and false if \n"
                "the parameter was not removed (return value possible only if \n"
                "throwIfExists is False). \n"
        )
        
    //   ParameterList& sublist(
    //     const string& name, bool mustAlreadyExist = false
    //     ,const string& docString = ""
    //     );
        .def("sublist", sublist1 ,
                ( args("name"),args("mustAlreadyExist")=false,args("docString")="" ),
                " plist.sublist( name ,mustAlreadyExist=False, docString="" ) -> plist \n"
                " \n"
                "reates an empty sublist and returns a reference to the sublist \n"
                "name. If the list already exists, returns reference to that \n"
                "sublist. If the name exists but is not a sublist, an exception is \n"
                "thrown.                 \n",
                return_internal_reference<>()
        )

    //   const std::string& name() const;
        .def("name", &get_name ,
                " plist.name( ) -> plist \n"
                )
    

    //   /*!  Query the existence of a parameter.  \return "true" if a
    //     parameter with this \c name exists, else "false".
    //   */
    //   bool isParameter(const string& name) const;
        .def("isParameter",&plist::isParameter,
                "plist.isParameter() -> bool \n\n"
                "Query the existence of a parameter return \"true\" if a\n"
                "parameter with this name exists, else \"false\"."
        )


    //   bool isSublist(const string& name) const;
        .def("isSublist",&plist::isSublist,
                "plist.isSublist() -> bool \n\n"
                "Query the existence of a parameter return \"true\" if a\n"
                "parameter with this name exists, else \"false\"."
        )    
        
    //   /*!  Printing method for parameter lists which takes an print options
    //   object.*/
    //   ostream& print(ostream& os, const PrintOptions &printOptions ) const;
        
    //   ostream& print(ostream& os, int indent = 0, bool showTypes = false, bool showFlags = true ) const;
        .def("_print", &_myprint ,
                ( args("pf")=NULL, args("indent")=0, args("showTypes")=false,
        	      args("showFlags")=true  ),
        
                "plist._print(os, indent=0, showTypes=False, showFlags=True ) -> os \n\n"
                "Printing method for parameter lists.  Indenting is used to indicate \n"
                "parameter list hierarchies. \n"
        )    

    //   void unused(ostream& os) const;
        // .def("unused", &plist::unused ,
        //         ( args("pf")=NULL, args("indent")=0, args("showTypes")=false,
        //         	      args("showFlags")=true  ),
        //         
        //         "plist.unused(os, indent=0, showTypes=False, showFlags=True ) -> os \n\n"
        //         "Print out unused parameters in the ParameterList \n"
        // 
        // )    

    //   std::string currentParametersString() const;
        .def("currentParametersString",&plist::currentParametersString,
                "Create a single formated string of all of the zero-level parameters in this list")

    //   void validateParameters(
    //     ParameterList const& validParamList,
    //     int const depth = 1000,
    //     EValidateUsed const validateUsed = VALIDATE_USED_ENABLED,
    //     EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED
    //     ) const;
        .def("validateParameters",&plist::validateParameters,
                ( args("validParamList"),
                  args("depth")=1000,
                  args("validateUsed")=Teuchos::VALIDATE_USED_ENABLED,
                  args("validateDefaults")=Teuchos::VALIDATE_DEFAULTS_ENABLED ),

                "plist.validateParameters(validParamList, \n"
                "   depth=1000, \n"
                "   validateUsed=VALIDATE_USED_ENABLED, \n"
                "   validateDefaults=VALIDATE_DEFAULTS_ENABLED) \n\n"
                "  \n"
                "Validate the parameters is this list given valid selections in                     \n"
                "the input list.                                                                    \n"
                "                       \n"
                "validParamList                                                                     \n"
                "          [in] This is the list that the parameters and sublist in this            \n"
                "          are compared against.                                                    \n"
                "depth                  \n"
                "          [in] Determines the number of levels of depth that the validation will   \n"
                "          recurse into.  A value of dpeth=0 means that only the top level          \n"
                "          parameters and sublists will be checked.  Default: depth = large number. \n"
                "validateUsed           \n"
                "          [in] Determines if parameters that have been used are checked            \n"
                "          against those in validParamList.  Default:                               \n"
                "          validateDefaults = VALIDATE_DEFAULTS_ENABLED.                            \n"
                "validateDefaults       \n"
                "          [in] Determines if parameters set at their default values                \n"
                "          using get(name,defaultVal)                                               \n"
                "          are checked against those in validParamList.  Default:                   \n"
                "          validateDefaults = VALIDATE_DEFAULTS_ENABLED.                            \n"
                "                       \n"
                "If a parameter in this is not found in validParamList                              \n"
                "then an exception of type Exceptions::InvalidParameterName will                    \n"
                "be thrown which will contain an excellent error message returned by                \n"
                "excpt.what().  If the parameter exists but has the wrong type,                     \n"
                "then an exception type Exceptions::InvalidParameterType will be                    \n"
                "thrown.  If the parameter exists and has the right type, but the value is          \n"
                "not valid then an exception type                                                   \n"
                "Exceptions::InvalidParameterValue will be thrown.                                  \n"
                "                       \n"
                "A breath-first search is performed to validate all of the parameters in            \n"
                "one sublist before moving into nested subslist.                                    \n"
    
        )
      // void validateParametersAndSetDefaults(
      //   ParameterList const& validParamList,
      //   int const depth = 1000
      //   );
        .def("validateParametersAndSetDefaults",&plist::validateParametersAndSetDefaults,
                ( args("validParamList"),
                  args("depth")=1000  ),

                 "Validate the parameters is this list given valid selections in    \n"
                 "the input list and set defaults for those not set.                \n"
                 "    \n"
                 "validParamList    \n"
                 "    [in] This is the list that the parameters and sublist in this    \n"
                 "    are compared against.                                            \n"
                 "depth    \n"
                 "    [in] Determines the number of levels of depth that the validation will    \n"
                 "    recurse into.  A value of dpeth=0 means that only the top level           \n"
                 "    parameters and sublists will be checked.  Default: depth = large number.  \n"
                 "    \n"
                 "If a parameter in this is not found in validParamList                 \n"
                 "then an exception of type Exceptions::InvalidParameterName will       \n"
                 "be thrown which will contain an excellent error message returned by   \n"
                 "excpt.what().  If the parameter exists but has the wrong type,        \n"
                 "then an exception type Exceptions::InvalidParameterType will be       \n"
                 "thrown.  If the parameter exists and has the right type, but the value is    \n"
                 "not valid then an exception type                              \n"
                 "Exceptions::InvalidParameterValue will be thrown.  If a       \n"
                 "parameter in validParamList does not exist in this,           \n"
                 "vthen it will be set at its default value as determined by    \n"
                 "validParamList.    \n"
                 "    \n"
                 "A breath-first search is performed to validate all of the parameters in    \n"
                 "one sublist before moving into nested subslist.    \n"
        )
    
    // ####################################################################################
    // ## Injected Methods
    // ####################################################################################

        //////////////////////////////////////////////
        // ** The following methods are added to ** //
        // ** give ParameterList a PyDict "feel" ** //
        //////////////////////////////////////////////
        .def("__setitem__", &python_plist_tools::setPythonParameter)
        .def("__getitem__", &python_plist_tools::getPythonParameter)
        
        .def("type",&type)
        .def("__cmp__",&__cmp__)
        .def("__contains__",&__contains__)
        .def("__eq__",&__eq__)
        .def("__iter__",&__iter__)
        .def("__len__",&__len__)
        .def("__ne__",&__ne__)
        .def("__repr__",&__repr__)
        .def("__str__",&__str__)
        .def("has_key",&has_key)
        .def("items",&items)
        .def("iteritems",&iteritems)
        .def("iterkeys",&iterkeys)
        .def("itervalues",&itervalues)
        .def("keys",&keys)
        .def("update",&update)
        .def("values",&values)
        .def("asDict",&asDict)
    ;
}