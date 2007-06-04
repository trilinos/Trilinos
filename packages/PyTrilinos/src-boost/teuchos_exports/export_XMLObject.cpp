
// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

#include "Teuchos_XMLObject.hpp"

typedef Teuchos::XMLObject xml;

void expose_xmlobject()
{
    
    class_<xml>( "XMLObject", "An object representation of a subset of XML data" ,
            init<>("Empty constructor") )


    // XMLObject(const string& tag);
        .def( init<string>( "Construct using a node labeled by tag" ) )
    
    //  * \brief Construct with a pointer to the low-level representation. 
    //  *
    //  * This is used to allow construction of an XMLObject from the
    //  * XMLObjectImplem* return value of ExceptionBase::toXML().
    //  */
    // XMLObject(XMLObjectImplem* ptr);
    
    // XMLObject deepCopy() const ;
        .def("deepCopy", &xml::deepCopy,                
                "obj.deepCopy() -> XMLObject \n"
                "Make a deep copy of this object" )

    // const string& getTag() const
        .def( "getTag" ,&xml::getTag,
                "obj.getTag() -> str \n"
                "Return the tag of the current node",
                return_value_policy<copy_const_reference>())
    
    // bool hasAttribute(const string& name) const 
        .def( "hasAttribute" ,&xml::hasAttribute,
                "obj.hasAttribute() -> bool \n"
                "Find out if the current node has an attribute of the specified name")
    
    // const string& getAttribute(const string& name) const 
        .def( "getAttribute" ,&xml::getAttribute,
                "obj.getAttribute(name) -> str \n"
                "Return the value of the attribute with the specified name",
                return_value_policy<copy_const_reference>())
                
    // const string& getRequired(const string& name) const ;
        .def( "getRequired" ,&xml::getRequired,
                args("name"),
                "obj.getRequired(name) -> str \n"
                "Get an attribute, throwing an exception if it is not found",
                return_value_policy<copy_const_reference>()
        )
        
    // double getRequiredDouble(const string& name) const 
        .def( "getRequiredDouble" ,&xml::getRequiredDouble,
                args("name"),
                "obj.getRequiredDouble(name) -> float \n"
                "Get a required attribute, returning it as a double"
                )
    
    // int getRequiredInt(const string& name) const
        .def( "getRequiredInt" ,&xml::getRequiredInt,
                args("name"),
                "obj.getRequiredInt(name) -> int \n"
                "Get a required attribute, returning it as an int"
        )
    
    // bool getRequiredBool(const string& name) const ;
        .def( "getRequiredBool" ,&xml::getRequiredBool,
                args("name"),
                "obj.getRequiredBool(name) -> bool \n"
                "Get a required attribute, returning it as a bool"
        )
        
    // string getWithDefault(const string& name, 
    //                       const string& defaultValue) const ;
        .def( "getWithDefault" ,&xml::getWithDefault,
                args("name","defaultValue"),
                "obj.getWithDefault(name,defaultValue) -> str \n"
                "Get an attribute, assigning a default value if the requested \n"
                "XMLObject does not exists"
        )

    // int numChildren() const 
        .def( "numChildren" ,&xml::numChildren,
                "obj.numChildren( ) -> int \n"
                "Return the number of child nodes owned by this node"
        )
        
    // const XMLObject& getChild(int i) const
        .def( "getChild" ,&xml::getChild,
                args("i"),
                "obj.getChild( i ) -> XMLObject \n"
                "Return the i-th child node",
                return_internal_reference<>()
        )
        
    // int numContentLines() const 
        .def( "numContentLines" ,&xml::numContentLines,
                "obj.numContentLines( ) -> int \n"
                "Return the number of lines of character content stored in this node"
        )

    // const string& getContentLine(int i) const 
        .def( "getContentLine" ,&xml::getContentLine,
                "obj.getContentLine( i ) -> str \n"
                "Return the i-th line of character content stored in this node",
                return_value_policy<copy_const_reference>()
        )

    // string toString() const 
        .def( "toString" ,&xml::toString,
                "obj.toString( ) -> str \n"
                "Represent this node and its children as a string"
        )

    // string toString() const 
        .def( "__str__" ,&xml::toString,
                "obj.toString( ) -> str \n"
                "Represent this node and its children as a string"
        )

    // void print(ostream& os, int indent) const

    // string header() const 
        .def( "header" ,&xml::header,
                "obj.header( ) -> str \n"
                "Write the header for this object to a string"
        )

    // string terminatedHeader() const 
        .def( "terminatedHeader" ,&xml::terminatedHeader,
                "obj.terminatedHeader( ) -> str \n"
                "Write the header for this object to a string"
        )
        
    // string footer() const 
        .def( "footer" ,&xml::footer,
                "obj.footer( ) -> str \n"
                "Write the footer for this object to a string"
        )
        
    // Find out if a node is empty
    // bool isEmpty() const 
        .def( "isEmpty" ,&xml::isEmpty,
                "obj.isEmpty( ) -> bool \n"
                "Find out if a node is empty"
        )

    // void checkTag(const string& expected) const ;
        .def( "checkTag" ,&xml::checkTag,
                args("expected"),
                "obj.checkTag( expected ) -> None \n"
                "Check that a tag is equal to an expected string"
        )

    // void addAttribute(const string& name, const string& value)
        .def( "addAttribute" ,&xml::addAttribute,
                args("name","value"),
                "obj.addAttribute( name, value ) -> None \n"
                "Add an attribute to the current node's atribute list"
        )

    // void addDouble(const string& name, double val)
        .def( "addDouble" ,&xml::addDouble,
                args("name","val"),
                "obj.addDouble( name, val ) -> None \n"
                "Add a double as an attribute"
        )

    // void addInt(const string& name, int val)
        .def( "addInt" ,&xml::addInt,
                args("name","val"),
                "obj.addInt( name, val ) -> None \n"
                "Add an int as an attribute"
        )

    // void addBool(const string& name, bool val)
        .def( "addBool" ,&xml::addBool,
                args("name","val"),
                "obj.addBool( name, val ) -> None \n"
                "Add an int as an attribute"
        )

    // void addChild(const XMLObject& child)
        .def( "addChild" ,&xml::addChild,
                args("name","val"),
                "obj.addChild( child ) -> None \n"
                "Add a child node to the node"
        )

    // void addContent(const string& contentLine)
        .def( "addContent" ,&xml::addContent,
                args("contentLine"),
                "obj.addContent( contentLine ) -> None \n"
                "Add a line of character content"
        )
    ;

}