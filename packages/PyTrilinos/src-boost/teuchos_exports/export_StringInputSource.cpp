

#include "Teuchos_StringInputSource.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;


void expose_str_inputsource()
{
    // StringInputSource(const string& text);
    class_< Teuchos::StringInputSource > ( "StringInputSource",
            "Definition of XMLInputSource derived class for reading XML from \n"
        	"a string    \n",
        	init< string >() )
        	

    // virtual RefCountPtr<XMLInputStream> stream() const;
        .def( "stream" , &Teuchos::StringInputSource::stream,
                "Create a StringInputStream" 
        )
        .def("getObject", &Teuchos::XMLInputSource::getObject,
            "source.getObject() -> XMLObject \n"
            "Get an object by invoking the TreeBuildingXMLHandler on the \n"
            "input data \n"
        )
    ;
}

