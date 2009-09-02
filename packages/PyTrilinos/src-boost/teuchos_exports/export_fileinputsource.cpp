

#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLInputSource.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;
#include "teuchos_call_policies.hpp"

class xisWrap : Teuchos::XMLInputSource, wrapper<Teuchos::XMLInputSource>
{
public:
    // virtual RefCountPtr<XMLInputStream> stream() const;
    Teuchos::RefCountPtr< Teuchos::XMLInputStream > stream( )
    {
        return this->get_override("stream")();
    }
};

void expose_fileinputsource()
{
    // class_< xisWrap >( "XMLInputSource" ,
    //     "XMLInputSource represents a source of XML input that can be parsed \n"
    //     "to produce an XMLObject.  \n"
    //     "Note: The source might be a file, a socket, a \n"
    //     "string. The XMLObject is created with a call to the getObject() method. \n"
    //     "The source gets its data from a XMLInputStream object that is \n"
    //     "created (internally) to work with this source. \n"
    //     " \n "
    //     "getObject() is implemented with EXPAT if Teuchos is configured with '--enable-expat' \n" )
    //     
    //     .def("stream", pure_virtual(&xisWrap::stream))
    //     
    //     // XMLObject getObject() const ;        
    //     .def("getObject", &xisWrap::getObject,
    //         "source.getObject() -> XMLObject \n"
    //         "Get an object by invoking the TreeBuildingXMLHandler on the \n"
    //         "input data \n"
    //     )
    // ;

      // FileInputSource(const string& filename);
    class_<Teuchos::FileInputSource /*, bases< Teuchos::XMLInputSource >*/ >("FileInputSource",
            "Definition of XMLInputSource derived class for reading XML from a file",
            init<string>() )
        // virtual RefCountPtr<XMLInputStream> stream() const;
        .def( "stream", &Teuchos::FileInputSource::stream ,
                        del_dangleing_rcp<>(),
                        "Create a FileInputStream"
        )
        .def("getObject", &Teuchos::XMLInputSource::getObject,
            "source.getObject() -> XMLObject \n"
            "Get an object by invoking the TreeBuildingXMLHandler on the \n"
            "input data \n"
        )
        

    ;
}
