
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_FileInputSource.hpp"

typedef Teuchos::XMLParameterListWriter writer;
typedef Teuchos::XMLParameterListReader reader;

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

// There is somthing wrong with the method ??
Teuchos::XMLObject my_toXML(writer const& self, const Teuchos::ParameterList& p)
{
    return self.toXML(p);
}

Teuchos::ParameterList my_toplist(reader const& self, const Teuchos::XMLObject& xml)
{
    return self.toParameterList(xml);
}



void expose_xml_r_w()
{
    
    class_<reader>( "XMLParameterListReader" ,
             "Writes an XML object to a parameter list \n",
             init<>("Construct a reader") )

        // ParameterList toParameterList(const XMLObject& xml) const ;
        .def( "toParameterList", &my_toplist,
                    "reader.toParameterList( xml ) -> plist  \n"
                    "Write the given XML object to a parameter list"
        )
    ;
    
    //########    //########
    //########    //########
    
    class_<writer>( "XMLParameterListWriter" ,
             "Writes a ParameterList to an XML object \n",
             init<>("Construct a writer") )

 		// XMLObject toXML(const ParameterList& p) const ;
        .def( "toXML", &my_toXML,
                    "writer.toXML( p ) -> XMLObject \n"
                    "Write the given list to an XML object"
        )
    ;    
}