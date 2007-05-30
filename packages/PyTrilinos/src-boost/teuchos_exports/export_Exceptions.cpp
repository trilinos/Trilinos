
#include <boost/python.hpp>
using namespace boost::python;

#include "Teuchos_XMLObject.hpp"

void EmptyXMLError_translator(Teuchos::EmptyXMLError const& x) 
{
    PyErr_SetString( PyExc_RuntimeError, x.what() );
}


void convert_exceptions()
{
    register_exception_translator< Teuchos::EmptyXMLError >(EmptyXMLError_translator);
}
