
#include "Teuchos_XMLInputSource.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

// ################################################################################
// ## General functions and classes
// ################################################################################

template< class obj >
class py_rcp
{
    typedef Teuchos::RefCountPtr<obj> rcpobj;
public:
    static void* extract_rcp(PyObject* o)
    {   
        object boostobj = object(handle<>(borrowed( o )) );
        obj *counted = extract<obj*>( boostobj );
        rcpobj *myObj = new rcpobj( counted , true );
        return myObj;
    }
    static void from_python()
    {
        converter::registry::insert( &extract_rcp , type_id< rcpobj >() );
    }
    
    struct rcp_to_object
    {
        static PyObject* convert(rcpobj const& x)
        {
        	object myObj = object( x.get() );
        	Py_XINCREF( myObj.ptr() );
            return myObj.ptr();
        }
    };

    static void to_python()
    {
        to_python_converter< rcpobj, rcp_to_object>();
    }
};

// ################################################################################
// ## Extracts and Inserts
// ################################################################################
void extract_teuchos_misc()
{
    py_rcp< Teuchos::XMLInputStream >::from_python();
    py_rcp< Teuchos::XMLInputStream >::to_python();
}
// ################################################################################
// ##
// ################################################################################
