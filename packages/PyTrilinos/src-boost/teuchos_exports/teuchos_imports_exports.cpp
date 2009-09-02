
#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FILEstream.hpp"
#include "python_ParameterEntry.hpp"
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
        object* boostobj = new object(handle<>( o ) );
        obj *counted = extract<obj*>( *boostobj );
        rcpobj *counter = new rcpobj( counted , false );
        Teuchos::set_extra_data( *boostobj, "boostpyobj", Teuchos::inOutArg(counter));
        return counter;
    }
    static void from_python()
    {
        converter::registry::insert( &extract_rcp , type_id< rcpobj >() );
    }
    
    struct rcp_to_object
    {
        static PyObject* convert(rcpobj const& x)
        {
            object *myObj = NULL;
            myObj = Teuchos::get_optional_extra_data< object >( const_cast< rcpobj& >(x),"boostpyobj" );
            if ( myObj == NULL)
            {
        	    myObj = new object( x.get() );
        	    Teuchos::set_extra_data( *myObj, "boostpyobj",
                      Teuchos::inOutArg(const_cast< rcpobj& >(x)) );
        	}
        	else
        	{
        	    Py_XINCREF( myObj->ptr() );
        	}
            return myObj->ptr();
        }
    };

    static void to_python()
    {
        to_python_converter< rcpobj, rcp_to_object>();
    }
};

// ################################################################################
// ## Property List
// ################################################################################

void* extract_plist(PyObject * o)
{
     if ( !PyDict_Check( o ) )
         return NULL;
     else
        return (void *) python_plist_tools::pyDictToNewParameterList(o,raiseError);
}


struct plist_to_object
{
    static PyObject* convert(const Teuchos::ParameterList & plist)
    {
        return python_plist_tools::parameterListToNewPyDict( plist,raiseError);
    }
};

// ################################################################################
// ## IO
// ################################################################################

void* extract_ostream(PyObject * pf)
{
    if ( !PyFile_Check( pf ) )
    {
        PyErr_SetString( PyExc_TypeError, "must be a file object" );
        throw_error_already_set();
        return NULL;
    }
    else
    {
        std::FILE *f = PyFile_AsFile(pf);
    	Teuchos::FILEstream buffer(f);
    	ostream *os = new ostream(&buffer);
    	return os;
    }
    
}

struct ostream_to_python
{
    static PyObject* convert(const ostream & plist)
    {
        return NULL;

    }
};


// ################################################################################
// ## Extracts and Inserts
// ################################################################################
void extract_teuchos_misc()
{
    py_rcp< Teuchos::XMLInputStream >::from_python();
    py_rcp< Teuchos::XMLInputStream >::to_python();
    
    // register_ptr_to_python< shared_ptr< Teuchos::XMLInputStream > >();
    
    converter::registry::insert( &extract_plist , type_id< Teuchos::ParameterList >() );
    to_python_converter< Teuchos::ParameterList, plist_to_object >();
    
    converter::registry::insert( &extract_ostream , type_id< std::ostream >() );
    //to_python_converter< std::ostream , ostream_to_python >();
}
// ################################################################################
// ##
// ################################################################################
