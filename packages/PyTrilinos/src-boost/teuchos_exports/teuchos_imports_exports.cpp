
#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_ParameterList.hpp"
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


void add_something_to_plist(Teuchos::ParameterList& plist,object o,string name)
{
    extract<string>            x1(o);
    extract<int>               x2(o);
    extract<double>            x3(o);
    extract<complex< double> > x4(o);
    
    if ( x1.check() )
        plist.set(name, x1() );
    else if ( x2.check() )
        plist.set(name, x2());
    else if ( x3.check() )
        plist.set(name, x3());
    else if ( x4.check() )
        plist.set(name, x4());
        
    return;
    
}

void* extract_plist(PyObject* o)
{   
    object keys,key,something;
    string name;
    dict pydict = dict(handle<>(borrowed( o )) );
    
    Teuchos::ParameterList * plist = new Teuchos::ParameterList();
    
    int len = extract<int>( pydict.attr("__len__")() );
    
    keys = pydict.keys();
    
    for(int i=0; i < len; i++ )
    {
        key = keys[i];
        something = pydict[ key ];
        name = extract<string>(key);
        add_something_to_plist(*plist,something,name);
    }
    
    return plist;
}




// ################################################################################
// ## Extracts and Inserts
// ################################################################################
void extract_teuchos_misc()
{
    py_rcp< Teuchos::XMLInputStream >::from_python();
    py_rcp< Teuchos::XMLInputStream >::to_python();
    
    converter::registry::insert( &extract_plist , type_id< Teuchos::ParameterList >() );
}
// ################################################################################
// ##
// ################################################################################
