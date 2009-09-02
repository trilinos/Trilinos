
#include <Teuchos_RefCountPtr.hpp>
// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;


template<class BasePolicy_ = default_call_policies>
struct del_dangleing_rcp : BasePolicy_
{
    void * my_rcp;
    static bool precall(PyObject*)
    {
        return true;
    }
    
    template <class ArgumentPackage>
    static PyObject* postcall(ArgumentPackage const& args_, PyObject* result)
    {
        result = BasePolicy_::postcall(args_, result);
        return result;
    }
};
