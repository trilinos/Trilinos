// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER


// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

// Teuchos initialization
#include "Teuchos_Version.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputSource.hpp"
// ################################################################################
// ## General functions and classes
// ################################################################################

template<class BasePolicy_ = default_call_policies>
struct del_dangleing_rcp : BasePolicy_
{
    void * my_rcp;
    template <class ArgumentPackage>
    static PyObject* postcall(ArgumentPackage const& args_, PyObject* result)
    {
        result = BasePolicy_::postcall(args_, result);
        return result;
    }
};

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

void expose_time()
{
  // Teuchos version support
  def("Teuchos_Version", Teuchos::Teuchos_Version);

  // Teuchos Time support
  class_<Teuchos::Time>("Time", "Time(name,start=False)\n"
			"Basic wall-clock timer class.",
			init<std::string, bool>( ( args("name" ),
						   args("start")=false ) )
		       )
    .def("start", &Teuchos::Time::start, ( args("reset")=false ),
	 "Starts the timer." )
    .def("stop", &Teuchos::Time::stop,
	 "Stops the timer and returns the total elapsed time.")
    .def("totalElapsedTime", &Teuchos::Time::totalElapsedTime,
	 ( args("readCurrentTime")=false ),
	 "Returns the total time accumulated by this timer.\n"
	 "This should only be called when the clock is stopped.")
    .def("reset", &Teuchos::Time::reset,
	 "Resets the cumulative time and number of times this timer has been called.")
    .def("isRunning", &Teuchos::Time::isRunning,
	 "Indicates if this timer is currently running.")
    .def("name", &Teuchos::Time::name, return_value_policy<copy_const_reference>(),
	 "Returns the name of this timer.")
    .def("incrementNumCalls", &Teuchos::Time::incrementNumCalls,
	 "Increment the number of times this timer has been called.")
    .def("numCalls", &Teuchos::Time::numCalls,
	 "Returns the number of times this timer has been called.")
	 
    .def("wallTime", &Teuchos::Time::wallTime,
	 "Returns the current wall-clock time in seconds.")
	.staticmethod("wallTime")
    
    ;
}



void expose_fileinputsource()
{
      // FileInputSource(const string& filename);
    class_<Teuchos::FileInputSource>("foo","Definition of XMLInputSource derived class for reading XML from a file",
                        init<string>() )
        // virtual RefCountPtr<XMLInputStream> stream() const;
        .def( "stream", &Teuchos::FileInputSource::stream ,
                        del_dangleing_rcp<>(),
                        "Create a FileInputStream")
    ;
}


// Define the Teuchos python module
BOOST_PYTHON_MODULE(_Teuchos)
{
    // Teuchos Time support
    expose_time();
    
    
    // dependant on: plist
    //               XMLInputSource
    // tests won't work yet !
    expose_fileinputsource();
    
    extract_teuchos_misc();
}
