
#include "Teuchos_Time.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;


void expose_time()
{

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
