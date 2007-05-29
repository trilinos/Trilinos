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

// Define the Teuchos python module
BOOST_PYTHON_MODULE(_Teuchos)
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
    .def("name", &Teuchos::Time::name, return_internal_reference<>(),
	 "Returns the name of this timer.")
    .def("incrementNumCalls", &Teuchos::Time::incrementNumCalls,
	 "Increment the number of times this timer has been called.")
    .def("numCalls", &Teuchos::Time::numCalls,
	 "Returns the number of times this timer has been called.")
    .def("wallTime", &Teuchos::Time::wallTime,
	 "Returns the current wall-clock time in seconds.")
    ;
}
