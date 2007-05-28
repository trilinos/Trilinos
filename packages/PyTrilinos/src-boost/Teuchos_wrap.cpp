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
  class_<Teuchos::Time>("Time", init<std::string, optional<bool> >())
    .def("start", &Teuchos::Time::start)
    .def("stop", &Teuchos::Time::stop)
    .def("totalElapsedTime", &Teuchos::Time::totalElapsedTime)
    .def("reset", &Teuchos::Time::reset)
    .def("isRunning", &Teuchos::Time::isRunning)
    //.def("name", &Teuchos::Time::name, return_value_policy<manage_new_object>())
    .def("incrementNumCalls", &Teuchos::Time::incrementNumCalls)
    .def("numCalls", &Teuchos::Time::numCalls)
    .def("wallTime", &Teuchos::Time::wallTime)
    ;
}
