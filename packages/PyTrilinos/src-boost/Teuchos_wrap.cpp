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
#include "Teuchos_ConfigDefs.hpp"

void extract_teuchos_misc();
void convert_exceptions();

void expose_time();
void expose_fileinputsource();
void expose_plist();
void expose_xmlobject();
void expose_xml_r_w();
void expose_str_inputsource();
void expose_scalartraits();


// Define the Teuchos python module
BOOST_PYTHON_MODULE(_Teuchos)
{
    // Teuchos version support
    def("Teuchos_Version", Teuchos::Teuchos_Version);
    
    // Teuchos Time support
    expose_time();
    
    // dependant on: plist
    //               XMLInputSource
    // tests won't work yet !
    expose_fileinputsource();
    //expose_plist();
    expose_xmlobject();
    expose_xml_r_w();
    expose_str_inputsource();
    expose_scalartraits();
    
    extract_teuchos_misc();
    convert_exceptions();
}
