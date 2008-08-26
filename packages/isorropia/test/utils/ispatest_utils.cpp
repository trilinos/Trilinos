//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#include <Isorropia_Exception.hpp>

#include <Teuchos_CommandLineProcessor.hpp>

/** ispatest is the namespace that contains isorropia's test-utilities
*/
namespace ispatest {

bool set_verbose(int localProc, int argc, char** argv)
{
  bool result = false;

  Teuchos::CommandLineProcessor clp(true,false);
  clp.setOption( "v", "q", &result, "Set if output is printed or not." );

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
    std::cout << "Proc "<< localProc <<", command-line parse_return "
            << parse_return << std::endl;
    throw Isorropia::Exception("Error parsing command-line.");
  }

  return(result);
}

}//namespace ispatest

