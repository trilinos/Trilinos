// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Phalanx_ConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_Generic.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Data Layout Testing
    // *********************************************************************
    {

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ctor
      cout << "\nTesting constructor...";
      
      RCP<DataLayout> node4 = 
	rcp(new Generic("Nodes", 4));
      
      RCP<DataLayout> quad4 = 
	rcp(new Generic("QuadPoints", 4));
      
      RCP<DataLayout> quad9 = 
	rcp(new Generic("QuadPoints", 9));
      
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // name()
      {
	cout << "Testing name() accessor...";
	RCP<Generic> g_node4 = rcp_dynamic_cast<Generic>(node4);
	TEST_FOR_EXCEPTION(is_null(g_node4), 
			   std::logic_error,
			   "dynamic cast from DataLayout to Generic failed!");
	
	TEST_FOR_EXCEPTION(g_node4->name() != std::string("Nodes"), 
			   std::logic_error,
			   "name() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      {
	cout << "Testing size() accessor...";
	TEST_FOR_EXCEPTION(node4->size() != 4, 
			   std::logic_error,
			   "size() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator==()
      {
	cout << "Testing operator==()...";
	
	RCP<DataLayout> unique_node4_copy = rcp(new Generic("Nodes", 4));
	
	// same data layout, different object
	TEST_FOR_EXCEPTION( !(*node4 == *unique_node4_copy), 
			    std::logic_error,
			    "operator==() failed for separate objects!");
	
	// different data layouts - name differentiation
	TEST_FOR_EXCEPTION( *node4 == *quad4, 
			    std::logic_error,
			    "operator==() failed for name differentiation!");
	
	// same name, different size 
	TEST_FOR_EXCEPTION( *quad4 == *quad9, 
			    std::logic_error,
			    "operator==() failed for size differentiation!");
	
	cout << "passed!" << endl;
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      cout << "Testing ostream...";
      ostringstream output;
      output << *node4;
      cout << "...passed: " << output.str() << endl;

    }

    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
