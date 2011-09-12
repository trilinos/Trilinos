// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include <map>

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Shards_Array.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// From test/Utilities directory
#include "Traits.hpp"

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Node)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Node)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(QP)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(QP)

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  Teuchos::GlobalMPISession mpi_session(&argc, &argv);
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Field Tag Testing
    // *********************************************************************
    {

      // Dummy data layouts (same size different name/type)
      RCP<DataLayout> node4 = 
	rcp(new MDALayout<Cell,Node>(100,4));
      RCP<DataLayout> quad4 = 
	rcp(new MDALayout<Cell,QP>(100,4));
      
      // Allocate tags with same name but different topology
      RCP<FieldTag> rcp_nodal_density = rcp(new Tag<double>("density", node4));
      RCP<FieldTag> rcp_qp_density = rcp(new Tag<double>("density", quad4));
      RCP<FieldTag> rcp_grad_qp_density = 
	rcp(new Tag<MyVector<double> >("density", quad4));

      // Get references to field tags
      FieldTag& nodal_density = *rcp_nodal_density;
      FieldTag& qp_density = *rcp_qp_density;
      FieldTag& grad_qp_density = *rcp_grad_qp_density;
    
      // test ostream
      cout << "Printing field tags" << endl;
      cout << nodal_density << endl;
      cout << qp_density << endl;
      cout << grad_qp_density << endl;
      cout << endl;
      
      // test operator ==
      cout << "Are nodal and qp fields equal (should be false)? = " 
	   << (nodal_density == qp_density) << endl;
      TEST_FOR_EXCEPTION(nodal_density == qp_density, std::logic_error,
			 "operator==() failed!");
      
      // New constructor that should be same as nodal_density
      RCP<FieldTag> rcp_nodal_density_copy = 
	rcp(new Tag<double>("density", node4));
      FieldTag& nodal_density_copy = *rcp_nodal_density_copy;

      cout << "Are nodal and nodal copy fields equal (should be true)? = " 
	   << (nodal_density == nodal_density_copy) << endl;
      TEST_FOR_EXCEPTION(!(nodal_density == nodal_density_copy), 
			 std::logic_error,
			 "operator==() failed for unique copy comparison!");
      
      cout << "Are scalar and vector fields "
	   << "equal (should be false)? = " 
	   << (qp_density == grad_qp_density) << endl;
      TEST_FOR_EXCEPTION(qp_density == grad_qp_density, 
			 std::logic_error,
			 "operator==() failed for data layout comparison !");
      
      // test operator =
      cout << "Testing operator=()...";
      Tag<double> op_eq("Garbage", node4);
      TEST_FOR_EXCEPTION(op_eq == nodal_density, std::logic_error, 
			 "Comparison failed.  Should be different!");
      op_eq = dynamic_cast< Tag<double>& >(nodal_density);
      TEST_FOR_EXCEPTION(op_eq != nodal_density, std::logic_error, 
			 "operator=() failed.  Tags should be the same!");
      cout << "Passed." << endl;

      // name() accessor
      cout << "Testing name() accessor...";
      TEST_FOR_EXCEPTION( (nodal_density.name() != std::string("density") ), 
			 std::logic_error,
			 "name() accessor failed!");
      cout << "Passed." << endl;
      
      // dataLayout() accessor
      const DataLayout& tmp = *node4;
      cout << "Testing dataLayout() accessor...";
      TEST_FOR_EXCEPTION(nodal_density.dataLayout() != tmp, 
			 std::logic_error,
			 "dataLayout() accessor failed!");
      cout << "Passed." << endl;
      
      
      // clone()
      cout << "Testing clone()...";
      RCP<FieldTag> my_copy = rcp_nodal_density->clone();
      TEST_FOR_EXCEPTION( *my_copy != *rcp_nodal_density , 
			 std::logic_error,
			 "name() accessor failed!");
      cout << "Passed." << endl;
      
      // Comparison for map key operations
      cout << "Testing stl std::map with key of RCP<FieldTag>...";
      
      map<RCP<FieldTag>, int, FTComp> my_map;
      my_map[rcp_nodal_density] = 0;
      my_map[rcp_qp_density] = 1;
      my_map[rcp_grad_qp_density] = 2;
      
      RCP<FieldTag> tmp_rcp_nodal_density = 
	rcp(new Tag<double>("density", node4));
      RCP<FieldTag> tmp_rcp_qp_density = 
	rcp(new Tag<double>("density", quad4));
      RCP<FieldTag> tmp_rcp_grad_qp_density = 
	rcp(new Tag<MyVector<double> >("density", quad4));
      
      // Make sure we can create a field tag and access matching map entry
      TEST_FOR_EXCEPTION(my_map[tmp_rcp_nodal_density] != 0,
			 std::logic_error,
			 "Failed to find correct FieldTag(0)!");
      TEST_FOR_EXCEPTION(my_map[tmp_rcp_qp_density] != 1,
			 std::logic_error,
			 "Failed to find correct FieldTag(1)!");
      TEST_FOR_EXCEPTION(my_map[tmp_rcp_grad_qp_density] != 2,
			 std::logic_error,
			 "Failed to find correct FieldTag(2)!");

      cout << "Passed." << endl;
      
      // Comparison for vector search operations
      cout << "Testing searching of std::vector of RCP<FieldTag> with find_if...";
      std::vector< Teuchos::RCP<FieldTag> > vec(0);
      vec.push_back(rcp_nodal_density);
      vec.push_back(rcp_qp_density);
      vec.push_back(rcp_grad_qp_density);
      
      { // reference version
	PHX::FTPredRef pred(*rcp_grad_qp_density);
	std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
	  std::find_if(vec.begin(), vec.end(), pred);
	
	TEST_FOR_EXCEPTION(*(*test) !=  *tmp_rcp_grad_qp_density,
			   std::logic_error, 
			   "FTPredRef failed to locate correct FieldTag!");
      }

      { // RCP version
	PHX::FTPred pred(rcp_grad_qp_density);
	std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
	  std::find_if(vec.begin(), vec.end(), pred);
	
	TEST_FOR_EXCEPTION(*(*test) !=  *tmp_rcp_grad_qp_density,
			   std::logic_error, 
			   "FTPred failed to locate correct FieldTag!");
      }

      cout << "Passed." << endl;
      
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
