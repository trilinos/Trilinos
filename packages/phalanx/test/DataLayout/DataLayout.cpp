// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
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
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Shards_Array.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Spatial : public shards::ArrayDimTag {
  Spatial(){};
  const char * name() const ;
  static const Spatial & tag();
};

struct Quadrature : public shards::ArrayDimTag {
  Quadrature(){};
  const char * name() const ;
  static const Quadrature & tag();
};

struct Node : public shards::ArrayDimTag {
  Node(){};
  const char * name() const ;
  static const Node & tag();
};

struct Cell : public shards::ArrayDimTag {
  Cell(){};
  const char * name() const ;
  static const Cell & tag();
};

const char * Spatial::name() const 
{ static const char n[] = "Spatial" ; return n ; }
const Spatial & Spatial::tag() 
{ static const Spatial myself ; return myself ; }

const char * Quadrature::name() const 
{ static const char n[] = "Quadrature" ; return n ; }
const Quadrature & Quadrature::tag() 
{ static const Quadrature myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::tag() 
{ static const Node myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::tag() 
{ static const Cell myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      cout << "\n************ Testing MDALayout *****************";

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ctor
      cout << "\nTesting constructor...";
      
      MDALayout<Cell,Cell,Cell,Cell,Cell,Cell,Cell,Cell> 
	rank8(1,2,3,4,5,6,7,8);
      MDALayout<Cell,Cell,Cell,Cell,Cell,Cell,Cell> rank7(1,2,3,4,5,6,7);
      MDALayout<Cell,Cell,Cell,Cell,Cell,Cell> rank6(1,2,3,4,5,6);
      MDALayout<Cell,Cell,Cell,Cell,Cell> rank5(1,2,3,4,5);
      MDALayout<Cell,Cell,Cell,Cell> rank4(1,2,3,4);
      MDALayout<Cell,Cell,Cell> rank3(1,2,3);
      MDALayout<Cell,Cell> rank2(1,2);
      MDALayout<Cell> rank1(1);
      MDALayout<Cell,Node,Spatial,Spatial> n_mat(100,4,2,2);
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // identifier()
      {
	cout << "Testing identifier() accessor...";
	TEST_FOR_EXCEPTION(n_mat.identifier() != 
			   std::string("MDA<Cell,Node,Spatial,Spatial>(100,4,2,2)"), 
			   std::logic_error,
			   "name() accessor failed!");
	cout << "passed!" << endl;
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // rank()
      {
	cout << "Testing rank() accessor...";
	TEST_FOR_EXCEPTION(n_mat.rank() != 4, 
			   std::logic_error,
			   "rank() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // dimensions()
      {
	cout << "Testing dimensions() accessor...";
	std::vector<PHX::DataLayout::size_type> dims;
	n_mat.dimensions(dims);
	TEST_FOR_EXCEPTION(dims[0] != 100, 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(dims[1] != 4, 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(dims[2] != 2, 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(dims[3] != 2, 
			   std::logic_error,
			   "dimensions() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      {
	cout << "Testing size() accessor...";
	TEST_FOR_EXCEPTION(n_mat.size() != 1600, 
			   std::logic_error,
			   "size() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator==()
      {
	cout << "Testing operator==()...";
	
	MDALayout<Cell,Node,Spatial> n_vec_a(100,4,2);
	MDALayout<Cell,Node,Spatial> n_vec_b(50,4,2);
	MDALayout<Cell,Node,Spatial> n_vec_c(50,4,2);
	MDALayout<Cell,Quadrature,Spatial> qp_vec(100,4,2);
	MDALayout<Cell,Quadrature,Spatial,Spatial> qp_mat(100,4,2,2);
	
	// same data layout, different object
	TEST_FOR_EXCEPTION( (n_vec_b != n_vec_c), 
			    std::logic_error,
			    "operator==() failed test 1!");
      
	// same data layout, different dim sizes
	TEST_FOR_EXCEPTION( (n_vec_a == n_vec_c), 
			    std::logic_error,
			    "operator==() failed test 2!");
      
	// different types, same rank, same dim sizes
	TEST_FOR_EXCEPTION( (n_vec_a == qp_vec), 
			    std::logic_error,
			    "operator==() failed test 3!");
      
	// different types
	TEST_FOR_EXCEPTION( (n_vec_a == qp_mat), 
			    std::logic_error,
			    "operator==() failed test 4!");

	cout << "passed!" << endl;
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // dimension() 
      {
	cout << "Testing dimension()...";
	TEST_FOR_EXCEPTION(rank1.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank2.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank2.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank3.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank3.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank3.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank4.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank4.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank4.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank4.dimension(3) != 4, std::logic_error,
			   "dimension() failed to return correct dimension!");
			   
	TEST_FOR_EXCEPTION(rank5.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank5.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank5.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank5.dimension(3) != 4, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank5.dimension(4) != 5, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank6.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank6.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank6.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank6.dimension(3) != 4, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank6.dimension(4) != 5, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank6.dimension(5) != 6, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank7.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(3) != 4, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(4) != 5, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(5) != 6, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank7.dimension(6) != 7, std::logic_error,
			   "dimension() failed to return correct dimension!");

	TEST_FOR_EXCEPTION(rank8.dimension(0) != 1, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(1) != 2, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(2) != 3, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(3) != 4, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(4) != 5, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(5) != 6, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(6) != 7, std::logic_error,
			   "dimension() failed to return correct dimension!");
	TEST_FOR_EXCEPTION(rank8.dimension(7) != 8, std::logic_error,
			   "dimension() failed to return correct dimension!");

	cout << "passed!" << endl;
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      {
	cout << "Testing ostream...";
	ostringstream output;
	output << rank8 << endl;
	output << rank7 << endl;
	output << rank6 << endl;   
	output << rank5 << endl;
	output << rank4 << endl;
	output << rank3 << endl;
	output << rank2 << endl;
	output << rank1 << endl;
	output << n_mat << endl;
	cout << "...passed:\n" << output.str() << endl;
      }

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
