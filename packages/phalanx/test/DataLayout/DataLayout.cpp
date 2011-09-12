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


#include "Phalanx_ConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

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

struct Ordinal1 : public shards::ArrayDimTag {
  Ordinal1(){};
  const char * name() const ;
  static const Ordinal1 & tag();
};

struct Ordinal2 : public shards::ArrayDimTag {
  Ordinal2(){};
  const char * name() const ;
  static const Ordinal2 & tag();
};

struct Ordinal3 : public shards::ArrayDimTag {
  Ordinal3(){};
  const char * name() const ;
  static const Ordinal3 & tag();
};

struct Ordinal4 : public shards::ArrayDimTag {
  Ordinal4(){};
  const char * name() const ;
  static const Ordinal4 & tag();
};

struct Ordinal5 : public shards::ArrayDimTag {
  Ordinal5(){};
  const char * name() const ;
  static const Ordinal5 & tag();
};

struct Ordinal6 : public shards::ArrayDimTag {
  Ordinal6(){};
  const char * name() const ;
  static const Ordinal6 & tag();
};

struct Ordinal7 : public shards::ArrayDimTag {
  Ordinal7(){};
  const char * name() const ;
  static const Ordinal7 & tag();
};

struct Ordinal8 : public shards::ArrayDimTag {
  Ordinal8(){};
  const char * name() const ;
  static const Ordinal8 & tag();
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

const char * Ordinal1::name() const 
{ static const char n[] = "Ordinal1" ; return n ; }
const Ordinal1 & Ordinal1::tag() 
{ static const Ordinal1 myself ; return myself ; }

const char * Ordinal2::name() const 
{ static const char n[] = "Ordinal2" ; return n ; }
const Ordinal2 & Ordinal2::tag() 
{ static const Ordinal2 myself ; return myself ; }

const char * Ordinal3::name() const 
{ static const char n[] = "Ordinal3" ; return n ; }
const Ordinal3 & Ordinal3::tag() 
{ static const Ordinal3 myself ; return myself ; }

const char * Ordinal4::name() const 
{ static const char n[] = "Ordinal4" ; return n ; }
const Ordinal4 & Ordinal4::tag() 
{ static const Ordinal4 myself ; return myself ; }

const char * Ordinal5::name() const 
{ static const char n[] = "Ordinal5" ; return n ; }
const Ordinal5 & Ordinal5::tag() 
{ static const Ordinal5 myself ; return myself ; }

const char * Ordinal6::name() const 
{ static const char n[] = "Ordinal6" ; return n ; }
const Ordinal6 & Ordinal6::tag() 
{ static const Ordinal6 myself ; return myself ; }

const char * Ordinal7::name() const 
{ static const char n[] = "Ordinal7" ; return n ; }
const Ordinal7 & Ordinal7::tag() 
{ static const Ordinal7 myself ; return myself ; }

const char * Ordinal8::name() const 
{ static const char n[] = "Ordinal8" ; return n ; }
const Ordinal8 & Ordinal8::tag() 
{ static const Ordinal8 myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

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
      
      MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7,Ordinal8> rank8(1,2,3,4,5,6,7,8);
      MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7> rank7(1,2,3,4,5,6,7);
      MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6> rank6(1,2,3,4,5,6);
      MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5> rank5(1,2,3,4,5);
      MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4> rank4(1,2,3,4);
      MDALayout<Ordinal1,Ordinal2,Ordinal3> rank3(1,2,3);
      MDALayout<Ordinal1,Ordinal2> rank2(1,2);
      MDALayout<Ordinal1> rank1(1);
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
      // names()
      {
	cout << "Testing names() accessor...";
	std::vector<std::string> names;
	n_mat.names(names);
	TEST_FOR_EXCEPTION(names[0] != "Cell", 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(names[1] != "Node", 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(names[2] != "Spatial", 
			   std::logic_error,
			   "dimensions() accessor failed!");
	TEST_FOR_EXCEPTION(names[3] != "Spatial", 
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
      // name() 
      {
	cout << "Testing name()...";


	cout << "\n" << rank1.name(0) << std::endl;



	TEST_FOR_EXCEPTION(rank1.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank2.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank2.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank3.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank3.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank3.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank4.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank4.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank4.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank4.name(3) != "Ordinal4", std::logic_error,
			   "name() failed to return correct name!");
			   
	TEST_FOR_EXCEPTION(rank5.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank5.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank5.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank5.name(3) != "Ordinal4", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank5.name(4) != "Ordinal5", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank6.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank6.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank6.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank6.name(3) != "Ordinal4", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank6.name(4) != "Ordinal5", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank6.name(5) != "Ordinal6", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank7.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(3) != "Ordinal4", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(4) != "Ordinal5", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(5) != "Ordinal6", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank7.name(6) != "Ordinal7", std::logic_error,
			   "name() failed to return correct name!");

	TEST_FOR_EXCEPTION(rank8.name(0) != "Ordinal1", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(1) != "Ordinal2", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(2) != "Ordinal3", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(3) != "Ordinal4", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(4) != "Ordinal5", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(5) != "Ordinal6", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(6) != "Ordinal7", std::logic_error,
			   "name() failed to return correct name!");
	TEST_FOR_EXCEPTION(rank8.name(7) != "Ordinal8", std::logic_error,
			   "name() failed to return correct name!");

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
