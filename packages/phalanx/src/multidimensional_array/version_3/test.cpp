/*
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

*/


#include <iostream>
#include "Array.hpp"

struct Spatial : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Spatial& descriptor();
};

struct Quadrature : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Quadrature& descriptor();
};

struct Node : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Node& descriptor();
};

struct Cell : public phdmesh::ArrayDimTag {
  const char * name() const ;
  static const Cell& descriptor();
};

const char * Spatial::name() const 
{ static const char n[] = "Spatial" ; return n ; }
const Spatial & Spatial::descriptor() 
{ static const Spatial myself ; return myself ; }

const char * Quadrature::name() const 
{ static const char n[] = "Quadrature" ; return n ; }
const Quadrature & Quadrature::descriptor() 
{ static const Quadrature myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::descriptor() 
{ static const Node myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::descriptor() 
{ static const Cell myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace phdmesh;
  
  try {
    
    {
      cout << "\nStarting MultiDimensionalArray Example!\n" << endl;

      std::vector<double> memory(72);
      Array<double,NaturalOrder,Cell,Node,Spatial,Spatial> 
	a(&memory[0], 3, 4, 2, 2);

      std::vector<Array<double,NaturalOrder,Cell,Node,Spatial,Spatial>::size_type> dimensions;

      std::size_t cell_dim = a.dimension(0);
      std::size_t node_dim = a.dimension(1);
      std::size_t row_dim = a.dimension(2);
      std::size_t col_dim = a.dimension(3);

      for (std::size_t cell = 0; cell < cell_dim; ++cell)
	for (std::size_t node = 0; node < node_dim; ++node)
	  for (std::size_t row = 0; row < row_dim; ++row)
	    for (std::size_t col = 0; col < col_dim; ++col)
	      {
		a(cell,node,row,col) = 2.0;
		cout << "a[" << cell <<"," << node << "," 
		     << row << "," << col << "] = " 
		     << a(cell,node,row,col) << endl; 
	      }
      
      cout << endl;

      for (std::size_t i = 0; i < a.size(); ++i)
	{
	  a[i] = 3.0;
	  cout << "a[" << i << "] = " << a[i] << endl;
	}

      // check truncating down to matrix
      Array<double,NaturalOrder,Spatial,Spatial> m = (a.truncate(0)).truncate(0);


      cout << "\nFinished MultiDimensionalArray Example!\n" << endl;
    }

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nRun has completed successfully!\n" << std::endl; 
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

  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
