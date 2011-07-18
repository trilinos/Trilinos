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

#include <vector>

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef PHX::MDField<double>::size_type size_type;

// ********************************************************
// Dimension tags for this problem
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dim)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Dim)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Quadrature)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Quadrature)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Node)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Node)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Point)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Point)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

// ********************************************************
template<typename VectorType>
void simulated_intrepid_integrate(VectorType& v) {
  
  if (v.rank() == 3)
    v(0,0,0) = 1.0;
  else
    v(0,1) = 1.0;
}


// ********************************************************
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
    // Start of MDField Testing
    // *********************************************************************
    {

      typedef MDField<double,Cell,Node>::size_type size_type;

      std::vector<size_type> dims(3);
      dims[0] = 10;
      dims[1] = 4;
      dims[2] = 3;

      RCP<DataLayout> quad_vector = 
	rcp(new MDALayout<Cell,Quadrature,Dim>(dims[0],dims[1],dims[2]));
      
      int size = quad_vector->size();

      TEST_FOR_EXCEPTION(size != dims[0]*dims[1]*dims[2], std::runtime_error, 
			 "Size mismatch on MDField!");

      ArrayRCP<double> a_mem = arcp<double>(size);
      ArrayRCP<double> b_mem = arcp<double>(size);

      for (int i=0; i < a_mem.size(); ++i)
	a_mem[i] = static_cast<double>(i);

      for (int i=0; i < b_mem.size(); ++i)
	b_mem[i] = static_cast<double>(i);

      MDField<double,Cell,Point,Dim> a("density",quad_vector);
      MDField<double> b("density",quad_vector);

      a.setFieldData(a_mem);
      b.setFieldData(b_mem);

      simulated_intrepid_integrate(a);     
      simulated_intrepid_integrate(b);     

      // ***********************
      // Shards tests
      // ***********************

      ArrayRCP<double> c_mem = arcp<double>(size);
      ArrayRCP<double> d_mem = arcp<double>(size);

      for (int i=0; i < c_mem.size(); ++i)
	c_mem[i] = static_cast<double>(i);

      for (int i=0; i < d_mem.size(); ++i)
	d_mem[i] = static_cast<double>(i);

      shards::Array<double,shards::NaturalOrder,Cell,Node,Dim> c(c_mem.get(),
								 dims[0], 
								 dims[1],
								 dims[2]);

      size_type rank = dims.size();

      const ArrayRCP<const shards::ArrayDimTag*> tags = 
	arcp<const shards::ArrayDimTag*>(rank);
      tags[0] = &Cell::tag();
      tags[1] = &Point::tag();
      tags[2] = &Dim::tag();
      
      shards::Array<double,shards::NaturalOrder> d(d_mem.get(),rank,
						   &dims[0],tags.get());
      
      simulated_intrepid_integrate(d); 
      simulated_intrepid_integrate((const shards::Array<double,shards::NaturalOrder>&)(c));    

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
