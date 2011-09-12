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
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
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
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Field Tag Testing
    // *********************************************************************
    {

      // Dummy data layouts
      RCP<DataLayout> node4 = rcp(new MDALayout<Cell,Node>(25,4));
      RCP<DataLayout> quad4 = rcp(new MDALayout<Cell,QP>(25,4));
      const int size = node4->size();
      TEST_FOR_EXCEPTION(node4->size() != quad4->size(), std::runtime_error,
			 "Array sizes fixed to 100 for this test!");
      
      // Tags with same name but different topology
      Tag<double> nodal_density("density", node4);
      Tag<double> qp_density("density", quad4);
      Tag< MyVector<double> > grad_qp_density("density", quad4);
      Tag< MyVector<MyTraits::FadType> > f_grad_qp_density("density", quad4);
      

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Ctors
      cout << "Testing ctor with FieldTag...";
      Field<double> a(nodal_density);
      Field< MyVector<double> > b(grad_qp_density);
      cout << "passed!" << endl;
      
      cout << "Testing ctor with individual data...";
      Field<MyTraits::FadType> c("density", node4);
      Field< MyVector<MyTraits::FadType> > d("density", quad4);
      cout << "passed!" << endl;
      
      cout << "Testing empty ctor...";
      Field<double> e;
      Field< MyVector<MyTraits::FadType> > f;
      cout << "passed!" << endl;
      
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // FieldTag accessor
      cout << "Testing fieldTag() accessor...";

      const FieldTag& test_a = a.fieldTag();
      TEST_FOR_EXCEPTION( !(test_a == nodal_density),
			  std::logic_error,
			  "fieldTag() accessor failed!");
      
      const FieldTag& test_b = b.fieldTag();
      TEST_FOR_EXCEPTION( !(test_b == grad_qp_density),
			  std::logic_error,
			  "fieldTag() accessor failed!");
      
      const FieldTag& test_d = d.fieldTag();
      TEST_FOR_EXCEPTION( !(test_d == f_grad_qp_density),
			  std::logic_error,
			  "fieldTag() accessor failed!");
      
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // setFieldTag()
      cout << "Testing setFieldTag()...";
      e.setFieldTag(nodal_density);
      f.setFieldTag(f_grad_qp_density);
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // setFieldData()
      cout << "Testing getFieldData()...";
      ArrayRCP<double> a_scalar_scalar = 
	arcp<double>(size);
      ArrayRCP< MyVector<double> > b_vector_scalar = 
	arcp< MyVector<double> >(size);
      ArrayRCP<MyTraits::FadType> c_scalar_fad = 
	arcp<MyTraits::FadType>(size);
      ArrayRCP< MyVector<MyTraits::FadType> > d_vector_fad = 
	arcp< MyVector<MyTraits::FadType> >(size);
      ArrayRCP<double> e_scalar_scalar = 
	arcp<double>(size);
      ArrayRCP< MyVector<MyTraits::FadType> > f_vector_fad = 
	arcp< MyVector<MyTraits::FadType> >(size);

      a.setFieldData(a_scalar_scalar);
      b.setFieldData(b_vector_scalar);
      c.setFieldData(c_scalar_fad);
      d.setFieldData(d_vector_fad);
      e.setFieldData(e_scalar_scalar);
      f.setFieldData(f_vector_fad);
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      cout << "Testing size() method...";
      TEST_FOR_EXCEPTION( a.size() != size , std::logic_error, "Size of array a is not equal to requested size.");
      TEST_FOR_EXCEPTION( b.size() != size , std::logic_error, "Size of array b is not equal to requested size.");
      TEST_FOR_EXCEPTION( c.size() != size , std::logic_error, "Size of array c is not equal to requested size.");
      TEST_FOR_EXCEPTION( d.size() != size , std::logic_error, "Size of array d is not equal to requested size.");
      TEST_FOR_EXCEPTION( e.size() != size , std::logic_error, "Size of array e is not equal to requested size.");
      TEST_FOR_EXCEPTION( f.size() != size , std::logic_error, "Size of array f is not equal to requested size.");

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator[]
      cout << "Testing operator[]() accessor...";
      
      for (int i = 0; i < a.size(); ++i) {
	a[i] = 5.0;
	b[i] = 5.0;
	c[i] = 5.0;
	d[i] = MyVector<MyTraits::FadType>(5.0);
	e[i] = 5.0;
	f[i] = MyVector<MyTraits::FadType>(5.0);
      }

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      cout << "Testing operator<<()...";
      ostringstream output;
      output << a << endl;
      cout << "passed!" << endl;
      //cout << output.str() << endl; 
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
