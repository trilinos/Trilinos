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
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

typedef PHX::MDField<double>::size_type size_type;

// Dimension tags for this problem
struct Dim : public shards::ArrayDimTag {
  Dim(){};
  const char * name() const ;
  static const Dim& tag();
};

struct Quadrature : public shards::ArrayDimTag {
  Quadrature(){};
  const char * name() const ;
  static const Quadrature& tag();
};

struct Node : public shards::ArrayDimTag {
  Node(){};
  const char * name() const ;
  static const Node& tag();
};

struct Cell : public shards::ArrayDimTag {
  Cell(){};
  const char * name() const ;
  static const Cell& tag();
};

struct Point : public shards::ArrayDimTag {
  Point(){};
  const char * name() const ;
  static const Point& tag();
};

const char * Dim::name() const 
{ static const char n[] = "Dim" ; return n ; }
const Dim & Dim::tag() 
{ static const Dim myself ; return myself ; }

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

const char * Point::name() const 
{ static const char n[] = "Point" ; return n ; }
const Point & Point::tag() 
{ static const Point myself ; return myself ; }

template<typename T>
void assign(T& y, const T& x)
{
  for (size_type i=0; i < x.size(); ++i) {
    y[i] = x[i];

    // The following line should cause compilation failure if const is
    // implemented correctly
    //x[i] = y[i];
  }
}


int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of MDField Testing
    // *********************************************************************
    {

      typedef MDField<double,Cell,Node>::size_type size_type;

      // Dummy data layouts
      size_type num_cells = 100;
      RCP<DataLayout> node_scalar = 
	rcp(new MDALayout<Cell,Node>(num_cells,4));
      RCP<DataLayout> quad_scalar = 
	rcp(new MDALayout<Cell,Quadrature>(num_cells,4));
      RCP<DataLayout> node_vector = 
	rcp(new MDALayout<Cell,Node,Dim>(num_cells,4,3));
      RCP<DataLayout> quad_vector = 
	rcp(new MDALayout<Cell,Quadrature,Dim>(num_cells,4,3));

      // Tags with same name but different topology
      Tag<double> nodal_density("density", node_scalar);
      Tag<double> qp_density("density", quad_scalar);
      Tag<double> grad_qp_density("density", quad_vector);
      Tag<MyTraits::FadType> f_grad_qp_density("density",quad_vector);

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Ctors
      cout << "Testing ctor with FieldTag...";
      MDField<double> a(nodal_density);
      MDField<double> b(grad_qp_density);
      cout << "passed!" << endl;
      
      cout << "Testing ctor with individual data...";
      MDField<MyTraits::FadType> c("density", node_scalar);
      MDField<MyTraits::FadType> d("density", quad_vector);
      cout << "passed!" << endl;
      
      cout << "Testing empty ctor...";
      MDField<double> e;
      MDField<MyTraits::FadType> f;
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
      cout << "Testing setFieldData()...";
      ArrayRCP<double> a_mem = 
	arcp<double>(a.fieldTag().dataLayout().size());
      ArrayRCP<double> b_mem = 
	arcp<double>(b.fieldTag().dataLayout().size());
      ArrayRCP<MyTraits::FadType> c_mem = 
	arcp<MyTraits::FadType>(c.fieldTag().dataLayout().size());
      ArrayRCP<MyTraits::FadType> d_mem = 
	arcp<MyTraits::FadType>(d.fieldTag().dataLayout().size());
      ArrayRCP<double> e_mem =
	arcp<double>(e.fieldTag().dataLayout().size());
      ArrayRCP<MyTraits::FadType> f_mem = 
	arcp<MyTraits::FadType>(f.fieldTag().dataLayout().size());

      a.setFieldData(a_mem);
      b.setFieldData(b_mem);
      c.setFieldData(c_mem);
      d.setFieldData(d_mem);
      e.setFieldData(e_mem);
      f.setFieldData(f_mem);

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // rank()
      cout << "Testing rank() method...";
      TEST_FOR_EXCEPTION(a.rank() != 2, std::logic_error,
			 "Rank in a is wrong!");
      TEST_FOR_EXCEPTION(b.rank() != 3, std::logic_error,
			 "Rank in b is wrong!");
      cout << "passed!" << endl;
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // dimension()
      cout << "Testing dimension() method...";
      TEST_FOR_EXCEPTION(b.dimension(0) != num_cells, std::logic_error,
			 "Cell dimesion is wrong!");
      TEST_FOR_EXCEPTION(b.dimension(1) != 4, std::logic_error,
			 "Quadrature dimesion is wrong!");
      TEST_FOR_EXCEPTION(b.dimension(2) != 3, std::logic_error,
			 "Dim dimesion is wrong!");
      cout << "passed!" << endl;
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // dimensions()
      cout << "Testing dimensions() method...";
      std::vector<size_type> dims;
      b.dimensions(dims);
      TEST_FOR_EXCEPTION(dims.size() != 3, std::logic_error,
			 "Number of dimesions is wrong!");
      TEST_FOR_EXCEPTION(dims[0] != 100, std::logic_error,
			 "Number of dimesions is wrong!");
      TEST_FOR_EXCEPTION(dims[1] != 4, std::logic_error,
			 "Number of dimesions is wrong!");
      TEST_FOR_EXCEPTION(dims[2] != 3, std::logic_error,
			 "Number of dimesions is wrong!");
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      cout << "Testing size() method...";
      TEST_FOR_EXCEPTION(a.size() != node_scalar->size(), 
			 std::logic_error, 
			 "Size of array a is not equal to requested size.");
      TEST_FOR_EXCEPTION(b.size() != quad_vector->size(), 
			 std::logic_error, 
			 "Size of array b is not equal to requested size.");
      TEST_FOR_EXCEPTION(c.size() != node_scalar->size(), 
			 std::logic_error, 
			 "Size of array c is not equal to requested size.");
      TEST_FOR_EXCEPTION(d.size() != quad_vector->size(), 
			 std::logic_error, 
			 "Size of array d is not equal to requested size.");
      TEST_FOR_EXCEPTION(e.size() != node_scalar->size(),
			 std::logic_error,
			 "Size of array e is not equal to requested size.");
      TEST_FOR_EXCEPTION(f.size() != quad_vector->size(),
			 std::logic_error,
			 "Size of array f is not equal to requested size.");
      cout << "passed!" << endl;


      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator()
      cout << "Testing operator()(...) accessors...";
      

      RCP<DataLayout> d1 = rcp(new MDALayout<Cell,Dim>(num_cells,2));
    
      RCP<DataLayout> d2 = rcp(new MDALayout<Cell,Dim,Dim>(num_cells,1,2));

      RCP<DataLayout> d3 = 
	rcp(new MDALayout<Cell,Dim,Dim,Dim>(num_cells,1,2,3));

      RCP<DataLayout> d4 = 
	rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4));

      RCP<DataLayout> d5 = 
	rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5));

      RCP<DataLayout> d6 = 
	rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5,6));
      
      RCP<DataLayout> d7 = 
	rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5,6,7));

      MDField<double> f1("Test1",d1);
      MDField<double> f2("Test7",d2);
      MDField<double> f3("Test7",d3);
      MDField<double> f4("Test7",d4);
      MDField<double> f5("Test7",d5);
      MDField<double> f6("Test7",d6);
      MDField<double> f7("Test7",d7);

      ArrayRCP<double> mem1 = arcp<double>(d1->size());
      ArrayRCP<double> mem2 = arcp<double>(d2->size());
      ArrayRCP<double> mem3 = arcp<double>(d3->size());
      ArrayRCP<double> mem4 = arcp<double>(d4->size());
      ArrayRCP<double> mem5 = arcp<double>(d5->size());
      ArrayRCP<double> mem6 = arcp<double>(d6->size());
      ArrayRCP<double> mem7 = arcp<double>(d7->size());
      
      f1.setFieldData(mem1);
      f2.setFieldData(mem2);
      f3.setFieldData(mem3);
      f4.setFieldData(mem4);
      f5.setFieldData(mem5);
      f6.setFieldData(mem6);
      f7.setFieldData(mem7);

      // Access last entry in contiguous array
      f1(99,0) = 1.0;
      f2(99,0,1) = 1.0;
      f3(99,0,1,2) = 1.0;
      f4(99,0,1,2,3) = 1.0;
      f5(99,0,1,2,3,4) = 1.0;
      f6(99,0,1,2,3,4,5) = 1.0;
      f7(99,0,1,2,3,4,5,6) = 1.0;

      // Test const/ non-const versions
      assign(f1,f1);
      assign(f2,f2);
      assign(f3,f3);
      assign(f4,f4);
      assign(f5,f5);
      assign(f6,f6);
      assign(f7,f7);

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator[]
      cout << "Testing operator[]() accessor...";
      
      for (size_type i = 0; i < a.size(); ++i)
	a[i] = 5.0;
      for (size_type i = 0; i < b.size(); ++i)
	b[i] = 5.0;
      for (size_type i = 0; i < c.size(); ++i)
	c[i] = 5.0;
      for (size_type i = 0; i < d.size(); ++i)
	d[i] = MyTraits::FadType(5.0);
      for (size_type i = 0; i < e.size(); ++i)
	e[i] = 5.0;
      for (size_type i = 0; i < f.size(); ++i)
	f[i] = MyTraits::FadType(5.0);

      cout << "passed!" << endl;

      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      cout << "Testing operator<<()...";
      ostringstream output;
      output << a;
      TEST_FOR_EXCEPTION(output.str() != "MDField(100,4): Tag: density, double, DataLayout: MDA<Cell,Node>(100,4)", std::logic_error, "String match failed!"); 
      cout << "passed!" << endl;
      cout << output.str() << endl;

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




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
