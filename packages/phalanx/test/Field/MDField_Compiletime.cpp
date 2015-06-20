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


#include "Phalanx_config.hpp"
#include "Phalanx.hpp"
#include "Phalanx_DimTag.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

typedef PHX::MDField<double>::size_type size_type;

// Dimension tags for this problem
struct Dim : public PHX::DimTag {
  Dim(){};
  const char * name() const ;
  static const Dim& tag();
};

struct Quadrature : public PHX::DimTag {
  Quadrature(){};
  const char * name() const ;
  static const Quadrature& tag();
};

struct Node : public PHX::DimTag {
  Node(){};
  const char * name() const ;
  static const Node& tag();
};

struct Cell : public PHX::DimTag {
  Cell(){};
  const char * name() const ;
  static const Cell& tag();
};

struct Point : public PHX::DimTag {
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

TEUCHOS_UNIT_TEST(mdfield, CompileTimeChecked)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
  TimeMonitor tm(*total_time);
  
  PHX::InitializeKokkosDevice();
  
  // *********************************************************************
  // Start of MDField Testing
  // *********************************************************************
  {
    
    typedef MDField<double,Cell,Node>::size_type size_type;
    
    // Dummy data layouts
    size_type num_cells = 100;
    RCP<DataLayout> node_scalar = rcp(new MDALayout<Cell,Node>(num_cells,4));
    RCP<DataLayout> quad_scalar = rcp(new MDALayout<Cell,Quadrature>(num_cells,4));
    RCP<DataLayout> node_vector = rcp(new MDALayout<Cell,Node,Dim>(num_cells,4,3));
    RCP<DataLayout> quad_vector = rcp(new MDALayout<Cell,Quadrature,Dim>(num_cells,4,3));
    
    // Tags with same name but different topology
    Tag<double> nodal_density("density", node_scalar);
    Tag<double> qp_density("density", quad_scalar);
    Tag<double> grad_qp_density("density", quad_vector);
    Tag<MyTraits::FadType> f_grad_qp_density("density",quad_vector);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Ctors
    cout << "Testing ctor with FieldTag...";
    MDField<double,Cell,Node> a(nodal_density);
    MDField<double,Cell,Quadrature,Dim> b(grad_qp_density);
    cout << "passed!" << endl;
    
    cout << "Testing ctor with individual data...";
    MDField<MyTraits::FadType,Cell,Node> 
      c("density", node_scalar);
    MDField<MyTraits::FadType,Cell,Quadrature,Dim> 
      d("density", quad_vector);
    cout << "passed!" << endl;
    
    cout << "Testing empty ctor...";
    MDField<double,Cell,Point> e;
    MDField<MyTraits::FadType,Cell,Point,Dim> f;
    cout << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FieldTag accessor
    cout << "Testing fieldTag() accessor...";
    
    const FieldTag& test_a = a.fieldTag();
    TEUCHOS_TEST_FOR_EXCEPTION( !(test_a == nodal_density),
				std::logic_error,
				"fieldTag() accessor failed!");
    
    const FieldTag& test_b = b.fieldTag();
    TEUCHOS_TEST_FOR_EXCEPTION( !(test_b == grad_qp_density),
				std::logic_error,
				"fieldTag() accessor failed!");
    
    const FieldTag& test_d = d.fieldTag();
    TEUCHOS_TEST_FOR_EXCEPTION( !(test_d == f_grad_qp_density),
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
    const size_type derivative_dim = 8;
    const std::vector<PHX::index_size_type> ddims(1,derivative_dim);
    boost::any a_mem = PHX::KokkosViewFactory<double,PHX::Device>::buildView(a.fieldTag());
    boost::any b_mem = PHX::KokkosViewFactory<double,PHX::Device>::buildView(b.fieldTag());
    boost::any c_mem = PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(c.fieldTag(),ddims);
    boost::any d_mem = PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(d.fieldTag(),ddims);
    boost::any e_mem = PHX::KokkosViewFactory<double,PHX::Device>::buildView(e.fieldTag());
    boost::any f_mem = PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f.fieldTag(),ddims);
    
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
    TEUCHOS_TEST_FOR_EXCEPTION(a.rank() != 2, std::logic_error,
			       "Rank in a is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.rank() != 3, std::logic_error,
			       "Rank in b is wrong!");
    cout << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimension()
    cout << "Testing dimension() method...";
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(0) != num_cells, std::logic_error,
			       "Cell dimesion is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(1) != 4, std::logic_error,
			       "Quadrature dimesion is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(2) != 3, std::logic_error,
			       "Dim dimesion is wrong!");
    cout << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimensions()
    cout << "Testing dimensions() method...";
    std::vector<size_type> dims;
    b.dimensions(dims);
    TEUCHOS_TEST_FOR_EXCEPTION(dims.size() != 3, std::logic_error,
			       "Number of dimesions is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(dims[0] != 100, std::logic_error,
			       "Number of dimesions is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(dims[1] != 4, std::logic_error,
			       "Number of dimesions is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(dims[2] != 3, std::logic_error,
			       "Number of dimesions is wrong!");
    cout << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // size()
    cout << "Testing size() method...";

    TEUCHOS_TEST_FOR_EXCEPTION(a.size() != node_scalar->size(), 
			       std::logic_error, 
			       "Size of array a is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(b.size() != quad_vector->size(), 
			       std::logic_error, 
			       "Size of array b is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(c.size() != node_scalar->size(), 
			       std::logic_error, 
			       "Size of array c is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(d.size() != quad_vector->size() , 
			       std::logic_error, 
			       "Size of array d is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(e.size() != node_scalar->size(),
			       std::logic_error,
			       "Size of array e is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(f.size() != quad_vector->size() ,
			       std::logic_error,
			       "Size of array f is not equal to requested size.");
    cout << "passed!" << endl;
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // operator()
    cout << "Testing operator()(...) accessors...";
    
    RCP<DataLayout> d1 = rcp(new MDALayout<Cell>(num_cells));
    RCP<DataLayout> d2 = rcp(new MDALayout<Cell,Dim>(num_cells,1));    
    RCP<DataLayout> d3 = rcp(new MDALayout<Cell,Dim,Dim>(num_cells,1,2));
    RCP<DataLayout> d4 = rcp(new MDALayout<Cell,Dim,Dim,Dim>(num_cells,1,2,3));
    RCP<DataLayout> d5 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4));
    RCP<DataLayout> d6 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5));
    RCP<DataLayout> d7 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5,6));
    
    MDField<double,Cell> f1("Test1",d1);
    MDField<double,Cell,Dim> f2("Test2",d2);
    MDField<double,Cell,Dim,Dim> f3("Test3",d3);
    MDField<double,Cell,Dim,Dim,Dim> f4("Test4",d4);
    MDField<double,Cell,Dim,Dim,Dim,Dim> f5("Test5",d5);
    MDField<double,Cell,Dim,Dim,Dim,Dim,Dim> f6("Test6",d6);
    MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim> f7("Test7",d7);
    
    f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f1.fieldTag()));
    f2.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag()));
    f3.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f3.fieldTag()));
    f4.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f4.fieldTag()));
    f5.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f5.fieldTag()));
    f6.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f6.fieldTag()));
    f7.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f7.fieldTag()));

    // Access last entry in contiguous array
    f1(99) = 1.0;
    f2(99,0) = 1.0;
    f3(99,0,1) = 1.0;
    f4(99,0,1,2) = 1.0;
    f5(99,0,1,2,3) = 1.0;
    f6(99,0,1,2,3,4) = 1.0;
    f7(99,0,1,2,3,4,5) = 1.0;
    
    // Test const/ non-const versions
    const MDField<double,Cell>& cf1 = f1;
    const MDField<double,Cell,Dim>& cf2 = f2;
    const MDField<double,Cell,Dim,Dim>& cf3 = f3;
    const MDField<double,Cell,Dim,Dim,Dim>& cf4 = f4;
    const MDField<double,Cell,Dim,Dim,Dim,Dim>& cf5 = f5;
    const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim>& cf6 = f6;
    const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim>& cf7 = f7;
    
    for (size_type i=0; i < f7.dimension(0); ++i) {
      f1(i) = cf1(i);      
      for (size_type j=0; j < f7.dimension(1); ++j) {
	f2(i,j) = cf2(i,j);
	for (size_type k=0; k < f7.dimension(2); ++k) {
	  f3(i,j,k) = cf3(i,j,k);
	  for (size_type l=0; l < f7.dimension(3); ++l) {
	    f4(i,j,k,l) = cf4(i,j,k,l);
	    for (size_type m=0; m < f7.dimension(4); ++m) {
	      f5(i,j,k,l,m) = cf5(i,j,k,l,m);
	      for (size_type n=0; n < f7.dimension(5); ++n) {
		f6(i,j,k,l,m,n) = cf6(i,j,k,l,m,n);
		for (size_type o=0; o < f7.dimension(6); ++o) {
		  f7(i,j,k,l,m,n,o) = cf7(i,j,k,l,m,n,o);
		}
	      }
	    }
	  }
	}
      }
    }
    
    cout << "passed!" << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for const mdfield assignment from non-const factory
    // boost::any.  the field manager always sotres the non-const
    // version.
    {
      MDField<const double,Cell> c_f1("CONST Test1",d1);
      MDField<const double,Cell,Dim> c_f2("CONST Test2",d2);
      MDField<const double,Cell,Dim,Dim> c_f3("CONST Test3",d3);
      MDField<const double,Cell,Dim,Dim,Dim> c_f4("CONST Test4",d4);
      MDField<const double,Cell,Dim,Dim,Dim,Dim> c_f5("CONST Test5",d5);
      MDField<const double,Cell,Dim,Dim,Dim,Dim,Dim> c_f6("CONST Test6",d6);
      MDField<const double,Cell,Dim,Dim,Dim,Dim,Dim,Dim> c_f7("CONST Test7",d7);
      
      // Note that the factory never uses a const scalar type
      c_f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f1.fieldTag()));
      c_f2.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f2.fieldTag()));
      c_f3.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f3.fieldTag()));
      c_f4.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f4.fieldTag()));
      c_f5.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f5.fieldTag()));
      c_f6.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f6.fieldTag()));
      c_f7.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c_f7.fieldTag()));
    }




    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for debug build array rank enforcement
    TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag())),boost::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    cout << "Testing operator<<()...";
    ostringstream output;
    output << a;
    TEUCHOS_TEST_FOR_EXCEPTION(output.str() != "MDField<Cell,Node>(100,4): Tag: density, double, DataLayout: <Cell,Node>(100,4)", std::logic_error, "String match failed!"); 
    cout << "passed!" << endl;
    cout << output.str() << endl;
  }

  PHX::FinalizeKokkosDevice();  

  TimeMonitor::summarize();
}
