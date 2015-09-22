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

TEUCHOS_UNIT_TEST(mdfield, RuntimeTimeChecked)
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
    out << "Testing ctor with FieldTag...";
    MDField<double> a(nodal_density);
    MDField<double> b(grad_qp_density);
    out << "passed!" << endl;
    
    out << "Testing ctor with individual data...";
    MDField<MyTraits::FadType> c("density", node_scalar);
    MDField<MyTraits::FadType> d("density", quad_vector);
    out << "passed!" << endl;
    
    out << "Testing empty ctor...";
    MDField<double> e;
    MDField<MyTraits::FadType> f;
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FieldTag accessor
    out << "Testing fieldTag() accessor...";
    
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
    
    out << "passed!" << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // setFieldTag()
    out << "Testing setFieldTag()...";
    e.setFieldTag(nodal_density);
    f.setFieldTag(f_grad_qp_density);
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // setFieldData()
    out << "Testing setFieldData()...";
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
    
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // rank()
    out << "Testing rank() method...";
    TEUCHOS_TEST_FOR_EXCEPTION(a.rank() != 2, std::logic_error,
			       "Rank in a is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.rank() != 3, std::logic_error,
			       "Rank in b is wrong!");
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimension()
    out << "Testing dimension() method...";
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(0) != num_cells, std::logic_error,
			       "Cell dimesion is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(1) != 4, std::logic_error,
			 "Quadrature dimesion is wrong!");
    TEUCHOS_TEST_FOR_EXCEPTION(b.dimension(2) != 3, std::logic_error,
			       "Dim dimesion is wrong!");
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimensions()
    out << "Testing dimensions() method...";
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
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // size()
    out << "Testing size() method...";

   TEUCHOS_TEST_FOR_EXCEPTION(a.size() != node_scalar->size(), 
			       std::logic_error, 
			       "Size of array a is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(b.size() != quad_vector->size(), 
			       std::logic_error, 
			       "Size of array b is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(c.size() != node_scalar->size(), 
			       std::logic_error, 
			       "Size of array c is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(d.size() != quad_vector->size(), 
			       std::logic_error, 
			       "Size of array d is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(e.size() != node_scalar->size(),
			       std::logic_error,
			       "Size of array e is not equal to requested size.");
    TEUCHOS_TEST_FOR_EXCEPTION(f.size() != quad_vector->size(),
			       std::logic_error,
			       "Size of array f is not equal to requested size.");
    out << "passed!" << endl;
 
   

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // operator()
    out << "Testing operator()(...) accessors...";
    
    RCP<DataLayout> d1 = rcp(new MDALayout<Cell>(num_cells));    
    RCP<DataLayout> d2 = rcp(new MDALayout<Cell,Dim>(num_cells,1));    
    RCP<DataLayout> d3 = rcp(new MDALayout<Cell,Dim,Dim>(num_cells,1,2));
    RCP<DataLayout> d4 = rcp(new MDALayout<Cell,Dim,Dim,Dim>(num_cells,1,2,3));
    RCP<DataLayout> d5 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4));
    RCP<DataLayout> d6 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5));
    RCP<DataLayout> d7 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5,6));
   
   //double
 
    MDField<double> f1("Test1",d1);
    MDField<double> f2("Test2",d2);
    MDField<double> f3("Test3",d3);
    MDField<double> f4("Test4",d4);
    MDField<double> f5("Test5",d5);
    MDField<double> f6("Test6",d6);
    MDField<double> f7("Test7",d7);
  
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
    const MDField<double>& cf1 = f1;
    const MDField<double>& cf2 = f2;
    const MDField<double>& cf3 = f3;
    const MDField<double>& cf4 = f4;
    const MDField<double>& cf5 = f5;
    const MDField<double>& cf6 = f6;
    const MDField<double>& cf7 = f7;
    
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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for const mdfield assignment from non-const factory
    // boost::any.  the field manager always sotres the non-const
    // version.
    {
      MDField<const double> c_f1("CONST Test1",d1);
      MDField<const double> c_f2("CONST Test2",d2);
      MDField<const double> c_f3("CONST Test3",d3);
      MDField<const double> c_f4("CONST Test4",d4);
      MDField<const double> c_f5("CONST Test5",d5);
      MDField<const double> c_f6("CONST Test6",d6);
      MDField<const double> c_f7("CONST Test7",d7);
      
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
    //Fad type
    MDField<MyTraits::FadType> f1_fad("FTest1",d1);
    MDField<MyTraits::FadType> f2_fad("FTest2",d2);
    MDField<MyTraits::FadType> f3_fad("FTest3",d3);
    MDField<MyTraits::FadType> f4_fad("FTest4",d4);
    MDField<MyTraits::FadType> f5_fad("FTest5",d5);
    MDField<MyTraits::FadType> f6_fad("FTest6",d6);
    MDField<MyTraits::FadType> f7_fad("FTest7",d7);

    f1_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f1.fieldTag(), ddims));
    f2_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f2.fieldTag(), ddims));
    f3_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f3.fieldTag(), ddims));
    f4_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f4.fieldTag(), ddims));
    f5_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f5.fieldTag(), ddims));
    f6_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f6.fieldTag(), ddims));
    f7_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,PHX::Device>::buildView(f7.fieldTag(), ddims));

    // Access last entry in contiguous array
    f1_fad(99) = MyTraits::FadType(1.0);
    f2_fad(99,0) = MyTraits::FadType(1.0);
    f3_fad(99,0,1) = MyTraits::FadType(1.0);
    f4_fad(99,0,1,2) = MyTraits::FadType(1.0);
    f5_fad(99,0,1,2,3) = MyTraits::FadType(1.0);
    f6_fad(99,0,1,2,3,4) = MyTraits::FadType(1.0);
    f7_fad(99,0,1,2,3,4,5) = MyTraits::FadType(1.0);

    // Test const/ non-const versions
    const MDField<MyTraits::FadType>& cf1_fad = f1_fad;
    const MDField<MyTraits::FadType>& cf2_fad = f2_fad;
    const MDField<MyTraits::FadType>& cf3_fad = f3_fad;
    const MDField<MyTraits::FadType>& cf4_fad = f4_fad;
    const MDField<MyTraits::FadType>& cf5_fad = f5_fad;
    const MDField<MyTraits::FadType>& cf6_fad = f6_fad;
    const MDField<MyTraits::FadType>& cf7_fad = f7_fad;

    for (size_type i=0; i < f7_fad.dimension(0); ++i) {
      f1_fad(i) = cf1_fad(i);
      for (size_type j=0; j < f7_fad.dimension(1); ++j) {
        f2_fad(i,j) = cf2_fad(i,j);
        for (size_type k=0; k < f7_fad.dimension(2); ++k) {
          f3_fad(i,j,k) = cf3_fad(i,j,k);
          for (size_type l=0; l < f7_fad.dimension(3); ++l) {
            f4_fad(i,j,k,l) = cf4_fad(i,j,k,l);
            for (size_type m=0; m < f7_fad.dimension(4); ++m) {
              f5_fad(i,j,k,l,m) = cf5_fad(i,j,k,l,m);
              for (size_type n=0; n < f7_fad.dimension(5); ++n) {
                f6_fad(i,j,k,l,m,n) = cf6_fad(i,j,k,l,m,n);
                for (size_type o=0; o < f7_fad.dimension(6); ++o) {
                  f7_fad(i,j,k,l,m,n,o) = cf7_fad(i,j,k,l,m,n,o);
                }
              }
            }
          }
        }
      }
    }

    
    out << "passed!" << endl;
      

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    //deep_copy
//    
//    out << "f4"<<f4 <<endl;
//    f4.deep_copy(double(0.0));
//    out << endl<<"f4"<<f4 <<endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // operator[]
    out << "Testing operator[](...) accessors...";

    f1[f1.size()-1] = 2.0;
    f2[f2.size()-1] = 3.0;
    f3[f3.size()-1] = 4.0;
    f4[f4.size()-1] = 5.0;
    f5[f5.size()-1] = 6.0;
    f6[f6.size()-1] = 7.0;

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY(f1[f1.size()-1], cf1[f1.size()-1], tol);
    TEST_FLOATING_EQUALITY(f2[f2.size()-1], cf2[f2.size()-1], tol);
    TEST_FLOATING_EQUALITY(f3[f3.size()-1], cf3[f3.size()-1], tol);
    TEST_FLOATING_EQUALITY(f4[f4.size()-1], cf4[f4.size()-1], tol);
    TEST_FLOATING_EQUALITY(f5[f5.size()-1], cf5[f5.size()-1], tol);
    TEST_FLOATING_EQUALITY(f6[f6.size()-1], cf6[f6.size()-1], tol);

   
    // fad checking
    f1_fad[f1_fad.size()-1] = MyTraits::FadType(2.0);
    f2_fad[f2_fad.size()-1] = MyTraits::FadType(3.0);
    f3_fad[f3_fad.size()-1] = MyTraits::FadType(4.0);
    f4_fad[f4_fad.size()-1] = MyTraits::FadType(5.0);
    f5_fad[f5_fad.size()-1] = MyTraits::FadType(6.0);
    f6_fad[f6_fad.size()-1] = MyTraits::FadType(7.0);
    
    out << "passed!" << endl;
      
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for array rank enforcement
    TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag())),boost::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    out << "Testing operator<<()...";
    ostringstream output;
    output << a;
    TEUCHOS_TEST_FOR_EXCEPTION(output.str() != "MDField(100,4): Tag: density, double, DataLayout: <Cell,Node>(100,4)", std::logic_error, "String match failed!"); 
    out << "passed!" << endl;
    out << output.str() << endl;
  }

  PHX::FinalizeKokkosDevice();  
  TimeMonitor::summarize();
}
