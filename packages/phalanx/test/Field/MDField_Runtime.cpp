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


#include "Kokkos_DynRankView_Fad.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_MDField_UnmanagedAllocator.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

// Dimension tags for this problem
// Dimension tags for this problem
struct Dim {};
struct Quadrature{};
struct Node{};
struct Cell{};
struct Point{};

namespace PHX {
  template<> struct is_extent<Dim> : std::true_type {};
  template<> struct is_extent<Quadrature> : std::true_type {};
  template<> struct is_extent<Node> : std::true_type {};
  template<> struct is_extent<Cell> : std::true_type {};
  template<> struct is_extent<Point> : std::true_type {};
}

// Demonstrate ability to override print string
namespace PHX {
  template<> std::string print<Dim>(){return "D";}
  template<> std::string print<Quadrature>(){return "Q";}
  template<> std::string print<Node>(){return "N";}
  template<> std::string print<Cell>(){return "C";}
  template<> std::string print<Point>(){return "P";}
}

TEUCHOS_UNIT_TEST(mdfield, RuntimeTimeChecked)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
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

    // Test create, set and get for RCP<FieldTag>
    {
      RCP<PHX::FieldTag> t1 = rcp(new PHX::Tag<double>("test",node_scalar));
      MDField<double> f1(t1); // rcp(tag) ctor
      MDField<double> f2; 
      f2.setFieldTag(t1); // set rcp(tag)
      auto t2 = f1.fieldTagPtr(); // accessor
      auto t3 = f2.fieldTagPtr();
      TEST_ASSERT(nonnull(t1));
      TEST_ASSERT(nonnull(t2));
      TEST_ASSERT(nonnull(t3));
      TEST_EQUALITY(*t1,*t2);
      TEST_EQUALITY(*t2,*t3);
    }

    // MDField ctor interoperability between const and non-const tags
    {
      Tag<double> nonconst_tag("non-const tag", node_scalar);
      Tag<const double> const_tag("const tag", node_scalar);

      // Create a const field from a non-const field tag
      MDField<const double> c_field1(nonconst_tag);

      // Create a non-const field from a const field tag
      MDField<double> nc_field1(const_tag);

      // Create a non-const field from a non-const field tag
      MDField<double> nc_field2(nonconst_tag);

      // Create a const field from a const field tag
      MDField<const double> c_field2(const_tag);
    }

    // Copy constructor from const/non-const MDFields. NOTE: this
    // tests for assignment from both Compiletime and DynRank
    // MDFields.
    {
      RCP<DataLayout> ctor_dl_p  = rcp(new MDALayout<Cell,Point>(10,4));
      MDField<double,Cell,Point> ctor_nonconst_p("ctor_nonconst_p",ctor_dl_p);
      MDField<const double,Cell,Point> ctor_const_p("ctor_const_p",ctor_dl_p);
      MDField<double> ctor_dyn_nonconst_p("ctor_nonconst_p",ctor_dl_p);
      MDField<const double> ctor_dyn_const_p("ctor_const_p",ctor_dl_p);

      MDField<double> cc1(ctor_nonconst_p);       // non-const from non-const
      MDField<const double> cc2(ctor_nonconst_p); // const from non-const
      MDField<const double> cc3(ctor_const_p);    // const from const
      MDField<const double> cc1dyn(ctor_dyn_nonconst_p);    // const from non-const
      MDField<const double> cc2dyn(ctor_dyn_nonconst_p);    // const from non-const
      MDField<const double> cc3dyn(ctor_dyn_const_p);    // const from const

      // NOTE: we allow the tag template types to be DIFFERENT as long
      // as the rank is the same! A field might use the "Point" DimTag
      // but another evaluator might reference the same field using
      // QuadraturePoint DimTag.
      RCP<DataLayout> ctor_dl_qp = rcp(new MDALayout<Cell,Quadrature>(10,4));
      MDField<double,Cell,Quadrature> ctor_nonconst_qp("ctor_nonconst",ctor_dl_qp);
      MDField<const double,Cell,Quadrature> ctor_const_qp("ctor_const_qp",ctor_dl_qp); 
      MDField<double,Cell,Quadrature> ctor_dyn_nonconst_qp("ctor_nonconst",ctor_dl_qp);
      MDField<const double,Cell,Quadrature> ctor_dyn_const_qp("ctor_const_qp",ctor_dl_qp); 

      // Repeat test above but with different tags for Quadrature --> Point
      MDField<double> cc4(ctor_nonconst_qp);       // non-const from non-const
      MDField<const double> cc5(ctor_nonconst_qp); // const from non-const
      MDField<const double> cc6(ctor_const_qp);    // const from const
      MDField<double> cc4dyn(ctor_dyn_nonconst_qp);       // non-const from non-const
      MDField<const double> cc5dyn(ctor_dyn_nonconst_qp); // const from non-const
      MDField<const double> cc6dyn(ctor_dyn_const_qp);    // const from const

      // While we have these objects, lets test the assignment operator as well
      MDField<double> cc7(ctor_nonconst_p);         // non-const from non-const
      MDField<const double> cc8(ctor_nonconst_p);   // const from non-const
      MDField<const double> cc9(ctor_const_p);      // const from const
      MDField<double> cc10(ctor_nonconst_qp);       // non-const from non-const
      MDField<const double> cc11(ctor_nonconst_qp); // const from non-const
      MDField<const double> cc12(ctor_const_qp);    // const from const

      cc7 = ctor_nonconst_p;
      cc8 = ctor_nonconst_p;
      cc9 = ctor_const_p;
      cc10 = ctor_nonconst_qp;
      cc11 = ctor_nonconst_qp;
      cc12 = ctor_const_qp;

      cc7 = ctor_dyn_nonconst_p;
      cc8 = ctor_dyn_nonconst_p;
      cc9 = ctor_dyn_const_p;
      cc10 = ctor_dyn_nonconst_qp;
      cc11 = ctor_dyn_nonconst_qp;
      cc12 = ctor_dyn_const_qp;
    }
    
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
    PHX::any a_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(a.fieldTag());
    PHX::any b_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(b.fieldTag());
    PHX::any c_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(c.fieldTag(),ddims);
    PHX::any d_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(d.fieldTag(),ddims);
    PHX::any e_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(e.fieldTag());
    PHX::any f_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f.fieldTag(),ddims);

    a.setFieldData(a_mem);
    b.setFieldData(b_mem);
    c.setFieldData(c_mem);
    d.setFieldData(d_mem);
    e.setFieldData(e_mem);
    f.setFieldData(f_mem);
    
    out << "passed!" << endl;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // rank()
    TEST_EQUALITY(a.rank(), 2);
    TEST_EQUALITY(b.rank(), 3);
    TEST_EQUALITY(c.rank(), 2);
    TEST_EQUALITY(d.rank(), 3);
    TEST_EQUALITY(e.rank(), 2);
    TEST_EQUALITY(f.rank(), 3);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // extent()
    TEST_EQUALITY(b.extent(0), num_cells);
    TEST_EQUALITY(b.extent(1), 4);
    TEST_EQUALITY(b.extent(2), 3);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // extent_int()
    TEST_EQUALITY(b.extent_int(0), static_cast<int>(num_cells));
    TEST_EQUALITY(b.extent_int(1), static_cast<int>(4));
    TEST_EQUALITY(b.extent_int(2), static_cast<int>(3));
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimension()
    TEST_EQUALITY(b.dimension(0), num_cells);
    TEST_EQUALITY(b.dimension(1), 4);
    TEST_EQUALITY(b.dimension(2), 3);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dimensions()
    std::vector<size_type> dims;
    b.dimensions(dims);
    TEST_EQUALITY(dims.size(), 3);
    TEST_EQUALITY(dims[0], 100);
    TEST_EQUALITY(dims[1], 4);
    TEST_EQUALITY(dims[2], 3);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // size()
    TEST_EQUALITY(a.size(), node_scalar->size());
    TEST_EQUALITY(b.size(), quad_vector->size());
    TEST_EQUALITY(c.size(), node_scalar->size());
    TEST_EQUALITY(d.size(), quad_vector->size());
    TEST_EQUALITY(e.size(), node_scalar->size());
    TEST_EQUALITY(f.size(), quad_vector->size());

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


    // Use unmanaged allocator
    MDField<double> f1 = PHX::allocateUnmanagedMDField<double>("Test1",d1);
    MDField<double> f2 = PHX::allocateUnmanagedMDField<double>("Test2",d2);
    MDField<double> f3 = PHX::allocateUnmanagedMDField<double>("Test3",d3);
    MDField<double> f4 = PHX::allocateUnmanagedMDField<double>("Test4",d4);
    MDField<double> f5 = PHX::allocateUnmanagedMDField<double>("Test5",d5);
    MDField<double> f6 = PHX::allocateUnmanagedMDField<double>("Test6",d6);
    MDField<double> f7 = PHX::allocateUnmanagedMDField<double>("Test7",d7);

    // Pre- Unmnaged allocator
    // MDField<double> f1("Test1",d1);
    // MDField<double> f2("Test2",d2);
    // MDField<double> f3("Test3",d3);
    // MDField<double> f4("Test4",d4);
    // MDField<double> f5("Test5",d5);
    // MDField<double> f6("Test6",d6);
    // MDField<double> f7("Test7",d7);
  
    // f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f1.fieldTag()));
    // f2.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag()));
    // f3.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f3.fieldTag()));
    // f4.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f4.fieldTag()));
    // f5.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f5.fieldTag()));
    // f6.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f6.fieldTag()));
    // f7.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f7.fieldTag()));

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
    // PHX::any.  the field manager always sotres the non-const
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
      c_f1.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f1.fieldTag()));
      c_f2.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f2.fieldTag()));
      c_f3.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f3.fieldTag()));
      c_f4.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f4.fieldTag()));
      c_f5.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f5.fieldTag()));
      c_f6.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f6.fieldTag()));
      c_f7.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(c_f7.fieldTag()));
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

    f1_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f1.fieldTag(), ddims));
    f2_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f2.fieldTag(), ddims));
    f3_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f3.fieldTag(), ddims));
    f4_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f4.fieldTag(), ddims));
    f5_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f5.fieldTag(), ddims));
    f6_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f6.fieldTag(), ddims));
    f7_fad.setFieldData(PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f7.fieldTag(), ddims));

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

    // The bracket operator use in Kokkos is normally liimited to a
    // rank-1 view. Intrepid requires access to the entire array for
    // views with rank greater than 1. This is a very dangerous hack
    // because it allows users to bypass the underlying layout. This
    // accessor will not work correctly with subviews or
    // padding. Kokkos will suppor this for now until we can remove
    // all use of bracket operator from intrepid. Teh below tests
    // verify that this capability works. We can remove the tests
    // below once the bracket op is no longer needed in our stack.

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
    // Do NOT check for array rank enforcement. DynRank MDField sets
    // the rank at runtime and we allow it to be changed!
    //TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag())),PHX::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // kokkos view accessors
    {
      // non-const view
      auto kva = a.get_view(); 
      kva(0,0) = 1.0;
      auto kvc = c.get_view(); 
      kvc(0,0) = MyTraits::FadType(1.0);
      // const view (view const, not const data)
      const auto const_kva = a.get_view(); 
      const_kva(0,0) = 1.0;
      const auto const_kvc = c.get_view(); 
      const_kvc(0,0) = MyTraits::FadType(1.0);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    out << "Testing operator<<()...";
    ostringstream output;
    output << a;
    TEST_EQUALITY(output.str(),"MDField(100,4): Tag: density, double, DataLayout: <C,N>(100,4)");
    out << "passed!" << endl;
    out << output.str() << endl;
  }

  TimeMonitor::summarize();
}
