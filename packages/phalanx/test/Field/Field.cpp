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

#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_Field_UnmanagedAllocator.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Field.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

TEUCHOS_UNIT_TEST(field, all)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
  TimeMonitor tm(*total_time);

  // *********************************************************************
  // Start of Field Testing
  // *********************************************************************
  {
    typedef Field<double,3>::size_type size_type;

    // Dummy data layouts
    size_type num_cells = 100;
    RCP<DataLayout> node_scalar = rcp(new Layout("NODE_SCALAR",num_cells,4));
    RCP<DataLayout> quad_scalar = rcp(new Layout("QUAD_SCALAR",num_cells,4));
    RCP<DataLayout> node_vector = rcp(new Layout("NODE_VECTOR",num_cells,4,3));
    RCP<DataLayout> quad_vector = rcp(new Layout("QUAD_VECTOR",num_cells,4,3));

    // Tags with same name but different topology
    Tag<double> nodal_density("density", node_scalar);
    Tag<double> qp_density("density", quad_scalar);
    Tag<double> grad_qp_density("density", quad_vector);
    Tag<MyTraits::FadType> f_grad_qp_density("density",quad_vector);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Ctors
    cout << "Testing ctor with FieldTag...";
    Field<double,2> a(nodal_density);
    Field<double,3> b(grad_qp_density);
    cout << "passed!" << endl;

    cout << "Testing ctor with individual data...";
    Field<MyTraits::FadType,2> c("density", node_scalar);
    Field<MyTraits::FadType,3> d("density", quad_vector);
    cout << "passed!" << endl;

    cout << "Testing empty ctor...";
    Field<double,2> e;
    Field<MyTraits::FadType,3> f;
    cout << "passed!" << endl;

    // Test create, set and get for RCP<FieldTag>
    {
      RCP<PHX::FieldTag> t1 = rcp(new PHX::Tag<double>("test",node_scalar));
      Field<double,2> f1(t1); // rcp(tag) ctor
      Field<double,2> f2; 
      f2.setFieldTag(t1); // set rcp(tag)
      auto t2 = f1.fieldTagPtr(); // accessor
      auto t3 = f2.fieldTagPtr();
      TEST_ASSERT(nonnull(t1));
      TEST_ASSERT(nonnull(t2));
      TEST_ASSERT(nonnull(t3));
      TEST_EQUALITY(*t1,*t2);
      TEST_EQUALITY(*t2,*t3);
    }

    // Field ctor interoperability between const and non-const tags
    {
      Tag<double> nonconst_tag("non-const tag", node_scalar);
      Tag<const double> const_tag("const tag", node_scalar);

      // Create a const field from a non-const field tag
      Field<const double,2> c_field1(nonconst_tag);

      // Create a non-const field from a const field tag
      Field<double,2> nc_field1(const_tag);

      // Create a non-const field from a non-const field tag
      Field<double,2> nc_field2(nonconst_tag);

      // Create a const field from a const field tag
      Field<const double,2> c_field2(const_tag);
    }

    // Copy constructor/assignment operator from const/non-const Fields
    {
      RCP<DataLayout> ctor_dl_p  = rcp(new Layout("test_layout",10,4));
      Field<double,2> ctor_nonconst_p("ctor_nonconst_p",ctor_dl_p);
      Field<const double,2> ctor_const_p("ctor_const_p",ctor_dl_p);

      Field<double,2> cc1(ctor_nonconst_p);       // non-const from non-const
      Field<const double,2> cc2(ctor_nonconst_p); // const from non-const
      Field<const double,2> cc3(ctor_const_p);    // const from const

      // While we have these objects, lets test the assignment operator as well
      Field<double,2> cc7(ctor_nonconst_p);         // non-const from non-const
      Field<const double,2> cc8(ctor_nonconst_p);   // const from non-const
      Field<const double,2> cc9(ctor_const_p);      // const from const
      cc7 = ctor_nonconst_p;
      cc8 = ctor_nonconst_p;
      cc9 = ctor_const_p;
    }

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

    cout << "passed!" << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // rank()
    TEST_EQUALITY(a.rank(),2);
    TEST_EQUALITY(b.rank(),3);
    TEST_EQUALITY(c.rank(),2);
    TEST_EQUALITY(d.rank(),3);
    TEST_EQUALITY(e.rank(),2);
    TEST_EQUALITY(f.rank(),3);

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
    // size()
    TEST_EQUALITY(a.size(), node_scalar->size());
    TEST_EQUALITY(b.size(), quad_vector->size());
    TEST_EQUALITY(c.size(), node_scalar->size());
    TEST_EQUALITY(d.size(), quad_vector->size());
    TEST_EQUALITY(e.size(), node_scalar->size());
    TEST_EQUALITY(f.size(), quad_vector->size());

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // operator()
    cout << "Testing operator()(...) accessors...";

    RCP<DataLayout> d1 = rcp(new Layout("d1",num_cells));
    RCP<DataLayout> d2 = rcp(new Layout("d2",num_cells,1));
    RCP<DataLayout> d3 = rcp(new Layout("d3",num_cells,1,2));
    RCP<DataLayout> d4 = rcp(new Layout("d4",num_cells,1,2,3));
    RCP<DataLayout> d5 = rcp(new Layout("d5",num_cells,1,2,3,4));
    RCP<DataLayout> d6 = rcp(new Layout("d6",num_cells,1,2,3,4,5));
    RCP<DataLayout> d7 = rcp(new Layout("d7",num_cells,1,2,3,4,5,6));

    // Use Unmanaged allocator
    Field<double,1> f1 = PHX::allocateUnmanagedField<double,1>("Test1",d1);
    Field<double,2> f2 = PHX::allocateUnmanagedField<double,2>("Test2",d2);
    Field<double,3> f3 = PHX::allocateUnmanagedField<double,3>("Test3",d3);
    Field<double,4> f4 = PHX::allocateUnmanagedField<double,4>("Test4",d4);
    Field<double,5> f5 = PHX::allocateUnmanagedField<double,5>("Test5",d5);
    Field<double,6> f6 = PHX::allocateUnmanagedField<double,6>("Test6",d6);
    Field<double,7> f7 = PHX::allocateUnmanagedField<double,7>("Test7",d7);

    // Access last entry in contiguous array
    f1(99) = 1.0;
    f2(99,0) = 1.0;
    f3(99,0,1) = 1.0;
    f4(99,0,1,2) = 1.0;
    f5(99,0,1,2,3) = 1.0;
    f6(99,0,1,2,3,4) = 1.0;
    f7(99,0,1,2,3,4,5) = 1.0;

    // Test const/ non-const versions
    const Field<double,1>& cf1 = f1;
    const Field<double,2>& cf2 = f2;
    const Field<double,3>& cf3 = f3;
    const Field<double,4>& cf4 = f4;
    const Field<double,5>& cf5 = f5;
    const Field<double,6>& cf6 = f6;
    const Field<double,7>& cf7 = f7;

    for (size_type i=0; i < f7.extent(0); ++i) {
      f1(i) = cf1(i);
      for (size_type j=0; j < f7.extent(1); ++j) {
	f2(i,j) = cf2(i,j);
	for (size_type k=0; k < f7.extent(2); ++k) {
	  f3(i,j,k) = cf3(i,j,k);
	  for (size_type l=0; l < f7.extent(3); ++l) {
	    f4(i,j,k,l) = cf4(i,j,k,l);
	    for (size_type m=0; m < f7.extent(4); ++m) {
	      f5(i,j,k,l,m) = cf5(i,j,k,l,m);
	      for (size_type n=0; n < f7.extent(5); ++n) {
		f6(i,j,k,l,m,n) = cf6(i,j,k,l,m,n);
		for (size_type o=0; o < f7.extent(6); ++o) {
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
    // PHX::any.  the field manager always stores the non-const
    // version.
    {
      Field<const double,1> c_f1("CONST Test1",d1);
      Field<const double,2> c_f2("CONST Test2",d2);
      Field<const double,3> c_f3("CONST Test3",d3);
      Field<const double,4> c_f4("CONST Test4",d4);
      Field<const double,5> c_f5("CONST Test5",d5);
      Field<const double,6> c_f6("CONST Test6",d6);
      Field<const double,7> c_f7("CONST Test7",d7);

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
    // Kokkos static View accessors
    {
      // non-const view
      auto kva = a.get_static_view();
      kva(0,0) = 1.0;
      auto kvc = c.get_static_view();
      kvc(0,0) = MyTraits::FadType(1.0);
      // const view (view const, not const data)
      const auto const_kva = a.get_static_view();
      const_kva(0,0) = 1.0;
      const auto const_kvc = c.get_static_view();
      const_kvc(0,0) = MyTraits::FadType(1.0);
    }
    // Kokkos DynRankView accessors
    {
      // non-const view
      auto kva = a.get_view();
      kva(0,0) = 1.0;
      auto kvc = c.get_view();
      kvc(0,0) = MyTraits::FadType(1.0);
      // const view (view const, not const data)
      const Kokkos::DynRankView<double,PHX::Device> const_kva = a.get_view();
      const_kva(0,0) = 1.0;
      const auto const_kvc = c.get_view();
      const_kvc(0,0) = MyTraits::FadType(1.0);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for debug build array rank enforcement
    TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(f2.fieldTag())),PHX::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    ostringstream output;
    output << a;
    TEST_EQUALITY(output.str(), "Field<2>(100,4): Tag: density, double, DataLayout: NODE_SCALAR(100,4)");

  }

  TimeMonitor::summarize();
}
