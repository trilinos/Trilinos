// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#include <limits> // for numeric_limits

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

// Functor to test accessors
template<typename F1,typename F2,typename F3,typename F4,typename F5,typename F6,typename F7,
	 typename CF1,typename CF2,typename CF3,typename CF4,typename CF5,typename CF6,typename CF7>
struct TestAssignmentFunctor {

  F1 f1_; F2 f2_; F3 f3_; F4 f4_; F5 f5_; F6 f6_; F7 f7_;
  CF1 cf1_; CF2 cf2_; CF3 cf3_; CF4 cf4_; CF5 cf5_; CF6 cf6_; CF7 cf7_;

  TestAssignmentFunctor(F1& f1,F2& f2,F3& f3,F4& f4,F5& f5, F6& f6, F7& f7,
			CF1& cf1,CF2& cf2,CF3& cf3,CF4& cf4,CF5& cf5, CF6& cf6, CF7& cf7)
    : f1_(f1),f2_(f2),f3_(f3),f4_(f4),f5_(f5),f6_(f6),f7_(f7),
      cf1_(cf1),cf2_(cf2),cf3_(cf3),cf4_(cf4),cf5_(cf5),cf6_(cf6),cf7_(cf7)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    using size_type = std::size_t;

    f1_(i) = cf1_(i);
    f1_.access(i) = cf1_.access(i);
    for (size_type j=0; j < f7_.dimension(1); ++j) {
      f2_(i,j) = cf2_(i,j);
      f2_(i,j) = cf2_.access(i,j);
      for (size_type k=0; k < f7_.dimension(2); ++k) {
	f3_(i,j,k) = cf3_(i,j,k);
	f3_(i,j,k) = cf3_.access(i,j,k);
	for (size_type l=0; l < f7_.dimension(3); ++l) {
	  f4_(i,j,k,l) = cf4_(i,j,k,l);
	  f4_(i,j,k,l) = cf4_.access(i,j,k,l);
	  for (size_type m=0; m < f7_.dimension(4); ++m) {
	    f5_(i,j,k,l,m) = cf5_(i,j,k,l,m);
	    f5_(i,j,k,l,m) = cf5_.access(i,j,k,l,m);
	    for (size_type n=0; n < f7_.dimension(5); ++n) {
	      f6_(i,j,k,l,m,n) = cf6_(i,j,k,l,m,n);
	      f6_(i,j,k,l,m,n) = cf6_.access(i,j,k,l,m,n);
	      for (size_type o=0; o < f7_.dimension(6); ++o) {
		f7_(i,j,k,l,m,n,o) = cf7_(i,j,k,l,m,n,o);
		f7_(i,j,k,l,m,n,o) = cf7_.access(i,j,k,l,m,n,o);
	      }
	    }
	  }
	}
      }
    }
  }
};

TEUCHOS_UNIT_TEST(mdfield, CompileTimeChecked)
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

    // Test create, set and get for RCP<FieldTag>
    {
      RCP<PHX::FieldTag> t1 = rcp(new PHX::Tag<double>("test",node_scalar));
      MDField<double,Cell,Node> f1(t1); // rcp(tag) ctor
      MDField<double,Cell,Node> f2;
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
      MDField<const double,Cell,Point> c_field1(nonconst_tag);

      // Create a non-const field from a const field tag
      MDField<double,Cell,Point> nc_field1(const_tag);

      // Create a non-const field from a non-const field tag
      MDField<double,Cell,Point> nc_field2(nonconst_tag);

      // Create a const field from a const field tag
      MDField<const double,Cell,Point> c_field2(const_tag);
    }

    // Copy constructor from const/non-const MDFields
    {
      RCP<DataLayout> ctor_dl_p  = rcp(new MDALayout<Cell,Point>(10,4));
      MDField<double,Cell,Point> ctor_nonconst_p("ctor_nonconst_p",ctor_dl_p);
      MDField<const double,Cell,Point> ctor_const_p("ctor_const_p",ctor_dl_p);

      MDField<double,Cell,Point> cc1(ctor_nonconst_p);       // non-const from non-const
      MDField<const double,Cell,Point> cc2(ctor_nonconst_p); // const from non-const
      MDField<const double,Cell,Point> cc3(ctor_const_p);    // const from const

      // NOTE: we allow the tag template types to be DIFFERENT as long
      // as the rank is the same! A field might use the "Point" DimTag
      // but another evaluator might reference the same field using
      // QuadraturePoint DimTag.
      RCP<DataLayout> ctor_dl_qp = rcp(new MDALayout<Cell,Quadrature>(10,4));
      MDField<double,Cell,Quadrature> ctor_nonconst_qp("ctor_nonconst",ctor_dl_qp);
      MDField<const double,Cell,Quadrature> ctor_const_qp("ctor_const_qp",ctor_dl_qp);

      // Repeat test above but with different tags for Quadrature --> Point
      MDField<double,Cell,Point> cc4(ctor_nonconst_qp);       // non-const from non-const
      MDField<const double,Cell,Point> cc5(ctor_nonconst_qp); // const from non-const
      MDField<const double,Cell,Point> cc6(ctor_const_qp);    // const from const

      // While we have these objects, lets test the assignment operator as well
      MDField<double,Cell,Point> cc7(ctor_nonconst_p);         // non-const from non-const
      MDField<const double,Cell,Point> cc8(ctor_nonconst_p);   // const from non-const
      MDField<const double,Cell,Point> cc9(ctor_const_p);      // const from const
      MDField<double,Cell,Point> cc10(ctor_nonconst_qp);       // non-const from non-const
      MDField<const double,Cell,Point> cc11(ctor_nonconst_qp); // const from non-const
      MDField<const double,Cell,Point> cc12(ctor_const_qp);    // const from const

      cc7 = ctor_nonconst_p;
      cc8 = ctor_nonconst_p;
      cc9 = ctor_const_p;
      cc10 = ctor_nonconst_qp;
      cc11 = ctor_nonconst_qp;
      cc12 = ctor_const_qp;
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
    std::any a_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(a.fieldTag());
    std::any b_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(b.fieldTag());
    std::any c_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(c.fieldTag(),ddims);
    std::any d_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(d.fieldTag(),ddims);
    std::any e_mem = PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(e.fieldTag());
    std::any f_mem = PHX::KokkosViewFactory<MyTraits::FadType,typename PHX::DevLayout<MyTraits::FadType>::type,PHX::Device>::buildView(f.fieldTag(),ddims);

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
    cout << "Testing operator()(...) accessors...";

    RCP<DataLayout> d1 = rcp(new MDALayout<Cell>(num_cells));
    RCP<DataLayout> d2 = rcp(new MDALayout<Cell,Dim>(num_cells,1));
    RCP<DataLayout> d3 = rcp(new MDALayout<Cell,Dim,Dim>(num_cells,1,2));
    RCP<DataLayout> d4 = rcp(new MDALayout<Cell,Dim,Dim,Dim>(num_cells,1,2,3));
    RCP<DataLayout> d5 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4));
    RCP<DataLayout> d6 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5));
    RCP<DataLayout> d7 = rcp(new MDALayout<Cell,Dim,Dim,Dim,Dim,Dim,Dim>(num_cells,1,2,3,4,5,6));

    // Use Unmanaged allocator
    MDField<double,Cell> f1 = PHX::allocateUnmanagedMDField<double,Cell>("Test1",d1);
    MDField<double,Cell,Dim> f2 = PHX::allocateUnmanagedMDField<double,Cell,Dim>("Test2",d2);
    MDField<double,Cell,Dim,Dim> f3 = PHX::allocateUnmanagedMDField<double,Cell,Dim,Dim>("Test3",d3);
    MDField<double,Cell,Dim,Dim,Dim> f4 = PHX::allocateUnmanagedMDField<double,Cell,Dim,Dim,Dim>("Test4",d4);
    MDField<double,Cell,Dim,Dim,Dim,Dim> f5 = PHX::allocateUnmanagedMDField<double,Cell,Dim,Dim,Dim,Dim>("Test5",d5);
    MDField<double,Cell,Dim,Dim,Dim,Dim,Dim> f6 = PHX::allocateUnmanagedMDField<double,Cell,Dim,Dim,Dim,Dim,Dim>("Test6",d6);
    MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim> f7 = PHX::allocateUnmanagedMDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim>("Test7",d7);

    // Test const/ non-const versions
    const MDField<double,Cell>& cf1 = f1;
    const MDField<double,Cell,Dim>& cf2 = f2;
    const MDField<double,Cell,Dim,Dim>& cf3 = f3;
    const MDField<double,Cell,Dim,Dim,Dim>& cf4 = f4;
    const MDField<double,Cell,Dim,Dim,Dim,Dim>& cf5 = f5;
    const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim>& cf6 = f6;
    const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim>& cf7 = f7;

    auto func = TestAssignmentFunctor<MDField<double,Cell>,
				      MDField<double,Cell,Dim>,
				      MDField<double,Cell,Dim,Dim>,
				      MDField<double,Cell,Dim,Dim,Dim>,
				      MDField<double,Cell,Dim,Dim,Dim,Dim>,
				      MDField<double,Cell,Dim,Dim,Dim,Dim,Dim>,
				      MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim>,
				      const MDField<double,Cell>,
				      const MDField<double,Cell,Dim>,
				      const MDField<double,Cell,Dim,Dim>,
				      const MDField<double,Cell,Dim,Dim,Dim>,
				      const MDField<double,Cell,Dim,Dim,Dim,Dim>,
				      const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim>,
				      const MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim>
				      >(f1,f2,f3,f4,f5,f6,f7,
					cf1,cf2,cf3,cf4,cf5,cf6,cf7);

    Kokkos::parallel_for("TestAssignmentFunctor",
			 Kokkos::RangePolicy<PHX::Device>(0,f7.extent(0)),
			 func);

    cout << "passed!" << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for const mdfield assignment from non-const factory
    // std::any. The field manager always stores the non-const
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
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){kva(0,0) = 1.0;});
      auto kvc = c.get_static_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){kvc(0,0) = MyTraits::FadType(1.0);});
      // const view (view const, not const data)
      const auto const_kva = a.get_static_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){const_kva(0,0) = 1.0;});
      const auto const_kvc = c.get_static_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){const_kvc(0,0) = MyTraits::FadType(1.0);});
    }
    // Kokkos DynRankView accessors
    {
      // non-const view
      auto kva = a.get_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){kva(0,0) = 1.0;});
      auto kvc = c.get_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){kvc(0,0) = MyTraits::FadType(1.0);});
      // const view (view const, not const data)
      const auto const_kva = a.get_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){const_kva(0,0) = 1.0;});
      const auto const_kvc = c.get_view();
      Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),KOKKOS_LAMBDA(const int ){const_kvc(0,0) = MyTraits::FadType(1.0);});
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Kokkos cast View accessor
    {
      decltype(a.get_static_view()) av;
      av = a;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PHX as_view accessor
    {
      Kokkos::deep_copy(a.get_static_view(), 5.);
      auto a_h = Kokkos::create_mirror_view(as_view(a));
      Kokkos::deep_copy(a_h, as_view(a));
      TEST_EQUALITY(a_h(0,0), 5);
      
      auto a_v = a.get_static_view();
      auto a_v_h = Kokkos::create_mirror_view(as_view(a_v));
      Kokkos::deep_copy(a_v_h, as_view(a_v));
      TEST_EQUALITY(a_v_h(0,0), 5);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for debug build array rank enforcement
    TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,typename PHX::DevLayout<double>::type,PHX::Device>::buildView(f2.fieldTag())),std::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    ostringstream output;
    output << a;
    // Disable below - name mangling not handled correctly on all platforms.
    // TEST_EQUALITY(output.str(), "MDField<C,N>(100,4): Tag: density, double, DataLayout: <C,N>(100,4)");
  }

  TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(mdfield, UnsafeCtor)
{
  using namespace PHX;

  int c = 5;
  int d = 2;
  std::string layout_name = "l";
  MDField<double,Cell> d1("d1",layout_name,c);
  MDField<double,Cell,Dim> d2("d2",layout_name,c,d);
  MDField<double,Cell,Dim,Dim> d3("d3",layout_name,c,d,d);
  MDField<double,Cell,Dim,Dim,Dim> d4("d4",layout_name,c,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim> d5("d5",layout_name,c,d,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim,Dim> d6("d6",layout_name,c,d,d,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim> d7("d7",layout_name,c,d,d,d,d,d,d);

  int derivative_dimension = 3;
  using FadType = MyTraits::FadType;
  MDField<FadType,Cell> f1("f1",layout_name,c,derivative_dimension);
  MDField<FadType,Cell,Dim> f2("f2",layout_name,c,d,derivative_dimension);
  MDField<FadType,Cell,Dim,Dim> f3("f3",layout_name,c,d,d,derivative_dimension);
  MDField<FadType,Cell,Dim,Dim,Dim> f4("f4",layout_name,c,d,d,d,derivative_dimension);
  MDField<FadType,Cell,Dim,Dim,Dim,Dim> f5("f5",layout_name,c,d,d,d,d,derivative_dimension);
  MDField<FadType,Cell,Dim,Dim,Dim,Dim,Dim> f6("f6",layout_name,c,d,d,d,d,d,derivative_dimension);
  MDField<FadType,Cell,Dim,Dim,Dim,Dim,Dim,Dim> f7("f7",layout_name,c,d,d,d,d,d,d,derivative_dimension);

  TEST_EQUALITY(d1.size(),std::size_t(c));
  TEST_EQUALITY(d2.size(),std::size_t(c*d));
  TEST_EQUALITY(d3.size(),std::size_t(c*d*d));
  TEST_EQUALITY(d4.size(),std::size_t(c*d*d*d));
  TEST_EQUALITY(d5.size(),std::size_t(c*d*d*d*d));
  TEST_EQUALITY(d6.size(),std::size_t(c*d*d*d*d*d));
  TEST_EQUALITY(d7.size(),std::size_t(c*d*d*d*d*d*d));
  TEST_EQUALITY(d1.rank(),std::size_t(1));
  TEST_EQUALITY(d2.rank(),std::size_t(2));
  TEST_EQUALITY(d3.rank(),std::size_t(3));
  TEST_EQUALITY(d4.rank(),std::size_t(4));
  TEST_EQUALITY(d5.rank(),std::size_t(5));
  TEST_EQUALITY(d6.rank(),std::size_t(6));
  TEST_EQUALITY(d7.rank(),std::size_t(7));
  TEST_EQUALITY(d1.size(),std::size_t(c));

  TEST_EQUALITY(f2.size(),std::size_t(c*d));
  TEST_EQUALITY(f3.size(),std::size_t(c*d*d));
  TEST_EQUALITY(f4.size(),std::size_t(c*d*d*d));
  TEST_EQUALITY(f5.size(),std::size_t(c*d*d*d*d));
  TEST_EQUALITY(f6.size(),std::size_t(c*d*d*d*d*d));
  TEST_EQUALITY(f7.size(),std::size_t(c*d*d*d*d*d*d));
  TEST_EQUALITY(f1.rank(),std::size_t(1));
  TEST_EQUALITY(f2.rank(),std::size_t(2));
  TEST_EQUALITY(f3.rank(),std::size_t(3));
  TEST_EQUALITY(f4.rank(),std::size_t(4));
  TEST_EQUALITY(f5.rank(),std::size_t(5));
  TEST_EQUALITY(f6.rank(),std::size_t(6));
  TEST_EQUALITY(f7.rank(),std::size_t(7));
}

TEUCHOS_UNIT_TEST(mdfield, releaseFieldData)
{
  using namespace PHX;

  const size_t c = 10;
  const size_t p = 9;
  const size_t d = 2;
  MDField<double,Cell,Point,Dim> a("a","",c,p,d);

  TEST_EQUALITY(a.extent(0),c);
  TEST_EQUALITY(a.extent(1),p);
  TEST_EQUALITY(a.extent(2),d);

  a.releaseFieldData();

  TEST_EQUALITY(a.extent(0),0);
  TEST_EQUALITY(a.extent(1),0);
  TEST_EQUALITY(a.extent(2),0);
}

TEUCHOS_UNIT_TEST(mdfield, AssignCompileTimeToRuntime)
{
  using namespace PHX;

  const size_t c = 3;
  const size_t p = 4;
  const size_t d = 2;
  MDField<double,Cell,Point,Dim> a("a","",c,p,d);

  MDField<double> a_dyn = a;

  Kokkos::deep_copy(a.get_static_view(),3.0);

  auto a_dyn_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),a_dyn.get_view());
  const auto tol = 100.0 * std::numeric_limits<double>::epsilon();

  for (size_t i = 0; i < a_dyn_host.extent(0); ++i) {
    for (size_t j = 0; j < a_dyn_host.extent(1); ++j) {
      for (size_t k = 0; k < a_dyn_host.extent(2); ++k) {
        TEST_FLOATING_EQUALITY(a_dyn_host(i,j,k),3.0,tol);
      }
    }
  }
}
