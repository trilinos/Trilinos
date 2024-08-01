// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

// Functor to test accessors
template<typename F,typename CF>
struct TestAssignmentFunctor {

  F f1_, f2_, f3_, f4_, f5_, f6_, f7_;
  CF cf1_, cf2_, cf3_, cf4_, cf5_, cf6_, cf7_;

  TestAssignmentFunctor(F& f1,F& f2,F& f3,F& f4,F& f5, F& f6, F& f7,
			CF& cf1,CF& cf2,CF& cf3,CF& cf4,CF& cf5, CF& cf6, CF& cf7)
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

#define PHX_TEST_FLOAT_EQUAL(A,B,tol,failed) if (std::fabs(A-B) > tol) failed += 1;

template<typename FieldType,typename ScalarType>
struct BracketAssignmentTest {
  FieldType f1_,f2_,f3_,f4_,f5_,f6_;
  BracketAssignmentTest(FieldType f1,FieldType f2,FieldType f3,FieldType f4,FieldType f5,FieldType f6)
    : f1_(f1),f2_(f2),f3_(f3),f4_(f4),f5_(f5),f6_(f6){}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int ) const {
    f1_[f1_.size()-1] = ScalarType(2.0);
    f2_[f2_.size()-1] = ScalarType(3.0);
    f3_[f3_.size()-1] = ScalarType(4.0);
    f4_[f4_.size()-1] = ScalarType(5.0);
    f5_[f5_.size()-1] = ScalarType(6.0);
    f6_[f6_.size()-1] = ScalarType(7.0);
  }
};

template<typename FieldType,typename ScalarType>
struct BracketAssignmentCheck {
  FieldType f1_,f2_,f3_,f4_,f5_,f6_;
  const FieldType cf1_,cf2_,cf3_,cf4_,cf5_,cf6_;
  BracketAssignmentCheck(FieldType f1,FieldType f2,FieldType f3,FieldType f4,FieldType f5,FieldType f6,
			 const FieldType cf1,const FieldType cf2,const FieldType cf3,const FieldType cf4,const FieldType cf5,const FieldType cf6)
    : f1_(f1),f2_(f2),f3_(f3),f4_(f4),f5_(f5),f6_(f6),
      cf1_(cf1),cf2_(cf2),cf3_(cf3),cf4_(cf4),cf5_(cf5),cf6_(cf6){}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int , int& lsum) const {
    const double tol = 1.0e-10;
    PHX_TEST_FLOAT_EQUAL(f1_[f1_.size()-1], cf1_[f1_.size()-1], tol, lsum);
    PHX_TEST_FLOAT_EQUAL(f2_[f2_.size()-1], cf2_[f2_.size()-1], tol, lsum);
    PHX_TEST_FLOAT_EQUAL(f3_[f3_.size()-1], cf3_[f3_.size()-1], tol, lsum);
    PHX_TEST_FLOAT_EQUAL(f4_[f4_.size()-1], cf4_[f4_.size()-1], tol, lsum);
    PHX_TEST_FLOAT_EQUAL(f5_[f5_.size()-1], cf5_[f5_.size()-1], tol, lsum);
    PHX_TEST_FLOAT_EQUAL(f6_[f6_.size()-1], cf6_[f6_.size()-1], tol, lsum);
  }
};

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

    // Use unmanaged allocator
    MDField<double> f1 = PHX::allocateUnmanagedMDField<double>("Test1",d1);
    MDField<double> f2 = PHX::allocateUnmanagedMDField<double>("Test2",d2);
    MDField<double> f3 = PHX::allocateUnmanagedMDField<double>("Test3",d3);
    MDField<double> f4 = PHX::allocateUnmanagedMDField<double>("Test4",d4);
    MDField<double> f5 = PHX::allocateUnmanagedMDField<double>("Test5",d5);
    MDField<double> f6 = PHX::allocateUnmanagedMDField<double>("Test6",d6);
    MDField<double> f7 = PHX::allocateUnmanagedMDField<double>("Test7",d7);

    // Test const/ non-const versions
    const MDField<double>& cf1 = f1;
    const MDField<double>& cf2 = f2;
    const MDField<double>& cf3 = f3;
    const MDField<double>& cf4 = f4;
    const MDField<double>& cf5 = f5;
    const MDField<double>& cf6 = f6;
    const MDField<double>& cf7 = f7;

    auto func = TestAssignmentFunctor<MDField<double>,const MDField<double>>
      (f1,f2,f3,f4,f5,f6,f7,cf1,cf2,cf3,cf4,cf5,cf6,cf7);

    Kokkos::parallel_for("TestAssignmentFunctor",
    			 Kokkos::RangePolicy<PHX::Device>(0,f7.extent(0)),
    			 func);
    PHX::Device().fence();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // check for const mdfield assignment from non-const factory
    // std::any. The field manager always stores the non-const
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

    // Test const/ non-const versions
    const MDField<MyTraits::FadType>& cf1_fad = f1_fad;
    const MDField<MyTraits::FadType>& cf2_fad = f2_fad;
    const MDField<MyTraits::FadType>& cf3_fad = f3_fad;
    const MDField<MyTraits::FadType>& cf4_fad = f4_fad;
    const MDField<MyTraits::FadType>& cf5_fad = f5_fad;
    const MDField<MyTraits::FadType>& cf6_fad = f6_fad;
    const MDField<MyTraits::FadType>& cf7_fad = f7_fad;

    auto func_fad = TestAssignmentFunctor<MDField<MyTraits::FadType>,const MDField<MyTraits::FadType>>
      (f1_fad,f2_fad,f3_fad,f4_fad,f5_fad,f6_fad,f7_fad,
       cf1_fad,cf2_fad,cf3_fad,cf4_fad,cf5_fad,cf6_fad,cf7_fad);

    Kokkos::parallel_for("TestAssignmentFunctor",
    			 Kokkos::RangePolicy<PHX::Device>(0,f7_fad.extent(0)),
    			 func_fad);
    PHX::Device().fence();

    out << "passed!" << endl;

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

    Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),BracketAssignmentTest<MDField<double>,double>(f1,f2,f3,f4,f5,f6));

    int num_failed = 0;
    Kokkos::parallel_reduce("t1",Kokkos::RangePolicy<PHX::Device>(0,1),BracketAssignmentCheck<MDField<double>,double>(f1,f2,f3,f4,f5,f6,cf1,cf2,cf3,cf4,cf5,cf6),num_failed);
    TEST_EQUALITY(num_failed,0);

    // fad checking
    Kokkos::parallel_for("t1",Kokkos::RangePolicy<PHX::Device>(0,1),BracketAssignmentTest<MDField<MyTraits::FadType>,MyTraits::FadType>(f1_fad,f2_fad,f3_fad,f4_fad,f5_fad,f6_fad));

    out << "passed!" << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Do NOT check for array rank enforcement. DynRank MDField sets
    // the rank at runtime and we allow it to be changed!
    //TEST_THROW(f1.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(f2.fieldTag())),std::bad_any_cast);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // kokkos view accessors
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
    // PHX as_view accessor
    {
      Kokkos::deep_copy(a.get_view(), 5.);
      auto a_h = Kokkos::create_mirror_view(as_view(a));
      Kokkos::deep_copy(a_h, as_view(a));
      TEST_EQUALITY(a_h(0,0), 5);
      
      auto a_v = a.get_static_view();
      auto a_v_h = Kokkos::create_mirror_view(as_view(a_v));
      Kokkos::deep_copy(a_v_h, as_view(a_v));
      TEST_EQUALITY(a_v_h(0,0), 5);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ostream
    out << "Testing operator<<()...";
    ostringstream output;
    output << a;
    // Disable below - name mangling not handled correctly on all platforms.
    // TEST_EQUALITY(output.str(),"MDField(100,4): Tag: density, double, DataLayout: <C,N>(100,4)");
    out << "passed!" << endl;
    out << output.str() << endl;
  }

  TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(mdfield, UnsafeCtor)
{
  using namespace PHX;

  int c = 5;
  int d = 2;
  std::string layout_name = "l";
  MDField<double> d1("d1",layout_name,c);
  MDField<double> d2("d2",layout_name,c,d);
  MDField<double> d3("d3",layout_name,c,d,d);
  MDField<double> d4("d4",layout_name,c,d,d,d);
  MDField<double> d5("d5",layout_name,c,d,d,d,d);
  MDField<double> d6("d6",layout_name,c,d,d,d,d,d);
  MDField<double> d7("d7",layout_name,c,d,d,d,d,d,d);

  int derivative_dimension = 3;
  using FadType = MyTraits::FadType;
  MDField<FadType> f1("f1",layout_name,c,derivative_dimension);
  MDField<FadType> f2("f2",layout_name,c,d,derivative_dimension);
  MDField<FadType> f3("f3",layout_name,c,d,d,derivative_dimension);
  MDField<FadType> f4("f4",layout_name,c,d,d,d,derivative_dimension);
  MDField<FadType> f5("f5",layout_name,c,d,d,d,d,derivative_dimension);
  MDField<FadType> f6("f6",layout_name,c,d,d,d,d,d,derivative_dimension);
  MDField<FadType> f7("f7",layout_name,c,d,d,d,d,d,d,derivative_dimension);

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
  MDField<double> a("a","",c,p,d);

  TEST_EQUALITY(a.extent(0),c);
  TEST_EQUALITY(a.extent(1),p);
  TEST_EQUALITY(a.extent(2),d);

  a.releaseFieldData();

  TEST_EQUALITY(a.extent(0),0);
  TEST_EQUALITY(a.extent(1),0);
  TEST_EQUALITY(a.extent(2),0);
}
