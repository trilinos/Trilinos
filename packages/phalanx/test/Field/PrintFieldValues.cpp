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
#include "Phalanx_MDField.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_PrintValues.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

#include <sstream>
#include <algorithm>

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

// ***************************
// Static MDFields
// ***************************
TEUCHOS_UNIT_TEST(PrintValues, MDField_Static)
{
  using namespace PHX;
  
  int d = 2;
  std::string layout_name = "layout";
  std::ostringstream os;
  
  // Static View
  MDField<double,Cell> v1("v1",layout_name,d);
  MDField<double,Cell,Dim> v2("v2",layout_name,d,d);
  MDField<double,Cell,Dim,Dim> v3("v3",layout_name,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim> v4("v4",layout_name,d,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim> v5("v5",layout_name,d,d,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim,Dim> v6("v6",layout_name,d,d,d,d,d,d);
  MDField<double,Cell,Dim,Dim,Dim,Dim,Dim,Dim> v7("v7",layout_name,d,d,d,d,d,d,d);

  v1.deep_copy(1);
  v2.deep_copy(2);
  v3.deep_copy(3);
  v4.deep_copy(4);
  v5.deep_copy(5);
  v6.deep_copy(6);
  v7.deep_copy(7);

  // Counts the number of output lines
  auto checkOutput = [&os, &out,&success](auto v) {
                       char end_line = '\n';
                       std::string output = os.str();
                       size_t num_lines = std::count(output.begin(),output.end(),end_line);
                       TEST_EQUALITY(num_lines,v.size());
                       os.str("");
                     };
  
  TEST_EQUALITY(os.str(),"");
  PHX::printValues(v1,os,"v1");
  checkOutput(v1);
  PHX::printValues(v2,os,"v2");
  checkOutput(v2);
  PHX::printValues(v3,os,"v3");
  checkOutput(v3);
  PHX::printValues(v4,os,"v4");
  checkOutput(v4);
  PHX::printValues(v5,os,"v5");
  checkOutput(v5);
  PHX::printValues(v6,os,"v6");
  checkOutput(v6);
  PHX::printValues(v7,os,"v7");
  checkOutput(v7);
}

// ***************************
// Dynamic MDFields
// ***************************
TEUCHOS_UNIT_TEST(PrintValues, MDField_Dynamic)
{
  using namespace PHX;
  
  int d = 2;
  std::string layout_name = "layout";
  std::ostringstream os;
  
  // Static View
  MDField<double> v1("v1",layout_name,d);
  MDField<double> v2("v2",layout_name,d,d);
  MDField<double> v3("v3",layout_name,d,d,d);
  MDField<double> v4("v4",layout_name,d,d,d,d);
  MDField<double> v5("v5",layout_name,d,d,d,d,d);
  MDField<double> v6("v6",layout_name,d,d,d,d,d,d);
  MDField<double> v7("v7",layout_name,d,d,d,d,d,d,d);

  v1.deep_copy(1);
  v2.deep_copy(2);
  v3.deep_copy(3);
  v4.deep_copy(4);
  v5.deep_copy(5);
  v6.deep_copy(6);
  v7.deep_copy(7);

  // Counts the number of output lines
  auto checkOutput = [&os, &out,&success](auto v) {
                       char end_line = '\n';
                       std::string output = os.str();
                       size_t num_lines = std::count(output.begin(),output.end(),end_line);
                       TEST_EQUALITY(num_lines,v.size());
                       os.str("");
                     };
  
  TEST_EQUALITY(os.str(),"");
  PHX::printValues(v1,os,"v1");
  checkOutput(v1);
  PHX::printValues(v2,os,"v2");
  checkOutput(v2);
  PHX::printValues(v3,os,"v3");
  checkOutput(v3);
  PHX::printValues(v4,os,"v4");
  checkOutput(v4);
  PHX::printValues(v5,os,"v5");
  checkOutput(v5);
  PHX::printValues(v6,os,"v6");
  checkOutput(v6);
  PHX::printValues(v7,os,"v7");
  checkOutput(v7);
}

// ***************************
// Fields
// ***************************
TEUCHOS_UNIT_TEST(PrintValues, Field)
{
  using namespace PHX;
  
  int d = 2;
  std::string layout_name = "layout";
  std::ostringstream os;
  
  // Static View
  Field<double,1> v1("v1",layout_name,d);
  Field<double,2> v2("v2",layout_name,d,d);
  Field<double,3> v3("v3",layout_name,d,d,d);
  Field<double,4> v4("v4",layout_name,d,d,d,d);
  Field<double,5> v5("v5",layout_name,d,d,d,d,d);
  Field<double,6> v6("v6",layout_name,d,d,d,d,d,d);
  Field<double,7> v7("v7",layout_name,d,d,d,d,d,d,d);

  v1.deep_copy(1);
  v2.deep_copy(2);
  v3.deep_copy(3);
  v4.deep_copy(4);
  v5.deep_copy(5);
  v6.deep_copy(6);
  v7.deep_copy(7);

  // Counts the number of output lines
  auto checkOutput = [&os, &out,&success](auto v) {
                       char end_line = '\n';
                       std::string output = os.str();
                       size_t num_lines = std::count(output.begin(),output.end(),end_line);
                       TEST_EQUALITY(num_lines,v.size());
                       os.str("");
                     };
  
  TEST_EQUALITY(os.str(),"");
  PHX::printValues(v1,os,"v1");
  checkOutput(v1);
  PHX::printValues(v2,os,"v2");
  checkOutput(v2);
  PHX::printValues(v3,os,"v3");
  checkOutput(v3);
  PHX::printValues(v4,os,"v4");
  checkOutput(v4);
  PHX::printValues(v5,os,"v5");
  checkOutput(v5);
  PHX::printValues(v6,os,"v6");
  checkOutput(v6);
  PHX::printValues(v7,os,"v7");
  checkOutput(v7);
}

// ***************************
// Check output string
// ***************************
TEUCHOS_UNIT_TEST(PrintValues, CheckString)
{
  PHX::MDField<int,Cell,Dim> f("My Field","layout",2,2);
  f.deep_copy(6);
  std::ostringstream os;
  
  PHX::printValues(f,os,f.fieldTag().name());

  std::string gold_value = "My Field(0,0) = 6\nMy Field(0,1) = 6\nMy Field(1,0) = 6\nMy Field(1,1) = 6\n";
  TEST_EQUALITY(os.str(),gold_value);
}

// ***************************
// Make sure type checking is working
// ***************************
template<typename T>
bool constexpr check_is_mdfield(T& mdfield)
{return PHX::is_mdfield<T>::value;}

template<typename T>
bool constexpr check_is_field(T& field)
{return PHX::is_field<T>::value;}

template<typename T>
bool constexpr check_is_view(T& view)
{return Kokkos::is_view<T>::value;}

template<typename T>
bool constexpr check_is_dyn_rank_view(T& dyn_rank_view)
{return PHX::is_dyn_rank_view<T>::value;}

TEUCHOS_UNIT_TEST(PrintValues, PrintMDFieldTypeChecks)
{
  int d = 3;
  std::string layout_name = "layout";


  // Static MDField
  {
    PHX::MDField<double,Cell,Dim> mdf("MDField",layout_name,d,d);
    const PHX::MDField<double,Cell,Dim> mdf_const("Const MDField",layout_name,d,d);
    static_assert(check_is_mdfield(mdf));
    static_assert(check_is_mdfield(mdf_const));
    static_assert(!check_is_field(mdf));
    static_assert(!check_is_field(mdf_const));
    static_assert(!check_is_view(mdf));
    static_assert(!check_is_view(mdf_const));
    static_assert(!check_is_dyn_rank_view(mdf));
    static_assert(!check_is_dyn_rank_view(mdf_const));
  }
  
  // Runtime MDField
  {
    PHX::MDField<double> mdf("MDField",layout_name,d,d);
    const PHX::MDField<double> mdf_const("Const MDField",layout_name,d,d);
    static_assert(check_is_mdfield(mdf));
    static_assert(check_is_mdfield(mdf_const));
    static_assert(!check_is_field(mdf));
    static_assert(!check_is_field(mdf_const));
    static_assert(!check_is_view(mdf));
    static_assert(!check_is_view(mdf_const));
    static_assert(!check_is_dyn_rank_view(mdf));
    static_assert(!check_is_dyn_rank_view(mdf_const));
  }
  
  // Kokkos View
  {
    Kokkos::View<double**> v("View",d,d);
    const Kokkos::View<double**> v_const("Const View",d,d);
    static_assert(!check_is_mdfield(v));
    static_assert(!check_is_mdfield(v_const));
    static_assert(!check_is_field(v));
    static_assert(!check_is_field(v_const));
    static_assert(check_is_view(v));
    static_assert(check_is_view(v_const));
    static_assert(!check_is_dyn_rank_view(v));
    static_assert(!check_is_dyn_rank_view(v_const));
  }
  
  // Kokkos DynRankView
  {
    Kokkos::DynRankView<double> v("DynRankView",d,d);
    const Kokkos::DynRankView<double> v_const("Const DynRankView",d,d);
    static_assert(!check_is_mdfield(v));
    static_assert(!check_is_mdfield(v_const));
    static_assert(!check_is_field(v));
    static_assert(!check_is_field(v_const));
    static_assert(!check_is_view(v));
    static_assert(!check_is_view(v_const));
    static_assert(check_is_dyn_rank_view(v));
    static_assert(check_is_dyn_rank_view(v_const));
  }
}
