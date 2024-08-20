// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_KokkosView_AllocationSize.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_Field.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"

PHX_EXTENT(D)

namespace phalanx_test {

  using real = double;
  using dfad = Sacado::Fad::DFad<double>;
  using elrcdfad = Sacado::ELRCacheFad::DFad<double>;
  using slfad = Sacado::Fad::SLFad<double,32>;
  using sfad = Sacado::Fad::SFad<double,32>;

  // ********************************************************
  template<typename ScalarT>
  void test_span(std::ostream& out,bool success)  {
    using namespace std;
    using namespace Teuchos;
    using namespace PHX;
    
    RCP<DataLayout> d1 = rcp(new MDALayout<D>(2));
    RCP<DataLayout> d2 = rcp(new MDALayout<D,D>(2,2));
    RCP<DataLayout> d3 = rcp(new MDALayout<D,D,D>(2,2,2));
    RCP<DataLayout> d4 = rcp(new MDALayout<D,D,D,D>(2,2,2,2));
    RCP<DataLayout> d5 = rcp(new MDALayout<D,D,D,D,D>(2,2,2,2,2));
    RCP<DataLayout> d6 = rcp(new MDALayout<D,D,D,D,D,D>(2,2,2,2,2,2));
    RCP<DataLayout> d7 = rcp(new MDALayout<D,D,D,D,D,D,D>(2,2,2,2,2,2,2));
    std::vector<PHX::index_size_type> derivative_dims;
    derivative_dims.push_back(10);
    
    using mem = PHX::KokkosViewFactory<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>;
    std::size_t alloc_size = 0;
    
    PHX::MDField<ScalarT,D> a1("a1",d1); // static rank
    PHX::MDField<ScalarT,D,D> a2("a2",d2); // static rank
    PHX::MDField<ScalarT,D,D,D> a3("a3",d3); // static rank
    PHX::MDField<ScalarT,D,D,D,D> a4("a4",d4); // static rank
    PHX::MDField<ScalarT,D,D,D,D,D> a5("a5",d5); // static rank
    PHX::MDField<ScalarT,D,D,D,D,D,D> a6("a6",d6); // static rank
    PHX::MDField<ScalarT,D,D,D,D,D,D,D> a7("a7",d7); // static rank

    a1.setFieldData(mem::buildView(a1.fieldTag(),derivative_dims));
    a2.setFieldData(mem::buildView(a2.fieldTag(),derivative_dims));
    a3.setFieldData(mem::buildView(a3.fieldTag(),derivative_dims));
    a4.setFieldData(mem::buildView(a4.fieldTag(),derivative_dims));
    a5.setFieldData(mem::buildView(a5.fieldTag(),derivative_dims));
    a6.setFieldData(mem::buildView(a6.fieldTag(),derivative_dims));
    a7.setFieldData(mem::buildView(a7.fieldTag(),derivative_dims));
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a1.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a1.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a2.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a2.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a3.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a3.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a4.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a4.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a5.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a5.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a6.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a6.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(a7.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,a7.span());
    
    PHX::MDField<ScalarT> b1("b",d1); // dynamic rank
    PHX::MDField<ScalarT> b2("b",d2); // dynamic rank
    PHX::MDField<ScalarT> b3("b",d3); // dynamic rank
    PHX::MDField<ScalarT> b4("b",d4); // dynamic rank
    PHX::MDField<ScalarT> b5("b",d5); // dynamic rank
    PHX::MDField<ScalarT> b6("b",d6); // dynamic rank
    PHX::MDField<ScalarT> b7("b",d7); // dynamic rank
    b1.setFieldData(mem::buildView(b1.fieldTag(),derivative_dims));
    b2.setFieldData(mem::buildView(b2.fieldTag(),derivative_dims));
    b3.setFieldData(mem::buildView(b3.fieldTag(),derivative_dims));
    b4.setFieldData(mem::buildView(b4.fieldTag(),derivative_dims));
    b5.setFieldData(mem::buildView(b5.fieldTag(),derivative_dims));
    b6.setFieldData(mem::buildView(b6.fieldTag(),derivative_dims));
    b7.setFieldData(mem::buildView(b7.fieldTag(),derivative_dims));
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b1.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b1.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b2.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b2.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b3.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b3.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b4.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b4.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b5.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b5.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b6.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b6.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(b7.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,b7.span());
    
    PHX::Field<ScalarT,1> c1("c1",d1);
    PHX::Field<ScalarT,2> c2("c2",d2);
    PHX::Field<ScalarT,3> c3("c3",d3);
    PHX::Field<ScalarT,4> c4("c4",d4);
    PHX::Field<ScalarT,5> c5("c5",d5);
    PHX::Field<ScalarT,6> c6("c6",d6);
    PHX::Field<ScalarT,7> c7("c7",d7);
    c1.setFieldData(mem::buildView(c1.fieldTag(),derivative_dims));
    c2.setFieldData(mem::buildView(c2.fieldTag(),derivative_dims));
    c3.setFieldData(mem::buildView(c3.fieldTag(),derivative_dims));
    c4.setFieldData(mem::buildView(c4.fieldTag(),derivative_dims));
    c5.setFieldData(mem::buildView(c5.fieldTag(),derivative_dims));
    c6.setFieldData(mem::buildView(c6.fieldTag(),derivative_dims));
    c7.setFieldData(mem::buildView(c7.fieldTag(),derivative_dims));
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c1.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c1.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c2.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c2.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c3.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c3.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c4.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c4.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c5.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c5.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c6.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c6.span());
    alloc_size = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(c7.fieldTag(),derivative_dims);
    TEST_EQUALITY(alloc_size,c7.span());
  }
  
  // ********************************************************
  TEUCHOS_UNIT_TEST(ViewAllocationSize, hidden_dimensions)
  {

    static_assert(!PHX::requires_dynamic_hidden_dimension<real>::value,
                  "ERROR: double does NOT require a dynamic hidden dimension");
    static_assert(!PHX::requires_dynamic_hidden_dimension<sfad>::value,
                  "ERROR: sfad does NOT require a dynamic hidden dimension");

    static_assert(PHX::requires_dynamic_hidden_dimension<dfad>::value,
                  "ERROR: dfad requires a dynamic hidden dimension");
    static_assert(PHX::requires_dynamic_hidden_dimension<elrcdfad>::value,
                  "ERROR: elrcdfad requires a dynamic hidden dimension");
    static_assert(PHX::requires_dynamic_hidden_dimension<slfad>::value,
                  "ERROR: slfad requires a dynamic hidden dimension");
  }

  TEUCHOS_UNIT_TEST(ViewAllocationSize, type_real)
  {
    test_span<real>(out, success);  
  }

  TEUCHOS_UNIT_TEST(ViewAllocationSize, type_dfad)
  {
    test_span<dfad>(out, success);
  }

  TEUCHOS_UNIT_TEST(ViewAllocationSize, type_elrcdfad)
  {
    test_span<elrcdfad>(out, success);
  }

  TEUCHOS_UNIT_TEST(ViewAllocationSize, type_slfad)
  {
    test_span<slfad>(out, success);
  }

  TEUCHOS_UNIT_TEST(ViewAllocationSize, type_sfad)
  {
    test_span<sfad>(out, success);
  }
  
} // namespace phalanx_test
