// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "KokkosExp_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ZeroSensitivities.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(zero_sensitivites, scalar)
  {

    panzer::Traits::RealType x = 1.0;
    TEST_FLOATING_EQUALITY(x,1.0,1e-12);
    panzer::zeroSensitivities(x);
    TEST_FLOATING_EQUALITY(x,1.0,1e-12);

    panzer::Traits::FadType y;
    y.val() = 1.0;
    y.resize(2);
    y.fastAccessDx(0) = 2.0;
    y.fastAccessDx(1) = 3.0;
    TEST_FLOATING_EQUALITY(y.val(),1.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(0),2.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(1),3.0,1e-12);
    panzer::zeroSensitivities(y);
    TEST_FLOATING_EQUALITY(y.val(),1.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(0),0.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(1),0.0,1e-12);
  }


  typedef PHX::index_size_type size_type;
  
  template <typename Scalar,typename Device,typename Array>
  class ComputeA {
    Array a_;
  public:
    typedef PHX::Device execution_space;
    
    ComputeA(Array& a)
      : a_(a)
    {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const size_type c) const
    {
      const size_type num_pts = a_.extent(1);
      for (size_type pt = 0; pt < num_pts; ++pt) {
	panzer::zeroSensitivities(a_(c,pt));
      }
    }
  };

  template <typename Device,typename Array>
  class ComputeB {
    Array a_;
  public:
    typedef PHX::Device execution_space;
    
    ComputeB(Array& a)
      : a_(a)
    {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const size_type c) const
    {
      const size_type num_pts = a_.extent(1);
      for (size_type pt = 0; pt < num_pts; ++pt) {
	panzer::zeroSensitivities(a_(c,pt));
      }
    }
  };

  TEUCHOS_UNIT_TEST(zero_sensitivites, mdfield)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using PHX::MDALayout;
    using PHX::MDField;
    using panzer::Cell;
    using panzer::BASIS;


    RCP<MDALayout<Cell,BASIS> > dl = rcp(new MDALayout<Cell,BASIS>(2,3));
    panzer::Traits::RealType tolerance = Teuchos::ScalarTraits<panzer::Traits::RealType>::eps()*100.0;

    // ***********************
    // RealType
    // ***********************
    {
      PHX::MDField<panzer::Traits::RealType,Cell,BASIS> a("rho",dl);
      a.setFieldData(PHX::KokkosViewFactory<panzer::Traits::RealType,typename PHX::DevLayout<panzer::Traits::RealType>::type,PHX::Device>::buildView(a.fieldTag()));
      
      // initialize
      a.deep_copy(2.0);
      // for (int cell = 0; cell < a.extent_int(0); ++cell) {
      //   for (int pt=0; pt < a.extent_int(1); ++pt) {
      //     a(cell,pt) = 2.0;
      //   }
      // }
      
      // Compute
      Kokkos::parallel_for(a.extent_int(0),ComputeA<panzer::Traits::RealType,PHX::Device,PHX::MDField<panzer::Traits::RealType,Cell,BASIS> > (a));
      typename PHX::Device().fence();

      auto a_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),a.get_static_view());
      
      // Check
      for (int cell = 0; cell < a.extent_int(0); ++cell)
	for (int pt=0; pt < a.extent_int(1); ++pt)
	  TEST_FLOATING_EQUALITY(a_host(cell,pt),2.0,tolerance);
    }

    // ***********************
    // FadType
    // ***********************
    {
      PHX::MDField<panzer::Traits::FadType,Cell,BASIS> a("rho",dl);
      std::vector<size_type> derivative_dimension;
      derivative_dimension.push_back(2);
      a.setFieldData(PHX::KokkosViewFactory<panzer::Traits::FadType,typename PHX::DevLayout<panzer::Traits::FadType>::type,PHX::Device>::buildView(a.fieldTag(),derivative_dimension));
      
      // Initialize (use raw kokkos to avoid compiler warning about MDField dtor with host function calls).
      auto a_dev = a.get_static_view();
      Kokkos::parallel_for("initialize zero sensitivity vec",a.extent(0),KOKKOS_LAMBDA (const int cell) {
          //for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a_dev.extent_int(1); ++pt) {
	  a_dev(cell,pt) = 1.0;
	  a_dev(cell,pt).fastAccessDx(0) = 2.0;
	  a_dev(cell,pt).fastAccessDx(1) = 3.0;
	}
      });

      // Check initial values
      auto a_host = Kokkos::create_mirror_view(a.get_static_view());
      Kokkos::deep_copy(a_host,a.get_static_view());
      for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a.extent_int(1); ++pt) {
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).val(),1.0,1e-12);
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).fastAccessDx(0),2.0,1e-12);
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).fastAccessDx(1),3.0,1e-12);
	}
      }

      // Compute
      Kokkos::parallel_for(a.extent(0),ComputeB<PHX::Device,PHX::MDField<panzer::Traits::FadType,Cell,BASIS> > (a));
      typename PHX::Device().fence();
      
      // Check
      Kokkos::deep_copy(a_host,a.get_static_view());
      for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a.extent_int(1); ++pt) {
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).val(),1.0,1e-12);
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).fastAccessDx(0),0.0,1e-12);
	  TEST_FLOATING_EQUALITY(a_host(cell,pt).fastAccessDx(1),0.0,1e-12);
	}
      }

    }    
  }
}
