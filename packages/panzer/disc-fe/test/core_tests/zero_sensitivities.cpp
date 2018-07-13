// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

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
      a.setFieldData(PHX::KokkosViewFactory<panzer::Traits::RealType,PHX::Device>::buildView(a.fieldTag()));
      
      // initialize
      for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a.extent_int(1); ++pt) {
	  a(cell,pt) = 2.0;
	}
      }
      
      // Compute
      Kokkos::parallel_for(a.extent_int(0),ComputeA<panzer::Traits::RealType,PHX::Device,PHX::MDField<panzer::Traits::RealType,Cell,BASIS> > (a));
      PHX::Device::fence();
      
      // Check
      for (int cell = 0; cell < a.extent_int(0); ++cell)
	for (int pt=0; pt < a.extent_int(1); ++pt)
	  TEST_FLOATING_EQUALITY(a(cell,pt),2.0,tolerance);
    }

    // ***********************
    // FadType
    // ***********************
    {
      PHX::MDField<panzer::Traits::FadType,Cell,BASIS> a("rho",dl);
      std::vector<size_type> derivative_dimension;
      derivative_dimension.push_back(2);
      a.setFieldData(PHX::KokkosViewFactory<panzer::Traits::FadType,PHX::Device>::buildView(a.fieldTag(),derivative_dimension));
      
      // initialize
      for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a.extent_int(1); ++pt) {
	  a(cell,pt) = 1.0;
	  a(cell,pt).fastAccessDx(0) = 2.0;
	  a(cell,pt).fastAccessDx(1) = 3.0;
	  TEST_FLOATING_EQUALITY(a(cell,pt).val(),1.0,1e-12);
	  TEST_FLOATING_EQUALITY(a(cell,pt).fastAccessDx(0),2.0,1e-12);
	  TEST_FLOATING_EQUALITY(a(cell,pt).fastAccessDx(1),3.0,1e-12);
	}
      }
      
      // Compute
      Kokkos::parallel_for(a.extent(0),ComputeB<PHX::Device,PHX::MDField<panzer::Traits::FadType,Cell,BASIS> > (a));
      PHX::Device::fence();
      
      // Check
      for (int cell = 0; cell < a.extent_int(0); ++cell) {
	for (int pt=0; pt < a.extent_int(1); ++pt) {
	  TEST_FLOATING_EQUALITY(a(cell,pt).val(),1.0,1e-12);
	  TEST_FLOATING_EQUALITY(a(cell,pt).fastAccessDx(0),0.0,1e-12);
	  TEST_FLOATING_EQUALITY(a(cell,pt).fastAccessDx(1),0.0,1e-12);
	}
      }

    }    
  }
}
