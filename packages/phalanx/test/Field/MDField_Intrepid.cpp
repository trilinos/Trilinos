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
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// ********************************************************
// Dimension tags for this problem
PHX_DIM_TAG_DECLARATION(Dim)
PHX_DIM_TAG_IMPLEMENTATION(Dim)

PHX_DIM_TAG_DECLARATION(Quadrature)
PHX_DIM_TAG_IMPLEMENTATION(Quadrature)

PHX_DIM_TAG_DECLARATION(Node)
PHX_DIM_TAG_IMPLEMENTATION(Node)

PHX_DIM_TAG_DECLARATION(Point)
PHX_DIM_TAG_IMPLEMENTATION(Point)

PHX_DIM_TAG_DECLARATION(Cell)
PHX_DIM_TAG_IMPLEMENTATION(Cell)

namespace phalanx_test {

  // ********************************************************
  TEUCHOS_UNIT_TEST(mdfield, IntrepidIssue)
  {
#ifdef HAVE_INTREPID_KOKKOSCORE
    using namespace std;
    using namespace Teuchos;
    using namespace PHX;
    
    PHX::InitializeKokkosDevice();

    RCP<DataLayout> dl = rcp(new MDALayout<Cell,Point,Dim>(10,4,3));
    typedef PHX::MDField<double,Cell,Node>::size_type size_type;

    // runtime non-const
    {
      PHX::MDField<double> a("density a",dl);
      a.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(a.fieldTag()));
      
      ArrayWrapper<double, PHX::MDField<double>, Rank<PHX::MDField<double> >::value, false> 
	a_wrap(a);
    }

    // runtime const
    {
      PHX::MDField<double> b("density b",dl);
      b.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(b.fieldTag()));
      
      ArrayWrapper<double, const PHX::MDField<double>, Rank<const PHX::MDField<double> >::value, false> 
	b_wrap(b);
    }

    // compiletime non-const
    {
      PHX::MDField<double,Cell,Point,Dim> c("density c",dl);
      c.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(c.fieldTag()));
      
      ArrayWrapper<double, PHX::MDField<double,Cell,Point,Dim>, Rank<PHX::MDField<double,Cell,Point,Dim> >::value, false> 
	c_wrap(c);
    }

    // compiletime const
    {
      PHX::MDField<double,Cell,Point,Dim> d("density d",dl);
      d.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(d.fieldTag()));
      
      ArrayWrapper<double, const PHX::MDField<double,Cell,Point,Dim>, Rank<const PHX::MDField<double,Cell,Point,Dim> >::value, false> 
	d_wrap(d);
    }


    // NOTE: Need to add a call to an intrepid function here.

    PHX::FinalizeKokkosDevice();
#else
   TEUCHOS_TEST_FOR_EXCEPTION( !(0==1),
                                std::logic_error,
                                "MDField_Intrepid has to have HAVE_INTREPID_KOKKOSCORE ON for this test");

#endif
  }
}
