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
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// ********************************************************
// Dimension tags for this problem
PHX_EXTENT(Dim)
PHX_EXTENT(Quadrature)
PHX_EXTENT(Node)
PHX_EXTENT(Point)
PHX_EXTENT(Cell)

/*
 This test simulates an issue with the implementation of Intrepid. To
 prevent a combinatorial explosion, Intrepid functions switch on rank
 in a single function and require dynamic rank views.
 */

using namespace std;
using namespace Teuchos;
using namespace PHX;
using namespace Kokkos;
using size_type = PHX::MDField<double>::size_type;

namespace phalanx_test {

  // ********************************************************
  template<typename VectorType>
  void simulated_intrepid_integrate(VectorType v)
  {
    parallel_for("simulated intrepid rank-3",
		 RangePolicy<PHX::Device>(0,10),
		 KOKKOS_LAMBDA (const int i)
    {
      if (v.rank() == 3) {
	for (size_type j=0; j < v.extent(1); ++j)
	  for (size_type k=0; k < v.extent(2); ++k)
	    v(i,j,k) = double(i+j+k);
      }
      else {
	for (size_type j=0; j < v.extent(1); ++j)
	  v(i,j) = double(i+j);
      }
    });
  }

  // ********************************************************
  TEUCHOS_UNIT_TEST(mdfield, IntrepidIssue)
  {
    PHX::MDField<double> b2("density 2","QP Layout",10,4);
    b2.deep_copy(2.0);

    PHX::MDField<double> b3("density 3","QP Layout",10,4,3);
    b3.deep_copy(3.0);

    simulated_intrepid_integrate(b2.get_view());
    simulated_intrepid_integrate(b3.get_view());
  }
}
