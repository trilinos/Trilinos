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

using namespace Teuchos;
using namespace PHX;
using size_type = PHX::MDField<double>::size_type;

namespace phalanx_test {

  // ********************************************************
  template<typename VectorType>
  void simulated_intrepid_integrate(VectorType v)
  {
    Kokkos::parallel_for("simulated intrepid rank-3",
		 Kokkos::RangePolicy<PHX::Device>(0,10),
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
