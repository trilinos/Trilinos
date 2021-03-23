// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_CellTools_Serial.hpp
    \brief  Definition file for the Intrepid2::Impl::CellTools class.
    \author Kyungjoo Kim
 */

#ifndef __INTREPID2_CELLTOOLS_SERIAL_HPP__
#define __INTREPID2_CELLTOOLS_SERIAL_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Kernels.hpp"

#include "Intrepid2_CellData.hpp"

namespace Intrepid2 {

namespace Impl {

/**
    \brief See Intrepid2::CellTools
 */

class CellTools {
public:

  struct Serial {

    // output:
    //   jacobian (D,sD) - jacobian matrix evaluated at a single point
    // input:
    //   grads    (N,sD) - hgrad basis grad values evaluated at a single point (C1/C2 element only)
    //   nodes    (N,D) - cell element-to-node connectivity
    template<typename jacobianViewType,
    typename basisGradViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    computeJacobian(const jacobianViewType  &jacobian, // D,sD
        const basisGradViewType &grads,    // N,sD
        const nodeViewType      &nodes) {  // N,D
      const auto N = nodes.extent(0);

      const auto  D = jacobian.extent(0);
      const auto sD = jacobian.extent(1);

      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE( N != grads.extent(0), "grad dimension_0 does not match to cardinality.");
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(sD != grads.extent(1), "grad dimension_1 does not match to space dim.");
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE( D != nodes.extent(1), "node dimension_1 does not match to space dim.");

      Kernels::Serial::gemm_trans_notrans(1.0, nodes, grads, 0.0, jacobian);
    }

    // output:
    //   point (D)   - mapped physical point
    // input:
    //   vals  (N)   - hgrad basis values evaluated at a single point (C1/C2 element only)
    //   nodes (N,D) - cell element-to-node connectivity
    template<typename PointViewType,
    typename basisValViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    mapToPhysicalFrame(const PointViewType    &point,    // D
        const basisValViewType &vals,     // N
        const nodeViewType     &nodes) {  // N,D
      const auto N = vals.extent(0);
      const auto D = point.extent(0);

      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(N != nodes.extent(0), "nodes dimension_0 does not match to vals dimension_0.");
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(D != nodes.extent(1), "node dimension_1 does not match to space dim.");

      Kernels::Serial::gemv_trans(1.0, nodes, vals, 0.0, point);
    }

    // template:
    //   implBasisType - impl basis function type e.g., Impl::Basis_HGRAD_QUAD_C1_FEM
    // output:
    //   xref (sD)   - point mapped to reference frame (subcell Dim)
    // input:
    //   xphy  (D)   - point in physical frame
    //   nodes (N,D) - cell element-to-node connectivity
    template<typename implBasisType,
    typename refPointViewType,
    typename phyPointViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    mapToReferenceFrame(const refPointViewType &xref, // sD
        const phyPointViewType &xphy, // D
        const nodeViewType &nodes) {  // N,D
      const ordinal_type sD = xref.extent(0);
      const ordinal_type D = xphy.extent(0);
      const ordinal_type N = nodes.extent(0);

      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(sD > D, "subcell dimension is greater than physical cell dimension.");
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(D != static_cast<ordinal_type>(nodes.extent(1)), "xphy dimension_0 does not match to space dim.");

      typedef typename refPointViewType::non_const_value_type value_type;

      // I want to use view instead of dynrankview
      // NMAX = 28, MAXDIM = 3
      value_type buf[27*3 + 27 + 9 + 9 + 9 + 9 + 3 + 3] = {}, *ptr = &buf[0];
      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> grads(ptr, N, sD); ptr += N*sD;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> vals(ptr, N); ptr += N;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> jac(ptr, D, sD); ptr += D*sD;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> metric(ptr, sD, sD); ptr += sD*sD;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> invmetric(ptr, sD, sD); ptr += sD*sD;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> invdf(ptr, sD, D); ptr += sD*D;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> xtmp(ptr, sD); ptr += sD;

      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> xold(ptr, sD); ptr += sD;

      // set initial guess
      for (ordinal_type j=0;j<D;++j) xold(j) = 0;

      const double tol = tolerence();
      for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
        // xtmp := F(xold);
        implBasisType::template Serial<OPERATOR_VALUE>::getValues(vals, xold);
        mapToPhysicalFrame(xtmp, vals, nodes);

        // DF^{-1}
        implBasisType::template Serial<OPERATOR_GRAD>::getValues(grads, xold);
        CellTools::Serial::computeJacobian(jac, grads, nodes);

        Kernels::Serial::gemm_trans_notrans(1.0, jac, jac, 0.0, metric);
        Kernels::Serial::inverse(invmetric, metric);
        Kernels::Serial::gemm_notrans_trans(1.0, invmetric, jac, 0.0, invdf);

        // Newton
        Kernels::Serial::z_is_axby(xtmp, 1.0, xphy, -1.0, xtmp);  // xtmp := xphy - F(xold);
        Kernels::Serial::gemv_notrans(1.0, invdf, xtmp, 0.0, xref); // xref := DF^{-1}( xphy - F(xold))
        Kernels::Serial::z_is_axby(xref, 1.0, xold,  1.0, xref); // xref += xold

        // l2 error
        Kernels::Serial::z_is_axby(xtmp, 1.0, xold, -1.0, xref);

        double err = Kernels::Serial::norm(xtmp, NORM_ONE);

        if (err < tol)
          break;

        Kernels::Serial::copy(xold, xref);
      }
    }

    template<typename refEdgeTanViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getReferenceEdgeTangent(const refEdgeTanViewType refEdgeTangent,
        const ParamViewType edgeParametrization,
        const ordinal_type edgeOrdinal ) {

      ordinal_type dim = edgeParametrization.extent(1);
      for(ordinal_type i = 0; i < dim; ++i) {
        refEdgeTangent(i) = edgeParametrization(edgeOrdinal, i, 1);
      }
    }

    template<typename refFaceTanViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getReferenceFaceTangent(const refFaceTanViewType refFaceTanU,
        const refFaceTanViewType refFaceTanV,
        const ParamViewType faceParametrization,
        const ordinal_type faceOrdinal) {

      // set refFaceTanU -> C_1(*)
      // set refFaceTanV -> C_2(*)
      ordinal_type dim = faceParametrization.extent(1);
      for(ordinal_type i = 0; i < dim; ++i) {
        refFaceTanU(i) = faceParametrization(faceOrdinal, i, 1);
        refFaceTanV(i) = faceParametrization(faceOrdinal, i, 2);
      }
    }

    template<typename edgeTangentViewType, typename ParamViewType, typename jacobianViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getPhysicalEdgeTangent(const edgeTangentViewType edgeTangent, // D
        const ParamViewType edgeParametrization,
        const jacobianViewType jacobian,         // D, D
        const ordinal_type edgeOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = edgeParametrization.extent(1);
      value_type buf[3];
      Kokkos::DynRankView<value_type, Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> refEdgeTangent(&buf[0], dim);

      getReferenceEdgeTangent(refEdgeTangent, edgeParametrization, edgeOrdinal);
      Kernels::Serial::matvec_product(edgeTangent, jacobian, refEdgeTangent);
    }

    template<typename faceTanViewType, typename ParamViewType, typename jacobianViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getPhysicalFaceTangents( faceTanViewType faceTanU, // D
        faceTanViewType faceTanV, // D
        const ParamViewType faceParametrization,
        const jacobianViewType jacobian,       // D, D
        const ordinal_type faceOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = faceParametrization.extent(1);
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(dim != 3,
          "computing face tangents requires dimension 3.");
      value_type buf[6];
      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> refFaceTanU(&buf[0], dim), refFaceTanV(&buf[3], dim);

      getReferenceFaceTangent(refFaceTanU,
          refFaceTanV,
          faceParametrization,
          faceOrdinal);

      Kernels::Serial::matvec_product_d3(faceTanU, jacobian, refFaceTanU);
      Kernels::Serial::matvec_product_d3(faceTanV, jacobian, refFaceTanV);
    }


    template<typename faceNormalViewType, typename ParamViewType, typename jacobianViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getPhysicalFaceNormal(const faceNormalViewType faceNormal, // D
        const ParamViewType faceParametrization,
        const jacobianViewType jacobian,       // D, D
        const ordinal_type faceOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = faceParametrization.extent(1);
      INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(dim != 3,
          "computing face normal requires dimension 3.");
      value_type buf[6];
      Kokkos::DynRankView<value_type,
      Kokkos::Impl::ActiveExecutionMemorySpace,
      Kokkos::MemoryUnmanaged> faceTanU(&buf[0], dim), faceTanV(&buf[3], dim);

      getPhysicalFaceTangents(faceTanU, faceTanV,
          faceParametrization,
          jacobian,
          faceOrdinal);
      Kernels::Serial::vector_product_d3(faceNormal, faceTanU, faceTanV);
    }

    template<typename sideNormalViewType, typename ParamViewType, typename jacobianViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getPhysicalSideNormal(const sideNormalViewType sideNormal, // D
        const ParamViewType sideParametrization,
        const jacobianViewType jacobian,       // D, D
        const ordinal_type sideOrdinal) {
      const ordinal_type dim = sideParametrization.extent(1);
      typedef typename ParamViewType::non_const_value_type value_type;
      switch (dim) {
      case 2: {
        value_type buf[3];
        Kokkos::DynRankView<value_type,
        Kokkos::Impl::ActiveExecutionMemorySpace,
        Kokkos::MemoryUnmanaged> edgeTangent(&buf[0], dim);
        getPhysicalEdgeTangent(edgeTangent, sideParametrization, jacobian, sideOrdinal);
        sideNormal(0) =  edgeTangent(1);
        sideNormal(1) = -edgeTangent(0);
        break;
      }
      case 3: {
        getPhysicalFaceNormal(sideNormal, sideParametrization, jacobian, sideOrdinal);
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT_DEVICE_SAFE(true, "cell dimension is out of range.");
        break;
      }
      }
    }
  };
};
}
}

#endif

