// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    //   jacobian (physDim,refDim) - jacobian matrix evaluated at a single point
    // input:
    //   grads    (numNodes,refDim) - hgrad basis grad values evaluated at a single point (C1/C2 element only)
    //   nodes    (numNodes,physDim) - cell element-to-node connectivity
    template<typename jacobianViewType,
    typename basisGradViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    computeJacobian(const jacobianViewType  &jacobian, // physDim,refDim
        const basisGradViewType &grads,    // numNodes,refDim
        const nodeViewType      &nodes) {  // numNodes,physDim
      const auto numNodes = nodes.extent(0);

      const auto  physDim = jacobian.extent(0);
      const auto refDim = jacobian.extent(1);

      INTREPID2_TEST_FOR_ABORT( numNodes != grads.extent(0), "grad dimension_0 does not match to cardinality.");
      INTREPID2_TEST_FOR_ABORT(refDim != grads.extent(1), "grad dimension_1 does not match to space dim.");
      INTREPID2_TEST_FOR_ABORT( physDim != nodes.extent(1), "node dimension_1 does not match to space dim.");

      Kernels::Serial::gemm_trans_notrans(1.0, nodes, grads, 0.0, jacobian);
    }

    // output:
    //   point (physDim)   - mapped physical point
    // input:
    //   vals  (numNodes)   - hgrad basis values evaluated at a single point (C1/C2 element only)
    //   nodes (numNodes,physDim) - cell element-to-node connectivity
    template<typename PointViewType,
    typename basisValViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    mapToPhysicalFrame(const PointViewType    &point,    // physDim
        const basisValViewType &vals,     // numNodes
        const nodeViewType     &nodes) {  // numNodes,physDim
      const auto numNodes = vals.extent(0);
      const auto physDim = point.extent(0);

      INTREPID2_TEST_FOR_ABORT(numNodes != nodes.extent(0), "nodes dimension_0 does not match to vals dimension_0.");
      INTREPID2_TEST_FOR_ABORT(physDim != nodes.extent(1), "node dimension_1 does not match to space dim.");

      Kernels::Serial::gemv_trans(1.0, nodes, vals, 0.0, point);
    }

    /** \brief  Computation of \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map.

        Applies \f$ F^{-1}_{c} \f$ for \b a cell workset to \b a single point (P,D) view.
        Computes a rank-2 (P,D) view such that
        \f[
        \mbox{refPoints}(p,d) = \Big(F^{-1}_c(physPoint(p,*)) \Big)_d
        \f]

        Requires pointer to HGrad basis that defines reference to physical cell mapping.
        See Section \ref sec_cell_topology_ref_map for definition of the mapping function.
        Note that the physical point can be in a space with larger dimension, e.g. if mapping a
        reference quadrilateral into a side in 3D space.

        \warning
        The array \c physPoints represents an arbitrary set (or sets) of points in the physical
        frame that are not required to belong in the physical cell (cells) that define(s) the reference
        to physical mapping. As a result, the images of these points in the reference frame
        are not necessarily contained in the reference cell corresponding to the specified
        cell topology.

        \param  refPoint        [out] - rank-2 view with dimension (refD) with the reference point
        \param  physPoint       [in]  - rank-2 view with dimension (physD) with the point in physical frame
        \param  cellNodes       [in]  - rank-2 array with dimensions (N,physD) with the nodes of the cell containing the physical point
        \return a boolean set to true if the algorithm converged, to false otherise.
    */
    template<typename implBasisType,
    typename refPointViewType,
    typename phyPointViewType,
    typename nodeViewType>
    KOKKOS_INLINE_FUNCTION
    static bool
    mapToReferenceFrame(const refPointViewType &refPoint, // refDim
        const phyPointViewType &physPoint, // physDim
        const nodeViewType &nodes,
        const double tol = tolerence()) {  // numNodes,physDim
      const ordinal_type refDim = refPoint.extent(0);
      const ordinal_type physDim = physPoint.extent(0);
      const ordinal_type numNodes = nodes.extent(0);

      INTREPID2_TEST_FOR_ABORT(refDim > physDim, "the dimension of the reference cell is greater than physical cell dimension.");
      INTREPID2_TEST_FOR_ABORT(physDim != static_cast<ordinal_type>(nodes.extent(1)), "physPoint dimension_0 does not match to space dim.");
      INTREPID2_TEST_FOR_ABORT(numNodes > 27, "function hard-coded to support at most mappings with 27 Dofs");

      typedef typename refPointViewType::non_const_value_type value_type;

      // I want to use view instead of dynrankview
      // NMAX = 27, MAXDIM = 3
      value_type buf[27*3 + 27 + 9 + 9 + 9 + 9 + 3 + 3] = {}, *ptr = &buf[0];
      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> grads(ptr, numNodes, refDim); ptr += numNodes*refDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> vals(ptr, numNodes); ptr += numNodes;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> jac(ptr, physDim, refDim); ptr += physDim*refDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> metric(ptr, refDim, refDim); ptr += refDim*refDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> invMetric(ptr, refDim, refDim); ptr += refDim*refDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> invDf(ptr, refDim, physDim); ptr += refDim*physDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> tmpPhysPoint(ptr, physDim); ptr += physDim;

      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> oldRefPoint(ptr, refDim); ptr += refDim;

      // set initial guess
      for (ordinal_type j=0;j<refDim;++j) oldRefPoint(j) = 0;

      double xphyNorm = Kernels::Serial::norm(physPoint, NORM_ONE);

      bool converged = false;
      for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
        // xtmp := F(oldRefPoint);
        implBasisType::template Serial<OPERATOR_VALUE>::getValues(vals, oldRefPoint);
        mapToPhysicalFrame(tmpPhysPoint, vals, nodes);
        Kernels::Serial::z_is_axby(tmpPhysPoint, 1.0, physPoint, -1.0, tmpPhysPoint);  //residual  xtmp := physPoint - F(oldRefPoint);
        double residualNorm = Kernels::Serial::norm(tmpPhysPoint, NORM_ONE);
        if (residualNorm < tol*xphyNorm) {
          converged = true;
          break;
        }

        // DF^{-1}
        implBasisType::template Serial<OPERATOR_GRAD>::getValues(grads, oldRefPoint);
        CellTools::Serial::computeJacobian(jac, grads, nodes);

        if(physDim == refDim) { //jacobian inverse
          Kernels::Serial::inverse(invDf, jac);
        } else { //jacobian is not square, we compute the pseudo-inverse (jac^T jac)^{-1} jac^T
          Kernels::Serial::gemm_trans_notrans(1.0, jac, jac, 0.0, metric);
          Kernels::Serial::inverse(invMetric, metric);
          Kernels::Serial::gemm_notrans_trans(1.0, invMetric, jac, 0.0, invDf);
        }

        // Newton
        Kernels::Serial::gemv_notrans(1.0, invDf, tmpPhysPoint, 0.0, refPoint); // refPoint := DF^{-1}( physPoint - F(oldRefPoint))
        Kernels::Serial::z_is_axby(refPoint, 1.0, oldRefPoint,  1.0, refPoint); // refPoint += oldRefPoint

        // temporarily overwriting oldRefPoint with oldRefPoint-refPoint
        Kernels::Serial::z_is_axby(oldRefPoint, 1.0, oldRefPoint, -1.0, refPoint);

        double err = Kernels::Serial::norm(oldRefPoint, NORM_ONE);
        if (err < tol) {
          converged = true;
          break;
        }

        Kernels::Serial::copy(oldRefPoint, refPoint);
      }
      return converged;
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
    getPhysicalEdgeTangent(const edgeTangentViewType edgeTangent, // physDim
        const ParamViewType edgeParametrization,
        const jacobianViewType jacobian,         // physDim, physDim
        const ordinal_type edgeOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = edgeParametrization.extent(1);
      value_type buf[3];
      Kokkos::DynRankView<value_type, Kokkos::AnonymousSpace,
      Kokkos::MemoryUnmanaged> refEdgeTangent(&buf[0], dim);

      getReferenceEdgeTangent(refEdgeTangent, edgeParametrization, edgeOrdinal);
      Kernels::Serial::matvec_product(edgeTangent, jacobian, refEdgeTangent);
    }

    template<typename faceTanViewType, typename ParamViewType, typename jacobianViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getPhysicalFaceTangents( faceTanViewType faceTanU, // physDim
        faceTanViewType faceTanV, // physDim
        const ParamViewType faceParametrization,
        const jacobianViewType jacobian,       // physDim, physDim
        const ordinal_type faceOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = faceParametrization.extent(1);
      INTREPID2_TEST_FOR_ABORT(dim != 3,
          "computing face tangents requires dimension 3.");
      value_type buf[6];
      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
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
    getPhysicalFaceNormal(const faceNormalViewType faceNormal, // physDim
        const ParamViewType faceParametrization,
        const jacobianViewType jacobian,       // physDim, physDim
        const ordinal_type faceOrdinal) {
      typedef typename ParamViewType::non_const_value_type value_type;
      const ordinal_type dim = faceParametrization.extent(1);
      INTREPID2_TEST_FOR_ABORT(dim != 3,
          "computing face normal requires dimension 3.");
      value_type buf[6];
      Kokkos::DynRankView<value_type,
      Kokkos::AnonymousSpace,
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
    getPhysicalSideNormal(const sideNormalViewType sideNormal, // physDim
        const ParamViewType sideParametrization,
        const jacobianViewType jacobian,       // physDim, physDim
        const ordinal_type sideOrdinal) {
      const ordinal_type dim = sideParametrization.extent(1);
      typedef typename ParamViewType::non_const_value_type value_type;
      switch (dim) {
      case 2: {
        value_type buf[3];
        Kokkos::DynRankView<value_type,
        Kokkos::AnonymousSpace,
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
        INTREPID2_TEST_FOR_ABORT(true, "cell dimension is out of range.");
        break;
      }
      }
    }
  };
};
}
}

#endif

