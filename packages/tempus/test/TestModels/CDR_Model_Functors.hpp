//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_CDR_MODEL_FUNCTORS_HPP
#define TEMPUS_CDR_MODEL_FUNCTORS_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

namespace Tempus_Test {

template <class TpetraVectorType>
struct CoordFiller {
  using SC       = typename TpetraVectorType::impl_scalar_type;
  using LO       = typename TpetraVectorType::local_ordinal_type;
  using GO       = typename TpetraVectorType::global_ordinal_type;
  using Map      = typename TpetraVectorType::map_type;
  using LocalMap = typename Map::local_map_type;
  using DV       = typename TpetraVectorType::dual_view_type;
  using View     = typename DV::t_dev;

  const View coordsView_;
  const SC zMin_;
  const SC dx_;
  const SC minGI_;

  CoordFiller(TpetraVectorType &coords, const SC zMin, const SC dx,
              const GO minGI)
    : coordsView_(coords.getLocalViewDevice(Tpetra::Access::ReadWrite)),
      zMin_(zMin),
      dx_(dx),
      minGI_(minGI)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i) const
  {
    coordsView_(i, 0) = zMin_ + dx_ * static_cast<SC>(minGI_ + i);
  }
};
// Finite Element Basis Object
template <class Scalar, class LO>
class TBasis {
 public:
  // Calculates the values of z and u at the specified Gauss point
  KOKKOS_INLINE_FUNCTION
  void computeBasis(LO gp, Scalar *z, Scalar *u, Scalar *u_dot)
  {
    if (gp == 0) {
      eta = -1.0 / sqrt(3.0);
      wt  = 1.0;
    }
    if (gp == 1) {
      eta = 1.0 / sqrt(3.0);
      wt  = 1.0;
    }

    // Calculate basis function and derivatives at Gauss point
    phi[0]    = (1.0 - eta) / 2.0;
    phi[1]    = (1.0 + eta) / 2.0;
    dphide[0] = -0.5;
    dphide[1] = 0.5;

    // Caculate function and derivative approximations at GP.
    dz      = 0.5 * (z[1] - z[0]);
    zz      = 0.0;
    uu      = 0.0;
    duu     = 0.0;
    uu_dot  = 0.0;
    duu_dot = 0.0;

    for (LO i = 0; i < 2; i++) {
      zz += z[i] * phi[i];
      uu += u[i] * phi[i];
      duu += u[i] * dphide[i];
      uu_dot += u_dot[i] * phi[i];
      duu_dot += u_dot[i] * dphide[i];
    }
  }

 public:
  // Variables that are calculated at the Gauss point
  Scalar phi[2];
  Scalar dphide[2];
  Scalar uu  = 0.0;
  Scalar zz  = 0.0;
  Scalar duu = 0.0;
  Scalar eta = 0.0;
  Scalar wt  = 0.0;
  Scalar dz  = 0.0;

  // These are only needed for transient
  Scalar uu_dot  = 0.0;
  Scalar duu_dot = 0.0;
};

//==================================================================

template <class TpetraVectorType, class TpetraMatrixType>
struct JacobianEvaluatorFunctor {
  using SC        = typename TpetraVectorType::impl_scalar_type;
  using LO        = typename TpetraVectorType::local_ordinal_type;
  using Map       = typename TpetraVectorType::map_type;
  using LocalMap  = typename Map::local_map_type;
  using DV        = typename TpetraVectorType::dual_view_type;
  using ConstView = typename DV::t_dev::const_type;
  using LocalMat  = typename TpetraMatrixType::local_matrix_device_type;

  const LocalMat jLocal_;
  const ConstView xView_;
  const ConstView uView_;
  const ConstView uDotView_;
  const LocalMap rowMap_;
  const LocalMap colMap_;
  const int myRank_;
  const SC a_;
  const SC k_;
  const SC alpha_;
  const SC beta_;

  JacobianEvaluatorFunctor(const TpetraMatrixType &J, const TpetraVectorType &x,
                           const TpetraVectorType &u,
                           const TpetraVectorType &uDot, const int &myRank,
                           const SC a, const SC k, const SC alpha,
                           const SC beta)
    : jLocal_(J.getLocalMatrixDevice()),
      xView_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uView_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uDotView_(uDot.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      rowMap_(J.getRowMap()->getLocalMap()),
      colMap_(J.getColMap()->getLocalMap()),
      myRank_(myRank),
      a_(a),
      k_(k),
      alpha_(alpha),
      beta_(beta)
  {
  }

  // Adds the contribution from element ne to the Jacobian matrix
  KOKKOS_INLINE_FUNCTION
  void operator()(const LO ne) const
  {
    const auto invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = xView_(ne, 0);
    xx[1] = xView_(ne + 1, 0);

    SC uu[2];
    uu[0] = uView_(ne, 0);
    uu[1] = uView_(ne + 1, 0);

    SC uuDot[2];
    uuDot[0] = uDotView_(ne, 0);
    uuDot[1] = uDotView_(ne + 1, 0);

    TBasis<SC, LO> basis;
    // Loop Over Gauss Points
    for (LO gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu, uuDot);

      // Loop over nodes in element
      for (LO i = 0; i < 2; ++i) {
        auto localRow =
            rowMap_.getLocalElement(colMap_.getGlobalElement(ne + i));

        if (localRow != invalid) {
          // Loop over trial functions
          for (LO j = 0; j < 2; ++j) {
            const auto localColumn = ne + j;
            double jac =
                basis.wt * basis.dz *
                (alpha_ * basis.phi[i] * basis.phi[j]  // transient
                 + beta_ * (+a_ / basis.dz * basis.dphide[j] *
                                basis.phi[i]  // convection
                            + (1.0 / (basis.dz * basis.dz)) * basis.dphide[j] *
                                  basis.dphide[i]  // diffusion
                            + 2.0 * k_ * basis.uu * basis.phi[j] *
                                  basis.phi[i]  // source
                            ));

            jLocal_.sumIntoValues(localRow, &localColumn, 1, &jac, false, true);
          }
        }
      }

      // Correct for Dirichlet BCs
      if ((myRank_ == 0) && (ne == 0)) {
        LO row = 0;
        LO col = 0;
        SC val = 1.0;
        jLocal_.replaceValues(row, &col, 1, &val);

        col = 1;
        val = 0.0;
        jLocal_.replaceValues(row, &col, 1, &val);
      }
    }
  }
};

//==================================================================

template <class TpetraVectorType, class TpetraMatrixType>
struct PreconditionerEvaluatorFunctor {
  using SC        = typename TpetraVectorType::impl_scalar_type;
  using LO        = typename TpetraVectorType::local_ordinal_type;
  using Map       = typename TpetraVectorType::map_type;
  using LocalMap  = typename Map::local_map_type;
  using DV        = typename TpetraVectorType::dual_view_type;
  using ConstView = typename DV::t_dev::const_type;
  using LocalMat  = typename TpetraMatrixType::local_matrix_device_type;

  const LocalMat mLocal_;
  const ConstView xView_;
  const ConstView uView_;
  const ConstView uDotView_;
  const LocalMap rowMap_;
  const LocalMap colMap_;
  const int myRank_;
  const SC a_;
  const SC k_;
  const SC alpha_;
  const SC beta_;

  PreconditionerEvaluatorFunctor(const TpetraMatrixType &M,
                                 const TpetraVectorType &x,
                                 const TpetraVectorType &u,
                                 const TpetraVectorType &uDot,
                                 const int &myRank, const SC a, const SC k,
                                 const SC alpha, const SC beta)
    : mLocal_(M.getLocalMatrixDevice()),
      xView_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uView_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uDotView_(uDot.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      rowMap_(M.getRowMap()->getLocalMap()),
      colMap_(M.getColMap()->getLocalMap()),
      myRank_(myRank),
      a_(a),
      k_(k),
      alpha_(alpha),
      beta_(beta)
  {
  }

  // Adds the contribution from element ne to the preconditioner matrix
  KOKKOS_INLINE_FUNCTION
  void operator()(const LO ne) const
  {
    const auto invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = xView_(ne, 0);
    xx[1] = xView_(ne + 1, 0);

    SC uu[2];
    uu[0] = uView_(ne, 0);
    uu[1] = uView_(ne + 1, 0);

    SC uuDot[2];
    uuDot[0] = uDotView_(ne, 0);
    uuDot[1] = uDotView_(ne + 1, 0);

    TBasis<SC, LO> basis;

    // Loop Over Gauss Points
    for (LO gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu, uuDot);

      // Loop over nodes in element
      for (LO i = 0; i < 2; ++i) {
        auto localRow =
            rowMap_.getLocalElement(colMap_.getGlobalElement(ne + i));

        // Loop over trial functions
        if (localRow != invalid) {
          for (LO j = 0; j < 2; ++j) {
            const auto localColumn = ne + j;
            if (rowMap_.getGlobalElement(localRow) ==
                colMap_.getGlobalElement(localColumn)) {
              auto value = basis.wt * basis.dz *
                           (alpha_ * basis.phi[i] * basis.phi[j]  // transient
                            + beta_ * (+a_ / basis.dz * basis.dphide[j] *
                                           basis.phi[i]  // convection
                                       + (1.0 / (basis.dz * basis.dz)) *
                                             basis.dphide[j] *
                                             basis.dphide[i]  // diffusion
                                       + 2.0 * k_ * basis.uu * basis.phi[j] *
                                             basis.phi[i]  // source
                                       ));
              mLocal_.sumIntoValues(localRow, &localColumn, 1, &value);
            }
          }
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((myRank_ == 0) && (ne == 0)) {
      LO row    = 0;
      LO column = 0;
      SC value  = 1.0;
      mLocal_.replaceValues(row, &column, 1, &value);
    }
  }
};

//==================================================================

template <class TpetraVectorType>
struct DfDp2EvaluatorFunctor {
  using SC        = typename TpetraVectorType::impl_scalar_type;
  using LO        = typename TpetraVectorType::local_ordinal_type;
  using Map       = typename TpetraVectorType::map_type;
  using LocalMap  = typename Map::local_map_type;
  using DV        = typename TpetraVectorType::dual_view_type;
  using View      = typename DV::t_dev;
  using ConstView = typename DV::t_dev::const_type;

  const View fView_;
  const ConstView xView_;
  const ConstView uView_;
  const ConstView uDotView_;
  const LocalMap rowMap_;
  const LocalMap colMap_;
  const int myRank_;
  const SC a_;
  const SC k_;

  DfDp2EvaluatorFunctor(TpetraVectorType &f, const TpetraVectorType &x,
                        const TpetraVectorType &u, const TpetraVectorType &uDot,
                        const int myRank, SC a, SC k)
    : fView_(f.getLocalViewDevice(Tpetra::Access::ReadWrite)),
      xView_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uView_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      uDotView_(uDot.getLocalViewDevice(Tpetra::Access::ReadOnly)),
      rowMap_(f.getMap()->getLocalMap()),
      colMap_(u.getMap()->getLocalMap()),
      myRank_(myRank),
      a_(a),
      k_(k)
  {
  }

  // Adds the contribution from element ne to the residual vector
  KOKKOS_INLINE_FUNCTION
  void operator()(const LO ne) const
  {
    const auto invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = xView_(ne, 0);
    xx[1] = xView_(ne + 1, 0);

    SC uu[2];
    uu[0] = uView_(ne, 0);
    uu[1] = uView_(ne + 1, 0);

    SC uuDot[2];
    uuDot[0] = uDotView_(ne, 0);
    uuDot[1] = uDotView_(ne + 1, 0);

    TBasis<SC, LO> basis;
    // Loop Over Gauss Points
    for (LO gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu, uuDot);

      // Loop over nodes in element
      for (LO i = 0; i < 2; ++i) {
        auto localRow =
            rowMap_.getLocalElement(colMap_.getGlobalElement(ne + i));
        if (localRow != invalid) {
          auto value =
              basis.wt * basis.dz *
              (basis.uu_dot * basis.phi[i]                  // transient
               + (a_ / basis.dz * basis.duu * basis.phi[i]  // convection
                  + 1.0 / (basis.dz * basis.dz)) *
                     basis.duu * basis.dphide[i]            // diffusion
               + k_ * basis.uu * basis.uu * basis.phi[i]);  // source
          Kokkos::atomic_add(&fView_(localRow, 0), value);
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((myRank_ == 0) && (ne == 0)) {
      SC value = 0.0;
      Kokkos::atomic_exchange(&fView_(0, 0), value);
    }
  }
};

}  // namespace Tempus_Test

#endif  // TEMPUS_CDR_MODEL_FUNCTORS_HPP
