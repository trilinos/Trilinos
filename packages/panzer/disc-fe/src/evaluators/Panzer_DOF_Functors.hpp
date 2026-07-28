// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_DOF_Functors_hpp__
#define __Panzer_DOF_Functors_hpp__

#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

//**********************************************************************

/**
 * \brief Device functors that interpolate DOF coefficients to values at integration points.
 *
 * Each functor here performs one Kokkos-parallel evaluation of
 * `dof_ip(cell,pt) = sum_bf dof_basis(cell,bf) * basis(cell,bf,pt)`,
 * i.e. it computes the value (or vector/gradient components) of a
 * field at every integration point of every cell in a workset, given
 * per-basis-function DOF coefficients. They are used by evaluators
 * such as panzer::DOF and panzer::DOFGradient.
 *
 * Two evaluation strategies are provided:
 *  - `*WithSens`: propagates derivatives with respect to all local DOF
 *    coefficients using whatever AD type `ScalarT` already carries.
 *    Correct in general, but loops over every derivative component.
 *  - `*FastSens`: exploits the fact that `dof_ip`'s derivative with
 *    respect to a given cell-local DOF depends only on that DOF's
 *    entry in `offsets`, so the functor sets exactly that one nonzero
 *    derivative component via `fastAccessDx()` instead of iterating
 *    over the whole FAD array. This is a significant performance
 *    optimization for the Jacobian evaluation type, but assumes the
 *    `basis` values themselves carry no derivative information (see
 *    the per-functor note about coordinate sensitivities).
 *
 * `*_Vector` variants compute all `spaceDim` components of a
 * vector-valued (or gradient) field; `*_Scalar` variants compute a
 * single scalar field.
 */
// This hides the EvaluateDOF functors outside of this file
namespace dof_functors {

/**
 * \brief Vector-valued DOF interpolation functor with full derivative propagation.
 * \tparam ScalarT the field scalar type (may carry AD derivatives).
 * \tparam Array the basis-value array type.
 * \tparam spaceDim number of spatial dimensions (vector components) to compute.
 * \sa dof_functors
 */
template <typename ScalarT,typename Array,int spaceDim>
class EvaluateDOFWithSens_Vector {
  PHX::View<const ScalarT**> dof_basis; // <C,P>
  PHX::View<ScalarT***> dof_ip; // <C,P,D>
  Array basis;
  const int numFields;
  const int numPoints;
  const int fadSize;
  const bool use_shared_memory;

public:
  using scratch_view = Kokkos::View<ScalarT* ,typename PHX::DevLayout<ScalarT>::type,typename PHX::exec_space::scratch_memory_space,Kokkos::MemoryUnmanaged>;

  /**
   * \brief Constructor.
   * \param in_dof_basis input field of per-cell, per-basis-function DOF coefficients.
   * \param in_dof_ip output field values at integration points.
   * \param in_basis basis values at integration points.
   * \param in_use_shared_memory if true, stage per-team data in Kokkos team scratch memory (see team_shmem_size()) to reduce redundant global-memory traffic across the contraction.
   */
  EvaluateDOFWithSens_Vector(PHX::View<const ScalarT**> in_dof_basis,
                             PHX::View<ScalarT***> in_dof_ip,
                             Array in_basis,
			     bool in_use_shared_memory = false)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), basis(in_basis),
      numFields(static_cast<int>(basis.extent(1))),
      numPoints(static_cast<int>(basis.extent(2))),
      fadSize(static_cast<int>(Kokkos::dimension_scalar(dof_basis))),
      use_shared_memory(in_use_shared_memory)
  {}

  /// \brief Computes dof_ip for one cell (team.league_rank()), optionally staging through team scratch memory.
  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    const int cell = team.league_rank();

    if (not use_shared_memory) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numPoints), [&] (const int& pt) {
	for (int d=0; d<spaceDim; ++d) {
	  // first initialize to the right thing (prevents over writing with 0)
	  // then loop over one less basis function
	  dof_ip(cell,pt,d) = dof_basis(cell, 0) * basis(cell, 0, pt, d);
	  // The start index is one, not zero since we used the zero index for initialization above.
	  for (int bf=1; bf<numFields; ++bf) {
	    dof_ip(cell,pt,d) += dof_basis(cell, bf) * basis(cell, bf, pt, d);
	  }
	}
      });
    }
    else {

      // Copy reused data into fast scratch space
      scratch_view dof_values;
      scratch_view point_values;
      if (Sacado::IsADType<ScalarT>::value) {
        dof_values = scratch_view(team.team_shmem(),numFields,fadSize);
        point_values = scratch_view(team.team_shmem(),numPoints,fadSize);
      }
      else {
        dof_values = scratch_view(team.team_shmem(),numFields);
        point_values = scratch_view(team.team_shmem(),numPoints);
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numFields), [&] (const int& dof) {
	dof_values(dof) = dof_basis(cell,dof);
      });

      team.team_barrier();

      for (int dim=0; dim < spaceDim; ++dim) {

	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numPoints), [&] (const int& pt) {
	  point_values(pt) = 0.0;
	});

	// Perform contraction
	for (int dof=0; dof<numFields; ++dof) {
	  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numPoints), [&] (const int& pt) {
	    point_values(pt) += dof_values(dof) * basis(cell,dof,pt,dim);
          });
	}

	// Copy to main memory
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numPoints), [&] (const int& pt) {
	  dof_ip(cell,pt,dim) = point_values(pt);
        });

      } // loop over dim
    } // if (use_shared_memory) {
  }

  /// \brief Returns the team scratch memory size (in bytes) required by operator() when use_shared_memory is enabled; 0 otherwise.
  size_t team_shmem_size(int /* team_size */ ) const
  {
    if (not use_shared_memory)
      return 0;

    size_t bytes;
    if (Sacado::IsADType<ScalarT>::value)
      bytes = scratch_view::shmem_size(numFields,fadSize) + scratch_view::shmem_size(numPoints,fadSize);
    else
      bytes = scratch_view::shmem_size(numFields) + scratch_view::shmem_size(numPoints);
    return bytes;
  }

};

/**
 * \brief Scalar-valued DOF interpolation functor with full derivative propagation.
 * \tparam ScalarT the field scalar type (may carry AD derivatives).
 * \tparam Array the basis-value array type.
 * \sa dof_functors
 */
template <typename ScalarT, typename Array>
class EvaluateDOFWithSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip;
  Array basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  /**
   * \brief Constructor.
   * \param in_dof_basis input field of per-cell, per-basis-function DOF coefficients.
   * \param in_dof_ip output field values at integration points.
   * \param in_basis basis values at integration points.
   */
  EvaluateDOFWithSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point> in_dof_ip,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
  /// \brief Computes dof_ip for all points of one cell.
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      // first initialize to the right thing (prevents over writing with 0)
      // then loop over one less basis function
      dof_ip(cell,pt) = dof_basis(cell, 0) * basis(cell, 0, pt);
      for (int bf=1; bf<numFields; bf++) {
        dof_ip(cell,pt) += dof_basis(cell, bf) * basis(cell, bf, pt);
      }
    }
  }
};

/**
 * \brief Vector-valued DOF interpolation functor using the fast-sensitivity optimization.
 *
 * Assumes each cell-local DOF's derivative contribution is confined to
 * a single entry of \c offsets, so only that one derivative component
 * is computed via \c fastAccessDx() rather than propagating a full AD
 * derivative array. Not valid if \c basis itself carries derivatives
 * (e.g. sensitivity to coordinates).
 * \tparam ScalarT the field scalar type (AD type carrying derivatives).
 * \tparam Array the basis-value array type.
 * \tparam spaceDim number of spatial dimensions (vector components) to compute.
 * \sa dof_functors
 */
template <typename ScalarT,typename Array,int spaceDim>
class EvaluateDOFFastSens_Vector {
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_ip;
  PHX::View<const int*> offsets;
  Array basis;

  const int numFields;
  const int numPoints;

public:
  typedef typename PHX::Device execution_space;

  /**
   * \brief Constructor.
   * \param in_dof_basis input field of per-cell, per-basis-function DOF coefficients.
   * \param in_dof_ip output field values at integration points.
   * \param in_offsets maps each basis function to the derivative array index holding its (only nonzero) sensitivity.
   * \param in_basis basis values at integration points.
   */
  EvaluateDOFFastSens_Vector(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point,Dim> in_dof_ip,
                             PHX::View<const int*> in_offsets,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis),
      numFields(in_basis.extent(1)),
      numPoints(in_basis.extent(2))
  {}

  /// \brief Computes dof_ip (value and its one nonzero sensitivity) for all points/dimensions of one cell.
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      for (int d=0; d<spaceDim; d++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function

        // This is a possible issue if you need sensitivity to coordinates (you will need to
        // change basis and then use the product rule!)
        dof_ip(cell,pt,d) = dof_basis(cell, 0).val() * basis(cell, 0, pt, d);
        dof_ip(cell,pt,d).fastAccessDx(offsets(0)) = dof_basis(cell, 0).fastAccessDx(offsets(0)) * Sacado::scalarValue(basis(cell, 0, pt, d));

        for (int bf=1; bf<numFields; bf++) {
          dof_ip(cell,pt,d).val() += dof_basis(cell, bf).val() * Sacado::scalarValue(basis(cell, bf, pt, d));
          dof_ip(cell,pt,d).fastAccessDx(offsets(bf)) += dof_basis(cell, bf).fastAccessDx(offsets(bf)) * Sacado::scalarValue(basis(cell, bf, pt, d));
        }
      }
    }
  }
};

/**
 * \brief Scalar-valued DOF interpolation functor using the fast-sensitivity optimization.
 *
 * Assumes each cell-local DOF's derivative contribution is confined to
 * a single entry of \c offsets, so only that one derivative component
 * is computed via \c fastAccessDx() rather than propagating a full AD
 * derivative array. Not valid if \c basis itself carries derivatives
 * (e.g. sensitivity to coordinates).
 * \tparam ScalarT the field scalar type (AD type carrying derivatives).
 * \tparam Array the basis-value array type.
 * \sa dof_functors
 */
template <typename ScalarT, typename Array>
class EvaluateDOFFastSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip;
  PHX::View<const int*> offsets;
  Array basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  /**
   * \brief Constructor.
   * \param in_dof_basis input field of per-cell, per-basis-function DOF coefficients.
   * \param in_dof_ip output field values at integration points.
   * \param in_offsets maps each basis function to the derivative array index holding its (only nonzero) sensitivity.
   * \param in_basis basis values at integration points.
   */
  EvaluateDOFFastSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point> in_dof_ip,
                             PHX::View<const int*> in_offsets,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
  /// \brief Computes dof_ip (value and its one nonzero sensitivity) for all points of one cell.
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      // first initialize to the right thing (prevents over writing with 0)
      // then loop over one less basis function

      // This is a possible issue if you need sensitivity to coordinates (you will need to
      // change basis and then use the product rule!)
      dof_ip(cell,pt) = dof_basis(cell, 0).val() * Sacado::scalarValue(basis(cell, 0, pt));
      dof_ip(cell,pt).fastAccessDx(offsets(0)) = dof_basis(cell, 0).fastAccessDx(offsets(0)) * Sacado::scalarValue(basis(cell, 0, pt));

      for (int bf=1; bf<numFields; bf++) {
        dof_ip(cell,pt).val() += dof_basis(cell, bf).val() * Sacado::scalarValue(basis(cell, bf, pt));
        dof_ip(cell,pt).fastAccessDx(offsets(bf)) += dof_basis(cell, bf).fastAccessDx(offsets(bf)) * Sacado::scalarValue(basis(cell, bf, pt));
      }
    }
  }
};

}

}

#endif
