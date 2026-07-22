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

// This hides the EvaluateDOF functors outside of this file
namespace dof_functors {

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

template <typename ScalarT, typename Array>
class EvaluateDOFWithSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip;
  Array basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateDOFWithSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point> in_dof_ip,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
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

  EvaluateDOFFastSens_Vector(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point,Dim> in_dof_ip,
                             PHX::View<const int*> in_offsets,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis),
      numFields(in_basis.extent(1)),
      numPoints(in_basis.extent(2))
  {}

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

  EvaluateDOFFastSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point> in_dof_ip,
                             PHX::View<const int*> in_offsets,
                             Array in_basis)
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
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
