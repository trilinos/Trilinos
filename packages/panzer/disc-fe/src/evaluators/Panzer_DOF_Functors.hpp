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

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateDOFWithSens_Vector(PHX::View<const ScalarT**> in_dof_basis,
                             PHX::View<ScalarT***> in_dof_ip,
                             Array in_basis) 
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), basis(in_basis)
  {
    numFields = static_cast<int>(basis.extent(1));
    numPoints = static_cast<int>(basis.extent(2));
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    const int cell = team.league_rank();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numPoints),KOKKOS_LAMBDA (const int& pt) {
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
  Kokkos::View<const int*,PHX::Device> offsets;
  Array basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateDOFFastSens_Vector(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point,Dim> in_dof_ip,
                             Kokkos::View<const int*,PHX::Device> in_offsets,
                             Array in_basis) 
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    typedef Sacado::ScalarValue<ScalarT> Value;

    for (int pt=0; pt<numPoints; pt++) {
      for (int d=0; d<spaceDim; d++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function

        // This is a possible issue if you need sensitivity to coordinates (you will need to
        // change basis and then use the product rule!)
        dof_ip(cell,pt,d) = dof_basis(cell, 0).val() * basis(cell, 0, pt, d);
        dof_ip(cell,pt,d).fastAccessDx(offsets(0)) = dof_basis(cell, 0).fastAccessDx(offsets(0)) * Value::eval(basis(cell, 0, pt, d));

        for (int bf=1; bf<numFields; bf++) {
          dof_ip(cell,pt,d).val() += dof_basis(cell, bf).val() * Value::eval(basis(cell, bf, pt, d));
          dof_ip(cell,pt,d).fastAccessDx(offsets(bf)) += dof_basis(cell, bf).fastAccessDx(offsets(bf)) * Value::eval(basis(cell, bf, pt, d));
        }
      }
    }
  }
};

template <typename ScalarT, typename Array>
class EvaluateDOFFastSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip;
  Kokkos::View<const int*,PHX::Device> offsets;
  Array basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateDOFFastSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_basis,
                             PHX::MDField<ScalarT,Cell,Point> in_dof_ip,
                             Kokkos::View<const int*,PHX::Device> in_offsets,
                             Array in_basis) 
    : dof_basis(in_dof_basis), dof_ip(in_dof_ip), offsets(in_offsets), basis(in_basis)
  {
    numFields = basis.extent(1);
    numPoints = basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    typedef Sacado::ScalarValue<ScalarT> Value;

    for (int pt=0; pt<numPoints; pt++) {
      // first initialize to the right thing (prevents over writing with 0)
      // then loop over one less basis function

      // This is a possible issue if you need sensitivity to coordinates (you will need to
      // change basis and then use the product rule!)
      dof_ip(cell,pt) = dof_basis(cell, 0).val() * Value::eval(basis(cell, 0, pt));
      dof_ip(cell,pt).fastAccessDx(offsets(0)) = dof_basis(cell, 0).fastAccessDx(offsets(0)) * Value::eval(basis(cell, 0, pt));

      for (int bf=1; bf<numFields; bf++) {
        dof_ip(cell,pt).val() += dof_basis(cell, bf).val() * Value::eval(basis(cell, bf, pt));
        dof_ip(cell,pt).fastAccessDx(offsets(bf)) += dof_basis(cell, bf).fastAccessDx(offsets(bf)) * Value::eval(basis(cell, bf, pt));
      }
    }
  }
};

}

}

#endif
