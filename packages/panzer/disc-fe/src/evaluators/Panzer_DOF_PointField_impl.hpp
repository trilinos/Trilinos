// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_POINT_FIELD_DECL_HPP
#define PANZER_DOF_POINT_FIELD_DECL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template <typename EvalT, typename TRAITST>
void DOF_PointField<EvalT,TRAITST>::initialize(const std::string & fieldName,
                                                const PureBasis & fieldBasis,
                                                const std::string & coordinateName,
                                                const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                                                const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                                                const std::string & postfixFieldName)
{
  intrepidBasis = fieldBasis.getIntrepid2Basis();

  int cellCount = fieldBasis.functional->extent(0);
  int coeffCount = fieldBasis.functional->extent(1);
  int pointCount = coordLayout->extent(0);
  int dimCount = coordLayout->extent(1);

  Teuchos::RCP<PHX::DataLayout> basisLayout = fieldBasis.functional;

  coordinates = PHX::MDField<const ScalarT,Point,Dim>(coordinateName,coordLayout);
  dof_coeff = PHX::MDField<const ScalarT>(fieldName,basisLayout);
  dof_field = PHX::MDField<ScalarT>(fieldName+postfixFieldName,quadLayout);

  this->addDependentField(coordinates);
  this->addDependentField(dof_coeff);
  this->addEvaluatedField(dof_field);

  // build data storage for temporary conversion
  basisRef    = Kokkos::DynRankView<double,PHX::Device>("basisRef",coeffCount,pointCount);
  basis       = Kokkos::DynRankView<double,PHX::Device>("basis",cellCount,coeffCount,pointCount);
  intrpCoords = Kokkos::DynRankView<double,PHX::Device>("intrpCoords",pointCount,dimCount);
  
  std::string n = "DOF_PointField: " + dof_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
void DOF_PointField<EvalT,TRAITST>::evaluateFields(typename TRAITST::EvalData workset)
{ 
  // Zero out arrays (intrepid does a sum! 1/17/2012)
  dof_field.deep_copy(ScalarT(0.0));

  // copy coordinates
  auto l_intrpCoords = PHX::as_view(intrpCoords);
  auto l_coordinates = coordinates.get_static_view();
  Kokkos::parallel_for("DOF PointFields", l_coordinates.extent_int(0), KOKKOS_LAMBDA (int i) {
     for (int j = 0; j < l_coordinates.extent_int(1); ++j)
      l_intrpCoords(i,j) = Sacado::scalarValue(l_coordinates(i,j));
    });

  if(workset.num_cells>0) {
    // evaluate at reference points
    intrepidBasis->getValues(basisRef, intrpCoords, Intrepid2::OPERATOR_VALUE);

    // transfer reference basis values to physical frame values
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::
      HGRADtransformVALUE(basis,basisRef);

    // evaluate function at specified points
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::
      evaluate(dof_field.get_view(),dof_coeff.get_view(),basis);
  }
  Kokkos::fence();
}

//**********************************************************************

}

#endif
