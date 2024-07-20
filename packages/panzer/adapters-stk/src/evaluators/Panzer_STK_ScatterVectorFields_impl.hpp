// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SCATTER_VECTOR_FIELDS_IMPL_HPP
#define PANZER_STK_SCATTER_VECTOR_FIELDS_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

template<typename EvalT, typename Traits>
ScatterVectorFields<EvalT, Traits>::
ScatterVectorFields(
  const Teuchos::ParameterList& p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
   TEUCHOS_ASSERT(false);
}

template <typename EvalT,typename TraitsT>
ScatterVectorFields<EvalT,TraitsT>::
ScatterVectorFields(const std::string & scatterName,
              const Teuchos::RCP<STK_Interface> mesh,
              const Teuchos::RCP<const panzer::PointRule> & pointRule,
              const std::vector<std::string> & names,
              const std::vector<double> & scaling)
   : mesh_(mesh)
   , scaling_(scaling)
{
  using panzer::Cell;
  using panzer::IP;
  using panzer::Dim;

  spatialDimension_ = pointRule->spatial_dimension;

  // this evaluator assumes you are evaluating at the cell centroid only
  TEUCHOS_ASSERT(pointRule->num_points==1);

  // build dependent fields
  names_ = names;
  scatterFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] =
      PHX::MDField<const ScalarT,Cell,IP,Dim>(names_[fd]+"_"+pointRule->getName(),pointRule->dl_vector);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Vector Fields");
}

template<typename EvalT, typename Traits>
void
ScatterVectorFields<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{
   TEUCHOS_ASSERT(false);
}

template < >
void ScatterVectorFields<panzer::Traits::Residual,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData workset)
{
  panzer::MDFieldArrayFactory af("",true);

  std::vector<std::string> dimStrings(3);
  dimStrings[0] = "X";
  dimStrings[1] = "Y";
  dimStrings[2] = "Z";

  // for convenience pull out some objects from workset
  const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;
  std::string blockId = this->wda(workset).block_id;

  for(int d=0;d<spatialDimension_;d++) {
    for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++) {
      PHX::MDField<const ScalarT,panzer::Cell,panzer::IP,panzer::Dim> & field = scatterFields_[fieldIndex];
      std::string fieldName = names_[fieldIndex]+dimStrings[d];

      PHX::MDField<double,panzer::Cell,panzer::NODE> cellValue
          = af.buildStaticArray<double,panzer::Cell,panzer::NODE>("",field.extent(0),1);

      // scaline field value only if the scaling parameter is specified, otherwise use 1.0
      double scaling = (scaling_.size()>0) ? scaling_[fieldIndex] : 1.0;

      auto cellValue_v = cellValue.get_static_view();
      auto field_v = field.get_static_view();
      Kokkos::parallel_for(field_v.extent(0), KOKKOS_LAMBDA (int i) {
	  cellValue_v(i,0) = field_v(i,0,d);
	});
      Kokkos::fence();

      // add in vector value at d^th dimension
      mesh_->setCellFieldData(fieldName,blockId,localCellIds,cellValue.get_view(),scaling);
    }
  }
}

} // end panzer_stk

#endif
