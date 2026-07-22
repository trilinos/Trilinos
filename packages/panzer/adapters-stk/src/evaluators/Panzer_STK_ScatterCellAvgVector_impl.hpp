// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SCATTER_CELL_AVG_VECTOR_IMPL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_VECTOR_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace panzer_stk {

template<typename EvalT, typename Traits>
ScatterCellAvgVector<EvalT, Traits>::
ScatterCellAvgVector(
  const Teuchos::ParameterList& p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

  std::string scatterName = p.get<std::string>("Scatter Name");

  if (p.isParameter("Variable Scale Factors Map"))
    varScaleFactors_ = p.get<Teuchos::RCP<std::map<std::string,double>>>("Variable Scale Factors Map");

  const std::vector<std::string> & names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::IntegrationRule> intRule =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // build dependent fields
  scatterFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd)
  {
    scatterFields_[fd] = PHX::MDField<const ScalarT,Cell,Point,Dim>(names[fd],intRule->dl_vector);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": PanzerSTK::ScatterCellAvgVectors");
}


template<typename EvalT, typename Traits>
void
ScatterCellAvgVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* d */,
  PHX::FieldManager<Traits>& /* fm */)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd)
  {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<double>(stk::topology::ELEMENT_RANK, fieldName);
  }
}


template<typename EvalT, typename Traits>
void
ScatterCellAvgVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  panzer::MDFieldArrayFactory af("",true);

  // for convenience pull out some objects from workset
  const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;
  std::string blockId = this->wda(workset).block_id;
  std::string d_mod[3] = {"X","Y","Z"};

  // loop over the number of vector fields requested for exodus output
  for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++)
  {
    PHX::MDField<const ScalarT,panzer::Cell,panzer::Point,panzer::Dim> & field = scatterFields_[fieldIndex];
    std::string fieldName = field.fieldTag().name();
    const int numCells = workset.num_cells;
    const int numPoints = field.extent(1);
    const int numDims = field.extent(2);

    for (int dim = 0; dim < numDims; dim++)
    {
      // std::vector<double> average(numCells,0.0);
      PHX::MDField<double,panzer::Cell,panzer::NODE> average = af.buildStaticArray<double,panzer::Cell,panzer::NODE>("",numCells,1);

      // write to double field
      Kokkos::parallel_for("ScatterCellAvgVector",numCells,KOKKOS_LAMBDA(const int i){
        average(i,0) = 0.0;
        for(int j = 0; j < numPoints; j++) { // loop over IPs
          average(i,0) += Sacado::scalarValue(field(i,j,dim));
        }
        average(i,0) /= numPoints;
      });
      double scalef = 1.0;

      if (!varScaleFactors_.is_null())
      {
        std::map<std::string,double> *tmp_sfs = varScaleFactors_.get();
        if(tmp_sfs->find(fieldName) != tmp_sfs->end())
          scalef = (*tmp_sfs)[fieldName];
      }

      PHX::Device().fence();
      mesh_->setCellFieldData(fieldName+d_mod[dim],blockId,localCellIds,average.get_view(),scalef);
    }
  }
}

} // end panzer_stk

#endif
