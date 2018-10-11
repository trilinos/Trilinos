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

#ifndef PANZER_DIRICHLET_RESIDUAL_FACEBASIS_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_FACEBASIS_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "Intrepid2_CellTools.hpp"

#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
DirichletResidual_FaceBasis<EvalT, Traits>::
DirichletResidual_FaceBasis(
  const Teuchos::ParameterList& p)
{
  std::string residual_name = p.get<std::string>("Residual Name");

  basis = p.get<Teuchos::RCP<const panzer::PureBasis> >("Basis");
  pointRule = p.get<Teuchos::RCP<const panzer::PointRule> >("Point Rule");

  std::string field_name = p.get<std::string>("DOF Name");
  std::string dof_name = field_name+"_"+pointRule->getName(); 
  std::string value_name = p.get<std::string>("Value Name");

  Teuchos::RCP<PHX::DataLayout> basis_layout = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout_dof    = pointRule->dl_vector;
  Teuchos::RCP<PHX::DataLayout> vector_layout_vector = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());
  TEUCHOS_ASSERT(basis_layout->extent(0)==vector_layout_dof->extent(0));
  TEUCHOS_ASSERT(basis_layout->extent(1)==vector_layout_dof->extent(1));
  TEUCHOS_ASSERT(basis->dimension()==vector_layout_dof->extent_int(2));
  TEUCHOS_ASSERT(vector_layout_vector->extent(0)==vector_layout_dof->extent(0));
  TEUCHOS_ASSERT(vector_layout_vector->extent(1)==vector_layout_dof->extent(1));
  TEUCHOS_ASSERT(vector_layout_vector->extent(2)==vector_layout_dof->extent(2));

  residual = PHX::MDField<ScalarT,Cell,BASIS>(residual_name, basis_layout);
  dof      = PHX::MDField<const ScalarT,Cell,Point,Dim>(dof_name, vector_layout_dof);
  value    = PHX::MDField<const ScalarT,Cell,Point,Dim>(value_name, vector_layout_vector);

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<double> (pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  constJac_ = pointValues.jac;
  this->addDependentField(constJac_);

  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Face Basis Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual_FaceBasis<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  orientations = sd.orientations_;
  this->utils.setFieldData(pointValues.jac,fm);
  faceNormal = Kokkos::createDynRankView(residual.get_static_view(),"faceNormal",dof.extent(0),dof.extent(1),dof.extent(2));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual_FaceBasis<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  // basic cell topology data
  const shards::CellTopology & parentCell = *basis->getCellTopology();
  const int cellDim = parentCell.getDimension();

  // face only, subcellOrd does not include edge count
  const int subcellDim = 2;
  const int subcellOrd = this->wda(workset).subcell_index;

  const int numFaces = parentCell.getSubcellCount(subcellDim);  
  const int numFaceDofs = dof.extent_int(1);

  TEUCHOS_ASSERT(cellDim == dof.extent_int(2));

  if(workset.num_cells<=0)
    return;
  else {
    Intrepid2::CellTools<PHX::exec_space>::getPhysicalFaceNormals(faceNormal,
                                                                  pointValues.jac.get_view(),
                                                                  subcellOrd,
                                                                  parentCell);

    const auto subcellTopo = shards::CellTopology(parentCell.getBaseCellTopologyData(subcellDim, subcellOrd));
    TEUCHOS_ASSERT(subcellTopo.getBaseKey() == shards::Triangle<>::key ||
                   subcellTopo.getBaseKey() == shards::Quadrilateral<>::key);

    const WorksetDetails & details = workset;

    int faceOrts[6] = {};
    for(index_t c=0;c<workset.num_cells;c++) {
      const auto ort = orientations->at(details.cell_local_ids[c]);
      ort.getFaceOrientation(faceOrts, numFaces); 

      // vertex count represent rotation count before it flips
      const double ortVal = faceOrts[subcellOrd] < static_cast<int>(subcellTopo.getVertexCount()) ? 1.0 : -1.0;
      for(int b=0;b<numFaceDofs;b++) {
        residual(c,b) = ScalarT(0.0);
        for(int d=0;d<cellDim;d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceNormal(c,b,d);
        residual(c,b) *= ortVal;
      } 
    }
  }
}

//**********************************************************************

}

#endif
