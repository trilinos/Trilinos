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

#ifndef PANZER_GATHER_TANGENTS_IMPL_HPP
#define PANZER_GATHER_TANGENTS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_PureBasis.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::GatherTangents<EvalT, Traits>::
GatherTangents(
  const Teuchos::ParameterList& p)
{ 
  dof_name = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  pointRule = p.get<Teuchos::RCP<const PointRule> >("Point Rule");

  Teuchos::RCP<PHX::DataLayout> basis_layout         = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout_vector = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());

  // setup the orientation field
  std::string orientationFieldName = basis->name() + " Orientation";
  // setup all fields to be evaluated and constructed
  pointValues = panzer::PointValues2<ScalarT> (pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  constJac_ = pointValues.jac;
  this->addDependentField(constJac_);

  gatherFieldTangents = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name+"_Tangents",vector_layout_vector);
  this->addEvaluatedField(gatherFieldTangents);

  this->setName("Gather Tangents");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  orientations = d.orientations_;

  // setup the field data object
  this->utils.setFieldData(gatherFieldTangents,fm);
  this->utils.setFieldData(pointValues.jac,fm);
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  if(workset.num_cells<=0)
    return;
  else {
    const shards::CellTopology & parentCell = *basis->getCellTopology();
    int cellDim = parentCell.getDimension();
    int numEdges = gatherFieldTangents.extent(1);

    auto workspace = Kokkos::createDynRankView(gatherFieldTangents.get_static_view(),"workspace", 4, cellDim);

    const WorksetDetails & details = workset;

    const auto worksetJacobians = pointValues.jac.get_view();

    // Loop over workset faces and edge points
    for(index_t c=0;c<workset.num_cells;c++) {

      int edgeOrts[12] = {};
      orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);

      for(int pt = 0; pt < numEdges; pt++) {
        auto phyEdgeTan = Kokkos::subview(gatherFieldTangents.get_static_view(), c, pt, Kokkos::ALL());
        auto ortEdgeTan = Kokkos::subview(workspace, 1, Kokkos::ALL());

        // Apply parent cell Jacobian to ref. edge tangent
        Intrepid2::Orientation::getReferenceEdgeTangent(ortEdgeTan,
                                                        pt,
                                                        parentCell,
                                                        edgeOrts[pt]);
        
        auto J = Kokkos::subview(worksetJacobians, c, pt, Kokkos::ALL(), Kokkos::ALL());
        Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);            

      }// for pt
    }// for pCell

  }

}

#endif
