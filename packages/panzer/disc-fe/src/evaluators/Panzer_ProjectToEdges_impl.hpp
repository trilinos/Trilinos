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

#ifndef PANZER_PROJECT_TO_EDGES_IMPL_HPP
#define PANZER_PROJECT_TO_EDGES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::ProjectToEdges<EvalT, Traits>::
ProjectToEdges(
  const Teuchos::ParameterList& p)
{ 
  dof_name = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  quad_degree = 0;
  if(p.isType<int>("Quadrature Order"))
    quad_degree = p.get<int>("Quadrature Order");
    
  Teuchos::RCP<PHX::DataLayout> basis_layout  = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());

  result = PHX::MDField<ScalarT,Cell,BASIS>(dof_name,basis_layout);
  this->addEvaluatedField(result);

  tangents = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Tangents",vector_layout);
  this->addDependentField(tangents);

  if(quad_degree > 0){
    const shards::CellTopology & parentCell = *basis->getCellTopology();                                                                                    
    Intrepid2::DefaultCubatureFactory quadFactory;                         
    Teuchos::RCP< Intrepid2::Cubature<PHX::exec_space,double,double> > quadRule                    
      = quadFactory.create<PHX::exec_space,double,double>(parentCell.getCellTopologyData(1,0), quad_degree);
    int numQPoints = quadRule->getNumPoints(); 
 
    vector_values.resize(numQPoints);
    for (int qp(0); qp < numQPoints; ++qp)
    {
      vector_values[qp] = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Vector"+"_qp_"+std::to_string(qp),vector_layout);
      this->addDependentField(vector_values[qp]);
    }

    // setup the orientation field
    std::string orientationFieldName = basis->name() + " Orientation";
    dof_orientation = PHX::MDField<const ScalarT,Cell,NODE>(orientationFieldName,basis_layout);
    this->addDependentField(dof_orientation);

    gatherFieldTangents = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name+"_Tangents",basis->functional_grad);
    this->addEvaluatedField(gatherFieldTangents);

  } else {
    vector_values.resize(1);
    vector_values[0] = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Vector",vector_layout);
    this->addDependentField(vector_values[0]);
  }

  this->setName("Project To Edges");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToEdges<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  d, 
		      PHX::FieldManager<Traits>& fm)
{
  orientations = d.orientations_;

  // setup the field data object
  this->utils.setFieldData(result,fm);
  for(unsigned qp = 0; qp < vector_values.size(); ++qp)
    this->utils.setFieldData(vector_values[qp],fm);
  this->utils.setFieldData(tangents,fm);

  if(quad_degree > 0){
    this->utils.setFieldData(dof_orientation,fm);
    this->utils.setFieldData(gatherFieldTangents,fm);
  }

  num_pts = vector_values[0].extent(1);
  num_dim = vector_values[0].extent(2);

  TEUCHOS_ASSERT(vector_values[0].extent(1) == tangents.extent(1));
  TEUCHOS_ASSERT(vector_values[0].extent(2) == tangents.extent(2));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToEdges<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const shards::CellTopology & parentCell = *basis->getCellTopology();
  const int intDegree = basis->order();
  TEUCHOS_ASSERT(intDegree == 1);
  Intrepid2::DefaultCubatureFactory quadFactory;
  Teuchos::RCP<Intrepid2::Cubature<PHX::exec_space,double,double>> edgeQuad;

  // One point quadrature if higher order quadrature not requested
  if (quad_degree == 0){

    // this should be edge cubature and collecting weights are always 2.
    // // Collect the reference edge information. For now, do nothing with the quadPts.
    // const unsigned num_edges = parentCell.getEdgeCount();
    // std::vector<double> refEdgeWt(num_edges, 0.0);
    // for (unsigned e=0; e<num_edges; e++) {
    //   edgeQuad = quadFactory.create<PHX::exec_space,double,double>(parentCell.getCellTopologyData(1,e), intDegree);
    //   const int numQPoints = edgeQuad->getNumPoints();
    //   Kokkos::DynRankView<double,PHX::Device> quadWts("quadWts",numQPoints);
    //   Kokkos::DynRankView<double,PHX::Device> quadPts("quadPts",numQPoints,num_dim);
    //   edgeQuad->getCubature(quadPts,quadWts);
    //   for (int q=0; q<numQPoints; q++)
    //     refEdgeWt[e] += quadWts(q);
    // }

    Kokkos::DynRankView<double,PHX::Device> v0("v0", num_dim), v1("v1", num_dim);
    const unsigned num_edges = parentCell.getEdgeCount();
    std::vector<double> refEdgeWt(num_edges, 0.0);
    for (unsigned e=0; e<num_edges; e++) {
      const auto v0_id = parentCell.getNodeMap(1, e, 0);
      const auto v1_id = parentCell.getNodeMap(1, e, 1);
      Intrepid2::CellTools<PHX::exec_space>::getReferenceVertex(v0, parentCell, v0_id);
      Intrepid2::CellTools<PHX::exec_space>::getReferenceVertex(v1, parentCell, v1_id);
      
      double norm = 0.0;
      for (int d=0;d<num_dim;++d)
        norm += (v0(d) - v1(d))*(v0(d) - v1(d));
      
      refEdgeWt[e] = sqrt(norm);
    }

    // Loop over the edges of the workset cells.
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int p = 0; p < num_pts; ++p) {
        result(cell,p) = ScalarT(0.0);
        for (int dim = 0; dim < num_dim; ++dim)
          result(cell,p) += vector_values[0](cell,p,dim) * tangents(cell,p,dim);
        result(cell,p) *= refEdgeWt[p];
      }
    }

  } else {

    TEUCHOS_ASSERT(false); // this doesn't work since we modified the way orientations are handled

    PHX::MDField<double,Cell,panzer::NODE,Dim> vertex_coords = workset.cell_vertex_coordinates;
    int subcell_dim = 1;

    // to compute tangents at qps (copied from GatherTangents)
    int numEdges = gatherFieldTangents.extent(1);
    Kokkos::DynRankView<ScalarT,PHX::Device> refEdgeTan = Kokkos::createDynRankView(gatherFieldTangents.get_static_view(),"refEdgeTan",numEdges,num_dim);
    for(int i=0;i<numEdges;i++) {
      Kokkos::DynRankView<double,PHX::Device> refEdgeTan_local("refEdgeTan_local",num_dim);
      Intrepid2::CellTools<PHX::exec_space>::getReferenceEdgeTangent(refEdgeTan_local, i, parentCell);

      for(int d=0;d<num_dim;d++)
        refEdgeTan(i,d) = refEdgeTan_local(d);
    }

    // Loop over the faces of the workset cells
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {

      // get nodal coordinates for this cell 
      Kokkos::DynRankView<double,PHX::Device> physicalNodes("physicalNodes",1,vertex_coords.extent(1),num_dim);
      for (int point(0); point < vertex_coords.extent_int(1); ++point)
      {
        for (int ict(0); ict < num_dim; ict++)
          physicalNodes(0, point, ict) = vertex_coords(cell, point, ict);
      }

      // loop over edges
      for (int p = 0; p < num_pts; ++p){
        result(cell,p) = ScalarT(0.0);

        // get quad weights/pts on reference 2d cell
        const shards::CellTopology & subcell = parentCell.getCellTopologyData(subcell_dim,p);     
        edgeQuad = quadFactory.create<PHX::exec_space,double,double>(subcell, quad_degree);
        TEUCHOS_ASSERT(
          edgeQuad->getNumPoints() == static_cast<int>(vector_values.size()));
        Kokkos::DynRankView<double,PHX::Device> quadWts("quadWts",edgeQuad->getNumPoints());
        Kokkos::DynRankView<double,PHX::Device> quadPts("quadPts",edgeQuad->getNumPoints(),subcell_dim);
        edgeQuad->getCubature(quadPts,quadWts);

        // map 1d quad pts to reference cell (3d)
        Kokkos::DynRankView<double,PHX::Device> refQuadPts("refQuadPts",edgeQuad->getNumPoints(),num_dim);
        Intrepid2::CellTools<PHX::exec_space>::mapToReferenceSubcell(refQuadPts, quadPts, subcell_dim, p, parentCell);


        // Calculate side jacobian
        Kokkos::DynRankView<double,PHX::Device> jacobianSide("jacobianSide", 1, edgeQuad->getNumPoints(), num_dim, num_dim);
        Intrepid2::CellTools<PHX::exec_space>::setJacobian(jacobianSide, refQuadPts, physicalNodes, parentCell);

        // Calculate weighted measure at quadrature points
        Kokkos::DynRankView<double,PHX::Device> weighted_measure("weighted_measure",1,edgeQuad->getNumPoints());
        Kokkos::DynRankView<double,PHX::Device> scratch_space("scratch_space",jacobianSide.span());
        Intrepid2::FunctionSpaceTools<PHX::exec_space>::computeEdgeMeasure(weighted_measure, jacobianSide, quadWts, p, parentCell,scratch_space);

        // loop over quadrature points
        for (int qp = 0; qp < edgeQuad->getNumPoints(); ++qp) {

          // get normal vector at quad points
          std::vector<ScalarT> edgeTan(num_dim);
          for(int i = 0; i < 3; i++) {
            edgeTan[i] = (Sacado::ScalarValue<ScalarT>::eval(jacobianSide(0,qp,i,0)*refEdgeTan(p,0))
                        + Sacado::ScalarValue<ScalarT>::eval(jacobianSide(0,qp,i,1)*refEdgeTan(p,1))
                        + Sacado::ScalarValue<ScalarT>::eval(jacobianSide(0,qp,i,2)*refEdgeTan(p,2)))
                        * dof_orientation(cell,p);
          }

          // compute the magnitude of the tangent vector
          ScalarT tnorm(0.0);
          for(int dim = 0; dim < num_dim; ++dim){
            tnorm += edgeTan[dim]*edgeTan[dim];
          }
          tnorm = std::sqrt(tnorm);

          // integrate vector dot t
          // normalize t since jacobian information is factored into both weighted measure and tangent
          for (int dim = 0; dim < num_dim; ++dim)
            result(cell,p) += weighted_measure(0,qp) * vector_values[qp](cell,p,dim) * edgeTan[dim] / tnorm;

        }
      }

    }

  }


}

#endif
