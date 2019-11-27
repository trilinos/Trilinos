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

#ifndef PANZER_PROJECT_TO_FACES_IMPL_HPP
#define PANZER_PROJECT_TO_FACES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <cstring>

template<typename EvalT,typename Traits>
panzer::ProjectToFaces<EvalT, Traits>::
ProjectToFaces(const Teuchos::ParameterList& p)
  : quad_degree(1),
    use_fast_method_on_rectangular_hex_mesh(false)
{ 
  dof_name = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
 
  if(p.isType<int>("Quadrature Order"))
    quad_degree = p.get<int>("Quadrature Order");

  if(p.isType<bool>("Use Fast Method for Rectangular Hex Mesh"))
    use_fast_method_on_rectangular_hex_mesh = p.get<bool>("Use Fast Method for Rectangular Hex Mesh");

  Teuchos::RCP<PHX::DataLayout> basis_layout  = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout  = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());

  result = PHX::MDField<ScalarT,Cell,BASIS>(dof_name,basis_layout);
  this->addEvaluatedField(result);

  normals = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Normals",vector_layout);
  this->addDependentField(normals);

  if(not use_fast_method_on_rectangular_hex_mesh){
    const shards::CellTopology & parentCell = *basis->getCellTopology();                                                                                    
    Intrepid2::DefaultCubatureFactory quadFactory;
    Teuchos::RCP< Intrepid2::Cubature<PHX::exec_space,double,double> > quadRule                    
      = quadFactory.create<PHX::exec_space,double,double>(parentCell.getCellTopologyData(2,0), quad_degree);
    int numQPoints = quadRule->getNumPoints(); 
 
    vector_values.resize(numQPoints);
    for (int qp(0); qp < numQPoints; ++qp)
    {
      vector_values[qp] = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Vector"+"_qp_"+std::to_string(qp),vector_layout);
      this->addDependentField(vector_values[qp]);
    }

    gatherFieldNormals = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name+"_Normals",basis->functional_grad);
    this->addEvaluatedField(gatherFieldNormals);

  } else {
    TEUCHOS_ASSERT(quad_degree == 1); // One pt quadrature for fast method
    TEUCHOS_ASSERT(std::strstr(basis->getCellTopology()->getBaseName(),"Hexahedron") != nullptr);

    vector_values.resize(1);
    vector_values[0] = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Vector",vector_layout);
    this->addDependentField(vector_values[0]);
  }

  if (use_fast_method_on_rectangular_hex_mesh)
    this->setName("Project To Faces (Fast Rectangular Hex)");
  else
    this->setName("Project To Faces");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  orientations = d.orientations_;

  num_pts  = result.extent(1);
  num_dim  = vector_values[0].extent(2);

  TEUCHOS_ASSERT(result.extent(1) == normals.extent(1));
  TEUCHOS_ASSERT(vector_values[0].extent(2) == normals.extent(2));


  // Cache off vectors and keep intrepid2 function calls on host (by
  // switching memory space to ProjectionSpace). This is done to
  // prevent many small UVM allocations of temporaries in the
  // intrepid2 functions.
  if (use_fast_method_on_rectangular_hex_mesh) {
    // Reminder to fix the fast method in the future 
    // quadWts = Kokkos::DynRankView<double,PHX::Device>("quadWts",numQPoints);
    // quadPts = Kokkos::DynRankView<double,PHX::Device>("quadPts",numQPoints,num_dim);
  } else {
    const shards::CellTopology & parentCell = *basis->getCellTopology();
    const int cellDim = parentCell.getDimension();
    using HostDRVType = Kokkos::DynRankView<ScalarT,ProjectionSpace>;
    refEdges = Kokkos::createDynRankViewWithType<HostDRVType>(result.get_static_view(),"ref_edges", 2, cellDim);
    phyEdges = Kokkos::createDynRankViewWithType<HostDRVType>(result.get_static_view(),"phy_edges", 2, cellDim);

    PHX::MDField<double,Cell,panzer::NODE,Dim> vertex_coords = (*d.worksets_)[0].cell_vertex_coordinates;
    physicalNodes = Kokkos::DynRankView<double,ProjectionSpace>("physicalNodes",1,vertex_coords.extent(1),num_dim);

    const int subcell_dim = 2;
    Intrepid2::DefaultCubatureFactory quadFactory;
    int maxNumFaceQP = 0;
    faceQuads.resize(num_pts);
    for (int p = 0; p < num_pts; ++p){
      const shards::CellTopology & subcell = parentCell.getCellTopologyData(subcell_dim,p);
      faceQuads[p] = quadFactory.create<ProjectionSpace::execution_space,double,double>(subcell, quad_degree);
      TEUCHOS_ASSERT(faceQuads[p]->getNumPoints() == static_cast<int>(vector_values.size()));
      maxNumFaceQP = std::max(maxNumFaceQP,faceQuads[p]->getNumPoints());
    }
    quadWts = Kokkos::DynRankView<double,ProjectionSpace>("quadWts",maxNumFaceQP);
    quadPts = Kokkos::DynRankView<double,ProjectionSpace>("quadPts",maxNumFaceQP,subcell_dim);
    refQuadPts = Kokkos::DynRankView<double,ProjectionSpace>("refQuadPts",maxNumFaceQP,num_dim);
    jacobianSide = Kokkos::DynRankView<double,ProjectionSpace>("jacobianSide",1,maxNumFaceQP,num_dim,num_dim);
    weighted_measure = Kokkos::DynRankView<double,ProjectionSpace>("weighted_measure",1,maxNumFaceQP);
    scratch_space = Kokkos::DynRankView<double,ProjectionSpace>("scratch_space",jacobianSide.span());
  }
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  // The coefficients being calculated here in the projection to the face basis
  // are defined as the integral over the face of the field dotted with the face
  // normal vector. For a first-order face basis, single point integration is
  // adequate, so the cubature here just provides the proper weighting.
  // For higher order, a distinction between "cell" and Gauss points will need
  // to be made so the field is appropriately projected.
  const shards::CellTopology & parentCell = *basis->getCellTopology();
  Intrepid2::DefaultCubatureFactory quadFactory;

  // Fast Method: One point quadrature on hex mesh. Assumes that the
  // face area can be computed from two edge normals. This is only
  // true if the mesh is square or rectangular. Will not work for
  // paved meshes.
  if (use_fast_method_on_rectangular_hex_mesh){
    Teuchos::RCP< Intrepid2::Cubature<PHX::exec_space,double,double> > faceQuad;

    // Collect the reference face information. For now, do nothing with the quadPts.
    const unsigned num_faces = parentCell.getFaceCount();
    std::vector<double> refFaceWt(num_faces, 0.0);
    for (unsigned f=0; f<num_faces; f++) {
      faceQuad = quadFactory.create<PHX::exec_space,double,double>(parentCell.getCellTopologyData(2,f), 1);
      const int numQPoints = faceQuad->getNumPoints();
      Kokkos::DynRankView<double,PHX::Device> quadWtsFast("quadWts",numQPoints);
      Kokkos::DynRankView<double,PHX::Device> quadPtsFast("quadPts",numQPoints,num_dim);
      faceQuad->getCubature(quadPtsFast,quadWtsFast);
      for (int q=0; q<numQPoints; q++)
        refFaceWt[f] += quadWtsFast(q);
    }
    
    // Loop over the faces of the workset cells.
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int p = 0; p < num_pts; ++p) {
        result(cell,p) = ScalarT(0.0);
        for (int dim = 0; dim < num_dim; ++dim)
          result(cell,p) += vector_values[0](cell,p,dim) * normals(cell,p,dim);
        result(cell,p) *= refFaceWt[p];
      }
    }

  } else {
    PHX::MDField<double,Cell,panzer::NODE,Dim> vertex_coords = workset.cell_vertex_coordinates;
    const int subcell_dim = 2;
    const int numFaces = Teuchos::as<int>(parentCell.getFaceCount());
    const WorksetDetails & details = workset;

    // Loop over the faces of the workset cells
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {

      // get nodal coordinates for this cell 
      for (int point(0); point < vertex_coords.extent_int(1); ++point)
      {
        for (int ict(0); ict < num_dim; ict++)
           physicalNodes(0,point,ict) = vertex_coords(cell,point,ict);
      }

      int faceOrts[6] = {};
      orientations->at(details.cell_local_ids[cell]).getFaceOrientation(faceOrts, numFaces);

      // loop over faces
      for (int p = 0; p < num_pts; ++p){
        result(cell,p) = ScalarT(0.0);

        auto ortEdgeTan_U = Kokkos::subview(refEdges, 0, Kokkos::ALL());
        auto ortEdgeTan_V = Kokkos::subview(refEdges, 1, Kokkos::ALL());

        // Apply parent cell Jacobian to ref. edge tangent
        Intrepid2::Orientation::getReferenceFaceTangents(ortEdgeTan_U,
                                                         ortEdgeTan_V,
                                                         p,
                                                         parentCell,
                                                         faceOrts[p]);

        // get quad weights/pts on reference 2d cell
        const shards::CellTopology & subcell = parentCell.getCellTopologyData(subcell_dim,p);     
        const auto& faceQuad = faceQuads[p];
        faceQuad->getCubature(quadPts,quadWts);

        // map 2d quad pts to reference cell (3d)
        Intrepid2::CellTools<ProjectionSpace::execution_space>::mapToReferenceSubcell(refQuadPts, quadPts, subcell_dim, p, parentCell);

        // Calculate side jacobian
        Intrepid2::CellTools<ProjectionSpace::execution_space>::setJacobian(jacobianSide, refQuadPts, physicalNodes, parentCell);

        // Calculate weighted measure at quadrature points
        Intrepid2::FunctionSpaceTools<ProjectionSpace::execution_space>::computeFaceMeasure(weighted_measure, jacobianSide, quadWts, p, parentCell, scratch_space);

        // loop over quadrature points
        for (int qp = 0; qp < faceQuad->getNumPoints(); ++qp) {

          auto phyEdgeTan_U = Kokkos::subview(phyEdges, 0, Kokkos::ALL());
          auto phyEdgeTan_V = Kokkos::subview(phyEdges, 1, Kokkos::ALL());
          auto J = Kokkos::subview(jacobianSide, 0, qp, Kokkos::ALL(), Kokkos::ALL());
  
          Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan_U, J, ortEdgeTan_U);            
          Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan_V, J, ortEdgeTan_V);            

          // normal = TanU x TanV
          ScalarT normal[3];
          normal[0] = (phyEdgeTan_U(1)*phyEdgeTan_V(2) - phyEdgeTan_U(2)*phyEdgeTan_V(1));
          normal[1] = (phyEdgeTan_U(2)*phyEdgeTan_V(0) - phyEdgeTan_U(0)*phyEdgeTan_V(2));
          normal[2] = (phyEdgeTan_U(0)*phyEdgeTan_V(1) - phyEdgeTan_U(1)*phyEdgeTan_V(0));

          // compute the magnitude of the normal vector
          ScalarT nnorm(0.0);
          for(int dim = 0; dim < num_dim; ++dim){
            nnorm += normal[dim]*normal[dim];
          }
          nnorm = std::sqrt(nnorm);

          // integrate vector dot n
          // normalize n since jacobian information is factored into both weighted measure and normal
          for (int dim = 0; dim < num_dim; ++dim)
            result(cell,p) += weighted_measure(0,qp) * vector_values[qp](cell,p,dim) * normal[dim] / nnorm;
        }
      }

    }

  } // end else (high order quad)

}

#endif
