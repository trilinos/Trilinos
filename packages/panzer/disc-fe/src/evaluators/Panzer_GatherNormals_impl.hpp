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

#ifndef PANZER_GATHER_NORMALS_IMPL_HPP
#define PANZER_GATHER_NORMALS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::GatherNormals<EvalT, Traits>::
GatherNormals(
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
  dof_orientation = PHX::MDField<ScalarT,Cell,NODE>(orientationFieldName,basis_layout);

  // setup all basis fields that are required
  MDFieldArrayFactory af_pv(pointRule->getName()+"_");

  // setup all fields to be evaluated and constructed
  pointValues.setupArrays(pointRule,af_pv);

  // the field manager will allocate all of these field
  this->addDependentField(dof_orientation);
  this->addDependentField(pointValues.jac);

  gatherFieldNormals = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name+"_Normals",vector_layout_vector);
  this->addEvaluatedField(gatherFieldNormals);

  this->setName("Gather Normals");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherNormals<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  // setup the field data object
  this->utils.setFieldData(gatherFieldNormals,fm);
  this->utils.setFieldData(dof_orientation,fm);
  this->utils.setFieldData(pointValues.jac,fm);

  faceNormal = Kokkos::createDynRankView(gatherFieldNormals.get_static_view(),
					 "faceNormal",
					 gatherFieldNormals.dimension(0),
					 gatherFieldNormals.dimension(1),
					 gatherFieldNormals.dimension(2));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherNormals<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  if(workset.num_cells<=0)
    return;

  const shards::CellTopology & parentCell = *basis->getCellTopology();
  int cellDim = parentCell.getDimension();
  int numFaces = gatherFieldNormals.dimension(1);

  // Collect the tangents for the element faces in reference space.
  // These are scaled such that U x V returns a unit normal,
  // **contrary to the Intrepid documentation**.
  Kokkos::DynRankView<ScalarT,PHX::Device> refFaceTanU = Kokkos::createDynRankView(gatherFieldNormals.get_static_view(),"refFaceTanU",numFaces,cellDim);
  Kokkos::DynRankView<ScalarT,PHX::Device> refFaceTanV = Kokkos::createDynRankView(gatherFieldNormals.get_static_view(),"refFaceTanV",numFaces,cellDim);
  for(int i=0;i<numFaces;i++) {
    Kokkos::DynRankView<double,PHX::Device> refTanU = Kokkos::DynRankView<double,PHX::Device>("refTanU",cellDim);
    Kokkos::DynRankView<double,PHX::Device> refTanV = Kokkos::DynRankView<double,PHX::Device>("refTanV",cellDim);
    Intrepid2::CellTools<double>::getReferenceFaceTangents(refTanU, refTanV, i, parentCell);
    for(int d=0;d<cellDim;d++) {
      refFaceTanU(i,d) = refTanU(d);
      refFaceTanV(i,d) = refTanV(d);
    }
  }

  // The conversion to physical space must be done by converting each face tangent,
  // then computing the normal: see Intrepid2_CellToolsDef.hpp.
  // This code duplicates Intrepid2::getPhysicalFaceNormals to avoid converting local
  // data structures to and from Intrepid data structures.
  // Note that the magnitude of the normal is related to the area of the physical face.
  for(index_t c=0;c<workset.num_cells;c++) {
    for(int f = 0; f < numFaces; f++) {

      std::vector<double> faceTanU(3);
      std::vector<double> faceTanV(3);
      for(int i = 0; i < cellDim; i++) {
        faceTanU[i] = Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,0)*refFaceTanU(f,0))
                    + Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,1)*refFaceTanU(f,1))
                    + Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,2)*refFaceTanU(f,2));
        faceTanV[i] = Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,0)*refFaceTanV(f,0))
                    + Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,1)*refFaceTanV(f,1))
                    + Sacado::ScalarValue<ScalarT>::eval(pointValues.jac(c,f,i,2)*refFaceTanV(f,2));
      }

      // normal = TanU x TanV
      gatherFieldNormals(c,f,0) = (faceTanU[1]*faceTanV[2] - faceTanU[2]*faceTanV[1])*dof_orientation(c,f);
      gatherFieldNormals(c,f,1) = (faceTanU[2]*faceTanV[0] - faceTanU[0]*faceTanV[2])*dof_orientation(c,f);
      gatherFieldNormals(c,f,2) = (faceTanU[0]*faceTanV[1] - faceTanU[1]*faceTanV[0])*dof_orientation(c,f);

    }
  }

}

#endif
