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
{ 
  dof_name = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  Teuchos::RCP<PHX::DataLayout> basis_layout  = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout  = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());

  result = PHX::MDField<ScalarT,Cell,BASIS>(dof_name,basis_layout);
  this->addEvaluatedField(result);

  normals = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Normals",vector_layout);
  this->addDependentField(normals);

  vector_values.resize(1);
  vector_values[0] = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name+"_Vector",vector_layout);
  this->addDependentField(vector_values[0]);
  
  this->setName("Project To Faces");

  num_faces  = result.extent(1);
  num_dim  = vector_values[0].extent(2);
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  num_faces  = result.extent(1);
  num_dim  = vector_values[0].extent(2);

  TEUCHOS_ASSERT(result.extent(1) == normals.extent(1));
  TEUCHOS_ASSERT(vector_values[0].extent(2) == normals.extent(2));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  // Restricting HDiv field, multiplied by the normal to the faces, into HVol on the faces.
  // This code assumes affine mapping and the projection into 1 quadrature point for each face,
  // which is identified with the face. This makes sense only for low order bases, for which
  // HVol is constant

  //TODO: make this work w/ high order basis
  const int intDegree = basis->order();
  TEUCHOS_ASSERT(intDegree == 1);
    
  // TODO: parallel for.
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int p = 0; p < num_faces; ++p) {
      result(cell,p) = ScalarT(0.0);
      for (int dim = 0; dim < num_dim; ++dim)
        result(cell,p) += vector_values[0](cell,p,dim) * normals(cell,p,dim);
    }
  }
}


#endif
