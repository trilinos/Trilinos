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

#ifndef PANZER_DOF_CURL_IMPL_HPP
#define PANZER_DOF_CURL_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DOFCurl,p) :
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the curl operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsCurl(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" does not support CURL");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" in DOF Curl should require orientations. So we are throwing.");

  // build dof_curl
  basis_dimension = basis->dimension();
  if(basis_dimension==2)
     dof_curl = PHX::MDField<ScalarT>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );
  else if(basis_dimension==3)
     dof_curl = PHX::MDField<ScalarT>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector );
  else
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFCurl only works for 2D and 3D basis functions"); } 

  // add to evaluation graph
  this->addEvaluatedField(dof_curl);
  this->addDependentField(dof_value);
  
  std::string n = "DOFCurl: " + dof_curl.fieldTag().name() + " ("+PHX::TypeString<EvalT>::value+")";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOFCurl,sd,fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_curl,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOFCurl,workset)
{ 
  if(workset.num_cells>0) {
    // evaluate at quadrature points
    if(basis_dimension==3) {
      const Intrepid::FieldContainer<double> & curl_basis = workset.bases[basis_index]->curl_basis;

      int numCells  = workset.num_cells;
      int numFields = curl_basis.dimension(1);
      int numPoints = curl_basis.dimension(2);
      int spaceDim  = curl_basis.dimension(3);

      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d=0; d<spaceDim; d++) {
            // first initialize to the right thing (prevents over writing with 0)
            // then loop over one less basis function
            ScalarT & curl = dof_curl(cell,pt,d);
            curl = dof_value(cell, 0) * curl_basis(cell, 0, pt, d);
            for (int bf=1; bf<numFields; bf++)
              curl += dof_value(cell, bf) * curl_basis(cell, bf, pt, d);
          }
        }
      }

    }
    else {
      // Zero out arrays
      for (int i = 0; i < dof_curl.size(); ++i)
        dof_curl[i] = 0.0;

      Intrepid::FunctionSpaceTools::evaluate<ScalarT>(dof_curl,dof_value,workset.bases[basis_index]->curl_basis);
    }
  }
}

//**********************************************************************

}

#endif
