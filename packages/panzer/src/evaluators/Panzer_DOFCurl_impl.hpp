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

namespace {

//**********************************************************************
template<typename ScalarT>                   
void evaluateCurl_withSens(int numCells,
                           int basis_dimension,
                           PHX::MDField<ScalarT> & dof_curl, 
                           PHX::MDField<ScalarT,Cell,Point> & dof_value,
                           const Intrepid::FieldContainer<double> & curl_basis)
{ 
  if(numCells>0) {
    // evaluate at quadrature points
    if(basis_dimension==3) {

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

      Intrepid::FunctionSpaceTools::evaluate<ScalarT>(dof_curl,dof_value,curl_basis);
    }
  }
}

}

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename Traits>                   
DOFCurl<EvalT, Traits>::
DOFCurl(const Teuchos::ParameterList & p) :
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
template<typename EvalT, typename Traits>                   
void DOFCurl<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_curl,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
template<typename EvalT, typename Traits>                   
void DOFCurl<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  evaluateCurl_withSens<ScalarT>(workset.num_cells,basis_dimension,dof_curl,dof_value,workset.bases[basis_index]->curl_basis);
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename Traits>                   
DOFCurl<panzer::Traits::Jacobian, Traits>::
DOFCurl(const Teuchos::ParameterList & p) :
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // do you specialize because you know where the basis functions are and can
  // skip a large number of AD calculations?
  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");
    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix

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
  
  std::string n = "DOFCurl: " + dof_curl.fieldTag().name() + " ("+PHX::TypeString<panzer::Traits::Jacobian>::value+")";
  this->setName(n);
}

//**********************************************************************
template<typename Traits>                   
void DOFCurl<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_curl,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

template<typename Traits>                   
void DOFCurl<panzer::Traits::Jacobian,Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  if(!accelerate_jacobian) {
    // do the case where we use the AD types to determine the derivatives
    evaluateCurl_withSens<ScalarT>(workset.num_cells,basis_dimension,dof_curl,dof_value,workset.bases[basis_index]->curl_basis);
    return;
  }

  if(workset.num_cells>0) {
    // evaluate at quadrature points
    if(basis_dimension==3) {
      const Intrepid::FieldContainer<double> & curl_basis = workset.bases[basis_index]->curl_basis;

      int numCells  = workset.num_cells;
      int numFields = curl_basis.dimension(1);
      int numPoints = curl_basis.dimension(2);
      int spaceDim  = curl_basis.dimension(3);

      int fadSize = dof_value(0,0).size(); // this is supposed to be fast
                                           // so assume that everything is the
                                           // same size!

      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d=0; d<spaceDim; d++) {
            // first initialize to the right thing (prevents over writing with 0)
            // then loop over one less basis function
            ScalarT & curl = dof_curl(cell,pt,d);
            curl = ScalarT(fadSize, dof_value(cell, 0).val() * curl_basis(cell, 0, pt, d));
            curl.fastAccessDx(offsets[0]) = dof_value(cell, 0).fastAccessDx(offsets[0]) * curl_basis(cell, 0, pt, d);
            for (int bf=1; bf<numFields; bf++) {
              curl.val() += dof_value(cell, bf).val() * curl_basis(cell, bf, pt, d);
              curl.fastAccessDx(offsets[bf]) += dof_value(cell, bf).fastAccessDx(offsets[bf]) * curl_basis(cell, bf, pt, d);
            }
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

}

#endif
