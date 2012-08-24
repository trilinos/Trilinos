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

#ifndef PANZER_DOF_IMPL_HPP
#define PANZER_DOF_IMPL_HPP

#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

namespace panzer {

//**********************************************************************

// This hides the evaluateDOF function outside of this file (if ETI is on)
namespace {

// A function for evaluating the DOFs. This is useful because 
// this code is needed twice, once for a DOF evaluator pulling from
// the workset and once for a DOF evaluator pulling from the field
// manager.
template <typename ScalarT,typename ArrayScalarT, typename ArrayT>
inline void evaluateDOF(PHX::MDField<ScalarT,Cell,Point> & dof_basis,PHX::MDField<ScalarT> & dof_ip,
                        PHX::MDField<ScalarT,Cell,BASIS> & dof_orientation,bool requires_orientation,
                        int num_cells,
                        panzer::BasisValues<ArrayScalarT,ArrayT> & basisValues)
{ 
  // Zero out arrays (intrepid does a sum! 1/17/2012)
  for (int i = 0; i < dof_ip.size(); ++i)
    dof_ip[i] = 0.0;

  if(num_cells>0) {
    if(requires_orientation) {
       // copy basis values, this will be modified by the orientation, so we don't
       // necessarily want to wipe them out
       Teuchos::Array<int> dimension; 
       for(int i=0;i<basisValues.basis.rank();i++)
          dimension.push_back(basisValues.basis.dimension(i));
       Intrepid::FieldContainer<ScalarT> bases(dimension);
       for(int i=0;i<bases.size();i++)
          bases[i] = basisValues.basis[i];

       // assign ScalarT "dof_orientation" to double "orientation"
       Intrepid::FieldContainer<double> orientation(dof_orientation.dimension(0),
                                                    dof_orientation.dimension(1));
       for(int i=0;i<dof_orientation.dimension(0);i++)
          for(int j=0;j<dof_orientation.dimension(1);j++)
             orientation(i,j) = Sacado::ScalarValue<ScalarT>::eval(dof_orientation(i,j));


       // make sure things are orientated correctly
       Intrepid::FunctionSpaceTools::
          applyFieldSigns<ScalarT>(bases,orientation);

       // evaluate at quadrature points
       Intrepid::FunctionSpaceTools::
         evaluate<ScalarT>(dof_ip,dof_basis,bases);
    }
    else // no orientation needed
       Intrepid::FunctionSpaceTools::
         evaluate<ScalarT>(dof_ip,dof_basis,basisValues.basis);
  }
}

}

//**********************************************************************
//* DOF evaluator
//**********************************************************************
PHX_EVALUATOR_CTOR(DOF,p) :
  dof_basis( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();
  requires_orientation = basis->requiresOrientations();

  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis())
     dof_ip = PHX::MDField<ScalarT>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
  else if(basis->isVectorBasis())
     dof_ip = PHX::MDField<ScalarT>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector);
  else
  { TEUCHOS_ASSERT(false); }

  this->addEvaluatedField(dof_ip);
  this->addDependentField(dof_basis);

  if(requires_orientation) {
     dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(p.get<std::string>("Name")+" Orientation", 
	                                                basis->functional);
     
     this->addDependentField(dof_orientation);
  }
  
  std::string n = "DOF: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOF,sd,fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  if(requires_orientation)
     this->utils.setFieldData(dof_orientation,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOF,workset)
{ 
  panzer::BasisValues<double,Intrepid::FieldContainer<double> > & basisValues = *workset.bases[basis_index];

  evaluateDOF(dof_basis,dof_ip,dof_orientation,requires_orientation,workset.num_cells,basisValues);
}

//**********************************************************************

//**********************************************************************
//* DOF_PointValues evaluator
//**********************************************************************
PHX_EVALUATOR_CTOR(DOF_PointValues,p)
{
  const std::string fieldName = p.get<std::string>("Name");
  basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  Teuchos::RCP<const PointRule> pointRule = p.get< Teuchos::RCP<const PointRule> >("Point Rule");
  requires_orientation = basis->requiresOrientations();

  const std::string evalName = fieldName+"_"+pointRule->getName();

  dof_basis = PHX::MDField<ScalarT,Cell,Point>(fieldName, basis->functional);

  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis())
     dof_ip = PHX::MDField<ScalarT>(
                evalName,
     	        pointRule->dl_scalar);
  else if(basis->isVectorBasis())
     dof_ip = PHX::MDField<ScalarT>(
                evalName,
     	        pointRule->dl_vector);
  else
  { TEUCHOS_ASSERT(false); }

  this->addEvaluatedField(dof_ip);
  this->addDependentField(dof_basis);

  if(requires_orientation) {
     dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(fieldName+" Orientation", 
	                                                basis->functional);
     
     this->addDependentField(dof_orientation);
  }

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  MDFieldArrayFactory af_bv(basis->name()+"_"+pointRule->getName()+"_");
  basisValues.setupArrays(layout,af_bv);

  // the field manager will allocate all of these field

  this->addDependentField(basisValues.basis_ref);      
  this->addDependentField(basisValues.basis);           

  if(basis->getElementSpace()==PureBasis::HGRAD) {
    this->addDependentField(basisValues.grad_basis_ref);   
    this->addDependentField(basisValues.grad_basis);        
  }

  if(basis->getElementSpace()==PureBasis::HCURL) {
    this->addDependentField(basisValues.curl_basis_ref);     
    this->addDependentField(basisValues.curl_basis);          
  }
  
  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOF_PointValues,sd,fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  if(requires_orientation)
     this->utils.setFieldData(dof_orientation,fm);

  // setup the pointers for the basis values data structure
  this->utils.setFieldData(basisValues.basis_ref,fm);      
  this->utils.setFieldData(basisValues.basis,fm);           

  if(basis->getElementSpace()==PureBasis::HGRAD) {
    this->utils.setFieldData(basisValues.grad_basis_ref,fm);   
    this->utils.setFieldData(basisValues.grad_basis,fm);        
  }

  if(basis->getElementSpace()==PureBasis::HCURL) {
    this->utils.setFieldData(basisValues.curl_basis_ref,fm);     
    this->utils.setFieldData(basisValues.curl_basis,fm);          
  }
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOF_PointValues,workset)
{ 
  evaluateDOF(dof_basis,dof_ip,dof_orientation,requires_orientation,workset.num_cells,basisValues);
}

}

#endif
