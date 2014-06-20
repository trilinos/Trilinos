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
template <typename ScalarT,typename ArrayT>
inline void evaluateDOF_withSens(PHX::MDField<ScalarT,Cell,Point> & dof_basis,PHX::MDField<ScalarT> & dof_ip,
                                 PHX::MDField<ScalarT,Cell,BASIS> & dof_orientation,
                                 bool is_vector_basis,
                                 int num_cells,
                                 ArrayT & basis)
{ 

  if(num_cells>0) {
    if(is_vector_basis) {
      int numCells  = basis.dimension(0);
      int numFields = basis.dimension(1);
      int numPoints = basis.dimension(2);
      int spaceDim  = basis.dimension(3);

      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d=0; d<spaceDim; d++) {
            // first initialize to the right thing (prevents over writing with 0)
            // then loop over one less basis function
            ScalarT & val = dof_ip(cell,pt,d);
            val = dof_basis(cell, 0) * basis(cell, 0, pt, d);
            for (int bf=1; bf<numFields; bf++)
              val += dof_basis(cell, bf) * basis(cell, bf, pt, d);
          }
        }
      } // for numCells

    }
    else { // no orientation needed
      // Zero out arrays (intrepid does a sum! 1/17/2012)
      for (int i = 0; i < dof_ip.size(); ++i)
        dof_ip[i] = 0.0;

      Intrepid::FunctionSpaceTools::
        evaluate<ScalarT>(dof_ip,dof_basis,basis);
    }
  }
}

template <typename ScalarT,typename ArrayT>
inline void evaluateDOF_fastSens(PHX::MDField<ScalarT,Cell,Point> & dof_basis,PHX::MDField<ScalarT> & dof_ip,
                                 PHX::MDField<ScalarT,Cell,BASIS> & dof_orientation,
                                 bool is_vector_basis,
                                 int num_cells,
                                 const std::vector<int> & offsets,
                                 ArrayT & basis)
{ 

  if(num_cells>0) {
    if(is_vector_basis) {
      int numCells  = basis.dimension(0);
      int numFields = basis.dimension(1);
      int numPoints = basis.dimension(2);
      int spaceDim  = basis.dimension(3);

      int fadSize = dof_basis(0,0).size(); // this is supposed to be fast
                                           // so assume that everything is the
                                           // same size!
      for (int cell=0; cell<numCells; cell++) {
        for (int pt=0; pt<numPoints; pt++) {
          for (int d=0; d<spaceDim; d++) {

            // first initialize to the right thing (prevents over writing with 0)
            // then loop over one less basis function
            ScalarT & val = dof_ip(cell,pt,d);

            // This is a possible issue if you need sensitivity to coordinates (you will need to
            // change basis and then use the product rule!)
            val = ScalarT(fadSize,dof_basis(cell, 0).val() * basis(cell, 0, pt, d).val());
            val.fastAccessDx(offsets[0]) = dof_basis(cell, 0).fastAccessDx(offsets[0]) * basis(cell, 0, pt, d).val();

            for (int bf=1; bf<numFields; bf++) {
              val.val() += dof_basis(cell, bf).val() * basis(cell, bf, pt, d).val();
              val.fastAccessDx(offsets[bf]) += dof_basis(cell, bf).fastAccessDx(offsets[bf]) * basis(cell, bf, pt, d).val();
            }
          }
        }
      } // for numCells

    }
    else { // no orientation needed
      // Zero out arrays (intrepid does a sum! 1/17/2012)
      for (int i = 0; i < dof_ip.size(); ++i)
        dof_ip[i] = 0.0;

      Intrepid::FunctionSpaceTools::
        evaluate<ScalarT>(dof_ip,dof_basis,basis);
    }
  }
}

}

//**********************************************************************
//* DOF_PointValues evaluator
//**********************************************************************

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename Traits>                   
DOF_PointValues<EvalT, Traits>::
DOF_PointValues(const Teuchos::ParameterList & p)
{
  const std::string fieldName = p.get<std::string>("Name");
  basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  Teuchos::RCP<const PointRule> pointRule = p.get< Teuchos::RCP<const PointRule> >("Point Rule");
  is_vector_basis = basis->isVectorBasis();

  std::string evalName = fieldName+"_"+pointRule->getName();
  if(p.isType<bool>("Use DOF Name")) {
    if(p.get<bool>("Use DOF Name"))
      evalName = fieldName;
  }

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

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  MDFieldArrayFactory af_bv(basis->name()+"_"+pointRule->getName()+"_");
  basisValues.setupArrays(layout,af_bv,false);

  // the field manager will allocate all of these field

  this->addDependentField(basisValues.basis_ref);      
  this->addDependentField(basisValues.basis);           

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name() + " ("+PHX::TypeString<EvalT>::value+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>                   
void DOF_PointValues<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  // setup the pointers for the basis values data structure
  this->utils.setFieldData(basisValues.basis_ref,fm);      
  this->utils.setFieldData(basisValues.basis,fm);           
}

//**********************************************************************
template<typename EvalT, typename Traits>                   
void DOF_PointValues<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  evaluateDOF_withSens(dof_basis,dof_ip,dof_orientation,is_vector_basis,workset.num_cells,basisValues.basis);
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename Traits>                   
DOF_PointValues<panzer::Traits::Jacobian, Traits>::
DOF_PointValues(const Teuchos::ParameterList & p)
{
  const std::string fieldName = p.get<std::string>("Name");
  basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  Teuchos::RCP<const PointRule> pointRule = p.get< Teuchos::RCP<const PointRule> >("Point Rule");
  is_vector_basis = basis->isVectorBasis();

  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");
    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix

  std::string evalName = fieldName+"_"+pointRule->getName();
  if(p.isType<bool>("Use DOF Name")) {
    if(p.get<bool>("Use DOF Name"))
      evalName = fieldName;
  }

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

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  MDFieldArrayFactory af_bv(basis->name()+"_"+pointRule->getName()+"_");
  basisValues.setupArrays(layout,af_bv,false);

  // the field manager will allocate all of these field

  this->addDependentField(basisValues.basis_ref);      
  this->addDependentField(basisValues.basis);           

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name() + " ("+PHX::TypeString<panzer::Traits::Jacobian>::value+")";
  this->setName(n);
}

//**********************************************************************
template<typename Traits>                   
void DOF_PointValues<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  // setup the pointers for the basis values data structure
  this->utils.setFieldData(basisValues.basis_ref,fm);      
  this->utils.setFieldData(basisValues.basis,fm);           
}

//**********************************************************************
template<typename Traits>                   
void DOF_PointValues<panzer::Traits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  if(accelerate_jacobian)
    evaluateDOF_fastSens(dof_basis,dof_ip,dof_orientation,is_vector_basis,workset.num_cells,offsets,basisValues.basis);
  else
    evaluateDOF_withSens(dof_basis,dof_ip,dof_orientation,is_vector_basis,workset.num_cells,basisValues.basis);
}

}

#endif
