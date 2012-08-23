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

#ifndef PANZER_DOF_POINT_FIELD_DECL_HPP
#define PANZER_DOF_POINT_FIELD_DECL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Intrepid_FunctionSpaceTools.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::initialize(const std::string & fieldName,
                                                const PureBasis & fieldBasis,
                                                const std::string & coordinateName,
                                                const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                                                const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                                                const std::string & postfixFieldName)
{
  intrepidBasis = fieldBasis.getIntrepidBasis();

  int cellCount = fieldBasis.functional->dimension(0);
  int coeffCount = fieldBasis.functional->dimension(1);
  int pointCount = coordLayout->dimension(0);
  int dimCount = coordLayout->dimension(1);

  Teuchos::RCP<PHX::DataLayout> basisLayout = fieldBasis.functional;

  coordinates = PHX::MDField<ScalarT,Point,Dim>(coordinateName,coordLayout);
  dof_coeff = PHX::MDField<ScalarT,Cell,Point>(fieldName,basisLayout);
  dof_field = PHX::MDField<ScalarT,Cell,Point>(fieldName+postfixFieldName,quadLayout);

  this->addDependentField(coordinates);
  this->addDependentField(dof_coeff);
  this->addEvaluatedField(dof_field);

  // build data storage for temporary conversion
  basisRef    = Intrepid::FieldContainer<double>(coeffCount,pointCount);
  basis       = Intrepid::FieldContainer<double>(cellCount,coeffCount,pointCount);
  intrpCoords = Intrepid::FieldContainer<double>(pointCount,dimCount);
  
  std::string n = "DOF_PointField: " + dof_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::postRegistrationSetup(typename TraitsT::SetupData d,
			                                  PHX::FieldManager<TraitsT>& fm)
{
  this->utils.setFieldData(coordinates,fm);
  this->utils.setFieldData(dof_coeff,fm);
  this->utils.setFieldData(dof_field,fm);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void DOF_PointField<EvalT,TraitsT>::evaluateFields(typename TraitsT::EvalData workset)
{ 
  // Zero out arrays (intrepid does a sum! 1/17/2012)
  for (int i = 0; i < dof_field.size(); ++i)
    dof_field[i] = 0.0;

  // copy coordinates
  for (int i = 0; i < coordinates.size(); ++i)
    intrpCoords[i] = Sacado::ScalarValue<ScalarT>::eval(coordinates[i]);

  if(workset.num_cells>0) {
    // evaluate at reference points
    intrepidBasis->getValues(basisRef, intrpCoords, Intrepid::OPERATOR_VALUE);

    // transfer reference basis values to physical frame values
    Intrepid::FunctionSpaceTools::
      HGRADtransformVALUE<ScalarT>(basis,
				  basisRef);

    // evaluate function at specified points
    Intrepid::FunctionSpaceTools::
      evaluate<ScalarT>(dof_field,dof_coeff,basis);
  }
}

//**********************************************************************

}

#endif
