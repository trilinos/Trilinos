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

#ifndef USER_APP_CONSTANT_MODEL_T_HPP
#define USER_APP_CONSTANT_MODEL_T_HPP

//**********************************************************************
PHX_EVALUATOR_CTOR_NAMESPACE(user_app,ConstantModel,p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);
  
  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(user_app::ConstantModel,worksets,fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

//**********************************************************************
PHX_EVALUATE_FIELDS(user_app::ConstantModel,d)
{ }

//**********************************************************************


#ifdef HAVE_STOKHOS

namespace user_app {

template <typename Traits>
ConstantModel<typename Traits::SGResidual,Traits>::ConstantModel(const Teuchos::ParameterList & p) :
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  if(p.isParameter("Expansion")) {
     expansion = p.get<Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >("Expansion");
     value.reset(expansion);
     value.copyForWrite();
     value.fastAccessCoeff(0) = p.get<double>("Value");
     value.fastAccessCoeff(1) = p.get<double>("UQ");
  }
  else {
     value = p.get<double>("Value");
  }

  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

template <typename Traits>
void ConstantModel<typename Traits::SGResidual,Traits>::postRegistrationSetup(typename Traits::SetupData d,
                                                                              PHX::FieldManager<Traits>& fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

template <typename Traits>
ConstantModel<typename Traits::SGJacobian,Traits>::ConstantModel(const Teuchos::ParameterList & p) :
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  if(p.isParameter("Expansion")) {
     expansion = p.get<Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >("Expansion");
     value.val().reset(expansion);
     value.val().copyForWrite();
     value.val().fastAccessCoeff(0) = p.get<double>("Value");
     value.val().fastAccessCoeff(1) = p.get<double>("UQ");
  }
  else {
     value.val() = p.get<double>("Value");
  }

  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

template <typename Traits>
void ConstantModel<typename Traits::SGJacobian,Traits>::postRegistrationSetup(typename Traits::SetupData d,
                                                                              PHX::FieldManager<Traits>& fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

}

#endif

#endif
