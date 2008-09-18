// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP

//**********************************************************************
PHX_EVALUATOR_CTOR(NonlinearSource,p) :
  source("Nonlinear Source", 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"))
{ 
  this->addEvaluatedField(source);
  this->addDependentField(density);
  this->addDependentField(temp);

  this->setName("NonlinearSource");
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(NonlinearSource,fm)
{
  this->utils.setFieldData(source, fm);
  this->utils.setFieldData(density, fm);
  this->utils.setFieldData(temp, fm);

  data_layout_size = source.fieldTag().dataLayout().size();
}

//**********************************************************************
PHX_EVALUATE_FIELDS(NonlinearSource,d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    source[i] = density[i] * temp[i] * temp[i];
}

//**********************************************************************
PHX_PRE_EVALUATE_FIELDS(NonlinearSource,d)
{ 
  std::cout << "In NonlinearSource Pre Op" << std::endl;
}

//**********************************************************************
PHX_POST_EVALUATE_FIELDS(NonlinearSource,d)
{ 
  std::cout << "In NonlinearSource Post Op" << std::endl;
}

//**********************************************************************

#endif
