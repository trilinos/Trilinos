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

//**********************************************************************
PHX_EVALUATOR_CTOR(Fourier,p) :
  flux("Energy_Flux", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  dc("Diffusion Coefficient", 
     p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  grad_temp("Temperature Gradient", 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{ 
  this->addEvaluatedField(flux);
  this->addDependentField(density);
  this->addDependentField(dc);
  this->addDependentField(grad_temp);

  this->setName("Fourier");
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Fourier,data,fm)
{
  this->utils.setFieldData(flux,fm);
  this->utils.setFieldData(density,fm);
  this->utils.setFieldData(dc,fm);
  this->utils.setFieldData(grad_temp,fm);

  cell_data_size = flux.fieldTag().dataLayout().size() / 
    flux.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Fourier,d)
{ 
  std::size_t size = d.size() * cell_data_size;

  for (std::size_t i = 0; i < size; ++i)
    flux[i] = - density[i] * dc[i] * grad_temp[i];
}

//**********************************************************************
