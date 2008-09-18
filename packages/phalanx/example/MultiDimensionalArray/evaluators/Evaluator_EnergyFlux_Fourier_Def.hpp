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
template<typename EvalT, typename Traits> 
Fourier<EvalT, Traits>::
Fourier(const Teuchos::ParameterList& p) :
  flux("Energy_Flux", 
       p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout")),
  density("Density", 
	  p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout")),
  dc("Heat Capacity", 
     p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout") ),
  grad_temp("Temperature Gradient", 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout") )
{ 
  this->addEvaluatedField(flux);
  this->addDependentField(density);
  this->addDependentField(dc);
  this->addDependentField(grad_temp);

  this->setName("Fourier");
}

//**********************************************************************
template<typename EvalT, typename Traits> 
void Fourier<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(flux,fm);
  this->utils.setFieldData(density,fm);
  this->utils.setFieldData(dc,fm);
  this->utils.setFieldData(grad_temp,fm);

  num_qp = flux.dimension(1);
  num_dim = flux.dimension(2);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Fourier<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t num_cells = d.num_cells;

  for (std::size_t cell = 0; cell < num_cells; ++cell)
    for (std::size_t qp = 0; qp < num_qp; ++qp)
      for (std::size_t dim = 0; dim < num_dim; ++dim)
	flux(cell,qp,dim) = 
	  - density(cell,qp) * dc(cell,qp) * grad_temp(cell,qp,dim);
}

//**********************************************************************
