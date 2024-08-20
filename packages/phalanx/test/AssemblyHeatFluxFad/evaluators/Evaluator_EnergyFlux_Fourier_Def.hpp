// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//**********************************************************************
template<typename EvalT, typename Traits> 
Fourier<EvalT, Traits>::
Fourier(const Teuchos::ParameterList& p) :
  flux("Energy_Flux", 
       p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout")),
  density("Density", 
	  p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout")),
  dc("Thermal Conductivity", 
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
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& /* fm */)
{  
  num_qp = static_cast<PHX::index_size_type>(flux.extent(1));
  num_dim = static_cast<PHX::index_size_type>(flux.extent(2));
}
//*********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void Fourier<EvalT, Traits>::operator () (const int i) const
{
  for (PHX::index_size_type qp = 0; qp < num_qp; ++qp)
    for (PHX::index_size_type dim = 0; dim < num_dim; ++dim)
      flux(i,qp,dim) =
	- density(i,qp) * dc(i,qp) * grad_temp(i,qp,dim);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Fourier<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  Kokkos::parallel_for(d.num_cells, *this);
}

//**********************************************************************
