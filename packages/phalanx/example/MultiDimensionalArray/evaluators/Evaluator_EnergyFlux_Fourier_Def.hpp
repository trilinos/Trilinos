// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
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
