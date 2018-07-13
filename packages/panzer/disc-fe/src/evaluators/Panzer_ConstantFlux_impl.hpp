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

#ifndef PANZER_CONSTANT_FLUX_IMPL_HPP
#define PANZER_CONSTANT_FLUX_IMPL_HPP

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
ConstantFlux<EvalT, Traits>::
ConstantFlux(
  const Teuchos::ParameterList& p) :
  flux( p.get<std::string>("Flux Field Name"), 
	p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  const Teuchos::ParameterList& flux_values = p.sublist("Flux Values");

  for (Teuchos::ParameterList::ConstIterator i = flux_values.begin(); i != flux_values.end(); ++i)
    values.push_back(Teuchos::getValue<double>(i->second));

  this->addEvaluatedField(flux);
  
  std::string n = "Constant: " + flux.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantFlux<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  using namespace PHX;
  this->utils.setFieldData(flux,fm);

  TEUCHOS_ASSERT(static_cast<std::size_t>(flux.extent(2)) == values.size());

  for (int cell = 0; cell < flux.extent_int(0); ++cell)
    for (int ip = 0; ip < flux.extent_int(1); ++ip)
      for (int dim = 0; dim < flux.extent_int(2); ++dim)
	flux(cell,ip,dim) = values[dim];
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantFlux<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* d */)
{ }

//**********************************************************************

}

#endif
