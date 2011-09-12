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
PHX_POST_REGISTRATION_SETUP(NonlinearSource,data,fm)
{
  this->utils.setFieldData(source, fm);
  this->utils.setFieldData(density, fm);
  this->utils.setFieldData(temp, fm);

  cell_data_size = source.fieldTag().dataLayout().size() /
    source.fieldTag().dataLayout().dimension(0) ;
}

//**********************************************************************
PHX_EVALUATE_FIELDS(NonlinearSource,d)
{ 
  std::size_t size = d.size() * cell_data_size;
  
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
