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


#include "Teuchos_ParameterList.hpp"

//**********************************************************************
PHX_EVALUATOR_CTOR(Constant,p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);
  
  std::string n = "Constant Provider: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Constant,data,fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);
  
  if (constant.rank() == 1) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
	   constant(i) = value;
  }
  else if (constant.rank() == 2) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	constant(i,j) = value;
  }
  else if (constant.rank() == 3) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	for (PHX::Device::size_type k = 0; k < constant.dimension(2); ++k)
	  constant(i,j,k) = value;
  }
  else if (constant.rank() == 4) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	for (PHX::Device::size_type k = 0; k < constant.dimension(2); ++k)
	  for (PHX::Device::size_type l = 0; l < constant.dimension(3); ++l)
	    constant(i,j,k,l) = value;
  }
  else if (constant.rank() == 5) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	for (PHX::Device::size_type k = 0; k < constant.dimension(2); ++k)
	  for (PHX::Device::size_type l = 0; l < constant.dimension(3); ++l)
	    for (PHX::Device::size_type m = 0; m < constant.dimension(4); ++m)
	      constant(i,j,k,l,m) = value;
  }
  else if (constant.rank() == 6) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	for (PHX::Device::size_type k = 0; k < constant.dimension(2); ++k)
	  for (PHX::Device::size_type l = 0; l < constant.dimension(3); ++l)
	    for (PHX::Device::size_type m = 0; m < constant.dimension(4); ++m)
	      for (PHX::Device::size_type n = 0; n < constant.dimension(5); ++n)
		constant(i,j,k,l,m,n) = value;
  }
  else if (constant.rank() == 7) {
    for (PHX::Device::size_type i = 0; i < constant.dimension(0); ++i)
      for (PHX::Device::size_type j = 0; j < constant.dimension(1); ++j)
	for (PHX::Device::size_type k = 0; k < constant.dimension(2); ++k)
	  for (PHX::Device::size_type l = 0; l < constant.dimension(3); ++l)
	    for (PHX::Device::size_type m = 0; m < constant.dimension(4); ++m)
	      for (PHX::Device::size_type n = 0; n < constant.dimension(5); ++n)
		for (PHX::Device::size_type o = 0; o < constant.dimension(6); ++o)
		  constant(i,j,k,l,m,n,o) = value;
  }
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Constant,d)
{ }

//**********************************************************************
